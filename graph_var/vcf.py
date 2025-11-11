"""
VCF file writing functionality for pangenome graphs.
"""
from typing import Optional
from tqdm import tqdm
from .utils import (
    edge_complement,
    nearly_identical_alleles,
    _node_recover,
    log_action,
)


def write_vcf_from_graph(
    graph: "PangenomeGraph",
    reference_edges: dict[tuple[str, str], tuple[str, str]],
    gfa_path: Optional[str],
    vcf_filename: str,
    chr_name: str,
    exclude_terminus: bool = True,
    size_threshold: int = None,
    check_degenerate: bool = False,
    log_path: str = None
) -> None:
    """
    Writes the variant call format (vcf) file from a PangenomeGraph.
    
    :param graph: PangenomeGraph instance
    :param gfa_path: the .gfa file, from which walks are read, or None to skip writing genotypes
    :param vcf_filename: the output vcf file path
    :param chr_name: the chromosome name in the first column of output vcf file
    :param exclude_terminus: whether to exclude terminus nodes
    :param size_threshold: the truncation length of ref and alt sequence
    :param check_degenerate: whether to exclude variants whose ref and alt alleles are identical
    :param log_path: path to log file
    :return:
    """
    if log_path:
        log_action(log_path, f"Start generating vcf")
    
    # Generate VCF metadata
    meta_info = _generate_vcf_metadata(chr_name)
    
    # Get allele counts and genotypes
    allele_count_dict = graph.allele_count()
    if gfa_path:
        sample_to_genotype = graph.genotypes_from_gfa(gfa_path, exclude_terminus)
    else:
        sample_to_genotype = {}
    # Preserve order from GFA file (dict insertion order in Python 3.7+)
    sample_ids = list(sample_to_genotype.keys())
    header_names = list(graph.vcf_attribute_names) + sample_ids

    print("Writing vcf file")
    with open(vcf_filename, 'w') as file:
        file.write(meta_info)
        file.write('#' + '\t'.join(header_names) + '\n')

        for idx, (u, v) in tqdm(enumerate(graph.sorted_variant_edges(exclude_terminus=exclude_terminus))):
            representative_variant_edge = (u, v)
            reference_edge = reference_edges[representative_variant_edge]

            if graph.direction(u) == -1 and graph.direction(v) == -1:
                u, v = edge_complement((u, v))
            edge = (u, v)

            ref_allele, alt_allele, last_letter_of_branch_point, branch_point = graph.ref_alt_alleles(edge)
            VT = graph.identify_variant_type(edge, ref_allele, alt_allele)

            motif = graph.annotate_repeat_motif(representative_variant_edge,
                                               ref_allele=ref_allele,
                                               alt_allele=alt_allele,
                                               branch_point=branch_point)
            motif = '.' if motif is None else motif

            if check_degenerate:
                if ref_allele == alt_allele:
                    continue

            prepend_letter_to_alleles = (len(ref_allele) == 0 or len(alt_allele) == 0)

            new_ref = '.'
            ref_allele_on_forward_reference_path = graph.direction(edge[0]) == 1 and graph.on_reference_path(edge)
            if not ref_allele_on_forward_reference_path:
                new_ref = ref_allele
                ref_allele = '.'
                prepend_letter_to_alleles = False
            
            if prepend_letter_to_alleles:
                ref_allele = last_letter_of_branch_point + ref_allele
                alt_allele = last_letter_of_branch_point + alt_allele
            
            if size_threshold:
                ref_allele = ref_allele[:size_threshold]
                alt_allele = alt_allele[:size_threshold]

            edge_vcf_position = graph.get_vcf_position(edge, prepend_letter_to_alleles)

            # Build genotype records for this edge
            genotype_records = _build_genotype_record(
                variant_edge=representative_variant_edge,
                reference_edge=reference_edge,
                sample_to_genotype=sample_to_genotype,
                is_inversion=graph.is_inversion(edge),
                sample_order=sample_ids
            )
            
            # Build VCF record
            allele_data_list = _build_vcf_record(
                chr_name=chr_name,
                edge_vcf_position=edge_vcf_position,
                edge=edge,
                ref_allele=ref_allele,
                alt_allele=alt_allele,
                sample_ids=sample_ids,
                genotype_records=genotype_records
            )

            RC = allele_count_dict[representative_variant_edge][0] if not graph.is_inversion(edge) else '.'
            AC = allele_count_dict[representative_variant_edge][1]
            AN = RC + AC if not graph.is_inversion(edge) else '.'

            # Get haplotype positions for this variant edge
            hap_positions = {hap_name: graph.get_vcf_position(edge, prepend_letter_to_alleles, hap_name)
                            for hap_name in graph.haplo_priorities.keys()}
                            
            # Format haplotype positions as "hap1:pos1,hap2:pos2,..." or "." if empty/root
            hp_str = ','.join([f'{hap}:{pos}' for hap, pos in sorted(hap_positions.items()) if pos is not None])
            
            INFO = _build_info_field(
                new_ref=new_ref,
                VT=VT,
                graph=graph,
                u=u,
                v=v,
                RC=RC,
                AC=AC,
                AN=AN,
                hp_str=hp_str,
                motif=motif,
                ref_allele=ref_allele,
                alt_allele=alt_allele
            )

            allele_data_list[7] = INFO

            file.write('\t'.join(allele_data_list) + '\n')
    
    if log_path:
        log_action(log_path, f"Writing vcf: {vcf_filename}")


def _generate_vcf_metadata(chr_name: str) -> str:
    """Generate VCF metadata header lines."""
    meta_info = f'##fileformat=VCFv4.2\n'
    meta_info += f'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    meta_info += f'##FORMAT=<ID=CR,Number=.,Type=Integer,Description="The reference allele count for each sample\'s haplotypes">\n'
    meta_info += f'##FORMAT=<ID=CA,Number=.,Type=Integer,Description="The alternative allele count for each sample\'s haplotypes">\n'
    meta_info += f'##INFO=<ID=NR,Number=1,Type=String,Description="Non-reference allele">\n'
    meta_info += f'##INFO=<ID=VT,Number=1,Type=String,Description="Variant type">\n'
    meta_info += f'##INFO=<ID=DR,Number=2,Type=Integer,Description="Distance from reference (variant edge\'s two nodes)">\n'
    meta_info += f'##INFO=<ID=RC,Number=1,Type=Integer,Description="The REF allele count">\n'
    meta_info += f'##INFO=<ID=AC,Number=A,Type=Integer,Description="The ALT allele count">\n'
    meta_info += f'##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">\n'
    meta_info += f'##INFO=<ID=PV,Number=2,Type=Integer,Description="Position of variant edge\'s two nodes">\n'
    meta_info += f'##INFO=<ID=HP,Number=.,Type=String,Description="Haplotype positions at reference tree edge (haplotype:position)">\n'
    meta_info += f'##INFO=<ID=TR_MOTIF,Number=1,Type=String,Description="Tandem repeat motif">\n'
    meta_info += f'##INFO=<ID=NIA,Number=1,Type=Integer,Description="Nearly identical alleles (1=yes, 0=no)">\n'
    meta_info += f'##contig=<ID={chr_name}>\n'
    return meta_info

def _build_genotype_record(variant_edge: tuple[str, str], 
                            reference_edge: tuple[str, str],
                            sample_to_genotype: dict[str, tuple],
                            is_inversion: bool,
                            sample_order: list[str]
                            ) -> list[str]:
    """Build the genotype record for a variant edge.
    
    Returns a list of genotype strings, one per sample in the order specified by sample_order.
    Format:
    - Haploid: "gt:cr:ca"
    - Diploid: "gt0|gt1:cr0,cr1:ca0,ca1"
    """
    result = []
    def _parse_gt(gt):
        return str(gt) if gt is not None else '.'
    
    for sample_name in sample_order:
        genotypes = sample_to_genotype[sample_name]
        
        if len(genotypes) == 1:
            # Haploid
            gt, cr, ca = genotypes[0].variant_record(variant_edge, reference_edge)
            if is_inversion:
                cr = '.'
            gt_str = _parse_gt(gt)
            result.append(f"{gt_str}:{cr}:{ca}")
        elif len(genotypes) == 2:
            # Diploid
            gt0, cr0, ca0 = genotypes[0].variant_record(variant_edge, reference_edge)
            gt1, cr1, ca1 = genotypes[1].variant_record(variant_edge, reference_edge)
            
            if is_inversion:
                cr0 = '.'
                cr1 = '.'
            
            gt0_str = _parse_gt(gt0)
            gt1_str = _parse_gt(gt1)
            
            result.append(f"{gt0_str}|{gt1_str}:{cr0},{cr1}:{ca0},{ca1}")
        else:
            raise ValueError(f"Sample {sample_name} has {len(genotypes)} haplotypes, expected 1 or 2")
    
    return result


def _build_vcf_record(
    chr_name: str,
    edge_vcf_position: int,
    edge: tuple,
    ref_allele: str,
    alt_allele: str,
    sample_ids: list,
    genotype_records: list
) -> list:
    """Build a VCF record (one line) as a list of fields."""
    allele_data_list = []
    # 'CHROM' 0
    allele_data_list.append(chr_name)
    # 'POS' 1
    allele_data_list.append(str(edge_vcf_position))
    # 'ID' 2
    allele_data_list.append(''.join(tuple(map(lambda x: _node_recover(x), edge))))
    # 'REF' 3
    allele_data_list.append(ref_allele if ref_allele else '.')
    # 'ALT' 4
    allele_data_list.append(alt_allele if alt_allele else '.')
    # 'QUAL' 5
    allele_data_list.append('60')
    # 'FILTER' 6
    allele_data_list.append('PASS')
    # 'INFO' 7
    allele_data_list.append(None)
    # 'FORMAT' 8
    allele_data_list.append('GT:CR:CA')

    # 'sample1', 'sample2', ... 9 - end
    allele_data_list.extend(genotype_records)

    return allele_data_list


def _build_info_field(
    new_ref: str,
    VT: str,
    graph,
    u: str,
    v: str,
    RC,
    AC,
    AN,
    hp_str: str,
    motif: str,
    ref_allele: str,
    alt_allele: str
) -> str:
    """Build the INFO field for a VCF record."""
    INFO = (f'NR={new_ref if new_ref else "."};'
            f'VT={VT};'
            f'DR={int(graph.nodes[u]["distance_from_reference"])},{int(graph.nodes[v]["distance_from_reference"])};'
            f'RC={RC};'
            f'AC={AC};'
            f'AN={AN};'
            f'PV={int(graph.nodes[u]["position"])},{int(graph.nodes[v]["position"])};'
            f'HP={hp_str};'
            f'TR_MOTIF={motif}')

    if nearly_identical_alleles(ref_allele, alt_allele):
        INFO += ';NIA=1'
    else:
        INFO += ';NIA=0'

    return INFO

"""
VCF file writing functionality for pangenome graphs.
"""
from typing import Optional, TYPE_CHECKING, Callable
from dataclasses import dataclass
from tqdm import tqdm
from .utils import (
    edge_complement,
    nearly_identical_alleles,
    _node_recover,
    log_action,
)

if TYPE_CHECKING:
    from .graph import PangenomeGraph


@dataclass
class VariantInfo:
    """Encapsulates all information about a variant for VCF writing."""
    chr_name: str
    edge_vcf_position: int
    ref_allele: str
    alt_allele: str
    # Node and edge data
    edge_data: dict
    node_u_data: dict
    node_v_data: dict
    branch_point_node_data: dict
    # Allele counts
    ref_allele_count: int
    alt_allele_count: int
    
    @classmethod
    def from_graph(cls, 
                   graph: "PangenomeGraph",
                   edge: tuple[str, str],
                   representative_variant_edge: tuple[str, str],
                   reference_edge: tuple[str, str],
                   chr_name: str,
                   size_threshold: Optional[int] = None) -> "VariantInfo":
        """Compute variant information from a graph and variant edge."""
        
        # Get ref/alt alleles and branch point
        ref_allele_raw, alt_allele_raw, last_letter_of_branch_point, branch_point = graph.ref_alt_alleles(edge)
        
        # Determine if we need to prepend letter to alleles
        prepend_letter_to_alleles = (len(ref_allele_raw) == 0 or len(alt_allele_raw) == 0)
        
        # Check if ref allele is on forward reference path
        ref_allele_on_forward_reference_path = graph.direction(edge[0]) == 1 and graph.on_reference_path(edge)
        
        # Store the original ref allele before modification (for new_ref field)
        original_ref_allele = ref_allele_raw
        
        ref_allele = ref_allele_raw
        alt_allele = alt_allele_raw
        
        if not ref_allele_on_forward_reference_path:
            ref_allele = '.'
            prepend_letter_to_alleles = False
        
        # Prepend letter if needed
        if prepend_letter_to_alleles:
            ref_allele = last_letter_of_branch_point + ref_allele
            alt_allele = last_letter_of_branch_point + alt_allele
        
        # Truncate if size threshold specified
        if size_threshold:
            ref_allele = ref_allele[:size_threshold]
            alt_allele = alt_allele[:size_threshold]
        
        # Get VCF position
        edge_vcf_position = graph.get_vcf_position(edge, prepend_letter_to_alleles)
        
        # Get node and edge data
        u, v = edge
        edge_data = dict(graph.edges[edge])
        # Store raw alleles and metadata for eval functions
        edge_data['ref_allele_raw'] = ref_allele_raw
        edge_data['alt_allele_raw'] = alt_allele_raw
        edge_data['original_ref_allele'] = original_ref_allele
        edge_data['ref_on_forward_path'] = ref_allele_on_forward_reference_path
        edge_data['size_threshold'] = size_threshold
        
        node_u_data = dict(graph.nodes[u])
        node_v_data = dict(graph.nodes[v])
        branch_point_node_data = dict(graph.nodes[branch_point])
        
        # Get allele counts
        allele_count_dict = graph.allele_count()
        ref_allele_count, alt_allele_count = allele_count_dict[representative_variant_edge]
        
        return cls(
            chr_name=chr_name,
            edge_vcf_position=edge_vcf_position,
            ref_allele=ref_allele,
            alt_allele=alt_allele,
            edge_data=edge_data,
            node_u_data=node_u_data,
            node_v_data=node_v_data,
            branch_point_node_data=branch_point_node_data,
            ref_allele_count=ref_allele_count,
            alt_allele_count=alt_allele_count
        )


@dataclass
class InfoField:
    """VCF INFO field definition."""
    id: str
    number: str
    type: str
    description: str
    evaluate: Callable[[VariantInfo], str]
    
    def get_header(self) -> str:
        """Return the VCF header line for this INFO field."""
        return f'##INFO=<ID={self.id},Number={self.number},Type={self.type},Description="{self.description}">'


def get_default_info_fields() -> list[InfoField]:
    """Returns the standard set of INFO fields for VCF output."""
    
    # Define evaluation functions for each INFO field
    def _eval_non_reference_allele(variant_info: VariantInfo) -> str:
        # If ref is not on forward reference path, return the original ref allele
        if not variant_info.edge_data.get('ref_on_forward_path', True):
            return variant_info.edge_data.get('original_ref_allele', '.')
        return '.'

    def _eval_variant_type(variant_info: VariantInfo) -> str:
        # Use raw alleles to determine variant type
        ref = variant_info.edge_data.get('ref_allele_raw', '')
        alt = variant_info.edge_data.get('alt_allele_raw', '')
        
        ref_len = len(ref)
        alt_len = len(alt)
        
        if ref_len == alt_len:
            if ref_len == 1:
                return 'SNP'
            else:
                return 'MNP'
        elif ref_len == 0:
            return 'INS'
        elif alt_len == 0:
            return 'DEL'
        elif ref_len < alt_len:
            return 'INS'
        else:
            return 'DEL'

    def _eval_distance_from_reference(variant_info: VariantInfo) -> str:
        dr_u = int(variant_info.node_u_data["distance_from_reference"])
        dr_v = int(variant_info.node_v_data["distance_from_reference"])
        return f"{dr_u},{dr_v}"

    def _eval_ref_allele_count(variant_info: VariantInfo) -> str:
        # Check if this is an inversion
        is_inversion = variant_info.edge_data.get('is_inversion', False)
        return '.' if is_inversion else str(variant_info.ref_allele_count)

    def _eval_alt_allele_count(variant_info: VariantInfo) -> str:
        return str(variant_info.alt_allele_count)

    def _eval_total_allele_count(variant_info: VariantInfo) -> str:
        an = variant_info.ref_allele_count + variant_info.alt_allele_count
        # Check if this is an inversion
        is_inversion = variant_info.edge_data.get('is_inversion', False)
        return '.' if is_inversion else str(an)

    def _eval_position_of_variant(variant_info: VariantInfo) -> str:
        pos_u = int(variant_info.node_u_data["position"])
        pos_v = int(variant_info.node_v_data["position"])
        return f"{pos_u},{pos_v}"

    def _eval_haplotype_positions(variant_info: VariantInfo) -> str:
        # Compute haplotype positions using branch point position
        # Logic: position = branch_point_position + 1 - int(prepend_letter_to_alleles)
        # Haplotype positions are stored as numeric fields in branch_point_node_data
        # Skip known non-position fields
        skip_fields = {'direction', 'sequence', 'position', 'right_position', 'distance_from_reference', 'on_reference_path'}
        
        branch_point_positions = {}
        for field_name, pos in variant_info.branch_point_node_data.items():
            if field_name not in skip_fields and pos is not None:
                # This should be a haplotype-specific position
                try:
                    # Determine prepend_letter based on raw alleles
                    ref_raw = variant_info.edge_data.get('ref_allele_raw', '')
                    alt_raw = variant_info.edge_data.get('alt_allele_raw', '')
                    prepend = (len(ref_raw) == 0 or len(alt_raw) == 0) and variant_info.edge_data.get('ref_on_forward_path', True)
                    vcf_pos = int(pos) + 1 - int(prepend)
                    branch_point_positions[field_name] = vcf_pos
                except (ValueError, TypeError):
                    # Skip non-numeric fields
                    pass
        
        hp_str = ','.join([f'{hap}:{pos}' for hap, pos in sorted(branch_point_positions.items()) if pos is not None])
        return hp_str if hp_str else '.'

    def _eval_tandem_repeat_motif(variant_info: VariantInfo) -> str:
        # Motif is stored in edge_data from graph
        motif = variant_info.edge_data.get('motif')
        return motif if motif is not None else '.'

    def _eval_nearly_identical_alleles(variant_info: VariantInfo) -> str:
        nia = int(nearly_identical_alleles(variant_info.ref_allele, variant_info.alt_allele))
        return str(nia)
    
    # Return list of InfoField instances
    return [
        InfoField("NR", "1", "String", "Non-reference allele", _eval_non_reference_allele),
        InfoField("VT", "1", "String", "Variant type", _eval_variant_type),
        InfoField("DR", "2", "Integer", "Distance from reference (variant edge's two nodes)", _eval_distance_from_reference),
        InfoField("RC", "1", "Integer", "The REF allele count", _eval_ref_allele_count),
        InfoField("AC", "A", "Integer", "The ALT allele count", _eval_alt_allele_count),
        InfoField("AN", "1", "Integer", "Total number of alleles in called genotypes", _eval_total_allele_count),
        InfoField("PV", "2", "Integer", "Position of variant edge's two nodes", _eval_position_of_variant),
        InfoField("HP", ".", "String", "Haplotype positions at reference tree edge (haplotype:position)", _eval_haplotype_positions),
        InfoField("TR_MOTIF", "1", "String", "Tandem repeat motif", _eval_tandem_repeat_motif),
        InfoField("NIA", "1", "Integer", "Nearly identical alleles (1=yes, 0=no)", _eval_nearly_identical_alleles)
    ]


def write_vcf_from_graph_legacy(
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


def write_vcf_from_graph(
    graph: "PangenomeGraph",
    reference_edges: dict[tuple[str, str], tuple[str, str]],
    gfa_path: Optional[str],
    vcf_filename: str,
    chr_name: str,
    exclude_terminus: bool = True,
    size_threshold: int = None,
    check_degenerate: bool = False,
    log_path: str = None,
    info_fields: Optional[list[InfoField]] = None
) -> None:
    """
    Refactored VCF writing function using class-based design.
    
    :param graph: PangenomeGraph instance
    :param reference_edges: Dictionary mapping variant edges to reference edges
    :param gfa_path: the .gfa file, from which walks are read, or None to skip writing genotypes
    :param vcf_filename: the output vcf file path
    :param chr_name: the chromosome name in the first column of output vcf file
    :param exclude_terminus: whether to exclude terminus nodes
    :param size_threshold: the truncation length of ref and alt sequence
    :param check_degenerate: whether to exclude variants whose ref and alt alleles are identical
    :param log_path: path to log file
    :param info_fields: List of InfoField instances to include (defaults to all standard fields)
    :return:
    """
    if log_path:
        log_action(log_path, f"Start generating vcf")
    
    # Default INFO fields if none specified
    if info_fields is None:
        info_fields = get_default_info_fields()
    
    # Generate VCF metadata with INFO field headers
    meta_info = '##fileformat=VCFv4.2\n'
    meta_info += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    meta_info += '##FORMAT=<ID=CR,Number=R,Type=Integer,Description="Number of reads supporting the REF allele">\n'
    meta_info += '##FORMAT=<ID=CA,Number=A,Type=Integer,Description="Number of reads supporting the ALT allele">\n'
    
    # Add INFO field headers
    for field in info_fields:
        meta_info += field.get_header() + '\n'
    
    meta_info += f'##contig=<ID={chr_name}>\n'
    
    # Get genotypes if GFA path provided
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
            
            # Create VariantInfo from graph
            variant_info = VariantInfo.from_graph(
                graph=graph,
                edge=edge,
                representative_variant_edge=representative_variant_edge,
                reference_edge=reference_edge,
                chr_name=chr_name,
                size_threshold=size_threshold
            )
            
            # Check for degenerate alleles if requested
            if check_degenerate:
                if variant_info.ref_allele == variant_info.alt_allele:
                    continue
            
            # Build genotype records for this edge if samples provided
            if sample_to_genotype:
                is_inversion = variant_info.edge_data.get('is_inversion', False)
                genotype_records = _build_genotype_record(
                    variant_edge=representative_variant_edge,
                    reference_edge=reference_edge,
                    sample_to_genotype=sample_to_genotype,
                    is_inversion=is_inversion,
                    sample_order=sample_ids
                )
            else:
                genotype_records = []
            
            # Build VCF record
            allele_data_list = _build_vcf_record(
                chr_name=variant_info.chr_name,
                edge_vcf_position=variant_info.edge_vcf_position,
                edge=edge,
                ref_allele=variant_info.ref_allele,
                alt_allele=variant_info.alt_allele,
                sample_ids=sample_ids,
                genotype_records=genotype_records
            )
            
            # Build INFO field using InfoField instances
            info_parts = []
            for field in info_fields:
                value = field.evaluate(variant_info)
                info_parts.append(f"{field.id}={value}")
            
            allele_data_list[7] = ';'.join(info_parts)
            
            file.write('\t'.join(allele_data_list) + '\n')
    
    if log_path:
        log_action(log_path, f"Writing vcf: {vcf_filename}")

"""
VCF file writing functionality for pangenome graphs.
"""
from typing import Optional, TYPE_CHECKING, Callable
from dataclasses import dataclass
from importlib.metadata import version, PackageNotFoundError
import logging
from .utils import (
    edge_complement,
    nearly_identical_alleles,
    node_recover,
    get_from_biedge_dict,
)
from .vcf_io import open_vcf_text_for_write

if TYPE_CHECKING:
    from .graph import PangenomeGraph

# Get package version
try:
    __version__ = version("pantree")
except PackageNotFoundError:
    __version__ = "unknown"


@dataclass
class _VariantData:
    """Encapsulates all information about a variant for VCF writing."""
    chr_name: str
    ref_allele_raw: str
    alt_allele_raw: str
    # Node and edge data
    edge_data: dict
    node_u_data: dict
    node_v_data: dict
    branch_point_node_data: dict
    # Allele counts (computed from genotype data)
    ref_allele_count: int
    alt_allele_count: int
    # Context for eval functions (parameters, etc.)
    context: dict

    @property
    def needs_prepend(self) -> bool:
        return (len(self.ref_allele_raw) == 0 or len(self.alt_allele_raw) == 0) \
            and self.node_v_data['on_reference_path']

    @property
    def is_degenerate(self) -> bool:
        """Check if ref and alt alleles are identical after VCF formatting."""
        return self.ref_allele_raw == self.alt_allele_raw

    def get_variant_position(self, position_field: str) -> int:
        branch_position = self.branch_point_node_data.get(position_field)
        if branch_position is None:
            return None
        return branch_position + 1 - int(self.needs_prepend)

    @staticmethod
    def compute_allele_counts(sample_to_genotype: dict,
                             variant_edge: tuple[str, str],
                             ref_edge_idx: int) -> tuple[int, int]:
        """Compute ref and alt allele counts from genotype data.

        Args:
            sample_to_genotype: Dict mapping sample names to tuples of Genotype objects
            variant_edge: The variant edge
            ref_edge_idx: The index of the reference edge

        Returns:
            Tuple of (ref_count, alt_count)
        """
        ref_count = 0
        alt_count = 0

        for genotypes in sample_to_genotype.values():
            for genotype in genotypes:
                # Get counts for this haplotype
                cr = genotype.get_ref_count(ref_edge_idx)
                ca = get_from_biedge_dict(genotype.alt_counts, variant_edge, 0)
                ref_count += cr
                alt_count += ca

        return ref_count, alt_count

    @classmethod
    def from_graph(cls,
                   graph: "PangenomeGraph",
                   edge: tuple[str, str],
                   reference_edge: tuple[str, str],
                   chr_name: str,
                   sample_to_genotype: dict,
                   size_threshold: float = float('inf'),
                   **kwargs) -> "_VariantData":
        """Compute variant information from a graph and variant edge."""

        # Get ref/alt alleles and branch point
        ref_allele_raw, alt_allele_raw, _, branch_point = graph.ref_alt_alleles(edge)

        # Check if ref allele is on forward reference path
        u, v = edge
        ref_allele_on_forward_reference_path = graph.direction(u) == 1 and graph.on_reference_path(edge)

        # Determine if we need to prepend letter to alleles
        prepend_letter_to_alleles = (len(ref_allele_raw) == 0 or len(alt_allele_raw) == 0) and ref_allele_on_forward_reference_path

        # Get node and edge data (store only what's in the graph)
        edge_data = dict(graph.edges[edge])
        node_u_data = dict(graph.nodes[u])
        node_v_data = dict(graph.nodes[v])
        branch_point_node_data = dict(graph.nodes[branch_point])

        # Compute and store repeat motif
        motif = graph.annotate_repeat_motif(edge, ref_allele_raw, alt_allele_raw, branch_point)
        edge_data['motif'] = motif

        # Compute allele counts from genotype data
        if sample_to_genotype:
            ref_edge_idx = int(graph.edges[reference_edge]['index'])
            ref_allele_count, alt_allele_count = cls.compute_allele_counts(
                sample_to_genotype, edge, ref_edge_idx
            )
        else:
            # No genotype data, set to zero
            ref_allele_count = 0
            alt_allele_count = 0

        # Create context dict for eval functions and VCF defaults
        # Store haplotype field names from graph's haplo_priorities
        haplotype_fields = set(graph.haplo_priorities.keys()) if hasattr(graph, 'haplo_priorities') else set()
        default_context = {
            'size_threshold': size_threshold,
            'haplotype_fields': haplotype_fields,
            'qual': '60',
            'filter_field': 'PASS',
            'format_field': 'GT:CR:CA'
        }
        # Merge with any user-provided context overrides
        context = {**default_context, **kwargs.get('context', {})}

        return cls(
            chr_name=chr_name,
            ref_allele_raw=ref_allele_raw,
            alt_allele_raw=alt_allele_raw,
            edge_data=edge_data,
            node_u_data=node_u_data,
            node_v_data=node_v_data,
            branch_point_node_data=branch_point_node_data,
            ref_allele_count=ref_allele_count,
            alt_allele_count=alt_allele_count,
            context=context,
        )

    @staticmethod
    def compute_vcf_alleles(ref_raw: str, alt_raw: str, last_letter: str,
                            ref_on_forward_path: bool, size_threshold: Optional[int]) -> tuple[str, str]:
        """Compute the final VCF ref and alt alleles from raw data."""
        ref_allele = ref_raw
        alt_allele = alt_raw

        if not ref_on_forward_path:
            ref_allele = '.'
        else:
            # Prepend letter if needed
            prepend = (len(ref_raw) == 0 or len(alt_raw) == 0)
            if prepend:
                ref_allele = last_letter + ref_allele
                alt_allele = last_letter + alt_allele

        # Truncate if size threshold specified (check for finite value)
        if size_threshold is not None and size_threshold != float('inf'):
            ref_allele = ref_allele[:int(size_threshold)]
            alt_allele = alt_allele[:int(size_threshold)]

        return ref_allele, alt_allele


@dataclass
class _VariantRecord:
    """VCF record built from _VariantData."""
    chr_name: str
    variant_id: str
    ref_allele: str
    alt_allele: str
    vcf_position: int
    qual: str
    filter_field: str
    info: str
    format_field: str
    genotype_records: list[str]

    @classmethod
    def from_variant_data(cls,
                         variant_info: _VariantData,
                         edge: tuple[str, str],
                         genotype_records: list[str],
                         info_fields: list["_InfoField"]) -> "_VariantRecord":
        """Build a VCF record from _VariantData."""
        # Compute VCF alleles
        ref_allele, alt_allele = _VariantData.compute_vcf_alleles(
            variant_info.ref_allele_raw, variant_info.alt_allele_raw,
            variant_info.branch_point_node_data.get('sequence', '')[-1] if variant_info.branch_point_node_data.get('sequence') else '',
            variant_info.node_u_data.get('direction', 1) == 1 and variant_info.node_v_data.get('on_reference_path', 1) == 1,
            variant_info.context.get('size_threshold')
        )

        # Build INFO field
        info_parts = []
        for field in info_fields:
            value = field.evaluate(variant_info)
            info_parts.append(f"{field.id}={value}")
        info = ';'.join(info_parts)

        # Build variant ID
        variant_id = ''.join(tuple(map(lambda x: node_recover(x), edge)))

        return cls(
            chr_name=variant_info.chr_name,
            vcf_position=variant_info.get_variant_position('position'),
            variant_id=variant_id,
            ref_allele=ref_allele if ref_allele else '.',
            alt_allele=alt_allele if alt_allele else '.',
            qual=variant_info.context.get('qual', '60'),
            filter_field=variant_info.context.get('filter_field', 'PASS'),
            info=info,
            format_field=variant_info.context.get('format_field', 'GT:CR:CA'),
            genotype_records=genotype_records
        )

    def to_vcf_line(self) -> str:
        """Convert record to VCF line."""
        fields = [
            self.chr_name,
            str(self.vcf_position),
            self.variant_id,
            self.ref_allele,
            self.alt_allele,
            self.qual,
            self.filter_field,
            self.info,
        ]
        if self.genotype_records:
            fields.append(self.format_field)
            fields.extend(self.genotype_records)
        return '\t'.join(fields) + '\n'

    def sort_key(self) -> tuple[str, int, int, str]:
        """Sort by VCF coordinates while preserving graph-order ties with UIDX."""
        uidx = -1
        for info_field in self.info.split(';'):
            key, _, value = info_field.partition('=')
            if key == 'UIDX':
                try:
                    uidx = int(value)
                except ValueError:
                    uidx = -1
                break
        return self.chr_name, int(self.vcf_position), uidx, self.variant_id


@dataclass
class _InfoField:
    """VCF INFO field definition."""
    id: str
    number: str
    type: str
    description: str
    evaluate: Callable[[_VariantData], str]

    def get_header(self) -> str:
        """Return the VCF header line for this INFO field."""
        return f'##INFO=<ID={self.id},Number={self.number},Type={self.type},Description="{self.description}">'


def _get_default_info_fields() -> list[_InfoField]:
    """Returns the standard set of INFO fields for VCF output."""

    info_fields = []

    # Helper to get raw alleles
    def _get_raw_alleles(variant_info: _VariantData) -> tuple[str, str]:
        """Get raw ref and alt alleles."""
        return variant_info.ref_allele_raw, variant_info.alt_allele_raw

    # Helper to compute VCF alleles from raw alleles
    def _compute_vcf_alleles(variant_info: _VariantData) -> tuple[str, str]:
        """Compute final VCF ref and alt alleles from raw data."""
        ref_raw, alt_raw = _get_raw_alleles(variant_info)
        ref_on_forward = _ref_on_forward_path(variant_info)
        last_letter = _get_last_letter_of_branch_point(variant_info)
        size_threshold = variant_info.context.get('size_threshold')

        ref_allele = ref_raw
        alt_allele = alt_raw

        if not ref_on_forward:
            ref_allele = '.'
        else:
            # Prepend letter if needed (when either allele is empty)
            if variant_info.needs_prepend:
                ref_allele = last_letter + ref_allele
                alt_allele = last_letter + alt_allele

        # Truncate if size threshold specified (check for finite value)
        if size_threshold is not None and size_threshold != float('inf'):
            ref_allele = ref_allele[:int(size_threshold)]
            alt_allele = alt_allele[:int(size_threshold)]

        return ref_allele, alt_allele

    # Helper to get last letter of branch point
    def _get_last_letter_of_branch_point(variant_info: _VariantData) -> str:
        """Get the last letter of the branch point sequence."""
        branch_seq = variant_info.branch_point_node_data.get('sequence', '')
        return branch_seq[-1] if branch_seq else 'N'

    # Helper to check if ref is on forward reference path
    def _ref_on_forward_path(variant_info: _VariantData) -> bool:
        """Check if ref allele is on forward reference path."""
        # Check direction of source node and if edge is on reference path
        direction = variant_info.node_u_data.get('direction', 1)
        on_ref_path = variant_info.node_v_data.get('on_reference_path', 1)
        return direction == 1 and on_ref_path == 1

    # Define evaluation functions and _InfoFields
    def _eval_non_reference_allele(variant_info: _VariantData) -> str:
        # If ref is not on forward reference path, return the raw ref allele
        if not _ref_on_forward_path(variant_info):
            ref_raw, _ = _get_raw_alleles(variant_info)
            return ref_raw if ref_raw else '.'
        return '.'
    info_fields.append(_InfoField("NR", "1", "String", "Non-reference allele", _eval_non_reference_allele))

    def _eval_variant_type(variant_info: _VariantData) -> str:
        if variant_info.node_u_data['direction'] != variant_info.node_v_data['direction']:
            return 'INV'

        if variant_info.edge_data['is_back_edge']:
            return 'DUP'

        # Use raw alleles to determine variant type
        ref, alt = _get_raw_alleles(variant_info)

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
        return 'REP'

    info_fields.append(_InfoField("VT", "1", "String", "Variant type", _eval_variant_type))

    def _eval_tree_position(variant_info: _VariantData) -> str:
        tp = variant_info.get_variant_position("tree_position")
        return '.' if tp is None else str(tp)
    info_fields.append(_InfoField("TP", "1", "Integer", "Tree position of the variant edge's branch point", _eval_tree_position))

    def _eval_ref_allele_count(variant_info: _VariantData) -> str:
        # Check if this is an inversion
        is_inversion = variant_info.edge_data.get('is_inversion', False)
        return '.' if is_inversion else str(variant_info.ref_allele_count)
    info_fields.append(_InfoField("RC", "1", "Integer", "The REF allele count", _eval_ref_allele_count))

    def _eval_alt_allele_count(variant_info: _VariantData) -> str:
        return str(variant_info.alt_allele_count)
    info_fields.append(_InfoField("AC", "A", "Integer", "The ALT allele count", _eval_alt_allele_count))

    def _eval_total_allele_count(variant_info: _VariantData) -> str:
        an = variant_info.ref_allele_count + variant_info.alt_allele_count
        # Check if this is an inversion
        is_inversion = variant_info.edge_data.get('is_inversion', False)
        return '.' if is_inversion else str(an)
    info_fields.append(_InfoField("AN", "1", "Integer", "Total number of alleles in called genotypes", _eval_total_allele_count))

    def _eval_haplotype_positions(variant_info: _VariantData) -> str:
        haplotype_fields = variant_info.context.get('haplotype_fields', set())
        positions = {}
        for field_name in haplotype_fields:
            positions[field_name] = variant_info.get_variant_position(field_name)

        hp_str = ','.join([f'{hap}:{pos}' for hap, pos in sorted(positions.items()) if pos is not None])
        return hp_str if hp_str else '.'
    info_fields.append(_InfoField("HP", ".", "String", "Haplotype positions at reference tree edge (haplotype:position)", _eval_haplotype_positions))

    def _eval_tandem_repeat_motif(variant_info: _VariantData) -> str:
        # Motif is stored in edge_data from graph
        motif = variant_info.edge_data.get('motif')
        return motif if motif is not None else '.'
    info_fields.append(_InfoField("TR_MOTIF", "1", "String", "Tandem repeat motif", _eval_tandem_repeat_motif))

    def _eval_nearly_identical_alleles(variant_info: _VariantData) -> str:
        ref_allele, alt_allele = _compute_vcf_alleles(variant_info)
        nia = int(nearly_identical_alleles(ref_allele, alt_allele))
        return str(nia)
    info_fields.append(_InfoField("NIA", "1", "Integer", "Nearly identical alleles (1=yes, 0=no)", _eval_nearly_identical_alleles))

    def _eval_u_index(variant_info: _VariantData) -> str:
        return str(variant_info.node_u_data['index'])
    info_fields.append(_InfoField("UIDX", "1", "Integer", "Index of node u", _eval_u_index))

    return info_fields


def _build_genotype_record(variant_edge: tuple[str, str],
                            reference_edge: tuple[str, str],
                            ref_edge_idx: int,
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
            gt, cr, ca = genotypes[0].variant_record(variant_edge, reference_edge, ref_edge_idx)
            if is_inversion:
                cr = '.'
            gt_str = _parse_gt(gt)
            result.append(f"{gt_str}:{cr}:{ca}")
        elif len(genotypes) == 2:
            # Diploid
            gt0, cr0, ca0 = genotypes[0].variant_record(variant_edge, reference_edge, ref_edge_idx)
            gt1, cr1, ca1 = genotypes[1].variant_record(variant_edge, reference_edge, ref_edge_idx)

            if is_inversion:
                cr0 = '.'
                cr1 = '.'

            gt0_str = _parse_gt(gt0)
            gt1_str = _parse_gt(gt1)

            result.append(f"{gt0_str}|{gt1_str}:{cr0},{cr1}:{ca0},{ca1}")
        else:
            raise ValueError(f"Sample {sample_name} has {len(genotypes)} haplotypes, expected 1 or 2")

    return result




def _build_vcf_header(chr_name: str,
                         info_fields: list[_InfoField],
                         sample_ids: list[str]) -> str:
    """Generate VCF header including metadata and column names.

    Args:
        chr_name: Chromosome name
        info_fields: List of INFO field definitions
        sample_ids: List of sample IDs for genotype columns

    Returns:
        Complete VCF header as string
    """
    # Generate VCF metadata with INFO field headers
    meta_info = '##fileformat=VCFv4.2\n'
    meta_info += f'##source=pantree v{__version__}\n'
    if sample_ids:
        meta_info += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype: 1 if ALT is present, 0 if absent, . if missing">\n'
        meta_info += '##FORMAT=<ID=CR,Number=R,Type=Integer,Description="Number of times visiting the REF allele">\n'
        meta_info += '##FORMAT=<ID=CA,Number=A,Type=Integer,Description="Number of times visiting the ALT allele">\n'

    # Add INFO field headers
    for field in info_fields:
        meta_info += field.get_header() + '\n'

    meta_info += f'##contig=<ID={chr_name}>\n'

    # Add column header line
    header_names = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    if sample_ids:
        header_names += ['FORMAT'] + sample_ids
    meta_info += '#' + '\t'.join(header_names) + '\n'

    return meta_info


def write_vcf_from_graph(
    graph: "PangenomeGraph",
    gfa_path: Optional[str],
    vcf_filename: str,
    chr_name: str,
    logger: logging.Logger,
    exclude_terminus: bool = True,
    size_threshold: float = float('inf'),
    check_degenerate: bool = False,
    info_fields: Optional[list[_InfoField]] = None,
    no_missingness: bool = False,
) -> None:
    """
    Refactored VCF writing function using class-based design.

    :param graph: PangenomeGraph instance
    :param gfa_path: the .gfa file, from which walks are read, or None to skip writing genotypes
    :param vcf_filename: the output vcf file path (use .vcf.gz extension for BGZF output)
    :param chr_name: the chromosome name in the first column of output vcf file
    :param logger: Logger instance for logging actions
    :param exclude_terminus: whether to exclude terminus nodes
    :param size_threshold: the truncation length of ref and alt sequence
    :param check_degenerate: whether to exclude variants whose ref and alt alleles are identical
    :param info_fields: List of _InfoField instances to include (defaults to all standard fields)
    :param no_missingness: whether to skip missingness computation for genotypes
    :return:
    """
    logger.info(f"Start generating vcf")

    if info_fields is None:
        info_fields = _get_default_info_fields()

    reference_edges = graph.get_reference_edges()

    if gfa_path:
        logger.info(f"Getting genotypes from GFA file: {gfa_path}")
        sample_to_genotype = graph.genotypes_from_gfa(gfa_path, exclude_terminus, skip_missing=no_missingness)
    else:
        sample_to_genotype = {}

    sample_ids = list(sample_to_genotype.keys())

    vcf_header = _build_vcf_header(chr_name, info_fields, sample_ids)

    logger.info(f"Writing vcf: {vcf_filename}")

    vcf_records = []
    for u, v in graph.sorted_variant_edges(exclude_terminus=exclude_terminus):
        reference_edge = reference_edges[(u, v)]

        if graph.direction(u) == -1 and graph.direction(v) == -1:
            u, v = edge_complement((u, v))
        edge = (u, v)

        # Create _VariantData from graph
        variant_info = _VariantData.from_graph(
            graph=graph,
            edge=edge,
            reference_edge=reference_edge,
            chr_name=chr_name,
            sample_to_genotype=sample_to_genotype,
            size_threshold=size_threshold
        )

        # Check for degenerate alleles if requested
        if check_degenerate and variant_info.is_degenerate:
            continue

        # Build genotype records for this edge if samples provided
        if sample_to_genotype:
            is_inversion = variant_info.edge_data.get('is_inversion', False)
            ref_edge_idx = int(graph.edges[reference_edge]['index'])
            genotype_records = _build_genotype_record(
                variant_edge=edge,
                reference_edge=reference_edge,
                ref_edge_idx=ref_edge_idx,
                sample_to_genotype=sample_to_genotype,
                is_inversion=is_inversion,
                sample_order=sample_ids
            )
        else:
            genotype_records = []

        # Build VCF record from variant data
        vcf_records.append(_VariantRecord.from_variant_data(
            variant_info=variant_info,
            edge=edge,
            genotype_records=genotype_records,
            info_fields=info_fields
        ))

    with open_vcf_text_for_write(vcf_filename) as file_handle:
        file_handle.write(vcf_header)
        for vcf_record in sorted(vcf_records, key=lambda record: record.sort_key()):
            file_handle.write(vcf_record.to_vcf_line())

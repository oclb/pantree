from __future__ import annotations
import polars as pl
from dataclasses import dataclass
from pantree.read_vcf import read_vcf_to_lazyframe

def overlap(a: Allele, b: Allele):
    ia = a.interval
    ib = b.interval
    start = max(ia[0], ib[0])
    end = min(ia[1], ib[1])
    return start, end

@dataclass
class Allele:
    sequence: list
    position: int

    def __getitem__(self, other: slice):
        return self.sequence[other.start - self.position : other.stop - self.position]

    def __len__(self):
        return len(self.sequence)

    @property
    def interval(self):
        return self.position, self.position + len(self.sequence)

    def is_disjoint(self, other: Allele):
        ovl = overlap(self, other)
        return ovl[1] < ovl[0]

    def matches(self, other: Allele):
        ovl = overlap(self, other)
        return self[ovl[0] : ovl[1]] == other[ovl[0] : ovl[1]]

    def prepend(self, other: Allele):
        if other.position >= self.position:
            return
        if self.is_disjoint(other):
            raise ValueError("Alleles are disjoint")
        ovl = overlap(self, other)
        other_len = len(other) - ovl[1] + ovl[0]
        self.position = other.position
        self.sequence = other.sequence[:other_len] + self.sequence

    def replace(self, ref: Allele, alt: Allele) -> bool:
        """
        Replace the ref region with alt in this allele.

        Returns:
            True if replacement succeeded, False if alleles don't match.
        """
        # Extend to cover ref if needed
        if ref.position < self.position:
            # Need to extend left with ref bases
            ovl = overlap(self, ref)
            if ovl[1] >= ovl[0]:  # Not disjoint or adjacent
                ref_len_to_prepend = ovl[0] - ref.position
                self.position = ref.position
                self.sequence = ref.sequence[:ref_len_to_prepend] + self.sequence

        # Validate that ref matches the existing sequence
        ovl = overlap(self, ref)
        if ovl[1] >= ovl[0]:  # Not disjoint or adjacent
            start_idx = ovl[0] - self.position
            end_idx = ovl[1] - self.position
            existing = self.sequence[start_idx:end_idx]
            expected = ref.sequence[:end_idx - start_idx]
            if existing != expected:
                return False

        # Now replace the ref region with alt
        start_idx = ovl[0] - self.position
        end_idx = ovl[1] - self.position
        self.sequence = self.sequence[:start_idx] + alt.sequence + self.sequence[end_idx:]
        return True

def _get_dup_ranges(df: pl.DataFrame) -> list[range]:
    """Get ranges covered by DUP variants (using TP, tree position)."""
    ranges = []
    for line in df.iter_rows(named=True):
        if line.get('VT') == 'DUP':
            ranges.append(range(line['TP'], line['TP'] + len(line['REF'])))
    return ranges


@dataclass
class ConsolidatedVariant:
    """A consolidated variant with status code."""
    ref: Allele
    alt: Allele
    status: str  # OK, INV, DUP, SKIP, ERR

def get_walk_variants(vcf: pl.LazyFrame) -> list[ConsolidatedVariant]:
    """
    Combine overlapping variants from a VCF into consolidated ref/alt sequences.

    Args:
        vcf: A Polars LazyFrame with columns UIDX (sort order), TP (position), REF, ALT, VT

    Returns:
        List of ConsolidatedVariant objects with status codes:
        - OK: Successfully consolidated
        - INV: Inversion variant (not aggregated with others)
        - DUP: Duplication variant (not aggregated with others)
        - SKIP: Skipped due to being in a DUP region
        - ERR: Allele mismatch detected (possibly from fragmented assembly)
    """
    df = vcf.sort("UIDX", descending=True).collect()
    dup_ranges: list[range] = _get_dup_ranges(df)

    results: list[ConsolidatedVariant] = []

    current_ref: Allele | None = None
    current_alt: Allele | None = None

    for line in df.iter_rows(named=True):
        ref = Allele(list(line["REF"]), line["TP"])
        alt = Allele(list(line["ALT"]), line["TP"])
        vt = line.get("VT", "")
        pos = line.get("POS")

        # Handle INV variants - add separately, never aggregate
        if vt == "INV":
            # Finalize current group first if exists
            if current_ref is not None:
                results.append(ConsolidatedVariant(current_ref, current_alt, "OK"))
                current_ref = None
                current_alt = None
            results.append(ConsolidatedVariant(ref, alt, "INV"))
            continue

        # Handle DUP variants - add separately, never aggregate
        if vt == "DUP":
            # Finalize current group first if exists
            if current_ref is not None:
                results.append(ConsolidatedVariant(current_ref, current_alt, "OK"))
                current_ref = None
                current_alt = None
            results.append(ConsolidatedVariant(ref, alt, "DUP"))
            continue

        # Check if this variant is in a DUP region
        in_dup_range = any(r.start <= line["TP"] < r.stop for r in dup_ranges)
        if in_dup_range:
            # Finalize current group first if exists
            if current_ref is not None:
                results.append(ConsolidatedVariant(current_ref, current_alt, "OK"))
                current_ref = None
                current_alt = None
            results.append(ConsolidatedVariant(ref, alt, "SKIP"))
            continue

        # Start a new group if none exists
        if current_ref is None:
            current_ref = ref
            current_alt = alt
            continue

        assert current_alt is not None

        # Check if disjoint - start new group
        if current_alt.is_disjoint(ref):
            results.append(ConsolidatedVariant(current_ref, current_alt, "OK"))
            current_ref = ref
            current_alt = alt
            continue

        # Try to aggregate
        current_ref.prepend(ref)
        success = current_alt.replace(ref, alt)

        if not success:
            # Allele mismatch - likely from fragmented assembly
            print(
                f"Warning: Allele mismatch at POS={pos}. "
                f"This may be caused by a fragmented assembly or by the presence of duplications or inversions. "
                f"Starting new variant group."
            )
            # Finalize current group as OK (it was fine until this variant)
            results.append(ConsolidatedVariant(current_ref, current_alt, "OK"))
            # Start new group with this variant marked as ERR
            results.append(ConsolidatedVariant(ref, alt, "ERR"))
            current_ref = None
            current_alt = None

    # Append the final group
    if current_ref is not None:
        results.append(ConsolidatedVariant(current_ref, current_alt, "OK"))

    return results


def process_haplotype_variants(
    vcf_path: str,
    sample_name: str,
    haplotype_number: int,
    output_vcf_path: str
):
    """
    Process a VCF file to produce consolidated ref/alt alleles for a specific haplotype.

    Args:
        vcf_path: Path to input VCF file
        sample_name: Name of the sample to analyze
        haplotype_number: Which haplotype to use (0 or 1 for diploid)
        output_vcf_path: Path to write output VCF file
    """
    # Read VCF
    vcf_df = read_vcf_to_lazyframe(vcf_path)

    # Get the GT column for this sample
    gt_col = f"{sample_name}_GT"

    # Filter to rows where this haplotype has the alt allele
    # GT format is "a|b" or "a/b", extract the appropriate allele
    vcf_filtered = vcf_df.with_columns(
        pl.col(gt_col).str.split_exact("|", 1).alias("gt_split")
    ).with_columns(
        pl.col("gt_split").struct.field(f"field_{haplotype_number}").alias("haplotype_gt")
    )

    # Build filter condition - always require haplotype_gt == "1"
    filter_cond = pl.col("haplotype_gt") == "1"
    vcf_filtered = vcf_filtered.filter(filter_cond)

    # Replace REF with NR when NR is not '.'
    # NR contains the actual reference allele for variants off the linear reference
    # Only do this if NR column exists
    if "NR" in vcf_filtered.collect_schema().names():
        vcf_filtered = vcf_filtered.with_columns(
            pl.when(pl.col("NR") != ".")
            .then(pl.col("NR"))
            .otherwise(pl.col("REF"))
            .alias("REF")
        )

    # Replace '.' with empty string in REF and ALT columns
    # '.' represents missing/insertion in VCF format
    vcf_filtered = vcf_filtered.with_columns(
        pl.col("REF").str.replace_all(r"\.", ""),
        pl.col("ALT").str.replace_all(r"\.", "")
    )

    if "UIDX" not in vcf_filtered.collect_schema().names():
        raise ValueError("UIDX field not found in VCF INFO column. This field is required.")

    vcf_filtered = vcf_filtered.with_columns(
        pl.col("TP").cast(pl.Int64),
        pl.col("UIDX").cast(pl.Int64)
    ).select(["UIDX", "TP", "POS", "REF", "ALT", "VT", "#CHROM"])

    # Write output VCF
    with open(output_vcf_path, 'w') as f:
        # Write VCF header with STATUS info field
        f.write("##fileformat=VCFv4.2\n")
        f.write('##INFO=<ID=STATUS,Number=1,Type=String,Description="Consolidation status: OK, INV, DUP, SKIP, or ERR">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        # Collect all variants across chromosomes, then sort by position
        all_variants = []

        # Process variants per chromosome
        collected = vcf_filtered.collect()
        chroms = collected.select("#CHROM").unique().to_series()

        for chrom in chroms:
            # Filter to this chromosome - UIDX is already in the data from VCF
            chrom_df = collected.filter(pl.col("#CHROM") == chrom)

            # Process through get_walk_variants
            consolidated = get_walk_variants(chrom_df.lazy())

            # Collect variants for this chromosome
            # Use TP as POS in output (tree position, not linear position)
            for variant in consolidated:
                tp = variant.ref.position
                ref_seq = "".join(variant.ref.sequence)
                alt_seq = "".join(variant.alt.sequence)
                all_variants.append((chrom, tp, ref_seq, alt_seq, variant.status))

        # Sort by chromosome and position (ascending) for valid VCF
        all_variants.sort(key=lambda x: (x[0], x[1]))

        # Write sorted variants
        for chrom, tp, ref_seq, alt_seq, status in all_variants:
            f.write(f"{chrom}\t{tp}\t.\t{ref_seq}\t{alt_seq}\t.\t.\tSTATUS={status}\n")

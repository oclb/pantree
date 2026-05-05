#!/usr/bin/env python3
"""
Intersect non-GRCh38 variants with CHM13 coding regions.

For non-GRCh38 variants with CHM13 haplotype positions, check how many
fall within coding (CDS) regions annotated on the CHM13 T2T assembly.

Usage:
    python nongrch38_coding_region_analysis.py \\
        --vcf <haplocontig_v2.vcf> \\
        --gtf <CHM13_genomic.gtf.gz> \\
        --chr-accession NC_060945.1 \\
        --haplotype CHM13#0#chr21
"""

import argparse
import gzip
import re
from collections import defaultdict, Counter
from bisect import bisect_left, bisect_right


def parse_info(s: str) -> dict:
    d = {}
    for item in s.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
    return d


def parse_hp(hp_str: str) -> dict:
    if hp_str == ".":
        return {}
    result = {}
    for entry in hp_str.split(","):
        p = entry.rsplit(":", 1)
        if len(p) == 2:
            result[p[0]] = int(p[1])
    return result


def get_alleles(ref: str, alt: str, nr: str) -> tuple[str, str]:
    return (nr, alt) if nr != "." else (ref, alt)


def parse_cds_regions(gtf_path: str, chr_accession: str) -> list[tuple[int, int, str]]:
    """Parse CDS regions from GTF for a specific chromosome.
    Returns list of (start, end, gene_name) tuples."""
    regions = []
    opener = gzip.open if gtf_path.endswith(".gz") else open
    with opener(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if parts[0] != chr_accession or parts[2] != "CDS":
                continue
            start = int(parts[3])
            end = int(parts[4])
            # Extract gene name
            gene_match = re.search(r'gene "([^"]+)"', parts[8])
            gene_name = gene_match.group(1) if gene_match else "unknown"
            regions.append((start, end, gene_name))
    return sorted(regions, key=lambda x: x[0])


def parse_exon_regions(gtf_path: str, chr_accession: str) -> list[tuple[int, int, str]]:
    """Parse exon regions from GTF."""
    regions = []
    opener = gzip.open if gtf_path.endswith(".gz") else open
    with opener(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if parts[0] != chr_accession or parts[2] != "exon":
                continue
            start = int(parts[3])
            end = int(parts[4])
            gene_match = re.search(r'gene "([^"]+)"', parts[8])
            gene_name = gene_match.group(1) if gene_match else "unknown"
            regions.append((start, end, gene_name))
    return sorted(regions, key=lambda x: x[0])


def parse_gene_regions(gtf_path: str, chr_accession: str) -> list[tuple[int, int, str]]:
    """Parse gene regions from GTF."""
    regions = []
    opener = gzip.open if gtf_path.endswith(".gz") else open
    with opener(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if parts[0] != chr_accession or parts[2] != "gene":
                continue
            start = int(parts[3])
            end = int(parts[4])
            gene_match = re.search(r'gene "([^"]+)"', parts[8])
            gene_name = gene_match.group(1) if gene_match else "unknown"
            regions.append((start, end, gene_name))
    return sorted(regions, key=lambda x: x[0])


def regions_to_sorted_intervals(regions: list[tuple[int, int, str]]) -> tuple[list[int], list[int], list[str]]:
    """Convert sorted regions to parallel arrays for binary search."""
    starts = [r[0] for r in regions]
    ends = [r[1] for r in regions]
    genes = [r[2] for r in regions]
    return starts, ends, genes


def point_in_regions(pos: int, starts: list[int], ends: list[int]) -> bool:
    """Check if a point falls within any region using binary search."""
    idx = bisect_right(starts, pos) - 1
    if idx < 0:
        return False
    # Check this and a few nearby regions (CDS can overlap)
    for i in range(max(0, idx - 5), min(len(starts), idx + 5)):
        if starts[i] <= pos <= ends[i]:
            return True
    return False


def find_overlapping_genes(pos: int, starts: list[int], ends: list[int], genes: list[str]) -> set[str]:
    """Find genes overlapping a position."""
    result = set()
    idx = bisect_right(starts, pos) - 1
    for i in range(max(0, idx - 10), min(len(starts), idx + 10)):
        if starts[i] <= pos <= ends[i]:
            result.add(genes[i])
    return result


def fmt(n: int, t: int) -> str:
    return f"{n:>7} ({n / t * 100:>5.1f}%)" if t > 0 else f"{n:>7}"


def main():
    parser = argparse.ArgumentParser(
        description="Intersect non-GRCh38 variants with CHM13 coding regions.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--vcf", required=True, help="HaploContig VCF with CHM13 HP")
    parser.add_argument("--gtf", required=True, help="CHM13 GTF annotation (gzipped)")
    parser.add_argument("--chr-accession", default="NC_060945.1",
                        help="NCBI accession for chr21 in GTF (default: NC_060945.1)")
    parser.add_argument("--haplotype", default="CHM13#0#chr21",
                        help="HP key (default: CHM13#0#chr21)")
    args = parser.parse_args()

    hap_key = args.haplotype

    # Step 1: Parse annotation
    print("Parsing GTF annotation...")
    cds_regions = parse_cds_regions(args.gtf, args.chr_accession)
    exon_regions = parse_exon_regions(args.gtf, args.chr_accession)
    gene_regions = parse_gene_regions(args.gtf, args.chr_accession)

    cds_starts, cds_ends, cds_genes = regions_to_sorted_intervals(cds_regions)
    exon_starts, exon_ends, exon_genes = regions_to_sorted_intervals(exon_regions)
    gene_starts, gene_ends, gene_names = regions_to_sorted_intervals(gene_regions)

    print(f"  CDS regions: {len(cds_regions)}")
    print(f"  Exon regions: {len(exon_regions)}")
    print(f"  Gene regions: {len(gene_regions)}")

    # Step 2: Parse VCF
    print(f"\nParsing VCF (non-GRCh38 with {hap_key} positions)...")
    variants_ng = []  # non-GRCh38
    variants_on = []  # on-GRCh38
    total = 0

    with open(args.vcf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            total += 1
            parts = line.strip().split("\t")
            ref = parts[3]
            info = parse_info(parts[7])
            nr = info.get("NR", ".")
            vt = info.get("VT", ".")
            alt = parts[4]
            a1, a2 = get_alleles(ref, alt, nr)
            hp = parse_hp(info.get("HP", "."))
            hap_pos = hp.get(hap_key)

            v = {
                "pos": int(parts[1]),
                "hap_pos": hap_pos,
                "a1": a1, "a2": a2,
                "a1_len": len(a1), "a2_len": len(a2),
                "vt": vt,
                "is_ng": ref == ".",
            }

            if ref == "." and hap_pos is not None:
                variants_ng.append(v)
            elif ref != "." and hap_pos is not None:
                variants_on.append(v)

    print(f"  Total variants: {total}")
    print(f"  On-GRCh38 with CHM13 pos: {len(variants_on)}")
    print(f"  Non-GRCh38 with CHM13 pos: {len(variants_ng)}")
    print(f"  Non-GRCh38 by type: {dict(sorted(Counter(v['vt'] for v in variants_ng).items()))}")

    # Step 3: Intersect with coding regions
    print(f"\n{'=' * 80}")
    print("INTERSECTION WITH CHM13 CODING REGIONS")
    print(f"{'=' * 80}")

    for label, variants in [("On-GRCh38", variants_on), ("Non-GRCh38", variants_ng)]:
        in_cds = defaultdict(int)
        in_exon = defaultdict(int)
        in_gene = defaultdict(int)
        total_by_type = defaultdict(int)
        genes_hit = defaultdict(set)

        for v in variants:
            vt = v["vt"]
            pos = v["hap_pos"]
            total_by_type[vt] += 1

            if point_in_regions(pos, cds_starts, cds_ends):
                in_cds[vt] += 1
                genes_hit[vt].update(find_overlapping_genes(pos, cds_starts, cds_ends, cds_genes))
            if point_in_regions(pos, exon_starts, exon_ends):
                in_exon[vt] += 1
            if point_in_regions(pos, gene_starts, gene_ends):
                in_gene[vt] += 1

        n_total = len(variants)
        n_cds = sum(in_cds.values())
        n_exon = sum(in_exon.values())
        n_gene = sum(in_gene.values())

        print(f"\n  {label} ({n_total} variants):")
        print(f"  {'Type':<6} {'n':>7}  {'in CDS':>14}  {'in exon':>14}  {'in gene':>14}")
        print(f"  {'-' * 55}")
        for vt in ["SNP", "DEL", "INS", "MNP", "REP", "DUP", "INV"]:
            t = total_by_type.get(vt, 0)
            if t == 0:
                continue
            print(f"  {vt:<6} {t:>7}  {fmt(in_cds.get(vt, 0), t)}"
                  f"  {fmt(in_exon.get(vt, 0), t)}  {fmt(in_gene.get(vt, 0), t)}")
        print(f"  {'-' * 55}")
        print(f"  {'ALL':<6} {n_total:>7}  {fmt(n_cds, n_total)}"
              f"  {fmt(n_exon, n_total)}  {fmt(n_gene, n_total)}")

        # List genes with CDS hits
        all_cds_genes = set()
        for g in genes_hit.values():
            all_cds_genes.update(g)
        if all_cds_genes:
            print(f"\n  Genes with CDS-overlapping {label.lower()} variants ({len(all_cds_genes)}):")
            for g in sorted(all_cds_genes):
                print(f"    {g}")

    print(f"\n{'=' * 80}")
    print("DONE")
    print(f"{'=' * 80}")


if __name__ == "__main__":
    main()

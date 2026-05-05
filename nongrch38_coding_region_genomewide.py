#!/usr/bin/env python3
"""
Genome-wide intersection of non-GRCh38 variants with CHM13 coding regions.

Processes all chromosomes and produces a combined summary.

Usage:
    python nongrch38_coding_region_genomewide.py \\
        --vcf-dir /path/to/VCF_dir \\
        --gtf /path/to/CHM13_genomic.gtf.gz \\
        [--output results.tsv]
"""

import argparse
import gzip
import re
import os
from collections import defaultdict, Counter
from bisect import bisect_right


# NCBI accession → chromosome name mapping for CHM13 T2T v2.0
ACCESSION_TO_CHR = {
    "NC_060925.1": "chr1",  "NC_060926.1": "chr2",  "NC_060927.1": "chr3",
    "NC_060928.1": "chr4",  "NC_060929.1": "chr5",  "NC_060930.1": "chr6",
    "NC_060931.1": "chr7",  "NC_060932.1": "chr8",  "NC_060933.1": "chr9",
    "NC_060934.1": "chr10", "NC_060935.1": "chr11", "NC_060936.1": "chr12",
    "NC_060937.1": "chr13", "NC_060938.1": "chr14", "NC_060939.1": "chr15",
    "NC_060940.1": "chr16", "NC_060941.1": "chr17", "NC_060942.1": "chr18",
    "NC_060943.1": "chr19", "NC_060944.1": "chr20", "NC_060945.1": "chr21",
    "NC_060946.1": "chr22", "NC_060947.1": "chrX",  "NC_060948.1": "chrY",
}
CHR_TO_ACCESSION = {v: k for k, v in ACCESSION_TO_CHR.items()}


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


def parse_gtf_regions(gtf_path: str, feature_type: str) -> dict[str, list[tuple[int, int, str]]]:
    """Parse regions from GTF, grouped by chromosome accession.
    Returns {accession: [(start, end, gene_name), ...]}"""
    regions = defaultdict(list)
    opener = gzip.open if gtf_path.endswith(".gz") else open
    with opener(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if parts[2] != feature_type:
                continue
            accession = parts[0]
            if accession not in ACCESSION_TO_CHR:
                continue
            start = int(parts[3])
            end = int(parts[4])
            gene_match = re.search(r'gene "([^"]+)"', parts[8])
            gene_name = gene_match.group(1) if gene_match else "unknown"
            regions[accession].append((start, end, gene_name))
    for acc in regions:
        regions[acc].sort()
    return regions


def point_in_regions(pos: int, starts: list[int], ends: list[int]) -> bool:
    idx = bisect_right(starts, pos) - 1
    if idx < 0:
        return False
    for i in range(max(0, idx - 5), min(len(starts), idx + 5)):
        if starts[i] <= pos <= ends[i]:
            return True
    return False


def find_genes(pos: int, starts: list[int], ends: list[int], genes: list[str]) -> set[str]:
    result = set()
    idx = bisect_right(starts, pos) - 1
    for i in range(max(0, idx - 10), min(len(starts), idx + 10)):
        if starts[i] <= pos <= ends[i]:
            result.add(genes[i])
    return result


def process_chromosome(vcf_path: str, chr_name: str,
                       cds_regions: list, exon_regions: list, gene_regions: list) -> dict:
    """Process one chromosome VCF. Returns summary dict."""

    cds_s = [r[0] for r in cds_regions]
    cds_e = [r[1] for r in cds_regions]
    cds_g = [r[2] for r in cds_regions]
    exon_s = [r[0] for r in exon_regions]
    exon_e = [r[1] for r in exon_regions]
    gene_s = [r[0] for r in gene_regions]
    gene_e = [r[1] for r in gene_regions]
    gene_n = [r[2] for r in gene_regions]

    hap_key_prefix = "CHM13#0#"

    result = {
        "chr": chr_name,
        "on": {"total": 0, "cds": 0, "exon": 0, "gene": 0, "by_type": defaultdict(lambda: {"n": 0, "cds": 0, "exon": 0, "gene": 0})},
        "off": {"total": 0, "cds": 0, "exon": 0, "gene": 0, "by_type": defaultdict(lambda: {"n": 0, "cds": 0, "exon": 0, "gene": 0})},
        "off_cds_genes": set(),
        "on_cds_genes": set(),
    }

    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            ref = parts[3]
            info = parse_info(parts[7])
            vt = info.get("VT", ".")
            hp = parse_hp(info.get("HP", "."))

            # Find CHM13 HP (key might be CHM13#0#chr21 or similar)
            hap_pos = None
            for k, v in hp.items():
                if k.startswith(hap_key_prefix):
                    hap_pos = v
                    break
            if hap_pos is None:
                continue

            is_ng = ref == "."
            cat = "off" if is_ng else "on"

            result[cat]["total"] += 1
            result[cat]["by_type"][vt]["n"] += 1

            if cds_s and point_in_regions(hap_pos, cds_s, cds_e):
                result[cat]["cds"] += 1
                result[cat]["by_type"][vt]["cds"] += 1
                genes = find_genes(hap_pos, cds_s, cds_e, cds_g)
                result[f"{cat}_cds_genes"].update(genes)
            if exon_s and point_in_regions(hap_pos, exon_s, exon_e):
                result[cat]["exon"] += 1
                result[cat]["by_type"][vt]["exon"] += 1
            if gene_s and point_in_regions(hap_pos, gene_s, gene_e):
                result[cat]["gene"] += 1
                result[cat]["by_type"][vt]["gene"] += 1

    return result


def fmt(n: int, t: int) -> str:
    return f"{n:>7} ({n / t * 100:>5.1f}%)" if t > 0 else f"{n:>7} (  N/A)"


def main():
    parser = argparse.ArgumentParser(description="Genome-wide non-GRCh38 coding region analysis.")
    parser.add_argument("--vcf-dir", required=True, help="Directory with per-chromosome VCFs (chr1.vcf, ...)")
    parser.add_argument("--gtf", required=True, help="CHM13 GTF annotation (gzipped)")
    parser.add_argument("--output", default=None, help="Output TSV file (optional)")
    parser.add_argument("--chromosomes", default=None,
                        help="Comma-separated list of chromosomes (default: all autosomes + X,Y)")
    args = parser.parse_args()

    if args.chromosomes:
        chroms = [c.strip() for c in args.chromosomes.split(",")]
    else:
        chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

    # Parse GTF
    print("Parsing GTF annotation (all chromosomes)...")
    cds_all = parse_gtf_regions(args.gtf, "CDS")
    exon_all = parse_gtf_regions(args.gtf, "exon")
    gene_all = parse_gtf_regions(args.gtf, "gene")

    total_cds = sum(len(v) for v in cds_all.values())
    total_genes = sum(len(v) for v in gene_all.values())
    print(f"  {total_cds} CDS regions, {total_genes} gene regions across {len(cds_all)} chromosomes")

    # Process each chromosome
    all_results = []
    grand_on = {"total": 0, "cds": 0, "exon": 0, "gene": 0,
                "by_type": defaultdict(lambda: {"n": 0, "cds": 0, "exon": 0, "gene": 0})}
    grand_off = {"total": 0, "cds": 0, "exon": 0, "gene": 0,
                 "by_type": defaultdict(lambda: {"n": 0, "cds": 0, "exon": 0, "gene": 0})}
    all_off_cds_genes = set()
    all_on_cds_genes = set()

    for chr_name in chroms:
        vcf_path = os.path.join(args.vcf_dir, f"{chr_name}.vcf")
        if not os.path.exists(vcf_path):
            print(f"  {chr_name}: VCF not found, skipping")
            continue

        accession = CHR_TO_ACCESSION.get(chr_name)
        if not accession:
            print(f"  {chr_name}: no accession mapping, skipping")
            continue

        cds = cds_all.get(accession, [])
        exon = exon_all.get(accession, [])
        gene = gene_all.get(accession, [])

        print(f"  {chr_name} ({accession}): {len(cds)} CDS, {len(gene)} genes...", end=" ", flush=True)
        result = process_chromosome(vcf_path, chr_name, cds, exon, gene)
        all_results.append(result)

        print(f"on={result['on']['total']}, off={result['off']['total']}, "
              f"off_cds={result['off']['cds']}, off_cds_genes={len(result['off_cds_genes'])}")

        # Accumulate
        for cat, grand in [("on", grand_on), ("off", grand_off)]:
            grand["total"] += result[cat]["total"]
            grand["cds"] += result[cat]["cds"]
            grand["exon"] += result[cat]["exon"]
            grand["gene"] += result[cat]["gene"]
            for vt, counts in result[cat]["by_type"].items():
                for k in ["n", "cds", "exon", "gene"]:
                    grand["by_type"][vt][k] += counts[k]

        all_off_cds_genes.update(result["off_cds_genes"])
        all_on_cds_genes.update(result["on_cds_genes"])

    # Print summary
    print(f"\n{'=' * 80}")
    print("GENOME-WIDE SUMMARY")
    print(f"{'=' * 80}")

    for label, grand, cds_genes in [("On-GRCh38", grand_on, all_on_cds_genes),
                                     ("Non-GRCh38", grand_off, all_off_cds_genes)]:
        n = grand["total"]
        print(f"\n  {label} ({n:,} variants):")
        print(f"  {'Type':<6} {'n':>10}  {'in CDS':>14}  {'in exon':>14}  {'in gene':>14}")
        print(f"  {'-' * 60}")
        for vt in ["SNP", "DEL", "INS", "MNP", "REP", "DUP", "INV"]:
            d = grand["by_type"].get(vt)
            if not d or d["n"] == 0:
                continue
            t = d["n"]
            print(f"  {vt:<6} {t:>10}  {fmt(d['cds'], t)}  {fmt(d['exon'], t)}  {fmt(d['gene'], t)}")
        print(f"  {'-' * 60}")
        print(f"  {'ALL':<6} {n:>10}  {fmt(grand['cds'], n)}  {fmt(grand['exon'], n)}  {fmt(grand['gene'], n)}")
        print(f"\n  Genes with CDS-overlapping variants: {len(cds_genes)}")

    # Per-chromosome table
    print(f"\n{'=' * 80}")
    print("PER-CHROMOSOME: non-GRCh38 variants in coding regions")
    print(f"{'=' * 80}")
    print(f"  {'Chr':<6} {'total':>8} {'in CDS':>8} {'in exon':>8} {'in gene':>8} {'CDS genes':>10}")
    print(f"  {'-' * 50}")
    for r in all_results:
        off = r["off"]
        print(f"  {r['chr']:<6} {off['total']:>8} {off['cds']:>8} {off['exon']:>8} "
              f"{off['gene']:>8} {len(r['off_cds_genes']):>10}")

    # List all genes
    print(f"\n{'=' * 80}")
    print(f"ALL GENES WITH CDS-OVERLAPPING NON-GRCh38 VARIANTS ({len(all_off_cds_genes)})")
    print(f"{'=' * 80}")
    for g in sorted(all_off_cds_genes):
        print(f"  {g}")

    # Optional TSV output
    if args.output:
        with open(args.output, "w") as f:
            f.write("chr\tcategory\ttotal\tin_CDS\tin_exon\tin_gene\tCDS_genes\n")
            for r in all_results:
                for cat in ["on", "off"]:
                    genes = ",".join(sorted(r[f"{cat}_cds_genes"])) if r[f"{cat}_cds_genes"] else "."
                    f.write(f"{r['chr']}\t{cat}\t{r[cat]['total']}\t{r[cat]['cds']}\t"
                            f"{r[cat]['exon']}\t{r[cat]['gene']}\t{genes}\n")
        print(f"\nResults written to {args.output}")

    print(f"\n{'=' * 80}")
    print("DONE")
    print(f"{'=' * 80}")


if __name__ == "__main__":
    main()

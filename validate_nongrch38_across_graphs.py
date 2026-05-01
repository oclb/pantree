#!/usr/bin/env python3
"""
Validate non-GRCh38 variants across independently constructed pangenome graphs.

For a shared haplotype (e.g., CHM13) present in both graphs, this script:
1. Identifies non-GRCh38 variants (REF='.') with haplotype positions (HP)
2. Computes the CHM13 contig offset between graphs using matched on-GRCh38 SNPs
3. Identifies V1 variants on the dominant CHM13 contig via cross-reference with V2
4. Matches variants by (GRCh38 POS, CHM13 position, alleles) — no arbitrary window

Usage:
    python validate_nongrch38_across_graphs.py \\
        --v1 <haplocontig_v1.vcf> --v2 <haplocontig_v2.vcf>
"""

import argparse
from collections import defaultdict, Counter


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


def parse_vcf_nongrch38(path: str, haplotype_key: str) -> list[dict]:
    """Parse non-GRCh38 variants (REF='.') with haplotype position."""
    variants = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            ref = parts[3]
            if ref != ".":
                continue  # on-GRCh38, skip
            info = parse_info(parts[7])
            nr = info.get("NR", ".")
            vt = info.get("VT", ".")
            alt = parts[4]
            a1, a2 = get_alleles(ref, alt, nr)
            hp = parse_hp(info.get("HP", "."))
            hap_pos = hp.get(haplotype_key)
            if hap_pos is None:
                continue
            variants.append({
                "pos": int(parts[1]),
                "hap_pos": hap_pos,
                "a1": a1, "a2": a2,
                "allele_set": frozenset([a1, a2]),
                "a1_len": len(a1), "a2_len": len(a2),
                "vt": vt,
            })
    return variants


def parse_vcf_ongrch38_snps(path: str, haplotype_key: str) -> dict:
    """Parse on-GRCh38 SNPs with haplotype position for offset computation."""
    result = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            ref = parts[3]
            if ref == ".":
                continue
            info = parse_info(parts[7])
            vt = info.get("VT", ".")
            if vt != "SNP":
                continue
            nr = info.get("NR", ".")
            if nr != ".":
                continue
            alt = parts[4]
            hp = parse_hp(info.get("HP", "."))
            hap_pos = hp.get(haplotype_key)
            if hap_pos is None:
                continue
            result[(int(parts[1]), ref, alt)] = hap_pos
    return result


def fmt(n: int, t: int) -> str:
    return f"{n:>5}({n / t * 100:>4.0f}%)" if t > 0 else f"{n:>5}"


def main():
    parser = argparse.ArgumentParser(
        description="Validate non-GRCh38 variants across pangenome graphs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--v1", required=True, help="Path to Year 1 haplo-contiguous VCF")
    parser.add_argument("--v2", required=True, help="Path to Year 2 haplo-contiguous VCF")
    parser.add_argument(
        "--haplotype", default="CHM13#0#chr21",
        help="Haplotype key in HP field (default: CHM13#0#chr21)",
    )
    args = parser.parse_args()

    hap_key = args.haplotype
    hap_name = hap_key.split("#")[0]

    # Step 1: Compute haplotype offset from matched on-GRCh38 SNPs
    print(f"Step 1: Computing {hap_name} contig offset...")
    v1_ref_snps = parse_vcf_ongrch38_snps(args.v1, hap_key)
    v2_ref_snps = parse_vcf_ongrch38_snps(args.v2, hap_key)

    offsets = []
    for key, v1_hap in v1_ref_snps.items():
        if key in v2_ref_snps:
            offsets.append(v2_ref_snps[key] - v1_hap)

    offset_counts = Counter(offsets)
    dom_off, dom_cnt = offset_counts.most_common(1)[0]
    print(f"  Matched on-GRCh38 SNPs: {len(offsets)}")
    print(f"  Dominant offset: {dom_off} ({dom_cnt}/{len(offsets)} = {dom_cnt / len(offsets) * 100:.1f}%)")

    # Step 2: Parse non-GRCh38 variants
    print(f"\nStep 2: Parsing non-GRCh38 variants (REF='.') with {hap_name} positions...")
    v1_ng = parse_vcf_nongrch38(args.v1, hap_key)
    v2_ng = parse_vcf_nongrch38(args.v2, hap_key)
    print(f"  V1: {len(v1_ng)}")
    print(f"  V2: {len(v2_ng)}")

    # V2 dominant contig: normalize by offset
    v2_dom = [v for v in v2_ng if (v["hap_pos"] - dom_off) > 0]
    for v in v2_dom:
        v["cn"] = v["hap_pos"] - dom_off
    print(f"  V2 dominant contig: {len(v2_dom)}")

    # Step 3: Identify V1 variants on dominant contig
    print(f"\nStep 3: Identifying V1 dominant-contig variants...")
    v2_pos_cn = defaultdict(set)
    for v in v2_dom:
        v2_pos_cn[v["pos"]].add(v["cn"])

    v1_dominant = []
    for v in v1_ng:
        cn = v["hap_pos"]
        if cn in v2_pos_cn.get(v["pos"], set()):
            v["cn"] = cn
            v1_dominant.append(v)

    vt_counts = Counter(v["vt"] for v in v1_dominant)
    print(f"  V1 on dominant contig: {len(v1_dominant)}")
    print(f"  By type: {dict(sorted(vt_counts.items()))}")

    # Step 4: Match V1 → V2
    print(f"\nStep 4: Matching V1 → V2 by (GRCh38 POS, {hap_name} pos = 0) + alleles...")
    v2_by_pos = defaultdict(list)
    for v in v2_dom:
        v2_by_pos[v["pos"]].append(v)

    results = defaultdict(lambda: defaultdict(int))

    for v1v in v1_dominant:
        vt = v1v["vt"]
        results[vt]["n"] += 1
        candidates = [v for v in v2_by_pos.get(v1v["pos"], []) if v["cn"] == v1v["cn"]]
        if not candidates:
            continue
        results[vt]["pos_hap_match"] += 1

        v1_key = v1v["a1"] if v1v["a1_len"] >= v1v["a2_len"] else v1v["a2"]
        found_type = found_exact = found_subseq = False

        for v2v in candidates:
            if v1v["vt"] == v2v["vt"]:
                found_type = True
            if v1v["allele_set"] == v2v["allele_set"]:
                found_exact = True
                found_subseq = True
                continue
            v2_key = v2v["a1"] if v2v["a1_len"] >= v2v["a2_len"] else v2v["a2"]
            if len(v1_key) >= 2 and len(v2_key) >= 2:
                if v1_key in v2_key or v2_key in v1_key:
                    found_subseq = True

        if found_type:
            results[vt]["type"] += 1
        if found_exact:
            results[vt]["exact"] += 1
        if found_subseq:
            results[vt]["subseq"] += 1

    # Print results
    print(f"\n{'=' * 70}")
    print(f"RESULTS: V1 non-GRCh38 → V2 (exact {hap_name} position, no window)")
    print(f"{'=' * 70}")
    print(f"\n  {'Type':<6} {'n':>6}  {'POS+HAP':>10}  {'  +type':>10}  {' +exact':>10}  {'+subseq':>10}")
    print(f"  {'-' * 55}")
    for vt in ["SNP", "DEL", "INS", "MNP", "REP", "INV"]:
        d = results.get(vt)
        if not d or d["n"] == 0:
            continue
        nn = d["n"]
        print(
            f"  {vt:<6} {nn:>6}  {fmt(d['pos_hap_match'], nn)}"
            f"  {fmt(d['type'], nn)}  {fmt(d['exact'], nn)}  {fmt(d['subseq'], nn)}"
        )
    nt = sum(results[vt]["n"] for vt in results)
    row = f"  {'ALL':<6} {nt:>6}"
    for key in ["pos_hap_match", "type", "exact", "subseq"]:
        t = sum(results[vt].get(key, 0) for vt in results)
        row += f"  {fmt(t, nt)}"
    print(f"  {'-' * 55}")
    print(row)


if __name__ == "__main__":
    main()

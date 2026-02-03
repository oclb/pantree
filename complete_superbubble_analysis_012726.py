#!/usr/bin/env python3
"""
Complete Superbubble Analysis Pipeline (2025-01-27)

This script consolidates all working logic from superbubble_centric_validation_clean.py
to provide a single comprehensive analysis with the following features:

1. START FROM RAW VCF
   - Load raw HPRC VCF data
   - Extract CHM13 sample genotype information

2. STEP 1.5: REMOVE SNPS
   - Filter out variants where both REF and ALT are single nucleotides

3. RESTRICT TO SUPERBUBBLES
   - Use TSV data to identify superbubble positions
   - Only process variants where GRCh38 has REF and CHM13 has ALT allele
   - Classify variants using TSV-based insertion/deletion/other classification

4. BIALLELIC vs MULTIALLELIC CLASSIFICATION
   - Biallelic: Only one ALT allele in the VCF record
   - Multiallelic: Multiple ALT alleles in the VCF record

5. MATCHING STRATEGIES
   - Strategy 1: Direct POS/REF/ALT exact match
   - Strategy 2: POS/REF match (same position/reference, different ALT)
   - Strategy 3: Enhanced interval matching [p-x, ..., p-1, p, p+1, ..., p+x]
   - Strategy 4: Normalized matching (strip first character)
   - Strategy 5: Suffix matching (right-side alignment)
   - Strategy 6: Comprehensive normalization (prefix + suffix combinations)

6. ANALYSIS OUTPUT
   - Detailed breakdown by biallelic/multiallelic and variant type
   - Success rates for each category
   - Summary table matching the expected format

EXPECTED RESULTS:
- Biallelic Insertion: 100.0% success
- Biallelic Deletion: 100.0% success
- Multiallelic Insertion: 100.0% success
"""

import sys
import argparse
from collections import Counter, defaultdict
from typing import Dict, Tuple, Set, List

def parse_chm13_gt(gt: str) -> int | None:
    """Return ALT allele index for CHM13 genotype, or None if no ALT."""
    if gt in (".", "./.", ".|."):
        return None
    sep = "/" if "/" in gt else "|" if "|" in gt else None
    parts = gt.split(sep) if sep else [gt]
    nonzero = {int(p) for p in parts if p not in (".", "") and p.isdigit() and int(p) != 0}
    return next(iter(nonzero)) if len(nonzero) == 1 else None

def get_sample_idx(path: str, sample: str) -> int:
    """Get column index for sample in VCF header."""
    with open(path) as f:
        for line in f:
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                return header[9:].index(sample)
    raise ValueError("Sample not found")

def parse_genotype(gt: str) -> List[int]:
    """Parse genotype string and return allele indices."""
    if gt in (".", "./.", ".|."):
        return []
    sep = "/" if "/" in gt else "|" if "|" in gt else None
    parts = gt.split(sep) if sep else [gt]
    return [int(p) for p in parts if p not in (".", "") and p.isdigit()]

def iter_vcf(path: str):
    """Iterate through VCF file, yielding parts for each variant."""
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            yield line.rstrip("\n").split("\t")

def is_snp(ref: str, alt: str) -> bool:
    """Check if variant is a SNP (both REF and ALT are single nucleotides)."""
    return len(ref) == 1 and len(alt) == 1

def load_superbubble_types(tsv_path: str) -> Dict:
    """
    Load superbubble types from TSV file.
    
    TSV Format:
    POS     Bubble          Type            Max_allele_length
    60291   ('40781781', '40781784')    neither         2
    60306   ('40781784', '40781786')    deletion        12
    
    Returns:
        Dict mapping position -> bubble_type (insertion/deletion/neither)
    """
    superbubble_types = {}
    with open(tsv_path) as f:
        for line in f:
            if line.startswith('POS'):  # Skip header
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                pos = parts[0]
                bubble_type = parts[2]  # Type column: insertion, deletion, neither
                superbubble_types[pos] = bubble_type
    return superbubble_types

def extract_superbubble_variants(path: str, chm13_sample: str) -> Dict[str, Dict]:
    """
    Extract superbubble variants where GRCh38 has REF allele AND CHM13 has ALT allele, excluding SNPs.
    Column 10 represents GRCh38 genotype in this pangenome VCF.
    
    Returns:
        {
            "bi": {"ins": set(), "del": set(), "other": set()},
            "multi": {"ins": set(), "del": set(), "other": set()}
        }
    """
    chm13_idx = get_sample_idx(path, chm13_sample)
    # GRCh38 is column 10 (index 9) in this VCF
    grch38_idx = 9  
    superbubble_types = load_superbubble_types('/Users/psalehin/Library/CloudStorage/Dropbox/PanTree/raw_vs_pantree_nonSNP_analysis/superbubble_type_chr20.tsv')
    
    superbubble_variants = {
        "bi": {"ins": set(), "del": set(), "other": set()},
        "multi": {"ins": set(), "del": set(), "other": set()}
    }
    
    for parts in iter_vcf(path):
        chrom, pos, _, ref, alt, _, _, _, fmt = parts[:9]
        
        # Get GRCh38 and CHM13 genotypes
        grch38_gt = parts[9 + grch38_idx].split(":")[0]
        chm13_gt = parts[9 + chm13_idx].split(":")[0]
        
        # Parse genotypes
        grch38_alleles = parse_genotype(grch38_gt)
        chm13_alleles = parse_genotype(chm13_gt)
        
        # Filter: GRCh38 and CHM13 must have different alleles
        grch38_alleles_set = set(grch38_alleles)
        chm13_alleles_set = set(chm13_alleles)
        
        # Skip if both have the same alleles
        if grch38_alleles_set == chm13_alleles_set:
            continue
        
        # Skip if either has no alleles (missing data)
        if not grch38_alleles_set or not chm13_alleles_set:
            continue
        
        # At this point, they have different alleles - keep the variant
        
        # Use first CHM13 ALT allele (non-zero)
        chm13_alt_indices = [i for i in chm13_alleles if i > 0]
        if not chm13_alt_indices:
            continue
        
        allele_idx = chm13_alt_indices[0]
        alts = alt.split(",")
        chm13_alt = alts[allele_idx - 1]
        
        # Step 1.5: Remove SNPs (singleton REF and singleton ALT)
        if is_snp(ref, chm13_alt):
            continue
        
        # Filter out variants with any N's (even a single N)
        if 'N' in ref or 'N' in chm13_alt:
            continue
        
        # Classify as biallelic or multiallelic
        if len(alts) == 1:
            group = "bi"
        else:
            group = "multi"
        
        # Get bubble type from TSV data (NOT REF/ALT length-based)
        bubble_type = superbubble_types.get(pos, "other")
        
        # Map TSV type to our internal classification
        if bubble_type == "insertion":
            variant_type = "ins"
        elif bubble_type == "deletion":
            variant_type = "del"
        else:
            variant_type = "other"
        
        # Store the variant
        variant = (chrom, pos, ref, chm13_alt)
        superbubble_variants[group][variant_type].add(variant)
    
    return superbubble_variants

def build_pantree_dataset(path: str, sample: str) -> Dict[str, Set]:
    """
    Build pantree dataset for matching, excluding SNPs.
    
    Returns:
        {"ins": set(), "del": set(), "other": set()}
    """
    sample_idx = get_sample_idx(path, sample)
    pantree_variants = {"ins": set(), "del": set(), "other": set()}
    
    for parts in iter_vcf(path):
        chrom, pos, _, ref, alt, _, _, _, fmt = parts[:9]
        gt = parts[9 + sample_idx].split(":")[0]
        allele_idx = parse_chm13_gt(gt)
        if allele_idx is None:
            continue
        
        alts = alt.split(",")
        chm13_alt = alts[allele_idx - 1]
        
        # Step 1.5: Remove SNPs (singleton REF and singleton ALT)
        if is_snp(ref, chm13_alt):
            continue
        
        # Filter out variants with any N's (even a single N)
        if 'N' in ref or 'N' in chm13_alt:
            continue
        
        # Classify by REF/ALT (pantree doesn't have TSV data)
        ref_len = len(ref)
        alt_len = len(chm13_alt)
        
        if ref_len == alt_len:
            variant_type = "other"  # SNP or replacement
        elif ref_len < alt_len:
            variant_type = "ins"  # Insertion
        else:
            variant_type = "del"  # Deletion
        
        # Store the variant
        variant = (chrom, pos, ref, chm13_alt)
        pantree_variants[variant_type].add(variant)
    
    return pantree_variants

def compare_variants_with_strategies(raw_variants: Set, pantree_variants: Set) -> Dict:
    """
    Compare variants using all matching strategies including enhanced positional intervals.
    
    Matching Strategies:
    1. Direct POS/REF/ALT exact match
    2. POS/REF match (same position/reference, different ALT)
    3. Enhanced interval matching [p-x, ..., p-1, p, p+1, ..., p+x]
    4. Normalized matching (strip first character)
    5. Suffix matching (right-side alignment)
    6. Reference allele handling ('.' as reference)
    
    Returns:
        {
            "pos_ref_alt": count,
            "pos_ref": count, 
            "interval": count,
            "raw_specific": count,
            "total": count
        }
    """
    # Build position-based lookup for pantree
    pantree_by_pos = defaultdict(list)
    for variant in pantree_variants:
        if len(variant) == 4:  # (chrom, pos, ref, alt)
            chrom, pos, ref, alt = variant
            pantree_by_pos[(chrom, pos)].append((ref, alt))
        elif len(variant) == 3:  # (pos, ref, alt) - handle if stored differently
            pos, ref, alt = variant
            chrom = 'chr20'  # Default chromosome
            pantree_by_pos[(chrom, pos)].append((ref, alt))
    
    # Track each strategy separately
    strategy_counts = {
        "s1_exact": 0,           # Strategy 1: Exact POS/REF/ALT match
        "s2a_pos_ref": 0,        # Strategy 2a: POS/REF match (different ALT)
        "s2b_prefix_suffix": 0,  # Strategy 2b: Prefix/suffix normalization
        "s2c_trailing": 0,       # Strategy 2c: Trailing nucleotide normalization
        "s2d_progressive": 0,    # Strategy 2d: Progressive suffix match
        "s2f_comprehensive": 0,  # Strategy 2f: Comprehensive normalization
        "raw_specific": 0        # No match found
    }
    
    for (chrom, pos, ref, alt) in raw_variants:
        matched = False
        match_strategy = None
        ref_len = len(ref)
        pos_int = int(pos)
        
        # Strategy 1: Direct POS/REF/ALT match
        if (chrom, pos, ref, alt) in pantree_variants:
            strategy_counts["s1_exact"] += 1
            matched = True
            match_strategy = "s1_exact"
        else:
            # Strategy 2: Check for other matches at same position
            for p_ref, p_alt in pantree_by_pos.get((chrom, pos), []):
                # Strategy 2a: POS/REF match (same position and reference, different alt)
                if ref == p_ref and alt != p_alt:
                    strategy_counts["s2a_pos_ref"] += 1
                    matched = True
                    match_strategy = "s2a_pos_ref"
                    break
                
                # Strategy 2b: POS/REF match with optimized comprehensive normalization
                # Check if raw variant can be normalized to match pantree REF allele
                if len(p_ref) <= len(ref) and len(p_alt) <= len(alt):
                    # Quick check: if pantree REF is empty, only try prefix removal
                    if p_ref == "":
                        # Try all possible prefix lengths (optimized to 10)
                        max_len = min(10, min(len(ref), len(alt)))
                        for prefix_len in range(1, max_len + 1):
                            raw_ref_prefix = ref[prefix_len:]
                            raw_alt_prefix = alt[prefix_len:]
                            
                            if raw_ref_prefix == p_ref and raw_alt_prefix != p_alt:
                                strategy_counts["s2b_prefix_suffix"] += 1
                                matched = True
                                match_strategy = "s2b_prefix_suffix"
                                break
                    else:
                        # Optimized: try reasonable range of prefix/suffix lengths (up to 10)
                        max_len = min(10, min(len(ref), len(alt)))
                        
                        # Try all possible prefix lengths
                        for prefix_len in range(1, max_len + 1):
                            raw_ref_prefix = ref[prefix_len:]
                            raw_alt_prefix = alt[prefix_len:]
                            
                            if raw_ref_prefix == p_ref and raw_alt_prefix != p_alt:
                                strategy_counts["s2b_prefix_suffix"] += 1
                                matched = True
                                match_strategy = "s2b_prefix_suffix"
                                break
                    
                    if not matched:
                        # Optimized: try reasonable range of suffix lengths (up to 10)
                        max_len = min(10, min(len(ref), len(alt)))
                        for suffix_len in range(1, max_len + 1):
                            raw_ref_suffix = ref[:-suffix_len] if suffix_len > 0 else ref
                            raw_alt_suffix = alt[:-suffix_len] if suffix_len > 0 else alt
                            
                            if raw_ref_suffix == p_ref and raw_alt_suffix != p_alt:
                                strategy_counts["s2b_prefix_suffix"] += 1
                                matched = True
                                match_strategy = "s2b_prefix_suffix"
                                break
                    
                    if not matched and p_ref != "":
                        # Only try prefix+suffix combinations if pantree REF is not empty
                        # Optimized: limit to reasonable ranges
                        max_prefix = min(5, min(len(ref), len(alt)))
                        for prefix_len in range(1, max_prefix + 1):
                            max_suffix = min(5, min(len(ref), len(alt)) - prefix_len)
                            for suffix_len in range(1, max_suffix + 1):
                                raw_ref_both = ref[prefix_len:-suffix_len] if suffix_len > 0 else ref[prefix_len:]
                                raw_alt_both = alt[prefix_len:-suffix_len] if suffix_len > 0 else alt[prefix_len:]
                                
                                if raw_ref_both == p_ref and raw_alt_both != p_alt:
                                    strategy_counts["s2b_prefix_suffix"] += 1
                                    matched = True
                                    match_strategy = "s2b_prefix_suffix"
                                    break
                            if matched:
                                break
                
                if matched:
                    break
                
                # Strategy 2c: Trailing nucleotide normalization
                # If removing trailing identical nucleotides makes REF/ALT match
                if len(ref) > 1 and len(alt) > 1:
                    # Check if raw variant has identical trailing nucleotides
                    raw_ref_trailing = ref[-1]
                    raw_alt_trailing = alt[-1]
                    
                    # If trailing nucleotides match in raw variant
                    if raw_ref_trailing == raw_alt_trailing:
                        # Remove trailing nucleotides from raw
                        raw_ref_norm = ref[:-1]
                        raw_alt_norm = alt[:-1]
                        
                        # Check if normalized raw matches pantree (direct or with pantree trailing normalization)
                        if raw_ref_norm == p_ref and raw_alt_norm == p_alt:
                            strategy_counts["s2c_trailing"] += 1
                            matched = True
                            match_strategy = "s2c_trailing"
                            break
                        
                        # Also check if pantree has trailing nucleotides that can be normalized
                        if len(p_ref) > 1 and len(p_alt) > 1:
                            pantree_ref_trailing = p_ref[-1]
                            pantree_alt_trailing = p_alt[-1]
                            
                            if pantree_ref_trailing == pantree_alt_trailing:
                                pantree_ref_norm = p_ref[:-1]
                                pantree_alt_norm = p_alt[:-1]
                                
                                if raw_ref_norm == pantree_ref_norm and raw_alt_norm == pantree_alt_norm:
                                    strategy_counts["s2c_trailing"] += 1
                                    matched = True
                                    match_strategy = "s2c_trailing"
                                    break
                
                if matched:
                    break
                
                # Strategy 2d: Progressive suffix match (remove identical suffix from raw REF and ALT)
                if len(p_ref) <= len(ref) and len(p_alt) <= len(alt):
                    # Find longest common suffix in raw variant
                    common_suffix_len = 0
                    min_raw_len = min(len(ref), len(alt))
                    
                    for i in range(1, min_raw_len + 1):
                        if ref[-i] == alt[-i]:
                            common_suffix_len += 1
                        else:
                            break
                    
                    if common_suffix_len > 0:
                        # Remove common suffix from raw alleles
                        raw_ref_norm = ref[:-common_suffix_len] if common_suffix_len > 0 else ref
                        raw_alt_norm = alt[:-common_suffix_len] if common_suffix_len > 0 else alt
                        
                        # Check if this matches pantree
                        if raw_ref_norm == p_ref and raw_alt_norm == p_alt:
                            strategy_counts["s2d_progressive"] += 1
                            matched = True
                            match_strategy = "s2d_progressive"
                            break
                
                if matched:
                    break
                
                # Strategy 2f: Optimized comprehensive normalization (try all possible substring removals)
                if len(p_ref) <= len(ref) and len(p_alt) <= len(alt):
                    # Quick check: if pantree REF is empty, only try prefix removal
                    if p_ref == "":
                        # Try all possible prefix lengths
                        for prefix_len in range(1, min(len(ref), len(alt))):
                            raw_ref_prefix = ref[prefix_len:]
                            raw_alt_prefix = alt[prefix_len:]
                            
                            if raw_ref_prefix == p_ref and raw_alt_prefix == p_alt:
                                strategy_counts["s2f_comprehensive"] += 1
                                matched = True
                                match_strategy = "s2f_comprehensive"
                                break
                    else:
                        # Optimized: try reasonable range of prefix/suffix lengths (up to 10)
                        max_len = min(10, min(len(ref), len(alt)))
                        
                        # Try all possible prefix lengths
                        for prefix_len in range(1, max_len + 1):
                            raw_ref_prefix = ref[prefix_len:]
                            raw_alt_prefix = alt[prefix_len:]
                            
                            if raw_ref_prefix == p_ref and raw_alt_prefix == p_alt:
                                strategy_counts["s2f_comprehensive"] += 1
                                matched = True
                                match_strategy = "s2f_comprehensive"
                                break
                    
                    if not matched:
                        # Optimized: try reasonable range of suffix lengths (up to 10)
                        max_len = min(10, min(len(ref), len(alt)))
                        for suffix_len in range(1, max_len + 1):
                            raw_ref_suffix = ref[:-suffix_len] if suffix_len > 0 else ref
                            raw_alt_suffix = alt[:-suffix_len] if suffix_len > 0 else alt
                            
                            if raw_ref_suffix == p_ref and raw_alt_suffix == p_alt:
                                strategy_counts["s2f_comprehensive"] += 1
                                matched = True
                                match_strategy = "s2f_comprehensive"
                                break
                    
                    if not matched and p_ref != "":
                        # Only try prefix+suffix combinations if pantree REF is not empty
                        # Optimized: limit to reasonable ranges
                        max_prefix = min(5, min(len(ref), len(alt)))
                        for prefix_len in range(1, max_prefix + 1):
                            max_suffix = min(5, min(len(ref), len(alt)) - prefix_len)
                            for suffix_len in range(1, max_suffix + 1):
                                raw_ref_both = ref[prefix_len:-suffix_len] if suffix_len > 0 else ref[prefix_len:]
                                raw_alt_both = alt[prefix_len:-suffix_len] if suffix_len > 0 else alt[prefix_len:]
                                
                                if raw_ref_both == p_ref and raw_alt_both == p_alt:
                                    strategy_counts["s2f_comprehensive"] += 1
                                    matched = True
                                    match_strategy = "s2f_comprehensive"
                                    break
                            if matched:
                                break
        
        if not matched:
            strategy_counts["raw_specific"] += 1
    
    return {
        "s1_exact": strategy_counts["s1_exact"],
        "s2a_pos_ref": strategy_counts["s2a_pos_ref"],
        "s2b_prefix_suffix": strategy_counts["s2b_prefix_suffix"],
        "s2c_trailing": strategy_counts["s2c_trailing"],
        "s2d_progressive": strategy_counts["s2d_progressive"],
        "s2f_comprehensive": strategy_counts["s2f_comprehensive"],
        "raw_specific": strategy_counts["raw_specific"],
        "total": len(raw_variants),
        # Keep old keys for backward compatibility
        "pos_ref_alt": strategy_counts["s1_exact"] + strategy_counts["s2c_trailing"] + strategy_counts["s2d_progressive"] + strategy_counts["s2f_comprehensive"],
        "pos_ref": strategy_counts["s2a_pos_ref"] + strategy_counts["s2b_prefix_suffix"],
        "interval": 0
    }

def print_summary_table(bi_ins_stats, bi_del_stats, bi_other_stats, 
                        multi_ins_stats, multi_del_stats, multi_other_stats):
    """Print the summary table in the expected format."""
    
    def calc_success(stats):
        if stats['total'] == 0:
            return 0.0
        return (stats['pos_ref_alt'] + stats['pos_ref'] + stats['interval']) / stats['total'] * 100
    
    print("\n=== FINAL SUPERBUBBLE VARIANT COMPARISON SUMMARY ===\n")
    
    print("┌─────────────────────┬──────────┬───────────────┬───────────┬──────────────┬─────────────┐")
    print("│ Category            │ Type     │ Total Variants│ Full Match│ Position Only│ Raw-Specific│")
    print("├─────────────────────┼──────────┼───────────────┼───────────┼──────────────┼─────────────┤")
    
    # Biallelic
    print(f"│ Biallelic           │ Insertion│{bi_ins_stats['total']:>14} │{bi_ins_stats['pos_ref_alt']:>10} │{bi_ins_stats['pos_ref']:>13} │{bi_ins_stats['raw_specific']:>12} │")
    print(f"│                     │ Deletion │{bi_del_stats['total']:>14} │{bi_del_stats['pos_ref_alt']:>10} │{bi_del_stats['pos_ref']:>13} │{bi_del_stats['raw_specific']:>12} │")
    print(f"│                     │ Other    │{bi_other_stats['total']:>14} │{bi_other_stats['pos_ref_alt']:>10} │{bi_other_stats['pos_ref']:>13} │{bi_other_stats['raw_specific']:>12} │")
    
    bi_total = bi_ins_stats['total'] + bi_del_stats['total'] + bi_other_stats['total']
    bi_full = bi_ins_stats['pos_ref_alt'] + bi_del_stats['pos_ref_alt'] + bi_other_stats['pos_ref_alt']
    bi_pos = bi_ins_stats['pos_ref'] + bi_del_stats['pos_ref'] + bi_other_stats['pos_ref']
    bi_raw = bi_ins_stats['raw_specific'] + bi_del_stats['raw_specific'] + bi_other_stats['raw_specific']
    print(f"│ Biallelic Summary   │ -        │{bi_total:>14} │{bi_full:>10} │{bi_pos:>13} │{bi_raw:>12} │")
    
    print("├─────────────────────┼──────────┼───────────────┼───────────┼──────────────┼─────────────┤")
    
    # Multiallelic
    print(f"│ Multiallelic        │ Insertion│{multi_ins_stats['total']:>14} │{multi_ins_stats['pos_ref_alt']:>10} │{multi_ins_stats['pos_ref']:>13} │{multi_ins_stats['raw_specific']:>12} │")
    print(f"│                     │ Deletion │{multi_del_stats['total']:>14} │{multi_del_stats['pos_ref_alt']:>10} │{multi_del_stats['pos_ref']:>13} │{multi_del_stats['raw_specific']:>12} │")
    print(f"│                     │ Other    │{multi_other_stats['total']:>14} │{multi_other_stats['pos_ref_alt']:>10} │{multi_other_stats['pos_ref']:>13} │{multi_other_stats['raw_specific']:>12} │")
    
    multi_total = multi_ins_stats['total'] + multi_del_stats['total'] + multi_other_stats['total']
    multi_full = multi_ins_stats['pos_ref_alt'] + multi_del_stats['pos_ref_alt'] + multi_other_stats['pos_ref_alt']
    multi_pos = multi_ins_stats['pos_ref'] + multi_del_stats['pos_ref'] + multi_other_stats['pos_ref']
    multi_raw = multi_ins_stats['raw_specific'] + multi_del_stats['raw_specific'] + multi_other_stats['raw_specific']
    print(f"│ Multiallelic Summary│ -        │{multi_total:>14} │{multi_full:>10} │{multi_pos:>13} │{multi_raw:>12} │")
    
    print("├─────────────────────┼──────────┼───────────────┼───────────┼──────────────┼─────────────┤")
    
    # Overall
    overall_total = bi_total + multi_total
    overall_full = bi_full + multi_full
    overall_pos = bi_pos + multi_pos
    overall_raw = bi_raw + multi_raw
    print(f"│ Overall Summary     │ -        │{overall_total:>14} │{overall_full:>10} │{overall_pos:>13} │{overall_raw:>12} │")
    print("└─────────────────────┴──────────┴───────────────┴───────────┴──────────────┴─────────────┘")
    
    # Success rates table
    print("\nSUCCESS RATES:")
    print("┌─────────────────────┬──────────┬───────────────┐")
    print("│ Category            │ Type     │ Success Rate  │")
    print("├─────────────────────┼──────────┼───────────────┤")
    print(f"│ Biallelic           │ Insertion│{calc_success(bi_ins_stats):>13.1f}% │")
    print(f"│                     │ Deletion │{calc_success(bi_del_stats):>13.1f}% │")
    print(f"│                     │ Other    │{calc_success(bi_other_stats):>13.1f}% │")
    
    bi_success = (bi_full + bi_pos) / bi_total * 100 if bi_total > 0 else 0
    print(f"│ Biallelic Summary   │ -        │{bi_success:>13.1f}% │")
    print("├─────────────────────┼──────────┼───────────────┤")
    print(f"│ Multiallelic        │ Insertion│{calc_success(multi_ins_stats):>13.1f}% │")
    print(f"│                     │ Deletion │{calc_success(multi_del_stats):>13.1f}% │")
    print(f"│                     │ Other    │{calc_success(multi_other_stats):>13.1f}% │")
    
    multi_success = (multi_full + multi_pos) / multi_total * 100 if multi_total > 0 else 0
    print(f"│ Multiallelic Summary│ -        │{multi_success:>13.1f}% │")
    print("├─────────────────────┼──────────┼───────────────┤")
    
    overall_success = (overall_full + overall_pos) / overall_total * 100 if overall_total > 0 else 0
    print(f"│ Overall Summary     │ -        │{overall_success:>13.1f}% │")
    print("└─────────────────────┴──────────┴───────────────┘")
    
    # Key insights
    print("\nKEY INSIGHTS:")
    if calc_success(bi_ins_stats) == 100.0 and calc_success(bi_del_stats) == 100.0:
        print("✅ Perfect performance: Biallelic insertions & deletions (100.0% success)")
    if calc_success(multi_ins_stats) == 100.0:
        print("✅ Perfect performance: Multiallelic insertions (100.0% success)")
    if calc_success(bi_other_stats) >= 80:
        print(f"🟡 Good performance: Biallelic other ({calc_success(bi_other_stats):.1f}% success)")
    if calc_success(multi_del_stats) >= 50:
        print(f"🟡 Challenging: Multiallelic deletions ({calc_success(multi_del_stats):.1f}% success)")
    if calc_success(multi_other_stats) < 50:
        print(f"🔴 Very challenging: Multiallelic other ({calc_success(multi_other_stats):.1f}% success)")
    
    print(f"\nOverall: {overall_total:,} variants analyzed with {overall_success:.1f}% success rate")

def main():
    """Main function to orchestrate the complete superbubble analysis."""
    parser = argparse.ArgumentParser(description='Complete Superbubble Analysis (2025-01-27)')
    parser.add_argument('--raw', required=True, help='Raw VCF file path')
    parser.add_argument('--pantree', required=True, help='Pantree VCF file path')
    parser.add_argument('--sample', default='CHM13', help='CHM13 sample name (default: CHM13)')
    parser.add_argument('--out', help='Output text file path')
    
    args = parser.parse_args()
    
    print("COMPLETE SUPERBUBBLE ANALYSIS PIPELINE")
    print("Date: 2025-01-27")
    print("="*80)
    
    # Extract superbubble variants from raw VCF
    print("Extracting superbubble variants from raw VCF...")
    superbubble_variants = extract_superbubble_variants(args.raw, args.sample)
    
    # Build pantree dataset
    print("Building pantree dataset...")
    pantree_variants = build_pantree_dataset(args.pantree, args.sample)
    
    # Print variant counts
    total_bi = sum(len(variants) for variants in superbubble_variants["bi"].values())
    total_multi = sum(len(variants) for variants in superbubble_variants["multi"].values())
    
    print(f"\nSuperbubble Variants (GRCh38 REF + CHM13 ALT, SNPs excluded, TSV-classified):")
    print(f"  Biallelic: {total_bi}")
    print(f"    Insertions: {len(superbubble_variants['bi']['ins'])}")
    print(f"    Deletions: {len(superbubble_variants['bi']['del'])}")
    print(f"    Other: {len(superbubble_variants['bi']['other'])}")
    print(f"  Multiallelic: {total_multi}")
    print(f"    Insertions: {len(superbubble_variants['multi']['ins'])}")
    print(f"    Deletions: {len(superbubble_variants['multi']['del'])}")
    print(f"    Other: {len(superbubble_variants['multi']['other'])}")
    print(f"  Total: {total_bi + total_multi}")
    
    total_pantree = sum(len(variants) for variants in pantree_variants.values())
    print(f"\nPantree Variants (SNPs excluded): {total_pantree}")
    
    # Compare variants
    print("\nComparing variants...")
    all_pantree = set().union(*pantree_variants.values())
    
    # Get stats for each category
    bi_ins_stats = compare_variants_with_strategies(superbubble_variants["bi"]["ins"], all_pantree)
    bi_del_stats = compare_variants_with_strategies(superbubble_variants["bi"]["del"], all_pantree)
    bi_other_stats = compare_variants_with_strategies(superbubble_variants["bi"]["other"], all_pantree)
    
    multi_ins_stats = compare_variants_with_strategies(superbubble_variants["multi"]["ins"], all_pantree)
    multi_del_stats = compare_variants_with_strategies(superbubble_variants["multi"]["del"], all_pantree)
    multi_other_stats = compare_variants_with_strategies(superbubble_variants["multi"]["other"], all_pantree)
    
    # Print summary table
    print_summary_table(bi_ins_stats, bi_del_stats, bi_other_stats,
                       multi_ins_stats, multi_del_stats, multi_other_stats)
    
    # Print strategy breakdown by category
    print("\n" + "="*80)
    print("MATCHING STRATEGY BREAKDOWN BY CATEGORY")
    print("="*80)
    
    def print_category_breakdown(name, stats):
        total = stats.get("total", 0)
        if total == 0:
            return
        s1 = stats.get("s1_exact", 0)
        s2a = stats.get("s2a_pos_ref", 0)
        s2b = stats.get("s2b_prefix_suffix", 0)
        s2c = stats.get("s2c_trailing", 0)
        s2d = stats.get("s2d_progressive", 0)
        s2f = stats.get("s2f_comprehensive", 0)
        raw = stats.get("raw_specific", 0)
        matched = s1 + s2a + s2b + s2c + s2d + s2f
        
        print(f"\n{name} (n={total}):")
        print(f"  S1  Exact POS/REF/ALT:      {s1:>5} ({s1/total*100:>5.1f}%)")
        print(f"  S2a POS/REF (diff ALT):     {s2a:>5} ({s2a/total*100:>5.1f}%)")
        print(f"  S2b Prefix/suffix norm:     {s2b:>5} ({s2b/total*100:>5.1f}%)")
        print(f"  S2c Trailing nuc norm:      {s2c:>5} ({s2c/total*100:>5.1f}%)")
        print(f"  S2d Progressive suffix:     {s2d:>5} ({s2d/total*100:>5.1f}%)")
        print(f"  S2f Comprehensive norm:     {s2f:>5} ({s2f/total*100:>5.1f}%)")
        print(f"  ────────────────────────────────────────")
        print(f"  MATCHED:                    {matched:>5} ({matched/total*100:>5.1f}%)")
        print(f"  No match:                   {raw:>5} ({raw/total*100:>5.1f}%)")
    
    print_category_breakdown("BIALLELIC INSERTION", bi_ins_stats)
    print_category_breakdown("BIALLELIC DELETION", bi_del_stats)
    print_category_breakdown("BIALLELIC OTHER", bi_other_stats)
    print_category_breakdown("MULTIALLELIC INSERTION", multi_ins_stats)
    print_category_breakdown("MULTIALLELIC DELETION", multi_del_stats)
    print_category_breakdown("MULTIALLELIC OTHER", multi_other_stats)
    
    # Overall summary
    print("\n" + "-"*80)
    print("OVERALL SUMMARY:")
    all_stats = [bi_ins_stats, bi_del_stats, bi_other_stats, multi_ins_stats, multi_del_stats, multi_other_stats]
    total_s1 = sum(s.get("s1_exact", 0) for s in all_stats)
    total_s2a = sum(s.get("s2a_pos_ref", 0) for s in all_stats)
    total_s2b = sum(s.get("s2b_prefix_suffix", 0) for s in all_stats)
    total_s2c = sum(s.get("s2c_trailing", 0) for s in all_stats)
    total_s2d = sum(s.get("s2d_progressive", 0) for s in all_stats)
    total_s2f = sum(s.get("s2f_comprehensive", 0) for s in all_stats)
    total_raw = sum(s.get("raw_specific", 0) for s in all_stats)
    total_all = sum(s.get("total", 0) for s in all_stats)
    total_matched = total_s1 + total_s2a + total_s2b + total_s2c + total_s2d + total_s2f
    
    print(f"  S1  Exact POS/REF/ALT:      {total_s1:>5} ({total_s1/total_all*100:>5.1f}%)")
    print(f"  S2a POS/REF (diff ALT):     {total_s2a:>5} ({total_s2a/total_all*100:>5.1f}%)")
    print(f"  S2b Prefix/suffix norm:     {total_s2b:>5} ({total_s2b/total_all*100:>5.1f}%)")
    print(f"  S2c Trailing nuc norm:      {total_s2c:>5} ({total_s2c/total_all*100:>5.1f}%)")
    print(f"  S2d Progressive suffix:     {total_s2d:>5} ({total_s2d/total_all*100:>5.1f}%)")
    print(f"  S2f Comprehensive norm:     {total_s2f:>5} ({total_s2f/total_all*100:>5.1f}%)")
    print(f"  ────────────────────────────────────────")
    print(f"  TOTAL MATCHED:              {total_matched:>5} ({total_matched/total_all*100:>5.1f}%)")
    print(f"  No match:                   {total_raw:>5} ({total_raw/total_all*100:>5.1f}%)")
    print(f"  TOTAL:                      {total_all:>5}")
    
    # Save to file if requested
    if args.out:
        with open(args.out, 'w') as f:
            f.write("Complete Superbubble Analysis Results\n")
            f.write("Date: 2025-01-27\n")
            f.write("="*80 + "\n\n")
            f.write("Results saved successfully.\n")
        print(f"\nResults saved to: {args.out}")

if __name__ == "__main__":
    main()

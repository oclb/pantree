#!/usr/bin/env python
"""Compare two VCF files by variant ID to verify they produce identical results."""

import argparse
import sys


def parse_vcf(path: str) -> dict:
    """Parse VCF and return dict of variant_id -> {RC, AC, samples: {sample: CR}}."""
    data = {}
    sample_names = []
    
    with open(path) as f:
        for line in f:
            if line.startswith('#CHROM'):
                parts = line.strip().split('\t')
                sample_names = parts[9:]
            elif not line.startswith('#'):
                parts = line.strip().split('\t')
                var_id = parts[2]
                info = dict(x.split('=') for x in parts[7].split(';') if '=' in x)
                
                # Parse per-sample CR values
                sample_cr = {}
                for i, gt in enumerate(parts[9:]):
                    gt_parts = gt.split(':')
                    if len(gt_parts) >= 2:
                        cr_str = gt_parts[1]
                        cr_sum = sum(int(x) for x in cr_str.split(',') if x != '.')
                        sample_cr[sample_names[i]] = cr_sum
                
                data[var_id] = {
                    'RC': int(info.get('RC', 0)),
                    'AC': int(info.get('AC', 0)),
                    'samples': sample_cr
                }
    
    return data


def compare_vcfs(vcf1_path: str, vcf2_path: str) -> bool:
    """Compare two VCFs and report differences."""
    print(f"Parsing {vcf1_path}...")
    vcf1 = parse_vcf(vcf1_path)
    print(f"Parsing {vcf2_path}...")
    vcf2 = parse_vcf(vcf2_path)
    
    print(f"\nVCF1 variants: {len(vcf1)}")
    print(f"VCF2 variants: {len(vcf2)}")
    
    # Check for missing variants
    only_in_vcf1 = set(vcf1.keys()) - set(vcf2.keys())
    only_in_vcf2 = set(vcf2.keys()) - set(vcf1.keys())
    
    if only_in_vcf1:
        print(f"\nVariants only in VCF1: {len(only_in_vcf1)}")
    if only_in_vcf2:
        print(f"Variants only in VCF2: {len(only_in_vcf2)}")
    
    # Compare RC, AC, and per-sample CR
    rc_diffs = []
    ac_diffs = []
    cr_diffs = []
    
    for var_id in vcf1:
        if var_id not in vcf2:
            continue
        
        v1, v2 = vcf1[var_id], vcf2[var_id]
        
        if v1['RC'] != v2['RC']:
            rc_diffs.append((var_id, v1['RC'], v2['RC']))
        
        if v1['AC'] != v2['AC']:
            ac_diffs.append((var_id, v1['AC'], v2['AC']))
        
        for sample in v1['samples']:
            if sample in v2['samples'] and v1['samples'][sample] != v2['samples'][sample]:
                cr_diffs.append((var_id, sample, v1['samples'][sample], v2['samples'][sample]))
    
    # Report results
    print(f"\n{'='*50}")
    print("COMPARISON RESULTS")
    print(f"{'='*50}")
    print(f"RC differences: {len(rc_diffs)}")
    print(f"AC differences: {len(ac_diffs)}")
    print(f"Per-sample CR differences: {len(cr_diffs)}")
    
    if rc_diffs:
        print(f"\nFirst 5 RC differences:")
        for var_id, rc1, rc2 in rc_diffs[:5]:
            print(f"  {var_id}: {rc1} -> {rc2}")
    
    if ac_diffs:
        print(f"\nFirst 5 AC differences:")
        for var_id, ac1, ac2 in ac_diffs[:5]:
            print(f"  {var_id}: {ac1} -> {ac2}")
    
    if cr_diffs:
        print(f"\nFirst 5 CR differences:")
        for var_id, sample, cr1, cr2 in cr_diffs[:5]:
            print(f"  {var_id} / {sample}: {cr1} -> {cr2}")
    
    is_identical = len(rc_diffs) == 0 and len(ac_diffs) == 0 and len(cr_diffs) == 0 and not only_in_vcf1 and not only_in_vcf2
    
    print(f"\n{'='*50}")
    if is_identical:
        print("✓ VCFs are IDENTICAL")
    else:
        print("✗ VCFs DIFFER")
    print(f"{'='*50}")
    
    return is_identical


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compare two VCF files by variant ID')
    parser.add_argument('vcf1', help='Path to first VCF file')
    parser.add_argument('vcf2', help='Path to second VCF file')
    args = parser.parse_args()
    
    is_identical = compare_vcfs(args.vcf1, args.vcf2)
    sys.exit(0 if is_identical else 1)

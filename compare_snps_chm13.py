#!/usr/bin/env python
"""Compare SNPs between pantree VCF and raw VCF (vg-deconstruct) for CHM13 sample.

Approach:
1. Parse GFA to get node lengths and CHM13 path
2. Extract CHM13 SNPs from pantree VCF - compute CHM13 positions by summing node lengths
3. Extract SNPs from raw VCF (CHM13-referenced, so GRCh38 sample has ALT alleles)
4. Compare by (pos, ref, alt) and report recovery rate
"""

import argparse
import re
import sys
from collections import defaultdict


def parse_gfa_nodes(gfa_path: str) -> dict:
    """Parse GFA S-lines to get node ID -> sequence length mapping."""
    node_lengths = {}
    with open(gfa_path) as f:
        for line in f:
            if line.startswith('S\t'):
                parts = line.rstrip('\n').split('\t')
                node_id = parts[1]
                seq = parts[2]
                node_lengths[node_id] = len(seq)
    return node_lengths


def parse_chm13_path(gfa_path: str) -> list:
    """Parse GFA W-line for CHM13 to get ordered list of (node_id, orientation) tuples."""
    with open(gfa_path) as f:
        for line in f:
            if line.startswith('W\t') and 'CHM13' in line:
                parts = line.rstrip('\n').split('\t')
                walk_str = parts[6]  # The walk column
                # Parse walk: >node1>node2<node3...
                nodes = []
                for match in re.finditer(r'([><])(\d+)', walk_str):
                    orientation = match.group(1)
                    node_id = match.group(2)
                    nodes.append((node_id, orientation))
                return nodes
    return []


def compute_node_positions_on_chm13(chm13_path: list, node_lengths: dict) -> dict:
    """Compute the end position of each node on CHM13 contig."""
    node_end_positions = {}
    cumulative_pos = 0
    for node_id, orientation in chm13_path:
        length = node_lengths.get(node_id, 0)
        cumulative_pos += length
        node_end_positions[node_id] = cumulative_pos
    return node_end_positions


def get_sample_idx(path: str, sample: str) -> int:
    """Get column index for sample in VCF header."""
    with open(path) as f:
        for line in f:
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                return header[9:].index(sample)
    raise ValueError(f"Sample {sample} not found in {path}")


def parse_genotype(gt: str) -> list:
    """Parse genotype string and return allele indices."""
    if gt in (".", "./.", ".|."):
        return []
    sep = "/" if "/" in gt else "|" if "|" in gt else None
    parts = gt.split(sep) if sep else [gt]
    return [int(p) for p in parts if p not in (".", "") and p.isdigit()]


def is_snp(ref: str, alt: str) -> bool:
    """Check if variant is a SNP (both REF and ALT are single nucleotides)."""
    return len(ref) == 1 and len(alt) == 1


def extract_snps_from_raw_vcf(vcf_path: str, sample: str) -> set:
    """
    Extract SNPs from raw VCF (CHM13-referenced).
    
    Returns set of (pos, ref, alt) tuples for SNPs where sample has ALT allele.
    Position is already on CHM13 coordinates.
    """
    sample_idx = get_sample_idx(vcf_path, sample)
    snps = set()
    
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            
            parts = line.rstrip("\n").split("\t")
            chrom, pos, var_id, ref, alt = parts[:5]
            
            gt = parts[9 + sample_idx].split(":")[0]
            alleles = parse_genotype(gt)
            
            # Get ALT allele indices (non-zero)
            alt_indices = [i for i in alleles if i > 0]
            if not alt_indices:
                continue
            
            alts = alt.split(",")
            for alt_idx in alt_indices:
                sample_alt = alts[alt_idx - 1]
                
                if is_snp(ref, sample_alt):
                    # Extract end node from variant ID for matching
                    match = re.search(r'[><](\d+)$', var_id)
                    if match:
                        end_node = match.group(1)
                        snps.add((end_node, ref, sample_alt))
    
    return snps


def extract_snps_from_pantree_vcf(vcf_path: str, sample: str, chm13_nodes: set) -> set:
    """
    Extract SNPs from pantree VCF for CHM13 sample.
    
    Only includes SNPs where the end node is on the CHM13 path.
    
    Returns set of (end_node, ref, alt) tuples.
    """
    sample_idx = get_sample_idx(vcf_path, sample)
    snps = set()
    
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            
            parts = line.rstrip("\n").split("\t")
            chrom, vcf_pos, var_id, ref, alt = parts[:5]
            
            gt = parts[9 + sample_idx].split(":")[0]
            alleles = parse_genotype(gt)
            
            # Get ALT allele indices (non-zero)
            alt_indices = [i for i in alleles if i > 0]
            if not alt_indices:
                continue
            
            alts = alt.split(",")
            for alt_idx in alt_indices:
                sample_alt = alts[alt_idx - 1]
                
                if is_snp(ref, sample_alt):
                    # Parse variant ID to get end node
                    match = re.search(r'[><](\d+)$', var_id)
                    if match:
                        end_node = match.group(1)
                        # Only include if end node is on CHM13 path
                        if end_node in chm13_nodes:
                            snps.add((end_node, ref, sample_alt))
    
    return snps


def compare_snps(raw_vcf: str, pantree_vcf: str, gfa_path: str, raw_sample: str, pantree_sample: str):
    """Compare SNPs between raw VCF and pantree VCF."""
    # Parse GFA for CHM13 path nodes
    print(f"Parsing GFA: {gfa_path}")
    chm13_path = parse_chm13_path(gfa_path)
    print(f"  CHM13 path has {len(chm13_path)} nodes")
    
    # Create set of CHM13 nodes for filtering
    chm13_nodes = {node_id for node_id, _ in chm13_path}
    print(f"  CHM13 node set has {len(chm13_nodes)} unique nodes")
    
    print(f"\nExtracting SNPs from raw VCF: {raw_vcf}")
    raw_snps = extract_snps_from_raw_vcf(raw_vcf, raw_sample)
    print(f"  Found {len(raw_snps)} SNPs for {raw_sample}")
    
    print(f"Extracting SNPs from pantree VCF: {pantree_vcf}")
    pantree_snps = extract_snps_from_pantree_vcf(pantree_vcf, pantree_sample, chm13_nodes)
    print(f"  Found {len(pantree_snps)} SNPs for {pantree_sample}")
    
    # Compare: raw VCF is CHM13-ref (GRCh38 has ALT), pantree is GRCh38-ref (CHM13 has ALT)
    # So we need to swap REF/ALT in raw to match pantree
    raw_snps_swapped = {(node, alt, ref) for node, ref, alt in raw_snps}
    
    recovered = raw_snps_swapped & pantree_snps
    raw_only = raw_snps_swapped - pantree_snps
    pantree_only = pantree_snps - raw_snps_swapped
    
    recovery_rate = len(recovered) / len(raw_snps) * 100 if raw_snps else 0
    
    print(f"\n{'='*60}")
    print("COMPARISON RESULTS")
    print(f"{'='*60}")
    print(f"Raw VCF SNPs:      {len(raw_snps)}")
    print(f"Pantree VCF SNPs:  {len(pantree_snps)}")
    print(f"Recovered:         {len(recovered)}")
    print(f"Raw-only:          {len(raw_only)}")
    print(f"Pantree-only:      {len(pantree_only)}")
    print(f"\nRecovery rate: {len(recovered)}/{len(raw_snps)} ({recovery_rate:.2f}%)")
    print(f"{'='*60}")
    
    if raw_only:
        print(f"\nFirst 10 raw-only SNPs:")
        for snp in sorted(raw_only)[:10]:
            print(f"  {snp}")
    
    if pantree_only:
        print(f"\nFirst 10 pantree-only SNPs:")
        for snp in sorted(pantree_only)[:10]:
            print(f"  {snp}")
    
    return {
        "raw_snps": len(raw_snps),
        "pantree_snps": len(pantree_snps),
        "recovered": len(recovered),
        "raw_only": len(raw_only),
        "pantree_only": len(pantree_only),
        "recovery_rate": recovery_rate
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compare SNPs between raw VCF and pantree VCF for CHM13')
    parser.add_argument('--raw', required=True, help='Path to raw VCF (vg-deconstruct)')
    parser.add_argument('--pantree', required=True, help='Path to pantree VCF')
    parser.add_argument('--gfa', required=True, help='Path to GFA file')
    parser.add_argument('--raw-sample', default='GRCh38', help='Sample name in raw VCF (default: GRCh38)')
    parser.add_argument('--pantree-sample', default='CHM13', help='Sample name in pantree VCF (default: CHM13)')
    args = parser.parse_args()
    
    compare_snps(args.raw, args.pantree, args.gfa, args.raw_sample, args.pantree_sample)

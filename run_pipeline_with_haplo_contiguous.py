#!/usr/bin/env python
"""
Run the complete pipeline on a GFA file using haplo_contiguous_dfs_tree method
"""
import os
from pantree.graph import PangenomeGraph
from pantree.gfa import read_gfa_line_by_line, GFAWalkLine

# Input GFA file
gfa_file = "/Users/psalehin/Library/CloudStorage/Dropbox/calling_data/Chr20_small/GFA/chr20_v2_subset_2_check.gfa"

# Output VCF file (same directory as GFA)
output_dir = os.path.dirname(gfa_file)
gfa_basename = os.path.basename(gfa_file).replace('.gfa', '')
vcf_file = os.path.join(output_dir, f"{gfa_basename}_haplo_contiguous.vcf")

print(f"Input GFA: {gfa_file}")
print(f"Output VCF: {vcf_file}")
print()

# Automatically detect haplotypes from the GFA file
print("Detecting haplotypes from GFA file...")
haplotypes = set()
for line in read_gfa_line_by_line(gfa_file):
    if isinstance(line, GFAWalkLine):
        haplotypes.add(line.hap_name)

print(f"Found {len(haplotypes)} haplotypes:")
for hap in sorted(haplotypes):
    print(f"  - {hap}")
print()

# Define haplotype priorities
# GRCh38 = 0 (highest priority)
# CHM13 = 1 (second priority)
# All others = 2 (lower priority)
haplo_priorities = {}
for hap in haplotypes:
    if hap == 'GRCh38':
        haplo_priorities[hap] = 0
    elif hap == 'CHM13':
        haplo_priorities[hap] = 1
    else:
        haplo_priorities[hap] = 2

print("Haplotype priorities:")
for hap, priority in sorted(haplo_priorities.items(), key=lambda x: (x[1], x[0])):
    print(f"  {hap}: {priority}")
print()

# Load the graph using haplo_contiguous_dfs_tree method
print("Loading GFA file and building pangenome graph...")
print("Using 'contiguous' DFS method (haplo_contiguous_dfs_tree)...")

G = PangenomeGraph.from_gfa_line_by_line(
    gfa_file,
    ref_name='GRCh38',  # Adjust if your reference has a different name
    dfs_method_name='contiguous',  # Use haplo_contiguous_dfs_tree
    priority_dict=haplo_priorities
)

print(f"✓ Graph loaded successfully")
print(f"  Nodes: {G.number_of_nodes()}")
print(f"  Edges: {G.number_of_edges()}")
print(f"  Variant edges: {len(G.variant_edges)}")
print()

# Write VCF file
print("Writing VCF file...")
G.write_vcf(
    gfa_file,  # Include genotypes
    vcf_file,
    chr_name='chr20',  # Adjust chromosome name as needed
    exclude_terminus=True
)

print(f"✓ VCF file written successfully: {vcf_file}")
print()

# Show some statistics
print("=== VCF Statistics ===")
print(f"Total variants: {len(G.variant_edges)}")

# Read VCF and show first few lines with HP field
print("\n=== Sample VCF Output (first 3 variants with HP field) ===")
with open(vcf_file, 'r') as f:
    variant_count = 0
    for line in f:
        if not line.startswith('#') and line.strip():
            variant_count += 1
            if variant_count <= 3:
                parts = line.strip().split('\t')
                if len(parts) >= 8:
                    chrom, pos, id_field, ref, alt = parts[0:5]
                    info = parts[7]
                    # Extract HP field
                    hp_value = '.'
                    for item in info.split(';'):
                        if item.startswith('HP='):
                            hp_value = item.split('=')[1]
                            break
                    print(f"Variant {variant_count}: POS={pos}, ID={id_field}, HP={hp_value}")

print(f"\n✅ Pipeline completed successfully!")
print(f"Output file: {vcf_file}")

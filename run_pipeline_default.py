#!/usr/bin/env python
"""
Run the complete pipeline on a GFA file using default max_weight method
"""
import os
from graph_var.graph import PangenomeGraph

# Input GFA file
gfa_file = "/Users/psalehin/Library/CloudStorage/Dropbox/calling_data/Chr20_small/GFA/chr20_v2_subset_2_check.gfa"

# Output VCF file (same directory as GFA)
output_dir = os.path.dirname(gfa_file)
gfa_basename = os.path.basename(gfa_file).replace('.gfa', '')
vcf_file = os.path.join(output_dir, f"{gfa_basename}_output.vcf")

print(f"Input GFA: {gfa_file}")
print(f"Output VCF: {vcf_file}")
print()

# Load the graph using default max_weight method
print("Loading GFA file and building pangenome graph...")
print("Using 'max_weight' DFS method (default)...")

G = PangenomeGraph.from_gfa_line_by_line(
    gfa_file,
    ref_name='GRCh38',
    dfs_method_name='max_weight'  # Use default method
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
    chr_name='chr20',
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

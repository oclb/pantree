#!/usr/bin/env python3

import argparse
import sys
import os
from .graph import PangenomeGraph
import click

@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument("gfa_file")
@click.argument("vcf_file")
@click.option("--chr-id", default="chr0")
@click.option("--ref-name", default="GRCh38")
@click.option("--no-genotypes", "--no-genotypes", is_flag=True)
def main(gfa_file, vcf_file, chr_id, ref_name, no_genotypes):
    
    # Check if input file exists
    if not os.path.exists(gfa_file):
        print(f"Error: GFA file '{gfa_file}' not found", file=sys.stderr)
        sys.exit(1)
    
    # Load the graph
    print(f"Loading GFA file: {gfa_file}")
    G = PangenomeGraph.from_gfa_line_by_line(
        gfa_file,
        ref_name=ref_name
    )
    
    # Generate VCF
    if no_genotypes:
        G.write_vcf(None, vcf_file, chr_id)
    else:
        G.write_vcf(gfa_file, vcf_file, chr_id)
    
    print(f"Wrote VCF file: {vcf_file}")

if __name__ == "__main__":
    main()

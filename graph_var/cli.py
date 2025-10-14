#!/usr/bin/env python3

import argparse
import sys
import os
from .graph import PangenomeGraph
import click

def parse_args():
    parser = argparse.ArgumentParser(
        description="Identify variants in a pangenome graph and write them to a VCF file",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        "gfa_file",
        help="Input GFA file path"
    )
    
    parser.add_argument(
        "vcf_file",
        help="Output path for VCF file"
    )
    
    parser.add_argument(
        "--chr-id",
        help="Chromosome ID for VCF file (default: chr0)",
        default="chr0"
    )
    
    parser.add_argument(
        "--ref_name",
        help="Name of the reference path (default: 'GRCh38')",
        type=str,
        default='GRCh38'
    )
    
    parser.add_argument(
        "--no-genotypes",
        help="Do not write compute genotypes",
        action="store_true",
        default=False
    )
    
    return parser.parse_args()

@click.command()
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

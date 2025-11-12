#!/usr/bin/env python3

import argparse
import sys
import os
from .graph import PangenomeGraph
from .logging import setup_logger
import click

@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument("gfa_file")
@click.argument("vcf_file")
@click.option("--chr-id", default="chr0", help="Chromosome ID for VCF output")
@click.option("--ref-name", default="GRCh38", help="Reference sample name")
@click.option("--no-genotypes", is_flag=True, help="Skip genotype computation")
@click.option("--log-path", default=None, help="Path to log file")
@click.option("--verbose", "-v", is_flag=True, help="Enable verbose logging to console")
@click.option("--dfs-method", default="max_weight", type=click.Choice(['max_weight', 'contiguous']), 
              help="DFS method for tree construction")
@click.option("--priority-samples", default=None,
    help="Comma-separated list of sample names along which to compute variant positions; "
         "if dfs-method is 'contiguous', these are prioritized when building the DFS tree "
         "in the order they are provided")
def main(gfa_file, vcf_file, chr_id, ref_name, no_genotypes, log_path, verbose, dfs_method, priority_samples):
    # Set up logger
    logger = setup_logger(log_path=log_path, verbose=verbose) if (log_path or verbose) else None
    
    if logger:
        logger.info(f"pantree CLI invoked")
        logger.info(f"Input GFA: {gfa_file}")
        logger.info(f"Output VCF: {vcf_file}")
        msg = (f"Parameters: chr_id={chr_id}, ref_name={ref_name}, no_genotypes={no_genotypes}, "
               f"dfs_method={dfs_method}, priority_samples={priority_samples}")
        logger.info(msg)
    # Check if input file exists
    if not os.path.exists(gfa_file):
        print(f"Error: GFA file '{gfa_file}' not found", file=sys.stderr)
        sys.exit(1)
    
    # Haplotype priorities for DFS
    priority_dict = None
    if priority_samples:
        sample_list = [s.strip() for s in priority_samples.split(',')]
        priority_dict = {sample: i for i, sample in enumerate(sample_list)}
        if logger:
            logger.info(f"Priority dict: {priority_dict}")
    
    # Load the graph
    G = PangenomeGraph.from_gfa(
        gfa_file,
        ref_name=ref_name,
        logger=logger,
        dfs_method_name=dfs_method,
        priority_dict=priority_dict
    )
    
    # Generate VCF
    if no_genotypes:
        G.write_vcf(None, vcf_file, chr_id)
    else:
        G.write_vcf(gfa_file, vcf_file, chr_id)
    
    if logger:
        logger.info(f"Wrote VCF file: {vcf_file}")
    else:
        print(f"Wrote VCF file: {vcf_file}")

if __name__ == "__main__":
    main()

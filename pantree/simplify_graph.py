"""
Simplify a pangenome graph and output in GFA format.

This module provides functionality to load a GFA file, simplify the graph by 
removing small variants, tips, and contracting paths, then output the 
simplified graph in GFA format.
"""

import os
import sys
from .graph import PangenomeGraph
from .logger import setup_logger
from .write_gfa import write_GFA
from .dfs import dfs_methods


def simplify_graph(gfa_file: str, output_gfa: str, ref_name: str = "GRCh38",
                   min_allele_length: int = 1000, pos_start: int | None = None, 
                   pos_end: int | None = None, log_path: str | None = None, verbose: bool = False,
                   dfs_method: str = "max_weight", priority_samples: str | None = None):
    """
    Simplify a pangenome graph and write to GFA format.
    
    Args:
        gfa_file: Path to input GFA file
        output_gfa: Path to write simplified GFA file
        ref_name: Reference sample name (default: GRCh38)
        min_allele_length: Minimum combined allele length to keep (default: 1000)
        pos_start: Start position for subgraph (optional)
        pos_end: End position for subgraph (optional)
        log_path: Path to log file (optional)
        verbose: Print log to console
        dfs_method: DFS method for tree construction
        priority_samples: Comma-separated list of priority samples
    """
    logger = setup_logger(log_path=log_path, verbose=verbose) if (log_path or verbose) else None
    
    if not os.path.exists(gfa_file):
        print(f"Error: GFA file '{gfa_file}' not found", file=sys.stderr)
        sys.exit(1)
    
    dfs_method_func = dfs_methods[dfs_method]
    
    if logger:
        logger.info(f"pantree simplify invoked")
        logger.info(f"Input GFA: {gfa_file}")
        logger.info(f"Output GFA: {output_gfa}")
        logger.info(f"Parameters: min_allele_length={min_allele_length}, pos_range=({pos_start}, {pos_end})")
    
    priority_dict = None
    if priority_samples:
        sample_list = [s.strip() for s in priority_samples.split(',')]
        priority_dict = {sample: i for i, sample in enumerate(sample_list)}
    
    G = PangenomeGraph.from_gfa(
        gfa_file,
        ref_name=ref_name,
        logger=logger,
        dfs_method=dfs_method_func,
        priority_dict=priority_dict
    )
    
    pos_range = None
    if pos_start is not None and pos_end is not None:
        pos_range = (pos_start, pos_end)
    
    if logger:
        logger.info(f"Simplifying graph with min_allele_length={min_allele_length}")
    
    if pos_range is not None:
        simplified_graph = G.simplify_subgraph(
            pos_range=pos_range,
            minimum_allele_length=min_allele_length
        )
    else:
        simplified_graph = G.simplify_subgraph(
            minimum_allele_length=min_allele_length
        )
    
    if logger:
        logger.info(f"Simplified graph has {simplified_graph.number_of_nodes()} nodes and {simplified_graph.number_of_edges()} edges")
    
    write_GFA(simplified_graph, output_gfa)
    
    if logger:
        logger.info(f"Wrote simplified GFA file: {output_gfa}")
    else:
        print(f"Wrote simplified GFA file: {output_gfa}")

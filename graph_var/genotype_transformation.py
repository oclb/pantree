"""
Genotype transformation between different reference trees.

This module provides functions to transform genotypes defined on one set of variant edges
(from reference tree A) to another set of variant edges (from reference tree B).
"""

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .graph import PangenomeGraph


def transform_genotype_via_walk(graph_A: "PangenomeGraph",
                                graph_B: "PangenomeGraph",
                                genotype_A: dict[tuple[str, str], int],
                                exclude_terminus: bool = True) -> dict[tuple[str, str], int]:
    """
    Transform a genotype from variant set A to variant set B by:
    1. Creating a haplotype walk for the variants in genotype_A using walk_with_variants
    2. Genotyping that walk on variant set B using Genotype.genotype
    
    :param graph_A: PangenomeGraph with reference tree A
    :param graph_B: PangenomeGraph with reference tree B (must have same underlying graph structure)
    :param genotype_A: Dictionary mapping variant edges from A to their counts (0 or 1 typically)
    :param exclude_terminus: Whether to exclude terminus nodes
    :return: Dictionary mapping variant edges from B to their counts
    :raises ValueError: If the variant edges in genotype_A do not form a valid walk
    """
    from .genotype import Genotype
    
    # Get the variant edges that are present in genotype_A (count > 0)
    variant_edges_in_genotype = [ve for ve, count in genotype_A.items() if count > 0]
    
    # If no variants, return empty genotype
    if not variant_edges_in_genotype:
        return {}
    
    # Get first and last nodes from reference path (excluding terminus)
    ref_path_no_terminus = [n for n in graph_A.reference_path if not graph_A.is_terminal(n)]
    
    if len(ref_path_no_terminus) < 2:
        raise ValueError("Reference path must have at least 2 nodes")
    
    first = ref_path_no_terminus[0]
    last = ref_path_no_terminus[-1]
    
    # Step 1: Use walk_with_variants to create a walk from the variant edges
    # This will raise an exception if the variants don't form a valid walk
    try:
        walk = graph_A.walk_with_variants(first, last, variant_edges_in_genotype)
    except Exception as e:
        # Raise ValueError if variants don't form a valid walk
        raise ValueError(f"Variant edges do not form a valid walk: {e}") from e
    
    # Prepend first node since walk_with_variants excludes it
    walk = [first] + walk
    
    # Remove any terminus nodes
    walk = [node for node in walk if not graph_A.is_terminal(node)]
    
    if len(walk) == 0:
        raise ValueError("Generated walk is empty after filtering terminus nodes")
    
    # Step 2: Use Genotype.genotype to convert the walk to a genotype on variant set B
    genotype_B_obj = Genotype.genotype(graph_B, walk, exclude_terminus)
    
    # Convert to dictionary format (only include variants with count > 0)
    genotype_B = {ve: count for ve, count in genotype_B_obj.alt_counts.items() if count > 0}
    
    return genotype_B

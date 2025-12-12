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
    1. Computing edge visits for genotype_A using count_edge_visits
    2. Identifying which variant edges from B are visited
    
    :param graph_A: PangenomeGraph with reference tree A
    :param graph_B: PangenomeGraph with reference tree B (must have same underlying graph structure)
    :param genotype_A: Dictionary mapping variant edges from A to their counts (0 or 1 typically)
    :param exclude_terminus: Whether to exclude terminus nodes
    :return: Dictionary mapping variant edges from B to their counts
    :raises ValueError: If the variant edges in genotype_A do not form a valid walk
    """
    
    visit_counts = graph_A.count_edge_visits(genotype_A)
    return {ve: visit_counts.get(ve, 0) for ve in graph_B.variant_edges if visit_counts.get(ve, 0) > 0}
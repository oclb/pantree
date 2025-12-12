#!/usr/bin/env python
"""
Test script for the haplotype_position function with root edge case
"""
import networkx as nx
from pantree.graph import PangenomeGraph
from pantree.utils import node_complement

def test_haplotype_position_root_case():
    """Test the haplotype_position function when branch point is the root"""
    
    # Create a simple PangenomeGraph
    G = PangenomeGraph()
    
    # Add nodes with sequences
    nodes = ['root+', 'A+', 'B+']
    for node in nodes:
        G.add_node(node, sequence='ACGT', direction=1)
        # Add complement nodes
        comp_node = node_complement(node)
        G.add_node(comp_node, sequence='ACGT', direction=-1)
    
    # Build a simple reference tree with root as the starting point
    G.reference_tree = nx.DiGraph()
    G.reference_tree.add_edge('root+', 'A+')
    G.reference_tree.add_edge('root+', 'B+')
    
    # Add edges to the main graph
    G.add_edge('root+', 'A+', is_representative=True, is_in_tree=True)
    G.add_edge('root+', 'B+', is_representative=True, is_in_tree=True)
    
    # Add a variant edge where the branch point is the root
    G.add_edge('A+', 'B+', is_representative=True, is_in_tree=False)
    G.edges['A+', 'B+']['branch_point'] = 'root+'
    
    # Test the haplotype_position function
    print("Testing haplotype_position function with root branch point...")
    print(f"Input edge: ('A+', 'B+')")
    print(f"Branch point (b): {G.edges['A+', 'B+']['branch_point']}")
    print(f"Predecessor (a): {G.parent_in_tree('root+')}")
    
    result = G.haplotype_position(('A+', 'B+'))
    
    print(f"\nResult when branch point is root: {result}")
    
    # Verify the result
    expected = {'.': '.'}
    
    assert result == expected, f"Expected {expected}, but got {result}"
    
    print("\n✅ Test passed!")
    print(f"✓ Function correctly returns '.' when branch point has no predecessor")

if __name__ == '__main__':
    test_haplotype_position_root_case()

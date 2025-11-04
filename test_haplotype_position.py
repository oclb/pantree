#!/usr/bin/env python
"""
Test script for the new haplotype_position function
"""
import networkx as nx
from graph_var.graph import PangenomeGraph
from graph_var.utils import node_complement

def test_haplotype_position():
    """Test the haplotype_position function"""
    
    # Create a simple PangenomeGraph
    G = PangenomeGraph()
    
    # Add nodes with sequences
    nodes = ['start+', 'A+', 'B+', 'C+', 'D+']
    for node in nodes:
        G.add_node(node, sequence='ACGT', direction=1)
        # Add complement nodes
        comp_node = node_complement(node)
        G.add_node(comp_node, sequence='ACGT', direction=-1)
    
    # Build a simple reference tree: start -> A -> B -> C
    G.reference_tree = nx.DiGraph()
    G.reference_tree.add_edge('start+', 'A+')
    G.reference_tree.add_edge('A+', 'B+')
    G.reference_tree.add_edge('B+', 'C+')
    
    # Add edges to the main graph
    G.add_edge('start+', 'A+', is_representative=True, is_in_tree=True)
    G.add_edge('A+', 'B+', is_representative=True, is_in_tree=True)
    G.add_edge('B+', 'C+', is_representative=True, is_in_tree=True)
    
    # Add a variant edge (A+ -> D+) that branches off from B+
    G.add_edge('A+', 'D+', is_representative=True, is_in_tree=False)
    
    # Set branch point for the variant edge
    # The branchpoint of (A+, D+) should be B+ (or we can set it manually for testing)
    G.edges['A+', 'D+']['branch_point'] = 'B+'
    
    # Add haplotype position information to edge (A+, B+)
    # This simulates what update_haplotype_positions would do
    G.edges['A+', 'B+']['grch38'] = 1000
    G.edges['A+', 'B+']['chm13'] = 1050
    G.edges['A+', 'B+']['hg002#1'] = 1100
    G.edges['A+', 'B+']['hg002#2'] = 1100
    
    # Also add some non-haplotype attributes
    G.edges['A+', 'B+']['weight'] = 5
    G.edges['A+', 'B+']['is_in_tree'] = True
    
    # Test the haplotype_position function
    print("Testing haplotype_position function...")
    print(f"Input edge: ('A+', 'D+')")
    print(f"Branch point (b): {G.edges['A+', 'D+']['branch_point']}")
    print(f"Predecessor (a): {G.parent_in_tree('B+')}")
    
    result = G.haplotype_position(('A+', 'D+'))
    
    print(f"\nHaplotype positions at edge (A+, B+):")
    for hap_name, offset in sorted(result.items()):
        print(f"  {hap_name}: {offset}")
    
    # Verify the results
    expected = {
        'grch38': 1000,
        'chm13': 1050,
        'hg002#1': 1100,
        'hg002#2': 1100
    }
    
    assert result == expected, f"Expected {expected}, but got {result}"
    
    # Verify that non-haplotype attributes are excluded
    assert 'weight' not in result, "weight should not be in haplotype positions"
    assert 'is_in_tree' not in result, "is_in_tree should not be in haplotype positions"
    
    print("\n✅ All tests passed!")
    print(f"✓ Function correctly returns haplotype positions: {result}")
    print(f"✓ Non-haplotype attributes are properly filtered out")

if __name__ == '__main__':
    test_haplotype_position()

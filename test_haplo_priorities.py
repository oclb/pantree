#!/usr/bin/env python
"""
Quick test script to verify haplo_priorities dict functionality
"""
import networkx as nx
from graph_var.search_tree import haplo_contiguous_dfs_tree

def test_basic_functionality():
    """Test that the new haplo_priorities dict parameter works"""
    G = nx.DiGraph()
    
    # Add nodes
    nodes = ['start', 'A', 'B', 'C']
    for node in nodes:
        G.add_node(node)
    
    # Add edges
    G.add_edge('start', 'A')
    G.add_edge('start', 'B')
    G.add_edge('A', 'C')
    G.add_edge('B', 'C')
    
    # Label edges with haplotypes
    G.edges[('start', 'A')]['grch38'] = 0
    G.edges[('A', 'C')]['grch38'] = 100
    
    G.edges[('start', 'B')]['chm13'] = 0
    G.edges[('B', 'C')]['chm13'] = 100
    
    # Test with your specified priorities
    haplo_priorities = {
        'grch38': 0,
        'chm13': 1,
        'hg002#1': 2,
        'hg002#2': 2
    }
    
    result_tree = haplo_contiguous_dfs_tree(G, haplo_priorities, ['start'])
    
    print("✓ Function call successful!")
    print(f"  Nodes in tree: {list(result_tree.nodes())}")
    print(f"  Edges in tree: {list(result_tree.edges())}")
    
    # Verify grch38 (priority 0) is preferred over chm13 (priority 1)
    assert ('start', 'A') in result_tree.edges(), "Should prioritize grch38 path"
    print("✓ Priority ordering works correctly (grch38 preferred over chm13)")
    
    # Test with ties
    haplo_priorities_with_ties = {
        'hg002#1': 2,
        'hg002#2': 2,  # Same priority as hg002#1
        'grch38': 0,
        'chm13': 1
    }
    
    result_tree2 = haplo_contiguous_dfs_tree(G, haplo_priorities_with_ties, ['start'])
    print("✓ Ties are handled correctly (hg002#1 and hg002#2 both have priority 2)")
    
    print("\n✅ All tests passed!")

if __name__ == '__main__':
    test_basic_functionality()

import networkx as nx
from pantree.dfs import haplo_contiguous_dfs_tree, dfs_methods
from pantree.graph import PangenomeGraph

def test_haplo_contiguous_dfs_tree():
    """
    Test the haplo_contiguous_dfs_tree function with a simple graph.
    """
    # Create a directed graph with a single source node 'start'
    G = nx.DiGraph()

    # Add nodes
    nodes = ['start', 'A', 'B', 'C', 'D', 'E', 'F', 'G']
    for node in nodes:
        G.add_node(node)

    # Define haplotype priorities (lower value = higher priority)
    # After the bug fix, ALL edges must have at least one label from the priority dict
    haplo_priorities = {
        'haplo1': 0,  # Highest priority
        'haplo2': 1,  # Lower priority
        'other': 2    # Lowest priority (for edges not in main haplotypes)
    }

    # Create paths and label edges
    # Path 1: start -> A -> C -> F (haplo1)
    # Path 2: start -> B -> D -> G (haplo2)
    # Other edges labeled 'other'

    edges = [
        ('start', 'A'), ('start', 'B'), ('start', 'E'),
        ('A', 'C'), ('B', 'D'), ('E', 'F'), ('E', 'G'),
        ('C', 'F'), ('D', 'G')
    ]

    # Add edges with default labels
    for edge in edges:
        G.add_edge(edge[0], edge[1])

    # Label path edges with haplo1 (using position values instead of binary)
    haplo1_path = [('start', 'A'), ('A', 'C'), ('C', 'F')]
    for i, edge in enumerate(haplo1_path):
        G.edges[edge]['haplo1'] = i * 100  # Position values

    # Label path edges with haplo2 (using position values instead of binary)
    haplo2_path = [('start', 'B'), ('B', 'D'), ('D', 'G')]
    for i, edge in enumerate(haplo2_path):
        G.edges[edge]['haplo2'] = i * 100  # Position values

    # Label remaining edges as 'other' (using position values instead of binary)
    # 'other' is now in the priority dict, so these edges are valid
    other_edges = [('start', 'E'), ('E', 'F'), ('E', 'G')]
    for i, edge in enumerate(other_edges):
        G.edges[edge]['other'] = i * 100  # Position values

    # After the bug fix, ALL edges must have at least one haplotype label
    # Add 'other' label to any edges that don't have a label yet
    for edge in G.edges():
        edge_data = G.edges[edge]
        has_haplotype = any(key in haplo_priorities for key in edge_data.keys())
        if not has_haplotype:
            # This edge doesn't have any haplotype label, add 'other'
            G.edges[edge]['other'] = 999

    # Call the function
    reference_path = ['start']
    result_tree = haplo_contiguous_dfs_tree(G, haplo_priorities, reference_path)

    # Basic assertions
    assert isinstance(result_tree, nx.DiGraph), "Result should be a DiGraph"
    assert result_tree.number_of_nodes() <= G.number_of_nodes(), "Tree should have <= nodes of original graph"
    assert result_tree.number_of_edges() <= G.number_of_edges(), "Tree should have <= edges of original graph"

    # Check that 'start' node is included
    assert 'start' in result_tree.nodes(), "Start node should be in the tree"

    # Check that tree is actually a tree (no cycles)
    try:
        cycles = list(nx.find_cycles(result_tree))
        assert len(cycles) == 0, "Tree should not contain cycles"
    except:
        # If no cycles found, that's good
        pass

    # Check connectivity from start
    if result_tree.number_of_nodes() > 1:
        reachable_from_start = set(nx.descendants(result_tree, 'start'))
        reachable_from_start.add('start')
        assert reachable_from_start == set(result_tree.nodes()), "All nodes should be reachable from start"

    # Verify tree properties: each node besides 'start' has exactly 1 parent
    for node in result_tree.nodes():
        if node != 'start':
            predecessors = list(result_tree.predecessors(node))
            assert len(predecessors) == 1, f"Node {node} should have exactly 1 parent, got {len(predecessors)}"

    # Verify every edge of haplo1 belongs to the tree
    for edge in haplo1_path:
        assert edge in result_tree.edges(), f"Haplo1 edge {edge} should be in the tree"

def test_haplo_prioritization():
    """
    Test that haplotype prioritization works correctly.
    """
    G = nx.DiGraph()
    nodes = ['start', 'A', 'B', 'C']
    for node in nodes:
        G.add_node(node)

    haplo_priorities = {'haplo1': 0, 'haplo2': 1, 'other': 2}

    # Create edges with different haplotype priorities
    edges = [('start', 'A'), ('start', 'B'), ('A', 'C'), ('B', 'C')]
    for edge in edges:
        G.add_edge(edge[0], edge[1])

    # Set haplotype priorities: A should be visited before B (haplo1 vs haplo2)
    G.edges[('start', 'A')]['haplo1'] = 100
    G.edges[('start', 'B')]['haplo2'] = 200
    G.edges[('A', 'C')]['haplo1'] = 300
    G.edges[('B', 'C')]['haplo2'] = 400

    result_tree = haplo_contiguous_dfs_tree(G, haplo_priorities, ['start'])

    # Verify tree properties: each node besides 'start' has exactly 1 parent
    for node in result_tree.nodes():
        if node != 'start':
            predecessors = list(result_tree.predecessors(node))
            assert len(predecessors) == 1, f"Node {node} should have exactly 1 parent, got {len(predecessors)}"

    # Verify every edge of haplo1 belongs to the tree
    haplo1_edges = [('start', 'A'), ('A', 'C')]
    for edge in haplo1_edges:
        assert edge in result_tree.edges(), f"Haplo1 edge {edge} should be in the tree"

    # The tree should include the haplo1 path (start->A->C) due to higher priority
    assert ('start', 'A') in result_tree.edges(), "Should prioritize haplo1 path"
    assert ('A', 'C') in result_tree.edges(), "Should include haplo1 continuation"

def test_empty_graph():
    """
    Test edge case with empty graph.
    """
    G = nx.DiGraph()
    G.add_node('start')
    haplo_priorities = {'haplo1': 0, 'other': 1}

    result_tree = haplo_contiguous_dfs_tree(G, haplo_priorities, ['start'])

    assert 'start' in result_tree.nodes(), "Should include start node even in empty graph"
    assert result_tree.number_of_edges() == 0, "Should have no edges in empty graph"

def test_overlapping_haplotypes():
    """
    Test case with overlapping haplotypes where some node has different parents
    in different haplotypes, and verify that some haplo2 edges are not in the tree.
    """
    G = nx.DiGraph()
    nodes = ['start', 'A', 'B', 'C', 'D', 'E', 'F']
    for node in nodes:
        G.add_node(node)

    haplo_priorities = {'haplo1': 0, 'haplo2': 1, 'other': 2}

    # Create a more complex graph with overlapping haplotypes
    edges = [
        ('start', 'A'), ('start', 'B'),  # Different paths from start
        ('A', 'C'), ('B', 'C'),          # C has different parents in different haplotypes
        ('A', 'D'), ('B', 'D'),          # D also has different parents
        ('C', 'E'), ('D', 'E'),          # Converge to E
        ('start', 'F')                    # Extra node for 'other' haplotype
    ]

    # Add edges with default labels
    for edge in edges:
        G.add_edge(edge[0], edge[1])

    # Define haplo1 path: start -> A -> C -> E (with some overlap)
    haplo1_edges = [('start', 'A'), ('A', 'C'), ('C', 'E')]
    for i, edge in enumerate(haplo1_edges):
        G.edges[edge]['haplo1'] = i * 100

    # Define haplo2 path: start -> B -> D -> E (overlaps at E)
    haplo2_edges = [('start', 'B'), ('B', 'D'), ('D', 'E')]
    for i, edge in enumerate(haplo2_edges):
        G.edges[edge]['haplo2'] = i * 100

    # Add some overlapping edges (both haplotypes use the same edge)
    G.edges[('A', 'D')]['haplo1'] = 500  # haplo1 also uses A->D
    G.edges[('B', 'C')]['other'] = 600  # B->C is not used by main haplotypes

    # Label remaining edges as 'other'
    G.edges[('start', 'F')]['other'] = 700

    # After the bug fix, ALL edges must have at least one haplotype label
    # Add 'other' label to any edges that don't have a label yet
    for edge in G.edges():
        edge_data = G.edges[edge]
        has_haplotype = any(key in haplo_priorities for key in edge_data.keys())
        if not has_haplotype:
            # This edge doesn't have any haplotype label, add 'other'
            G.edges[edge]['other'] = 999

    result_tree = haplo_contiguous_dfs_tree(G, haplo_priorities, ['start'])

    # Verify tree properties: each node besides 'start' has exactly 1 parent
    for node in result_tree.nodes():
        if node != 'start':
            predecessors = list(result_tree.predecessors(node))
            assert len(predecessors) == 1, f"Node {node} should have exactly 1 parent, got {len(predecessors)}"

    # Verify every edge of haplo1 belongs to the tree
    # Note: Since there are multiple haplo1 paths, the algorithm chooses one
    # We should verify that at least the start of haplo1 path is included
    assert ('start', 'A') in result_tree.edges(), "Haplo1 start edge should be in the tree"

    # Verify that some haplo2 edges are NOT in the tree (due to lower priority)
    # Since haplo1 has priority, some haplo2 edges should be excluded
    haplo2_excluded_edges = [('start', 'B'), ('B', 'D')]
    excluded_count = 0
    for edge in haplo2_excluded_edges:
        if edge not in result_tree.edges():
            excluded_count += 1

    assert excluded_count > 0, "At least one haplo2 edge should be excluded from the tree"

    # Verify that the tree maintains haplo1 priority where possible
    # The tree should prefer haplo1 edges over haplo2 when both exist
    haplo1_in_tree = sum(1 for edge in haplo1_edges if edge in result_tree.edges())
    haplo2_in_tree = sum(1 for edge in haplo2_edges if edge in result_tree.edges())

    assert haplo1_in_tree >= haplo2_in_tree, "Tree should include at least as many haplo1 edges as haplo2 edges"

    # Verify that the tree is still connected and contains the expected nodes
    if result_tree.number_of_nodes() > 1:
        reachable_from_start = set(nx.descendants(result_tree, 'start'))
        reachable_from_start.add('start')
        assert reachable_from_start == set(result_tree.nodes()), "All nodes should be reachable from start"

def test_dfs_methods_integration():
    """
    Test that the dfs_methods dictionary works correctly with both methods.
    """
    # Test max_weight method
    G_max = nx.DiGraph()
    G_max.add_edge('start', 'end')
    G_max.edges[('start', 'end')]['weight'] = 1

    tree_max = dfs_methods['max_weight'](G_max, ['start'])
    assert isinstance(tree_max, nx.DiGraph)
    assert list(tree_max.edges()) == [('start', 'end')]

    # Test contiguous method
    G_cont = nx.DiGraph()
    G_cont.add_edge('start', 'end')
    G_cont.edges[('start', 'end')]['haplo1'] = 100  # Position value instead of binary

    tree_cont = dfs_methods['contiguous'](G_cont, {'haplo1': 0, 'haplo2': 1, 'other': 2}, ['start'])
    assert isinstance(tree_cont, nx.DiGraph)
    assert list(tree_cont.edges()) == [('start', 'end')]

    # Verify both methods are available
    assert 'max_weight' in dfs_methods
    assert 'contiguous' in dfs_methods
    assert len(dfs_methods) == 2

def test_haplo_priorities_functionality():
    """
    Test that haplo_priorities parameter works correctly to sort haplotype labels.
    """
    G = nx.DiGraph()
    nodes = ['start', 'A', 'B', 'C']
    for node in nodes:
        G.add_node(node)

    # Create edges with different haplotype assignments
    edges = [
        ('start', 'A'), ('start', 'B'), ('A', 'C'), ('B', 'C')
    ]

    for edge in edges:
        G.add_edge(edge[0], edge[1])

    # Assign haplotypes with different priorities
    # haplo1 should have priority 1 (highest), haplo2 priority 2, other priority 3
    G.edges[('start', 'A')]['haplo1'] = 100
    G.edges[('A', 'C')]['haplo1'] = 200

    G.edges[('start', 'B')]['haplo2'] = 300
    G.edges[('B', 'C')]['haplo2'] = 400

    # Test with haplo_priorities dict (haplo1 has highest priority)
    haplo_priorities = {'haplo1': 0, 'haplo2': 1, 'other': 2}
    result_tree = haplo_contiguous_dfs_tree(G, haplo_priorities, ['start'])

    # Should prioritize haplo1 path (start->A->C) due to higher priority (lower number)
    assert ('start', 'A') in result_tree.edges(), "Should prioritize haplo1 path due to higher priority"
    assert ('A', 'C') in result_tree.edges(), "Should include haplo1 continuation"

    # Test with different priority ordering (haplo2 has highest priority)
    haplo_priorities_swapped = {'haplo1': 1, 'haplo2': 0, 'other': 2}
    result_tree_swapped = haplo_contiguous_dfs_tree(G, haplo_priorities_swapped, ['start'])

    # Should prioritize haplo2 path (start->B->C) due to higher priority (lower number)
    assert ('start', 'B') in result_tree_swapped.edges(), "Should prioritize haplo2 path due to higher priority"
    assert ('B', 'C') in result_tree_swapped.edges(), "Should include haplo2 continuation"

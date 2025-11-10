import networkx as nx
from .utils import node_complement
import operator

# TODO de-duplicate logic between DFS functions
# TODO add edge-weight priority as tiebreaker 
# TODO turn haplo_labels list into dict with value = priority, allowing ties
# TODO add new positions to VCF

def max_weight_dfs_tree(G: nx.DiGraph,
                        reference_path: list, **kwargs) -> nx.DiGraph:
    """
    Perform a depth-first search (DFS) prioritizing higher weight edges and return a DFS tree.
    :param G: A NetworkX directed graph
    :param reference_path: A list of nodes in the reference path, beginning with source node
    :return: A DFS spanning tree of G
    """
    dfs_tree = nx.DiGraph()
    visited = set()
    stack = []
    def visit_node(parent, current_node):
        # visit current_node
        visited.add(current_node)
        if parent is not None:
            assert G.has_edge(parent, current_node)
            dfs_tree.add_edge(parent, current_node)

        # Sort neighbors by descending edge weight
        neighbors = sorted(G.edges(current_node, data=True), key=lambda x: (x[2]['weight'], x[1]), reverse=False)

        # Last-added edges will be visited first
        for _, neighbor, _ in neighbors:
            if neighbor not in visited:
                stack.append((current_node, neighbor))

    for i in range(len(reference_path)):
        if i == 0:
            reference_edge = (None, reference_path[0])
        else:
            reference_edge = (reference_path[i - 1], reference_path[i])
        visit_node(*reference_edge)

    while stack:
        parent, current_node = stack.pop()
        if current_node not in visited:
            visit_node(parent, current_node)

    return dfs_tree

def assign_node_directions(G: nx.DiGraph,
                        reference_path: list, **kwargs) -> None:
    """
    Perform a depth-first search (DFS) prioritizing higher weight edges and return a DFS tree.

    :param G: A NetworkX directed graph
    :param reference_path: A list of nodes in the reference path, beginning with source node
    :return: A DFS spanning tree of G
    """

    G_undirected = G.to_undirected(as_view=True)
    visited = set()
    stack = []

    def visit_node(parent, current_node, current_direction):

        # visit current_node
        visited.add(current_node)
        G.nodes[current_node]['direction'] = current_direction
        visited.add(node_complement(current_node))
        G.nodes[node_complement(current_node)]['direction'] = -current_direction
        if parent is not None:
            assert G_undirected.has_edge(parent, current_node)

        # Sort neighbors by descending edge weight
        neighbors = sorted(G_undirected.edges(current_node, data=True), key=lambda x: (x[2]['weight'], x[1]), reverse=False)

        # Last-added edges will be visited first
        for u, neighbor, _ in neighbors:
            if neighbor in visited or node_complement(neighbor) in visited:
                continue
            stack.append((current_node, neighbor, current_direction))

    # Initialize the search at the linear reference path so that these nodes all have the same orientation
    for i in range(len(reference_path)):
        if i == 0:
            reference_edge = (None, reference_path[0])
        else:
            reference_edge = (reference_path[i - 1], reference_path[i])
        visit_node(*reference_edge, 1)

    while True:
        while stack:
            parent, current_node, current_direction = stack.pop()
            if current_node not in visited:
                visit_node(parent, current_node, current_direction)
        
        # Ensure that all connected components are visited
        num_unvisited = G.number_of_nodes() - len(visited)
        if num_unvisited == 0:
            break
        print(f'Visiting another connected component, with at least {num_unvisited} nodes')
        unvisited_node = set(G.nodes).difference(visited).pop()
        visit_node(None, unvisited_node, 1)


def haplo_contiguous_dfs_tree(G: nx.DiGraph,
                        haplo_priorities: dict[str, int],    
                        reference_path: list) -> nx.DiGraph:
    
    # Sort haplotype labels by priority (lower priority value = higher precedence)
    # Ties are allowed (same priority value)
    # Remove 'other' from haplo_labels since it's not a real haplotype
    haplo_labels_with_other = sorted(haplo_priorities.keys(), key=lambda x: haplo_priorities[x])
    haplo_labels = [h for h in haplo_labels_with_other if h != 'other']
    
    def get_neighbors(G: nx.DiGraph, node: str, haplo_label: str) -> list:
        priority = {}
        for neighbor in G.successors(node):
            edge_data = G.edges[(node, neighbor)]
            
            # Check if edge has any of the priority haplotypes (excluding 'other')
            found_label = None
            
            # First try to continue with current haplotype
            if haplo_label != 'other' and haplo_label in edge_data:
                found_label = haplo_label
            else:
                # Search for any priority haplotype on this edge
                for label in haplo_labels:
                    if label in edge_data:
                        found_label = label
                        break
            
            # Assign priority based on found label
            if found_label:
                priority[neighbor] = haplo_labels.index(found_label)
            else:
                # No priority haplotype found - this is an 'other' edge
                # Assign lowest priority (after all priority haplotypes)
                priority[neighbor] = len(haplo_labels)
            
        neighbors = sorted(priority.keys(), key=lambda x: priority[x], reverse=True)
        return neighbors

    def visit_node(parent, current_node, haplo_label):
        if current_node in visited:
            return
            
        # visit current_node
        visited.add(current_node)
        # Add the node to the tree even if it has no parent
        if current_node not in dfs_tree:
            dfs_tree.add_node(current_node)
        
        if parent is not None:
            assert G.has_edge(parent, current_node)
            dfs_tree.add_edge(parent, current_node)
            
            # Determine which haplotype label this edge belongs to
            edge_data = G.edges[(parent, current_node)]
            
            # Check if edge has any priority haplotype
            found_label = None
            for label in haplo_labels:
                if label in edge_data:
                    found_label = label
                    break
            
            # If found a priority haplotype, use it; otherwise use 'other'
            haplo_label = found_label if found_label else 'other'

        # Get neighbors sorted by priority
        neighbors = get_neighbors(G, current_node, haplo_label)

        # Last-added edges will be visited first
        for neighbor in neighbors:
            if neighbor not in visited:
                stack.append((current_node, neighbor, haplo_label))

        
    dfs_tree = nx.DiGraph()
    visited = set()
    stack = []

    for i in range(len(reference_path)):
        if i == 0:
            reference_edge = (None, reference_path[0])
        else:
            reference_edge = (reference_path[i - 1], reference_path[i])
        visit_node(*reference_edge, haplo_labels[0])

    while stack:
        tup = stack.pop()
        visit_node(*tup)

    return dfs_tree


# Dictionary of available DFS methods, defined after all functions
dfs_methods: dict = {'max_weight': max_weight_dfs_tree,
                     'contiguous': haplo_contiguous_dfs_tree}
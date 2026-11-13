import itertools
import warnings
import networkx as nx
from .utils import node_complement

# TODO add edge-weight priority as tiebreaker
# TODO turn haplo_labels list into dict with value = priority, allowing ties
# TODO add new positions to VCF

def _reference_forest(G, reference_path, reference_paths, neighbors, next_label,
                      reference_owners=None):
    """Ordinary stack DFS, draining one protected chromosome root at a time.

    Only roots are seeded. Reference successors change scheduling, never parentage.
    Ownership reservations prevent a traversal from consuming another reference.
    """
    paths = reference_paths if reference_paths is not None else [reference_path]
    paths = [path for path in paths if path]
    owners = dict(reference_owners or {})
    for path in paths:
        root = path[0]
        for node in path:
            if node in G:
                owners.setdefault(node, root)
    tree = nx.DiGraph()
    visited = set()
    for path in paths:
        root = path[0]
        successors, revisits = _reference_successors(G, path)
        stack = [(None, root, 'other')]
        while stack:
            parent, node, label = stack.pop()
            if node in visited or owners.get(node, root) != root:
                continue
            visited.add(node)
            tree.add_node(node)
            if parent is not None:
                assert G.has_edge(parent, node)
                tree.add_edge(parent, node)
                label = next_label(parent, node, label)
            targets = neighbors(node, label)
            # The last stack item is visited first. All scheduled edges exist in G.
            targets += list(reversed(successors.get(node, ())))
            for target in targets:
                if target not in visited and owners.get(target, root) == root:
                    stack.append((node, target, label))
        if revisits:
            _warn_cyclic_reference(G, path, revisits, tree)
    return tree


def max_weight_dfs_tree(G: nx.DiGraph, reference_path: list,
                        reference_paths: list[list] | None = None,
                        reference_owners=None, **kwargs) -> nx.DiGraph:
    """DFS forest prioritizing reference successors, then higher edge weights."""
    def neighbors(node, label):
        return [v for _, v, _ in sorted(G.edges(node, data=True),
                                       key=lambda x: (x[2]['weight'], x[1]))]
    return _reference_forest(G, reference_path, reference_paths, neighbors,
                             lambda u, v, label: label, reference_owners)

def _reference_successors(G: nx.DiGraph, reference_path: list) -> tuple[dict, list]:
    """
    Maps each reference node in G to the nodes that follow it on the reference walk, in walk order. A cyclic walk
    revisits a node or its complement; after such steps the walk resumes at its next new node, which is therefore
    also a successor when G has that edge. Returns the map and the indices of revisiting steps.
    """
    seen = set()
    revisit_steps = []
    for i, node in enumerate(reference_path):
        oriented = node[-2:] in ('_+', '_-')
        if node in seen or (oriented and node_complement(node) in seen):
            revisit_steps.append(i)
        seen.add(node)

    revisits = set(revisit_steps)
    next_new_node = [None] * len(reference_path)
    upcoming = None
    for i in range(len(reference_path) - 1, -1, -1):
        next_new_node[i] = upcoming
        if i not in revisits:
            upcoming = reference_path[i]

    successors = {}
    for i, node in enumerate(reference_path[:-1]):
        if node not in G:
            continue
        for candidate in (reference_path[i + 1], next_new_node[i]):
            if candidate is not None and G.has_edge(node, candidate):
                targets = successors.setdefault(node, [])
                if candidate not in targets:
                    targets.append(candidate)
    return successors, revisit_steps

def _warn_cyclic_reference(G: nx.DiGraph, reference_path: list, revisit_steps: list, dfs_tree: nx.DiGraph,
                           merge_distance: int = 100) -> None:
    def length(node):
        if node not in G and node[-2:] in ('_+', '_-'):
            node = node_complement(node)
        return len(G.nodes[node].get('sequence', '')) if node in G else 0

    ends = list(itertools.accumulate(length(u) for u in reference_path))
    regions = []
    for i in revisit_steps:
        if regions and i - regions[-1][1] <= merge_distance:
            regions[-1][1] = i
        else:
            regions.append([i, i])
    spans = ', '.join(f'{ends[a] - length(reference_path[a]) + 1}-{ends[b]}' for a, b in regions[:10])
    if len(regions) > 10:
        spans += ', ...'

    revisits = set(revisit_steps)
    departures = sum(1 for i in range(1, len(reference_path)) if i not in revisits
                     and reference_path[i - 1] not in dfs_tree.pred.get(reference_path[i], ()))
    num_inverted = sum(reference_path[i] not in G for i in revisit_steps)
    warnings.warn(
        f'Reference walk is cyclic: {len(revisit_steps)} steps revisit a node or its complement '
        f'({num_inverted} in the opposite orientation). The edges they use are treated as variant edges, and the '
        f'reference tree departs from the reference walk at {departures} nodes. Affected regions '
        f'({len(regions)}; bp along the reference walk, 1-based): {spans}'
    )

def assign_node_directions(G: nx.DiGraph,
                        reference_path: list, reference_paths: list[list] | None = None) -> None:
    """
    Perform a depth-first search (DFS) prioritizing higher weight edges and return a DFS tree.

    :param G: A NetworkX directed graph
    :param reference_path: A list of nodes in the reference path, beginning with source node
    :return: A DFS spanning tree of G
    """

    G_undirected = G.to_undirected(as_view=True)
    unvisited = set(G.nodes())
    stack = []

    def visit_node(parent, current_node, current_direction):

        if current_node not in unvisited:
            return

        # visit current_node
        unvisited.remove(current_node)
        G.nodes[current_node]['direction'] = current_direction
        unvisited.remove(node_complement(current_node))
        G.nodes[node_complement(current_node)]['direction'] = -current_direction
        if parent is not None:
            assert G_undirected.has_edge(parent, current_node)

        # Sort neighbors by descending edge weight
        neighbors = sorted(G_undirected.edges(current_node, data=True), key=lambda x: (x[2]['weight'], x[1]), reverse=False)

        # Last-added edges will be visited first
        for u, neighbor, _ in neighbors:
            stack.append((current_node, neighbor, current_direction))

    # Initialize the search at the linear reference path so that these nodes all have the same orientation
    for path in reference_paths or [reference_path]:
        for i, node in enumerate(path):
            visit_node(None if i == 0 else path[i - 1], node, 1)

    component_seeds = iter(sorted(unvisited))
    while True:
        while stack:
            parent, current_node, current_direction = stack.pop()
            visit_node(parent, current_node, current_direction)

        # Ensure that all connected components are visited
        num_unvisited = len(unvisited)
        if num_unvisited == 0:
            break
        print(f'Visiting another connected component, with at most {num_unvisited} nodes')
        # Scan the stable order once; visit_node removes the seed and its complement.
        unvisited_node = next(node for node in component_seeds if node in unvisited)
        visit_node(None, unvisited_node, 1)


def haplo_contiguous_dfs_tree(G: nx.DiGraph,
                        haplo_priorities: dict[str, int],
                        reference_path: list, reference_paths: list[list] | None = None,
                        reference_owners=None) -> nx.DiGraph:

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

        neighbors = sorted(priority.keys(), key=lambda x: (priority[x], x), reverse=True)
        return neighbors

    def next_label(parent, node, current):
        edge_data = G.edges[parent, node]
        if current in edge_data and current != 'other':
            return current
        return next((label for label in haplo_labels if label in edge_data), 'other')

    return _reference_forest(G, reference_path, reference_paths,
                             lambda node, label: get_neighbors(G, node, label),
                             next_label, reference_owners)


# Dictionary of available DFS methods, defined after all functions
dfs_methods: dict = {'max_weight': max_weight_dfs_tree,
                     'contiguous': haplo_contiguous_dfs_tree}

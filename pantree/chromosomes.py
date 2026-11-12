"""Reference-tree construction and coordinates for a whole-genome graph.

Each chromosome has its own terminal pair and rooted tree. The trees are stored
in one DiGraph (a forest); graph edges between trees remain variant edges.
"""
from collections import deque

import networkx as nx

from .dfs import assign_node_directions
from .utils import edge_complement, node_complement


def initialize_chromosomes(graph, walk_starts, walk_ends, dfs_method):
    paths = dict(sorted(graph.reference_paths.items()))
    assign_node_directions(graph, [], reference_paths=list(paths.values()))
    graph.terminal_binodes = set()
    graph.chromosome_termini = {
        chrom: (f'+_terminus_{i}', f'-_terminus_{i}')
        for i, chrom in enumerate(paths)
    }
    for pair in graph.chromosome_termini.values():
        graph.terminal_binodes.update(pair)
    if any(name + '_+' in graph for name in graph.terminal_binodes):
        raise ValueError('Input segment ID collides with a chromosome terminal ID')

    # First reference occurrence owns a shared segment and its coordinate.
    reference_owner = {}
    lengths = {}
    for chrom, path in paths.items():
        pos = graph.reference_starts[chrom]
        seen = set()
        for node in path:
            segment = node[:-2]
            if segment in seen:
                raise ValueError(f'Reference chromosome {chrom!r} contains a repeated segment')
            seen.add(segment)
            pos += len(graph.nodes[node]['sequence'])
            positive = graph.positive_node(node)
            if positive not in reference_owner:
                reference_owner[positive] = chrom
                for oriented in (positive, node_complement(positive)):
                    graph.nodes[oriented].update(
                        chromosome=chrom, position=pos, tree_position=pos,
                        on_reference_path=1)
        lengths[chrom] = pos

    # Hints only choose roots for branches that cannot be reached forwards from
    # a reference. Final ownership always comes from the selected tree parent.
    hints = dict(reference_owner)
    queue = deque(hints)
    while queue:
        node = queue.popleft()
        for neighbor in sorted(set(graph.successors(node)) | set(graph.predecessors(node))):
            neighbor = graph.positive_node(neighbor)
            if neighbor not in hints:
                hints[neighbor] = hints[node]
                queue.append(neighbor)
    fallback = next(iter(paths))

    def add_link(u, v):
        if not graph.has_edge(u, v):
            graph.add_biedge(u, v)

    seed_paths = []
    for chrom, path in paths.items():
        start, end = graph.chromosome_termini[chrom]
        for name, position in ((start, graph.reference_starts[chrom]), (end, lengths[chrom])):
            graph.add_binode(name)
            for suffix in ('_+', '_-'):
                graph.nodes[name + suffix].update(
                    chromosome=chrom, position=position, tree_position=position,
                    on_reference_path=1)
        # Shared reference segments belong to the first reference that uses
        # them. Each remaining contiguous stretch can still seed this tree.
        chunk = [start + '_+']
        for node in path:
            if reference_owner[graph.positive_node(node)] != chrom:
                if len(chunk) > 1:
                    add_link(chunk[-1], end + '_+')
                    seed_paths.append(chunk + [end + '_+'])
                    chunk = [start + '_+']
                continue
            add_link(chunk[-1], node)
            chunk.append(node)
        add_link(chunk[-1], end + '_+')
        seed_paths.append(chunk + [end + '_+'])

    positive = graph.subgraph([n for n in graph if graph.direction(n) == 1])
    sources = {n for n, degree in positive.in_degree() if degree == 0 and not graph.is_terminal(n)}
    sinks = {n for n, degree in positive.out_degree() if degree == 0 and not graph.is_terminal(n)}
    for node in walk_starts:
        (sources if graph.direction(node) == 1 else sinks).add(graph.positive_node(node))
    for node in walk_ends:
        (sinks if graph.direction(node) == 1 else sources).add(graph.positive_node(node))
    for nodes, is_source in ((sources, True), (sinks, False)):
        for node in sorted(nodes):
            chrom = hints.get(node, fallback)
            start, end = graph.chromosome_termini[chrom]
            add_link(start + '_+', node) if is_source else add_link(node, end + '_+')

    def build_tree():
        return dfs_method(positive, reference_path=[], reference_paths=seed_paths,
                          haplo_priorities=graph.haplo_priorities)

    graph.reference_tree = build_tree()
    missing = set(positive) - set(graph.reference_tree)
    if missing:
        # A component with no reference and no walk can contain source cycles.
        components = nx.condensation(positive.subgraph(missing))
        for component, degree in components.in_degree():
            if degree == 0:
                node = min(components.nodes[component]['members'])
                start, _ = graph.chromosome_termini[hints.get(node, fallback)]
                add_link(start + '_+', node)
        graph.reference_tree = build_tree()

    for node in nx.topological_sort(graph.reference_tree):
        parent = graph.parent_in_tree(node)
        if parent is not None:
            chrom = graph.chromosome(parent)
            graph.nodes[node]['chromosome'] = chrom
            graph.nodes[node_complement(node)]['chromosome'] = chrom

    # Walk endpoints use the final tree assignment, including walks crossing
    # chromosomes. Additional terminal edges are excluded from exported sites.
    for node in walk_starts:
        start, end = graph.termini_for(node)
        add_link(start + '_+' if graph.direction(node) == 1 else end + '_-', node)
    for node in walk_ends:
        start, end = graph.termini_for(node)
        add_link(node, end + '_+' if graph.direction(node) == 1 else start + '_-')

    for u, v, data in graph.edges(data=True):
        data['is_in_tree'] = (graph.reference_tree.has_edge(u, v)
                              or graph.reference_tree.has_edge(*edge_complement((u, v))))
        if graph.chromosome(u) != graph.chromosome(v):
            assert not data['is_in_tree']
            data['trans_chrom'] = graph.chromosome(v)
    graph.variant_edges = {(u, v) for u, v, data in graph.edges(data=True)
                           if data['is_representative'] and not data['is_in_tree']}

    # Compute LCAs within a chromosome. A cross-chromosome edge has no common
    # ancestor in the forest; its source is its annotation anchor.
    pairs_by_chrom = {chrom: [] for chrom in paths}
    for u, v in graph.variant_edges:
        if not graph.is_trans_chromosomal((u, v)):
            pairs_by_chrom[graph.chromosome(u)].append((graph.positive_node(u), graph.positive_node(v)))
    branch_points = {}
    for chrom, pairs in pairs_by_chrom.items():
        if pairs:
            root = graph.chromosome_termini[chrom][0] + '_+'
            branch_points.update(nx.tree_all_pairs_lowest_common_ancestor(
                graph.reference_tree, root=root, pairs=pairs))
    for edge in graph.variant_edges:
        u, v = edge
        complement = edge_complement(edge)
        if graph.is_trans_chromosomal(edge):
            graph.edges[edge]['branch_point'] = graph.positive_node(u)
            graph.edges[complement]['branch_point'] = graph.positive_node(complement[0])
            continue
        branch = branch_points[graph.positive_node(u), graph.positive_node(v)]
        if branch == v or branch == node_complement(u):
            graph.edges[edge]['is_back_edge'] = True
            graph.edges[complement]['is_back_edge'] = True
            branch = graph.parent_in_tree(branch)
        graph.edges[edge]['branch_point'] = branch
        graph.edges[complement]['branch_point'] = branch

    graph.compute_binode_indices()
    for node in graph.sorted_positive_nodes:
        if graph.on_reference_path(node):
            continue
        parent = graph.parent_in_tree(node)
        for oriented in (node, node_complement(node)):
            graph.nodes[oriented]['position'] = graph.nodes[parent]['position']
            graph.nodes[oriented]['tree_position'] = graph.nodes[parent]['tree_position'] + len(graph.nodes[node]['sequence'])
    for node in reversed(graph.sorted_positive_nodes):
        chrom = graph.chromosome(node)
        if graph.on_reference_path(node):
            right = graph.nodes[node]['position']
        else:
            successors = [v for v in graph.successors(node)
                          if graph.chromosome(v) == chrom and not graph.is_back_edge((node, v))
                          and graph.direction(v) == 1]
            right = min((graph.nodes[v]['right_position'] for v in successors), default=lengths[chrom])
        for oriented in (node, node_complement(node)):
            graph.nodes[oriented]['right_position'] = right
    graph.reference_paths = {chrom: [graph.chromosome_termini[chrom][0] + '_+'] + path
                             + [graph.chromosome_termini[chrom][1] + '_+'] for chrom, path in paths.items()}
    graph.reference_path = next(iter(graph.reference_paths.values()))
    graph.position.cache_clear()

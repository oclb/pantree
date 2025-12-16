from __future__ import annotations
from functools import lru_cache, cached_property
from math import inf
import networkx as nx
import numpy as np
from .utils import (
    node_complement,
    edge_complement,
    sequence_complement,
    walk_complement,
    node_recover,
)
from .gfa import read_gfa_line_by_line, GFAWalkLine
from .logger import setup_logger
import logging
from .dfs import assign_node_directions, dfs_methods
from .genotype import Genotype
from .vcf import write_vcf_from_graph
import os
import time
from collections import defaultdict
from tqdm import tqdm
from typing import Any, Optional, Callable
import warnings
from functools import lru_cache

class PangenomeGraph(nx.DiGraph):
    reference_tree: nx.classes.digraph.DiGraph
    reference_path: list[str]
    variant_edges: set
    haplo_priorities: dict[str, int]
    number_of_biedges: int  # Needed because nx.number_of_edges() runs in O(number of edges)
    logger: logging.Logger
    walks: dict[str, list[str]]

    # A valid walk proceeds from either +_terminus_+ to -_terminus_+ or from -_terminus_- to +_terminus_-
    @property
    def termini(self) -> tuple[str, str]:
        return '+_terminus', '-_terminus'

    def is_terminal(self, node_or_edge) -> bool:
        if isinstance(node_or_edge, str):
            return node_or_edge[:-2] in self.termini
        elif isinstance(node_or_edge, tuple):
            return self.is_terminal(self.edges[node_or_edge]['branch_point']) or self.is_terminal(node_or_edge[1])
        else:
            raise TypeError

    # Each biedge has a representative edge, whichever is in the .gfa file
    @property
    def sorted_biedge_representatives(self) -> list[str]:
        edges = [edge_with_data for edge_with_data in self.edges(data=True) if edge_with_data[2]['is_representative']]
        return sorted(edges, key=lambda edge: edge[2]['index'])

    @property
    def sorted_positive_nodes(self) -> list[str]:
        nodes = [node for node, val in self.nodes.items() if val['direction'] == 1]
        return list(sorted(nodes, key=lambda u: self.nodes[u]['index']))

    @lru_cache(maxsize=2)
    def sorted_variant_edges(self, exclude_terminus=True) -> list[str]:
        if exclude_terminus:
            sorted_vars = [edge for edge in self.variant_edges if not self.is_terminal(edge)]
        else:
            sorted_vars = self.variant_edges
        sorted_vars = sorted(sorted_vars, key=lambda x: self.position(self.edges[x]['branch_point']))
        return sorted_vars

    @property
    def biedge_attribute_names(self) -> tuple:
        return 'index', 'weight', 'is_in_tree', 'branch_point', 'is_back_edge'

    @property
    def node_attribute_names(self) -> tuple:
        return 'direction', 'sequence', 'position', 'right_position', 'tree_position', 'on_reference_path', 'index'

    @property
    def vcf_attribute_names(self) -> tuple:
        return 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'

    @property
    def number_of_binodes(self) -> int:
        return self.number_of_nodes() // 2

    def on_reference_path(self, node_or_edge):
        if type(node_or_edge) is tuple:
            if not self.has_edge(*node_or_edge):
                msg = f"Graph does not have edge {node_or_edge}"
                self.logger.error(msg)
                raise ValueError(msg)
            return self.nodes[node_or_edge[1]]['on_reference_path']
        elif type(node_or_edge) is str:
            if not self.has_node(node_or_edge):
                msg = f"Graph does not have node {node_or_edge}"
                self.logger.error(msg)
                raise ValueError(msg)
            return self.nodes[node_or_edge]['on_reference_path']

    def parent_in_tree(self, node: str) -> str:
        if self.reference_tree.in_degree(node) == 0:
            return None
        return next(self.reference_tree.predecessors(node))

    @lru_cache(maxsize=None)
    def position(self, node_or_edge: Any) -> Any:
        if isinstance(node_or_edge, str):
            return self.nodes[node_or_edge]['position']
        elif isinstance(node_or_edge, tuple):
            return self.position(node_or_edge[0]), self.position(node_or_edge[1])
        else:
            raise TypeError

    def right_position(self, node_or_edge: Any) -> Any:
        if isinstance(node_or_edge, str):
            return self.nodes[node_or_edge]['right_position']
        elif isinstance(node_or_edge, tuple):
            return self.right_position(node_or_edge[0]), self.right_position(node_or_edge[1])
        else:
            raise TypeError

    @classmethod
    def from_gfa(cls,
                 gfa_file: str,
                 ref_name: str = 'GRCh38',
                 logger: Optional[logging.Logger] = None,
                 return_walks: bool = False,
                 dfs_method: Callable[[dict[str, Any]], list[str]] = dfs_methods['max_weight'],
                 priority_dict: dict[str, int] | None = None,
                 ) -> PangenomeGraph:
        """
        Reads a .gfa file into a PangenomeGraph object.
        :param gfa_file: path to a file name ending in .gfa
        :param ref_name: name of the reference sample
        :param logger: Logger instance (if None, uses a null logger that does nothing)
        :param return_walks: if True, return (graph, walks) tuple
        :param dfs_method: DFS method for tree construction (function that takes a graph and returns a list of nodes)
        :param priority_dict: dict mapping sample names to priority integers (lower = higher priority)
        """
        # Use null logger if none provided
        if logger is None:
            logger = logging.getLogger('pantree')
            logger.addHandler(logging.NullHandler())
            logger.setLevel(logging.CRITICAL + 1)  # Disable all logging

        if priority_dict is None:
            priority_dict = {ref_name: 0}

        logger.info(f"Loading GFA file: {gfa_file}")
        logger.info(f"Parameters: ref_name={ref_name}, priority_dict={priority_dict}")

        if not os.path.exists(gfa_file):
            msg = f"GFA file not found: {gfa_file}"
            logger.error(msg)
            raise FileNotFoundError(msg)


        G = cls()
        G.haplo_priorities = {}
        G.logger = logger

        gfa_basename = os.path.basename(gfa_file)
        walk_start_nodes = []
        walk_end_nodes = []

        logger.info(f"Reading GFA file: {gfa_basename}")
        for line in read_gfa_line_by_line(gfa_file, line_types=['S']):
            G.add_binode(line.node_id, line.sequence)

        for line in read_gfa_line_by_line(gfa_file, line_types=['L']):
            G.add_biedge(line.u, line.v)


        if return_walks:
            G.walks = {}

        for line in read_gfa_line_by_line(gfa_file, line_types=['W', 'P']):
            # Check if haplotype name starts with ref_name (e.g., 'GRCh38#0#chr20' starts with 'GRCh38')
            hit_reference = line.sample_name == ref_name
            if hit_reference:
                if G.reference_path:
                    warnings.warn("Reference path already exists. Skipping additional reference walks.")
                else:
                    G.add_reference_path(line.walk)

            # Extract sample#haplotype from hap_name (format: sample#haplotype#contig)
            if line.sample_name in priority_dict:
                G.haplo_priorities[line.hap_name] = priority_dict[line.sample_name]
                G.update_haplotype_positions(line)

            if return_walks:
                assert G.walks is not None
                G.walks[line.hap_name] = line.walk

            walk_start_nodes.append(line.walk[0])
            walk_end_nodes.append(line.walk[-1])

            G.compute_edge_weights([line.walk])

        logger.info(f"Loaded {len(G.nodes) / 2:.0f} binodes, {len(G.edges) / 2:.0f} biedges")

        assign_node_directions(G, G.reference_path)

        G.add_terminal_nodes(walk_start_nodes=walk_start_nodes, walk_end_nodes=walk_end_nodes)
        logger.info("Computing reference tree")
        G.compute_reference_tree(dfs_method, haplo_priorities=priority_dict)

        logger.info("Computing branch points")
        G.annotate_branch_points()

        logger.info("Computing positions")
        G.compute_binode_indices()
        G.compute_binode_positions()
        G.compute_binode_right_positions()

        return G

    @classmethod
    def from_gfa_line_by_line(cls, *args, **kwargs):
        """Backward compatibility alias for from_gfa."""
        return cls.from_gfa(*args, **kwargs)

    def __init__(self,
                 directed_graph: nx.classes.digraph.DiGraph | None = None,
                 reference_tree: nx.classes.digraph.DiGraph | None = None,
                 reference_path: list[str] | None = None,
                 variant_edges: set | None = None,
                 haplo_priorities: dict[str, int] | None = None,
                 logger: logging.Logger | None = None,
                 walks: dict[str, list[str]] | None = None
                 ):
        if directed_graph is None:
            directed_graph = nx.DiGraph()
        if reference_tree is None:
            reference_tree = nx.DiGraph()
        super().__init__(directed_graph)
        self.reference_path = reference_path if reference_path is not None else []
        self.variant_edges = variant_edges if variant_edges is not None else set()
        self.haplo_priorities = haplo_priorities if haplo_priorities is not None else {}
        self.logger = logger if logger is not None else logging.getLogger('pantree')
        self.number_of_biedges = np.sum(
            [count_or_not for _, _, count_or_not in directed_graph.edges(data='is_representative')]
        )
        self.walks = walks if walks is not None else {}

    def identify_variant_type(self,
                              edge: tuple[str, str],
                              ref: str | None = None,
                              alt: str | None = None,
                              ) -> str:
        var_type = None
        if self.is_inversion(edge):
            if var_type is not None:
                raise KeyError(f'Variant, {edge}, has dual type.')
            var_type = 'INV'
        if self.is_replacement(edge, ref=ref, alt=alt):
            if var_type is not None:
                raise KeyError(f'Variant, {edge}, has dual type.')
            var_type = 'REP'
        if self.is_insertion(edge, ref=ref, alt=alt):
            if var_type is not None:
                raise KeyError(f'Variant, {edge}, has dual type.')
            var_type = 'INS'
        if self.is_snp(edge, ref=ref, alt=alt):
            if var_type is not None:
                raise KeyError(f'Variant, {edge}, has dual type.')
            var_type = 'SNP'
        if self.is_mnp(edge, ref=ref, alt=alt):
            if var_type is not None:
                raise KeyError(f'Variant, {edge}, has dual type.')
            var_type = 'MNP'
        if self.is_back_edge(edge):
            if var_type is not None:
                raise KeyError(f'Variant, {edge}, has dual type.')
            var_type = 'DUP'
        if self.is_forward_edge(edge):
            if var_type is not None:
                raise KeyError(f'Variant, {edge}, has dual type.')
            var_type = 'DEL'
        assert var_type is not None

        return var_type

    def is_inversion(self, edge: tuple[str, str]) -> bool:
        u, v = edge
        return self.direction(u) != self.direction(v)

    def is_in_tree(self, edge: tuple[str, str]) -> bool:
        return (edge in self.reference_tree or edge_complement(edge) in self.reference_tree)

    def is_back_edge(self, edge: tuple[str, str]) -> bool:
        if self.is_inversion(edge):
            return False
        return self.edges[edge]['is_back_edge']

    def is_forward_edge(self, edge: tuple[str, str]) -> bool:
        if self.is_inversion(edge):
            return False
        if self.is_in_tree(edge):
            return False
        positive_edge = edge if self.direction(edge[0]) == 1 else edge_complement(edge)
        branch_point = self.edges[positive_edge]['branch_point']
        return branch_point == positive_edge[0]

    def is_crossing_edge(self, edge: tuple[str, str]) -> bool:
        if self.is_inversion(edge):
            return False
        if self.is_in_tree(edge):
            return False
        if self.is_forward_edge(edge):
            return False
        if self.is_back_edge(edge):
            return False
        return True

    def is_insertion(self,
                     edge: tuple[str, str],
                     ref: str | None = None,
                     alt: str | None = None,
                     ) -> bool:
        if not self.is_crossing_edge(edge):
            return False
        if ref is None or alt is None:
            ref, _, _, _ = self.ref_alt_alleles(edge)
        return len(ref) == 0

    def is_replacement(self,
                       edge: tuple[str, str],
                       ref: str | None = None,
                       alt: str | None = None,
                       ) -> bool:
        if not self.is_crossing_edge(edge):
            return False
        if ref is None or alt is None:
            ref, alt, _, _ = self.ref_alt_alleles(edge)
        return len(ref) != 0 and len(alt) != 0 and len(ref) != len(alt)

    def is_snp(self,
               edge: tuple[str, str],
               ref: str | None = None,
               alt: str | None = None,
               ) -> bool:
        if not self.is_crossing_edge(edge):
            return False
        if ref is None or alt is None:
            ref, alt, _, _ = self.ref_alt_alleles(edge)
        return len(ref) == len(alt) == 1

    def is_mnp(self,
               edge: tuple[str, str],
               ref: str | None = None,
               alt: str | None = None,
               ) -> bool:
        if not self.is_crossing_edge(edge):
            return False
        if ref is None or alt is None:
            ref, alt, _, _ = self.ref_alt_alleles(edge)
        return len(ref) == len(alt) > 1

    def add_reference_path(self, reference_walk):
        self.reference_path = reference_walk
        for node in reference_walk:
            self.nodes[node]['on_reference_path'] = 1
            self.nodes[node_complement(node)]['on_reference_path'] = 1


    def add_terminal_nodes(self, walk_start_nodes: list[str]=None, walk_end_nodes: list[str]=None):
        """Add two terminal binodes, +_terminus and -_terminus, to the graph. A valid walk proceeds
        from +_terminus_+ to -_terminus_+ or from -_terminus_- to +_terminus_-."""

        positive_subgraph = self.subgraph([n for n, direction in self.nodes(data="direction") if direction == 1])

        # source and sink nodes for positive-direction walks; complementary nodes are sink and source nodes for
        # negative-direction walks, respectively
        source_nodes = {node for node, degree in positive_subgraph.in_degree() if degree == 0}
        sink_nodes = {node for node, degree in positive_subgraph.out_degree() if degree == 0}

        if walk_start_nodes:
            source_nodes = source_nodes.union(
                {node for node in walk_start_nodes if self.direction(node) == 1}
            )
            sink_nodes = sink_nodes.union(
                {node_complement(node) for node in walk_start_nodes if self.direction(node) == -1}
            )
        if walk_end_nodes:
            sink_nodes = sink_nodes.union(
                {node for node in walk_end_nodes if self.direction(node) == 1}
            )
            source_nodes = source_nodes.union(
                {node_complement(node) for node in walk_end_nodes if self.direction(node) == -1}
            )

        plus_terminus, minus_terminus = self.termini
        self.add_binode(plus_terminus)
        self.add_binode(minus_terminus)

        for source_node in source_nodes:
            self.add_biedge(plus_terminus + '_+', source_node)

        for sink_node in sink_nodes:
            self.add_biedge(sink_node, minus_terminus + '_+')

        if self.reference_path == []:
            self.add_biedge(plus_terminus + '_+', minus_terminus + '_+')

        # reference path is assumed to be in the positive direction
        self.reference_path = [plus_terminus + '_+'] + self.reference_path + [minus_terminus + '_+']
        self.nodes[plus_terminus + '_+']['on_reference_path'] = 1
        self.nodes[plus_terminus + '_-']['on_reference_path'] = 1
        self.nodes[minus_terminus + '_+']['on_reference_path'] = 1
        self.nodes[minus_terminus + '_-']['on_reference_path'] = 1

        chromosome_length = self.nodes[self.reference_path[-2]]['position']
        self.nodes[plus_terminus + '_+']['position'] = 0
        self.nodes[plus_terminus + '_-']['position'] = 0
        self.nodes[minus_terminus + '_+']['position'] = chromosome_length
        self.nodes[minus_terminus + '_-']['position'] = chromosome_length

    def compute_reference_tree(self, dfs_method: Callable = dfs_methods['max_weight'], **kwargs):
        """
        Computes the reference tree, a DFS spanning tree of the positively-oriented subgraph; defines variant edges
        as those that are not in the reference tree or its complement.
        """
        # reference tree contains positive-direction nodes only, and no inversion edges
        # self.add_biedge(self.reference_path[0], self.reference_path[1])
        positive_subgraph = self.subgraph([n for n, direction in self.nodes(data="direction") if direction == 1])
        self.reference_tree = dfs_method(positive_subgraph, reference_path=self.reference_path, **kwargs)

        for edge in positive_subgraph.edges():
            edge_in_tree = self.reference_tree.has_edge(*edge)
            self.edges[edge]['is_in_tree'] = edge_in_tree
            self.edges[edge_complement(edge)]['is_in_tree'] = edge_in_tree

        # Variant edges are those not in the tree
        self.variant_edges = {(u, v) for u, v, data in self.edges(data=True)
                              if data['is_representative'] and not data['is_in_tree']}

    def get_reference_edges(self) -> dict[tuple[str, str], tuple[str, str]]:
        """Returns a dictionary mapping variant edges to their references edges."""
        return {e: self.representative_edge(self.reference_tree_edge(e)) for e in self.variant_edges}

    def write_vcf(self,
                  gfa_path: Optional[str],
                  vcf_filename: str,
                  chr_name: str,
                  exclude_terminus: bool = True,
                  size_threshold: float = float('inf'),
                  check_degenerate: bool = False,
                  ) -> None:
        """
        Writes the variant call format (vcf) file.
        :param gfa_path: the .gfa file, from which walks are read, or None to skip writing genotypes
        :param vcf_filename: the output vcf file path
        :param chr_name: the chromosome name in the first column of output vcf file
        :param size_threshold: the truncation length of ref and alt sequence
        :param check_degenerate: whether to exclude variants whose ref and alt alleles are identical
        :param log_path: path to log file
        :return:
        """
        # Use graph's logger, or create one if log_path provided
        write_vcf_from_graph(
            graph=self,
            gfa_path=gfa_path,
            vcf_filename=vcf_filename,
            chr_name=chr_name,
            logger=self.logger,
            exclude_terminus=exclude_terminus,
            size_threshold=size_threshold,
            check_degenerate=check_degenerate,
        )

    def genotype(self, walk: list[str], exclude_terminus: bool=True) -> Genotype:
        """
        Compute the phased genotype data for a single walk.
        """
        return Genotype.genotype(self, walk, exclude_terminus)

    def genotypes_from_gfa(self, gfa_path: str, exclude_terminus: bool=True) -> dict[str, tuple[Genotype]]:
        """
        Compute the phased genotype data for each sample in a GFA file. There can be 1 or 2 haplotypes
        per sample, which are stored in a tuple. There can be any positive number of walks per haplotype,
        which are aggregated together.
        """
        self.logger.info("Computing genotype for haplotypes")
        genotype_dict = {}
        names = defaultdict(set)
        for line in read_gfa_line_by_line(gfa_path, line_types=['W', 'P']):
            walk = line.walk
            haplotype_name = '#'.join(line.hap_name.split('#')[:2])
            names[line.sample_name].add(haplotype_name)
            genotype: Genotype = self.genotype(walk, exclude_terminus)
            if haplotype_name in genotype_dict:
                genotype_dict[haplotype_name].update(genotype)
            else:
                genotype_dict[haplotype_name] = genotype

        for _, genotype in genotype_dict.items():
            genotype.compute_missing_variants(self)

        result = {}
        for sample_name, haplo_names in names.items():
            assert len(haplo_names) in (1,2), "Samples should be haploid or diploid"
            result[sample_name] = tuple(genotype_dict[name] for name in haplo_names)

        return result

    def add_binode(self, binode: str, seq: str = ''):
        """
        Adds a binode to the bidirected graph, comprising a node and its complement.
        """
        node_data = {key: 0 for key in self.node_attribute_names}
        node_data['sequence'] = seq
        node_data['direction'] = 1

        self.add_node(str(binode) + '_+', **node_data)

        node_data['sequence'] = sequence_complement(seq)
        node_data['direction'] = -1

        self.add_node(str(binode) + '_-', **node_data)

    def add_biedge(self, node1: str, node2: str, *, weight: int = 0):
        """
        Adds a biedge to the bidirected graph, comprising an edge and its complement.
        """
        if self.has_edge(node1, node2):
            msg = f'Attempted to add duplicate biedge: {node1}, {node2}'
            self.logger.error(msg)
            raise ValueError(msg)

        edge_data = {key: 0 for key in self.biedge_attribute_names}
        edge_data['weight'] = weight
        edge_data['index'] = self.number_of_biedges
        edge_data['is_representative'] = True
        self.add_edge(node1, node2, **edge_data)

        edge_data['is_representative'] = False
        self.add_edge(node_complement(node2), node_complement(node1), **edge_data)
        self.number_of_biedges += 1

    def representative_edge(self, edge: tuple):
        """
        Returns the edge associated with a biedge which is in the .gfa file.
        :param edge: one of the two edges in the biedge
        :return: the one which is in the .gfa file
        """
        return edge if self.edges[edge]['is_representative'] else edge_complement(edge)

    def reference_tree_edge(self, variant_edge):
        _, v = self.positive_variant_edge(variant_edge)
        if self.is_inversion(variant_edge):
            v = self.positive_node(v)
        w = self.parent_in_tree(v)
        return w, v

    def positive_variant_edge(self, edge: tuple):
        return edge if self.direction(edge[0]) == 1 or self.is_inversion(edge) else edge_complement(edge)

    def update_haplotype_positions(self, line: GFAWalkLine) -> None:
        """
        For each edge (u->v) in the walks, store:
        primary_edge_pos = ((hap_id, contig_name, walk_id), position_bp)
        Keep only the tuple from the highest-priority hap (smaller hap_rank wins).
        position_bp is the cumulative bp offset at the START of (u->v) within that walk.
        """

        walk = line.walk
        hap_name = line.hap_name
        offset = line.contig_start

        for u in walk:
            offset += len(self.nodes[u]['sequence'])
            self.nodes[u][hap_name] = offset
            self.nodes[node_complement(u)][hap_name] = offset

    def compute_edge_weights(self, walks: list[list[str]]):
        """
        Computes the number of times that each edge is visited by some walk.
        """
        for walk in walks:
            for u, v in zip(walk[:-1], walk[1:]):
                self.edges[u, v]['weight'] += 1
                self.edges[node_complement(v), node_complement(u)]['weight'] += 1

    def annotate_branch_points(self) -> None:
        """
        Computes the branch point, i.e. the lowest common ancestor in the reference tree, of each variant biedge. The
        branch point is a positive orientation node. The edges (u,v), (u complement, v), (v,u), etc., all have the same
        branch point.
        """

        positive_direction_variants = [(self.positive_node(u), self.positive_node(v)) for u,v in self.variant_edges]

        # For each variant edge (u,v), yield the lowest common ancestor of u and v in the tree
        branch_point_tuples = nx.tree_all_pairs_lowest_common_ancestor(
            self.reference_tree,
            root=self.termini[0]+'_+',
            pairs=positive_direction_variants
        )
        branch_points = dict(branch_point_tuples)

        for edge in self.variant_edges:
            u, v = edge
            branch_point = branch_points[self.positive_node(u), self.positive_node(v)]
            if branch_point == v or branch_point == node_complement(u):
                assert self.reference_tree.in_degree(branch_point) == 1
                self.edges[edge]['is_back_edge'] = True
                self.edges[edge_complement(edge)]['is_back_edge'] = True
                branch_point = self.parent_in_tree(branch_point)
            self.edges[edge]['branch_point'] = branch_point
            self.edges[edge_complement(edge)]['branch_point'] = branch_point

    def walk_up_tree(self, ancestor, descendant) -> list:
        u = descendant
        result = []
        while u != ancestor:
            result.append(u)
            u = self.parent_in_tree(u)
        result.append(ancestor)
        return result

    def positive_node(self, node: str) -> str:
        return node if self.direction(node) == 1 else node_complement(node)

    def direction(self, node: str) -> int:
        return self.nodes[node]['direction']


    def walk_sequence(self, walk: list[str]) -> str:
        seq = ''
        for node in walk:
            seq += self.nodes[node]['sequence']
        return seq

    def walk_with_variants(self, first: str, last: str, variant_edges: list) -> list:
        """
        Computes a walk from first node to last node, including all of the variant edges in the list.
        The walk proceeds from the 'end' of the first node to the 'start' of the last node, where the sequence associated
        with each node lies between its 'start' and 'end'. Thus, the first node is excluded from the walk, and the
        last node is usually excluded. However, the last node is included if the walk goes from + to -, because the
        'start' of the - node is the 'end' of its complementary + node.

        :param first: first node of the walk, either + or -
        :param last: last node
        :param variant_edges: variant edges in the order they are encountered, with the correct orientations
        :return: the nodes on the walk in order, excluding first and last
        """
        walk_direction = self.direction(first)
        source_nodes = [first] + [v for _, v in variant_edges]
        sink_nodes = [u for u, _ in variant_edges] + [last]
        result = []
        for source, sink in zip(source_nodes, sink_nodes):
            include_source = self.direction(source) == walk_direction
            source = source if self.direction(source) == walk_direction else node_complement(source)
            sink = sink if self.direction(source) == walk_direction else node_complement(source)
            pair = (source, self.positive_node(sink)) if walk_direction == 1 else \
                (self.positive_node(sink), self.positive_node(source))
            segment = self.walk_up_tree(*pair)
            if not include_source:
                segment = segment[:-1] if walk_direction == 1 else segment[1:]
            segment = segment[::-1]
            result += segment if walk_direction == 1 else walk_complement(segment)

        # exclude first and sometimes last nodes
        if self.direction(last) == -1 and walk_direction == 1:
            return result[1:]
        return result[1:-1]

    def ref_alt_alleles(self,
                        variant_edge: tuple
                        ):
        """
        Computes the reference allele and alternative allele of the branch point for each variant edge.
        :param variant_edges: list of tuples (u,v).
        :return: dict mapping variant edges to tuples (ref, alt)
        """

        if self.is_in_tree(variant_edge):
            msg = "Ref and alt alleles are only defined for variant edges"
            self.logger.error(msg)
            raise ValueError(msg)

        u, v = variant_edge
        branch_point = self.edges[u, v]['branch_point']
        first, last = (branch_point, v) if self.direction(u) == 1 else (u, branch_point)

        ref_path = self.walk_with_variants(first, last, [])
        alt_path = self.walk_with_variants(first, last, [variant_edge])

        alt_allele = self.walk_sequence(alt_path)
        ref_allele = self.walk_sequence(ref_path)

        branch_sequence = self.nodes[branch_point]['sequence']
        if not branch_sequence:
            branch_sequence = 'N'
        last_letter_of_branch_point = branch_sequence[-1]

        return ref_allele, alt_allele, last_letter_of_branch_point, branch_point


    def count_edge_visits(self, genotype: dict) -> dict:
        """
        Computes the number of time that a walk visits every edge, given the number of times it visits every
        variant edge.
        :param genotype: for every variant edge that is visited, the number of visits
        :return: for every edge that is visited, the number of visits
        """

        sinks: list[str] = []
        sources: dict[str, int] = {}

        # Add variant edge endpoints as sources or sinks depending on their respective directions
        for variant_edge, visit_count in genotype.items():
            if variant_edge not in self.variant_edges:
                msg = "geno dictionary contains a key which is not a variant edge"
                self.logger.error(msg)
                raise ValueError(msg)

            u, v = variant_edge
            new_sinks = []
            new_sources = []
            if self.direction(u) == 1:
                new_sinks.append(u)
            else:
                new_sources.append(node_complement(u))
            if self.direction(v) == 1:
                new_sources.append(v)
            else:
                new_sinks.append(node_complement(v))

            for w in new_sources:
                sources[w] = sources.get(w, 0) + visit_count

            sinks += visit_count * new_sinks

        # Add sink and source nodes for the beginning and end of the walk, depending on the number of inversions
        num_sources_minus_sinks = np.sum([val for _, val in sources.items()]) - len(sinks)

        # Equal number of + to - and - to + inversions: walk from + terminus to - terminus
        if num_sources_minus_sinks == 0:
            sources[self.termini[0] + '_+'] = 1
            sinks += [self.termini[1] + '_+']

        # Odd number of inversions, with one more from + to - strand: walk from - to -
        elif num_sources_minus_sinks == 2:
            sinks += 2 * [self.termini[1] + '_+']

        # Odd number of inversions, with one more from - to + strand: walk from + to +
        elif num_sources_minus_sinks == -2:
            sources[self.termini[0] + '_+'] = 2

        else:
            msg = "The input genotype does not correspond to any valid walk"
            self.logger.error(msg)
            raise ValueError(msg)

        edge_visits = genotype.copy()
        for sink in sinks:
            # Walk up the tree (in the only possible direction) until reaching a source
            current_node = sink
            while sources.get(current_node, 0) == 0:
                # If current_node is the root, it means that the input genotype was invalid
                if self.reference_tree.in_degree(current_node) == 0:
                    msg = "The input genotype does not correspond to any valid walk"
                    self.logger.error(msg)
                    raise ValueError(msg)

                previous_node = current_node
                current_node = self.parent_in_tree(current_node)
                edge_representative = self.representative_edge((current_node, previous_node))
                edge_visits[edge_representative] = edge_visits.get(edge_representative, 0) + 1

            # when reaching the source, it is "used up"
            sources[current_node] -= 1

        return edge_visits

    def allele_length(self) -> dict:
        """
        Computes the allele length for each variant edge.
        :return: dictionary mapping variant edges to its own allele length, len(ref_allele) + len(alt_allele)
        """
        def _allele_length(variant_edge):
            ref_allele, alt_allele, _, _ = self.ref_alt_alleles(variant_edge)
            return len(ref_allele) + len(alt_allele)

        return {e: _allele_length(e) for e in self.sorted_variant_edges(exclude_terminus=False)}

    def compute_binode_indices(self):
        """
        Computes the topological sort of the binodes according to the DAG whose nodes are positive-direction
        nodes and whose edges are not back edges nor inversions.
        """
        positive_subgraph = self.subgraph([n for n, direction in self.nodes(data="direction") if direction == 1])
        positive_subgraph_no_cycle = self.edge_subgraph([edge for edge in positive_subgraph.edges() if not self.is_back_edge(edge)])
        order = nx.topological_sort(positive_subgraph_no_cycle)
        for i, node in enumerate(order):
            self.nodes[node]['index'] = i
            self.nodes[node_complement(node)]['index'] = i



    def compute_binode_positions(self):
        """
        Computes the position of each binode along the linear reference path, as well as the tree position
        (distance from reference plus the position of the node's ancestor on the reference path), in basepairs.
        """
        current_position = 0
        for u in self.reference_path:
            current_position += len(self.nodes[u]['sequence'])
            self.nodes[u]['position'] = current_position
            self.nodes[u]['tree_position'] = current_position
            self.nodes[node_complement(u)]['position'] = current_position
            self.nodes[node_complement(u)]['tree_position'] = current_position

        for u in self.sorted_positive_nodes:
            if self.nodes[u]['on_reference_path']:
                continue
            predecessor = self.parent_in_tree(u)
            assert predecessor is not None
            self.nodes[u]['position'] = self.nodes[predecessor]['position']
            self.nodes[node_complement(u)]['position'] = self.nodes[u]['position']

            self.nodes[u]['tree_position'] = (self.nodes[predecessor]['tree_position'] +
                                              len(self.nodes[u]['sequence']))
            self.nodes[node_complement(u)]['tree_position'] = self.nodes[u]['tree_position']

    def compute_binode_right_positions(self):
        """Computes the right position of each binode, defined as the minimum position of its successors in the
        positive subgraph minus back edges."""
        order = self.sorted_positive_nodes
        for u in reversed(list(order)):
            if self.on_reference_path(u):
                self.nodes[u]['right_position'] = self.nodes[u]['position']
                self.nodes[node_complement(u)]['right_position'] = self.nodes[node_complement(u)]['position']
                continue
            successor_positions = [self.right_position(v) for v in self.successors(u)]
            self.nodes[u]['right_position'] = np.min(successor_positions)
            self.nodes[node_complement(u)]['right_position'] = self.nodes[u]['right_position']

    def get_variants_at_interval(self, half_open_interval: tuple[int, int], exclude_root_edges=True) -> list:
        """
        Computes variant edges whose position intersects some interval, by iterating over all variant edges.
        :param half_open_interval: [start, end) of the interval
        :param exclude_root_edges: if True, edges beginning or ending at a terminus are excluded
        :return: list of edges
        """
        start, end = half_open_interval
        result = []

        for edge in self.variant_edges:
            u, v = edge
            positions = [self.nodes[u]['position'], self.nodes[v]['position']]
            x, y = min(positions), max(positions)
            if x < end and y >= start:
                result.append(edge)

        if exclude_root_edges:
            result = [(u,v) for u,v in result if u != '+_terminus_+' and v != '+_terminus_-']

        return result

    def classify_triallelic_bubble(self, endpoint_binodes: list, variants: list) -> str:
        """
        Classifies a top-level triallelic superbubble containing 2 variants as either:
        overlapping: one pair of the three alleles shares one subsequence
        interlocking: two different pairs of alleles share one subsequence each
        nested: one pair of alleles shares two noncontiguous subsequences
        properly_triallelic: no pair of alleles shares a subsequence
        # TODO unsure if this handles inversions properly

        :param endpoint_binodes: starting and ending binodes of the superbubble
        :param variants: list of variants within the superbubble; there must be two of them
        :return: the class of superbubble
        """

        if len(variants) != 2:
            msg = "Expected exactly 2 variants in the variant list"
            self.logger.error(msg)
            raise ValueError(msg)

        if any([self.is_back_edge(edge) for edge in variants]):
            return "not_triallelic"

        if self.is_inversion(variants[0]) != self.is_inversion(variants[1]):
            msg = "Found an inversion and a non-inversion in the variant list"
            self.logger.error(msg)
            raise ValueError(msg)

        endpoint_nodes = [node + '_+' for node in endpoint_binodes]
        endpoint_nodes += [node + '_-' for node in endpoint_binodes]
        if not all([self.has_node(node) for node in endpoint_nodes]):
            msg = "One or more of the endpoint nodes is not in the graph"
            self.logger.error(msg)
            raise ValueError(msg)

        # Detect which node of each binode demarcates the superbubble
        variant_branch_points = [self.edges[e]['branch_point'] for e in variants]
        assert not any([brach_point == 0 for brach_point in variant_branch_points]), "Branch point of variant not found"
        start_nodes = [node for node in endpoint_nodes if node in variant_branch_points]
        if len(start_nodes) == 0:
            return "not_triallelic"
        variant_end_points = []
        for u, v in variants:
            if self.direction(u) == 1 and self.direction(v) == 1:
                variant_end_points.append(v)
            elif self.direction(u) == -1 and self.direction(v) == -1:
                variant_end_points.append(node_complement(u))
            elif self.direction(u) == 1 and self.direction(v) == -1:
                variant_end_points.append(u)
                variant_end_points.append(node_complement(v))
            elif self.direction(u) == -1 and self.direction(v) == 1:
                variant_end_points.append(node_complement(u))
                variant_end_points.append(v)
            else:
                raise ValueError("Invalid direction for variant edge nodes.")

        end_nodes = [node for node in endpoint_nodes if node in variant_end_points]
        if len(end_nodes) == 0:
            return "not_triallelic"

        assert len(start_nodes) == 1 and len(end_nodes) == 1, \
            f"Found {len(start_nodes)} possible start nodes, and {len(end_nodes)} possible end nodes"

        start_node = start_nodes[0]
        end_node = end_nodes[0]

        # Get in-neighbors and out-neighbors excluding specific nodes
        start_degree = len({successor for successor in self.successors(start_node) if not self.is_terminal(successor)})
        end_degree = len({predecessor for predecessor in self.predecessors(end_node) if not self.is_terminal(predecessor)})

        assert start_degree <= 3 and end_degree <= 3, \
            f"Starting and ending nodes ({start_node} and {end_node}) of the bubble had degree {start_degree} and {end_degree}"

        if start_degree == 3 and end_degree == 3:
            return 'properly_triallelic'

        if start_degree == 3 or end_degree == 3:
            return 'overlapping'

        # e.g., [(0,1), (0,2), (1,2), (1,3), (2,3)] with variant edges [(1,2), (1,3)]
        if all([node == start_node for node in variant_branch_points]):
            return 'interlocking'

        for branch_point, end_point in zip(variant_branch_points, variant_end_points):
            if branch_point == start_node and end_point == end_node:
                return 'nested'

        # e.g., [(0,1), (0,2), (1,2), (1,3), (2,3)] with variant edges [(0,2), (1,3)]
        return 'interlocking'



    def _match_sequence_up_tree(self, sequence: str, node: str) -> bool:
        """Returns True if `sequence` is a suffix of the unique sequence beginning at a terminus and
        ending at `node` within the tree.
        """
        node_sequence = self.nodes[node]['sequence']
        putative_match_length = min(len(sequence), len(node_sequence))
        if self.nodes[node]['sequence'][-putative_match_length:] != sequence[-putative_match_length:]:
            return False
        if putative_match_length == len(sequence):
            return True
        remaining_sequence = sequence[:putative_match_length]
        return self._match_sequence_up_tree(remaining_sequence, self.parent_in_tree(node))

    def _match_sequence_down_tree(self, sequence: str, node: str) -> bool:
        """Returns True if `sequence' is a prefix of some sequence beginning at `node` in
        thre tree."""
        node_sequence = self.nodes[node]['sequence']
        putative_match_length = min(len(sequence), len(node_sequence))
        if self.nodes[node]['sequence'][:putative_match_length] != sequence[:putative_match_length]:
            return False
        if putative_match_length == len(sequence):
            return True
        remaining_sequence = sequence[putative_match_length:]
        tree_successors = self.reference_tree.successors(node)
        return any([self._match_sequence_down_tree(remaining_sequence, successor) \
            for successor in tree_successors])

    def _match_sequence_down_graph(self, sequence: str, node: str, end_node: str) -> bool:
        """Returns True if `sequence' is a prefix of some sequence beginning at `node` in
        the graph."""
        node_sequence = self.nodes[node]['sequence']
        putative_match_length = min(len(sequence), len(node_sequence))
        if self.nodes[node]['sequence'][:putative_match_length] != sequence[:putative_match_length]:
            return False
        if putative_match_length == len(sequence):
            return (node == end_node and len(sequence) == len(node_sequence))
        remaining_sequence = sequence[putative_match_length:]
        tree_successors = self.successors(node)
        return any([self._match_sequence_down_graph(remaining_sequence, successor, end_node) \
            for successor in tree_successors])

    def annotate_repeat_motif(self,
                              variant_edge: tuple[str, str],
                              ref_allele: str = None,
                              alt_allele: str = None,
                              branch_point: str = None) -> Optional[str]:
        """
        Returns the repeat motif of a variant edge if it is a repeat.
        Otherwise, returns None.
        """
        import sys
        recursion_limit = max(sys.getrecursionlimit(), len(self.reference_tree.nodes))
        sys.setrecursionlimit(recursion_limit)

        if self.is_inversion(variant_edge):
            return None
        variant_edge = self.positive_variant_edge(variant_edge)
        _, v = variant_edge
        if ref_allele is None or alt_allele is None or branch_point is None:
            ref_allele, alt_allele, _, branch_point = self.ref_alt_alleles(variant_edge)

        def get_repeat_motif(allele: str) -> Optional[str]:
            non_basepair_character = 'N'
            if any(letter == non_basepair_character for letter in allele):
                return None

            allele_length = len(allele)
            for repeat_length in range(1,allele_length+1):
                if allele_length % repeat_length != 0:
                    continue
                motif = allele[:repeat_length]
                if allele == motif * (allele_length // repeat_length):
                    break

            if self._match_sequence_up_tree(motif, branch_point):
                return motif
            if not self.is_back_edge(variant_edge):
                if self._match_sequence_down_tree(motif, v):
                    return motif
            else:
                return motif
            return None


        if len(alt_allele) == 0:
            return get_repeat_motif(ref_allele)
        if len(ref_allele) == 0:
            return get_repeat_motif(alt_allele)

    def missing_inversion_allele(self, variant_edge: tuple[str, str], minimum_alt_length: int=10) -> Optional[str]:
        """Checks if the alt allele of a variant edge matches the reverse complement of some other path from branch
        point to v, whether that be the reference allele or a different alt path. Returns the matching allele
        or None."""
        if not self.is_crossing_edge(variant_edge):
            return None
        variant_edge = self.positive_variant_edge(variant_edge)
        ref, alt, _, branch_point = self.ref_alt_alleles(variant_edge)
        if len(alt) < minimum_alt_length:
            return None

        end_node = variant_edge[1]
        alt_complement = self.nodes[branch_point]['sequence'] \
            + sequence_complement(alt) + self.nodes[end_node]['sequence']

        if self._match_sequence_down_graph(alt_complement, branch_point, end_node):
            return alt

        return None

    def simplify_subgraph(self,
                          pos_range: tuple = None,
                          endpoints: set[set] = None,
                          minimum_allele_length: int = 1000) -> nx.DiGraph:
        """Returns a simplified minor of the underlying directed graph, with variant
        edges remvoed if their combined allele length is less than the minimum, with
        'tips' removed, and with 'paths' contracted."""

        subgraph = self.delete_small_variants(pos_range, minimum_allele_length)
        assert not nx.is_empty(subgraph), "Invalid position range: Empty subgraph."
        if endpoints is None:
            node_in_degrees = {node: degree for node, degree in subgraph.in_degree()}
            node_out_degrees = {node: degree for node, degree in subgraph.out_degree()}
            endpoints = {node for node in subgraph.nodes
                         if self.nodes[node]['on_reference_path'] == 1
                         and (node_in_degrees[node] == 0 or node_out_degrees[node] == 0)}
            assert len(endpoints) == 4
        self.delete_tips(subgraph, endpoints)
        self.contract_paths(subgraph, endpoints)

        return subgraph

    def delete_small_variants(self, pos_range: tuple = None, minimum_allele_length: int = 5) -> nx.DiGraph:
        """Returns an edge subgraph of the underlying directed graph with variant edges removed
        if their combined allele length is less than the minimum."""
        if pos_range == None:
            simplified_graph = nx.DiGraph(self)
            variant_edge_set = self.variant_edges
        else:
            nodes_to_include = [node for node, data in self.nodes(data=True)
                                if (data['position'] >= pos_range[0] and data['position'] < pos_range[1]) or
                                (data['position'] == pos_range[1] and data['on_reference_path'] == 1)]
            simplified_graph = nx.DiGraph(self.subgraph(nodes_to_include))
            variant_edge_set = [variant for variant in self.variant_edges
                                 if min(self.position(variant)) >= pos_range[0] and max(self.position(variant)) <= pos_range[1]]

        small_variant_edges = []
        for variant_edge in variant_edge_set:
            ref, alt, _, _ = self.ref_alt_alleles(variant_edge)
            if len(ref) + len(alt) < minimum_allele_length:
                small_variant_edges.append(variant_edge)
                small_variant_edges.append(edge_complement(variant_edge))

        simplified_graph.remove_edges_from(small_variant_edges)

        return simplified_graph

    @staticmethod
    def delete_tips(simplified_graph: nx.DiGraph, end_points: set[str]) -> None:
        node_in_degrees = {node: degree for node, degree in simplified_graph.in_degree()}
        node_out_degrees = {node: degree for node, degree in simplified_graph.out_degree()}
        tips = []
        for node, in_degree in node_in_degrees.items():
            out_degree = node_out_degrees[node]
            if in_degree == 1 and out_degree == 0 and node not in end_points:
                tips.append(node)

        nodes_to_delete = []
        for tip in tips:
            node = tip
            while node_in_degrees[node] == 1 and node_out_degrees[node] == 0:
                parent = next(simplified_graph.predecessors(node))
                node_out_degrees[parent] -= 1
                nodes_to_delete.append(node)
                nodes_to_delete.append(node_complement(node))
                node = parent

        simplified_graph.remove_nodes_from(nodes_to_delete)


    @staticmethod
    def contract_paths(simplified_graph: nx.DiGraph, end_points: set[str]) -> dict[str, str]:
        """
        Contracts paths in a graph by combining nodes.
        :param simplified_graph: A NetworkX directed graph with tips deleted and small variants removed
        :param end_points: A set of end points in the graph, which won't be combined with other nodes
        :return: A dictionary mapping each node that was combined with another node to the node
        with which it was combined
        """

        def combine_nodes(parent: str, child: str) -> None:
            neighbors = list(simplified_graph.successors(child))
            sequence = simplified_graph.nodes[child]['sequence']
            simplified_graph.remove_node(child)
            simplified_graph.remove_node(node_complement(child))
            for neighbor in neighbors:
                simplified_graph.add_edge(parent, neighbor)
                simplified_graph.add_edge(*edge_complement((parent, neighbor)))
            simplified_graph.nodes[parent]['sequence'] += sequence
            simplified_graph.nodes[node_complement(parent)]['sequence'] = sequence + \
                simplified_graph.nodes[node_complement(parent)]['sequence']


        node_in_degrees = {node: degree for node, degree in simplified_graph.in_degree()}
        node_out_degrees = {node: degree for node, degree in simplified_graph.out_degree()}

        starting_points = list(end_points)
        for node, out_degree in node_out_degrees.items():
            in_degree = node_in_degrees[node]
            if simplified_graph.nodes[node]['direction'] == 1 and (out_degree > 1 or in_degree > 1):
                starting_points.append(node)

        node_mapping: dict[str, str] = {}
        for start in starting_points:
            for node in list(simplified_graph.successors(start)):
                parent = start
                combined_nodes = []
                while node_in_degrees[node] == 1 and node_out_degrees[node] == 1:
                    if node_out_degrees[parent] == 1:
                        combine_nodes(parent, node) # TODO replace with combine_path
                        combined_nodes.append(node)
                        node = next(simplified_graph.successors(parent))
                    else:
                        parent = node
                        node = next(simplified_graph.successors(node))

                for combined_node in combined_nodes:
                    node_mapping[combined_node] = parent

        return node_mapping

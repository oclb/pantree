from __future__ import annotations
from dataclasses import dataclass, field
from bisect import bisect_left
from math import inf
import logging
import sys
import numpy as np
from .utils import edge_complement, get_from_biedge_dict


def _deep_getsizeof(obj, seen=None):
    """Recursively calculate size of objects including contents."""
    size = sys.getsizeof(obj)
    if seen is None:
        seen = set()
    obj_id = id(obj)
    if obj_id in seen:
        return 0
    seen.add(obj_id)
    if isinstance(obj, dict):
        size += sum([_deep_getsizeof(k, seen) + _deep_getsizeof(v, seen) for k, v in obj.items()])
    elif isinstance(obj, (list, tuple, set, frozenset)):
        size += sum([_deep_getsizeof(i, seen) for i in obj])
    return size


@dataclass
class Genotype:
    ref_visited: np.ndarray  # bool array indexed by edge index
    ref_counts_multi: dict[int, int]  # edge_index -> count for edges visited >= 2 times
    alt_counts: dict[tuple[str, str], int]
    linear_coverage: list[tuple[int, int]]
    exclude_terminus: bool
    num_edges: int  # total number of edges in graph
    missing_variants: set[tuple[str, str]] | None = None
    coverage_by_chrom: dict[str, list[tuple[int, int]]] = field(default_factory=dict)

    @classmethod
    def genotype(cls, graph: "PangenomeGraph", walk: list[str], exclude_terminus: bool) -> "Genotype":
        """
        Computes the number of time that a walk visits each variant edge.
        :param graph: PangenomeGraph object
        :param walk: list of nodes
        :return: Genotype object
        """

        # Append start and end nodes to walk
        input_walk = walk
        walk = graph.terminal_walk(walk)

        num_edges = graph.number_of_edges() // 2 + 1  # +1 for safety
        ref_visited = np.zeros(num_edges, dtype=np.bool_)
        ref_counts_multi = {}  # only for edges visited >= 2 times
        count_alt = {}
        min_pos = inf
        max_pos = -inf
        for e in zip(walk[:-1], walk[1:]):
            if not graph.has_edge(*e):
                msg = f"Specified list contains edge {e} which is not present in the graph"
                graph.logger.error(msg)
                raise ValueError(msg)

            if not graph.is_terminal(e[0]):
                min_pos = min(min_pos, graph.position(e[0]))
                max_pos = max(max_pos, graph.right_position(e[0]))

            if graph.edges[e]['is_in_tree']:
                edge_idx = int(graph.edges[e]['index'])
                if ref_visited[edge_idx]:
                    # Already visited, increment multi count
                    ref_counts_multi[edge_idx] = ref_counts_multi.get(edge_idx, 1) + 1
                else:
                    ref_visited[edge_idx] = True
            else:
                count_alt[e] = count_alt.get(e, 0) + 1

        coverage = {}
        run_chrom, low, high = None, inf, -inf
        for node in input_walk:
            chrom = graph.nodes[node].get('chromosome', 'chr0')
            if chrom != run_chrom and run_chrom is not None:
                coverage.setdefault(run_chrom, []).append((low, high))
                low, high = inf, -inf
            run_chrom = chrom
            low = min(low, graph.position(node))
            high = max(high, graph.right_position(node))
        if run_chrom is not None:
            coverage.setdefault(run_chrom, []).append((low, high))
        return cls(ref_visited, ref_counts_multi, count_alt, [(min_pos, max_pos)],
                   exclude_terminus, num_edges, coverage_by_chrom=coverage)

    def update(self, other: 'Genotype'):
        # For ref_visited bitmap: edges visited in either become multi if visited in both
        both_visited = self.ref_visited & other.ref_visited
        for edge_idx in np.where(both_visited)[0]:
            # This edge was visited in both, so total count >= 2
            count_self = self.ref_counts_multi.get(edge_idx, 1)
            count_other = other.ref_counts_multi.get(edge_idx, 1)
            self.ref_counts_multi[edge_idx] = count_self + count_other
        # Merge ref_visited
        self.ref_visited = self.ref_visited | other.ref_visited
        # Merge ref_counts_multi for edges only in other
        for edge_idx, count in other.ref_counts_multi.items():
            if edge_idx not in self.ref_counts_multi:
                self.ref_counts_multi[edge_idx] = count
        for key, val in other.alt_counts.items():
            self.alt_counts[key] = self.alt_counts.get(key, 0) + val
        self.linear_coverage += other.linear_coverage
        for chrom, intervals in other.coverage_by_chrom.items():
            self.coverage_by_chrom.setdefault(chrom, []).extend(intervals)

    def memory_breakdown(self) -> dict[str, int]:
        """Return memory usage in bytes for each component."""
        return {
            'ref_visited': self.ref_visited.nbytes,
            'ref_counts_multi': _deep_getsizeof(self.ref_counts_multi),
            'alt_counts': _deep_getsizeof(self.alt_counts),
            'linear_coverage': _deep_getsizeof(self.linear_coverage),
            'missing_variants': _deep_getsizeof(self.missing_variants) if self.missing_variants else 0,
            'ref_visited_count': int(np.sum(self.ref_visited)),
            'ref_counts_multi_len': len(self.ref_counts_multi),
            'alt_counts_len': len(self.alt_counts),
        }

    def compute_missing_variants(self,
                             graph: "PangenomeGraph"):
        """
        Computes variant edges that are missing from a haplotype.
        :param graph: PangenomeGraph object"""

        if len(graph.chromosome_termini) > 1:
            self._compute_missing_by_chromosome(graph)
            return

        # order walks and variants by position
        source_positions = np.sort([x[1] for x in self.linear_coverage] + [graph.position('+_terminus_+')])
        sink_positions = np.sort([x[0] for x in self.linear_coverage] + [graph.right_position('-_terminus_+')])
        sorted_variant_edges = graph.sorted_variant_edges(exclude_terminus=self.exclude_terminus)
        sorted_variant_positions = [min(*graph.position(e)) for e in sorted_variant_edges]

        result = []
        for source, sink in zip(source_positions, sink_positions):
            # indices of first and last variant edges u,v s.t. position of u in between source and sink
            first = np.searchsorted(sorted_variant_positions, source + 1, side='left')
            last = np.searchsorted(sorted_variant_positions, sink - 1, side='right')
            for i in range(first, last):
                if max(*graph.right_position(sorted_variant_edges[i])) <= sink:
                    result.append(sorted_variant_edges[i])

        self.missing_variants = set(result)

    def _compute_missing_by_chromosome(self, graph):
        gaps_by_chrom = {}
        for chrom, (start, end) in graph.chromosome_termini.items():
            cursor = graph.position(start + '_+')
            gaps = []
            for low, high in sorted(self.coverage_by_chrom.get(chrom, [])):
                if low > cursor:
                    gaps.append((cursor, low))
                cursor = max(cursor, high)
            last = graph.right_position(end + '_+')
            if cursor < last:
                gaps.append((cursor, last))
            gaps_by_chrom[chrom] = ([left for left, _ in gaps], gaps)
        self.missing_variants = set()
        for edge in graph.sorted_variant_edges(exclude_terminus=self.exclude_terminus):
            u, v = edge
            chrom = graph.chromosome(u)
            if graph.is_trans_chromosomal(edge):
                low = high = graph.position(u)
            else:
                low = min(graph.position(edge)) + 1
                high = max(graph.right_position(edge))
            starts, gaps = gaps_by_chrom[chrom]
            index = bisect_left(starts, low) - 1
            if not self.coverage_by_chrom.get(chrom) or (index >= 0 and high <= gaps[index][1]):
                self.missing_variants.add(edge)

    def get_ref_count(self, edge_idx: int | list[int]) -> int:
        """Get the count of one tree edge or a set of source continuations."""
        if isinstance(edge_idx, list):
            return sum(self.get_ref_count(index) for index in edge_idx)
        if edge_idx in self.ref_counts_multi:
            return self.ref_counts_multi[edge_idx]
        elif self.ref_visited[edge_idx]:
            return 1
        else:
            return 0

    def variant_record(self, variant_edge: tuple[str, str], reference_edge: tuple[str, str] | None, ref_edge_idx: int | list[int]) -> tuple[int|None, int, int|None]:
        """
        For variant edge e, returns (GT, CR, CA) where GT is the genotype,
        CR is the reference allele count, and CA is the alternative allele count.
        """
        cr = self.get_ref_count(ref_edge_idx)
        ca = get_from_biedge_dict(self.alt_counts, variant_edge, 0)
        gt = None if (self.missing_variants is not None and variant_edge in self.missing_variants) else int(ca > 0)
        if cr + ca > 0 and gt is None:
            logging.warning(
                f"Variant edge {variant_edge} has allele counts (CR={cr}, CA={ca}) but genotype is missing. "
                "This may occur with vg-generated graphs. Consider using --no-missingness."
            )
            gt = int(ca > 0)
        return gt, cr, ca

    @staticmethod
    def _ref_edge_for_index(graph: "PangenomeGraph", edge_idx: int) -> tuple[str, str] | None:
        for u, v, data in graph.edges(data=True):
            if data['is_in_tree'] and int(data['index']) == edge_idx:
                return u, v
        return None

    @staticmethod
    def _check_extra_ref_edges(
        graph: "PangenomeGraph",
        genotype: "Genotype",
        expected_ref_indices: set[int],
        errors: list[str],
    ) -> None:
        for edge_idx in range(genotype.num_edges):
            if edge_idx in expected_ref_indices or genotype.get_ref_count(edge_idx) == 0:
                continue
            edge = Genotype._ref_edge_for_index(graph, edge_idx)
            edge_label = edge if edge is not None else f"index {edge_idx}"
            errors.append(f"Genotype has ref edge {edge_label} not traversed by any walk")

    @staticmethod
    def verify_genotype_matches_walk(
        graph: "PangenomeGraph",
        walk: list[str],
        genotype: "Genotype"
    ) -> tuple[bool, list[str]]:
        """
        Verify that a genotype correctly represents the edges traversed by a walk.
        
        For each edge in the walk:
        - If it's a tree edge (ref), it should be counted in ref_counts
        - If it's a variant edge (alt), it should be counted in alt_counts
        
        Args:
            graph: PangenomeGraph object
            walk: list of nodes (without terminal nodes)
            genotype: Genotype object to verify
            
        Returns:
            Tuple of (is_valid, list of error messages)
        """
        errors = []
        
        # Build expected counts by traversing the walk
        expected_ref = {}
        expected_alt = {}
        
        # Add terminal nodes to walk (same logic as in genotype() method)
        full_walk = graph.terminal_walk(walk)
        
        for e in zip(full_walk[:-1], full_walk[1:]):
            if not graph.has_edge(*e):
                errors.append(f"Walk contains edge {e} which is not in the graph")
                continue
                
            if graph.edges[e]['is_in_tree']:
                expected_ref[e] = expected_ref.get(e, 0) + 1
            else:
                expected_alt[e] = expected_alt.get(e, 0) + 1
        
        # Compare expected vs actual ref counts
        expected_ref_indices = set()
        for edge, expected_count in expected_ref.items():
            edge_idx = int(graph.edges[edge]['index'])
            expected_ref_indices.add(edge_idx)
            actual_count = genotype.get_ref_count(edge_idx)
            if actual_count < expected_count:
                errors.append(
                    f"Ref edge {edge}: expected count >= {expected_count}, got {actual_count}"
                )
        Genotype._check_extra_ref_edges(graph, genotype, expected_ref_indices, errors)
        
        # Compare expected vs actual alt counts
        for edge, expected_count in expected_alt.items():
            actual_count = get_from_biedge_dict(genotype.alt_counts, edge, 0)
            if actual_count < expected_count:
                errors.append(
                    f"Alt edge {edge}: expected count >= {expected_count}, got {actual_count}"
                )
        
        return len(errors) == 0, errors

    @staticmethod
    def verify_genotype_matches_walks(
        graph: "PangenomeGraph",
        walks: list[list[str]],
        genotype: "Genotype"
    ) -> tuple[bool, list[str]]:
        """
        Verify that a genotype correctly represents the edges traversed by multiple walks
        (which may be aggregated together for a single haplotype).
        
        Args:
            graph: PangenomeGraph object
            walks: list of walks (each walk is a list of nodes without terminal nodes)
            genotype: Genotype object to verify
            
        Returns:
            Tuple of (is_valid, list of error messages)
        """
        errors = []
        
        # Build expected counts by traversing all walks
        expected_ref = {}
        expected_alt = {}
        
        for walk in walks:
            # Add terminal nodes to walk
            full_walk = graph.terminal_walk(walk)
            
            for e in zip(full_walk[:-1], full_walk[1:]):
                if not graph.has_edge(*e):
                    errors.append(f"Walk contains edge {e} which is not in the graph")
                    continue
                    
                if graph.edges[e]['is_in_tree']:
                    expected_ref[e] = expected_ref.get(e, 0) + 1
                else:
                    expected_alt[e] = expected_alt.get(e, 0) + 1
        
        # Compare expected vs actual ref counts
        expected_ref_indices = set()
        for edge, expected_count in expected_ref.items():
            edge_idx = int(graph.edges[edge]['index'])
            expected_ref_indices.add(edge_idx)
            actual_count = genotype.get_ref_count(edge_idx)
            if actual_count != expected_count:
                errors.append(
                    f"Ref edge {edge}: expected count {expected_count}, got {actual_count}"
                )
        Genotype._check_extra_ref_edges(graph, genotype, expected_ref_indices, errors)
        
        # Compare expected vs actual alt counts
        for edge, expected_count in expected_alt.items():
            actual_count = get_from_biedge_dict(genotype.alt_counts, edge, 0)
            if actual_count != expected_count:
                errors.append(
                    f"Alt edge {edge}: expected count {expected_count}, got {actual_count}"
                )
        
        for edge in genotype.alt_counts:
            if edge not in expected_alt and edge_complement(edge) not in expected_alt:
                errors.append(f"Genotype has alt edge {edge} not traversed by any walk")
        
        return len(errors) == 0, errors

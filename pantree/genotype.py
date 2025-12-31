from __future__ import annotations
from dataclasses import dataclass
from math import inf
import logging
import numpy as np
from .utils import edge_complement, get_from_biedge_dict

@dataclass
class Genotype:
    ref_counts: dict[str, int]
    alt_counts: dict[str, int]
    linear_coverage: list[tuple[int, int]]
    exclude_terminus: bool
    missing_variants: set[tuple[str, str]] | None = None

    @classmethod
    def genotype(cls, graph: "PangenomeGraph", walk: list[str], exclude_terminus: bool) -> "Genotype":
        """
        Computes the number of time that a walk visits each variant edge.
        :param graph: PangenomeGraph object
        :param walk: list of nodes
        :return: Genotype object
        """

        # Append start and end nodes to walk
        start = [graph.termini[0] + '_+' if graph.direction(walk[0]) == 1 else graph.termini[1] + '_-']
        end = [graph.termini[1] + '_+' if graph.direction(walk[-1]) == 1 else graph.termini[0] + '_-']
        walk = start + walk + end

        count_ref = {}
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
                count_ref[e] = count_ref.get(e, 0) + 1
            else:
                count_alt[e] = count_alt.get(e, 0) + 1

        return cls(count_ref, count_alt, [(min_pos, max_pos)], exclude_terminus)

    def update(self, other: 'Genotype'):
        for key, val in other.ref_counts.items():
            self.ref_counts[key] = self.ref_counts.get(key, 0) + val
        for key, val in other.alt_counts.items():
            self.alt_counts[key] = self.alt_counts.get(key, 0) + val
        self.linear_coverage += other.linear_coverage

    def compute_missing_variants(self,
                             graph: "PangenomeGraph"):
        """
        Computes variant edges that are missing from a haplotype.
        :param graph: PangenomeGraph object"""

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

    def variant_record(self, variant_edge: tuple[str, str], reference_edge: tuple[str, str]) -> tuple[int|None, int, int|None]:
        """
        For variant edge e, returns (GT, CR, CA) where GT is the genotype,
        CR is the reference allele count, and CA is the alternative allele count.
        """
        cr = get_from_biedge_dict(self.ref_counts, reference_edge, 0)
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
        start = [graph.termini[0] + '_+' if graph.direction(walk[0]) == 1 else graph.termini[1] + '_-']
        end = [graph.termini[1] + '_+' if graph.direction(walk[-1]) == 1 else graph.termini[0] + '_-']
        full_walk = start + walk + end
        
        for e in zip(full_walk[:-1], full_walk[1:]):
            if not graph.has_edge(*e):
                errors.append(f"Walk contains edge {e} which is not in the graph")
                continue
                
            if graph.edges[e]['is_in_tree']:
                expected_ref[e] = expected_ref.get(e, 0) + 1
            else:
                expected_alt[e] = expected_alt.get(e, 0) + 1
        
        # Compare expected vs actual ref counts
        for edge, expected_count in expected_ref.items():
            actual_count = get_from_biedge_dict(genotype.ref_counts, edge, 0)
            if actual_count < expected_count:
                errors.append(
                    f"Ref edge {edge}: expected count >= {expected_count}, got {actual_count}"
                )
        
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
            start = [graph.termini[0] + '_+' if graph.direction(walk[0]) == 1 else graph.termini[1] + '_-']
            end = [graph.termini[1] + '_+' if graph.direction(walk[-1]) == 1 else graph.termini[0] + '_-']
            full_walk = start + walk + end
            
            for e in zip(full_walk[:-1], full_walk[1:]):
                if not graph.has_edge(*e):
                    errors.append(f"Walk contains edge {e} which is not in the graph")
                    continue
                    
                if graph.edges[e]['is_in_tree']:
                    expected_ref[e] = expected_ref.get(e, 0) + 1
                else:
                    expected_alt[e] = expected_alt.get(e, 0) + 1
        
        # Compare expected vs actual ref counts
        for edge, expected_count in expected_ref.items():
            actual_count = get_from_biedge_dict(genotype.ref_counts, edge, 0)
            if actual_count != expected_count:
                errors.append(
                    f"Ref edge {edge}: expected count {expected_count}, got {actual_count}"
                )
        
        # Compare expected vs actual alt counts
        for edge, expected_count in expected_alt.items():
            actual_count = get_from_biedge_dict(genotype.alt_counts, edge, 0)
            if actual_count != expected_count:
                errors.append(
                    f"Alt edge {edge}: expected count {expected_count}, got {actual_count}"
                )
        
        # Check for extra edges in genotype that weren't in walks
        for edge in genotype.ref_counts:
            if edge not in expected_ref and edge_complement(edge) not in expected_ref:
                errors.append(f"Genotype has ref edge {edge} not traversed by any walk")
        
        for edge in genotype.alt_counts:
            if edge not in expected_alt and edge_complement(edge) not in expected_alt:
                errors.append(f"Genotype has alt edge {edge} not traversed by any walk")
        
        return len(errors) == 0, errors

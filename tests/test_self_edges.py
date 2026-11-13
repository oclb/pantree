"""Self-complementary hairpins and ordinary loops must survive every graph seam."""

from copy import deepcopy
from collections import Counter
import os
from pathlib import Path
import subprocess
import sys
from unittest.mock import patch

import networkx as nx
import pytest

from pantree.graph import PangenomeGraph
from pantree.genotype import Genotype
from pantree.gfa import read_gfa_line_by_line
from pantree.utils import edge_complement, node_complement, sequence_complement
from pantree.write_gfa import write_GFA


DATA = Path(__file__).parent / 'data'


def assert_bidirected(graph):
    for node, data in graph.nodes(data=True):
        complement = graph.nodes[node_complement(node)]
        assert set(data) >= set(PangenomeGraph().node_attribute_names)
        assert complement['sequence'] == sequence_complement(data['sequence'])
        assert complement['original_ids'] == list(reversed(data['original_ids']))
        assert complement['direction'] == -data['direction']
    for edge in graph.edges:
        assert graph.has_edge(*edge_complement(edge))


@pytest.mark.parametrize('edge', [('a_+', 'a_-'), ('a_-', 'a_+'),
                                  ('a_+', 'a_+'), ('a_-', 'a_-')])
def test_insert_and_count_self_edges(edge):
    graph = PangenomeGraph()
    graph.add_binode('a', 'ACG')
    graph.add_biedge(*edge, weight=3)
    complement = edge_complement(edge)
    assert set(graph.edges) == {edge, complement}
    assert graph.edges[edge]['is_representative'] is True
    assert graph.sorted_biedge_representatives == [(*edge, graph.edges[edge])]
    assert graph.representative_edge(complement) == edge
    graph.compute_edge_weights([[*edge], [*edge], [*complement]])
    assert graph.edges[edge]['weight'] == graph.edges[complement]['weight'] == 6
    assert graph.edges[edge]['index'] == graph.edges[complement]['index'] == 0
    assert type(graph.edges[edge]['index']) is int
    assert graph.number_of_biedges == 1
    with pytest.raises(ValueError, match='duplicate biedge'):
        graph.add_biedge(*complement)


def test_self_edge_genotype_bitmap_bounds():
    graph = PangenomeGraph.from_gfa(str(DATA / 'self_edge_vcf.gfa'), ref_name='REF', return_walks=True)
    genotypes = graph.genotypes_from_gfa(str(DATA / 'self_edge_vcf.gfa'))
    assert genotypes['ALT'][0].num_edges > max(data['index'] for *_, data in graph.edges(data=True))
    assert {('1_+', '1_-'), ('1_-', '1_+'), ('2_+', '2_-'), ('1_+', '1_+')} <= graph.variant_edges
    for genotype in (genotypes['REF'][0], genotypes['ALT'][0]):
        assert genotype.num_edges == graph.biedge_index_capacity == 8
    # The forward REF exit uses the highest-index tree edge, formerly out of bounds.
    assert genotypes['REF'][0].get_ref_count(7) == 1
    assert genotypes['ALT'][0].get_ref_count(7) == 0
    representatives = graph.sorted_biedge_representatives
    assert len(representatives) == graph.number_of_biedges == 8
    assert {data['index'] for *_, data in representatives} == set(range(8))
    for edge in [('1_+', '1_-'), ('1_-', '1_+'), ('2_+', '2_-'), ('1_+', '1_+')]:
        reference = graph.reference_tree_edge(edge)
        index = graph.edges[reference]['index']
        assert genotypes['REF'][0].variant_record(edge, reference, index) == (0, 1, 0)
        assert genotypes['ALT'][0].variant_record(edge, reference, index) == (1, 1, 1)


def test_sparse_indices_reinitialize_capacity_and_allocate_unique_index():
    graph = PangenomeGraph.from_gfa(str(DATA / 'self_edge_vcf.gfa'), ref_name='REF')
    # A copied/filtered graph retains its original assigned index space.
    graph.remove_edge('1_+', '1_-')
    wrapped = PangenomeGraph(nx.DiGraph(graph), reference_tree=graph.reference_tree,
                             reference_path=graph.reference_path)
    assert wrapped.number_of_biedges == 7
    assert wrapped.biedge_index_capacity == 8
    genotype = wrapped.genotype(['1_+', '2_+'])
    assert genotype.num_edges == 8
    assert genotype.get_ref_count(7) == 1
    wrapped.add_biedge('1_+', '1_-')
    assert wrapped.edges['1_+', '1_-']['index'] == 8
    assert wrapped.biedge_index_capacity == 9
    assert len({data['index'] for *_, data in wrapped.sorted_biedge_representatives}) == 8


@pytest.mark.parametrize('operation', ['copy', 'view', 'subgraph', 'edge_subgraph'])
@pytest.mark.parametrize('sparse', [False, True])
def test_networkx_copies_and_views_preserve_genotype_index_capacity(operation, sparse):
    graph = PangenomeGraph.from_gfa(str(DATA / 'no_variants.gfa'), ref_name='ref')
    if sparse:
        for _, _, data in graph.edges(data=True):
            data['index'] = 100 + 7 * data['index']
        graph = PangenomeGraph(nx.DiGraph(graph))
    walk = ['1_+', '2_+', '3_+']
    expected = graph.genotype(walk)
    if operation == 'copy':
        copied = graph.copy()
    elif operation == 'view':
        copied = graph.copy(as_view=True)
    elif operation == 'subgraph':
        copied = graph.subgraph(list(graph))
    else:
        copied = graph.edge_subgraph(list(graph.edges))
    assert copied.biedge_index_capacity == graph.biedge_index_capacity
    # Repeated haplotypes must not rescan the graph's edge indices.
    with patch.object(type(copied.edges), '__call__', side_effect=AssertionError('edge scan')):
        for _ in range(3):
            genotype = copied.genotype(walk)
            assert genotype.num_edges == expected.num_edges
            assert genotype.ref_visited.tolist() == expected.ref_visited.tolist()
            assert genotype.ref_counts_multi == expected.ref_counts_multi
            assert genotype.alt_counts == expected.alt_counts


def test_index_capacity_tracks_source_growth_in_views_but_not_independent_copies():
    graph = PangenomeGraph.from_gfa(str(DATA / 'no_variants.gfa'), ref_name='ref')
    view = graph.copy(as_view=True)
    copied = graph.copy()
    original_capacity = graph.biedge_index_capacity
    graph.add_biedge('1_+', '3_+')
    assert view.biedge_index_capacity == graph.biedge_index_capacity == original_capacity + 1
    assert copied.biedge_index_capacity == original_capacity
    copied.add_biedge('1_+', '3_+')
    assert copied.edges['1_+', '3_+']['index'] == original_capacity
    assert copied.biedge_index_capacity == original_capacity + 1


def test_index_capacity_recovers_after_clearing_graph_metadata():
    graph = PangenomeGraph.from_gfa(str(DATA / 'no_variants.gfa'), ref_name='ref')
    expected = graph.genotype(['1_+', '2_+', '3_+'])
    graph.graph.clear()
    assert graph.biedge_index_capacity == expected.num_edges
    with patch.object(type(graph.edges), '__call__', side_effect=AssertionError('edge scan')):
        assert graph.genotype(['1_+', '2_+', '3_+']).ref_visited.tolist() == expected.ref_visited.tolist()


@pytest.mark.parametrize('operation', ['copy', 'view_copy', 'deepcopy', 'constructor'])
@pytest.mark.parametrize('sparse', [False, True])
def test_metadata_free_copy_recovers_from_populated_indices(operation, sparse):
    graph = PangenomeGraph.from_gfa(str(DATA / 'no_variants.gfa'), ref_name='ref')
    if sparse:
        for _, _, data in graph.edges(data=True):
            data['index'] = 100 + 7 * data['index']
        graph = PangenomeGraph(nx.DiGraph(graph))
    walk = ['1_+', '2_+', '3_+']
    expected = graph.genotype(walk)
    graph.graph.clear()
    if operation == 'copy':
        copied = graph.copy()
    elif operation == 'view_copy':
        copied = graph.subgraph([node for node in graph if node.endswith('_+')]).copy()
    elif operation == 'deepcopy':
        copied = deepcopy(graph)
    else:
        copied = PangenomeGraph(nx.DiGraph(graph))
    # The source getter has not run since metadata was cleared.
    assert copied.biedge_index_capacity == expected.num_edges
    with patch.object(type(copied.edges), '__call__', side_effect=AssertionError('edge scan')):
        for _ in range(3):
            assert copied.genotype(walk).ref_visited.tolist() == expected.ref_visited.tolist()
    copied.add_biedge('1_+', '3_+')
    assert copied.edges['1_+', '3_+']['index'] == expected.num_edges


@pytest.mark.parametrize('operation', ['nodes', 'edges', 'chained'])
@pytest.mark.parametrize('sparse', [False, True])
def test_metadata_free_filtered_view_recovers_the_owning_index_space(operation, sparse):
    graph = PangenomeGraph.from_gfa(str(DATA / 'no_variants.gfa'), ref_name='ref')
    if sparse:
        for _, _, data in graph.edges(data=True):
            data['index'] = 100 + 7 * data['index']
        graph = PangenomeGraph(nx.DiGraph(graph))
    walk = ['1_+', '2_+', '3_+']
    expected = graph.genotype(walk)
    if operation == 'nodes':
        view = graph.subgraph(['1_+', '2_+'])
    elif operation == 'edges':
        view = graph.edge_subgraph([('1_+', '2_+')])
    else:
        view = graph.subgraph(['1_+', '2_+']).edge_subgraph([('1_+', '2_+')]).copy(as_view=True)
    occupied_indices = {data['index'] for _, _, data in graph.edges(data=True)}
    graph.graph.clear()
    # Access the filtered view first: its single edge does not define the source's capacity.
    assert view.biedge_index_capacity == expected.num_edges
    assert graph.genotype(walk).ref_visited.tolist() == expected.ref_visited.tolist()
    graph.add_biedge('1_+', '3_+')
    assert graph.edges['1_+', '3_+']['index'] not in occupied_indices
    assert graph.edges['1_+', '3_+']['index'] == expected.num_edges
    with patch.object(type(graph.edges), '__call__', side_effect=AssertionError('edge scan')):
        for _ in range(3):
            assert view.biedge_index_capacity == expected.num_edges + 1
            assert graph.genotype(walk).ref_visited.tolist() == expected.ref_visited.tolist() + [False]


@pytest.mark.parametrize(('fixture', 'hidden', 'walk'), [
    ('no_variants.gfa', ('1_+', '2_+'), ['3_-', '2_-', '1_-']),
    ('self_edge_repeated.gfa', ('a_+', 'a_+'), ['z_-', 'a_-', 'a_-', 'a_-', 's_-']),
])
def test_verifiers_accept_filtered_views_without_the_representative(fixture, hidden, walk):
    graph = PangenomeGraph.from_gfa(str(DATA / fixture), ref_name='ref')
    view = graph.edge_subgraph([edge for edge in graph.edges if edge != hidden])
    assert not view.has_edge(*hidden)
    genotype = view.genotype(walk)
    assert Genotype.verify_genotype_matches_walk(view, walk, genotype) == (True, [])
    assert Genotype.verify_genotype_matches_walks(view, [walk], genotype) == (True, [])
    combined = deepcopy(genotype)
    combined.update(genotype)
    assert Genotype.verify_genotype_matches_walks(view, [walk, walk], combined) == (True, [])

    # The helpers still reject an undercount rather than simply accepting filtered views.
    bad = deepcopy(genotype)
    if hidden in bad.alt_counts:
        bad.alt_counts[hidden] = 0
    else:
        index = view.edges[walk[0], walk[1]]['index']
        bad.ref_visited[index] = False
    assert Genotype.verify_genotype_matches_walk(view, walk, bad)[0] is False
    assert Genotype.verify_genotype_matches_walks(view, [walk], bad)[0] is False


@pytest.mark.parametrize('reverse_loop', [False, True])
def test_repeated_self_edges_exact_genotypes_and_emitted_alleles(tmp_path, reverse_loop):
    source = tmp_path / 'repeated.gfa'
    content = (DATA / 'self_edge_repeated.gfa').read_text()
    if reverse_loop:
        content = content.replace('L\ta\t+\ta\t+\t0M', 'L\ta\t-\ta\t-\t0M')
    source.write_text(content)
    graph = PangenomeGraph.from_gfa(str(source), ref_name='ref', return_walks=True)
    counts = {('a_+', 'a_-'): 2, ('a_-', 'a_+'): 2,
              (('a_-', 'a_-') if reverse_loop else ('a_+', 'a_+')): 3}
    assert graph.variant_edges == set(counts)
    assert set(graph.sorted_variant_edges()) == set(counts)
    genotypes = graph.genotypes_from_gfa(str(source))
    reference, alternate = genotypes['ref'][0], genotypes['ALT'][0]
    assert alternate.alt_counts == counts
    for edge, count in counts.items():
        ref_edge = graph.reference_tree_edge(edge)
        index = graph.edges[ref_edge]['index']
        assert reference.variant_record(edge, ref_edge, index) == (0, 1, 0)
        assert alternate.variant_record(edge, ref_edge, index) == (1, 1, count)
        assert alternate.variant_record(edge_complement(edge), ref_edge, index) == (1, 1, count)
        assert graph.edges[edge]['weight'] == count
        assert graph.edges[edge_complement(edge)]['weight'] == count
    assert graph.ref_alt_alleles(('a_+', 'a_+'))[:2] == ('', 'ACG')
    for edge in [('a_+', 'a_-'), ('a_-', 'a_+')]:
        # A hairpin changes strand at a segment boundary, with no intervening DNA.
        assert graph.ref_alt_alleles(edge)[:2] == ('', '')
    valid, errors = alternate.verify_genotype_matches_walks(graph, [graph.walks['ALT#0#chr1']], alternate)
    assert valid, errors
    vcf = tmp_path / 'self-edges.vcf'
    graph.write_vcf(str(source), str(vcf), 'chr1')
    records = {fields[2]: fields for row in vcf.read_text().splitlines()
               if not row.startswith('#') and (fields := row.split('\t'))}
    assert set(records) == {'>a<a', '<a>a', '>a>a'}
    assert records['>a>a'][1:5] == ['2', '>a>a', 'A', 'AACG']
    assert records['>a<a'][1:5] == ['5', '>a<a', 'G', 'G']
    assert records['<a>a'][1:5] == ['2', '<a>a', '.', '.']
    for variant_id, count, variant_type in [('>a>a', 3, 'DUP'), ('>a<a', 2, 'INV'), ('<a>a', 2, 'INV')]:
        row = records[variant_id]
        info = dict(item.split('=', 1) for item in row[7].split(';'))
        assert info['VT'] == variant_type
        assert (info['RC'], info['AC']) == ('2', str(count))
        assert row[8:] == ['GT:CR:CA', '0:1:0', f'1:1:{count}']
    alternate.update(graph.genotype(graph.walks['ALT#0#chr1']))
    for edge, count in counts.items():
        ref_edge = graph.reference_tree_edge(edge)
        assert alternate.variant_record(edge, ref_edge, graph.edges[ref_edge]['index']) == (1, 2, 2 * count)
    valid, errors = alternate.verify_genotype_matches_walks(graph, [graph.walks['ALT#0#chr1']] * 2, alternate)
    assert valid, errors


@pytest.mark.parametrize(('fixture', 'ref_name'), [
    ('self_edge_simplify.gfa', 'REF'), ('disconnected_self_edge.gfa', 'ref'),
])
@pytest.mark.parametrize('threshold', [0, 5, 50])
def test_self_edge_simplify_attributes(fixture, ref_name, threshold):
    graph = PangenomeGraph.from_gfa(str(DATA / fixture), ref_name=ref_name, return_walks=True)
    original = deepcopy((dict(graph.nodes(data=True)), dict(((u, v), data) for u, v, data in graph.edges(data=True)),
                         graph.reference_path, graph.walks))
    if fixture == 'disconnected_self_edge.gfa':
        # The corrected change-2 direction assignment is necessary to expose the P2.
        assert graph.direction('c_+') == graph.direction('f_-') == 1
    simplified = graph.simplify_subgraph(minimum_allele_length=threshold)
    assert_bidirected(simplified)
    assert original == (dict(graph.nodes(data=True)), dict(((u, v), data) for u, v, data in graph.edges(data=True)),
                        graph.reference_path, graph.walks)


@pytest.mark.parametrize('orientation', ['+', '-'])
def test_contraction_retains_hairpin_with_deleted_rewiring_target(orientation):
    graph = PangenomeGraph()
    graph.add_binode('b', 'T')
    graph.add_binode('a', 'ACG')
    node = f'a_{orientation}'
    graph.add_biedge('b_+', node)
    graph.add_biedge(node, node_complement(node))
    before = deepcopy(nx.node_link_data(graph, edges='edges'))
    # Contracting b -> a would reattach a -> complement(a) to a deleted endpoint.
    PangenomeGraph.contract_paths(graph, {'b_+', 'b_-'})
    assert nx.node_link_data(graph, edges='edges') == before
    assert_bidirected(graph)


@pytest.mark.parametrize('orientation', ['+', '-'])
def test_contraction_retains_ordinary_self_loop(orientation):
    graph = PangenomeGraph()
    graph.add_binode('a', 'ACG')
    node = f'a_{orientation}'
    graph.add_biedge(node, node)
    before = deepcopy(nx.node_link_data(graph, edges='edges'))
    PangenomeGraph.contract_paths(graph, set(graph))
    assert nx.node_link_data(graph, edges='edges') == before
    assert_bidirected(graph)


def test_contraction_closes_cycle_as_paired_loops_without_merging_itself():
    graph = PangenomeGraph()
    graph.add_binode('b', 'T')
    graph.add_binode('a', 'ACG')
    graph.add_biedge('b_+', 'a_+')
    graph.add_biedge('a_+', 'b_+', weight=4)
    expected_data = dict(graph.edges['a_+', 'b_+'])
    expected_complement_data = dict(graph.edges['b_-', 'a_-'])
    PangenomeGraph.contract_paths(graph, {'b_+', 'b_-'})
    assert set(graph) == {'b_+', 'b_-'}
    assert set(graph.edges) == {('b_+', 'b_+'), ('b_-', 'b_-')}
    assert graph.nodes['b_+']['sequence'] == 'TACG'
    assert graph.nodes['b_+']['original_ids'] == ['b', 'a']
    assert graph.edges['b_+', 'b_+'] == expected_data
    assert graph.edges['b_-', 'b_-'] == expected_complement_data
    assert_bidirected(graph)


def reload_gfa(path):
    """The simplified-GFA contract has S/L records, without a reference W line."""
    graph = PangenomeGraph()
    for record in read_gfa_line_by_line(str(path), ['S']):
        graph.add_binode(record.node_id, record.sequence)
        graph.nodes[f'{record.node_id}_+']['original_ids'] = record.original_ids
        graph.nodes[f'{record.node_id}_-']['original_ids'] = list(reversed(record.original_ids))
    for record in read_gfa_line_by_line(str(path), ['L']):
        graph.add_biedge(record.u, record.v)
    return graph


def sequence_topology(graph):
    """Compare the complete oriented sequence/provenance graph despite segment renaming."""
    result = nx.DiGraph()
    for node, data in graph.nodes(data=True):
        result.add_node(node, dna=data['sequence'], ids=data['original_ids'])
    result.add_edges_from(graph.edges)
    return result


@pytest.mark.parametrize(('fixture', 'ref_name'), [
    ('self_edge_vcf.gfa', 'REF'), ('self_edge_simplify.gfa', 'REF'),
    ('disconnected_self_edge.gfa', 'ref'), ('self_edge_repeated.gfa', 'ref'),
    ('self_edge_repeated.gfa', 'ALT'),
])
@pytest.mark.parametrize('threshold', [0, 1, 3, 4, 5, 50])
def test_self_edge_thresholds_source_dna_and_roundtrip(tmp_path, fixture, ref_name, threshold):
    source = DATA / fixture
    graph = PangenomeGraph.from_gfa(str(source), ref_name=ref_name)
    filtered = graph.delete_small_variants(minimum_allele_length=threshold)
    reference_edges = set(zip(graph.reference_path, graph.reference_path[1:]))
    assert reference_edges <= set(filtered.edges)
    if fixture == 'self_edge_repeated.gfa':
        for edge, length in [(('a_+', 'a_-'), 0), (('a_-', 'a_+'), 0), (('a_+', 'a_+'), 3)]:
            assert filtered.has_edge(*edge) == (ref_name == 'ALT' or threshold <= length)
    simplified = graph.simplify_subgraph(minimum_allele_length=threshold)
    assert_bidirected(simplified)
    # Each original segment occurs once per retained orientation, even after contractions.
    source_ids = Counter(segment for _, data in simplified.nodes(data=True) for segment in data['original_ids'])
    assert all(count == 2 for count in source_ids.values())
    if threshold == 0 or ref_name == 'ALT':
        assert set(source_ids) == {node[:-2] for node in graph if not graph.is_terminal(node)}
    # Every merged segment spells a real oriented path in the source graph.
    for _, data in simplified.nodes(data=True):
        paths = [[]]
        for segment in data['original_ids']:
            paths = [path + [node] for path in paths for node in (f'{segment}_+', f'{segment}_-')
                     if not path or graph.has_edge(path[-1], node)]
        assert any(graph.walk_sequence(path) == data['sequence'] for path in paths)
    if fixture == 'disconnected_self_edge.gfa':
        if threshold == 0:
            assert simplified.has_edge('f_-', 'f_+')
            assert simplified.nodes['f_-']['sequence'] == 'ACG'
            assert simplified.nodes['f_+']['sequence'] == 'CGT'
            assert ''.join(simplified.nodes[n]['sequence'] for n in ['c_+', 'f_-', 'f_+', 'c_-']) == 'TACGCGTA'
        if threshold == 50:
            assert set(source_ids) == {'r'}
            assert any(data['sequence'] == 'AC' and data['original_ids'] == ['r']
                       for _, data in simplified.nodes(data=True))
    output = tmp_path / 'simplified.gfa'
    write_GFA(simplified, output, source_gfa=str(source))
    assert output.read_text().splitlines()[0] == f'H\tVN:Z:1.0\tog:Z:{source}'
    reloaded = reload_gfa(output)
    assert_bidirected(reloaded)
    assert nx.is_isomorphic(sequence_topology(simplified), sequence_topology(reloaded),
                            node_match=lambda left, right: left == right)
    # Self-complementary links are serialized once; ordinary loops have two orientations internally.
    assert len(list(read_gfa_line_by_line(str(output), ['L']))) == len(reloaded.sorted_biedge_representatives)


@pytest.mark.parametrize(('fixture', 'ref_name', 'threshold'), [
    ('disconnected_self_edge.gfa', 'ref', 50), ('disconnected_self_edge.gfa', 'ref', 0),
    ('self_edge_simplify.gfa', 'REF', 5), ('self_edge_simplify.gfa', 'REF', 0),
    ('self_edge_repeated.gfa', 'ALT', 50),
])
def test_self_edge_cli_simplification_reproducible_and_terminates(tmp_path, fixture, ref_name, threshold):
    outputs = []
    for seed in (0, 2, 42):
        output = tmp_path / f'seed-{seed}.gfa'
        completed = subprocess.run(
            [sys.executable, '-m', 'pantree.cli', 'simplify', str(DATA / fixture), str(output),
             '--ref-name', ref_name, '--min-allele-length', str(threshold)],
            env={**os.environ, 'PYTHONHASHSEED': str(seed)},
            capture_output=True, text=True, timeout=30,
        )
        assert completed.returncode == 0, completed.stdout + completed.stderr
        assert_bidirected(reload_gfa(output))
        outputs.append(output.read_bytes())
    assert outputs[0] == outputs[1] == outputs[2]

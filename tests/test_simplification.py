"""Sequence, provenance, and reproducibility contracts for simplification."""

import json
import os
from pathlib import Path
import subprocess
import sys

from Bio.Seq import Seq
import networkx as nx
import pytest

from pantree.dfs import assign_node_directions
from pantree.gfa import read_gfa_line_by_line
from pantree.graph import PangenomeGraph
from pantree.utils import node_complement
from pantree.write_gfa import write_GFA


DATA = Path(__file__).parent / "data"


def assert_complement_pairs(graph):
    for node, data in graph.nodes(data=True):
        complement = graph.nodes[node_complement(node)]
        assert complement['sequence'] == str(Seq(data['sequence']).reverse_complement())
        assert complement['original_ids'] == list(reversed(data['original_ids']))


def emitted_segments(path):
    return [
        (fields[2], json.loads(fields[3][5:]))
        for line in path.read_text().splitlines()
        if (fields := line.split('\t'))[0] == 'S'
    ]


@pytest.mark.parametrize('contract_from', ['forward', 'reverse'])
@pytest.mark.parametrize(('orientations', 'expected_dna'), [
    ('+++', 'ACGTTGAGCTA'),
    ('---', 'CGTTCAATAGC'),
    ('+-+', 'ACGTCAAGCTA'),
    ('-+-', 'CGTTTGATAGC'),
])
def test_contract_paths_oriented_sequence_and_provenance(
    tmp_path, contract_from, orientations, expected_dna
):
    """Force each contraction direction with ordinary, non-palindromic DNA."""
    graph = PangenomeGraph()
    start, end = ('a-start', 'z-end') if contract_from == 'forward' else ('z-start', 'a-end')
    graph.add_binode(start)
    graph.add_binode(end)
    for boundary in (start, end):
        for orientation in '+-':
            graph.nodes[f'{boundary}_{orientation}']['original_ids'] = []
    for segment, sequence in [('m', 'ACG'), ('n', 'TTGA'), ('p', 'GCTA')]:
        graph.add_binode(segment, sequence)
    walk = [f'{segment}_{orientation}' for segment, orientation in zip('mnp', orientations)]
    path = [f'{start}_+', *walk, f'{end}_+']
    for u, v in zip(path, path[1:]):
        graph.add_biedge(u, v)
    endpoints = {f'{boundary}_{orientation}' for boundary in (start, end) for orientation in '+-'}

    PangenomeGraph.contract_paths(graph, endpoints)

    retained = start if contract_from == 'forward' else end
    assert graph.nodes[f'{retained}_+']['sequence'] == expected_dna
    assert graph.nodes[f'{retained}_+']['original_ids'] == ['m', 'n', 'p']
    assert set(graph.nodes) == endpoints
    assert set(graph.edges) == {(f'{start}_+', f'{end}_+'), (f'{end}_-', f'{start}_-')}
    assert_complement_pairs(graph)
    output = tmp_path / 'contracted.gfa'
    write_GFA(graph, output, source_gfa='oriented-source.gfa')
    assert (expected_dna, ['m', 'n', 'p']) in emitted_segments(output)
    assert output.read_text().splitlines()[0] == 'H\tVN:Z:1.0\tog:Z:oriented-source.gfa'


@pytest.mark.parametrize('minimum_allele_length', [0, 5, 50])
def test_nested_simplification_matches_source_walks(tmp_path, minimum_allele_length):
    source = DATA / 'simple_nested.gfa'
    graph = PangenomeGraph.from_gfa(str(source), ref_name='sample1', return_walks=True)
    original = {node: (data['sequence'], list(data['original_ids']))
                for node, data in graph.nodes(data=True)}
    simplified = graph.simplify_subgraph(minimum_allele_length=minimum_allele_length)
    assert_complement_pairs(simplified)
    # Every retained segment spells a contiguous, known forward source walk.
    walks = [[node[:-2] for node in walk] for walk in graph.walks.values()]
    for node, data in simplified.nodes(data=True):
        if not node.endswith('_+'):
            continue
        ids = data['original_ids']
        assert data['sequence'] == ''.join(graph.nodes[f'{segment}_+']['sequence'] for segment in ids)
        assert not ids or any(walk[i:i + len(ids)] == ids
                              for walk in walks for i in range(len(walk)))
    output = tmp_path / 'nested.gfa'
    write_GFA(simplified, output, source_gfa=str(source))
    for sequence, ids in emitted_segments(output):
        assert sequence == ''.join(graph.nodes[f'{segment}_+']['sequence'] for segment in ids)
    if minimum_allele_length == 50:
        assert ('AGGATGGCACAAGATCCGGTAGGTACTGCATCATG',
                ['1', '2', '4', '5', '6', '8', '9', '10']) in emitted_segments(output)
    assert original == {node: (data['sequence'], data['original_ids'])
                        for node, data in graph.nodes(data=True)}


@pytest.mark.parametrize('reverse_insertion', [False, True])
def test_disconnected_directions_visit_seed_and_follow_edges(reverse_insertion):
    graph = PangenomeGraph()
    segments = ['ref', 'a', 'b', 'c']
    for segment in reversed(segments) if reverse_insertion else segments:
        graph.add_binode(segment, 'ACG')
    # The second component requires flipping b; the isolated c is a third component.
    graph.add_biedge('a_+', 'b_-')
    assign_node_directions(graph, ['ref_+'])
    assert graph.direction('ref_+') == 1
    assert graph.direction('a_+') == 1
    assert graph.direction('b_-') == 1
    assert graph.direction('c_+') == 1
    for segment in segments:
        assert graph.direction(f'{segment}_+') == -graph.direction(f'{segment}_-')


@pytest.fixture
def disconnected_gfa(tmp_path):
    """A reference, an inverted non-reference path, and a separate isolated segment."""
    source = tmp_path / 'disconnected.gfa'
    source.write_text(
        'H\tVN:Z:1.1\n'
        'S\tr1\tACG\nS\tr2\tTTGA\n'
        'S\ta\tGCTA\nS\tb\tAACG\nS\tc\tCAT\n'
        'L\tr1\t+\tr2\t+\t0M\nL\ta\t+\tb\t-\t0M\n'
        'W\tref\t0\t0\t0\t7\t>r1>r2\n'
        'W\tsample\t1\t0\t0\t8\t>a<b\n'
    )
    return source


def test_disconnected_simplification_sequence_and_topology(disconnected_gfa, tmp_path):
    graph = PangenomeGraph.from_gfa(str(disconnected_gfa), ref_name='ref')
    assert graph.direction('a_+') == graph.direction('b_-') == 1
    simplified = graph.simplify_subgraph(minimum_allele_length=0)
    assert_complement_pairs(simplified)
    # Three source-to-sink paths survive, with the inverted component spelled correctly.
    plus, minus = '+_terminus_+', '-_terminus_+'
    spelled = {''.join(simplified.nodes[node]['sequence'] for node in path)
               for path in nx.all_simple_paths(simplified, plus, minus)}
    assert spelled == {'ACGTTGA', 'GCTACGTT', 'CAT'}
    output = tmp_path / 'disconnected-output.gfa'
    write_GFA(simplified, output, source_gfa=str(disconnected_gfa))
    assert ('GCTACGTT', ['a', 'b']) in emitted_segments(output)


@pytest.mark.parametrize(('fixture', 'ref_name', 'minimum_allele_length', 'interval', 'method'), [
    ('simple_nested.gfa', 'sample1', 50, None, 'max_weight'),
    ('simple_nested.gfa', 'sample1', 0, None, 'max_weight'),
    ('simple_nested.gfa', 'sample1', 5, (5, 25), 'max_weight'),
    ('simple_nested.gfa', 'sample1', 5, None, 'contiguous'),
    ('disconnected', 'ref', 0, None, 'max_weight'),
])
def test_simplified_gfa_identical_across_hash_seeds(
    tmp_path, disconnected_gfa, fixture, ref_name, minimum_allele_length, interval, method
):
    source = disconnected_gfa if fixture == 'disconnected' else DATA / fixture
    script = (
        'import sys\n'
        'from pantree.simplify_graph import simplify_graph\n'
        f'simplify_graph(sys.argv[1], sys.argv[2], ref_name={ref_name!r}, '
        f'min_allele_length={minimum_allele_length}, dfs_method={method!r}, '
        f'pos_start={interval[0] if interval else None!r}, '
        f'pos_end={interval[1] if interval else None!r})\n'
    )
    outputs = []
    for seed in (0, 2, 42):
        output = tmp_path / f'seed-{seed}.gfa'
        result = subprocess.run(
            [sys.executable, '-c', script, str(source), str(output)],
            env={**os.environ, 'PYTHONHASHSEED': str(seed)},
            capture_output=True, text=True, timeout=30,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        outputs.append(output.read_bytes())
    assert outputs[0] == outputs[1] == outputs[2]


def test_gfa_serialization_independent_of_insertion_order(tmp_path):
    """Stable numbering and complementary-link selection on an already simplified graph."""
    graph = PangenomeGraph()
    for segment, sequence in [('z', 'ACG'), ('a', 'TTGA'), ('m', 'GCTA')]:
        graph.add_binode(segment, sequence)
    graph.add_biedge('z_+', 'a_-')
    graph.add_biedge('a_-', 'm_+')
    reordered = nx.DiGraph()
    reordered.add_nodes_from(reversed(list(graph.nodes(data=True))))
    reordered.add_edges_from(reversed(list(graph.edges(data=True))))
    original_output, reordered_output = tmp_path / 'original.gfa', tmp_path / 'reordered.gfa'
    write_GFA(graph, original_output, source_gfa='input.gfa')
    write_GFA(reordered, reordered_output, source_gfa='input.gfa')
    assert original_output.read_bytes() == reordered_output.read_bytes()
    # Re-reading verifies the serialized orientation of both links, not only their count.
    reread = PangenomeGraph()
    for line in read_gfa_line_by_line(str(original_output), line_types=['S']):
        reread.add_binode(line.node_id, line.sequence)
    for line in read_gfa_line_by_line(str(original_output), line_types=['L']):
        reread.add_biedge(line.u, line.v)
    assert set(reread.edges) == {('s3_+', 's1_-'), ('s1_+', 's3_-'),
                                 ('s1_-', 's2_+'), ('s2_-', 's1_+')}

"""Regression tests for reference walks that revisit segments (issue #22)."""
from pathlib import Path
import networkx as nx
import pytest
from pantree.graph import PangenomeGraph
from pantree.dfs import max_weight_dfs_tree

DATA = Path(__file__).parent / 'data'

@pytest.mark.parametrize('fixture', ['hairpin', 'tandem', 'hairpin_tip', 'acyclic_unreachable'])
def test_reference_tree_spans_positive_nodes_and_commands_complete(fixture, tmp_path):
    path = DATA / f'{fixture}.gfa'
    with pytest.warns(UserWarning):
        graph = PangenomeGraph.from_gfa(str(path), ref_name='REF')
    assert nx.is_arborescence(graph.reference_tree)
    assert set(graph.reference_tree) == {n for n, d in graph.nodes(data='direction') if d == 1}
    filtered = graph.delete_small_variants(minimum_allele_length=50)
    assert all(filtered.has_edge(u, v) for u, v in zip(graph.reference_path, graph.reference_path[1:]))
    simplified = graph.simplify_subgraph(minimum_allele_length=50)
    assert simplified.number_of_nodes() > 0
    output = tmp_path / 'out.vcf'
    graph.write_vcf(str(path), str(output), chr_name='chr1')
    content = output.read_text()
    assert '#CHROM' in content
    if fixture == 'hairpin':
        records = [line.split('\t') for line in content.splitlines() if not line.startswith('#')]
        assert any('VT=INV' in row[7] and row[9:] == ['1:1:1', '0:1:0']
                   for row in records)

@pytest.mark.parametrize('sample', ['HG02257', 'HG02622', 'HG03486', 'HG03492', 'HG03540'])
def test_existing_c4_cyclic_references(sample, tmp_path):
    path = DATA / 'c4a_with_inversion_and_sequences.gfa'
    with pytest.warns(UserWarning, match='Reference walk is cyclic'):
        graph = PangenomeGraph.from_gfa(str(path), ref_name=sample)
    assert nx.is_arborescence(graph.reference_tree)
    graph.simplify_subgraph(minimum_allele_length=50)
    graph.write_vcf(str(path), str(tmp_path / 'out.vcf'), chr_name='chr1')


def test_acyclic_reference_keeps_priority_over_edge_weight():
    graph = nx.DiGraph()
    graph.add_weighted_edges_from([('a', 'b', 1), ('a', 'c', 100), ('b', 'd', 1), ('c', 'd', 100)])
    tree = max_weight_dfs_tree(graph, ['a', 'b', 'd'])
    assert set(tree.edges) == {('a', 'b'), ('b', 'd'), ('a', 'c')}


def test_cyclic_priority_remains_depth_first():
    graph = nx.DiGraph()
    graph.add_weighted_edges_from([('a', 'b', 1), ('b', 'c', 1), ('c', 'b', 1),
                                   ('b', 'd', 1), ('c', 'd', 1), ('d', 'a', 1)])
    with pytest.warns(UserWarning, match='Reference walk is cyclic'):
        tree = max_weight_dfs_tree(graph, ['a', 'b', 'c', 'b', 'd'])
    assert nx.is_arborescence(tree)
    assert ('c', 'd') in tree.edges  # Do not jump out of the active DFS branch.


def test_contraction_skips_removed_complements(tmp_path):
    path = tmp_path / 'contraction.gfa'
    path.write_text((DATA / 'hairpin.gfa').read_text()
                    .replace('L\t3\t+\t2\t-', 'L\t3\t+\t1\t-')
                    .replace('L\t2\t-\t4', 'L\t1\t-\t4')
                    .replace('>1>2>3<2>4', '>1>2>3<1>4'))
    with pytest.warns(UserWarning, match='Reference walk is cyclic'):
        graph = PangenomeGraph.from_gfa(str(path), ref_name='REF')
    simplified = graph.simplify_subgraph(minimum_allele_length=50)
    assert simplified.number_of_nodes() > 0
    assert all('sequence' in data for _, data in simplified.nodes(data=True))

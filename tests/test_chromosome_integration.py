"""Combined PR24 / cyclic-reference and self-link integration contracts."""
from copy import deepcopy
import os
from pathlib import Path
import re
import subprocess
import sys
from unittest.mock import patch

import networkx as nx
import pytest
from click.testing import CliRunner

from pantree.cli import cli
from pantree.dfs import dfs_methods
from pantree.genotype import Genotype
from pantree.gfa import read_gfa_line_by_line
from pantree.graph import PangenomeGraph
from pantree.utils import edge_complement, node_complement

DATA = Path(__file__).parent / 'data'


@pytest.mark.parametrize('method', list(dfs_methods))
@pytest.mark.parametrize('clear_metadata', [False, True])
@pytest.mark.parametrize('operation', ['copy', 'live', 'nodes', 'edges', 'filtered_copy', 'deepcopy', 'import'])
def test_multichromosome_copy_view_terminal_metadata(method, clear_metadata, operation):
    graph = PangenomeGraph.from_gfa(str(DATA / 'multi_chrom.gfa'), ref_name='ref',
                                    dfs_method=dfs_methods[method])
    walk = ['14_-', '12_-', '11_-']
    expected = graph.genotype(walk)
    expected_termini = dict(graph.chromosome_termini)
    kept_edges = [edge for edge in graph.edges if edge != ('11_+', '12_+')]
    if clear_metadata:
        graph.graph.clear()
    if operation == 'copy':
        copied = graph.copy()
    elif operation == 'live':
        copied = graph.copy(as_view=True)
    elif operation == 'nodes':
        copied = graph.subgraph(list(graph))
    elif operation == 'edges':
        copied = graph.edge_subgraph(kept_edges)
    elif operation == 'filtered_copy':
        copied = graph.edge_subgraph(kept_edges).copy()
    elif operation == 'deepcopy':
        copied = deepcopy(graph)
    else:
        copied = PangenomeGraph(nx.DiGraph(graph))
    assert copied.chromosome_termini == expected_termini
    assert copied.terminal_walk(walk) == ['-_terminus_1_-', *walk, '+_terminus_1_-']
    assert all(copied.is_terminal(name + sign) for pair in expected_termini.values()
               for name in pair for sign in ('_+', '_-'))
    assert copied.biedge_index_capacity == expected.num_edges
    with patch.object(type(copied.edges), '__call__', side_effect=AssertionError('edge scan')):
        for _ in range(3):
            genotype = copied.genotype(walk)
            assert genotype.ref_visited.tolist() == expected.ref_visited.tolist()
            assert genotype.alt_counts == expected.alt_counts
            assert genotype.coverage_by_chrom == expected.coverage_by_chrom
            assert Genotype.verify_genotype_matches_walk(copied, walk, genotype) == (True, [])
            assert Genotype.verify_genotype_matches_walks(copied, [walk], genotype) == (True, [])
    bad = deepcopy(genotype)
    index = copied.edges['14_-', '12_-']['index']
    bad.ref_visited[index] = False
    assert Genotype.verify_genotype_matches_walk(copied, walk, bad)[0] is False
    assert Genotype.verify_genotype_matches_walks(copied, [walk], bad)[0] is False


@pytest.mark.parametrize('method', list(dfs_methods))
@pytest.mark.parametrize('multiple_chromosomes', [False, True])
def test_ignored_reference_haplotype_fragments_are_genotyped(tmp_path, method, multiple_chromosomes):
    text = ('S\t1\tA\nS\t2\tC\nS\t3\tG\nS\t4\tT\n'
            'L\t1\t+\t2\t+\t0M\nL\t2\t+\t4\t+\t0M\n'
            'L\t1\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\n'
            'W\tref\t0\tchr1\t0\t3\t>1>2>4\n')
    if multiple_chromosomes:
        text += 'S\t5\tA\nW\tref\t0\tchr2\t0\t1\t>5\n'
    text += ('W\tref\t1\tassembly\t0\t1\t>1\n'
             'W\tref\t1\tassembly\t1\t3\t>3>4\n')
    path = tmp_path / 'fragments.gfa'
    path.write_text(text)
    with pytest.warns(UserWarning, match='Multiple reference haplotypes'):
        graph = PangenomeGraph.from_gfa(str(path), ref_name='ref', dfs_method=dfs_methods[method])
    assert set(graph.reference_paths) == ({'chr1', 'chr2'} if multiple_chromosomes else {'chr1'})
    rows = records(graph, path, tmp_path / 'fragments.vcf')
    row = rows['chr1', '>3>4']
    assert row[3:5] == ['C', 'G']
    assert row[8:] == ['GT:CR:CA', '0|1:1,0:0,1']
    path.write_text(text + 'W\tref\t0\tchr1\t3\t4\t>4\n')
    with pytest.raises(ValueError, match='Multiple reference walks'):
        PangenomeGraph.from_gfa(str(path), ref_name='ref', dfs_method=dfs_methods[method])


@pytest.mark.parametrize('method', list(dfs_methods))
@pytest.mark.parametrize('independent_copy', [False, True])
def test_cleared_metadata_on_a_filtered_cross_chromosome_walk(tmp_path, method, independent_copy):
    path = tmp_path / 'cross.gfa'
    path.write_text((DATA / 'multi_chrom.gfa').read_text() + 'L\t2\t+\t12\t+\t0M\n')
    graph = PangenomeGraph.from_gfa(str(path), ref_name='ref', dfs_method=dfs_methods[method])
    walk = ['1_+', '2_+', '12_+', '14_+']
    expected = graph.genotype(walk)
    terminal_walk = graph.terminal_walk(walk)
    edges = list(zip(terminal_walk, terminal_walk[1:]))
    graph.graph.clear()
    copied = graph.edge_subgraph(edges)
    if independent_copy:
        copied = copied.copy()
    # Only chr1's start and chr2's end are present; the other terminal nodes are hidden.
    assert copied.terminal_walk(walk) == terminal_walk
    genotype = copied.genotype(walk)
    assert genotype.alt_counts == expected.alt_counts
    assert genotype.coverage_by_chrom == expected.coverage_by_chrom
    for edge in edges:
        index = copied.edges[edge]['index']
        assert genotype.get_ref_count(index) == expected.get_ref_count(index)
    assert Genotype.verify_genotype_matches_walk(copied, walk, genotype) == (True, [])
    assert Genotype.verify_genotype_matches_walks(copied, [walk], genotype) == (True, [])


def renamed(fixture, chrom, offset=0):
    rows = []
    for line in (DATA / f'{fixture}.gfa').read_text().splitlines():
        p = line.split('\t')
        if p[0] == 'S':
            p[1] = chrom + '_' + p[1]
        elif p[0] == 'L':
            p[1], p[3] = chrom + '_' + p[1], chrom + '_' + p[3]
        elif p[0] == 'W':
            p[3] = chrom
            p[4:6] = [str(offset), str(offset + int(p[5]))]
            p[6] = re.sub(r'([<>])', lambda m: m[1] + chrom + '_', p[6])
        rows.append('\t'.join(p))
    return '\n'.join(rows) + '\n'


def forest_contract(graph):
    assert nx.is_branching(graph.reference_tree)
    roots = {start + '_+' for start, _ in graph.chromosome_termini.values()}
    assert {n for n, degree in graph.reference_tree.in_degree() if degree == 0} == roots
    assert set(graph.reference_tree) == set(graph.sorted_positive_nodes)
    for root in roots:
        chrom = graph.chromosome(root)
        assert all(graph.chromosome(n) == chrom for n in nx.descendants(graph.reference_tree, root))
    for u, v in graph.reference_tree.edges:
        assert graph.has_edge(u, v)
        assert graph.chromosome(u) == graph.chromosome(v)
    for node in graph:
        assert graph.chromosome(node) == graph.chromosome(node_complement(node))
    for u, v, data in graph.edges(data=True):
        if graph.chromosome(u) != graph.chromosome(v):
            assert not data['is_in_tree']
            assert data['trans_chrom'] == graph.chromosome(v)
            assert graph.edges[edge_complement((u, v))]['trans_chrom'] == graph.chromosome(u)
        else:
            assert 'trans_chrom' not in data


def records(graph, path, output):
    graph.write_vcf(str(path), str(output))
    rows = [line.split('\t') for line in output.read_text().splitlines() if not line.startswith('#')]
    return {(row[0], row[2]): row for row in rows}


@pytest.mark.parametrize('method', list(dfs_methods))
@pytest.mark.parametrize('fixture,length', [('hairpin', 40), ('tandem', 48), ('hairpin_tip', 40)])
def test_two_cyclic_references_keep_steps_coordinates_counts_and_vcf(tmp_path, method, fixture, length):
    path = tmp_path / 'cyclic.gfa'
    path.write_text(renamed(fixture, 'chr1', 100) + renamed(fixture, 'chr2', 1000))
    with pytest.warns(UserWarning, match='Reference walk is cyclic'):
        graph = PangenomeGraph.from_gfa(str(path), ref_name='REF', dfs_method=dfs_methods[method])
    forest_contract(graph)
    rows = records(graph, path, tmp_path / 'cyclic.vcf')
    for chrom, offset in [('chr1', 100), ('chr2', 1000)]:
        start, end = graph.chromosome_termini[chrom]
        expected = {start: offset, end: offset + length,
                    f'{chrom}_1': offset + 8, f'{chrom}_2': offset + 16,
                    f'{chrom}_3': offset + 24, f'{chrom}_4': offset + length}
        for segment, position in expected.items():
            for sign in '+-':
                data = graph.nodes[f'{segment}_{sign}']
                assert (data['position'], data['tree_position'], data['right_position']) == (position,) * 3
        walk = next(w.walk for w in read_gfa_line_by_line(str(path), ['W'])
                    if w.sample_name == 'REF' and w.sequence_name == chrom)
        assert graph.reference_paths[chrom][1:-1] == walk
        assert all(graph.has_edge(u, v) for u, v in zip(walk, walk[1:]))
        genotype = graph.genotype(walk)
        assert Genotype.verify_genotype_matches_walk(graph, walk, genotype) == (True, [])
        assert Genotype.verify_genotype_matches_walks(graph, [walk], genotype) == (True, [])
        if fixture == 'tandem':
            assert graph.is_back_edge((f'{chrom}_3_+', f'{chrom}_2_+'))
            idx = graph.edges[f'{chrom}_2_+', f'{chrom}_3_+']['index']
            assert genotype.get_ref_count(idx) == 2
            row = rows[chrom, f'>{chrom}_3>{chrom}_2']
            assert (int(row[1]), row[3], row[4]) == (offset + 8, 'T', 'TTTGACCAAGGCATTCA')
        else:
            row = rows[chrom, f'>{chrom}_3<{chrom}_2']
            assert (int(row[1]), row[3], row[4]) == (offset + 16, 'A', 'AGGCATTCA')
        assert row[9:] == ['1:1:1', '0:1:0']
        assert 'TRANS_CHROM=' not in row[7]
        if fixture == 'hairpin_tip':
            # There is no biological 3 -> 4 link. The inverted continuation is
            # represented faithfully, and terminal repair spans the positive 4.
            assert not graph.has_edge(f'{chrom}_3_+', f'{chrom}_4_+')


@pytest.mark.parametrize('method', list(dfs_methods))
@pytest.mark.parametrize('fixture', ['acyclic_unreachable', 'disconnected_self_edge', 'self_edge_vcf'])
def test_multichromosome_unreachable_regions_and_self_links(tmp_path, method, fixture):
    ref = 'REF' if fixture != 'disconnected_self_edge' else 'ref'
    path = tmp_path / 'self.gfa'
    path.write_text(renamed(fixture, 'chr1') + renamed(fixture, 'chr2'))
    graph = PangenomeGraph.from_gfa(str(path), ref_name=ref, dfs_method=dfs_methods[method])
    forest_contract(graph)
    # Forest construction may add terminal attachments only.
    biological = {(w.u, w.v) for w in read_gfa_line_by_line(str(path), ['L'])}
    biological |= {edge_complement(e) for e in biological}
    assert {e for e in graph.edges if not any(graph.is_terminal(n) for n in e)} == biological
    genotypes = graph.genotypes_from_gfa(str(path))
    for line in read_gfa_line_by_line(str(path), ['W']):
        gt = graph.genotype(line.walk)
        assert gt.num_edges == graph.biedge_index_capacity
        assert Genotype.verify_genotype_matches_walks(graph, [line.walk], gt) == (True, [])
    if fixture == 'self_edge_vcf':
        rows = records(graph, path, tmp_path / 'self.vcf')
        for chrom in ('chr1', 'chr2'):
            # Links anchored at the artificial root keep the existing terminal
            # exclusion; interior hairpins are exported in ordinary mode.
            for suffix in ('>1<1', '>2<2'):
                variant_id = re.sub(r'([<>])', lambda m: m[1] + chrom + '_', suffix)
                row = rows[chrom, variant_id]
                assert row[9:] == ['0:1:0', '1:1:1']
            edge = (f'{chrom}_1_+', f'{chrom}_1_-')
            assert graph.edges[edge]['is_representative']
            assert graph.edges[edge]['weight'] == 1
            for link in [('1_+', '1_-'), ('1_-', '1_+'), ('1_+', '1_+'), ('2_+', '2_-')]:
                link = tuple(chrom + '_' + n for n in link)
                indices = [graph.edges[e]['index'] for e in graph.reference_tree_edges(link)]
                assert genotypes['REF'][0].variant_record(link, None, indices) == (0, 1, 0)
                assert genotypes['ALT'][0].variant_record(link, None, indices) == (1, 1, 1)
        assert len(graph.sorted_biedge_representatives) == graph.number_of_biedges
        assert genotypes['ALT'][0].num_edges > max(d['index'] for *_, d in graph.edges(data=True))


@pytest.mark.parametrize('method', list(dfs_methods))
@pytest.mark.parametrize('shared', ['>2', '<2'])
def test_shared_segments_keep_lexicographic_owner_without_new_biological_links(tmp_path, method, shared):
    text = (DATA / 'multi_chrom.gfa').read_text()
    # Present chr2 first to check lexical rather than input ownership.
    rows = text.splitlines()
    rows.sort(key=lambda s: 0 if s.startswith('W\tref\t0\tchr2') else 1)
    text = '\n'.join(rows) + '\n'
    text = text.replace('>11>12>14', '>11' + shared + '>14')
    sign = '+' if shared == '>2' else '-'
    text += f'L\t11\t+\t2\t{sign}\t0M\nL\t2\t{sign}\t14\t+\t0M\n'
    path = tmp_path / 'shared.gfa'
    path.write_text(text)
    graph = PangenomeGraph.from_gfa(str(path), ref_name='ref', dfs_method=dfs_methods[method])
    forest_contract(graph)
    assert graph.chromosome('2_+') == 'chr1'
    assert graph.chromosome('14_+') == graph.chromosome('11_+') == 'chr2'
    assert not graph.has_edge('11_+', '14_+')
    assert graph.reference_paths['chr2'][1:-1] == ['11_+', '2_' + sign, '14_+']
    rows = records(graph, path, tmp_path / 'shared.vcf')
    assert any(row[0] == 'chr2' and 'TRANS_CHROM=chr1' in row[7] for row in rows.values())


@pytest.mark.parametrize('method', list(dfs_methods))
def test_forest_dfs_is_not_reference_replay(tmp_path, method):
    graph = nx.DiGraph()
    graph.add_weighted_edges_from([('a', 'b', 1), ('b', 'c', 1), ('c', 'b', 1),
                                   ('b', 'd', 1), ('c', 'd', 1), ('d', 'a', 1),
                                   ('d', 'x', 99), ('x', 'y', 1), ('y', 'b', 99)])
    with pytest.warns(UserWarning, match='Reference walk is cyclic'):
        tree = dfs_methods[method](graph, reference_path=[],
                                   reference_paths=[['a', 'b', 'c', 'b', 'd'], ['x', 'y']],
                                   haplo_priorities={})
    assert set(tree.edges) == {('a', 'b'), ('b', 'c'), ('c', 'd'), ('x', 'y')}
    assert nx.is_branching(tree)


@pytest.mark.parametrize('method', list(dfs_methods))
def test_multichromosome_contiguous_priority_is_retained(tmp_path, method):
    path = tmp_path / 'priority.gfa'
    path.write_text(renamed('contiguous_priority', 'chr1') + renamed('contiguous_priority', 'chr2'))
    for priorities, expected in [({'sample1': 0, 'sample2': 1}, '2'),
                                 ({'sample2': 0, 'sample1': 1}, '3')]:
        graph = PangenomeGraph.from_gfa(str(path), ref_name='ref',
                                       priority_dict=priorities, dfs_method=dfs_methods[method])
        forest_contract(graph)
        if method == 'contiguous':
            for chrom in ('chr1', 'chr2'):
                assert graph.parent_in_tree(f'{chrom}_5_+') == f'{chrom}_{expected}_+'


@pytest.mark.parametrize('fixture,ref', [('simple_nested', 'sample1'),
                                         ('c4a_with_inversion_and_sequences', 'HG02257')])
def test_diploid_reference_selects_first_haplotype_with_warning(fixture, ref):
    path = DATA / f'{fixture}.gfa'
    walks = [w for w in read_gfa_line_by_line(str(path), ['W']) if w.sample_name == ref]
    with pytest.warns(UserWarning, match='Multiple reference haplotypes'):
        graph = PangenomeGraph.from_gfa(str(path), ref_name=ref)
    assert len(graph.reference_paths) == 1
    assert graph.reference_path[1:-1] == walks[0].walk
    assert nx.is_arborescence(graph.reference_tree)


def test_alternate_haplotype_selection_and_duplicate_fragment_errors(tmp_path):
    path = tmp_path / 'ref.gfa'
    text = (DATA / 'multi_chrom.gfa').read_text()
    alternate = 'W\tref\t1\tassembly1\t0\t6\t>1>2>4\nW\tref\t1\tassembly2\t0\t6\t>11>12>14\n'
    path.write_text(text + alternate)
    with pytest.warns(UserWarning, match='Multiple reference haplotypes'):
        graph = PangenomeGraph.from_gfa(str(path), ref_name='ref')
    assert list(graph.reference_paths) == ['chr1', 'chr2']
    path.write_text(text + 'W\tref\t0\tchr2\t6\t9\t>14\n')
    with pytest.raises(ValueError, match='Multiple reference walks'):
        PangenomeGraph.from_gfa(str(path), ref_name='ref')
    # A repeated record is equally ambiguous even if the coordinates are equal.
    path.write_text(text + 'W\tref\t0\tchr2\t0\t6\t>11>12>14\n')
    with pytest.raises(ValueError, match='Multiple reference walks'):
        PangenomeGraph.from_gfa(str(path), ref_name='ref')


@pytest.mark.parametrize('method', list(dfs_methods))
def test_single_chromosome_offset_cli_override_and_cyclic_contiguous(tmp_path, method):
    path = tmp_path / 'offset.gfa'
    path.write_text(renamed('hairpin', 'chr1', 100))
    output = tmp_path / 'override.vcf'
    result = CliRunner().invoke(cli, ['gfa2vcf', str(path), str(output), '--ref-name', 'REF',
                                      '--chr-id', 'override', '--dfs-method', method])
    assert result.exit_code == 0, result.output
    rows = [r.split('\t') for r in output.read_text().splitlines() if not r.startswith('#')]
    row = next(row for row in rows if row[2] == '>chr1_3<chr1_2')
    assert (row[0], row[1], row[3], row[4], row[9:]) == ('override', '116', 'A', 'AGGCATTCA', ['1:1:1', '0:1:0'])


@pytest.mark.parametrize('method', list(dfs_methods))
def test_cyclic_multichromosome_output_is_hash_seed_deterministic(tmp_path, method):
    path = tmp_path / 'deterministic.gfa'
    path.write_text(renamed('hairpin_tip', 'chr1') + renamed('tandem', 'chr2') +
                    'S\tunreachable\tAG\nL\tunreachable\t+\tunreachable\t+\t0M\n'
                    'L\tchr1_2\t+\tchr2_2\t-\t0M\n')
    code = '''from pantree.graph import PangenomeGraph
from pantree.dfs import dfs_methods
import sys
G = PangenomeGraph.from_gfa(sys.argv[1], ref_name='REF', dfs_method=dfs_methods[sys.argv[3]])
G.write_vcf(sys.argv[1], sys.argv[2])
print(sorted(G.reference_tree.edges))
'''
    results = []
    for seed in ('0', '2'):
        output = tmp_path / (seed + '.vcf')
        run = subprocess.run([sys.executable, '-c', code, str(path), str(output), method],
                             env={**os.environ, 'PYTHONHASHSEED': seed}, capture_output=True, text=True, check=True)
        results.append((run.stdout, output.read_bytes()))
    assert results[0] == results[1]


@pytest.mark.parametrize('method', list(dfs_methods))
def test_repeated_self_links_export_normalized_counts_on_two_chromosomes(tmp_path, method):
    path = tmp_path / 'repeated.gfa'
    path.write_text(renamed('self_edge_repeated', 'chr1') + renamed('self_edge_repeated', 'chr2'))
    graph = PangenomeGraph.from_gfa(str(path), ref_name='ref', dfs_method=dfs_methods[method])
    forest_contract(graph)
    rows = records(graph, path, tmp_path / 'repeated.vcf')
    assert len(rows) == 6
    for chrom in ('chr1', 'chr2'):
        for prefix, count, pos, alleles in [('>a>a', 3, 2, ('A', 'AACG')),
                                            ('>a<a', 2, 5, ('G', 'G')),
                                            ('<a>a', 2, 2, ('.', '.'))]:
            variant_id = re.sub(r'([<>])', lambda m: m[1] + chrom + '_', prefix)
            row = rows[chrom, variant_id]
            assert (int(row[1]), row[3], row[4]) == (pos, *alleles)
            assert row[9:] == ['0:1:0', f'1:1:{count}']
            assert f'AC={count}' in row[7].split(';')


@pytest.mark.parametrize('method', list(dfs_methods))
def test_reachability_repair_cannot_steal_an_existing_chromosome_branch(tmp_path, method):
    # The first chromosome's 4 is reachable only via inversion or its self-loop.
    # The second chromosome reaches its 4 via an off-reference branch 'a'.
    # Repairing chr1 must preserve a's initial chr2 ancestry, even with 4 -> a.
    text = ''
    for chrom in ('chr1', 'chr2'):
        part = renamed('hairpin_tip', chrom)
        part = '\n'.join(line for line in part.splitlines()
                         if not line.startswith('W\tALT') and
                         line != f'L\t{chrom}_2\t+\t{chrom}_4\t+\t0M') + '\n'
        text += part + f'L\t{chrom}_4\t+\t{chrom}_4\t+\t0M\n'
    text += ('S\ta\tG\nL\tchr2_2\t+\ta\t+\t0M\n'
             'L\ta\t+\tchr2_4\t+\t0M\nL\tchr1_4\t+\ta\t+\t0M\n')
    path = tmp_path / 'repair.gfa'
    path.write_text(text)
    with pytest.warns(UserWarning, match='unreachable'):
        graph = PangenomeGraph.from_gfa(str(path), ref_name='REF', dfs_method=dfs_methods[method])
    forest_contract(graph)
    assert graph.chromosome('a_+') == 'chr2'
    assert graph.parent_in_tree('chr2_4_+') == 'a_+'
    assert graph.edges['chr1_4_+', 'a_+']['trans_chrom'] == 'chr2'


@pytest.mark.parametrize('method', list(dfs_methods))
@pytest.mark.parametrize('fixture', ['hairpin', 'tandem'])
def test_shared_segments_in_two_cyclic_references(tmp_path, method, fixture):
    second = renamed(fixture, 'chr2').replace('chr2_2', 'chr1_2')
    second = '\n'.join(line for line in second.splitlines() if not line.startswith('S\tchr1_2\t')) + '\n'
    path = tmp_path / 'shared_cyclic.gfa'
    path.write_text(renamed(fixture, 'chr1') + second)
    with pytest.warns(UserWarning, match='Reference walk is cyclic'):
        graph = PangenomeGraph.from_gfa(str(path), ref_name='REF', dfs_method=dfs_methods[method])
    forest_contract(graph)
    assert graph.chromosome('chr1_2_+') == 'chr1'
    assert graph.chromosome('chr2_3_+') == graph.chromosome('chr2_4_+') == 'chr2'
    assert not graph.has_edge('chr2_1_+', 'chr2_3_+')
    assert graph.nodes['chr1_2_+']['position'] == 16
    end = graph.chromosome_termini['chr2'][1] + '_+'
    assert graph.position(end) == (40 if fixture == 'hairpin' else 48)
    rows = records(graph, path, tmp_path / 'shared_cyclic.vcf')
    row = rows['chr2', '>chr2_1>chr1_2']
    assert 'TRANS_CHROM=chr1' in row[7].split(';')
    assert row[9:] == ['1:0:1', '1:0:1']


@pytest.mark.parametrize('method', list(dfs_methods))
def test_branched_cyclic_missingness_remains_chromosome_specific(tmp_path, method):
    path = tmp_path / 'missing.gfa'
    path.write_text(renamed('cyclic_missingness', 'chr1', 100) +
                    renamed('cyclic_missingness', 'chr2', 1000) +
                    'W\tonly1\t0\tcontig\t0\t6\t>chr1_1>chr1_2>chr1_3>chr1_2>chr1_4>chr1_5\n')
    with pytest.warns(UserWarning, match='Reference walk is cyclic'):
        graph = PangenomeGraph.from_gfa(str(path), ref_name='REF', dfs_method=dfs_methods[method])
    forest_contract(graph)
    output = tmp_path / 'missing.vcf'
    rows = records(graph, path, output)
    header = next(line.split('\t') for line in output.read_text().splitlines() if line.startswith('#CHROM'))
    for chrom in ('chr1', 'chr2'):
        row = rows[chrom, f'>{chrom}_4>{chrom}_5']
        assert row[header.index('PARTIAL')] == '.:0:0'
        # The cyclic reference itself traverses this non-tree link.
        assert row[header.index('only1')] == ('1:0:1' if chrom == 'chr1' else '.:0:0')

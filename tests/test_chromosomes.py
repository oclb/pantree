from pathlib import Path
import gzip
import re

import networkx as nx
import pytest
from click.testing import CliRunner

from pantree.cli import cli
from pantree.dfs import dfs_methods
from pantree.gfa import GFAWalkLine
from pantree.graph import PangenomeGraph
from pantree.read_vcf import read_vcf_to_lazyframe
from pantree.utils import edge_complement, node_complement
from pantree.walk_variants import process_haplotype_variants


FIXTURE = Path('tests/data/multi_chrom.gfa')


def load(tmp_path, extra='', transform=lambda s: s, method='max_weight'):
    path = tmp_path / 'graph.gfa'
    path.write_text(transform(FIXTURE.read_text()) + extra)
    graph = PangenomeGraph.from_gfa(str(path), ref_name='ref', dfs_method=dfs_methods[method])
    return graph, path


def export(graph, path, tmp_path, compressed=False):
    output = tmp_path / ('out.vcf.gz' if compressed else 'out.vcf')
    graph.write_vcf(str(path), str(output))
    text = gzip.open(output, 'rt').read() if compressed else output.read_text()
    return output, text, read_vcf_to_lazyframe(str(output)).collect()


@pytest.mark.parametrize('method', list(dfs_methods))
def test_disconnected_chromosomes_match_separate_runs(tmp_path, method):
    graph, path = load(tmp_path, method=method)
    _, text, combined = export(graph, path, tmp_path, compressed=True)
    assert set(combined['#CHROM']) == {'chr1', 'chr2'}
    assert combined['POS'].to_list() == [4, 4]
    assert combined['REF'].to_list() == ['C', 'A']
    assert combined['ALT'].to_list() == ['T', 'G']
    assert combined['sample_GT'].to_list() == ['1|0', '1|0']
    assert 'TRANS_CHROM=' not in '\n'.join(l for l in text.splitlines() if not l.startswith('#'))
    assert len(graph.terminal_binodes) == 4
    assert nx.is_forest(graph.reference_tree)
    for node in graph:
        assert graph.chromosome(node) == graph.chromosome(node_complement(node))
    for u, v in graph.reference_tree.edges:
        assert graph.chromosome(u) == graph.chromosome(v)

    for chrom in ('chr1', 'chr2'):
        rows = []
        for line in FIXTURE.read_text().splitlines():
            p = line.split('\t')
            if p[0] == 'H' or (p[0] in ('S', 'L') and (int(p[1]) < 10) == (chrom == 'chr1')):
                rows.append(line)
            elif p[0] == 'W' and (int(p[6].split('>')[1]) < 10) == (chrom == 'chr1'):
                rows.append(line)
        separate_path = tmp_path / (chrom + '.gfa')
        separate_path.write_text('\n'.join(rows) + '\n')
        separate = PangenomeGraph.from_gfa(str(separate_path), ref_name='ref', dfs_method=dfs_methods[method])
        _, _, result = export(separate, separate_path, tmp_path)
        columns = ['POS', 'ID', 'REF', 'ALT', 'VT', 'TP', 'sample_GT', 'sample_CR', 'sample_CA']
        assert result.select(columns).rows() == combined.filter(combined['#CHROM'] == chrom).select(columns).rows()


@pytest.mark.parametrize('orientation', ['++', '+-', '-+', '--'])
@pytest.mark.parametrize('method', list(dfs_methods))
def test_trans_chrom_is_scalar_and_follows_source(tmp_path, orientation, method):
    a, b = orientation
    edge = ('2_' + a, '12_' + b)
    walk = ('>1>2' if a == '+' else '<4<2') + ('>12>14' if b == '+' else '<12<11')
    extra = f'L\t2\t{a}\t12\t{b}\t0M\nW\tcross\t0\tcontig\t0\t8\t{walk}\n'
    graph, path = load(tmp_path, extra, method=method)
    assert edge in graph.variant_edges
    assert not graph.edges[edge]['is_in_tree']
    assert graph.edges[edge]['trans_chrom'] == 'chr2'
    assert graph.edges[edge_complement(edge)]['trans_chrom'] == 'chr1'
    output, text, data = export(graph, path, tmp_path)
    row = data.filter(data['TRANS_CHROM'].is_not_null()).row(0, named=True)
    assert row['#CHROM'] == 'chr1'
    assert row['TRANS_CHROM'] == 'chr2'
    assert row['VT'] == row['SVTYPE'] == 'BND'
    assert row['POS'] == 4
    assert row['REF'] == 'C'
    assert row['ALT'] == ('C.' if a == '+' else '.C')
    assert row['cross_GT'] == '1'
    assert '##INFO=<ID=TRANS_CHROM,Number=1,Type=String' in text
    with pytest.raises(ValueError, match='cross-chromosome'):
        process_haplotype_variants(str(output), 'cross', 0, str(tmp_path / 'consolidated.vcf'))


def test_cross_chromosome_cycle_and_tree_dependent_assignment(tmp_path):
    extra = 'S\t20\tAC\nL\t2\t+\t20\t+\t0M\nL\t12\t+\t20\t+\t0M\nL\t20\t+\t4\t+\t0M\nL\t20\t+\t14\t+\t0M\nL\t14\t+\t1\t+\t0M\n'
    graph, path = load(tmp_path, extra)
    assert graph.chromosome('20_+') in {'chr1', 'chr2'}
    for u, v in graph.edges:
        if graph.chromosome(u) != graph.chromosome(v):
            assert not graph.edges[u, v]['is_in_tree']
            assert graph.edges[u, v]['trans_chrom'] == graph.chromosome(v)
    export(graph, path, tmp_path)


def test_missing_chromosome_is_not_covered_by_other_chromosome(tmp_path):
    graph, path = load(tmp_path, 'W\tonly1\t0\tcontig\t0\t6\t>1>2>4\n')
    _, _, data = export(graph, path, tmp_path)
    assert data['only1_GT'].to_list() == ['0', '.']


@pytest.mark.parametrize('source_direction', ['+', '-'])
def test_trans_reference_counts_use_source_chromosome(tmp_path, source_direction):
    extra = (f'L\t2\t{source_direction}\t12\t+\t0M\n'
             'W\tonly2\t0\tcontig\t0\t6\t>11>12>14\n'
             'W\tonly1\t0\tcontig\t0\t6\t>1>2>4\n')
    graph, path = load(tmp_path, extra)
    _, _, data = export(graph, path, tmp_path)
    row = data.filter(data['TRANS_CHROM'].is_not_null()).row(0, named=True)
    assert (row['only2_GT'], row['only2_CR'], row['only2_CA']) == ('.', '0', '0')
    assert (row['only1_GT'], row['only1_CR'], row['only1_CA']) == ('0', '1', '0')


def test_trans_source_counts_sum_multiple_tree_continuations(tmp_path):
    extra = ('S\t20\tAG\nL\t2\t+\t20\t+\t0M\nL\t20\t+\t4\t+\t0M\n'
             'L\t2\t+\t12\t+\t0M\nW\tbranch\t0\tx\t0\t8\t>1>2>20>4\n')
    graph, path = load(tmp_path, extra)
    assert len(graph.reference_tree_edges(('2_+', '12_+'))) == 2
    _, _, data = export(graph, path, tmp_path)
    row = data.filter(data['TRANS_CHROM'].is_not_null()).row(0, named=True)
    assert row['branch_CR'] == '1'
    assert row['sample_CR'] == '0,1'


def test_mc_fragmented_walks_keep_haplotype_and_missingness(tmp_path):
    def clipped(text):
        return text.replace('W\tsample\t1\tcontigB\t0\t6\t>11>13>14',
                            'W\tsample\t1\tcontigB\t0\t3\t>11\nW\tsample\t1\tcontigB\t4\t6\t>14')
    graph, path = load(tmp_path, transform=clipped)
    _, _, data = export(graph, path, tmp_path)
    assert data['sample_GT'].to_list() == ['1|0', '.|0']


def test_shared_reference_segment_has_one_owner(tmp_path):
    def shared(text):
        return text.replace('>11>12>14', '>11>2>14') + 'L\t11\t+\t2\t+\t0M\nL\t2\t+\t14\t+\t0M\n'
    graph, path = load(tmp_path, transform=shared)
    assert graph.chromosome('2_+') == 'chr1'
    assert graph.edges['11_+', '2_+']['trans_chrom'] == 'chr1'
    export(graph, path, tmp_path)


def test_cli_infers_chromosomes_and_rejects_global_override(tmp_path):
    output = tmp_path / 'all.vcf.gz'
    runner = CliRunner()
    args = ['gfa2vcf', str(FIXTURE), str(output), '--ref-name', 'ref']
    assert runner.invoke(cli, args).exit_code == 0
    result = runner.invoke(cli, args + ['--chr-id', 'wrong'])
    assert result.exit_code != 0
    assert 'single-chromosome' in str(result.exception)


def test_pansn_paths_and_w_sequence_names():
    path = GFAWalkLine.from_parts_P(['P', 'ref#0#chr2', '1+,2+', '*'])
    assert path.sample_name == 'ref'
    assert path.sequence_name == 'chr2'
    walk = GFAWalkLine.from_parts_W(['W', 'sample', '1', 'assembly_contig', '100', '102', '>1>2'])
    assert walk.sequence_name == 'assembly_contig'
    assert walk.contig_start == 100


def test_optional_info_fields_can_be_reordered(tmp_path):
    path = tmp_path / 'optional.vcf'
    path.write_text('##fileformat=VCFv4.2\n##INFO=<ID=VT,Number=1,Type=String,Description="Type">\n'
                    '##INFO=<ID=TRANS_CHROM,Number=1,Type=String,Description="Other chromosome">\n'
                    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                    'chr1\t1\ta\tA\tT\t.\tPASS\tVT=SNP\n'
                    'chr1\t2\tb\tC\tC.\t.\tPASS\tTRANS_CHROM=chr2;VT=BND\n')
    data = read_vcf_to_lazyframe(str(path)).collect()
    assert data['VT'].to_list() == ['SNP', 'BND']
    assert data['TRANS_CHROM'].to_list() == [None, 'chr2']


@pytest.mark.parametrize('method', list(dfs_methods))
def test_existing_inversion_and_duplication_graph_as_two_chromosomes(tmp_path, method):
    original = Path('tests/data/c4a_with_inversion_and_sequences.gfa').read_text()
    def renamed(chrom):
        result = []
        for line in original.splitlines():
            parts = line.split('\t')
            if parts[0] == 'S':
                parts[1] = chrom + '_' + parts[1]
            elif parts[0] == 'L':
                parts[1] = chrom + '_' + parts[1]
                parts[3] = chrom + '_' + parts[3]
            elif parts[0] == 'W':
                parts[3] = chrom
                parts[6] = re.sub(r'([<>])', lambda m: m[1] + chrom + '_', parts[6])
            result.append('\t'.join(parts))
        return '\n'.join(result) + '\n'
    path = tmp_path / 'two.gfa'
    path.write_text(renamed('chr1') + renamed('chr2'))
    graph = PangenomeGraph.from_gfa(str(path), ref_name='GRCh38', dfs_method=dfs_methods[method])
    _, _, combined = export(graph, path, tmp_path)
    for chrom in ('chr1', 'chr2'):
        path.write_text(renamed(chrom))
        single = PangenomeGraph.from_gfa(str(path), ref_name='GRCh38', dfs_method=dfs_methods[method])
        _, _, separate = export(single, path, tmp_path)
        cols = ['POS', 'ID', 'REF', 'ALT', 'VT', 'TP'] + [c for c in separate.columns if c.endswith(('_GT', '_CR', '_CA'))]
        assert separate.select(cols).sort('ID').rows() == combined.filter(combined['#CHROM'] == chrom).select(cols).sort('ID').rows()


def test_reference_offsets_and_invalid_reference_inputs(tmp_path):
    graph, path = load(tmp_path, transform=lambda s: s.replace('ref\t0\tchr2\t0\t6', 'ref\t0\tchr2\t100\t106'))
    _, _, data = export(graph, path, tmp_path)
    assert data['POS'].to_list() == [4, 104]
    result = CliRunner().invoke(cli, ['gfa2vcf', str(path), str(tmp_path / 'missing.vcf'), '--ref-name', 'absent'])
    assert result.exit_code != 0
    assert 'No reference walks' in result.output
    path.write_text(FIXTURE.read_text() + 'W\tref\t0\tchr1\t6\t9\t>1\n')
    with pytest.raises(ValueError, match='Multiple reference walks'):
        PangenomeGraph.from_gfa(str(path), ref_name='ref')


def test_multi_chromosome_pansn_paths(tmp_path):
    def paths(text):
        rows = []
        for line in text.splitlines():
            p = line.split('\t')
            if p[0] == 'W':
                path = ','.join(node + '+' for node in p[6].split('>')[1:])
                line = f'P\t{p[1]}#{p[2]}#{p[3]}\t{path}\t*'
            rows.append(line)
        return '\n'.join(rows) + '\n'
    graph, path = load(tmp_path, transform=paths)
    _, _, data = export(graph, path, tmp_path)
    assert data['sample_GT'].to_list() == ['1|0', '1|0']
    assert data['#CHROM'].to_list() == ['chr1', 'chr2']

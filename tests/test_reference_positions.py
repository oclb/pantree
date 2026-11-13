"""First-occurrence reference coordinates and their downstream projections."""
from pathlib import Path

import pytest

from pantree.graph import PangenomeGraph
from pantree.gfa import read_gfa_line_by_line


DATA = Path(__file__).parent / 'data'


@pytest.mark.parametrize('fixture, length, last_third', [
    ('hairpin', 40, 24),
    ('tandem', 48, 40),
])
@pytest.mark.parametrize('contig_start', [0, 100])
def test_revisited_segments_keep_first_positions(fixture, length, last_third,
                                                contig_start, tmp_path):
    lines = (DATA / f'{fixture}.gfa').read_text().splitlines()
    for i, line in enumerate(lines):
        fields = line.split('\t')
        if fields[:2] == ['W', 'REF']:
            fields[4:6] = [str(contig_start), str(contig_start + length)]
            lines[i] = '\t'.join(fields)
    path = tmp_path / 'reference.gfa'
    path.write_text('\n'.join(lines) + '\n')
    with pytest.warns(UserWarning, match='Reference walk is cyclic'):
        graph = PangenomeGraph.from_gfa(str(path), ref_name='REF')

    # Positions are segment ends; revisits count toward later nodes and the terminus.
    expected = {'+_terminus': 0, '1': 8, '2': 16, '3': 24,
                '4': length, '-_terminus': length}
    for segment, position in expected.items():
        for orientation in ['+', '-']:
            data = graph.nodes[f'{segment}_{orientation}']
            assert data['position'] == position
            assert data['tree_position'] == position
            assert data['right_position'] == position

    # This base uses W offsets only for haplotype annotations, which retain
    # their existing last-occurrence convention independently of reference positions.
    for segment, last_position in [('1', 8), ('2', 32),
                                   ('3', last_third), ('4', length)]:
        for orientation in ['+', '-']:
            assert graph.nodes[f'{segment}_{orientation}']['REF#0#chr1'] == (
                contig_start + last_position
            )


def test_end_terminus_keeps_full_length_when_last_segment_is_revisited(tmp_path):
    path = tmp_path / 'revisited_end.gfa'
    path.write_text((DATA / 'tandem.gfa').read_text().replace(
        'W\tREF\t0\tchr1\t0\t48\t>1>2>3>2>3>4',
        'W\tREF\t0\tchr1\t0\t40\t>1>2>3>2>3',
    ))
    with pytest.warns(UserWarning, match='Reference walk is cyclic'):
        graph = PangenomeGraph.from_gfa(str(path), ref_name='REF')
    for orientation in ['+', '-']:
        for key in ['position', 'tree_position', 'right_position']:
            assert graph.nodes[f'3_{orientation}'][key] == 24
            assert graph.nodes[f'-_terminus_{orientation}'][key] == 40


@pytest.mark.parametrize('fixture', ['simple_nested', 'simple_gfa10'])
def test_acyclic_reference_coordinates_and_vcf_are_unchanged(fixture, tmp_path):
    path = DATA / f'{fixture}.gfa'
    graph = PangenomeGraph.from_gfa(str(path), ref_name='ref')
    # Both W and P paths keep the same reference and off-reference projections.
    expected = {
        '+_terminus': (0, 0, 0), '1': (5, 5, 5), '2': (10, 10, 10),
        '3': (5, 10, 16), '4': (16, 16, 16), '5': (16, 21, 20),
        '6': (16, 22, 20), '7': (16, 22, 20), '8': (16, 26, 20),
        '9': (20, 20, 20), '10': (25, 25, 25), '11': (29, 29, 29),
        '-_terminus': (29, 29, 29),
    }
    for segment, positions in expected.items():
        for orientation in ['+', '-']:
            data = graph.nodes[f'{segment}_{orientation}']
            assert tuple(data[key] for key in ['position', 'tree_position',
                                               'right_position']) == positions

    output = tmp_path / 'acyclic.vcf'
    graph.write_vcf(str(path), str(output), chr_name='chr1')
    rows = [line.split('\t') for line in output.read_text().splitlines()
            if not line.startswith('#')]
    assert [(row[2], int(row[1]), row[3], row[4]) for row in rows] == [
        ('>3>4', 6, 'GGCAC', 'CTTAG'),
        ('>8>9', 16, 'C', 'CCGGTAGGTAC'),
        ('>7>8', 17, '.', 'A'),
        ('>9>11', 20, 'ATCATG', 'A'),
    ]


@pytest.mark.parametrize('fixture, length, variant_id, variant_alleles', [
    ('hairpin', 40, '>3<2', (16, 'A', 'AGGCATTCA')),
    ('tandem', 48, '>3>2', (8, 'T', 'TTTGACCAAGGCATTCA')),
])
def test_revisit_vcf_alleles_and_suffix_missingness(fixture, length, variant_id,
                                                  variant_alleles, tmp_path):
    path = tmp_path / 'with_branch.gfa'
    path.write_text((DATA / f'{fixture}.gfa').read_text() +
                    'S\t5\tG\n'
                    'L\t2\t+\t5\t+\t0M\n'
                    'L\t5\t+\t3\t+\t0M\n'
                    'W\tBUBBLE\t0\tchr1\t0\t33\t>1>2>5>3>4\n'
                    'W\tTAIL\t0\tchr1\t0\t16\t>3>4\n')
    with pytest.warns(UserWarning, match='Reference walk is cyclic'):
        graph = PangenomeGraph.from_gfa(str(path), ref_name='REF')

    # Off-reference nodes inherit the first position of their reference ancestor.
    for node in ['5_+', '5_-']:
        assert graph.nodes[node]['position'] == 16
        assert graph.nodes[node]['tree_position'] == 17
        assert graph.nodes[node]['right_position'] == 24

    genotypes = graph.genotypes_from_gfa(str(path))
    for sample in ['REF', 'ALT', 'BUBBLE']:
        assert genotypes[sample][0].linear_coverage == [(8, length)]
        assert genotypes[sample][0].missing_variants == set()
    tail = genotypes['TAIL'][0]
    assert tail.linear_coverage == [(24, length)]
    revisit_edge = ('3_+', '2_-' if fixture == 'hairpin' else '2_+')
    assert tail.missing_variants == {revisit_edge, ('5_+', '3_+')}

    output = tmp_path / 'cyclic.vcf'
    graph.write_vcf(str(path), str(output), chr_name='chr1')
    rows = {row[2]: row for line in output.read_text().splitlines()
            if not line.startswith('#') for row in [line.split('\t')]}
    expected_alleles = {variant_id: variant_alleles, '>5>3': (16, 'A', 'AG')}
    if fixture == 'hairpin':
        expected_alleles['<2>4'] = (8, '.', 'TGAATGCCTTGGTCAA')
        assert rows['<2>4'][9:] == ['1:0:1', '0:1:0', '0:1:0', '0:1:0']
    assert {key: (int(row[1]), row[3], row[4])
            for key, row in rows.items()} == expected_alleles
    for key, (position, _, _) in expected_alleles.items():
        assert f'TP={position}' in rows[key][7].split(';')
    assert rows[variant_id][9:] == ['1:1:1', '0:1:0', '0:1:0', '.:0:0']
    assert rows['>5>3'][9:] == ['0:1:0' if fixture == 'hairpin' else '0:2:0',
                              '0:1:0', '1:0:1', '.:0:0']

    # Disabling inference changes only the absent suffix genotypes to called zeros.
    unmasked = tmp_path / 'unmasked.vcf'
    graph.write_vcf(str(path), str(unmasked), chr_name='chr1', no_missingness=True)
    unmasked_rows = {row[2]: row for line in unmasked.read_text().splitlines()
                     if not line.startswith('#') for row in [line.split('\t')]}
    for key, row in rows.items():
        assert unmasked_rows[key][:-1] == row[:-1]
        assert unmasked_rows[key][-1] == ('0:0:0' if key in [variant_id, '>5>3']
                                         else row[-1])


@pytest.mark.parametrize('exclude_terminus', [True, False])
def test_branched_cyclic_missingness_matches_coverage_gap_scan(exclude_terminus,
                                                            tmp_path):
    path = DATA / 'cyclic_missingness.gfa'
    with pytest.warns(UserWarning, match='Reference walk is cyclic'):
        graph = PangenomeGraph.from_gfa(str(path), ref_name='REF')
    # Branch-point order need not sort the endpoints used for missingness.
    endpoint_positions = [min(graph.position(edge))
                          for edge in graph.sorted_variant_edges()]
    assert endpoint_positions != sorted(endpoint_positions)

    for line in read_gfa_line_by_line(str(path), line_types=['W']):
        genotype = graph.genotype(line.walk, exclude_terminus=exclude_terminus)
        genotype.compute_missing_variants(graph)
        gap_starts = sorted([end for _, end in genotype.linear_coverage]
                            + [graph.position('+_terminus_+')])
        gap_ends = sorted([start for start, _ in genotype.linear_coverage]
                          + [graph.right_position('-_terminus_+')])
        gaps = list(zip(gap_starts, gap_ends))
        expected = {
            edge for edge in graph.variant_edges
            if (not exclude_terminus or not graph.is_terminal(edge))
            and any(start < min(graph.position(edge)) < end
                    and max(graph.right_position(edge)) <= end
                    for start, end in gaps)
        }
        assert genotype.missing_variants == expected
        if line.sample_name == 'PARTIAL':
            assert ('4_+', '5_+') in genotype.missing_variants

    output = tmp_path / 'partial.vcf'
    graph.write_vcf(str(path), str(output), chr_name='chr1',
                    exclude_terminus=exclude_terminus)
    rows = [line.split('\t') for line in output.read_text().splitlines()]
    header = next(row for row in rows if row[0] == '#CHROM')
    uncovered = next(row for row in rows if len(row) > 2 and row[2] == '>4>5')
    assert uncovered[header.index('PARTIAL')] == '.:0:0'

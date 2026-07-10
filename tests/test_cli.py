from pantree.cli import cli
from pantree.graph import PangenomeGraph
from click.testing import CliRunner
import gzip
import json
import os
import tempfile


def test_main_simple_nested():
    """Test main function with simple_nested.gfa"""
    runner = CliRunner()
    
    # Get the path to the test data file
    test_dir = os.path.dirname(__file__)
    gfa_file = os.path.join(test_dir, "data", "simple_nested.gfa")
    
    # Create a temporary output file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as tmp:
        vcf_file = tmp.name
    
    try:
        # Run the main function
        result = runner.invoke(cli, ['gfa2vcf', gfa_file, vcf_file, '--chr-id', 'chr1', '--ref-name', 'ref'])
        
        # Print output for inspection
        print(f"Exit code: {result.exit_code}")
        print(f"Output: {result.output}")
        if result.exception:
            print(f"Exception: {result.exception}")
            import traceback
            traceback.print_exception(type(result.exception), result.exception, result.exception.__traceback__)
        
        # Check if VCF file was created
        assert os.path.exists(vcf_file), "VCF file was not created"
        
        with open(vcf_file, 'r') as f:
            vcf_contents = f.read()
            print(f"\nVCF file created at: {vcf_file}")
            print(f"VCF contents:\n{vcf_contents}")
            
            lines = vcf_contents.strip().split('\n')
            
            # Find the header line
            header_line = None
            variant_lines = []
            for line in lines:
                if line.startswith('#CHROM'):
                    header_line = line
                elif not line.startswith('#') and line.strip():
                    variant_lines.append(line)
            
            assert header_line is not None, "Header line not found in VCF"
            
            # Parse header to get sample names
            header_fields = header_line.split('\t')
            sample_names = header_fields[9:]  # Samples start after FORMAT column
            
            # Assert 1: Verify we have all expected samples (now aggregated by sample name)
            # With the new genotype system, samples are aggregated: 'ref', 'sample1', 'sample2'
            assert any('ref' in s for s in sample_names), f"ref sample not found in VCF. Samples: {sample_names}"
            assert any('sample1' in s for s in sample_names), f"sample1 not found in VCF. Samples: {sample_names}"
            assert any('sample2' in s for s in sample_names), f"sample2 not found in VCF. Samples: {sample_names}"
            assert len(sample_names) == 3, f"Expected 3 samples (ref, sample1, sample2), got {len(sample_names)}: {sample_names}"
            print(f"\n✓ All 3 samples present: {sample_names}")
            
            # Assert 2: Verify we have the expected number of variants (4 variants in the reference file)
            assert len(variant_lines) == 4, f"Expected 4 variants, got {len(variant_lines)}"
            print(f"✓ Found 4 variants as expected")
            
            # Assert 3: Verify that the 'ref' sample has no alt alleles
            ref_sample_idx = [i for i, s in enumerate(sample_names) if 'ref' in s][0]
            sample_column_idx = 9 + ref_sample_idx  # 9 is the FORMAT column index
            
            for i, variant_line in enumerate(variant_lines):
                fields = variant_line.split('\t')
                ref_genotype_info = fields[sample_column_idx]
                
                # Parse the genotype info (format: GT:CR:CA)
                gt_parts = ref_genotype_info.split(':')
                if len(gt_parts) >= 3:
                    ca_values = gt_parts[2]  # CA (alternative allele count)
                    # CA should be all zeros or dots for ref sample
                    assert all(c in '0.,|' for c in ca_values), \
                        f"Variant {i+1}: ref sample has non-zero alt alleles: {ca_values}"
            
            print(f"✓ ref sample has no alt alleles in any variant")
            
    finally:
        # Clean up
        if os.path.exists(vcf_file):
            os.unlink(vcf_file)


def test_cli_gfa2vcf_and_consolidate_with_bgzf_vcf(tmp_path):
    """Test CLI can write BGZF VCF and use it as consolidate input."""
    runner = CliRunner()
    gfa_file = os.path.join(os.path.dirname(__file__), "data", "simple_nested.gfa")
    vcf_file = tmp_path / "simple_nested.vcf.gz"
    output_vcf = tmp_path / "sample2.hap0.vcf"

    result = runner.invoke(
        cli,
        ['gfa2vcf', gfa_file, str(vcf_file), '--chr-id', 'chr1', '--ref-name', 'ref']
    )

    assert result.exit_code == 0, result.output
    assert vcf_file.exists()
    with gzip.open(vcf_file, 'rt') as f:
        assert f.readline().strip() == '##fileformat=VCFv4.2'

    result = runner.invoke(
        cli,
        ['consolidate', str(vcf_file), 'sample2', '0', str(output_vcf)]
    )

    assert result.exit_code == 0, result.output
    assert output_vcf.exists()
    assert output_vcf.read_text().startswith('##fileformat=VCFv4.2')


def test_cli_simplify_writes_origin_metadata(tmp_path):
    """Test simplify output maps simplified segments to original GFA segment IDs."""
    runner = CliRunner()
    gfa_file = os.path.join(os.path.dirname(__file__), "data", "simple_nested.gfa")
    simplified_gfa = tmp_path / "simple_nested.simplified.gfa"

    result = runner.invoke(
        cli,
        ['simplify', gfa_file, str(simplified_gfa), '--ref-name', 'ref']
    )

    assert result.exit_code == 0, result.output
    lines = simplified_gfa.read_text().splitlines()
    assert lines[0].split('\t') == ['H', 'VN:Z:1.0', f'og:Z:{gfa_file}']

    segment_origin_lists = []
    origins_by_segment = {}
    for line in lines:
        if not line.startswith('S\t'):
            continue
        fields = line.split('\t')
        origin_fields = [field for field in fields[3:] if field.startswith('oi:J:')]
        assert len(origin_fields) == 1
        original_ids = json.loads(origin_fields[0][5:])
        assert all(isinstance(original_id, str) for original_id in original_ids)
        assert all(not original_id.endswith(('_+', '_-')) for original_id in original_ids)
        assert all('terminus' not in original_id for original_id in original_ids)
        origins_by_segment[fields[1]] = original_ids
        segment_origin_lists.append(original_ids)

    assert segment_origin_lists
    assert any(len(original_ids) > 1 for original_ids in segment_origin_lists)
    assert any(len(original_ids) == 0 for original_ids in segment_origin_lists)

    reread = PangenomeGraph.from_gfa_line_by_line(str(simplified_gfa), ref_name='ref')
    assert reread.number_of_nodes() > 0
    for segment_id, original_ids in origins_by_segment.items():
        assert reread.nodes[f'{segment_id}_+']['original_ids'] == original_ids
        assert reread.nodes[f'{segment_id}_-']['original_ids'] == list(reversed(original_ids))

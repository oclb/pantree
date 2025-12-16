from pantree.cli import cli
from click.testing import CliRunner
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


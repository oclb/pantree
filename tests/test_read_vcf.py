import polars as pl
from pantree.read_vcf import read_vcf_to_lazyframe


def test_read_vcf_to_lazyframe():
    """Test that read_vcf_to_lazyframe returns a LazyFrame with parsed fields."""
    vcf_path = 'tests/data/tmp5.vcf'
    
    # Get LazyFrame
    lf = read_vcf_to_lazyframe(vcf_path)
    
    # Check it's a LazyFrame
    assert isinstance(lf, pl.LazyFrame)
    
    # Collect and check results
    df = lf.collect()
    
    # Check shape
    assert df.shape[0] == 9  # 9 rows
    assert df.shape[1] == 167  # All columns including parsed fields
    
    # Check standard VCF columns exist
    assert '#CHROM' in df.columns
    assert 'POS' in df.columns
    assert 'ID' in df.columns
    assert 'REF' in df.columns
    assert 'ALT' in df.columns
    
    # Check INFO fields were parsed
    assert 'NR' in df.columns
    assert 'VT' in df.columns
    assert 'DR' in df.columns
    assert 'RC' in df.columns
    assert 'AC' in df.columns
    assert 'AN' in df.columns
    assert 'HP' in df.columns
    assert 'TR_MOTIF' in df.columns
    assert 'NIA' in df.columns
    
    # Check genotype fields were parsed for samples
    assert 'CHM13_GT' in df.columns
    assert 'CHM13_CR' in df.columns
    assert 'CHM13_CA' in df.columns
    assert 'GRCh38_GT' in df.columns
    assert 'NA21309_GT' in df.columns
    
    # Check some data values
    assert df['#CHROM'][0] == 'chr0'
    assert df['POS'][0] == 9
    assert df['VT'][0] == 'DUP'

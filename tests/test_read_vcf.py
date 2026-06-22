import polars as pl
import gzip
from Bio import bgzf
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


def test_read_vcf_to_lazyframe_bgzf_matches_uncompressed(tmp_path):
    """Test compressed BGZF input parses the same as uncompressed input."""
    vcf_path = 'tests/data/tmp5.vcf'
    compressed_path = tmp_path / 'tmp5.vcf.gz'

    with open(vcf_path, 'r') as source, bgzf.open(compressed_path, 'wt') as target:
        target.write(source.read())

    uncompressed = read_vcf_to_lazyframe(vcf_path).collect()
    compressed = read_vcf_to_lazyframe(str(compressed_path)).collect()

    assert isinstance(read_vcf_to_lazyframe(str(compressed_path)), pl.LazyFrame)
    assert compressed.equals(uncompressed)


def test_read_vcf_to_lazyframe_gzip_matches_uncompressed(tmp_path):
    """Test legacy gzip-compressed input parses the same as uncompressed input."""
    vcf_path = 'tests/data/tmp5.vcf'
    compressed_path = tmp_path / 'tmp5.legacy.vcf.gz'

    with open(vcf_path, 'rt') as source, gzip.open(compressed_path, 'wt') as target:
        target.write(source.read())

    uncompressed = read_vcf_to_lazyframe(vcf_path).collect()
    compressed = read_vcf_to_lazyframe(str(compressed_path)).collect()

    assert isinstance(read_vcf_to_lazyframe(str(compressed_path)), pl.LazyFrame)
    assert compressed.equals(uncompressed)

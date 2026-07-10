import polars as pl
import gzip
from Bio import bgzf
from pantree.graph import PangenomeGraph
from pantree.read_vcf import read_vcf_to_lazyframe


def write_current_vcf(tmp_path):
    gfa_path = 'tests/data/simple_nested.gfa'
    vcf_path = tmp_path / 'simple_nested.current.vcf'
    graph = PangenomeGraph.from_gfa_line_by_line(gfa_path, ref_name='ref')
    graph.write_vcf(gfa_path, str(vcf_path), chr_name='chr0')
    return vcf_path


def write_current_no_genotype_vcf(tmp_path):
    gfa_path = 'tests/data/simple_nested.gfa'
    vcf_path = tmp_path / 'simple_nested.no_genotypes.vcf'
    graph = PangenomeGraph.from_gfa_line_by_line(gfa_path, ref_name='ref')
    graph.write_vcf(None, str(vcf_path), chr_name='chr0')
    return vcf_path


def test_read_vcf_to_lazyframe(tmp_path):
    """Test that read_vcf_to_lazyframe returns a LazyFrame with parsed fields."""
    vcf_path = write_current_vcf(tmp_path)
    
    # Get LazyFrame
    lf = read_vcf_to_lazyframe(str(vcf_path))
    
    # Check it's a LazyFrame
    assert isinstance(lf, pl.LazyFrame)
    
    # Collect and check results
    df = lf.collect()
    
    # Check shape
    assert df.shape[0] == 4
    
    # Check standard VCF columns exist
    assert '#CHROM' in df.columns
    assert 'POS' in df.columns
    assert 'ID' in df.columns
    assert 'REF' in df.columns
    assert 'ALT' in df.columns
    
    # Check INFO fields were parsed
    assert 'NR' in df.columns
    assert 'VT' in df.columns
    assert 'TP' in df.columns
    assert 'RC' in df.columns
    assert 'AC' in df.columns
    assert 'AN' in df.columns
    assert 'HP' in df.columns
    assert 'TR_MOTIF' in df.columns
    assert 'NIA' in df.columns
    assert 'UIDX' in df.columns
    
    # Check genotype fields were parsed for samples
    assert 'ref_GT' in df.columns
    assert 'ref_CR' in df.columns
    assert 'ref_CA' in df.columns
    assert 'sample1_GT' in df.columns
    assert 'sample2_GT' in df.columns
    
    # Check some data values
    assert df['#CHROM'][0] == 'chr0'
    assert df['POS'][0] > 0
    assert df['VT'][0] in {'SNP', 'MNP', 'INS', 'DEL', 'DUP', 'INV', 'COMPLEX'}
    assert df['UIDX'][0].isdigit()


def test_read_vcf_to_lazyframe_without_format_column(tmp_path):
    """Test no-genotype VCFs parse without FORMAT or sample-derived columns."""
    vcf_path = write_current_no_genotype_vcf(tmp_path)

    df = read_vcf_to_lazyframe(str(vcf_path)).collect()

    assert df.shape[0] == 4
    assert '#CHROM' in df.columns
    assert 'INFO' not in df.columns
    assert 'FORMAT' not in df.columns
    assert 'VT' in df.columns
    assert 'TP' in df.columns
    assert 'UIDX' in df.columns
    assert not any(column.endswith('_GT') for column in df.columns)
    assert not any(column.endswith('_CR') for column in df.columns)
    assert not any(column.endswith('_CA') for column in df.columns)


def test_read_vcf_to_lazyframe_bgzf_matches_uncompressed(tmp_path):
    """Test compressed BGZF input parses the same as uncompressed input."""
    vcf_path = write_current_vcf(tmp_path)
    compressed_path = tmp_path / 'simple_nested.current.vcf.gz'

    with open(vcf_path, 'r') as source, bgzf.open(compressed_path, 'wt') as target:
        target.write(source.read())

    uncompressed = read_vcf_to_lazyframe(str(vcf_path)).collect()
    compressed = read_vcf_to_lazyframe(str(compressed_path)).collect()

    assert isinstance(read_vcf_to_lazyframe(str(compressed_path)), pl.LazyFrame)
    assert compressed.equals(uncompressed)


def test_read_vcf_to_lazyframe_gzip_matches_uncompressed(tmp_path):
    """Test legacy gzip-compressed input parses the same as uncompressed input."""
    vcf_path = write_current_vcf(tmp_path)
    compressed_path = tmp_path / 'tmp5.legacy.vcf.gz'

    with open(vcf_path, 'rt') as source, gzip.open(compressed_path, 'wt') as target:
        target.write(source.read())

    uncompressed = read_vcf_to_lazyframe(str(vcf_path)).collect()
    compressed = read_vcf_to_lazyframe(str(compressed_path)).collect()

    assert isinstance(read_vcf_to_lazyframe(str(compressed_path)), pl.LazyFrame)
    assert compressed.equals(uncompressed)

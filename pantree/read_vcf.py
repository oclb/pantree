import polars as pl
import re

from .vcf_io import open_vcf_text_for_read

def read_vcf_to_lazyframe(vcf_path: str) -> pl.LazyFrame:
    # Parse header to extract INFO field IDs and FORMAT fields
    info_fields = []
    format_fields = []
    sample_names = []
    with open_vcf_text_for_read(vcf_path) as f:
        for line in f:
            if line.startswith('##INFO='):
                match = re.search(r'ID=([^,]+)', line)
                if match:
                    info_fields.append(match.group(1))
            elif line.startswith('##FORMAT='):
                match = re.search(r'ID=([^,]+)', line)
                if match:
                    format_fields.append(match.group(1))
            elif line.startswith('#CHROM'):
                # Extract sample names from header line
                parts = line.strip().split('\t')
                if 'FORMAT' in parts:
                    format_idx = parts.index('FORMAT')
                    sample_names = parts[format_idx + 1:]
                break

    # Read VCF data as LazyFrame
    result = pl.scan_csv(vcf_path, separator='\t', has_header=True, comment_prefix='##')

    # Build regex pattern dynamically from INFO fields
    pattern_parts = [f'{field}=([^;]+)' for field in info_fields]
    pattern = ';'.join(pattern_parts)

    # Parse INFO field into separate columns using generated regex
    result = result.with_columns(
        pl.col('INFO').str.extract_groups(pattern).struct.rename_fields(info_fields)
    ).unnest('INFO')

    # Parse genotype fields for each sample
    if format_fields:
        format_pattern = ':'.join([r'([^:]+)' for _ in format_fields])
        for sample in sample_names:
            result = result.with_columns(
                pl.col(sample).str.extract_groups(format_pattern).struct.rename_fields(
                    [f'{sample}_{field}' for field in format_fields]
                )
            ).unnest(sample)

    return result

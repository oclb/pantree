# `pantree consolidate`

Collapse overlapping and nested records for one sample haplotype:

```bash
uv run pantree consolidate input.vcf.gz HG002 0 HG002.hap0.vcf
```

Arguments are the input VCF, sample name, zero-based haplotype number, and output VCF. Input may be uncompressed VCF or compressed `.vcf.gz`.

Use this for a direct representation of one haplotype relative to the linear reference, rather than Pantree's graph-edge variant decomposition.

Run `uv run pantree consolidate --help` for the authoritative argument list.

Ordinary variants are consolidated separately within each chromosome. If the
selected haplotype carries a `TRANS_CHROM` variant, consolidation raises an error
before writing output: cross-chromosome breakends cannot be combined into a
within-chromosome REF/ALT sequence.

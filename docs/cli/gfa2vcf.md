# `pantree gfa2vcf`

Convert a pangenome graph to VCF:

```bash
uv run pantree gfa2vcf input.gfa output.vcf.gz \
  --chr-id chr21 \
  --ref-name GRCh38
```

An output ending in `.vcf.gz` is BGZF-compressed.

## Tree construction

`--dfs-method max_weight` is the default. To favor contiguous traversal along selected samples:

```bash
uv run pantree gfa2vcf input.gfa output.vcf \
  --dfs-method contiguous \
  --priority-samples HG002,CHM13
```

Priority samples are ordered and receive haplotype-relative position annotations.

## Genotype controls

- `--no-genotypes` writes sites without computing sample genotypes.
- `--no-missingness` retains genotypes but disables reference-span missingness inference. This is useful for some `vg`-generated graphs.

Use `--verbose` for console progress or `--log-path FILE` for a persistent log. Run `uv run pantree gfa2vcf --help` for the authoritative option list.

# `pantree gfa2vcf`

Convert a pangenome graph to VCF:

```bash
uv run pantree gfa2vcf input.gfa output.vcf.gz \
  --chr-id chr21 \
  --ref-name GRCh38
```

An output ending in `.vcf.gz` is BGZF-compressed.

For a multi-chromosome graph, omit `--chr-id`:

```bash
uv run pantree gfa2vcf whole-genome.gfa.gz whole-genome.vcf.gz --ref-name GRCh38
```

Chromosome names are inferred from the selected reference sample's `W` sequence
names, or PanSN `P` path names (`sample#haplotype#chromosome`). Sample contig names
need not match the reference chromosome names. `--chr-id` overrides the output
label only for a single-chromosome input; supplying it for multiple chromosomes
raises an error. The previous default label `chr0` is replaced by the inferred
reference sequence name.

Use a reference with one intact walk per chromosome. Multiple reference walks for
the same chromosome are rejected; sample haplotypes may have multiple fragments.
Whole-genome conversion keeps the complete graph in memory, so separate runs
remain an option when memory is limited.

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

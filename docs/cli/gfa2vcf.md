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

Pantree selects the first matching reference haplotype in input order and uses
all its chromosome walks. Additional reference haplotypes produce a warning and
remain available for genotyping; their assembly-contig names do not create extra
reference chromosomes. Supply one intact walk per chromosome in the selected
haplotype. Duplicate or fragmented reference walks with the same chromosome and
haplotype are rejected, including in alternate reference haplotypes. Other sample
haplotypes may have multiple fragments. Reference coordinates include the selected
`W` start offset; `P` coordinates start at zero.
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

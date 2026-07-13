# Pantree documentation

Pantree identifies variants in a pangenome graph by constructing a reference tree and representing non-tree edges as variants.

## Guides

- [Concepts](concepts.md): graph representation, reference trees, variants, genotypes, and consolidation.
- [`gfa2vcf`](cli/gfa2vcf.md): convert a GFA graph to VCF.
- [`simplify`](cli/simplify.md): remove small variation and contract paths while retaining provenance.
- [`consolidate`](cli/consolidate.md): produce one haplotype-versus-reference VCF.
- [File formats](file-formats.md): supported GFA records and Pantree output contracts.

Install with `uv sync`, inspect the commands with `uv run pantree --help`, and run tests with `uv run pytest`.

# Pantree

Converts pangenome graphs (GFA format) into VCF files. Variants are defined as edges outside a reference tree constructed via DFS on the pangenome graph.

**Paper**: Nowbandegani PS, Zhang S, Hu H, Li H, O'Connor LJ. "Defining and cataloging variants in pangenome graphs." *bioRxiv* 2025. https://doi.org/10.1101/2025.08.04.668502

## Tech Stack

- Python 3.10+ (target: 3.13)
- Package manager: `uv`
- Build system: hatchling / pyproject.toml
- Key libs: networkx, numpy, pandas, polars, click, scipy

## Commands

```bash
# Setup
uv venv && uv sync

# Run
uv run pantree gfa2vcf <gfa_file> <vcf_file> [options]
uv run pantree consolidate <vcf_file> <sample> <haplotype> <output>

# Test
uv run pytest
uv run pytest --cov=pantree
```

## Project Structure

```
pantree/              # Main package
  graph.py            # Core PangenomeGraph class (extends nx.DiGraph)
  cli.py              # Click CLI entry point
  dfs.py              # DFS tree construction (max_weight, haplo_contiguous)
  vcf.py              # VCF file writing
  gfa.py              # GFA file parsing
  genotype.py         # Genotype counting per walk
  walk_variants.py    # Haplotype variant consolidation
  evaluating_functions.py  # VCF INFO field definitions
  utils.py            # Node/edge complement helpers
  logger.py           # Custom logging with memory tracking
tests/                # pytest suite with GFA/VCF fixtures in tests/data/
```

## Key Concepts

- **Bidirected graph**: Nodes are `id_direction` (e.g., "1_+", "2_-"). Complement flips both ID and direction.
- **Terminal nodes**: `+_terminus_+`, `-_terminus_-`
- **Reference tree**: Built via DFS from reference path. Edges NOT in tree = variants.
- **Variant types**: SNP, INS, DEL, INV, DUP, MNP, COMPLEX
- **Walks**: Sample haplotypes as node paths through the graph
- **Edge attributes**: `index`, `weight`, `is_in_tree`, `branch_point`, `is_back_edge`

## Data Flow

GFA → parse segments/links/walks → build reference tree (DFS) → identify variant edges → compute genotypes per walk → determine missingness → extract REF/ALT alleles → write VCF

## Notebook

This project uses a separate notebook repository for analysis logs. See `notebook/INDEX.md` for a summary of past work.

# Pantree Agent Context

Pantree converts pangenome graphs in GFA format into VCF records. It builds a reference tree over a bidirected graph and treats graph edges outside that tree as variants.

## Project map

- `pantree/graph.py`: central `PangenomeGraph` type and graph-to-variant orchestration.
- `pantree/gfa.py`: streaming parser for GFA `S`, `L`, `P`, and `W` records.
- `pantree/dfs.py`: `max_weight` and `contiguous` reference-tree construction.
- `pantree/genotype.py` and `pantree/genotype_transformation.py`: allele counting, missingness, and genotype transformations.
- `pantree/vcf.py`, `pantree/vcf_io.py`, and `pantree/read_vcf.py`: VCF construction and compressed/uncompressed I/O.
- `pantree/walk_variants.py`: collapse nested or overlapping variants for one haplotype.
- `pantree/simplify_graph.py` and `pantree/write_gfa.py`: simplify a graph and retain source-segment provenance.
- `pantree/cli.py`: `pantree gfa2vcf`, `pantree consolidate`, and `pantree simplify`.
- `tests/`: pytest suite; small GFA fixtures live in `tests/data/`.
- `notebook/`: separate git repository containing durable project history and todos.

## Commands

```bash
uv sync
uv run pantree --help
uv run pytest
uv run pytest --cov=pantree
```

Run focused tests while developing, then the full suite before handing off a behavioral change.

## Agent skills

- Use `pantree-setup` for installation, dependency, or environment validation work.
- Use `pantree-usage` for CLI/API workflows, scientific concepts, and file-format contracts.

## Terminology

- **Bidirected graph:** each input segment has forward and reverse nodes, represented internally as `<segment>_+` and `<segment>_-`.
- **Complement:** the reverse-orientation counterpart of a node or edge.
- **Reference tree:** DFS tree used as the baseline representation of the graph.
- **Variant edge:** an edge outside the reference tree.
- **Walk:** an input sample haplotype represented as a path through graph nodes.
- **Priority sample:** a sample whose haplotypes influence contiguous DFS traversal and haplotype-position annotations.
- **Consolidation:** collapsing nested or overlapping variants into one haplotype-versus-reference representation.
- **Origin metadata:** simplified-GFA tags linking output segments to the input GFA and its segment IDs.

## Gotchas

- Preserve both orientations and complement edges when changing graph topology.
- GFA `P` paths and `W` walks use different syntax and naming conventions.
- Missingness inference may be inappropriate for graphs produced by `vg`; use `--no-missingness` for this case.
- A `.vcf.gz` output must be BGZF-compatible and sorted for tabix use, not merely gzip-compressed.
- Simplified GFA output must preserve `og:Z` header provenance and `oi:J` segment provenance.
- Treat zero-byte Dropbox files as potentially unsynced; ask the user to sync rather than overwrite them.
- Check `git status` in both the main repository and `notebook/` before modifying or committing files.

---
name: pantree-usage
description: Run, debug, or explain Pantree analyses, CLI commands, Python APIs, scientific concepts, and GFA/VCF contracts. Use for gfa2vcf conversion, graph simplification, haplotype consolidation, reference-tree behavior, genotype interpretation, or Pantree file formats.
---

# Pantree usage

Start by inspecting the installed command surface:

```bash
uv run pantree --help
uv run pantree <command> --help
```

Select the relevant reference:

- Read `references/concepts.md` before explaining or changing reference-tree, variant, coordinate, genotype, or consolidation behavior.
- Read `references/gfa2vcf.md` for graph-to-VCF commands and option selection.
- Read `references/simplify.md` for simplification and origin-metadata behavior.
- Read `references/consolidate.md` for single-haplotype consolidation.
- Read `references/file-formats.md` before changing or relying on GFA, VCF, compression, or provenance contracts.

Prefer small fixtures in `tests/data/` for verification. Check that a fixture is present and non-empty. Preserve user data and never overwrite an input GFA or VCF with generated output.

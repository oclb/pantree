# `pantree simplify`

Simplify a graph while retaining source-segment provenance:

```bash
uv run pantree simplify input.gfa simplified.gfa \
  --ref-name GRCh38 \
  --min-allele-length 1000
```

`--min-allele-length` is the minimum combined allele length for retaining a variant. Use `--pos-start` and `--pos-end` to restrict simplification to a reference interval.

Output uses new sequential segment IDs. The `og:Z` header tag identifies the source GFA, and every segment has an `oi:J` JSON array of represented input segment IDs. Artificial terminal segments may have an empty array.

Run `uv run pantree simplify --help` for the authoritative option list.

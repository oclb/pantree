# File formats

## Input GFA

Pantree reads plain-text or gzip-compressed GFA and recognizes:

- `S` segment records, including optional `oi:J` origin metadata.
- `L` link records.
- GFA 1.0 `P` path records.
- GFA 1.1 `W` walk records.

Internally, oriented segment identifiers are `<segment>_+` and `<segment>_-`. `P` records use comma-separated segment orientations. `W` records use `>` and `<` traversal markers. Pantree constructs a `W` haplotype name from its sample, haplotype, and sequence fields.

For whole-genome inputs, chromosome names come from the selected reference's
`W` sequence field or the final field of a PanSN `P` name. For plain `P` names,
the entire name is the sample and sequence name. Segment IDs must be unique
across the graph, and the selected reference must have one intact walk per
chromosome. Reference-walk start offsets are retained.

## Output VCF

`gfa2vcf` writes `.vcf` or BGZF-compressed `.vcf.gz`. Compressed records are coordinate-sorted for tabix compatibility. Coordinates and alleles follow VCF anchoring rules, including a preceding base for otherwise empty indel alleles.

The INFO and FORMAT declarations written by the installed version are authoritative. Inspect a generated header or the VCF writer before building downstream parsing around optional annotations.

`TRANS_CHROM` is an optional `Number=1,Type=String` INFO tag containing the
destination node's chromosome when it differs from the source chromosome in
`CHROM`. It is omitted on ordinary records. INFO tags may be absent or reordered;
Pantree's reader represents absent tags as null.

Cross-chromosome edges have no within-chromosome REF/ALT path pair. They are
exported as single breakends (`VT=BND;SVTYPE=BND`), anchored at the source's
reference-tree projection, with ALT `A.` or `.A` for a reference base `A` according
to source orientation. `TRANS_CHROM` retains the destination chromosome; no precise
partner position or reconstructed rearrangement is claimed. The graph edge ID
retains both oriented segment IDs. See the [VCF 4.2 single-breakend specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf).

## Simplified GFA

Simplified output is GFA 1.0 with `S` and `L` records. Links use `0M` overlap because simplification does not retain original overlap information.

```gfa
H	VN:Z:1.0	og:Z:input.gfa
S	s1	ACGT	oi:J:["12","13","14"]
S	s2	*	oi:J:[]
```

- `og:Z` records the input GFA path.
- `oi:J` is a JSON array of original segment IDs represented by an output segment.
- An empty `oi:J` denotes an artificial segment with no corresponding input segment.

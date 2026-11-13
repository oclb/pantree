# Concepts

## Bidirected graph

Pantree represents each GFA segment in both orientations. Segment `12` becomes `12_+` and `12_-`; reversing an edge reverses its direction and complements both endpoints. Graph transformations must maintain this symmetry.

A self-complementary link such as `12_+ -> 12_-` has one directed representation. An ordinary loop such as `12_+ -> 12_+` has a distinct complementary loop, `12_- -> 12_-`. Each link has one representative and one shared index; each walk traversal contributes one visit to that link, including when a walk visits both complementary orientations.

Artificial terminal nodes connect walk ends so tree construction can operate over a common graph. They are implementation helpers rather than biological segments.

## Reference tree and variants

Pantree constructs a directed DFS tree rooted along the requested reference. Tree edges define the baseline traversal; edges outside the tree define variants. The tree path between a variant edge's endpoints supplies the reference allele, while the non-tree traversal supplies the alternate allele.

The default `max_weight` DFS favors heavily used edges. The `contiguous` method attempts to remain on individual priority-sample haplotypes. Priority-sample order controls which haplotypes are preferred when a traversal switches.

## Chromosome assignment

A whole-genome graph has two terminal binodes per reference chromosome. Its
reference tree is stored as a forest with one root per chromosome. Both DFS
methods start an ordinary stack DFS at each chromosome root, prioritizing
existing reference-successor edges. Reference nodes and roots are reserved for
their owning chromosome before traversal. Each node and
its complement inherit one chromosome from their tree ancestry. An edge between
different chromosomes is therefore always a variant edge.

For a reported variant edge `(u,v)`, `CHROM` is the chromosome of `u`. If the
chromosome of `v` differs, the directed edge has a `trans_chrom` property containing
that single chromosome name. The VCF writer exports it as `TRANS_CHROM`. The
complement edge has the opposite source and destination chromosomes. Cross-chromosome
records retain the input link's representative orientation.

Chromosome assignment depends on the tree. Shared reference segments are owned
by the first reference chromosome in lexicographic order. Branches with no forward
route from a reference are attached through terminal edges, using a nearby
reference to choose the root. Components with no connection to any reference use
the first chromosome's root; their assigned coordinates are projections and do
not establish a biological chromosome location.

`reference_paths` and `chromosome_termini` expose the per-chromosome paths and
terminal pairs, retaining every reference step including shared segments and
same- or opposite-orientation revisits. Only artificial terminal links are added
to span regions unreachable within the owning chromosome; reference walks are
never shortened to bridge a shared segment. The legacy `reference_path` remains the first path. Use
`chromosome(node)` and `termini_for(node)` for chromosome-aware code. Low-level
`count_edge_visits` and graph simplification still require a single chromosome.

## Coordinates and genotypes

VCF coordinates derive from positions along the named reference walk. Pantree counts traversal of variant edges by each input walk to construct genotypes. Missingness is inferred from the portions of the reference spanned by a walk; disable this inference for graphs whose reference representation violates those assumptions.

A segment on the reference walk and its complement share the position at the end of its first occurrence. Every visit contributes to the cumulative offset, so later segments and the end terminus include the lengths of repeated segments. Reference coordinates include the selected GFA `W` start offset for both single- and multi-chromosome inputs; `P` paths start at zero. The end terminus is the start offset plus the full walk length. Priority-haplotype annotations also include `W` offsets and retain their existing last-occurrence rule.

When priority samples are supplied, Pantree can also record positions relative to those haplotypes. Indels follow the VCF requirement that REF and ALT be non-empty, so Pantree may prepend an anchoring base and decrement the position.

Whole-genome missingness is evaluated separately for each chromosome. Walk
fragments are combined within a sample haplotype, and a walk crossing chromosomes
contributes a separate coverage interval for each contiguous chromosome run.
Coverage on another chromosome cannot turn a missing genotype into reference.

Both DFS methods support cyclic reference walks, including inverted revisits.
First-occurrence coordinates describe a graph projection: a cyclic reference can
still yield nonmonotonic projections or linear alleles such as `REF=.` in existing
inversion cases. Loading and counting such walks does not guarantee a VCF suitable
for every linear-reference consumer. Inspect node IDs and alleles for that use;
`--no-missingness` disables coordinate-based inference. Consolidation retains its
existing linear-allele limitations for cyclic references.

## Consolidation

Graph-to-VCF conversion can emit nested variants. Consolidation selects one sample haplotype and collapses overlapping or nested records into a direct haplotype-versus-linear-reference representation.

## Simplification and provenance

Simplification removes variation below a requested combined allele length, trims optional reference-coordinate bounds, and contracts compatible paths. Output segment IDs are new, so output records the source GFA in `og:Z` and represented input segment IDs in `oi:J`.

Links traversed by the reference walk remain protected at every length threshold. A hairpin at a segment boundary has no intervening allele sequence, so a positive threshold can remove it when it is off the reference walk. Retained hairpins and loops remain explicit wherever contraction would delete an endpoint needed to preserve their topology.

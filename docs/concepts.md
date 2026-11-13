# Concepts

## Bidirected graph

Pantree represents each GFA segment in both orientations. Segment `12` becomes `12_+` and `12_-`; reversing an edge reverses its direction and complements both endpoints. Graph transformations must maintain this symmetry.

A self-complementary link such as `12_+ -> 12_-` has one directed representation. An ordinary loop such as `12_+ -> 12_+` has a distinct complementary loop, `12_- -> 12_-`. Each link has one representative and one shared index; each walk traversal contributes one visit to that link, including when a walk visits both complementary orientations.

Artificial terminal nodes connect walk ends so tree construction can operate over a common graph. They are implementation helpers rather than biological segments.

## Reference tree and variants

Pantree constructs a directed DFS tree rooted along the requested reference. Tree edges define the baseline traversal; edges outside the tree define variants. The tree path between a variant edge's endpoints supplies the reference allele, while the non-tree traversal supplies the alternate allele.

The default `max_weight` DFS favors heavily used edges. The `contiguous` method attempts to remain on individual priority-sample haplotypes. Priority-sample order controls which haplotypes are preferred when a traversal switches.

## Coordinates and genotypes

VCF coordinates derive from positions along the named reference walk. Pantree counts traversal of variant edges by each input walk to construct genotypes. Missingness is inferred from the portions of the reference spanned by a walk; disable this inference for graphs whose reference representation violates those assumptions.

A segment on the reference walk and its complement share the position at the end of its first occurrence. Every visit contributes to the cumulative offset, so later segments and the end terminus include the lengths of repeated segments. Reference coordinates start at zero for the loaded walk; GFA `W` start offsets apply to haplotype position annotations.

When priority samples are supplied, Pantree can also record positions relative to those haplotypes. Indels follow the VCF requirement that REF and ALT be non-empty, so Pantree may prepend an anchoring base and decrement the position.

## Consolidation

Graph-to-VCF conversion can emit nested variants. Consolidation selects one sample haplotype and collapses overlapping or nested records into a direct haplotype-versus-linear-reference representation.

## Simplification and provenance

Simplification removes variation below a requested combined allele length, trims optional reference-coordinate bounds, and contracts compatible paths. Output segment IDs are new, so output records the source GFA in `og:Z` and represented input segment IDs in `oi:J`.

Links traversed by the reference walk remain protected at every length threshold. A hairpin at a segment boundary has no intervening allele sequence, so a positive threshold can remove it when it is off the reference walk. Retained hairpins and loops remain explicit wherever contraction would delete an endpoint needed to preserve their topology.

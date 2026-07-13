# Concepts

## Bidirected graph

Pantree represents each GFA segment in both orientations. Segment `12` becomes `12_+` and `12_-`; reversing an edge reverses its direction and complements both endpoints. Graph transformations must maintain this symmetry.

Artificial terminal nodes connect walk ends so tree construction can operate over a common graph. They are implementation helpers rather than biological segments.

## Reference tree and variants

Pantree constructs a directed DFS tree rooted along the requested reference. Tree edges define the baseline traversal; edges outside the tree define variants. The tree path between a variant edge's endpoints supplies the reference allele, while the non-tree traversal supplies the alternate allele.

The default `max_weight` DFS favors heavily used edges. The `contiguous` method attempts to remain on individual priority-sample haplotypes. Priority-sample order controls which haplotypes are preferred when a traversal switches.

## Coordinates and genotypes

VCF coordinates derive from positions along the named reference walk. Pantree counts traversal of variant edges by each input walk to construct genotypes. Missingness is inferred from the portions of the reference spanned by a walk; disable this inference for graphs whose reference representation violates those assumptions.

When priority samples are supplied, Pantree can also record positions relative to those haplotypes. Indels follow the VCF requirement that REF and ALT be non-empty, so Pantree may prepend an anchoring base and decrement the position.

## Consolidation

Graph-to-VCF conversion can emit nested variants. Consolidation selects one sample haplotype and collapses overlapping or nested records into a direct haplotype-versus-linear-reference representation.

## Simplification and provenance

Simplification removes variation below a requested combined allele length, trims optional reference-coordinate bounds, and contracts compatible paths. Output segment IDs are new, so output records the source GFA in `og:Z` and represented input segment IDs in `oi:J`.

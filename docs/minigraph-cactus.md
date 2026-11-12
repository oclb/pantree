# Minigraph-Cactus inputs

Use a full or clipped base-level Minigraph-Cactus GFA with embedded sample walks.
The separate `.sv.gfa.gz` minigraph output lacks those walks. Select the primary
reference with `--ref-name`; the ordinary MC pipeline preserves its paths through
clipping. `--collapse` can introduce reference cycles and is outside this contract.

Chromosome names come from reference `W` sequence names. Sample assembly contigs
can have different names, and clipped haplotypes can contain multiple `W` fragments
with the same sample, haplotype, and sequence name but different start coordinates.
Pantree aggregates these fragments for genotyping.

The full graph retains material removed by clipping. Choose full when that
sequence matters to the analysis; clipped input limits analysis to the retained
graph. Neither input choice restores variants already removed upstream.

## Compatibility evidence

Checked against Cactus commit
[`4c4d95b321d1c62fc9ad87207e2e8d2f8f3de065`](https://github.com/ComparativeGenomicsToolkit/cactus/tree/4c4d95b321d1c62fc9ad87207e2e8d2f8f3de065):

1. The [output specification](https://github.com/ComparativeGenomicsToolkit/cactus/blob/4c4d95b321d1c62fc9ad87207e2e8d2f8f3de065/doc/pangenome.md#output)
   specifies GFA 1.1 walks; its examples distinguish sample, haplotype, sequence,
   and offset fields.
2. [`join_vg` and `make_vg_indexes`](https://github.com/ComparativeGenomicsToolkit/cactus/blob/4c4d95b321d1c62fc9ad87207e2e8d2f8f3de065/src/cactus/refmap/cactus_graphmap_join.py)
   make node IDs unique across chromosomes with `vg ids -j`, then concatenate
   chromosome GFA records under a shared header.
3. The [clipping documentation](https://github.com/ComparativeGenomicsToolkit/cactus/blob/4c4d95b321d1c62fc9ad87207e2e8d2f8f3de065/doc/pangenome.md#clipping-filtering-and-indexing)
   describes full and clipped graphs and the sequence removed by clipping.

`tests/data/multi_chrom.gfa` and `tests/test_chromosomes.py` exercise this output
contract, including different assembly contig names and fragmented sample walks.
The compatibility check is against upstream code/documentation and synthetic
fixtures; it is not a real-data whole-genome benchmark.

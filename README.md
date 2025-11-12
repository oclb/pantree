# Identifying genetic variants in the pangenome using a reference tree

## Table of Contents
- [Introduction](#introduction)
- [Installation](#installation)
- [Command Line Interface](#command-line-interface)
- [Usage](#usage)

## Introduction
`pantree` converts a pangenome graph `.gfa` file into a `.vcf` file containing variants identified in the graph. It creates a reference tree and defines variants as edges that are not in the reference tree. For more information, please see our [preprint](https://www.biorxiv.org/content/10.1101/2025.08.04.668502v1).

`pantree` is a work in progress. If you use our code, please do reach out with questions and feedback.

## Installation

You can install `pantree` using `uv`:

```bash
git clone https://github.com/oclb/pantree.git
cd pantree
uv venv
uv sync
```

## Command Line Interface

```bash
pantree <gfa_file> <vcf_file> [options]
```

### Required Arguments
- `gfa_file`: Path to the input GFA file containing the pangenome graph
- `vcf_file`: Output path for VCF file

### Optional Arguments
- `--chr-id TEXT`: Chromosome ID for VCF output (default: "chr0")
- `--ref-name TEXT`: Reference sample name (default: "GRCh38")
- `--no-genotypes`: Skip genotype computation
- `--log-path TEXT`: Path to log file for tracking progress and memory usage
- `--verbose, -v`: Enable verbose logging to console
- `--dfs-method [max_weight|contiguous]`: DFS method for reference tree construction (default: "max_weight")
  - `max_weight`: Prioritize edges with higher weights (more walks)
  - `contiguous`: Prioritize contiguous haplotype paths
- `--priority-samples TEXT`: Comma-separated list of sample names. For haplotypes belonging to samples in this list, haplotype positions are computed (for variant edges whose branch point belongs to the haplotype). With the 'contiguous' DFS method, these haplotypes are prioritized when building the DFS tree, in the order they are specified. 

### Example Usage
```bash
pantree input.gfa output.vcf

pantree input.gfa output.vcf \
  --chr-id chr20 \
  --ref-name GRCh38 \
  --log-path analysis.log \
  --verbose \
  --dfs-method contiguous \
  --priority-samples "GRCh38,CHM13,HG002"
```

## Python API Usage
```python
from graph_var import PangenomeGraph, Genotype
from graph_var.logging import setup_logger

gfa_path = "/path/to/graph.gfa"

logger = setup_logger(log_path="analysis.log", verbose=True)
G: PangenomeGraph = PangenomeGraph.from_gfa(
    gfa_path, 
    ref_name="GRCh38",
    logger=logger,
    dfs_method_name="max_weight",
    priority_dict={"GRCh38": 0, "CHM13": 1, "HG002": 2}
)

# Also return walks; causes increased memory requirements
walks: list[list[str]]
G, walks = PangenomeGraph.from_gfa(gfa_path, return_walks=True)

# Get the genotype of some walk
genotype: Genotype = G.genotype(walks[0])

# Generate VCF file with genotypes
vcf_path = "/path/to/output.vcf"
chr_id = "chr1"
G.write_vcf(gfa_path, vcf_path, chr_id)

# Generate VCF without genotypes
G.write_vcf(None, vcf_path, chr_id)

# Get genotypes from GFA
sample_to_genotype = G.genotypes_from_gfa(gfa_path)

# Access graph properties
print(f"Number of nodes: {G.number_of_nodes()}")
print(f"Number of variant edges: {len(G.variant_edges)}")
print(f"Reference path length: {len(G.reference_path)}")
```

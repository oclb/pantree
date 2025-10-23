# Identifying genetic variants in the pangenome using a reference tree

## Table of Contents
- [Introduction](#introduction)
- [Installation](#installation)
- [Command Line Interface](#command-line-interface)
- [Usage](#usage)

## Introduction
`pantree` converts a pangenome graph `.gfa` file into a `.vcf` file containing variants identified in the graph. It creates a reference tree and defines variants as edges that are not in the reference tree. For more information, please see our [preprint](https://www.biorxiv.org/content/10.1101/2025.08.04.668502v1).

`pantree` supports both **GFA1.0** (P-lines) and **GFA1.1** (W-lines) formats, making it compatible with graphs from various tools including PGGB and Minigraph-Cactus.

`pantree` is a work in progress. If you use our code, please do reach out with questions and feedback.

## Installation

You can install `pantree` using `uv`:

```bash
git clone https://github.com/oclb/graph_var.git
cd graph_var
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
- `--chr-id`: Chromosome ID for VCF file (default: "chr0")
- `--ref-name`: Reference name to use (default: "GRCh38")
- `--no-genotypes`: Skip writing genotype information to VCF

### Example Usage
```bash
# Basic usage - analyze a GFA file and output VCF
pantree input.gfa output.vcf

# Specify chromosome ID
pantree input.gfa output.vcf --chr-id chr20

# Use a different reference name
pantree input.gfa output.vcf --ref-name hg19

# Generate VCF without genotype information
pantree input.gfa output.vcf --no-genotypes

```

## Usage
```python
from graph_var import PangenomeGraph

# Read a .gfa file
gfa_path = "/path/to/graph.gfa"
reference_path_index = 1
G, walks, walk_sample_names = PangenomeGraph.from_gfa(gfa_path, 
                                                return_walks=True, compressed=False, 
                                                reference_path_index=reference_path_index)

# Generate vcf file
vcf_path = "/path/to/output.vcf"
chr_id = "chr1"
G.write_vcf(gfa_path, vcf_path, chr_id)

# Enumerate variants of different types
edge_type_count: dict = G.variant_edges_summary()

# Get the genotype of a walk, then reconstruct edge visit counts
genotype: dict = G.genotype(walks[0])
edge_visit_counts: dict = G.count_edge_visits(genotype)

```

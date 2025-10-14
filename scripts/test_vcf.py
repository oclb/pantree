from graph_var import PangenomeGraph
import matplotlib.pyplot as plt
import networkx as nx
import os

def main():
    # Read a .gfa file
    gfa_path = "/Users/lukeoconnor/Downloads/grch38_chr20.gfa"
    gfa_path = "data/c4a_with_inversion_and_sequences.gfa"
    G, walks, walk_sample_names = PangenomeGraph.from_gfa(gfa_path,
                                                    return_walks=True, compressed=False)
    G.write_vcf(gfa_path, "output/test.vcf", chr_name='chr6', check_degenerate=True)
if __name__ == "__main__":
    main()
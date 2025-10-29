import re
import gzip
import pickle

import time
import os
import psutil
from datetime import datetime
from dataclasses import dataclass

def sequence_complement(s: str) -> str:
    def _base_complement(letter: str) -> str:
        base_dict: dict = {'A': 'T', 'C': 'G', 'T' : 'A', 'G': 'C'}
        return base_dict[letter] if letter in base_dict else letter

    return ''.join(_base_complement(c) for c in reversed(s))

def node_complement(s: str):
    return s[:-1] + _flip(s[-1])

def edge_complement(e: tuple[str, str]):
    return node_complement(e[1]), node_complement(e[0])

def walk_complement(w: list[str]) -> list[str]:
    return [node_complement(v) for v in w[::-1]]

def _flip(s):
    if s == '+':
        return '-'
    elif s == '-':
        return '+'
    else:
        raise ValueError()

def _node_recover(node_id):
    if len(node_id.split('_')) == 2:
        node, direction = node_id.split('_')
    else:
        terminus, node, direction = node_id.split('_')
        node = 'start' + node.title() if terminus == '+' else 'end' + node.title()
    if direction == '+':
        return '>'+node
    elif direction == '-':
        return '<'+node
    else:
        raise ValueError()

def _node_convert(node_id):
    direction, node = node_id[0], node_id[1:]
    if direction == '>':
        return node+"_+"
    elif direction == '<':
        return node+"_-"
    else:
        raise ValueError()

def read_gfa(filename, ref_name='GRCh38'):
  nodes = []
  sequences = []
  edges = []
  walks = []
  walk_sample_names = []
  pattern = r'\w+|[<>]'
  reference_index = 0
  hit_reference = False

  data_dict = {}

  if filename.endswith('.gz'):
      file = gzip.open(filename, 'rt')
  else:
      file = open(filename, 'r')
  for line in file:
      parts = line.strip().split('\t')
      if parts[0] == 'S':
          nodes.append(parts[1])
          sequences.append(parts[2])
      elif parts[0] == 'L':
          edge = (parts[1], parts[3], parts[2], parts[4])
          edges.append(edge)
      elif parts[0] == 'W':
          if not hit_reference:
              if parts[1] == ref_name:
                  hit_reference = True
              else:
                  reference_index += 1
          sample_name = parts[1]+'_'+parts[2]
          p = parts[6]
          matches = re.findall(pattern, p)
          # List to store node IDs
          node_ids = []

          # Iterate through the matches to generate node IDs
          for i, match in enumerate(matches):
              if match not in ['<', '>']:  # If the match is a word
                  if matches[i - 1] == '<':
                      # For '<', generate IDs with '-'
                      node_ids.append(f'{match}_-')
                  elif matches[i - 1] == '>':
                      # For '>', generate IDs with '+'
                      node_ids.append(f'{match}_+')
          #node_ids.append('end_node')
          walks.append(node_ids)
          walk_sample_names.append(sample_name)

  data_dict['nodes'] = nodes
  data_dict['edges'] = edges
  data_dict['walks'] = walks
  data_dict['walk_sample_names'] = walk_sample_names
  data_dict['sequences'] = sequences
  data_dict['reference_index'] = reference_index

  return data_dict

@dataclass
class GFAWalkLine:
    hap_name: str
    contig_start: int
    walk: list[str]

    @staticmethod
    def _get_hap_name(sample_name: str, hap_index: int) -> str:
        if int(hap_index) == 0:
            return sample_name
        
        return sample_name + '_' + str(hap_index)

    @classmethod
    def from_parts_P(cls, parts: list[str]):
        hap_name = parts[1]
        contig_start = 0
        walk = cls._parse_P_walk(parts[2])
        return cls(hap_name, contig_start, walk)

    @classmethod
    def from_parts_W(cls, parts: list[str]):
        hap_name = cls._get_hap_name(parts[1], parts[2])
        contig_start = int(parts[4]) if parts[4] != '*' else 0
        walk = cls._parse_W_walk(parts[6])
        return cls(hap_name, contig_start, walk)

    @staticmethod
    def _parse_P_walk(segment_names: str) -> list[str]:
        node_ids = []
        for segment in segment_names.split(','):
            segment = segment.strip()
            if not segment:
                continue
            # Last character is orientation (+/-)
            node = segment[:-1]
            orientation = segment[-1]
            if orientation == '+':
                node_ids.append(f'{node}_+')
            elif orientation == '-':
                node_ids.append(f'{node}_-')
        return node_ids

    @staticmethod
    def _parse_W_walk(walk: str) -> list[str]:
        pattern = r'\w+|[<>]'
        matches = re.findall(pattern, walk)
        node_ids = []
        for i, match in enumerate(matches):
            if match not in ['<', '>']:  # If the match is a word
                if matches[i - 1] == '<':
                    # For '<', generate IDs with '-'
                    node_ids.append(f'{match}_-')
                elif matches[i - 1] == '>':
                    # For '>', generate IDs with '+'
                    node_ids.append(f'{match}_+')
        return node_ids

@dataclass
class GFANodeLine:
    node_id: str
    sequence: str

    @classmethod
    def from_parts(cls, parts: list[str]):
        return cls(parts[1], parts[2])

@dataclass
class GFAEdgeLine:
    u: str
    v: str
    
    @classmethod
    def from_parts(cls, parts: list[str]):
        u = parts[1] + '_' + parts[2]
        v = parts[3] + '_' + parts[4]
        return cls(u, v)

line_getters = {
    'S': GFANodeLine.from_parts,
    'L': GFAEdgeLine.from_parts,
    'W': GFAWalkLine.from_parts_W,
    'P': GFAWalkLine.from_parts_P,
}

def read_gfa_line_by_line(filename: str):
    if filename.endswith('.gz'):
        opener = gzip.open
    else:
        opener = open

    with opener(filename, 'rt') as file:
        for line in file:
            parts = line.strip().split('\t')
            # Skip empty lines and header lines
            if not parts[0] or parts[0] == 'H' or parts[0] == 'C':
                continue
            line_getter = line_getters[parts[0]]
            yield line_getter(parts)

def save_graph_to_pkl(G, path, compressed=False):
    if compressed:
        with gzip.open(path, 'wb') as file:
            pickle.dump(G, file, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        with open(path, 'wb') as file:
            pickle.dump(G, file, protocol=pickle.HIGHEST_PROTOCOL)

def load_graph_from_pkl(path, compressed=False):
    if compressed:
        with gzip.open(path, 'rb') as file:
            G = pickle.load(file)
    else:
        with open(path, 'rb') as file:
            G = pickle.load(file)
    return G

def merge_dicts(dicts: list[dict]):
    merged = {}
    for d in dicts:
        for k, v in d.items():
            if k in merged:
                merged[k] += v
            else:
                merged[k] = v
    return merged

def group_walks_by_name(walks: list, names: list) -> dict:
    """
    Group walks by their names into sublists.
    
    Args:
        walks: List of walks (each walk is a list of nodes)
        names: List of names corresponding to each walk
    
    Returns:
        List of lists, where each sublist contains all walks with the same name
    """
    # Create a dictionary to store walks by name
    from collections import defaultdict
    walks_by_name = defaultdict(list)
    
    # Group walks by their names
    for walk, name in zip(walks, names):
        walks_by_name[name].append(walk)
    
    return walks_by_name

def nearly_identical_alleles(allele1: str, allele2: str, threshold: int = 10):
    """
    Returns True if two sequences differ by at most one base.
    
    Args:
        allele1: String representation of an allele
        allele2: String representation of an allele
        threshold: Minimum length threshold - only applies the check if both alleles are >= threshold
    
    Returns:
        True if the alleles are nearly identical, False otherwise
    """
    # Only apply threshold check if both sequences are long enough
    # This prevents filtering out short sequences that are legitimately nearly identical
    if len(allele1) >= threshold and len(allele2) >= threshold:
        if len(allele1) + len(allele2) < 2 * threshold:
            return False

    length_difference = len(allele1) - len(allele2)
    if abs(length_difference) > 1:
        return False

    mismatches = 0
    idx1, idx2 = 0, 0
    while idx1 < len(allele1) and idx2 < len(allele2):
        base1 = allele1[idx1]
        base2 = allele2[idx2]
        if base1 == base2:
            idx1 += 1
            idx2 += 1
            continue
        mismatches += 1
        if mismatches > 1:
            return False
        if length_difference >= 0: # SNP or deletion
            idx1 += 1
        if length_difference <= 0: # SNP or insertion
            idx2 += 1
        
    return True

def log_action(log_path: str, action: str):
    # Get timestamp
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Measure memory in MB
    memory_mb = psutil.Process().memory_info().rss / 1024 / 1024

    # Build the header string
    log_entry = f"{timestamp},{memory_mb:.2f} MB,{action}\n"

    # Append to the end of the log file
    with open(log_path, "a") as log_file:
        log_file.write(log_entry)
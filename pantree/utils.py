import gzip
import pickle
import functools

def sequence_complement(s: str) -> str:
    def _base_complement(letter: str) -> str:
        base_dict: dict = {'A': 'T', 'C': 'G', 'T' : 'A', 'G': 'C'}
        return base_dict[letter] if letter in base_dict else letter

    return ''.join(_base_complement(c) for c in reversed(s))

@functools.lru_cache(maxsize=None)
def node_complement(s: str):
    return s[:-1] + _flip(s[-1])

@functools.lru_cache(maxsize=None)
def edge_complement(e: tuple[str, str]):
    return node_complement(e[1]), node_complement(e[0])

def walk_complement(w: list[str]) -> list[str]:
    return [node_complement(v) for v in w[::-1]]

@functools.lru_cache(maxsize=None)
def _flip(s):
    if s == '+':
        return '-'
    elif s == '-':
        return '+'
    else:
        raise ValueError()

def node_recover(node_id):
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

def get_from_biedge_dict(biedge_dict: dict, edge: tuple[str, str], default: any) -> any:
    if edge in biedge_dict:
        return biedge_dict[edge]
    return biedge_dict.get(edge_complement(edge), default)
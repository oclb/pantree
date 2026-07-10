"""
GFA file parsing utilities.
"""
import gzip
import json
import re
from dataclasses import dataclass


@dataclass
class GFAWalkLine:
    hap_name: str
    sample_name: str
    contig_start: int
    walk: list[str]

    @staticmethod
    def _get_hap_name(parts: list[str]) -> str:
        return '#'.join(parts[1:4])

    @classmethod
    def from_parts_P(cls, parts: list[str]):
        hap_name = parts[1]
        sample_name = parts[1]
        contig_start = 0
        walk = cls._parse_P_walk(parts[2])
        return cls(hap_name, sample_name, contig_start, walk)

    @classmethod
    def from_parts_W(cls, parts: list[str]):
        hap_name = cls._get_hap_name(parts)
        sample_name = parts[1]
        contig_start = int(parts[4]) if parts[4] != '*' else 0 # TODO
        walk = cls._parse_W_walk(parts[6])
        return cls(hap_name, sample_name, contig_start, walk)

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
    original_ids: list[str] | None = None

    @classmethod
    def from_parts(cls, parts: list[str]):
        return cls(parts[1], parts[2], cls._parse_original_ids(parts[3:]))

    @staticmethod
    def _parse_original_ids(optional_fields: list[str]) -> list[str] | None:
        for field in optional_fields:
            if not field.startswith('oi:J:'):
                continue
            try:
                original_ids = json.loads(field[5:])
            except json.JSONDecodeError:
                return None
            if isinstance(original_ids, list):
                return [str(original_id) for original_id in original_ids]
        return None


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


def read_gfa_line_by_line(filename: str, line_types: list[str]|None = None):
    """
    Read a GFA file line by line, yielding parsed line objects.

    Args:
        filename: Path to GFA file (can be .gz compressed)
        line_types: List of line types to parse (S, L, W, P). If None, parses all types.

    Yields:
        GFANodeLine, GFAEdgeLine, or GFAWalkLine objects depending on line type
    """
    if line_types is None:
        line_types = ['S', 'L', 'W', 'P']

    if filename.endswith('.gz'):
        opener = gzip.open
    else:
        opener = open

    with opener(filename, 'rt') as file:
        for line in file:
            parts = line.strip().split('\t')
            # Skip empty lines and header lines
            if not parts[0] or parts[0] not in line_types:
                continue
            line_getter = line_getters[parts[0]]
            yield line_getter(parts)

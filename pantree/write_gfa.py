"""
Write simplified NetworkX graph to GFA format.
"""
import json

import networkx as nx


def _optional_field(tag: str, field_type: str, value: str) -> str:
    return f"{tag}:{field_type}:{value}"


def _flip_orientation(orientation: str) -> str:
    return '-' if orientation == '+' else '+'


def _edge_complement_id(edge_id: tuple[str, str, str, str]) -> tuple[str, str, str, str]:
    u_segment, u_orient, v_segment, v_orient = edge_id
    return v_segment, _flip_orientation(v_orient), u_segment, _flip_orientation(u_orient)


def write_GFA(simplified_graph: nx.DiGraph, output_path: str, source_gfa: str | None = None) -> None:
    """
    Write a simplified NetworkX graph to GFA format.
    
    The function handles nodes that have been concatenated during simplification,
    assigning new node IDs and concatenating their sequences.
    
    Args:
        simplified_graph: NetworkX DiGraph with node attributes 'sequence' and 'direction'
        output_path: Path to write the GFA file
        source_gfa: Optional path to the source GFA file
    """
    
    # Create mapping from bidirectional node pairs to unique segment IDs
    # Nodes come in pairs like "node_+" and "node_-"
    segment_mapping = {}
    segment_counter = 1
    
    # Get all unique segments (each bidirectional pair gets one segment ID)
    processed_bases = set()
    for node in simplified_graph.nodes():
        # Extract base name (remove orientation suffix _+ or _-)
        if node.endswith('_+') or node.endswith('_-'):
            base_name = node[:-2]
        else:
            # Handle nodes without orientation suffix
            base_name = node
            
        if base_name not in processed_bases:
            processed_bases.add(base_name)
            segment_mapping[base_name] = f"s{segment_counter}"
            segment_counter += 1
    
    with open(output_path, 'w') as f:
        # Write header
        header_fields = ["H", _optional_field("VN", "Z", "1.0")]
        if source_gfa is not None:
            header_fields.append(_optional_field("og", "Z", source_gfa))
        f.write("\t".join(header_fields) + "\n")
        
        # Write S (Segment) lines
        # For each segment, use the forward orientation node to get the sequence
        for base_name, segment_id in sorted(segment_mapping.items(), key=lambda x: int(x[1][1:])):
            forward_node = f"{base_name}_+"
            
            if forward_node in simplified_graph.nodes():
                sequence = simplified_graph.nodes[forward_node].get('sequence', '*')
                origin_node = forward_node
            else:
                # Try without orientation suffix
                if base_name in simplified_graph.nodes():
                    sequence = simplified_graph.nodes[base_name].get('sequence', '*')
                    origin_node = base_name
                else:
                    sequence = '*'
                    origin_node = None

            if origin_node is not None:
                original_ids = simplified_graph.nodes[origin_node].get('original_ids', [base_name])
            else:
                original_ids = [base_name]
            origin_tag = _optional_field("oi", "J", json.dumps(original_ids, separators=(',', ':')))
            
            f.write(f"S\t{segment_id}\t{sequence}\t{origin_tag}\n")
        
        # Write L (Link) lines
        # Track written edges to avoid duplicates
        written_edges = set()
        
        for u, v in simplified_graph.edges():
            # Extract base names and orientations
            if u.endswith('_+'):
                u_base = u[:-2]
                u_orient = '+'
            elif u.endswith('_-'):
                u_base = u[:-2]
                u_orient = '-'
            else:
                u_base = u
                u_orient = '+'
                
            if v.endswith('_+'):
                v_base = v[:-2]
                v_orient = '+'
            elif v.endswith('_-'):
                v_base = v[:-2]
                v_orient = '-'
            else:
                v_base = v
                v_orient = '+'
            
            # Get segment IDs
            u_segment = segment_mapping.get(u_base, u_base)
            v_segment = segment_mapping.get(v_base, v_base)
            
            # Create edge identifier to avoid duplicates
            edge_id = (u_segment, u_orient, v_segment, v_orient)
            
            if edge_id not in written_edges and _edge_complement_id(edge_id) not in written_edges:
                written_edges.add(edge_id)
                # GFA format: L from_segment from_orient to_segment to_orient overlap
                # Using 0M for overlap (no overlap information)
                f.write(f"L\t{u_segment}\t{u_orient}\t{v_segment}\t{v_orient}\t0M\n")

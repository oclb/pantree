"""
Unit tests for haplotype position tracking functionality.
Tests the update_haplotype_positions method which tracks edge positions
across different haplotypes with priority-based updates.
"""
import unittest
import os
from graph_var.graph import PangenomeGraph
from graph_var.utils import read_gfa_line_by_line, GFAWalkLine


class TestHaplotypePositions(unittest.TestCase):
    """Test haplotype position tracking"""
    
    @classmethod
    def setUpClass(cls):
        """Load the graph once for all tests"""
        test_dir = os.path.dirname(__file__)
        cls.gfa_file = os.path.join(test_dir, "data", "simple_nested.gfa")
        cls.ref_name = 'ref'
        cls.G = PangenomeGraph.from_gfa_line_by_line(
            cls.gfa_file,
            ref_name=cls.ref_name
        )
    
    def test_reference_edge_positions_match_node_positions(self):
        """
        Test that for edges in the reference path, the primary_edge_pos matches
        the node position.
        
        The edge position represents the cumulative position after traversing node u,
        which should match the node's position (computed as cumulative length up to
        and including that node).
        """
        ref_name = self.ref_name
        
        # Check edges along the reference path
        ref_path = self.G.reference_path
        for i in range(len(ref_path) - 1):
            u, v = ref_path[i], ref_path[i + 1]
            
            if (u, v) not in self.G.edges:
                continue
            
            edge_data = self.G.edges[u, v]
            primary_edge_pos = edge_data.get('primary_edge_pos')
            
            # Skip if no position (e.g., terminal edges)
            if primary_edge_pos is None:
                continue
            
            hap_name, position = primary_edge_pos
            node_position = self.G.nodes[u].get('position')
            
            # For reference path edges, position should match
            self.assertIsNotNone(
                node_position,
                f"Node {u} on reference path should have a position"
            )
            self.assertEqual(
                node_position,
                position,
                f"Reference edge ({u}, {v}) has primary_edge_pos=({hap_name}, {position}), "
                f"but node {u} has position={node_position}"
            )
    
    def test_all_non_terminal_edges_have_haplotype_positions(self):
        """Test that all non-terminal edges have haplotype position data assigned"""
        edges_without_pos = []
        
        for (u, v), edge_data in self.G.edges.items():
            # Skip terminal edges (edges involving terminus nodes)
            if 'terminus' in u or 'terminus' in v:
                continue
            
            # Check if edge has any haplotype position data
            # Haplotype keys are in format: sample#haplotype#contig (e.g., 'ref#0#0', 'sample1#1#0')
            # Exclude known edge attributes
            excluded_keys = {'weight', 'is_in_tree', 'branch_point', 'is_back_edge', 'is_representative', 'index', 'is_inversion'}
            haplotype_keys = [key for key in edge_data.keys() if key not in excluded_keys and isinstance(edge_data[key], (int, float))]
            if not haplotype_keys:
                edges_without_pos.append((u, v))
        
        self.assertEqual(
            len(edges_without_pos),
            0,
            f"Found {len(edges_without_pos)} non-terminal edges without haplotype position data: {edges_without_pos[:5]}"
        )
    
    def test_haplotype_position_format(self):
        """Test that haplotype position data has the correct format (position values)"""
        for (u, v), edge_data in self.G.edges.items():
            # Check for haplotype keys
            haplotype_keys = [key for key in edge_data.keys() if key in ['ref', 'sample1_1', 'sample1_2', 'sample2_1', 'sample2_2']]
            
            for haplo_key in haplotype_keys:
                position = edge_data[haplo_key]
                
                # Should be an integer position
                self.assertIsInstance(
                    position,
                    int,
                    f"haplotype position for {haplo_key} should be an int, got {type(position)}"
                )
                self.assertGreaterEqual(
                    position,
                    0,
                    f"haplotype position for {haplo_key} should be >= 0, got {position}"
                )
    
    def test_complementary_edges_have_same_position(self):
        """Test that complementary edges have the same haplotype position data"""
        from graph_var.utils import edge_complement
        
        for (u, v), edge_data in self.G.edges.items():
            # Check for haplotype keys
            haplotype_keys = [key for key in edge_data.keys() if key in ['ref', 'sample1_1', 'sample1_2', 'sample2_1', 'sample2_2']]
            
            if not haplotype_keys:
                continue
            
            # Get the complementary edge
            comp_edge = edge_complement((u, v))
            if comp_edge not in self.G.edges:
                continue
            
            comp_edge_data = self.G.edges[comp_edge]
            
            # Check that both edges have the same haplotype keys and positions
            for haplo_key in haplotype_keys:
                self.assertIn(
                    haplo_key,
                    comp_edge_data,
                    f"Complementary edge {comp_edge} missing haplotype key {haplo_key}"
                )
                self.assertEqual(
                    edge_data[haplo_key],
                    comp_edge_data[haplo_key],
                    f"Complementary edges should have same position for {haplo_key}"
                )
    
    def test_haplotype_positions_are_assigned(self):
        """Test that haplotype positions are assigned to edges"""
        # Create a minimal graph to test position assignment
        test_dir = os.path.dirname(__file__)
        gfa_file = os.path.join(test_dir, "data", "simple_nested.gfa")
        
        # Load the graph
        G = PangenomeGraph.from_gfa_line_by_line(
            gfa_file,
            ref_name='ref'
        )
        
        # Check that edges have haplotype position data
        edges_with_positions = 0
        for (u, v), edge_data in G.edges.items():
            # Haplotype keys are in format: sample#haplotype#contig (e.g., 'ref#0#0', 'sample1#1#0')
            excluded_keys = {'weight', 'is_in_tree', 'branch_point', 'is_back_edge', 'is_representative', 'index', 'is_inversion'}
            haplotype_keys = [key for key in edge_data.keys() if key not in excluded_keys and isinstance(edge_data[key], (int, float))]
            if haplotype_keys:
                edges_with_positions += 1
        
        # Should have some edges with haplotype positions
        self.assertGreater(
            edges_with_positions,
            0,
            "Should have at least some edges with haplotype position data"
        )
    
    def test_position_increases_along_walk(self):
        """Test that positions increase monotonically along a walk"""
        # Read walks from the GFA file
        test_dir = os.path.dirname(__file__)
        gfa_file = os.path.join(test_dir, "data", "simple_nested.gfa")
        
        for line in read_gfa_line_by_line(gfa_file):
            if not isinstance(line, GFAWalkLine):
                continue
            
            walk = line.walk
            hap_name = line.hap_name
            
            # Track positions along the walk
            positions = []
            for u, v in zip(walk[:-1], walk[1:]):
                if (u, v) in self.G.edges:
                    edge_data = self.G.edges[u, v]
                    if hap_name in edge_data:
                        positions.append(edge_data[hap_name])
            
            # Positions should be monotonically increasing
            if len(positions) > 1:
                for i in range(len(positions) - 1):
                    self.assertLessEqual(
                        positions[i],
                        positions[i + 1],
                        f"Positions should increase along walk {hap_name}: "
                        f"{positions[i]} > {positions[i + 1]} at index {i}"
                    )


class TestHaplotypePositionsC4A(unittest.TestCase):
    """Test haplotype positions on the C4A inversion graph"""
    
    @classmethod
    def setUpClass(cls):
        """Load the C4A graph once for all tests"""
        test_dir = os.path.dirname(__file__)
        cls.gfa_file = os.path.join(test_dir, "data", "c4a_with_inversion_and_sequences.gfa")
        cls.ref_name = 'GRCh38'
        cls.G = PangenomeGraph.from_gfa_line_by_line(
            cls.gfa_file,
            ref_name=cls.ref_name
        )
    
    def test_reference_edge_positions_match_node_positions(self):
        """
        Test that for edges in the reference path, the primary_edge_pos matches
        the node position.
        """
        ref_name = self.ref_name
        
        # Check edges along the reference path
        ref_path = self.G.reference_path
        mismatches = []
        
        for i in range(len(ref_path) - 1):
            u, v = ref_path[i], ref_path[i + 1]
            
            if (u, v) not in self.G.edges:
                continue
            
            edge_data = self.G.edges[u, v]
            primary_edge_pos = edge_data.get('primary_edge_pos')
            
            if primary_edge_pos is None:
                continue
            
            hap_name, position = primary_edge_pos
            node_position = self.G.nodes[u].get('position')
            
            if node_position != position:
                mismatches.append((u, v, position, node_position))
        
        self.assertEqual(
            len(mismatches),
            0,
            f"Found {len(mismatches)} reference path edges with position mismatches: {mismatches[:5]}"
        )
    
    def test_reference_path_edges_have_positions(self):
        """Test that all non-terminal edges in the reference path have haplotype position data"""
        ref_path = self.G.reference_path
        edges_without_pos = []
        
        for i in range(len(ref_path) - 1):
            u, v = ref_path[i], ref_path[i + 1]
            
            # Skip terminal edges
            if 'terminus' in u or 'terminus' in v:
                continue
            
            if (u, v) not in self.G.edges:
                continue
            
            edge_data = self.G.edges[u, v]
            # Check for any haplotype position data (exclude standard graph attributes)
            haplotype_keys = [key for key in edge_data.keys() if key not in ['weight', 'direction', 'sequence', 'position', 'branch_point', 'index', 'is_back_edge', 'is_in_tree', 'is_representative']]
            if not haplotype_keys:
                edges_without_pos.append((u, v))
        
        self.assertEqual(
            len(edges_without_pos),
            0,
            f"Found {len(edges_without_pos)} non-terminal reference path edges without positions: {edges_without_pos}"
        )


if __name__ == '__main__':
    unittest.main()

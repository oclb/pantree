"""
Unit tests for haplotype position tracking functionality.
Tests the update_haplotype_positions method which tracks edge positions
across different haplotypes with priority-based updates.
"""
import unittest
import os
from graph_var.graph import PangenomeGraph
from graph_var.gfa import read_gfa_line_by_line, GFAWalkLine


class TestHaplotypePositions(unittest.TestCase):
    """Test haplotype position tracking"""
    
    @classmethod
    def setUpClass(cls):
        """Load the graph once for all tests"""
        test_dir = os.path.dirname(__file__)
        cls.gfa_file = os.path.join(test_dir, "data", "simple_nested.gfa")
        cls.ref_name = 'ref'
        priority_dict = {
            'ref': 0,
            'sample1': 1,
            'sample2': 2,
        }
        cls.G = PangenomeGraph.from_gfa_line_by_line(
            cls.gfa_file,
            ref_name=cls.ref_name,
            priority_dict=priority_dict
        )
    
    def test_reference_nodes_have_haplotype_positions(self):
        """
        Test that nodes in the reference path have haplotype position data.
        
        The node position represents the cumulative bp offset after traversing that node.
        """
        ref_name = self.ref_name
        
        # Check nodes along the reference path
        ref_path = self.G.reference_path
        for node in ref_path:
            # Skip terminus nodes
            if 'terminus' in node:
                continue
            
            node_data = self.G.nodes[node]
            
            # Node should have the reference haplotype position (format: ref#0#0)
            # Find any key that starts with ref_name
            ref_keys = [k for k in node_data.keys() if k.startswith(ref_name + '#')]
            self.assertGreater(
                len(ref_keys),
                0,
                f"Node {node} on reference path should have haplotype position for {ref_name}"
            )
            
            # Position should be a valid integer
            position = node_data[ref_keys[0]]
            self.assertIsInstance(
                position,
                int,
                f"Node {node} haplotype position should be int, got {type(position)}"
            )
            self.assertGreaterEqual(
                position,
                0,
                f"Node {node} haplotype position should be >= 0, got {position}"
            )
    
    def test_all_non_terminal_nodes_have_haplotype_positions(self):
        """Test that all non-terminal nodes have haplotype position data assigned"""
        nodes_without_pos = []
        
        for node, node_data in self.G.nodes.items():
            # Skip terminal nodes
            if 'terminus' in node:
                continue
            
            # Check if node has any haplotype position data
            # Haplotype keys are in format: sample#haplotype#contig (e.g., 'ref#0#0', 'sample1#1#0')
            # Exclude known node attributes
            excluded_keys = {'direction', 'sequence', 'position', 'right_position', 'distance_from_reference', 'on_reference_path'}
            haplotype_keys = [key for key in node_data.keys() if key not in excluded_keys and isinstance(node_data[key], (int, float))]
            if not haplotype_keys:
                nodes_without_pos.append(node)
        
        self.assertEqual(
            len(nodes_without_pos),
            0,
            f"Found {len(nodes_without_pos)} non-terminal nodes without haplotype position data: {nodes_without_pos[:5]}"
        )
    
    def test_haplotype_position_format(self):
        """Test that haplotype position data has the correct format (position values)"""
        for node, node_data in self.G.nodes.items():
            # Skip terminus nodes
            if 'terminus' in node:
                continue
            
            # Check for haplotype keys (exclude known node attributes)
            excluded_keys = {'direction', 'sequence', 'position', 'right_position', 'distance_from_reference', 'on_reference_path'}
            haplotype_keys = [key for key in node_data.keys() if key not in excluded_keys]
            
            for haplo_key in haplotype_keys:
                position = node_data[haplo_key]
                
                # Should be an integer position
                self.assertIsInstance(
                    position,
                    int,
                    f"Node {node} haplotype position for {haplo_key} should be an int, got {type(position)}"
                )
                self.assertGreaterEqual(
                    position,
                    0,
                    f"Node {node} haplotype position for {haplo_key} should be >= 0, got {position}"
                )
    
    def test_complementary_nodes_have_same_position(self):
        """Test that complementary nodes have the same haplotype position data"""
        from graph_var.utils import node_complement
        
        for node, node_data in self.G.nodes.items():
            # Skip terminus nodes
            if 'terminus' in node:
                continue
            
            # Check for haplotype keys (exclude known node attributes)
            excluded_keys = {'direction', 'sequence', 'position', 'right_position', 'distance_from_reference', 'on_reference_path'}
            haplotype_keys = [key for key in node_data.keys() if key not in excluded_keys]
            
            if not haplotype_keys:
                continue
            
            # Get the complementary node
            comp_node = node_complement(node)
            if comp_node not in self.G.nodes:
                continue
            
            comp_node_data = self.G.nodes[comp_node]
            
            # Check that both nodes have the same haplotype keys and positions
            for haplo_key in haplotype_keys:
                self.assertIn(
                    haplo_key,
                    comp_node_data,
                    f"Complementary node {comp_node} missing haplotype key {haplo_key}"
                )
                self.assertEqual(
                    node_data[haplo_key],
                    comp_node_data[haplo_key],
                    f"Complementary nodes {node} and {comp_node} should have same position for {haplo_key}"
                )
    
    def test_haplotype_positions_are_assigned(self):
        """Test that haplotype positions are assigned to nodes"""
        # Create a minimal graph to test position assignment
        test_dir = os.path.dirname(__file__)
        gfa_file = os.path.join(test_dir, "data", "simple_nested.gfa")
        
        priority_dict = {
            'ref': 0,
            'sample1': 1,
            'sample2': 2,
        }
        
        # Load the graph
        G = PangenomeGraph.from_gfa_line_by_line(
            gfa_file,
            ref_name='ref',
            priority_dict=priority_dict
        )
        
        # Check that nodes have haplotype position data
        nodes_with_positions = 0
        for node, node_data in G.nodes.items():
            # Skip terminus nodes
            if 'terminus' in node:
                continue
            
            # Haplotype keys are in format: sample#haplotype#contig (e.g., 'ref#0#0', 'sample1#1#0')
            excluded_keys = {'direction', 'sequence', 'position', 'right_position', 'distance_from_reference', 'on_reference_path'}
            haplotype_keys = [key for key in node_data.keys() if key not in excluded_keys and isinstance(node_data[key], (int, float))]
            if haplotype_keys:
                nodes_with_positions += 1
        
        # Should have some nodes with haplotype positions
        self.assertGreater(
            nodes_with_positions,
            0,
            "Should have at least some nodes with haplotype position data"
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
            for node in walk:
                if node in self.G.nodes:
                    node_data = self.G.nodes[node]
                    if hap_name in node_data:
                        positions.append(node_data[hap_name])
            
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
    
    def test_reference_nodes_have_haplotype_positions(self):
        """
        Test that nodes in the reference path have haplotype position data.
        """
        ref_name = self.ref_name
        
        # Check nodes along the reference path
        ref_path = self.G.reference_path
        nodes_without_pos = []
        
        for node in ref_path:
            # Skip terminus nodes
            if 'terminus' in node:
                continue
            
            node_data = self.G.nodes[node]
            
            # Node should have the reference haplotype position (format: ref#0#0)
            ref_keys = [k for k in node_data.keys() if k.startswith(ref_name + '#')]
            if len(ref_keys) == 0:
                nodes_without_pos.append(node)
        
        self.assertEqual(
            len(nodes_without_pos),
            0,
            f"Found {len(nodes_without_pos)} reference path nodes without haplotype position: {nodes_without_pos[:5]}"
        )
    
    def test_reference_path_nodes_have_positions(self):
        """Test that all non-terminal nodes in the reference path have haplotype position data"""
        ref_path = self.G.reference_path
        nodes_without_pos = []
        
        for node in ref_path:
            # Skip terminal nodes
            if 'terminus' in node:
                continue
            
            node_data = self.G.nodes[node]
            # Check for any haplotype position data (exclude standard node attributes)
            excluded_keys = {'direction', 'sequence', 'position', 'right_position', 'distance_from_reference', 'on_reference_path'}
            haplotype_keys = [key for key in node_data.keys() if key not in excluded_keys and isinstance(node_data[key], (int, float))]
            if not haplotype_keys:
                nodes_without_pos.append(node)
        
        self.assertEqual(
            len(nodes_without_pos),
            0,
            f"Found {len(nodes_without_pos)} non-terminal reference path nodes without positions: {nodes_without_pos}"
        )


if __name__ == '__main__':
    unittest.main()

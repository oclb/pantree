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
    
    def test_all_non_terminal_edges_have_primary_edge_pos(self):
        """Test that all non-terminal edges have primary_edge_pos assigned"""
        edges_without_pos = []
        
        for (u, v), edge_data in self.G.edges.items():
            # Skip terminal edges (edges involving terminus nodes)
            if 'terminus' in u or 'terminus' in v:
                continue
            
            if 'primary_edge_pos' not in edge_data:
                edges_without_pos.append((u, v))
        
        self.assertEqual(
            len(edges_without_pos),
            0,
            f"Found {len(edges_without_pos)} non-terminal edges without primary_edge_pos: {edges_without_pos[:5]}"
        )
    
    def test_primary_edge_pos_format(self):
        """Test that primary_edge_pos has the correct format (hap_name, position)"""
        for (u, v), edge_data in self.G.edges.items():
            primary_edge_pos = edge_data.get('primary_edge_pos')
            
            if primary_edge_pos is None:
                continue
            
            # Should be a tuple of (str, int)
            self.assertIsInstance(
                primary_edge_pos,
                tuple,
                f"primary_edge_pos should be a tuple, got {type(primary_edge_pos)}"
            )
            self.assertEqual(
                len(primary_edge_pos),
                2,
                f"primary_edge_pos should have 2 elements, got {len(primary_edge_pos)}"
            )
            
            hap_name, position = primary_edge_pos
            self.assertIsInstance(
                hap_name,
                str,
                f"Haplotype name should be str, got {type(hap_name)}"
            )
            self.assertIsInstance(
                position,
                int,
                f"Position should be int, got {type(position)}"
            )
            self.assertGreaterEqual(
                position,
                0,
                f"Position should be non-negative, got {position}"
            )
    
    def test_complementary_edges_have_same_position(self):
        """Test that complementary edges have the same primary_edge_pos"""
        from graph_var.utils import edge_complement
        
        for (u, v), edge_data in self.G.edges.items():
            primary_edge_pos = edge_data.get('primary_edge_pos')
            
            if primary_edge_pos is None:
                continue
            
            # Get the complementary edge
            comp_edge = edge_complement((u, v))
            if comp_edge not in self.G.edges:
                continue
            
            comp_edge_data = self.G.edges[comp_edge]
            comp_primary_edge_pos = comp_edge_data.get('primary_edge_pos')
            
            self.assertEqual(
                primary_edge_pos,
                comp_primary_edge_pos,
                f"Edge ({u}, {v}) has primary_edge_pos={primary_edge_pos}, "
                f"but its complement {comp_edge} has {comp_primary_edge_pos}"
            )
    
    def test_priority_dict_affects_position_assignment(self):
        """Test that priority_dict correctly prioritizes haplotypes"""
        # Create a minimal graph to test priority
        test_dir = os.path.dirname(__file__)
        gfa_file = os.path.join(test_dir, "data", "simple_nested.gfa")
        
        # Load with different priority
        priority_dict = {'ref': 0, 'sample1': 1, 'sample2': 2}
        G_with_priority = PangenomeGraph.from_gfa_line_by_line(
            gfa_file,
            ref_name='ref',
            priority_dict=priority_dict
        )
        
        # Check that reference positions are preserved
        ref_edge_count = 0
        for (u, v), edge_data in G_with_priority.edges.items():
            primary_edge_pos = edge_data.get('primary_edge_pos')
            if primary_edge_pos and primary_edge_pos[0] == 'ref':
                ref_edge_count += 1
        
        # Should have some edges with reference as primary
        self.assertGreater(
            ref_edge_count,
            0,
            "Should have at least some edges with reference as primary haplotype"
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
                    primary_edge_pos = edge_data.get('primary_edge_pos')
                    if primary_edge_pos and primary_edge_pos[0] == hap_name:
                        positions.append(primary_edge_pos[1])
            
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
        """Test that all non-terminal edges in the reference path have primary_edge_pos"""
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
            if 'primary_edge_pos' not in edge_data:
                edges_without_pos.append((u, v))
        
        self.assertEqual(
            len(edges_without_pos),
            0,
            f"Found {len(edges_without_pos)} non-terminal reference path edges without positions: {edges_without_pos}"
        )


if __name__ == '__main__':
    unittest.main()

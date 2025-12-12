"""
Test for GFA file ordering issue where W-lines appear before L-lines.
This tests the fix for the issue where compute_edge_weights is called
before edges are added to the graph.
"""

import unittest
import os
from pantree import PangenomeGraph


class TestGFAOrdering(unittest.TestCase):
    """Test that GFA files with different line orderings are handled correctly."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_data_dir = os.path.join(os.path.dirname(__file__), 'data')
        self.correct_order_gfa = os.path.join(self.test_data_dir, 'simple_nested.gfa')
        self.wrong_order_gfa = os.path.join(self.test_data_dir, 'simple_nested_wrong_order.gfa')
    
    def test_gfa_with_correct_order(self):
        """Test that GFA with standard ordering (S, L, W) works correctly."""
        G = PangenomeGraph.from_gfa_line_by_line(
            self.correct_order_gfa,
            ref_name='ref'
        )
        
        # Verify graph was constructed properly
        self.assertGreater(len(G.nodes), 0, "Graph should have nodes")
        self.assertGreater(len(G.edges), 0, "Graph should have edges")
        self.assertIsNotNone(G.reference_path, "Reference path should be set")
        self.assertGreater(len(G.reference_path), 0, "Reference path should not be empty")
    
    def test_gfa_with_walks_before_links(self):
        """Test that GFA with W-lines before L-lines (S, W, L) works correctly.
        
        This reproduces the issue where compute_edge_weights is called on walks
        before edges are added to the graph, causing an error.
        """
        # This should not raise an error
        G = PangenomeGraph.from_gfa_line_by_line(
            self.wrong_order_gfa,
            ref_name='ref'
        )
        
        # Verify graph was constructed properly
        self.assertGreater(len(G.nodes), 0, "Graph should have nodes")
        self.assertGreater(len(G.edges), 0, "Graph should have edges")
        self.assertIsNotNone(G.reference_path, "Reference path should be set")
        self.assertGreater(len(G.reference_path), 0, "Reference path should not be empty")
    
    def test_both_orderings_produce_same_graph(self):
        """Test that both orderings produce equivalent graphs."""
        G_correct = PangenomeGraph.from_gfa_line_by_line(
            self.correct_order_gfa,
            ref_name='ref'
        )
        
        G_wrong = PangenomeGraph.from_gfa_line_by_line(
            self.wrong_order_gfa,
            ref_name='ref'
        )
        
        # Both should have the same number of nodes and edges
        self.assertEqual(len(G_correct.nodes), len(G_wrong.nodes),
                        "Both graphs should have the same number of nodes")
        self.assertEqual(len(G_correct.edges), len(G_wrong.edges),
                        "Both graphs should have the same number of edges")
        self.assertEqual(len(G_correct.reference_path), len(G_wrong.reference_path),
                        "Both graphs should have the same reference path length")


if __name__ == '__main__':
    unittest.main()

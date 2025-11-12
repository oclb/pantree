"""
Extended unit tests for utils module
"""
import unittest
import os
import tempfile
from graph_var.utils import (
    sequence_complement,
    node_complement,
    edge_complement,
    walk_complement,
)
from graph_var.gfa import read_gfa_line_by_line, GFANodeLine, GFAEdgeLine, GFAWalkLine


class TestSequenceOperations(unittest.TestCase):
    """Test sequence complement operations"""
    
    def test_sequence_complement(self):
        """Test DNA sequence complement"""
        self.assertEqual(sequence_complement('ACGT'), 'ACGT')
        self.assertEqual(sequence_complement('A'), 'T')
        self.assertEqual(sequence_complement('T'), 'A')
        self.assertEqual(sequence_complement('C'), 'G')
        self.assertEqual(sequence_complement('G'), 'C')
        self.assertEqual(sequence_complement('AAAA'), 'TTTT')
        self.assertEqual(sequence_complement('ATCG'), 'CGAT')
    
    def test_node_complement(self):
        """Test node complement operation"""
        self.assertEqual(node_complement('1_+'), '1_-')
        self.assertEqual(node_complement('1_-'), '1_+')
        self.assertEqual(node_complement('node_+'), 'node_-')
        self.assertEqual(node_complement('node_-'), 'node_+')
    
    def test_edge_complement(self):
        """Test edge complement operation"""
        edge = ('1_+', '2_+')
        comp_edge = edge_complement(edge)
        self.assertEqual(comp_edge, ('2_-', '1_-'))
        
        # Double complement should give original
        double_comp = edge_complement(comp_edge)
        self.assertEqual(double_comp, edge)
    
    def test_walk_complement(self):
        """Test walk complement operation"""
        walk = ['1_+', '2_+', '3_+']
        comp_walk = walk_complement(walk)
        self.assertEqual(comp_walk, ['3_-', '2_-', '1_-'])
        
        # Single node walk
        single_walk = ['1_+']
        self.assertEqual(walk_complement(single_walk), ['1_-'])


class TestGFAReading(unittest.TestCase):
    """Test GFA file reading"""
    
    def test_read_gfa_line_by_line(self):
        """Test reading GFA file line by line"""
        
        test_dir = os.path.dirname(__file__)
        gfa_file = os.path.join(test_dir, "data", "simple_nested.gfa")
        
        segments = []
        links = []
        walks = []
        
        for line_obj in read_gfa_line_by_line(gfa_file):
            if isinstance(line_obj, GFANodeLine):
                segments.append(line_obj)
            elif isinstance(line_obj, GFAEdgeLine):
                links.append(line_obj)
            elif isinstance(line_obj, GFAWalkLine):
                walks.append(line_obj)
        
        # Should have segments, links, and walks
        self.assertGreater(len(segments), 0)
        self.assertGreater(len(links), 0)
        self.assertGreater(len(walks), 0)
        
        # Check segment format (GFANodeLine dataclass)
        for seg in segments:
            self.assertIsInstance(seg, GFANodeLine)
            self.assertIsInstance(seg.node_id, str)
            self.assertIsInstance(seg.sequence, str)
        
        # Check link format (GFAEdgeLine dataclass)
        for link in links:
            self.assertIsInstance(link, GFAEdgeLine)
            self.assertIsInstance(link.u, str)
            self.assertIsInstance(link.v, str)
        
        # Check walk format (GFAWalkLine dataclass)
        for walk in walks:
            self.assertIsInstance(walk, GFAWalkLine)
            self.assertIsInstance(walk.hap_name, str)
            self.assertIsInstance(walk.walk, list)
            # Check that walk contains node IDs in the expected format
            for node_id in walk.walk:
                self.assertIsInstance(node_id, str)
                self.assertTrue(node_id.endswith('_+') or node_id.endswith('_-'))


if __name__ == '__main__':
    unittest.main()

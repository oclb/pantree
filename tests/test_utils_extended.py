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
    merge_dicts,
    group_walks_by_name,
    read_gfa_line_by_line,
)


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


class TestDictOperations(unittest.TestCase):
    """Test dictionary operations"""
    
    def test_merge_dicts_empty(self):
        """Test merging empty dictionaries"""
        result = merge_dicts([])
        self.assertEqual(result, {})
    
    def test_merge_dicts_single(self):
        """Test merging single dictionary"""
        d = {'a': 1, 'b': 2}
        result = merge_dicts([d])
        self.assertEqual(result, d)
    
    def test_merge_dicts_multiple(self):
        """Test merging multiple dictionaries"""
        d1 = {'a': 1, 'b': 2}
        d2 = {'b': 3, 'c': 4}
        d3 = {'a': 5, 'd': 6}
        result = merge_dicts([d1, d2, d3])
        self.assertEqual(result, {'a': 6, 'b': 5, 'c': 4, 'd': 6})
    
    def test_merge_dicts_no_overlap(self):
        """Test merging dictionaries with no overlapping keys"""
        d1 = {'a': 1}
        d2 = {'b': 2}
        d3 = {'c': 3}
        result = merge_dicts([d1, d2, d3])
        self.assertEqual(result, {'a': 1, 'b': 2, 'c': 3})


class TestWalkGrouping(unittest.TestCase):
    """Test walk grouping operations"""
    
    def test_group_walks_by_name_empty(self):
        """Test grouping empty walks"""
        result = group_walks_by_name([], [])
        self.assertEqual(result, {})
    
    def test_group_walks_by_name_single(self):
        """Test grouping single walk"""
        walks = [['1_+', '2_+']]
        names = ['sample1']
        result = group_walks_by_name(walks, names)
        self.assertEqual(result, {'sample1': [['1_+', '2_+']]})
    
    def test_group_walks_by_name_multiple_same(self):
        """Test grouping multiple walks with same name"""
        walks = [['1_+', '2_+'], ['3_+', '4_+']]
        names = ['sample1', 'sample1']
        result = group_walks_by_name(walks, names)
        self.assertEqual(result, {'sample1': [['1_+', '2_+'], ['3_+', '4_+']]})
    
    def test_group_walks_by_name_multiple_different(self):
        """Test grouping multiple walks with different names"""
        walks = [['1_+', '2_+'], ['3_+', '4_+'], ['5_+', '6_+']]
        names = ['sample1', 'sample2', 'sample1']
        result = group_walks_by_name(walks, names)
        self.assertEqual(result, {
            'sample1': [['1_+', '2_+'], ['5_+', '6_+']],
            'sample2': [['3_+', '4_+']]
        })


class TestGFAReading(unittest.TestCase):
    """Test GFA file reading"""
    
    def test_read_gfa_line_by_line(self):
        """Test reading GFA file line by line"""
        from graph_var.utils import GFANodeLine, GFAEdgeLine, GFAWalkLine
        
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

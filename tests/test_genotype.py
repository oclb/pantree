"""
Unit tests for genotype-related functionality
"""
import unittest
import os
from graph_var.graph import PangenomeGraph


class TestGenotype(unittest.TestCase):
    """Test genotype computation and related methods"""
    
    @classmethod
    def setUpClass(cls):
        """Load the test GFA file once for all tests"""
        test_dir = os.path.dirname(__file__)
        cls.gfa_file = os.path.join(test_dir, "data", "simple_nested.gfa")
        cls.G = PangenomeGraph.from_gfa_line_by_line(cls.gfa_file, ref_name='ref')
    
    def test_genotype_with_walk(self):
        """Test genotype computation from a walk"""
        # Use a walk in the correct format (list of node_ids)
        walk = ['1_+', '2_+', '4_+', '9_+', '10_+', '11_+']
        result = self.G.genotype(walk, return_linear_coverage=True)
        
        # Should return a tuple of (cr_dict, ca_dict, linear_coverage)
        self.assertIsInstance(result, tuple)
        self.assertEqual(len(result), 3)
        cr_dict, ca_dict, linear_coverage = result
        self.assertIsInstance(cr_dict, dict)
        self.assertIsInstance(ca_dict, dict)
        # linear_coverage is a tuple of (min_pos, max_pos)
        self.assertIsInstance(linear_coverage, tuple)
    
    def test_get_missing_variants(self):
        """Test missing variant detection from linear coverage"""
        walk = ['1_+', '2_+', '4_+', '9_+', '10_+', '11_+']
        _, _, linear_coverage = self.G.genotype(walk, return_linear_coverage=True)
        
        missing = self.G.get_missing_variants([linear_coverage], exclude_terminus=True)
        self.assertIsInstance(missing, list)
    
    def test_count_edge_visits(self):
        """Test edge visit counting from genotype"""
        # Test with reference walk
        walk1 = ['1_+', '2_+', '4_+', '9_+', '10_+', '11_+']
        cr_dict1, ca_dict1, _ = self.G.genotype(walk1, return_linear_coverage=True)
        
        # count_edge_visits expects only variant edges (alt alleles)
        genotype1 = ca_dict1
        
        # Count edge visits
        edge_visits1 = self.G.count_edge_visits(genotype1)
        self.assertIsInstance(edge_visits1, dict)
        # Should have at least the variant edges from genotype
        for edge in genotype1:
            self.assertIn(edge, edge_visits1)
        
        # Test with alternative walk that takes different path
        walk2 = ['1_+', '3_+', '4_+', '5_+', '7_+', '8_+', '9_+', '11_+']
        cr_dict2, ca_dict2, _ = self.G.genotype(walk2, return_linear_coverage=True)
        genotype2 = ca_dict2
        
        edge_visits2 = self.G.count_edge_visits(genotype2)
        self.assertIsInstance(edge_visits2, dict)
        
        # Different walks should produce different edge visits
        self.assertNotEqual(edge_visits1, edge_visits2)
    
    def test_count_edge_visits_invalid_genotype(self):
        """Test that invalid genotype raises ValueError"""
        # Get a valid genotype first
        walk = ['1_+', '2_+', '4_+', '9_+', '10_+', '11_+']
        _, ca_dict, _ = self.G.genotype(walk, return_linear_coverage=True)
        genotype = ca_dict.copy()
        
        # Make it invalid by setting a visit count to 2 (graph is acyclic, so this is impossible)
        if genotype:
            first_edge = list(genotype.keys())[0]
            genotype[first_edge] = 2
            
            # Should raise ValueError for invalid genotype
            with self.assertRaises(ValueError) as context:
                self.G.count_edge_visits(genotype)
            
            self.assertIn("does not correspond to any valid walk", str(context.exception))


class TestPositionAndDistance(unittest.TestCase):
    """Test position and distance computation"""
    
    @classmethod
    def setUpClass(cls):
        """Load the test GFA file once for all tests"""
        test_dir = os.path.dirname(__file__)
        cls.gfa_file = os.path.join(test_dir, "data", "simple_nested.gfa")
        cls.G = PangenomeGraph.from_gfa_line_by_line(cls.gfa_file, ref_name='ref')
    
    def test_nodes_have_positions(self):
        """Test that nodes have position attributes"""
        for node in self.G.nodes():
            if self.G.is_terminal(node):
                continue
            self.assertIn('position', self.G.nodes[node])
            self.assertIsInstance(self.G.nodes[node]['position'], (int, float))
    
    def test_nodes_have_distance_from_reference(self):
        """Test that nodes have distance_from_reference attributes"""
        for node in self.G.nodes():
            if self.G.is_terminal(node):
                continue
            self.assertIn('distance_from_reference', self.G.nodes[node])
            self.assertIsInstance(self.G.nodes[node]['distance_from_reference'], (int, float))
    
    def test_get_vcf_position(self):
        """Test VCF position calculation"""
        for edge in self.G.variant_edges:
            pos = self.G.get_vcf_position(edge, prepend_letter_to_alleles=False)
            self.assertIsInstance(pos, int)
            self.assertGreater(pos, 0)


class TestEdgeProperties(unittest.TestCase):
    """Test edge property methods"""
    
    @classmethod
    def setUpClass(cls):
        """Load the test GFA file once for all tests"""
        test_dir = os.path.dirname(__file__)
        cls.gfa_file = os.path.join(test_dir, "data", "simple_nested.gfa")
        cls.G = PangenomeGraph.from_gfa_line_by_line(cls.gfa_file, ref_name='ref')
    
    def test_is_in_tree(self):
        """Test tree membership detection"""
        for edge in self.G.edges():
            result = self.G.is_in_tree(edge)
            self.assertIsInstance(result, bool)
    
    def test_is_terminal(self):
        """Test terminal node/edge detection"""
        # Test with nodes
        self.assertTrue(self.G.is_terminal('+_terminus_+'))
        self.assertTrue(self.G.is_terminal('-_terminus_+'))
        
        # Test with regular nodes
        for node in self.G.nodes():
            if node.startswith('+_terminus') or node.startswith('-_terminus'):
                self.assertTrue(self.G.is_terminal(node))
            else:
                self.assertFalse(self.G.is_terminal(node))
    
    def test_parent_in_tree(self):
        """Test parent node retrieval in tree"""
        for node in self.G.reference_path[1:]:  # Skip first node
            if self.G.is_terminal(node):
                continue
            parent = self.G.parent_in_tree(node)
            # Parent should be a string or None
            self.assertTrue(parent is None or isinstance(parent, str))
    
    def test_sorted_variant_edges(self):
        """Test that variant edges can be sorted"""
        sorted_edges = list(self.G.sorted_variant_edges(exclude_terminus=True))
        self.assertIsInstance(sorted_edges, list)
        self.assertGreater(len(sorted_edges), 0)


class TestSNPAndMNP(unittest.TestCase):
    """Test SNP and MNP detection"""
    
    @classmethod
    def setUpClass(cls):
        """Load the test GFA file once for all tests"""
        test_dir = os.path.dirname(__file__)
        cls.gfa_file = os.path.join(test_dir, "data", "simple_nested.gfa")
        cls.G = PangenomeGraph.from_gfa_line_by_line(cls.gfa_file, ref_name='ref')
    
    def test_is_snp(self):
        """Test SNP detection"""
        snps = [e for e in self.G.variant_edges if self.G.is_snp(e)]
        self.assertIsInstance(snps, list)
    
    def test_is_mnp(self):
        """Test MNP detection"""
        mnps = [e for e in self.G.variant_edges if self.G.is_mnp(e)]
        self.assertIsInstance(mnps, list)


if __name__ == '__main__':
    unittest.main()

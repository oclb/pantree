"""
Tests for genotype transformation functions.
"""

import unittest
from pantree import PangenomeGraph
from pantree.genotype_transformation import transform_genotype_via_walk


class TestGenotypeTransformation(unittest.TestCase):
    """Test genotype transformation via walk construction"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.G = PangenomeGraph.from_gfa('tests/data/simple_nested.gfa', ref_name='ref')
    
    def test_transform_empty_genotype(self):
        """Test transformation with empty genotype"""
        genotype_A = {}
        genotype_B = transform_genotype_via_walk(self.G, self.G, genotype_A)
        
        self.assertIsInstance(genotype_B, dict)
        self.assertEqual(len(genotype_B), 0)
    
    def test_transform_single_variant_success(self):
        """Test transformation with a single variant that forms a valid walk"""
        # Use an actual variant edge from the graph
        variant_edge = ('3_+', '4_+')
        genotype_A = {variant_edge: 1}
        
        genotype_B = transform_genotype_via_walk(self.G, self.G, genotype_A)
        
        self.assertIsInstance(genotype_B, dict)
        # When A == B, the same variant should appear in B
        self.assertIn(variant_edge, genotype_B)
        self.assertEqual(genotype_B[variant_edge], 1)
    

    def test_transform_with_single_variant(self):
        """Test transformation with a single variant"""
        # Use an actual variant edge from the graph
        genotype_A = {('3_+', '4_+'): 1}
        
        genotype_B = transform_genotype_via_walk(self.G, self.G, genotype_A)
        
        self.assertIsInstance(genotype_B, dict)
        # The variant should be in the result
        self.assertIn(('3_+', '4_+'), genotype_B)
    
    def test_transform_returns_dict(self):
        """Test that transformation returns a dictionary for valid variants"""
        # Use an actual variant edge from the graph that forms a valid walk
        variant_edge = ('3_+', '4_+')
        genotype_A = {variant_edge: 1}
        
        genotype_B = transform_genotype_via_walk(self.G, self.G, genotype_A)
        
        self.assertIsInstance(genotype_B, dict)
        for key in genotype_B.keys():
            self.assertIsInstance(key, tuple)
            self.assertEqual(len(key), 2)
        for value in genotype_B.values():
            self.assertIsInstance(value, int)
            self.assertGreater(value, 0)
    
    def test_transform_raises_value_error_for_invalid_variants(self):
        """Test that ValueError is raised for invalid variant combinations"""
        # Use a non-existent edge that is not a variant edge
        genotype_A = {('99_+', '100_+'): 1}
        
        # This should raise ValueError because this edge is not a variant edge
        with self.assertRaises(ValueError) as context:
            genotype_B = transform_genotype_via_walk(self.G, self.G, genotype_A)
        
        # Check that the error message mentions variant edge
        self.assertIn("variant edge", str(context.exception).lower())
    
    def test_transform_with_different_graphs(self):
        """Test that transformation works between different graphs"""
        # For this test, we use the same graph but the function should work
        # with different graphs that have the same underlying structure
        # Use an actual variant edge from the graph: ('3_+', '4_+')
        variant_edge = ('3_+', '4_+')
        genotype_A = {variant_edge: 1}
        
        # Transform from G to G (simulating A to B)
        genotype_B = transform_genotype_via_walk(self.G, self.G, genotype_A)
        
        self.assertIsInstance(genotype_B, dict)
        self.assertGreater(len(genotype_B), 0)


class TestGenotypeTransformationEdgeCases(unittest.TestCase):
    """Test edge cases for genotype transformation"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.G = PangenomeGraph.from_gfa('tests/data/simple_nested.gfa', ref_name='ref')
    
    def test_transform_all_variants_raises_error(self):
        """Test that using all variants raises ValueError if they don't form valid walk"""
        genotype_A = {ve: 1 for ve in self.G.variant_edges}
        
        # All variants together likely don't form a valid walk
        # This should raise ValueError
        with self.assertRaises(ValueError):
            genotype_B = transform_genotype_via_walk(self.G, self.G, genotype_A)
    

    def test_empty_genotype_returns_empty(self):
        """Test that empty genotype returns empty result"""
        genotype_A = {}
        
        genotype_B = transform_genotype_via_walk(self.G, self.G, genotype_A)
        
        # Empty genotype should return empty result
        self.assertIsInstance(genotype_B, dict)
        self.assertEqual(len(genotype_B), 0)





if __name__ == '__main__':
    unittest.main()

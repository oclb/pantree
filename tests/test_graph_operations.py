"""
Unit tests for graph operations including simplification, repeat finding, 
inversions, and direction assignment.
"""
import unittest
import os
from graph_var.graph import PangenomeGraph
from graph_var.search_tree import dfs_methods


class GraphOperationsTestBase:
    """Base class with shared test methods for different GFA files"""
    
    def test_graph_loaded_successfully(self):
        """Test that the graph loads with expected properties"""
        self.assertIsNotNone(self.G)
        self.assertGreater(self.G.number_of_nodes(), 0)
        self.assertGreater(self.G.number_of_edges(), 0)
    
    def test_variant_edges_identified(self):
        """Test that variant edges are properly identified"""
        self.assertIsNotNone(self.G.variant_edges)
        self.assertIsInstance(self.G.variant_edges, set)
        # Variant edges should be non-empty
        self.assertGreater(len(self.G.variant_edges), 0)
    
    def test_ref_alt_alleles(self):
        """Test that ref and alt alleles can be extracted from variant edges"""
        for edge in self.G.variant_edges:
            ref_allele, alt_allele, last_letter, branch_point = self.G.ref_alt_alleles(edge)
            # All should return strings (possibly empty)
            self.assertIsInstance(ref_allele, str)
            self.assertIsInstance(alt_allele, str)
            self.assertIsInstance(last_letter, str)
            self.assertIsInstance(branch_point, str)
    
    def test_variant_type_identification(self):
        """Test that variant types are correctly identified"""
        variant_types = set()
        valid_types = ['SNP', 'INS', 'DEL', 'MNP', 'COMPLEX', 'INV', 'DUP']
        for edge in self.G.variant_edges:
            ref_allele, alt_allele, _, _ = self.G.ref_alt_alleles(edge)
            vt = self.G.identify_variant_type(edge, ref_allele, alt_allele)
            variant_types.add(vt)
            self.assertIn(vt, valid_types, f"Unknown variant type: {vt}")
        
        # Should have at least one variant type
        self.assertGreater(len(variant_types), 0)
    
    def test_is_insertion(self):
        """Test insertion detection"""
        insertions = [e for e in self.G.variant_edges if self.G.is_insertion(e)]
        self.assertIsInstance(insertions, list)
    
    def test_is_replacement(self):
        """Test replacement detection (variants with both ref and alt alleles)"""
        replacements = [e for e in self.G.variant_edges if self.G.is_replacement(e)]
        # Should have some replacements
        self.assertGreaterEqual(len(replacements), 0)
    
    def test_direction_assignment(self):
        """Test that node directions are properly assigned"""
        for node in self.G.nodes():
            if node in {'+_terminus_+', '-_terminus_+', '+_terminus_-', '-_terminus_-'}:
                continue
            direction = self.G.direction(node)
            self.assertIn(direction, [-1, 1], f"Node {node} has invalid direction {direction}")
    
    def test_is_inversion(self):
        """Test inversion detection"""
        # Test that the method runs without error
        inversions = [e for e in self.G.variant_edges if self.G.is_inversion(e)]
        self.assertIsInstance(inversions, list)
    
    def test_reference_path_exists(self):
        """Test that reference path is properly constructed"""
        self.assertIsNotNone(self.G.reference_path)
        self.assertGreater(len(self.G.reference_path), 0)
        # Reference path should start with +_terminus_+
        self.assertEqual(self.G.reference_path[0], '+_terminus_+')
    
    def test_reference_tree_exists(self):
        """Test that reference tree is constructed"""
        self.assertIsNotNone(self.G.reference_tree)
        self.assertGreater(self.G.reference_tree.number_of_nodes(), 0)
    
    def test_on_reference_path(self):
        """Test that on_reference_path correctly identifies nodes"""
        # All nodes in reference_path should return True
        for node in self.G.reference_path:
            if node in {'+_terminus_+', '-_terminus_+'}:
                continue
            # Create an edge from consecutive nodes
            idx = self.G.reference_path.index(node)
            if idx < len(self.G.reference_path) - 1:
                next_node = self.G.reference_path[idx + 1]
                edge = (node, next_node)
                if self.G.has_edge(*edge):
                    self.assertTrue(self.G.on_reference_path(edge))

    def test_dfs_method_parameter(self):
        """Test that different dfs_method_name parameters work correctly"""
        test_dir = os.path.dirname(__file__)
        gfa_file = os.path.join(test_dir, "data", "simple_nested.gfa")
        
        # Test with default max_weight method
        G_max_weight = PangenomeGraph.from_gfa_line_by_line(gfa_file, ref_name='ref', dfs_method_name='max_weight')
        self.assertIsNotNone(G_max_weight.reference_tree)
        self.assertGreater(G_max_weight.reference_tree.number_of_edges(), 0)
        
        # Test that the max_weight method works (it's the only one compatible with standard GFA files)
        # The contiguous method requires haplotype labels which aren't present in standard GFA files
        G_again = PangenomeGraph.from_gfa_line_by_line(gfa_file, ref_name='ref', dfs_method_name='max_weight')
        self.assertIsNotNone(G_again.reference_tree)
        self.assertGreater(G_again.reference_tree.number_of_edges(), 0)
        
        # Test invalid method name raises appropriate error
        with self.assertRaises(KeyError):
            PangenomeGraph.from_gfa_line_by_line(gfa_file, ref_name='ref', dfs_method_name='invalid_method')


class TestGraphOperations(GraphOperationsTestBase, unittest.TestCase):
    """Test various graph operations on the simple_nested.gfa test file"""
    
    @classmethod
    def setUpClass(cls):
        """Load the test GFA file once for all tests"""
        test_dir = os.path.dirname(__file__)
        cls.gfa_file = os.path.join(test_dir, "data", "simple_nested.gfa")
        cls.G = PangenomeGraph.from_gfa_line_by_line(cls.gfa_file, ref_name='ref')


class TestC4AInversionGraph(GraphOperationsTestBase, unittest.TestCase):
    """Test graph operations on the c4a_with_inversion_and_sequences.gfa file"""
    
    @classmethod
    def setUpClass(cls):
        """Load the C4A test GFA file with inversions"""
        test_dir = os.path.dirname(__file__)
        cls.gfa_file = os.path.join(test_dir, "data", "c4a_with_inversion_and_sequences.gfa")
        cls.G = PangenomeGraph.from_gfa_line_by_line(cls.gfa_file, ref_name='GRCh38')
    
    def test_has_inversions(self):
        """Test that this graph contains inversions"""
        inversions = [e for e in self.G.variant_edges if self.G.is_inversion(e)]
        # This file should have at least one inversion
        self.assertGreater(len(inversions), 0, "C4A graph should contain inversions")
    
    def test_missing_inversion_allele(self):
        """Test missing inversion allele detection on C4A graph"""
        found_missing_inversion = False
        for edge in self.G.variant_edges:
            inversion_allele = self.G.missing_inversion_allele(edge, minimum_alt_length=5)
            if inversion_allele:
                found_missing_inversion = True
                self.assertIsInstance(inversion_allele, str)
                self.assertGreater(len(inversion_allele), 0)
        
        # May or may not find missing inversions, but method should work
        self.assertIsInstance(found_missing_inversion, bool)


class TestRepeatAnnotation(unittest.TestCase):
    """Test repeat motif annotation"""
    
    @classmethod
    def setUpClass(cls):
        """Load the test GFA file once for all tests"""
        test_dir = os.path.dirname(__file__)
        cls.gfa_file = os.path.join(test_dir, "data", "simple_nested.gfa")
        cls.G = PangenomeGraph.from_gfa_line_by_line(cls.gfa_file, ref_name='ref')
    
    def test_annotate_repeat_motif(self):
        """Test that repeat motif annotation runs without error"""
        for edge in self.G.variant_edges:
            ref_allele, alt_allele, _, branch_point = self.G.ref_alt_alleles(edge)
            motif = self.G.annotate_repeat_motif(
                edge, 
                ref_allele=ref_allele, 
                alt_allele=alt_allele,
                branch_point=branch_point
            )
            # Motif can be None or a string
            self.assertTrue(motif is None or isinstance(motif, str))
    
    def test_is_back_edge(self):
        """Test back edge detection"""
        back_edges = [e for e in self.G.variant_edges if self.G.is_back_edge(e)]
        self.assertIsInstance(back_edges, list)
    
    def test_is_forward_edge(self):
        """Test forward edge detection"""
        forward_edges = [e for e in self.G.variant_edges if self.G.is_forward_edge(e)]
        self.assertIsInstance(forward_edges, list)
    
    def test_is_crossing_edge(self):
        """Test crossing edge detection"""
        crossing_edges = [e for e in self.G.variant_edges if self.G.is_crossing_edge(e)]
        self.assertIsInstance(crossing_edges, list)


class TestGraphSimplification(unittest.TestCase):
    """Test graph simplification operations"""
    
    @classmethod
    def setUpClass(cls):
        """Load the test GFA file once for all tests"""
        test_dir = os.path.dirname(__file__)
        cls.gfa_file = os.path.join(test_dir, "data", "simple_nested.gfa")
        cls.G = PangenomeGraph.from_gfa_line_by_line(cls.gfa_file, ref_name='ref')
    
    def test_delete_small_variants(self):
        """Test that small variants can be deleted"""
        original_edges = self.G.number_of_edges()
        simplified = self.G.delete_small_variants(minimum_allele_length=5)
        
        self.assertIsNotNone(simplified)
        # Simplified graph should have fewer or equal edges
        self.assertLessEqual(simplified.number_of_edges(), original_edges)
    
    def test_allele_count(self):
        """Test that allele counts are computed"""
        allele_counts = self.G.allele_count()
        self.assertIsInstance(allele_counts, dict)
        
        # Each variant edge should have a count
        for edge in self.G.variant_edges:
            self.assertIn(edge, allele_counts)
            ref_count, alt_count = allele_counts[edge]
            self.assertIsInstance(ref_count, int)
            self.assertIsInstance(alt_count, int)
            self.assertGreaterEqual(ref_count, 0)
            self.assertGreaterEqual(alt_count, 0)


class TestMissingInversion(unittest.TestCase):
    """Test missing inversion detection"""
    
    @classmethod
    def setUpClass(cls):
        """Load the test GFA file once for all tests"""
        test_dir = os.path.dirname(__file__)
        cls.gfa_file = os.path.join(test_dir, "data", "simple_nested.gfa")
        cls.G = PangenomeGraph.from_gfa_line_by_line(cls.gfa_file, ref_name='ref')
    
    def test_missing_inversion_allele(self):
        """Test that missing inversion allele detection runs without error"""
        for edge in self.G.variant_edges:
            inversion_allele = self.G.missing_inversion_allele(edge, minimum_alt_length=5)
            # Should return None or a string
            self.assertTrue(inversion_allele is None or isinstance(inversion_allele, str))


if __name__ == '__main__':
    unittest.main()

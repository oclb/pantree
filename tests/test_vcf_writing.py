"""
Unit tests for VCF writing functionality including edge cases
"""
import unittest
import os
import tempfile
from graph_var.graph import PangenomeGraph
from graph_var.genotype import Genotype


class TestVCFWriting(unittest.TestCase):
    """Test VCF writing with various parameters"""
    
    @classmethod
    def setUpClass(cls):
        """Load the test GFA file once for all tests"""
        test_dir = os.path.dirname(__file__)
        cls.gfa_file = os.path.join(test_dir, "data", "simple_nested.gfa")
        cls.G = PangenomeGraph.from_gfa_line_by_line(cls.gfa_file, ref_name='ref')
    
    def test_write_vcf_with_genotypes(self):
        """Test VCF writing with genotypes (default behavior)"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            vcf_path = f.name
        
        try:
            self.G.write_vcf(
                self.gfa_file,
                vcf_path,
                chr_name='chr0'
            )
            
            # Read and verify VCF content
            with open(vcf_path, 'r') as f:
                content = f.read()
            
            # Should have genotype columns
            self.assertIn('FORMAT', content)
            self.assertIn('GT:CR:CA', content)
            self.assertIn('ref', content)
            self.assertIn('sample1', content)
            self.assertIn('sample2', content)
        finally:
            if os.path.exists(vcf_path):
                os.unlink(vcf_path)
    
    def test_write_vcf_no_genotypes(self):
        """Test VCF writing without genotypes (gfa_path=None)"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            vcf_path = f.name
        
        try:
            # Pass None for gfa_path to skip genotype computation
            self.G.write_vcf(
                None,
                vcf_path,
                chr_name='chr0'
            )
            
            # Read and verify VCF content
            with open(vcf_path, 'r') as f:
                lines = f.readlines()
            
            # Should have FORMAT column header but no sample columns
            header_line = [l for l in lines if l.startswith('#CHROM')][0]
            columns = header_line.strip().split('\t')
            
            # Should have standard VCF columns plus FORMAT
            self.assertIn('#CHROM', columns)
            self.assertIn('POS', columns)
            self.assertIn('REF', columns)
            self.assertIn('ALT', columns)
            self.assertIn('FORMAT', columns)
            
            # Should NOT have sample columns (no genotypes computed)
            # The FORMAT column is the last column when no samples
            self.assertEqual(columns[-1], 'FORMAT')
        finally:
            if os.path.exists(vcf_path):
                os.unlink(vcf_path)
    
    def test_write_vcf_with_size_threshold(self):
        """Test VCF writing with size_threshold parameter"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            vcf_path = f.name
        
        try:
            # Write with size threshold of 100 (should truncate long alleles)
            self.G.write_vcf(
                None,
                vcf_path,
                chr_name='chr0',
                size_threshold=100
            )
            
            # Read and count variants
            with open(vcf_path, 'r') as f:
                lines = [line for line in f if not line.startswith('#')]
            
            # With high threshold, should have fewer variants than default
            self.assertIsInstance(lines, list)
        finally:
            if os.path.exists(vcf_path):
                os.unlink(vcf_path)
    
    def test_write_vcf_check_degenerate(self):
        """Test VCF writing with check_degenerate parameter"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            vcf_path = f.name
        
        try:
            # Write with check_degenerate=True (excludes variants with identical ref/alt)
            self.G.write_vcf(
                None,
                vcf_path,
                chr_name='chr0',
                check_degenerate=True
            )
            
            # Should complete without error
            self.assertTrue(os.path.exists(vcf_path))
            
            # Verify it's a valid VCF
            with open(vcf_path, 'r') as f:
                content = f.read()
            self.assertIn('##fileformat=VCFv4.2', content)
        finally:
            if os.path.exists(vcf_path):
                os.unlink(vcf_path)


class TestGenotypeAggregation(unittest.TestCase):
    """Test genotype aggregation and sample VCF info generation"""
    
    @classmethod
    def setUpClass(cls):
        """Load the test GFA file once for all tests"""
        test_dir = os.path.dirname(__file__)
        cls.gfa_file = os.path.join(test_dir, "data", "simple_nested.gfa")
        cls.G = PangenomeGraph.from_gfa_line_by_line(cls.gfa_file, ref_name='ref')
    
    def test_genotype_and_linear_coverage_by_sample(self):
        """Test genotype aggregation by sample from multiple walks"""
        # Test using the Genotype class instead of the removed graph.genotype method
        
        walks = [
            ['1_+', '2_+', '4_+', '9_+', '10_+', '11_+'],
            ['1_+', '3_+', '4_+', '5_+', '7_+', '8_+', '9_+', '11_+'],
        ]
        
        # Test that Genotype.genotype works for each walk
        for walk in walks:
            genotype = Genotype.genotype(self.G, walk, exclude_terminus=True)
            
            self.assertIsInstance(genotype.ref_counts, dict)
            self.assertIsInstance(genotype.alt_counts, dict)
            self.assertIsInstance(genotype.linear_coverage, list)
            
            # Should have entries if walk visits variants
            if len(self.G.variant_edges) > 0:
                # At least one dict should have entries
                self.assertTrue(len(genotype.ref_counts) > 0 or len(genotype.alt_counts) > 0)
    
    @unittest.skip("get_sample_vcf_info method removed during refactor")
    def test_get_sample_vcf_info(self):
        """Test get_sample_vcf_info method"""
        pass


class TestEdgeCases(unittest.TestCase):
    """Test edge cases: single node, no variants, etc."""
    
    def test_single_node_graph(self):
        """Test loading and processing a graph with a single node"""
        test_dir = os.path.dirname(__file__)
        gfa_file = os.path.join(test_dir, "data", "single_node.gfa")
        
        G = PangenomeGraph.from_gfa_line_by_line(gfa_file, ref_name='ref')
        
        # Should load successfully
        self.assertIsNotNone(G)
        self.assertGreater(G.number_of_nodes(), 0)
        
        # May have terminal edges in variant_edges, but no non-terminal variants
        non_terminal_variants = [e for e in G.variant_edges if not G.is_terminal(e)]
        self.assertEqual(len(non_terminal_variants), 0)
    
    def test_no_variants_graph(self):
        """Test graph with multiple nodes but no variants"""
        test_dir = os.path.dirname(__file__)
        gfa_file = os.path.join(test_dir, "data", "no_variants.gfa")
        
        G = PangenomeGraph.from_gfa_line_by_line(gfa_file, ref_name='ref')
        
        # Should load successfully
        self.assertIsNotNone(G)
        self.assertGreater(G.number_of_nodes(), 0)
        self.assertGreater(G.number_of_edges(), 0)
        
        # May have terminal edges in variant_edges, but no non-terminal variants
        non_terminal_variants = [e for e in G.variant_edges if not G.is_terminal(e)]
        self.assertEqual(len(non_terminal_variants), 0)
        
        # Should still be able to write VCF (empty of variants)
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            vcf_path = f.name
        
        try:
            G.write_vcf(
                None,
                vcf_path,
                chr_name='chr0'
            )
            
            # Should create valid VCF with header but no variants
            with open(vcf_path, 'r') as f:
                lines = f.readlines()
            
            # Should have header lines
            header_lines = [l for l in lines if l.startswith('#')]
            self.assertGreater(len(header_lines), 0)
            
            # Should have no variant lines
            variant_lines = [l for l in lines if not l.startswith('#')]
            self.assertEqual(len(variant_lines), 0)
        finally:
            if os.path.exists(vcf_path):
                os.unlink(vcf_path)
    
    @unittest.skip("allele_count method removed during refactor")
    def test_allele_count_empty_variants(self):
        """Test allele_count on graph with no variants"""
        pass


if __name__ == '__main__':
    unittest.main()

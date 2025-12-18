"""
Test for GFA1.0 support (P-lines instead of W-lines).
This tests that pantree can handle both GFA1.0 and GFA1.1 formats.
"""

import unittest
import os
import tempfile
from pantree import PangenomeGraph


class TestGFA10Support(unittest.TestCase):
    """Test that GFA1.0 files with P-lines are handled correctly."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_data_dir = os.path.join(os.path.dirname(__file__), 'data')
        self.gfa10_file = os.path.join(self.test_data_dir, 'simple_gfa10.gfa')
        self.gfa11_file = os.path.join(self.test_data_dir, 'simple_nested.gfa')
    
    def test_load_gfa10_graph(self):
        """Test that GFA1.0 file can be loaded into a graph."""
        G = PangenomeGraph.from_gfa_line_by_line(
            self.gfa10_file,
            ref_name='ref'
        )
        
        # Verify graph was constructed properly
        self.assertGreater(len(G.nodes), 0, "Graph should have nodes")
        self.assertGreater(len(G.edges), 0, "Graph should have edges")
        self.assertIsNotNone(G.reference_path, "Reference path should be set")
        self.assertGreater(len(G.reference_path), 0, "Reference path should not be empty")
    
    def test_write_vcf_from_gfa10(self):
        """Test that VCF can be written from GFA1.0 file with genotypes.
        
        This reproduces the issue where sample_name is undefined when
        processing GFA1.0 files because P-lines are skipped.
        """
        G = PangenomeGraph.from_gfa_line_by_line(
            self.gfa10_file,
            ref_name='ref'
        )
        
        # Try to write VCF with genotypes - this should not raise an error
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            vcf_path = f.name
        
        try:
            # This should work without errors
            G.write_vcf(self.gfa10_file, vcf_path, 'chr1')
            
            # Verify VCF was created
            self.assertTrue(os.path.exists(vcf_path), "VCF file should be created")
            
            # Verify VCF has content
            with open(vcf_path, 'r') as f:
                content = f.read()
                self.assertIn('##fileformat=VCFv4.2', content, "VCF should have proper header")
                self.assertIn('#CHROM', content, "VCF should have column headers")
        finally:
            # Clean up
            if os.path.exists(vcf_path):
                os.remove(vcf_path)
    
    def test_gfa10_and_gfa11_produce_similar_graphs(self):
        """Test that GFA1.0 and GFA1.1 produce graphs with similar structure."""
        G10 = PangenomeGraph.from_gfa_line_by_line(
            self.gfa10_file,
            ref_name='ref'
        )
        
        G11 = PangenomeGraph.from_gfa_line_by_line(
            self.gfa11_file,
            ref_name='ref'
        )
        
        # Both should have the same number of nodes
        self.assertEqual(len(G10.nodes), len(G11.nodes),
                        "Both graphs should have the same number of nodes")
        
        # Edge counts may differ slightly due to terminal edge handling differences
        # between GFA 1.0 and 1.1 formats after the DFS bug fix
        # The difference should only be in terminal edges (at most 2 edges)
        edge_diff = abs(len(G10.edges) - len(G11.edges))
        self.assertLessEqual(edge_diff, 2,
                            f"Edge count difference ({edge_diff}) is too large. "
                            f"GFA 1.0: {len(G10.edges)}, GFA 1.1: {len(G11.edges)}")


if __name__ == '__main__':
    unittest.main()

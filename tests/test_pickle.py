"""
Unit tests for pickle save/load functions.
"""
import unittest
import os
import tempfile
from graph_var.utils import save_graph_to_pkl, load_graph_from_pkl
from graph_var.graph import PangenomeGraph


class TestPickleFunctions(unittest.TestCase):
    """Test graph pickle save and load functions"""
    
    @classmethod
    def setUpClass(cls):
        """Load a test graph once for all tests"""
        test_dir = os.path.dirname(__file__)
        gfa_file = os.path.join(test_dir, "data", "simple_nested.gfa")
        cls.G = PangenomeGraph.from_gfa_line_by_line(gfa_file, ref_name='ref')
    
    def test_save_and_load_uncompressed(self):
        """Test saving and loading graph without compression"""
        with tempfile.NamedTemporaryFile(suffix='.pkl', delete=False) as f:
            pkl_path = f.name
        
        try:
            # Save graph
            save_graph_to_pkl(self.G, pkl_path, compressed=False)
            self.assertTrue(os.path.exists(pkl_path))
            
            # Load graph
            loaded_G = load_graph_from_pkl(pkl_path, compressed=False)
            
            # Verify basic properties match
            self.assertEqual(loaded_G.number_of_nodes(), self.G.number_of_nodes())
            self.assertEqual(loaded_G.number_of_edges(), self.G.number_of_edges())
            self.assertEqual(len(loaded_G.variant_edges), len(self.G.variant_edges))
            
        finally:
            if os.path.exists(pkl_path):
                os.remove(pkl_path)
    
    def test_save_and_load_compressed(self):
        """Test saving and loading graph with gzip compression"""
        with tempfile.NamedTemporaryFile(suffix='.pkl.gz', delete=False) as f:
            pkl_path = f.name
        
        try:
            # Save graph compressed
            save_graph_to_pkl(self.G, pkl_path, compressed=True)
            self.assertTrue(os.path.exists(pkl_path))
            
            # Load graph compressed
            loaded_G = load_graph_from_pkl(pkl_path, compressed=True)
            
            # Verify basic properties match
            self.assertEqual(loaded_G.number_of_nodes(), self.G.number_of_nodes())
            self.assertEqual(loaded_G.number_of_edges(), self.G.number_of_edges())
            self.assertEqual(len(loaded_G.variant_edges), len(self.G.variant_edges))
            
        finally:
            if os.path.exists(pkl_path):
                os.remove(pkl_path)
    
    def test_compressed_file_smaller(self):
        """Test that compressed files are smaller than uncompressed"""
        with tempfile.NamedTemporaryFile(suffix='.pkl', delete=False) as f:
            uncompressed_path = f.name
        with tempfile.NamedTemporaryFile(suffix='.pkl.gz', delete=False) as f:
            compressed_path = f.name
        
        try:
            # Save both versions
            save_graph_to_pkl(self.G, uncompressed_path, compressed=False)
            save_graph_to_pkl(self.G, compressed_path, compressed=True)
            
            # Check file sizes
            uncompressed_size = os.path.getsize(uncompressed_path)
            compressed_size = os.path.getsize(compressed_path)
            
            # Compressed should be smaller
            self.assertLess(compressed_size, uncompressed_size)
            
        finally:
            if os.path.exists(uncompressed_path):
                os.remove(uncompressed_path)
            if os.path.exists(compressed_path):
                os.remove(compressed_path)


if __name__ == '__main__':
    unittest.main()

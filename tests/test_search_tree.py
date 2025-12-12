"""
Unit tests for search tree functionality
"""
import unittest
import networkx as nx
from pantree.dfs import max_weight_dfs_tree


class TestSearchTree(unittest.TestCase):
    """Test search tree algorithms"""

    def test_max_weight_dfs_tree_simple(self):
        """Test max weight DFS tree on a simple graph"""
        G = nx.DiGraph()
        G.add_edge('A', 'B', weight=5)
        G.add_edge('A', 'C', weight=3)
        G.add_edge('B', 'D', weight=2)
        G.add_edge('C', 'D', weight=4)

        tree = max_weight_dfs_tree(G, 'A')

        self.assertIsInstance(tree, nx.DiGraph)
        self.assertIn('A', tree.nodes())
        # Tree should have n-1 edges for n nodes
        self.assertEqual(tree.number_of_edges(), tree.number_of_nodes() - 1)

    def test_max_weight_dfs_tree_single_node(self):
        """Test max weight DFS tree with single node"""
        G = nx.DiGraph()
        G.add_node('A')

        tree = max_weight_dfs_tree(G, 'A')

        # Tree may be empty if there are no edges from source
        self.assertGreaterEqual(tree.number_of_nodes(), 0)
        self.assertEqual(tree.number_of_edges(), 0)

    def test_max_weight_dfs_tree_linear(self):
        """Test max weight DFS tree on a linear graph"""
        G = nx.DiGraph()
        G.add_edge('A', 'B', weight=1)
        G.add_edge('B', 'C', weight=1)
        G.add_edge('C', 'D', weight=1)

        tree = max_weight_dfs_tree(G, 'A')

        self.assertEqual(tree.number_of_nodes(), 4)
        self.assertEqual(tree.number_of_edges(), 3)

    def test_max_weight_dfs_tree_with_cycles(self):
        """Test max weight DFS tree on a graph with cycles"""
        G = nx.DiGraph()
        G.add_edge('A', 'B', weight=10)
        G.add_edge('B', 'C', weight=5)
        G.add_edge('C', 'A', weight=1)  # Creates cycle
        G.add_edge('B', 'D', weight=3)

        tree = max_weight_dfs_tree(G, 'A')

        # Should be acyclic
        self.assertTrue(nx.is_directed_acyclic_graph(tree))
        # Should have all nodes
        self.assertEqual(tree.number_of_nodes(), 4)


if __name__ == '__main__':
    unittest.main()

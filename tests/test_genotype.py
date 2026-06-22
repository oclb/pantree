"""
Unit tests for genotype-related functionality
"""
import unittest
import os
import json
import subprocess
import sys
import textwrap
from pantree.graph import PangenomeGraph
from pantree.genotype import Genotype
from pantree.gfa import read_gfa_line_by_line


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

    def test_nodes_have_tree_position(self):
        """Test that nodes have tree_position attributes"""
        for node in self.G.nodes():
            if self.G.is_terminal(node):
                continue
            self.assertIn('tree_position', self.G.nodes[node])
            self.assertIsInstance(self.G.nodes[node]['tree_position'], (int, float))

    @unittest.skip("get_vcf_position method removed - now in vcf.py")
    def test_get_vcf_position(self):
        """Test VCF position calculation"""
        pass


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


class TestGenotypeClass(unittest.TestCase):
    """Test the Genotype dataclass and its methods"""

    @classmethod
    def setUpClass(cls):
        """Load the test GFA file once for all tests"""
        test_dir = os.path.dirname(__file__)
        cls.gfa_file = os.path.join(test_dir, "data", "simple_nested.gfa")
        cls.G = PangenomeGraph.from_gfa_line_by_line(cls.gfa_file, ref_name='ref')

    def test_genotype_creation(self):
        """Test creating a Genotype from a walk"""
        walk = ['1_+', '2_+', '4_+', '9_+', '10_+', '11_+']
        genotype = Genotype.genotype(self.G, walk, exclude_terminus=True)

        # Check that genotype is a Genotype instance
        self.assertIsInstance(genotype, Genotype)

        # Check that it has the expected attributes
        self.assertIsInstance(genotype.ref_counts, dict)
        self.assertIsInstance(genotype.alt_counts, dict)
        self.assertIsInstance(genotype.linear_coverage, list)
        self.assertEqual(genotype.exclude_terminus, True)
        self.assertIsNone(genotype.missing_variants)

    def test_genotype_linear_coverage(self):
        """Test that linear coverage is computed correctly"""
        walk = ['1_+', '2_+', '4_+', '9_+', '10_+', '11_+']
        genotype = Genotype.genotype(self.G, walk, exclude_terminus=True)

        # Should have exactly one coverage interval
        self.assertEqual(len(genotype.linear_coverage), 1)
        min_pos, max_pos = genotype.linear_coverage[0]

        # Min should be less than max
        self.assertLessEqual(min_pos, max_pos)

        # Both should be non-negative
        self.assertGreaterEqual(min_pos, 0)
        self.assertGreaterEqual(max_pos, 0)

    def test_genotype_update(self):
        """Test updating a genotype with another genotype"""
        walk1 = ['1_+', '2_+', '4_+', '9_+', '10_+', '11_+']
        walk2 = ['1_+', '3_+', '4_+', '9_+', '11_+']

        genotype1 = Genotype.genotype(self.G, walk1, exclude_terminus=True)
        genotype2 = Genotype.genotype(self.G, walk2, exclude_terminus=True)

        # Store original counts
        orig_ref_count = sum(genotype1.ref_counts.values())
        orig_alt_count = sum(genotype1.alt_counts.values())
        orig_coverage_len = len(genotype1.linear_coverage)

        # Update genotype1 with genotype2
        genotype1.update(genotype2)

        # Check that counts increased
        new_ref_count = sum(genotype1.ref_counts.values())
        new_alt_count = sum(genotype1.alt_counts.values())
        new_coverage_len = len(genotype1.linear_coverage)

        # At least one should have increased
        self.assertTrue(
            new_ref_count >= orig_ref_count or
            new_alt_count >= orig_alt_count
        )

        # Coverage should have been appended
        self.assertEqual(new_coverage_len, orig_coverage_len + 1)

    def test_compute_missing_variants(self):
        """Test computing missing variants for a genotype"""
        walk = ['1_+', '2_+', '4_+', '9_+', '10_+', '11_+']
        genotype = Genotype.genotype(self.G, walk, exclude_terminus=True)

        # Initially missing_variants should be None
        self.assertIsNone(genotype.missing_variants)

        # Compute missing variants
        genotype.compute_missing_variants(self.G)

        # Now it should be a set
        self.assertIsInstance(genotype.missing_variants, set)

        # All elements should be tuples (edges)
        for variant in genotype.missing_variants:
            self.assertIsInstance(variant, tuple)
            self.assertEqual(len(variant), 2)

    def test_variant_record(self):
        """Test getting variant record for an edge"""
        walk = ['1_+', '2_+', '4_+', '9_+', '10_+', '11_+']
        genotype = Genotype.genotype(self.G, walk, exclude_terminus=True)
        genotype.compute_missing_variants(self.G)

        # Get a variant edge and its reference edge
        variant_edges = list(self.G.sorted_variant_edges(exclude_terminus=True))
        if variant_edges:
            edge = variant_edges[0]
            reference_edge = self.G.reference_tree_edge(edge)
            gt, cr, ca = genotype.variant_record(edge, reference_edge)

            # GT should be None or int
            self.assertTrue(gt is None or isinstance(gt, int))

            # CR and CA should be integers
            self.assertIsInstance(cr, int)
            self.assertIsInstance(ca, int)

            # Both should be non-negative
            self.assertGreaterEqual(cr, 0)
            self.assertGreaterEqual(ca, 0)


class TestGenotypesFromGFA(unittest.TestCase):
    """Test the genotypes_from_gfa method"""

    @classmethod
    def setUpClass(cls):
        """Load the test GFA file once for all tests"""
        test_dir = os.path.dirname(__file__)
        cls.gfa_file = os.path.join(test_dir, "data", "simple_nested.gfa")
        cls.G = PangenomeGraph.from_gfa_line_by_line(cls.gfa_file, ref_name='ref')

    def test_genotypes_from_gfa_returns_dict(self):
        """Test that genotypes_from_gfa returns a dictionary"""
        result = self.G.genotypes_from_gfa(self.gfa_file, exclude_terminus=True)

        # Should return a dictionary
        self.assertIsInstance(result, dict)

    def test_genotypes_from_gfa_has_samples(self):
        """Test that genotypes_from_gfa returns genotypes for each sample"""
        result = self.G.genotypes_from_gfa(self.gfa_file, exclude_terminus=True)

        # Should have sample1 and sample2 based on the GFA file
        self.assertIn('sample1', result)
        self.assertIn('sample2', result)

    def test_genotypes_from_gfa_tuple_structure(self):
        """Test that each sample has a tuple of genotypes"""
        result = self.G.genotypes_from_gfa(self.gfa_file, exclude_terminus=True)

        for sample_name, genotypes in result.items():
            # Each value should be a tuple
            self.assertIsInstance(genotypes, tuple)

            # Should have 1 or 2 haplotypes
            self.assertIn(len(genotypes), (1, 2))

            # Each genotype should be a Genotype instance
            for genotype in genotypes:
                self.assertIsInstance(genotype, Genotype)

    def test_genotypes_from_gfa_has_missing_variants(self):
        """Test that genotypes have missing_variants computed"""
        result = self.G.genotypes_from_gfa(self.gfa_file, exclude_terminus=True)

        for sample_name, genotypes in result.items():
            for genotype in genotypes:
                # missing_variants should be computed (not None)
                self.assertIsNotNone(genotype.missing_variants)
                self.assertIsInstance(genotype.missing_variants, set)

    def test_genotypes_from_gfa_diploid_samples(self):
        """Test that diploid samples have 2 haplotypes"""
        result = self.G.genotypes_from_gfa(self.gfa_file, exclude_terminus=True)

        # Based on simple_nested.gfa, both samples should be diploid
        for sample_name in ['sample1', 'sample2']:
            if sample_name in result:
                genotypes = result[sample_name]
                self.assertEqual(len(genotypes), 2,
                               f"{sample_name} should have 2 haplotypes")

    def test_genotypes_from_gfa_orders_haplotypes_by_gfa_index(self):
        """Test that diploid haplotypes are ordered by explicit GFA haplotype index"""
        result = self.G.genotypes_from_gfa(self.gfa_file, exclude_terminus=True)
        sample2_lines = [
            line for line in read_gfa_line_by_line(self.gfa_file, line_types=['W', 'P'])
            if line.sample_name == 'sample2'
        ]
        self.assertEqual(
            ['#'.join(line.hap_name.split('#')[:2]) for line in sample2_lines],
            ['sample2#1', 'sample2#2'],
        )

        expected_genotypes = []
        for line in sample2_lines:
            genotype = self.G.genotype(line.walk, exclude_terminus=True)
            genotype.compute_missing_variants(self.G)
            expected_genotypes.append(genotype)

        def genotype_signature(genotype):
            signature = []
            for edge in self.G.sorted_variant_edges(exclude_terminus=True):
                reference_edge = self.G.representative_edge(self.G.reference_tree_edge(edge))
                ref_edge_idx = int(self.G.edges[reference_edge]['index'])
                signature.append(genotype.variant_record(edge, reference_edge, ref_edge_idx))
            return signature

        self.assertEqual(
            [genotype_signature(genotype) for genotype in result['sample2']],
            [genotype_signature(genotype) for genotype in expected_genotypes],
        )

    def test_haplotype_tuple_stable_across_python_hash_seeds(self):
        """Test that phased haplotype tuple order is deterministic across hash seeds"""
        code = textwrap.dedent(f"""
            import json
            from pantree.graph import PangenomeGraph

            gfa_file = {self.gfa_file!r}
            graph = PangenomeGraph.from_gfa_line_by_line(gfa_file, ref_name='ref')
            genotypes = graph.genotypes_from_gfa(gfa_file, exclude_terminus=True)
            variant_edge = ('7_+', '8_+')
            reference_edge = ('6_+', '8_+')
            ref_edge_idx = int(graph.edges[reference_edge]['index'])
            records = [
                genotype.variant_record(variant_edge, reference_edge, ref_edge_idx)
                for genotype in genotypes['sample2']
            ]
            print(json.dumps(records))
        """)

        outputs = []
        for seed in ('1', '2', '3'):
            env = os.environ.copy()
            env['PYTHONHASHSEED'] = seed
            completed = subprocess.run(
                [sys.executable, '-c', code],
                check=True,
                capture_output=True,
                text=True,
                env=env,
            )
            outputs.append(json.loads(completed.stdout))

        self.assertEqual(outputs[0], outputs[1])
        self.assertEqual(outputs[0], outputs[2])
        self.assertEqual(outputs[0], [[1, 0, 1], [0, 1, 0]])


class TestGenotypeMatchesWalk(unittest.TestCase):
    """Test that genotypes correctly match the walks they were computed from"""

    def test_verify_genotype_matches_walk_simple_nested(self):
        """Test verification on simple_nested.gfa"""
        test_dir = os.path.dirname(__file__)
        gfa_file = os.path.join(test_dir, "data", "simple_nested.gfa")
        G = PangenomeGraph.from_gfa(gfa_file, ref_name='ref')

        # Test with a simple walk
        walk = ['1_+', '2_+', '4_+', '9_+', '10_+', '11_+']
        genotype = Genotype.genotype(G, walk, exclude_terminus=True)

        is_valid, errors = Genotype.verify_genotype_matches_walk(G, walk, genotype)
        self.assertTrue(is_valid, f"Genotype should match walk. Errors: {errors}")

    def test_verify_all_haplotypes_c4a_with_inversion(self):
        """Test that all haplotypes in c4a_with_inversion_and_sequences.gfa have matching genotypes"""
        test_dir = os.path.dirname(__file__)
        gfa_file = os.path.join(test_dir, "data", "c4a_with_inversion_and_sequences.gfa")

        if not os.path.exists(gfa_file):
            self.skipTest(f"GFA file not found: {gfa_file}")

        G = PangenomeGraph.from_gfa(gfa_file, ref_name='ref', return_walks=True)

        from collections import defaultdict
        from pantree.gfa import read_gfa_line_by_line

        haplotype_walks = defaultdict(list)
        for line in read_gfa_line_by_line(gfa_file, line_types=['W', 'P']):
            haplotype_name = '#'.join(line.hap_name.split('#')[:2])
            haplotype_walks[haplotype_name].append(line.walk)

        for hap_name, walks in haplotype_walks.items():
            if walks:
                # Compute genotype directly from walks
                genotype = Genotype.genotype(G, walks[0], exclude_terminus=True)
                for walk in walks[1:]:
                    genotype.update(Genotype.genotype(G, walk, exclude_terminus=True))

                is_valid, errors = Genotype.verify_genotype_matches_walks(G, walks, genotype)
                self.assertTrue(
                    is_valid,
                    f"Genotype for {hap_name} should match walks. Errors: {errors[:5]}"
                )



if __name__ == '__main__':
    unittest.main()

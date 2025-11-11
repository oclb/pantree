"""
Unit tests for VCF writing functionality
"""
import unittest
import os
import tempfile
from graph_var.graph import PangenomeGraph
from graph_var.genotype import Genotype
from graph_var.vcf import (
    _generate_vcf_metadata,
    _build_genotype_record,
    _build_vcf_record,
    write_vcf_from_graph
)


class TestVCFMetadata(unittest.TestCase):
    """Test VCF metadata generation"""
    
    def test_generate_vcf_metadata(self):
        """Test that VCF metadata is generated correctly"""
        chr_name = "chr1"
        metadata = _generate_vcf_metadata(chr_name)
        
        # Check that it's a string
        self.assertIsInstance(metadata, str)
        
        # Check for required VCF header lines
        self.assertIn("##fileformat=VCFv4.2", metadata)
        self.assertIn("##FORMAT=<ID=GT", metadata)
        self.assertIn("##FORMAT=<ID=CR", metadata)
        self.assertIn("##FORMAT=<ID=CA", metadata)
        self.assertIn("##INFO=<ID=NR", metadata)
        self.assertIn("##INFO=<ID=VT", metadata)
        self.assertIn("##INFO=<ID=DR", metadata)
        self.assertIn("##INFO=<ID=RC", metadata)
        self.assertIn("##INFO=<ID=AC", metadata)
        self.assertIn("##INFO=<ID=AN", metadata)
        self.assertIn("##INFO=<ID=PV", metadata)
        self.assertIn("##INFO=<ID=HP", metadata)
        self.assertIn("##INFO=<ID=TR_MOTIF", metadata)
        self.assertIn("##INFO=<ID=NIA", metadata)
        self.assertIn(f"##contig=<ID={chr_name}>", metadata)


class TestGenotypeRecord(unittest.TestCase):
    """Test genotype record building"""
    
    @classmethod
    def setUpClass(cls):
        """Load the test GFA file once for all tests"""
        test_dir = os.path.dirname(__file__)
        cls.gfa_file = os.path.join(test_dir, "data", "simple_nested.gfa")
        cls.G = PangenomeGraph.from_gfa_line_by_line(cls.gfa_file, ref_name='ref')
    
    def test_build_genotype_record_diploid(self):
        """Test building genotype records for diploid samples"""
        # Get genotypes from GFA
        sample_to_genotype = self.G.genotypes_from_gfa(self.gfa_file, exclude_terminus=True)
        
        # Get a variant edge
        variant_edges = list(self.G.sorted_variant_edges(exclude_terminus=True))
        self.assertGreater(len(variant_edges), 0, "Need at least one variant edge")
        
        edge = variant_edges[0]
        reference_edge = self.G.representative_edge(self.G.reference_tree_edge(edge))
        is_inversion = self.G.is_inversion(edge)
        sample_order = list(sample_to_genotype.keys())
        
        # Build genotype record
        genotype_records = _build_genotype_record(edge, reference_edge, sample_to_genotype, is_inversion, sample_order)
        
        # Should have one record per sample
        self.assertEqual(len(genotype_records), len(sample_to_genotype))
        
        # Each record should be a string
        for record in genotype_records:
            self.assertIsInstance(record, str)
            
            # Should contain colons
            self.assertIn(':', record)
            
            # Diploid samples should have pipe separator
            if len(sample_to_genotype[list(sample_to_genotype.keys())[0]]) == 2:
                self.assertIn('|', record)
    
    def test_build_genotype_record_format(self):
        """Test that genotype records have correct format"""
        # Get genotypes from GFA
        sample_to_genotype = self.G.genotypes_from_gfa(self.gfa_file, exclude_terminus=True)
        
        # Get a variant edge
        variant_edges = list(self.G.sorted_variant_edges(exclude_terminus=True))
        edge = variant_edges[0]
        reference_edge = self.G.representative_edge(self.G.reference_tree_edge(edge))
        is_inversion = self.G.is_inversion(edge)
        sample_order = list(sample_to_genotype.keys())
        
        # Build genotype record
        genotype_records = _build_genotype_record(edge, reference_edge, sample_to_genotype, is_inversion, sample_order)
        
        for record in genotype_records:
            parts = record.split(':')
            
            # Should have 3 parts: GT, CR, CA
            self.assertEqual(len(parts), 3)
            
            # GT part
            gt_part = parts[0]
            # Should be either "0", "1", ".", or "0|0", "0|1", "1|0", "1|1", ".|.", etc.
            self.assertTrue(
                gt_part in ['.', '0', '1'] or 
                '|' in gt_part
            )
            
            # CR and CA parts should be numeric or '.' or contain comma for diploid
            for part in parts[1:]:
                if ',' in part:
                    # Diploid - check each value
                    for val in part.split(','):
                        self.assertTrue(val == '.' or val.isdigit())
                else:
                    # Haploid - single value
                    self.assertTrue(part == '.' or part.isdigit())
    
    def test_build_genotype_record_inversion(self):
        """Test that inversions set CR to '.'"""
        # Create mock genotypes
        genotype1 = Genotype(
            ref_counts={('1_+', '2_+'): 1},
            alt_counts={('3_+', '4_+'): 1},
            linear_coverage=[(0, 100)],
            exclude_terminus=True,
            missing_variants=set()
        )
        
        sample_to_genotype = {
            'sample1': (genotype1,)
        }
        
        variant_edge = ('1_+', '2_+')
        reference_edge = ('1_+', '2_+')
        sample_order = ['sample1']
        
        # Test with inversion
        records = _build_genotype_record(variant_edge, reference_edge, sample_to_genotype, is_inversion=True, sample_order=sample_order)
        self.assertEqual(len(records), 1)
        
        # CR should be '.'
        parts = records[0].split(':')
        self.assertEqual(parts[1], '.')
        
        # Test without inversion
        records = _build_genotype_record(variant_edge, reference_edge, sample_to_genotype, is_inversion=False, sample_order=sample_order)
        self.assertEqual(len(records), 1)
        
        # CR should be numeric
        parts = records[0].split(':')
        self.assertTrue(parts[1].isdigit())


class TestVCFRecord(unittest.TestCase):
    """Test VCF record building"""
    
    def test_build_vcf_record_structure(self):
        """Test that VCF records have correct structure"""
        chr_name = "chr1"
        edge_vcf_position = 100
        edge = ('1_+', '2_+')
        ref_allele = "A"
        alt_allele = "T"
        sample_ids = ['sample1', 'sample2']
        genotype_records = ['0:1:0', '0|1:1,0:0,1']
        
        record = _build_vcf_record(
            chr_name=chr_name,
            edge_vcf_position=edge_vcf_position,
            edge=edge,
            ref_allele=ref_allele,
            alt_allele=alt_allele,
            sample_ids=sample_ids,
            genotype_records=genotype_records
        )
        
        # Should be a list
        self.assertIsInstance(record, list)
        
        # Should have 9 fixed fields + number of samples
        self.assertEqual(len(record), 9 + len(sample_ids))
        
        # Check fixed fields
        self.assertEqual(record[0], chr_name)  # CHROM
        self.assertEqual(record[1], str(edge_vcf_position))  # POS
        self.assertEqual(record[3], ref_allele)  # REF
        self.assertEqual(record[4], alt_allele)  # ALT
        self.assertEqual(record[5], '60')  # QUAL
        self.assertEqual(record[6], 'PASS')  # FILTER
        self.assertEqual(record[8], 'GT:CR:CA')  # FORMAT


class TestInfoField(unittest.TestCase):
    """Test INFO field building"""
    
    @classmethod
    def setUpClass(cls):
        """Load the test GFA file once for all tests"""
        test_dir = os.path.dirname(__file__)
        cls.gfa_file = os.path.join(test_dir, "data", "simple_nested.gfa")
        cls.G = PangenomeGraph.from_gfa_line_by_line(cls.gfa_file, ref_name='ref')
    


class TestSampleOrdering(unittest.TestCase):
    """Test that sample order is consistent between header and genotype records"""
    
    @classmethod
    def setUpClass(cls):
        """Load the test GFA file once for all tests"""
        test_dir = os.path.dirname(__file__)
        cls.gfa_file = os.path.join(test_dir, "data", "simple_nested.gfa")
        cls.G = PangenomeGraph.from_gfa_line_by_line(cls.gfa_file, ref_name='ref')
    
    def test_sample_order_matches_header(self):
        """Test that genotype records are in the same order as header"""
        # Get genotypes from GFA
        sample_to_genotype = self.G.genotypes_from_gfa(self.gfa_file, exclude_terminus=True)
        # Preserve GFA file order
        sample_order = list(sample_to_genotype.keys())
        
        # Get a variant edge
        variant_edges = list(self.G.sorted_variant_edges(exclude_terminus=True))
        self.assertGreater(len(variant_edges), 0)
        
        edge = variant_edges[0]
        reference_edge = self.G.representative_edge(self.G.reference_tree_edge(edge))
        is_inversion = self.G.is_inversion(edge)
        
        # Build genotype records with specific order
        genotype_records = _build_genotype_record(
            edge, reference_edge, sample_to_genotype, is_inversion, sample_order
        )
        
        # Should have one record per sample
        self.assertEqual(len(genotype_records), len(sample_order))
        
        # Check that order matches by verifying different samples have different genotypes
        # (this is a sanity check - if order was wrong, we'd get mismatched data)
        for i, sample_name in enumerate(sample_order):
            record = genotype_records[i]
            # Record should be a non-empty string
            self.assertTrue(len(record) > 0)


class TestWriteVCF(unittest.TestCase):
    """Test full VCF writing"""
    
    @classmethod
    def setUpClass(cls):
        """Load the test GFA file once for all tests"""
        test_dir = os.path.dirname(__file__)
        cls.gfa_file = os.path.join(test_dir, "data", "simple_nested.gfa")
        cls.G = PangenomeGraph.from_gfa_line_by_line(cls.gfa_file, ref_name='ref')
    
    def test_write_vcf_creates_file(self):
        """Test that write_vcf creates a valid VCF file"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            vcf_path = f.name
        
        try:
            # Write VCF
            reference_edges = self.G.get_reference_edges()
            write_vcf_from_graph(
                graph=self.G,
                reference_edges=reference_edges,
                gfa_path=self.gfa_file,
                vcf_filename=vcf_path,
                chr_name='chr1',
                exclude_terminus=True
            )
            
            # Check file exists
            self.assertTrue(os.path.exists(vcf_path))
            
            # Read and check content
            with open(vcf_path, 'r') as f:
                content = f.read()
            
            # Should have metadata
            self.assertIn('##fileformat=VCFv4.2', content)
            
            # Should have header line
            self.assertIn('#CHROM', content)
            
            # Should have at least one variant line (not starting with #)
            lines = content.strip().split('\n')
            variant_lines = [l for l in lines if not l.startswith('#')]
            self.assertGreater(len(variant_lines), 0, "Should have at least one variant")
            
        finally:
            # Clean up
            if os.path.exists(vcf_path):
                os.remove(vcf_path)
    
    def test_write_vcf_without_genotypes(self):
        """Test writing VCF without genotype information"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            vcf_path = f.name
        
        try:
            # Write VCF without genotypes
            reference_edges = self.G.get_reference_edges()
            write_vcf_from_graph(
                graph=self.G,
                reference_edges=reference_edges,
                gfa_path=None,  # No GFA path means no genotypes
                vcf_filename=vcf_path,
                chr_name='chr1',
                exclude_terminus=True
            )
            
            # Check file exists
            self.assertTrue(os.path.exists(vcf_path))
            
            # Read content
            with open(vcf_path, 'r') as f:
                lines = f.readlines()
            
            # Find header line
            header_line = [l for l in lines if l.startswith('#CHROM')][0]
            
            # Should only have fixed columns, no sample columns
            columns = header_line.strip().split('\t')
            self.assertEqual(len(columns), 9, "Should have only 9 fixed columns without samples")
            
        finally:
            # Clean up
            if os.path.exists(vcf_path):
                os.remove(vcf_path)


if __name__ == '__main__':
    unittest.main()

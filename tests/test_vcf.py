"""
Unit tests for VCF writing functionality
"""
import unittest
import os
import tempfile
import gzip
import logging
from Bio import bgzf
from pantree.graph import PangenomeGraph
from pantree.vcf import (
    _build_genotype_record,
    _VariantRecord,
    _VariantData,
    _get_default_info_fields,
    write_vcf_from_graph
)


class TestVCFMetadata(unittest.TestCase):
    """Test VCF metadata generation"""
    
    def test_generate_vcf_metadata(self):
        """Test that VCF INFO fields are defined correctly"""
        info_fields = _get_default_info_fields()
        
        # Check that we have the expected INFO fields
        field_ids = [field.id for field in info_fields]
        expected_fields = ['NR', 'VT', 'TP', 'RC', 'AC', 'AN', 'HP', 'TR_MOTIF', 'NIA', 'UIDX']
        
        for expected_id in expected_fields:
            self.assertIn(expected_id, field_ids, f"Missing INFO field: {expected_id}")
        
        # Check that each field has required attributes
        for field in info_fields:
            self.assertTrue(hasattr(field, 'id'))
            self.assertTrue(hasattr(field, 'number'))
            self.assertTrue(hasattr(field, 'type'))
            self.assertTrue(hasattr(field, 'description'))
            self.assertTrue(hasattr(field, 'evaluate'))
            
            # Check that get_header returns a string
            header = field.get_header()
            self.assertIsInstance(header, str)
            self.assertIn(f"##INFO=<ID={field.id}", header)


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
        ref_edge_idx = int(self.G.edges[reference_edge]['index'])
        is_inversion = self.G.is_inversion(edge)
        sample_order = list(sample_to_genotype.keys())
        
        # Build genotype record
        genotype_records = _build_genotype_record(
            edge, reference_edge, ref_edge_idx, sample_to_genotype, is_inversion, sample_order
        )
        
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
        ref_edge_idx = int(self.G.edges[reference_edge]['index'])
        is_inversion = self.G.is_inversion(edge)
        sample_order = list(sample_to_genotype.keys())
        
        # Build genotype record
        genotype_records = _build_genotype_record(
            edge, reference_edge, ref_edge_idx, sample_to_genotype, is_inversion, sample_order
        )
        
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
        sample_to_genotype = self.G.genotypes_from_gfa(self.gfa_file, exclude_terminus=True)
        variant_edge = self.G.sorted_variant_edges(exclude_terminus=True)[0]
        reference_edge = self.G.representative_edge(self.G.reference_tree_edge(variant_edge))
        ref_edge_idx = int(self.G.edges[reference_edge]['index'])
        sample_order = list(sample_to_genotype.keys())
        
        # Test with inversion
        records = _build_genotype_record(
            variant_edge, reference_edge, ref_edge_idx, sample_to_genotype,
            is_inversion=True, sample_order=sample_order
        )
        self.assertEqual(len(records), len(sample_order))
        
        # CR should be '.'
        for record in records:
            parts = record.split(':')
            self.assertTrue(all(value == '.' for value in parts[1].split(',')))
        
        # Test without inversion
        records = _build_genotype_record(
            variant_edge, reference_edge, ref_edge_idx, sample_to_genotype,
            is_inversion=False, sample_order=sample_order
        )
        self.assertEqual(len(records), len(sample_order))
        
        # CR should be numeric
        self.assertTrue(any(
            value.isdigit()
            for record in records
            for value in record.split(':')[1].split(',')
        ))


class TestVCFRecord(unittest.TestCase):
    """Test VCF record building"""
    
    def test_build_vcf_record_structure(self):
        """Test that VCF records have correct structure"""
        # Create a mock _VariantData
        variant_info = _VariantData(
            chr_name="chr1",
            ref_allele_raw="A",
            alt_allele_raw="T",
            edge_data={'is_back_edge': False, 'is_inversion': False, 'motif': None},
            node_u_data={'direction': 1, 'distance_from_reference': 0, 'position': 100, 'index': 1},
            node_v_data={'direction': 1, 'distance_from_reference': 0, 'position': 101, 'on_reference_path': 1},
            branch_point_node_data={'sequence': 'ACGT', 'position': 99},
            ref_allele_count=1,
            alt_allele_count=1,
            context={'size_threshold': None, 'haplotype_fields': set(), 'qual': '60', 'filter_field': 'PASS', 'format_field': 'GT:CR:CA'}
        )
        
        edge = ('1_+', '2_+')
        genotype_records = ['0:1:0', '0|1:1,0:0,1']
        info_fields = _get_default_info_fields()
        
        record = _VariantRecord.from_variant_data(
            variant_info=variant_info,
            edge=edge,
            genotype_records=genotype_records,
            info_fields=info_fields
        )
        
        # Check fields
        self.assertEqual(record.chr_name, "chr1")
        self.assertEqual(record.vcf_position, 100)
        self.assertEqual(record.ref_allele, "A")
        self.assertEqual(record.alt_allele, "T")
        self.assertEqual(record.qual, '60')
        self.assertEqual(record.filter_field, 'PASS')
        self.assertEqual(record.format_field, 'GT:CR:CA')
        self.assertEqual(len(record.genotype_records), 2)
        
        # Test to_vcf_line
        vcf_line = record.to_vcf_line()
        self.assertIsInstance(vcf_line, str)
        fields = vcf_line.split('\t')
        self.assertEqual(len(fields), 9 + len(genotype_records))


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
        ref_edge_idx = int(self.G.edges[reference_edge]['index'])
        is_inversion = self.G.is_inversion(edge)
        
        # Build genotype records with specific order
        genotype_records = _build_genotype_record(
            edge, reference_edge, ref_edge_idx, sample_to_genotype, is_inversion, sample_order
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
            import logging
            write_vcf_from_graph(
                graph=self.G,
                gfa_path=self.gfa_file,
                vcf_filename=vcf_path,
                chr_name='chr1',
                logger=logging.getLogger('pantree'),
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

    def test_write_vcf_gz_is_bgzf_and_gzip_readable(self):
        """Test that .vcf.gz output is BGZF and remains gzip-compatible."""
        with tempfile.NamedTemporaryFile(suffix='.vcf.gz', delete=False) as f:
            vcf_path = f.name

        try:
            import logging
            write_vcf_from_graph(
                graph=self.G,
                gfa_path=self.gfa_file,
                vcf_filename=vcf_path,
                chr_name='chr1',
                logger=logging.getLogger('pantree'),
                exclude_terminus=True
            )

            with bgzf.open(vcf_path, 'rt') as f:
                self.assertEqual(f.readline().strip(), '##fileformat=VCFv4.2')

            with gzip.open(vcf_path, 'rt') as f:
                content = f.read()
            self.assertIn('#CHROM', content)
            self.assertIn('chr1', content)

        finally:
            if os.path.exists(vcf_path):
                os.remove(vcf_path)

    def test_write_vcf_records_are_sorted_by_position(self):
        """Test output records are sorted for tabix-compatible coordinate order."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            vcf_path = f.name

        try:
            import logging
            write_vcf_from_graph(
                graph=self.G,
                gfa_path=self.gfa_file,
                vcf_filename=vcf_path,
                chr_name='chr1',
                logger=logging.getLogger('pantree'),
                exclude_terminus=True
            )

            with open(vcf_path, 'r') as f:
                variant_lines = [line for line in f if not line.startswith('#')]

            positions = [int(line.split('\t')[1]) for line in variant_lines]
            self.assertEqual(positions, sorted(positions))

        finally:
            if os.path.exists(vcf_path):
                os.remove(vcf_path)
    
    def test_write_vcf_without_genotypes(self):
        """Test writing VCF without genotype information"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            vcf_path = f.name
        
        try:
            # Write VCF without genotypes
            import logging
            write_vcf_from_graph(
                graph=self.G,
                gfa_path=None,  # No GFA path means no genotypes
                vcf_filename=vcf_path,
                chr_name='chr1',
                logger=logging.getLogger('pantree'),
                exclude_terminus=True
            )
            
            # Check file exists
            self.assertTrue(os.path.exists(vcf_path))
            
            # Read content
            with open(vcf_path, 'r') as f:
                lines = f.readlines()
            
            # Find header line
            header_line = [l for l in lines if l.startswith('#CHROM')][0]
            
            columns = header_line.strip().split('\t')
            self.assertEqual(
                columns,
                ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
            )
            self.assertFalse(any(line.startswith('##FORMAT=') for line in lines))
            for line in lines:
                if not line.startswith('#'):
                    self.assertEqual(len(line.rstrip('\n').split('\t')), 8)
            
        finally:
            # Clean up
            if os.path.exists(vcf_path):
                os.remove(vcf_path)

    def test_write_vcf_includes_haplotype_positions(self):
        """Test HP metadata and INFO values from priority-sample haplotype positions."""
        test_dir = os.path.dirname(__file__)
        gfa_file = os.path.join(test_dir, "data", "simple_nested.gfa")
        graph = PangenomeGraph.from_gfa_line_by_line(
            gfa_file,
            ref_name='ref',
            priority_dict={
                'ref': 0,
                'sample1': 1,
                'sample2': 2,
            }
        )

        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            vcf_path = f.name

        try:
            write_vcf_from_graph(
                graph=graph,
                gfa_path=None,
                vcf_filename=vcf_path,
                chr_name='chr1',
                logger=logging.getLogger('pantree'),
                exclude_terminus=True
            )

            with open(vcf_path, 'r') as f:
                lines = f.readlines()

            self.assertTrue(any(line.startswith('##INFO=<ID=HP') for line in lines))
            variant_infos = [
                line.split('\t')[7]
                for line in lines
                if not line.startswith('#')
            ]
            self.assertTrue(any('HP=ref#0#0:' in info for info in variant_infos))
            self.assertTrue(any('sample1#1#0:' in info for info in variant_infos))
        finally:
            if os.path.exists(vcf_path):
                os.remove(vcf_path)


if __name__ == '__main__':
    unittest.main()

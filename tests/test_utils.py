import unittest
from pantree.utils import nearly_identical_alleles

class TestUtils(unittest.TestCase):
    def test_nearly_identical_alleles(self):
        # Test identical sequences
        self.assertTrue(nearly_identical_alleles("ACGT", "ACGT"))
        
        # Test single SNP
        self.assertTrue(nearly_identical_alleles("ACGT", "ACTT"))
        self.assertTrue(nearly_identical_alleles("ACGT", "ACAT"))
        
        # Test single deletion
        self.assertTrue(nearly_identical_alleles("ACGT", "ACT"))
        self.assertTrue(nearly_identical_alleles("ACGT", "AGT"))
        
        # Test single insertion
        self.assertTrue(nearly_identical_alleles("ACT", "ACGT"))
        self.assertTrue(nearly_identical_alleles("AGT", "AGCT"))
        
        # Test sequences that should not be nearly identical
        # Multiple SNPs
        self.assertFalse(nearly_identical_alleles("ACGT", "AATT"))
        
        # Length difference > 1
        self.assertFalse(nearly_identical_alleles("ACGT", "AC"))
        self.assertFalse(nearly_identical_alleles("AC", "ACGT"))
        
        # Empty sequences
        self.assertTrue(nearly_identical_alleles("A", ""))
        self.assertTrue(nearly_identical_alleles("", "A"))
        self.assertTrue(nearly_identical_alleles("", ""))

if __name__ == '__main__':
    unittest.main()

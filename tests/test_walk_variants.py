"""
Test suite for pantree.walk_variants module.

This module tests the Allele class and related functions, which provide data structures
for working with genomic sequences at specific positions on a reference genome.

## Overview of walk_variants.py

The module implements:

1. **Allele class**: Represents a genomic sequence at a specific position
   - `sequence`: A list of bases (or any sequence type)
   - `position`: The starting genomic coordinate
   - `interval`: Property returning (start, end) genomic coordinates

2. **Key Methods**:
   - `__getitem__`: Slice using genomic coordinates (not list indices)
   - `is_disjoint`: Check if two alleles don't overlap
   - `matches`: Check if overlapping regions have the same sequence
   - `prepend`: Extend allele leftward by prepending another overlapping allele
   - `replace`: Apply a variant by replacing a ref allele with an alt allele

3. **overlap() function**: Calculate the overlapping genomic interval between two alleles

4. **get_walk_variants() function**: Combine overlapping variants from a VCF into
   consolidated ref/alt sequences

## Use Case

This is designed for variant graph traversal, where you need to:
- Track sequences at genomic positions
- Extend sequences as you walk through a graph
- Apply variants (substitutions, insertions, deletions) to sequences
- Combine multiple overlapping variants into single consolidated sequences
"""

import pytest
import polars as pl
from pantree.graph import PangenomeGraph
from pantree.walk_variants import Allele, overlap, get_walk_variants, process_haplotype_variants


class TestOverlapFunction:
    """Tests for the overlap() function."""

    def test_overlap_complete_overlap(self):
        """Test when one allele completely contains another."""
        a = Allele(sequence=['A', 'T', 'G', 'C'], position=10)
        b = Allele(sequence=['T', 'G'], position=11)
        start, end = overlap(a, b)
        assert start == 11
        assert end == 13

    def test_overlap_partial_overlap(self):
        """Test when two alleles partially overlap."""
        a = Allele(sequence=['A', 'T', 'G', 'C'], position=10)
        b = Allele(sequence=['G', 'C', 'A'], position=12)
        start, end = overlap(a, b)
        assert start == 12
        assert end == 14

    def test_overlap_no_overlap(self):
        """Test when two alleles don't overlap."""
        a = Allele(sequence=['A', 'T'], position=10)
        b = Allele(sequence=['G', 'C'], position=15)
        start, end = overlap(a, b)
        assert start == 15
        assert end == 12  # end < start indicates no overlap

    def test_overlap_adjacent(self):
        """Test when two alleles are adjacent but don't overlap."""
        a = Allele(sequence=['A', 'T'], position=10)
        b = Allele(sequence=['G', 'C'], position=12)
        start, end = overlap(a, b)
        assert start == 12
        assert end == 12  # Adjacent, no overlap

    def test_overlap_exact_same(self):
        """Test when two alleles have the same position and length."""
        a = Allele(sequence=['A', 'T', 'G'], position=10)
        b = Allele(sequence=['A', 'T', 'G'], position=10)
        start, end = overlap(a, b)
        assert start == 10
        assert end == 13


class TestAlleleBasics:
    """Tests for basic Allele properties and methods."""

    def test_allele_creation(self):
        """Test basic allele creation."""
        allele = Allele(sequence=['A', 'T', 'G', 'C'], position=100)
        assert allele.sequence == ['A', 'T', 'G', 'C']
        assert allele.position == 100

    def test_allele_length(self):
        """Test __len__ method."""
        allele = Allele(sequence=['A', 'T', 'G', 'C'], position=100)
        assert len(allele) == 4

    def test_allele_interval(self):
        """Test interval property."""
        allele = Allele(sequence=['A', 'T', 'G', 'C'], position=100)
        assert allele.interval == (100, 104)

    def test_allele_interval_empty(self):
        """Test interval property with empty sequence."""
        allele = Allele(sequence=[], position=50)
        assert allele.interval == (50, 50)


class TestAlleleGetItem:
    """Tests for the __getitem__ slicing method."""

    def test_getitem_full_range(self):
        """Test slicing the entire allele."""
        allele = Allele(sequence=['A', 'T', 'G', 'C'], position=10)
        result = allele[10:14]
        assert result == ['A', 'T', 'G', 'C']

    def test_getitem_partial_range(self):
        """Test slicing a portion of the allele."""
        allele = Allele(sequence=['A', 'T', 'G', 'C'], position=10)
        result = allele[11:13]
        assert result == ['T', 'G']

    def test_getitem_single_base(self):
        """Test slicing a single position."""
        allele = Allele(sequence=['A', 'T', 'G', 'C'], position=10)
        result = allele[10:11]
        assert result == ['A']

    def test_getitem_end_portion(self):
        """Test slicing the end portion."""
        allele = Allele(sequence=['A', 'T', 'G', 'C'], position=10)
        result = allele[12:14]
        assert result == ['G', 'C']

    def test_getitem_start_portion(self):
        """Test slicing the start portion."""
        allele = Allele(sequence=['A', 'T', 'G', 'C'], position=10)
        result = allele[10:12]
        assert result == ['A', 'T']


class TestAlleleIsDisjoint:
    """Tests for the is_disjoint method."""

    def test_is_disjoint_no_overlap(self):
        """Test alleles that don't overlap."""
        a = Allele(sequence=['A', 'T'], position=10)
        b = Allele(sequence=['G', 'C'], position=15)
        assert a.is_disjoint(b)
        assert b.is_disjoint(a)

    def test_is_disjoint_adjacent(self):
        """Test alleles that are adjacent (touching) are NOT disjoint."""
        a = Allele(sequence=['A', 'T'], position=10)
        b = Allele(sequence=['G', 'C'], position=12)
        assert not a.is_disjoint(b)
        assert not b.is_disjoint(a)

    def test_is_disjoint_with_overlap(self):
        """Test alleles that overlap."""
        a = Allele(sequence=['A', 'T', 'G', 'C'], position=10)
        b = Allele(sequence=['G', 'C', 'A'], position=12)
        assert not a.is_disjoint(b)
        assert not b.is_disjoint(a)

    def test_is_disjoint_complete_overlap(self):
        """Test when one allele is completely inside another."""
        a = Allele(sequence=['A', 'T', 'G', 'C'], position=10)
        b = Allele(sequence=['T', 'G'], position=11)
        assert not a.is_disjoint(b)
        assert not b.is_disjoint(a)


class TestAlleleMatches:
    """Tests for the matches method."""

    def test_matches_identical_alleles(self):
        """Test matching with identical alleles."""
        a = Allele(sequence=['A', 'T', 'G', 'C'], position=10)
        b = Allele(sequence=['A', 'T', 'G', 'C'], position=10)
        assert a.matches(b)
        assert b.matches(a)

    def test_matches_overlapping_same_sequence(self):
        """Test matching alleles with same sequence in overlap."""
        a = Allele(sequence=['A', 'T', 'G', 'C'], position=10)
        b = Allele(sequence=['G', 'C', 'A', 'T'], position=12)
        assert a.matches(b)
        assert b.matches(a)

    def test_matches_overlapping_different_sequence(self):
        """Test non-matching alleles with different sequence in overlap."""
        a = Allele(sequence=['A', 'T', 'G', 'C'], position=10)
        b = Allele(sequence=['A', 'A', 'T', 'T'], position=12)
        assert not a.matches(b)
        assert not b.matches(a)

    def test_matches_partial_overlap_matching(self):
        """Test matching with partial overlap."""
        a = Allele(sequence=['A', 'T', 'G', 'C', 'G'], position=10)
        b = Allele(sequence=['T', 'G'], position=11)
        assert a.matches(b)
        assert b.matches(a)

    def test_matches_single_base_overlap(self):
        """Test matching with single base overlap."""
        a = Allele(sequence=['A', 'T'], position=10)
        b = Allele(sequence=['T', 'G'], position=11)
        assert a.matches(b)
        assert b.matches(a)


class TestAllelePrepend:
    """Tests for the prepend method."""

    def test_prepend_extends_sequence(self):
        """Test prepending an allele that extends to the left."""
        a = Allele(sequence=['G', 'C', 'A'], position=12)
        b = Allele(sequence=['A', 'T', 'G', 'C'], position=10)
        a.prepend(b)
        assert a.position == 10
        assert a.sequence == ['A', 'T', 'G', 'C', 'A']

    def test_prepend_with_disjoint_raises_error(self):
        """Test prepending disjoint alleles raises ValueError."""
        a = Allele(sequence=['G', 'C'], position=15)
        b = Allele(sequence=['A', 'T'], position=10)
        with pytest.raises(ValueError, match="Alleles are disjoint"):
            a.prepend(b)

    def test_prepend_with_mismatch_no_longer_requires_match(self):
        """Test prepending no longer requires matching sequences."""
        a = Allele(sequence=['G', 'C', 'A'], position=12)
        b = Allele(sequence=['A', 'T', 'A', 'A'], position=10)
        a.prepend(b)
        # Prepends non-overlapping part from b
        assert a.position == 10
        assert a.sequence == ['A', 'T', 'G', 'C', 'A']

    def test_prepend_no_op_when_other_after(self):
        """Test prepending when other starts at or after self (no-op)."""
        a = Allele(sequence=['A', 'T', 'G', 'C'], position=10)
        original_position = a.position
        original_sequence = a.sequence.copy()
        b = Allele(sequence=['G', 'C'], position=12)
        a.prepend(b)
        assert a.position == original_position
        assert a.sequence == original_sequence

    def test_prepend_exact_overlap(self):
        """Test prepending when alleles exactly overlap (no extension)."""
        a = Allele(sequence=['A', 'T', 'G'], position=10)
        b = Allele(sequence=['A', 'T', 'G'], position=10)
        a.prepend(b)
        assert a.position == 10
        assert a.sequence == ['A', 'T', 'G']

    def test_prepend_partial_overlap(self):
        """Test prepending with partial overlap."""
        a = Allele(sequence=['T', 'G', 'C', 'A'], position=11)
        b = Allele(sequence=['A', 'T', 'G'], position=10)
        a.prepend(b)
        assert a.position == 10
        assert a.sequence == ['A', 'T', 'G', 'C', 'A']

    def test_prepend_modifies_in_place(self):
        """Test that prepend modifies the allele in place."""
        a = Allele(sequence=['G', 'C'], position=12)
        b = Allele(sequence=['A', 'T', 'G'], position=10)
        a.prepend(b)
        assert a.position == 10
        assert a.sequence == ['A', 'T', 'G', 'C']


class TestAlleleReplace:
    """Tests for the replace method."""

    def test_replace_simple_substitution(self):
        """Test simple base substitution."""
        allele = Allele(sequence=['A', 'T', 'G', 'C'], position=10)
        ref = Allele(sequence=['T'], position=11)
        alt = Allele(sequence=['C'], position=11)
        allele.replace(ref, alt)
        assert allele.sequence == ['A', 'C', 'G', 'C']
        assert allele.position == 10

    def test_replace_deletion(self):
        """Test deletion (ref longer than alt)."""
        allele = Allele(sequence=['A', 'T', 'G', 'C'], position=10)
        ref = Allele(sequence=['T', 'G'], position=11)
        alt = Allele(sequence=[], position=11)
        allele.replace(ref, alt)
        assert allele.sequence == ['A', 'C']
        assert allele.position == 10

    def test_replace_insertion(self):
        """Test insertion (alt longer than ref)."""
        allele = Allele(sequence=['A', 'T', 'G', 'C'], position=10)
        ref = Allele(sequence=['T'], position=11)
        alt = Allele(sequence=['T', 'T', 'T'], position=11)
        allele.replace(ref, alt)
        assert allele.sequence == ['A', 'T', 'T', 'T', 'G', 'C']
        assert allele.position == 10

    def test_replace_entire_sequence(self):
        """Test replacing the entire sequence."""
        allele = Allele(sequence=['A', 'T', 'G', 'C'], position=10)
        ref = Allele(sequence=['A', 'T', 'G', 'C'], position=10)
        alt = Allele(sequence=['G', 'G'], position=10)
        allele.replace(ref, alt)
        assert allele.sequence == ['G', 'G']
        assert allele.position == 10

    def test_replace_at_start(self):
        """Test replacing at the start of the sequence."""
        allele = Allele(sequence=['A', 'T', 'G', 'C'], position=10)
        ref = Allele(sequence=['A'], position=10)
        alt = Allele(sequence=['G'], position=10)
        allele.replace(ref, alt)
        assert allele.sequence == ['G', 'T', 'G', 'C']
        assert allele.position == 10

    def test_replace_at_end(self):
        """Test replacing at the end of the sequence."""
        allele = Allele(sequence=['A', 'T', 'G', 'C'], position=10)
        ref = Allele(sequence=['C'], position=13)
        alt = Allele(sequence=['A'], position=13)
        allele.replace(ref, alt)
        assert allele.sequence == ['A', 'T', 'G', 'A']
        assert allele.position == 10

    def test_replace_with_prepend(self):
        """Test replace when ref extends before allele start."""
        allele = Allele(sequence=['T', 'G', 'C'], position=11)
        ref = Allele(sequence=['A', 'T'], position=10)
        alt = Allele(sequence=['G', 'G'], position=10)
        allele.replace(ref, alt)
        assert allele.position == 10
        assert allele.sequence == ['G', 'G', 'G', 'C']

    def test_replace_multiple_bases(self):
        """Test replacing multiple bases in the middle."""
        allele = Allele(sequence=['A', 'T', 'G', 'C', 'A'], position=10)
        ref = Allele(sequence=['T', 'G', 'C'], position=11)
        alt = Allele(sequence=['C', 'C'], position=11)
        allele.replace(ref, alt)
        assert allele.sequence == ['A', 'C', 'C', 'A']
        assert allele.position == 10

    def test_replace_with_position_zero(self):
        """Test replace works correctly with position 0."""
        allele = Allele(sequence=['A', 'T', 'G', 'C'], position=0)
        ref = Allele(sequence=['T'], position=1)
        alt = Allele(sequence=['C'], position=1)
        allele.replace(ref, alt)
        assert allele.sequence == ['A', 'C', 'G', 'C']
        assert allele.position == 0


class TestAllelePrependLeftOnly:
    """Tests confirming prepend only extends to the left."""

    def test_prepend_no_op_when_other_at_same_or_after(self):
        """Test that prepend is no-op when other doesn't extend left of allele."""
        allele = Allele(sequence=['A', 'T', 'G'], position=10)
        ref = Allele(sequence=['A', 'T'], position=10)
        allele.prepend(ref)
        # Should be no-op since ref doesn't extend left
        assert allele.position == 10
        assert allele.sequence == ['A', 'T', 'G']

    def test_prepend_left_extension_only(self):
        """Test prepend only extends to the left."""
        allele = Allele(sequence=['T', 'G'], position=11)
        ref = Allele(sequence=['A', 'T'], position=10)
        allele.prepend(ref)
        assert allele.position == 10
        assert allele.sequence == ['A', 'T', 'G']


class TestAlleleEdgeCases:
    """Tests for edge cases and special scenarios."""

    def test_empty_sequence_allele(self):
        """Test alleles with empty sequences."""
        allele = Allele(sequence=[], position=10)
        assert len(allele) == 0
        assert allele.interval == (10, 10)

    def test_single_base_allele(self):
        """Test alleles with single base."""
        allele = Allele(sequence=['A'], position=10)
        assert len(allele) == 1
        assert allele.interval == (10, 11)
        assert allele[10:11] == ['A']

    def test_allele_with_string_sequences(self):
        """Test that alleles work with string sequences (not just lists)."""
        allele = Allele(sequence="ATGC", position=10)
        assert len(allele) == 4
        assert allele[10:12] == "AT"


class TestGetWalkVariants:
    """Tests for the get_walk_variants function."""

    def test_single_variant(self):
        """Test with a single variant."""
        df = pl.LazyFrame({
            "UIDX": [0],
            "TP": [10],
            "REF": ["A"],
            "ALT": ["T"]
        })
        results = get_walk_variants(df)
        assert len(results) == 1
        assert results[0].ref.sequence == ['A']
        assert results[0].ref.position == 10
        assert results[0].alt.sequence == ['T']
        assert results[0].alt.position == 10
        assert results[0].status == "OK"

    def test_two_disjoint_variants(self):
        """Test with two non-overlapping variants."""
        df = pl.LazyFrame({
            "UIDX": [0, 1],
            "TP": [10, 20],
            "REF": ["A", "C"],
            "ALT": ["T", "G"]
        })
        results = get_walk_variants(df)
        assert len(results) == 2
        assert all(r.status == "OK" for r in results)

    def test_empty_vcf(self):
        """Test with empty VCF."""
        df = pl.LazyFrame({
            "UIDX": [],
            "TP": [],
            "REF": [],
            "ALT": []
        }, schema={"UIDX": pl.Int64, "TP": pl.Int64, "REF": pl.Utf8, "ALT": pl.Utf8})
        results = get_walk_variants(df)
        assert len(results) == 0

    def test_example1(self):
        """Test example 1: A->AAAA at 10 (UIDX=4), A->T at 12 (UIDX=2).

        Expected output: 'A' 'AATA' 10
        """
        df = pl.LazyFrame({
            "UIDX": [4, 2],
            "TP": [10, 12],
            "REF": ["A", "A"],
            "ALT": ["AAAA", "T"]
        })
        results = get_walk_variants(df)

        assert len(results) == 1
        assert results[0].ref.position == 10
        assert results[0].ref.sequence == ['A']
        assert results[0].alt.position == 10
        assert results[0].alt.sequence == ['A', 'A', 'T', 'A']
        assert results[0].status == "OK"

    def test_example2(self):
        """Test example 2: A->AAAA at 10 (UIDX=4), T->A at 12 (UIDX=2).

        Expected output: ERR status (ref mismatch) - no longer raises exception
        """
        df = pl.LazyFrame({
            "UIDX": [4, 2],
            "TP": [10, 12],
            "REF": ["A", "T"],
            "ALT": ["AAAA", "A"]
        })
        results = get_walk_variants(df)
        # Should have 2 results: the first variant as OK, and the mismatched one as ERR
        assert len(results) == 2
        assert results[0].status == "OK"
        assert results[1].status == "ERR"

    def test_example3(self):
        """Test example 3: T->A at 12 (UIDX=2), AAAAA->TTT at 10 (UIDX=4).

        Expected output: 'AAAAA' 'TTA' 10
        """
        df = pl.LazyFrame({
            "UIDX": [2, 4],
            "TP": [12, 10],
            "REF": ["T", "AAAAA"],
            "ALT": ["A", "TTT"]
        })
        results = get_walk_variants(df)

        assert len(results) == 1
        assert results[0].ref.position == 10
        assert results[0].ref.sequence == ['A', 'A', 'A', 'A', 'A']
        assert results[0].alt.position == 10
        assert results[0].alt.sequence == ['T', 'T', 'A']
        assert results[0].status == "OK"

    def test_example4(self):
        """Test example 4: TT->AA at 12 (UIDX=2), CA->'' at 11 (UIDX=1).

        Expected output: 'CTT' 'A' 11
        """
        df = pl.LazyFrame({
            "UIDX": [2, 1],
            "TP": [12, 11],
            "REF": ["TT", "CA"],
            "ALT": ["AA", ""]
        })
        results = get_walk_variants(df)

        assert len(results) == 1
        assert results[0].ref.position == 11
        assert results[0].ref.sequence == ['C', 'T', 'T']
        assert results[0].alt.position == 11
        assert results[0].alt.sequence == ['A']
        assert results[0].status == "OK"

    def test_example5(self):
        """Test example 5: T->A at 12 (UIDX=2), C->G at 11 (UIDX=1).

        Expected output: 'TC' 'AG' 11
        """
        df = pl.LazyFrame({
            "UIDX": [2, 1],
            "TP": [12, 11],
            "REF": ["T", "C"],
            "ALT": ["A", "G"]
        })
        results = get_walk_variants(df)

        assert len(results) == 1
        assert results[0].ref.position == 11
        assert results[0].ref.sequence == ['C', 'T']
        assert results[0].alt.position == 11
        assert results[0].alt.sequence == ['G', 'A']
        assert results[0].status == "OK"

    def test_inv_variant_separate(self):
        """Test that INV variants are returned separately with INV status."""
        df = pl.LazyFrame({
            "UIDX": [0, 1],
            "TP": [10, 15],
            "REF": ["A", "C"],
            "ALT": ["T", "G"],
            "VT": ["SNP", "INV"]
        })
        results = get_walk_variants(df)
        assert len(results) == 2
        # INV variant should have INV status
        inv_results = [r for r in results if r.status == "INV"]
        assert len(inv_results) == 1

    def test_dup_variant_separate(self):
        """Test that DUP variants are returned separately with DUP status."""
        df = pl.LazyFrame({
            "UIDX": [0, 1],
            "TP": [10, 15],
            "REF": ["A", "CC"],
            "ALT": ["T", "CCCC"],
            "VT": ["SNP", "DUP"]
        })
        results = get_walk_variants(df)
        assert len(results) == 2
        # DUP variant should have DUP status
        dup_results = [r for r in results if r.status == "DUP"]
        assert len(dup_results) == 1


def test_process_haplotype_variants_accepts_bgzf_input(tmp_path):
    """Test consolidate input can be a Pantree-generated BGZF VCF."""
    gfa_path = "tests/data/simple_nested.gfa"
    uncompressed_vcf = tmp_path / "input.vcf"
    compressed_vcf = tmp_path / "input.vcf.gz"
    uncompressed_output = tmp_path / "uncompressed.out.vcf"
    compressed_output = tmp_path / "compressed.out.vcf"

    graph = PangenomeGraph.from_gfa_line_by_line(gfa_path, ref_name='ref')
    graph.write_vcf(gfa_path, str(uncompressed_vcf), chr_name='chr0')
    graph.write_vcf(gfa_path, str(compressed_vcf), chr_name='chr0')

    process_haplotype_variants(str(uncompressed_vcf), 'sample2', 0, str(uncompressed_output))
    process_haplotype_variants(str(compressed_vcf), 'sample2', 0, str(compressed_output))

    assert compressed_output.read_text() == uncompressed_output.read_text()

"""
Tests for overlap detection and Tm calculation.
"""

import pytest

from src.core.overlap import (
    find_overlap,
    get_Tm,
    design_overlap,
    gc_content,
    check_homopolymer,
)


class TestFindOverlap:
    """Tests for find_overlap function."""

    def test_exact_overlap(self):
        """Test finding exact overlap between sequences."""
        seq1 = "ATGCGATCGATCGATCGATCGTAGCTAGCT"
        seq2 = "GCTAGCTATGCATGCATGC"

        length, overlap = find_overlap(seq1, seq2, min_overlap=8)

        assert length == 8
        assert overlap == "GCTAGCTA"

    def test_no_overlap(self):
        """Test when sequences don't overlap."""
        seq1 = "AAAAAAAAAA"
        seq2 = "TTTTTTTTTT"

        length, overlap = find_overlap(seq1, seq2, min_overlap=5)

        assert length == 0
        assert overlap == ""

    def test_full_overlap(self):
        """Test when one sequence is entirely an overlap."""
        seq1 = "ATGCGATCGATCG"
        seq2 = "ATCGATCGATCGAAAAA"

        length, overlap = find_overlap(seq1, seq2, min_overlap=10)

        assert length >= 10
        assert seq1.endswith(overlap)
        assert seq2.startswith(overlap)

    def test_case_insensitive(self):
        """Test that overlap search is case insensitive."""
        seq1 = "atgcgatcg"
        seq2 = "GATCGAAAA"

        length, overlap = find_overlap(seq1, seq2, min_overlap=5)

        assert length == 5
        assert overlap == "GATCG"


class TestGetTm:
    """Tests for Tm calculation."""

    def test_basic_tm(self):
        """Test basic Tm calculation."""
        seq = "ATGCGATCGATCGATCGATCGTAG"  # 24 bp

        tm = get_Tm(seq)

        # Should be in reasonable range for this sequence
        assert 40 < tm < 70

    def test_high_gc_higher_tm(self):
        """Test that high GC content gives higher Tm."""
        low_gc = "AAATTTAAATTTAAATTT"
        high_gc = "GGGCCCGGGCCCGGGCCC"

        tm_low = get_Tm(low_gc)
        tm_high = get_Tm(high_gc)

        assert tm_high > tm_low

    def test_gc_method(self):
        """Test GC method gives different result."""
        seq = "ATGCGATCGATCGATCGATCGTAG"

        tm_nn = get_Tm(seq, method="nearest_neighbor")
        tm_gc = get_Tm(seq, method="gc")

        # Methods should give different (but similar) results
        assert abs(tm_nn - tm_gc) < 20  # Within 20°C


class TestGCContent:
    """Tests for GC content calculation."""

    def test_all_gc(self):
        """Test 100% GC content."""
        seq = "GGGGCCCC"

        assert gc_content(seq) == 1.0

    def test_all_at(self):
        """Test 0% GC content."""
        seq = "AAAATTTT"

        assert gc_content(seq) == 0.0

    def test_50_percent(self):
        """Test 50% GC content."""
        seq = "ATGC"

        assert gc_content(seq) == 0.5

    def test_empty_sequence(self):
        """Test empty sequence returns 0."""
        assert gc_content("") == 0.0


class TestCheckHomopolymer:
    """Tests for homopolymer detection."""

    def test_find_long_run(self):
        """Test finding a long homopolymer run."""
        seq = "ATGCAAAAAAAAAATGC"

        violations = check_homopolymer(seq, max_run=6)

        assert len(violations) == 1
        assert violations[0][1] == "A"  # Nucleotide
        assert violations[0][2] == 10  # Length

    def test_no_violations(self):
        """Test sequence with no long runs."""
        seq = "ATGCATGCATGC"

        violations = check_homopolymer(seq, max_run=6)

        assert len(violations) == 0

    def test_multiple_violations(self):
        """Test finding multiple homopolymer runs."""
        seq = "AAAAAAAAATTTTTTTTT"

        violations = check_homopolymer(seq, max_run=6)

        assert len(violations) == 2


class TestDesignOverlap:
    """Tests for overlap design."""

    def test_design_with_matching_sequences(self):
        """Test designing overlap with matching ends."""
        seq1_end = "NNNNNNATGCGATCGATCGATCGATCGATCGATCG"
        seq2_start = "ATGCGATCGATCGATCGATCGATCGATCGNNNNN"

        length, overlap = design_overlap(
            target_tm=52.0,
            seq1_end=seq1_end,
            seq2_start=seq2_start
        )

        assert 15 <= length <= 70
        assert overlap in seq1_end
        assert overlap in seq2_start

    def test_design_no_match_raises(self):
        """Test that non-matching sequences raise error."""
        seq1_end = "AAAAAAAAAAAAAAAA"
        seq2_start = "TTTTTTTTTTTTTTTT"

        with pytest.raises(ValueError, match="Could not design"):
            design_overlap(
                target_tm=52.0,
                seq1_end=seq1_end,
                seq2_start=seq2_start
            )

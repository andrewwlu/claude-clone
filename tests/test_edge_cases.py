"""
Edge case and regression tests.

These tests cover bugs that were fixed in past versions.
"""

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from src.core.backbone import cut_backbone_with_enzymes, find_cut_sites
from src.core.assembly import rotate_assembly
from src.core.overlap import find_overlap, check_homopolymer
from src.core.validation import validate_sequence


class TestBackboneCutEdgeCases:
    """Edge cases for backbone cutting."""

    def test_asymmetric_overhang_enzyme(self):
        """Test cutting with enzymes that have asymmetric overhangs.

        Regression test for v0.2.1 bug where asymmetric overhangs
        caused incorrect cut positions.
        """
        # Sequence with EcoRI site (G/AATTC)
        seq = "ATGC" * 10 + "GAATTC" + "ATGC" * 10

        sites = find_cut_sites(seq, "EcoRI")

        assert len(sites) == 1
        # Cut should be after G (position 40 + 1 = 41)
        assert sites[0] == 41

    def test_multiple_sites_raises(self):
        """Test that multiple cut sites raises appropriate error."""
        # Sequence with two XhoI sites
        seq = "CTCGAG" + "ATGC" * 20 + "CTCGAG"

        backbone = SeqRecord(Seq(seq), id="test")

        with pytest.raises(ValueError, match="multiple times"):
            cut_backbone_with_enzymes(backbone, "XhoI", "HindIII")

    def test_no_sites_raises(self):
        """Test that missing cut sites raises appropriate error."""
        seq = "ATGCATGCATGC" * 10

        backbone = SeqRecord(Seq(seq), id="test")

        with pytest.raises(ValueError, match="does not cut"):
            cut_backbone_with_enzymes(backbone, "XhoI", "HindIII")

    def test_sites_too_close(self):
        """Test that cut sites too close together raises error."""
        # XhoI and HindIII only 10 bp apart
        seq = "ATGC" * 10 + "CTCGAG" + "ATGC" * 2 + "AAGCTT" + "ATGC" * 10

        backbone = SeqRecord(Seq(seq), id="test")

        with pytest.raises(ValueError, match="too close"):
            cut_backbone_with_enzymes(backbone, "XhoI", "HindIII", min_overlap=20)


class TestRotationEdgeCases:
    """Edge cases for sequence rotation."""

    def test_feature_across_origin_after_rotation(self):
        """Test feature that spans origin after rotation.

        Regression test for issue #12 fixed in v0.3.2.
        """
        seq = "NNNNNSTARTSEQ" + "ATGCGATCG"
        record = SeqRecord(Seq(seq), id="test")

        # Add feature near end that will wrap after rotation
        record.features.append(SeqFeature(
            FeatureLocation(18, 22),  # Near end
            type="misc_feature"
        ))

        rotated = rotate_assembly(record, "STARTSEQ")

        # After rotation, original position 18 should be at different position
        assert len(rotated.features) > 0

    def test_rotation_preserves_feature_count(self):
        """Test that rotation doesn't lose features."""
        seq = "NNNSTARTNNNN"
        record = SeqRecord(Seq(seq), id="test")

        # Add multiple features
        for i in range(5):
            record.features.append(SeqFeature(
                FeatureLocation(i, i + 2),
                type="misc_feature"
            ))

        rotated = rotate_assembly(record, "START")

        assert len(rotated.features) == len(record.features)


class TestHomopolymerBoundaryCase:
    """Test homopolymer detection at boundaries."""

    def test_homopolymer_at_fragment_boundary(self):
        """Test detection of homopolymer at sequence boundaries.

        Regression test for v0.3.2 bug where homopolymers at
        fragment boundaries were missed.
        """
        # Homopolymer starts at position 0
        seq = "AAAAAAAAATGCGATCG"

        violations = check_homopolymer(seq, max_run=6)

        assert len(violations) == 1
        assert violations[0][0] == 0  # Position
        assert violations[0][1] == "A"  # Nucleotide

    def test_homopolymer_at_end(self):
        """Test detection of homopolymer at sequence end."""
        seq = "ATGCGATCGAAAAAAAA"

        violations = check_homopolymer(seq, max_run=6)

        assert len(violations) == 1
        assert violations[0][1] == "A"


class TestIdentical3PrimeEnds:
    """Test handling of fragments with identical 3' ends."""

    def test_identical_ends_overlap(self):
        """Test finding overlap when fragments have identical 3' ends.

        Regression test for v0.4.1 bug where identical 3' ends
        caused incorrect assembly order.
        """
        # Two fragments with same 3' end
        end_seq = "GCTAGCTAGCTAGCTAGCTA"
        frag1 = "AAAA" + end_seq
        frag2 = "TTTT" + end_seq

        # Should find the overlap correctly
        length1, overlap1 = find_overlap(frag1, end_seq, min_overlap=15)
        length2, overlap2 = find_overlap(frag2, end_seq, min_overlap=15)

        assert length1 == len(end_seq)
        assert length2 == len(end_seq)
        assert overlap1 == overlap2


class TestSequenceValidationEdgeCases:
    """Edge cases for sequence validation."""

    def test_very_short_sequence(self):
        """Test validation of very short sequence."""
        result = validate_sequence("ATGC", min_length=50)

        assert not result
        assert any("too short" in e for e in result.errors)

    def test_sequence_with_ambiguous_bases(self):
        """Test validation rejects ambiguous bases."""
        result = validate_sequence("ATGCNNNNATGC")

        assert not result
        assert any("Invalid" in e for e in result.errors)

    def test_internal_stop_codon(self):
        """Test detection of internal stop codons."""
        # Start codon, then stop codon in the middle, then more codons
        seq = "ATGTAAGCGATCGATCG"  # TAA at position 3-5

        result = validate_sequence(seq, check_expression=True)

        assert any("Internal stop" in e for e in result.errors)

    def test_extreme_gc_content(self):
        """Test detection of extreme GC content."""
        high_gc = "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"

        result = validate_sequence(high_gc)

        # Should have warning or error about high GC
        assert any("GC" in str(w) for w in result.warnings + result.errors)

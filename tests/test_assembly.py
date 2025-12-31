"""
Tests for core assembly functions.
"""

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from src.core.assembly import (
    assemble_plasmid,
    rotate_assembly,
    validate_overlaps,
    AssemblyError,
    OverlapError,
)


class TestAssemblePlasmid:
    """Tests for assemble_plasmid function."""

    def test_basic_assembly(self):
        """Test basic two-fragment assembly."""
        # Create backbone with overlap regions
        backbone_seq = "ATGCGATCGATCGATCGATCG" + "TAGCTAGCTAGCTAGCTAGCTAGC"
        backbone = SeqRecord(Seq(backbone_seq), id="backbone")

        # Create fragment that overlaps with backbone ends
        fragment_seq = "TAGCTAGCTAGCTAGCTAGCTAGC" + "AAAAAAAAAAAAAAAAA" + "ATGCGATCGATCGATCGATCG"
        fragment = SeqRecord(Seq(fragment_seq), id="fragment1")

        # Assemble
        result = assemble_plasmid(backbone, [fragment], "P-TEST", validate=False)

        assert result.id == "P-TEST"
        assert len(result.seq) > 0

    def test_empty_fragments_raises(self):
        """Test that empty fragment list raises error."""
        backbone = SeqRecord(Seq("ATGCGATCGATCG"), id="backbone")

        with pytest.raises(AssemblyError, match="No fragments"):
            assemble_plasmid(backbone, [], "P-TEST")

    def test_many_fragments_warning(self):
        """Test warning for >6 fragments."""
        backbone = SeqRecord(Seq("ATGC" * 100), id="backbone")
        fragments = [SeqRecord(Seq("ATGC" * 50), id=f"frag{i}") for i in range(7)]

        with pytest.warns(UserWarning, match="efficiency drops"):
            # This will fail validation but should still warn
            try:
                assemble_plasmid(backbone, fragments, "P-TEST", validate=False)
            except Exception:
                pass


class TestValidateOverlaps:
    """Tests for overlap validation."""

    def test_short_overlap_fails(self):
        """Test that overlaps shorter than min fail."""
        backbone = SeqRecord(Seq("ATGCGATCGATCG"), id="backbone")
        fragment = SeqRecord(Seq("GATCGATCGATCGATCG"), id="fragment")

        with pytest.raises(OverlapError, match="too short"):
            validate_overlaps([fragment], backbone, circular=False)

    def test_long_overlap_fails(self):
        """Test that overlaps longer than max fail."""
        # Create 80bp overlap (exceeds 70bp max)
        overlap = "ATGC" * 20  # 80 bp

        backbone = SeqRecord(Seq(overlap + "NNNNN"), id="backbone")
        fragment = SeqRecord(Seq("NNNNN" + overlap), id="fragment")

        with pytest.raises(OverlapError, match="too long"):
            validate_overlaps([fragment], backbone, circular=False)


class TestRotateAssembly:
    """Tests for rotate_assembly function."""

    def test_basic_rotation(self):
        """Test rotation to canonical start."""
        seq = "NNNNNATGCDEFGHIJKLMNOP"
        record = SeqRecord(Seq(seq), id="test")

        rotated = rotate_assembly(record, "ATGC")

        assert str(rotated.seq).startswith("ATGC")

    def test_rotation_not_found_raises(self):
        """Test that missing start sequence raises error."""
        record = SeqRecord(Seq("ATGCGATCG"), id="test")

        with pytest.raises(ValueError, match="not found"):
            rotate_assembly(record, "ZZZZZZ")

    def test_feature_positions_adjusted(self):
        """Test that feature positions are updated after rotation.

        This was a bug fixed in v0.3.2 (issue #12).
        """
        from Bio.SeqFeature import SeqFeature, FeatureLocation

        seq = "NNNNNATGCGATCGATCG"
        record = SeqRecord(Seq(seq), id="test")
        record.features.append(SeqFeature(
            FeatureLocation(5, 10),
            type="CDS",
            qualifiers={"gene": ["test"]}
        ))

        rotated = rotate_assembly(record, "ATGC")

        # Feature should now be at position 0
        assert len(rotated.features) > 0
        assert rotated.features[0].location.start == 0

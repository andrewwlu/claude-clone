"""
Tests for GenBank I/O functions.
"""

import pytest
import tempfile
from pathlib import Path

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from src.io.genbank import (
    read_seq_file,
    write_gb_file,
    create_dna_seqrecord,
    create_seqfeature,
    validate_genbank,
)


class TestCreateDNASeqRecord:
    """Tests for create_dna_seqrecord function."""

    def test_basic_creation(self):
        """Test basic SeqRecord creation."""
        record = create_dna_seqrecord(
            seq="ATGCGATCG",
            sid="TEST001",
            sdescr="Test plasmid"
        )

        assert record.id == "TEST001"
        assert str(record.seq) == "ATGCGATCG"
        assert record.description == "Test plasmid"

    def test_circular_topology(self):
        """Test circular topology annotation."""
        record = create_dna_seqrecord(
            seq="ATGC",
            sid="TEST",
            sdescr="Test",
            topology="circular"
        )

        assert record.annotations["topology"] == "circular"

    def test_linear_topology(self):
        """Test linear topology annotation."""
        record = create_dna_seqrecord(
            seq="ATGC",
            sid="TEST",
            sdescr="Test",
            topology="linear"
        )

        assert record.annotations["topology"] == "linear"

    def test_with_features(self):
        """Test adding features."""
        features = [
            {"start": 0, "end": 9, "type": "CDS", "label": "gene1"},
        ]

        record = create_dna_seqrecord(
            seq="ATGCGATCG",
            sid="TEST",
            sdescr="Test",
            sfeaturelist=features
        )

        # Should have source + 1 custom feature
        assert len(record.features) == 2
        assert record.features[1].type == "CDS"

    def test_source_feature_always_added(self):
        """Test that source feature is always added."""
        record = create_dna_seqrecord(
            seq="ATGC",
            sid="TEST",
            sdescr="Test"
        )

        source_features = [f for f in record.features if f.type == "source"]
        assert len(source_features) == 1


class TestCreateSeqFeature:
    """Tests for create_seqfeature function."""

    def test_basic_feature(self):
        """Test basic feature creation."""
        feature = create_seqfeature(
            start=0,
            end=100,
            strand=1,
            ftype="CDS",
            flabel="MyGene"
        )

        assert feature.type == "CDS"
        assert feature.location.start == 0
        assert feature.location.end == 100
        assert feature.location.strand == 1

    def test_reverse_strand(self):
        """Test reverse strand feature."""
        feature = create_seqfeature(
            start=0,
            end=100,
            strand=-1,
            ftype="CDS"
        )

        assert feature.location.strand == -1

    def test_cds_qualifiers(self):
        """Test CDS feature qualifiers."""
        feature = create_seqfeature(
            start=0,
            end=100,
            ftype="CDS",
            flabel="TestGene",
            ftranslation="MVLSPADKTNV"
        )

        assert "gene" in feature.qualifiers
        assert "product" in feature.qualifiers
        assert "translation" in feature.qualifiers


class TestWriteAndReadGenbankFile:
    """Tests for writing and reading GenBank files."""

    def test_write_and_read(self):
        """Test writing and reading back a GenBank file."""
        record = create_dna_seqrecord(
            seq="ATGCGATCGATCGATCGATCG",
            sid="TESTFILE",
            sdescr="Test file"
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            path = write_gb_file(record, tmpdir, "test.gb")

            assert path.exists()

            # Read back
            read_record = read_seq_file(path)

            assert read_record.id == "TESTFILE"
            assert str(read_record.seq) == "ATGCGATCGATCGATCGATCG"

    def test_auto_filename(self):
        """Test automatic filename generation."""
        record = SeqRecord(Seq("ATGC"), id="AUTONAME")

        with tempfile.TemporaryDirectory() as tmpdir:
            path = write_gb_file(record, tmpdir)

            assert path.name == "AUTONAME.gb"


class TestValidateGenbank:
    """Tests for GenBank validation."""

    def test_valid_file(self):
        """Test validation of valid GenBank file."""
        record = create_dna_seqrecord(
            seq="ATGCGATCG",
            sid="VALID",
            sdescr="Valid file"
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            path = write_gb_file(record, tmpdir)
            is_valid, errors, warnings = validate_genbank(path)

            assert is_valid
            assert len(errors) == 0

    def test_missing_file(self):
        """Test validation of non-existent file."""
        is_valid, errors, warnings = validate_genbank("/nonexistent/file.gb")

        assert not is_valid
        assert len(errors) > 0


class TestReadSeqFile:
    """Tests for reading sequence files."""

    def test_read_genbank(self):
        """Test reading GenBank file."""
        record = create_dna_seqrecord(
            seq="ATGC",
            sid="GB",
            sdescr="GenBank"
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            path = write_gb_file(record, tmpdir, "test.gb")
            read_record = read_seq_file(path)

            assert read_record.id == "GB"

    def test_read_nonexistent_raises(self):
        """Test reading non-existent file raises error."""
        with pytest.raises(FileNotFoundError):
            read_seq_file("/nonexistent/file.gb")

    def test_unknown_extension_raises(self):
        """Test unknown file extension raises error."""
        with tempfile.NamedTemporaryFile(suffix=".xyz") as f:
            with pytest.raises(ValueError, match="Unrecognized"):
                read_seq_file(f.name)

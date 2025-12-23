"""
GenBank file I/O operations.

This module handles reading and writing GenBank format sequence files,
as well as creating properly annotated SeqRecord objects.
"""

import os
from pathlib import Path
from datetime import datetime
from typing import List, Optional, Dict, Any, Union

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


def read_seq_file(file_path: Union[str, Path]) -> SeqRecord:
    """
    Read a sequence file (.dna, .gb, .gbk, .fasta).

    Args:
        file_path: Path to sequence file

    Returns:
        SeqRecord object

    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file format not recognized

    Example:
        >>> record = read_seq_file("data/backbones/pUC19.gb")
        >>> print(f"Loaded {record.id}: {len(record.seq)} bp")
    """
    file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")

    # Determine format from extension
    ext = file_path.suffix.lower()
    format_map = {
        ".gb": "genbank",
        ".gbk": "genbank",
        ".genbank": "genbank",
        ".dna": "genbank",  # SnapGene uses .dna but it's often GenBank format
        ".fasta": "fasta",
        ".fa": "fasta",
        ".fna": "fasta",
    }

    file_format = format_map.get(ext)
    if not file_format:
        raise ValueError(f"Unrecognized file extension: {ext}")

    # Handle .dna files (might be binary SnapGene format)
    if ext == ".dna":
        try:
            return SeqIO.read(file_path, "genbank")
        except Exception:
            # Try reading as SnapGene binary (requires snapgene_reader)
            try:
                from snapgene_reader import snapgene_file_to_seqrecord
                return snapgene_file_to_seqrecord(str(file_path))
            except ImportError:
                raise ValueError(
                    f"Cannot read SnapGene .dna file. "
                    "Install snapgene_reader: pip install snapgene_reader"
                )

    return SeqIO.read(file_path, file_format)


def write_gb_file(
    seqrecord: SeqRecord,
    out_dir: Union[str, Path],
    file_name: Optional[str] = None
) -> Path:
    """
    Write SeqRecord to GenBank file.

    Args:
        seqrecord: SeqRecord to write
        out_dir: Output directory
        file_name: Output filename (default: seqrecord.id + ".gb")

    Returns:
        Path to written file

    Example:
        >>> path = write_gb_file(plasmid, "data/plasmids/")
        >>> print(f"Wrote {path}")
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if file_name is None:
        file_name = f"{seqrecord.id}.gb"

    if not file_name.endswith((".gb", ".gbk")):
        file_name += ".gb"

    out_path = out_dir / file_name

    # Ensure required annotations
    if "molecule_type" not in seqrecord.annotations:
        seqrecord.annotations["molecule_type"] = "DNA"

    SeqIO.write(seqrecord, out_path, "genbank")

    return out_path


def create_dna_seqrecord(
    seq: str,
    sid: str,
    sdescr: str,
    sfeaturelist: Optional[List[Dict[str, Any]]] = None,
    topology: str = "circular"
) -> SeqRecord:
    """
    Create SeqRecord for GenBank output.

    Args:
        seq: DNA sequence string
        sid: Sequence ID
        sdescr: Sequence description
        sfeaturelist: List of feature dictionaries
        topology: "circular" or "linear"

    Returns:
        SeqRecord ready for GenBank output

    Example:
        >>> features = [
        ...     {"start": 0, "end": 100, "type": "promoter", "label": "T7"}
        ... ]
        >>> record = create_dna_seqrecord("ATGC...", "P-0654", "My plasmid", features)
    """
    # Create base record
    record = SeqRecord(
        Seq(seq.upper()),
        id=sid,
        name=sid[:16],  # GenBank name field limited to 16 chars
        description=sdescr,
        annotations={
            "molecule_type": "DNA",
            "topology": topology,
            "organism": "synthetic construct",
            "source": "synthetic construct",
            "date": datetime.now().strftime("%d-%b-%Y").upper(),
        }
    )

    # Add source feature (required for valid GenBank)
    source_feature = SeqFeature(
        FeatureLocation(0, len(seq)),
        type="source",
        qualifiers={
            "organism": ["synthetic construct"],
            "mol_type": ["other DNA"],
        }
    )
    record.features.append(source_feature)

    # Add provided features
    if sfeaturelist:
        for feat in sfeaturelist:
            feature = create_seqfeature(
                start=feat.get("start", 0),
                end=feat.get("end", len(seq)),
                strand=feat.get("strand", 1),
                fid=feat.get("id", ""),
                ftype=feat.get("type", "misc_feature"),
                flabel=feat.get("label", feat.get("id", "")),
                ftranslation=feat.get("translation")
            )
            record.features.append(feature)

    return record


def create_seqfeature(
    start: int,
    end: int,
    strand: int = 1,
    fid: str = "",
    ftype: Optional[str] = None,
    flabel: Optional[str] = None,
    ftranslation: Optional[str] = None,
    fnote: Optional[str] = None
) -> SeqFeature:
    """
    Create SeqFeature for GenBank annotations.

    Args:
        start: Start position (0-indexed)
        end: End position (exclusive)
        strand: 1 for forward, -1 for reverse
        fid: Feature ID
        ftype: Feature type (CDS, gene, etc.)
        flabel: Feature label/name
        ftranslation: Protein translation (for CDS)
        fnote: Additional note

    Returns:
        SeqFeature object

    Example:
        >>> cds = create_seqfeature(100, 400, 1, "gene1", "CDS", "MyGene")
    """
    if ftype is None:
        ftype = "misc_feature"

    qualifiers = {}

    if flabel:
        qualifiers["label"] = [flabel]
        if ftype == "CDS":
            qualifiers["gene"] = [flabel]
            qualifiers["product"] = [flabel]

    if ftranslation:
        qualifiers["translation"] = [ftranslation]

    if fnote:
        qualifiers["note"] = [fnote]

    if fid:
        qualifiers["locus_tag"] = [fid]

    return SeqFeature(
        FeatureLocation(start, end, strand=strand),
        type=ftype,
        qualifiers=qualifiers
    )


def validate_genbank(file_path: Union[str, Path]) -> tuple:
    """
    Validate a GenBank file.

    Args:
        file_path: Path to GenBank file

    Returns:
        Tuple of (is_valid, errors, warnings)
    """
    errors = []
    warnings = []

    try:
        record = SeqIO.read(file_path, "genbank")
    except Exception as e:
        return (False, [f"Failed to parse: {e}"], [])

    # Check required elements
    if not record.id:
        errors.append("Missing sequence ID")

    if not record.seq:
        errors.append("Missing sequence")

    if "molecule_type" not in record.annotations:
        errors.append("Missing molecule_type annotation")

    # Check for source feature
    has_source = any(f.type == "source" for f in record.features)
    if not has_source:
        errors.append("Missing source feature")

    # Check CDS features
    for feature in record.features:
        if feature.type == "CDS":
            if "translation" not in feature.qualifiers:
                warnings.append(f"CDS at {feature.location} missing translation")
            if "product" not in feature.qualifiers:
                warnings.append(f"CDS at {feature.location} missing product")

    return (len(errors) == 0, errors, warnings)

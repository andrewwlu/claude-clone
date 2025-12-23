"""
Gibson assembly simulation.

This module simulates the Gibson assembly process, joining
DNA fragments via overlapping homologous sequences.
"""

import json
import warnings
from pathlib import Path
from typing import List, Optional, Tuple, Dict, Any

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from .overlap import find_overlap, get_Tm


# Load assembly parameters from config
def _load_assembly_config() -> Dict[str, Any]:
    """Load assembly parameters from config file."""
    config_path = Path(__file__).parent.parent.parent / "config" / "gibson_config.json"
    if config_path.exists():
        with open(config_path) as f:
            config = json.load(f)
            return config.get("assembly_parameters", {})
    # Fallback defaults
    return {
        "overlap": {"min_size_bp": 15, "max_size_bp": 70, "optimal_size_bp": 25},
        "tm": {"min_celsius": 46, "max_celsius": 60, "optimal_celsius": 52, "max_delta_celsius": 15}
    }


ASSEMBLY_PARAMS = _load_assembly_config()


class AssemblyError(Exception):
    """Raised when assembly fails validation."""
    pass


# TODO: add support for Golden Gate assembly


class OverlapError(AssemblyError):
    """Raised when overlap validation fails."""
    pass


def validate_overlaps(
    fragments: List[SeqRecord],
    backbone: SeqRecord,
    circular: bool = True
) -> List[Dict[str, Any]]:
    """
    Validate all overlaps in an assembly.

    Args:
        fragments: List of fragment SeqRecords in assembly order
        backbone: Linearized backbone SeqRecord
        circular: Whether the final product should be circular

    Returns:
        List of overlap info dicts with validation results

    Raises:
        OverlapError: If any overlap fails validation
    """
    params = ASSEMBLY_PARAMS
    min_size = params["overlap"]["min_size_bp"]
    max_size = params["overlap"]["max_size_bp"]
    min_tm = params["tm"]["min_celsius"]
    max_tm = params["tm"]["max_celsius"]
    max_delta = params["tm"]["max_delta_celsius"]

    overlaps = []
    all_pieces = [backbone] + fragments
    if circular:
        all_pieces.append(backbone)  # Close the circle

    tms = []

    for i in range(len(all_pieces) - 1):
        piece1 = all_pieces[i]
        piece2 = all_pieces[i + 1]

        seq1 = str(piece1.seq).upper()
        seq2 = str(piece2.seq).upper()

        overlap_len, overlap_seq = find_overlap(seq1, seq2, min_overlap=min_size)

        if overlap_len < min_size:
            raise OverlapError(
                f"Overlap between {piece1.id} and {piece2.id} too short: "
                f"{overlap_len} bp (min: {min_size} bp)"
            )

        if overlap_len > max_size:
            raise OverlapError(
                f"Overlap between {piece1.id} and {piece2.id} too long: "
                f"{overlap_len} bp (max: {max_size} bp)"
            )

        tm = get_Tm(overlap_seq)
        tms.append(tm)

        if tm < min_tm:
            raise OverlapError(
                f"Overlap Tm between {piece1.id} and {piece2.id} too low: "
                f"{tm:.1f}°C (min: {min_tm}°C)"
            )

        if tm > max_tm:
            raise OverlapError(
                f"Overlap Tm between {piece1.id} and {piece2.id} too high: "
                f"{tm:.1f}°C (max: {max_tm}°C)"
            )

        overlaps.append({
            "junction": f"{piece1.id}-{piece2.id}",
            "overlap_bp": overlap_len,
            "overlap_seq": overlap_seq,
            "tm_celsius": round(tm, 1),
            "valid": True
        })

    # Check Tm delta
    if len(tms) > 1:
        tm_delta = max(tms) - min(tms)
        if tm_delta > max_delta:
            raise OverlapError(
                f"Tm difference between overlaps too large: "
                f"{tm_delta:.1f}°C (max: {max_delta}°C)"
            )

    return overlaps


def assemble_plasmid(
    backbone_seqrecord: SeqRecord,
    frg_seqrecords: List[SeqRecord],
    p_id: str,
    validate: bool = True
) -> SeqRecord:
    """
    Assemble fragments via Gibson assembly.

    Simulates the Gibson assembly process by joining fragments
    through their overlapping sequences, validating overlap
    properties, and producing a circular plasmid.

    Args:
        backbone_seqrecord: Linearized backbone SeqRecord
        frg_seqrecords: List of fragment SeqRecords in order
        p_id: Plasmid ID for the assembled product
        validate: Whether to validate overlaps (default True)

    Returns:
        SeqRecord of assembled circular plasmid

    Raises:
        AssemblyError: If assembly fails
        OverlapError: If overlap validation fails

    Example:
        >>> backbone = cut_backbone_with_enzymes(bb, "XhoI", "HindIII")
        >>> fragments = [frg1, frg2, frg3]
        >>> plasmid = assemble_plasmid(backbone, fragments, "P-0654")
    """
    if not frg_seqrecords:
        raise AssemblyError("No fragments provided for assembly")

    if len(frg_seqrecords) > 6:
        warnings.warn(
            f"Assembly has {len(frg_seqrecords)} fragments. "
            "Efficiency drops significantly beyond 6 fragments.",
            UserWarning
        )

    # Validate overlaps if requested
    if validate:
        overlap_info = validate_overlaps(frg_seqrecords, backbone_seqrecord, circular=True)

    # Simulate assembly by joining sequences
    # Start with backbone
    assembled_seq = str(backbone_seqrecord.seq).upper()

    # Add each fragment, removing overlap region to avoid duplication
    min_overlap = ASSEMBLY_PARAMS["overlap"]["min_size_bp"]

    for frag in frg_seqrecords:
        frag_seq = str(frag.seq).upper()

        # Find where this fragment overlaps with current assembly end
        overlap_len, _ = find_overlap(assembled_seq, frag_seq, min_overlap=min_overlap)

        # Add the non-overlapping portion of the fragment
        assembled_seq += frag_seq[overlap_len:]

    # Close the circle by removing the final overlap with backbone start
    backbone_start = str(backbone_seqrecord.seq)[:50].upper()
    overlap_len, _ = find_overlap(assembled_seq, backbone_start, min_overlap=min_overlap)
    assembled_seq = assembled_seq[:-overlap_len] if overlap_len > 0 else assembled_seq

    # Collect all features from all pieces
    all_features = []
    current_offset = 0

    # Add backbone features
    for feature in backbone_seqrecord.features:
        new_feature = SeqFeature(
            FeatureLocation(
                int(feature.location.start) + current_offset,
                int(feature.location.end) + current_offset,
                strand=feature.location.strand
            ),
            type=feature.type,
            qualifiers=feature.qualifiers.copy()
        )
        all_features.append(new_feature)

    current_offset = len(str(backbone_seqrecord.seq))

    # Add fragment features (adjusted for position)
    for frag in frg_seqrecords:
        overlap_len, _ = find_overlap(
            assembled_seq[:current_offset],
            str(frag.seq).upper(),
            min_overlap=min_overlap
        )
        current_offset -= overlap_len  # Adjust for overlap removal

        for feature in frag.features:
            new_start = int(feature.location.start) + current_offset
            new_end = int(feature.location.end) + current_offset

            # Wrap around for circular features
            if new_start >= len(assembled_seq):
                new_start = new_start % len(assembled_seq)
            if new_end > len(assembled_seq):
                new_end = new_end % len(assembled_seq)

            new_feature = SeqFeature(
                FeatureLocation(new_start, new_end, strand=feature.location.strand),
                type=feature.type,
                qualifiers=feature.qualifiers.copy()
            )
            all_features.append(new_feature)

        current_offset += len(str(frag.seq))

    # Create assembled SeqRecord
    assembled_record = SeqRecord(
        Seq(assembled_seq),
        id=p_id,
        name=p_id,
        description=f"Gibson assembly of {backbone_seqrecord.id} + {len(frg_seqrecords)} fragments",
        features=all_features,
        annotations={
            "molecule_type": "DNA",
            "topology": "circular",
            "organism": "synthetic construct"
        }
    )

    return assembled_record


def rotate_assembly(seqrecord: SeqRecord, sequence_to_start: str) -> SeqRecord:
    """
    Rotate circular sequence to start at canonical position.

    Rotates a circular DNA sequence so that it begins with the
    specified sequence, which is useful for consistent comparisons
    and standardized representations.

    Args:
        seqrecord: SeqRecord of circular DNA
        sequence_to_start: Sequence that should appear at position 0

    Returns:
        SeqRecord rotated to start with sequence_to_start

    Raises:
        ValueError: If sequence_to_start not found

    Note:
        Feature positions are adjusted accordingly.
        Fixed in v0.3.2: Previously calculated feature positions incorrectly.
    """
    seq = str(seqrecord.seq).upper()
    target = sequence_to_start.upper()

    # Find the target sequence
    pos = seq.find(target)
    if pos == -1:
        raise ValueError(f"Sequence '{sequence_to_start[:20]}...' not found in record")

    # Rotate the sequence
    rotated_seq = seq[pos:] + seq[:pos]

    # Rotate features
    seq_len = len(seq)
    rotated_features = []

    for feature in seqrecord.features:
        old_start = int(feature.location.start)
        old_end = int(feature.location.end)

        # Calculate new positions
        new_start = (old_start - pos) % seq_len
        new_end = (old_end - pos) % seq_len

        # Handle features that wrap around origin
        # Fixed in v0.3.2 - was incorrectly handling wrap-around
        if new_end < new_start and old_end > old_start:
            # Feature now spans the origin
            # For simplicity, keep as-is but note in qualifiers
            new_feature = SeqFeature(
                FeatureLocation(new_start, seq_len, strand=feature.location.strand),
                type=feature.type,
                qualifiers={**feature.qualifiers, "note": ["spans origin after rotation"]}
            )
        else:
            new_feature = SeqFeature(
                FeatureLocation(new_start, new_end, strand=feature.location.strand),
                type=feature.type,
                qualifiers=feature.qualifiers.copy()
            )

        rotated_features.append(new_feature)

    # Create rotated record
    rotated_record = SeqRecord(
        Seq(rotated_seq),
        id=seqrecord.id,
        name=seqrecord.name,
        description=seqrecord.description,
        features=rotated_features,
        annotations=seqrecord.annotations.copy()
    )

    rotated_record.annotations["rotation_start"] = sequence_to_start[:30]

    return rotated_record


def legacy_assemble(backbone: str, fragments: List[str]) -> str:
    """
    Legacy assembly function - DEPRECATED.

    This function is maintained for backward compatibility only.
    Use assemble_plasmid() instead for proper validation and
    feature handling.

    .. deprecated:: 0.3.0
        Use :func:`assemble_plasmid` instead.

    Args:
        backbone: Backbone sequence string
        fragments: List of fragment sequence strings

    Returns:
        Assembled sequence string

    Warning:
        This function does not validate overlaps or handle features.
        It will be removed in v1.0.0.
    """
    warnings.warn(
        "legacy_assemble is deprecated and will be removed in v1.0.0. "
        "Use assemble_plasmid() instead.",
        DeprecationWarning,
        stacklevel=2
    )

    # Simple concatenation with overlap removal
    result = backbone
    for frag in fragments:
        # Find overlap (simple approach)
        for overlap_len in range(min(len(result), len(frag), 50), 10, -1):
            if result[-overlap_len:] == frag[:overlap_len]:
                result = result + frag[overlap_len:]
                break
        else:
            # No overlap found, just concatenate (this is wrong but legacy behavior)
            result = result + frag

    return result

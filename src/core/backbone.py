"""
Backbone cutting with restriction enzymes.

This module handles linearization of circular backbone plasmids
using restriction enzyme digestion for Gibson assembly.
"""

import json
from pathlib import Path
from typing import Optional, Tuple, List, Dict, Any

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


# Load enzyme definitions from config
def _load_enzyme_config() -> Dict[str, Dict[str, str]]:
    """Load restriction enzyme definitions from config file."""
    config_path = Path(__file__).parent.parent.parent / "config" / "gibson_config.json"
    if config_path.exists():
        with open(config_path) as f:
            config = json.load(f)
            return config.get("restriction_enzymes", {})
    # Fallback to hardcoded defaults
    return {
        "XhoI": {"recognition": "CTCGAG", "cut": "C/TCGAG"},
        "HindIII": {"recognition": "AAGCTT", "cut": "A/AGCTT"},
        "EcoRI": {"recognition": "GAATTC", "cut": "G/AATTC"},
        "BamHI": {"recognition": "GGATCC", "cut": "G/GATCC"},
    }


RESTRICTION_ENZYMES = _load_enzyme_config()


def find_cut_sites(sequence: str, enzyme: str) -> List[int]:
    """
    Find all cut sites for a restriction enzyme in a sequence.

    Args:
        sequence: DNA sequence string
        enzyme: Enzyme name (e.g., 'XhoI')

    Returns:
        List of cut positions (0-indexed)

    Raises:
        ValueError: If enzyme is not recognized
    """
    if enzyme not in RESTRICTION_ENZYMES:
        raise ValueError(f"Unknown enzyme: {enzyme}. Available: {list(RESTRICTION_ENZYMES.keys())}")

    recognition = RESTRICTION_ENZYMES[enzyme]["recognition"]
    cut_pattern = RESTRICTION_ENZYMES[enzyme]["cut"]

    # Find cut position within recognition sequence
    cut_offset = cut_pattern.index("/")

    # Find all occurrences
    sites = []
    pos = 0
    seq_upper = sequence.upper()
    while True:
        pos = seq_upper.find(recognition, pos)
        if pos == -1:
            break
        sites.append(pos + cut_offset)
        pos += 1

    return sites


def cut_backbone_with_enzymes(
    backbone_seqrecord: SeqRecord,
    enzyme1: str = "XhoI",
    enzyme2: str = "HindIII",
    min_overlap: int = 20
) -> SeqRecord:
    """
    Cut backbone with two restriction enzymes.

    Linearizes a circular backbone by cutting with two enzymes,
    removing the region between them to create an insertion site
    for Gibson assembly.

    Args:
        backbone_seqrecord: SeqRecord of circular backbone
        enzyme1: First restriction enzyme name
        enzyme2: Second restriction enzyme name
        min_overlap: Minimum distance between cut sites

    Returns:
        SeqRecord of linearized backbone

    Raises:
        ValueError: If enzymes don't cut exactly once each

    Note:
        Handles palindromic enzymes correctly.
        Feature positions are adjusted for the linear product.

    Example:
        >>> backbone = SeqIO.read("pUC19.gb", "genbank")
        >>> linear = cut_backbone_with_enzymes(backbone, "XhoI", "HindIII")
    """
    sequence = str(backbone_seqrecord.seq).upper()

    # Find cut sites
    sites1 = find_cut_sites(sequence, enzyme1)
    sites2 = find_cut_sites(sequence, enzyme2)

    # Validate single cuts
    if len(sites1) == 0:
        raise ValueError(f"{enzyme1} does not cut this backbone")
    if len(sites1) > 1:
        raise ValueError(f"{enzyme1} cuts multiple times at positions: {sites1}")
    if len(sites2) == 0:
        raise ValueError(f"{enzyme2} does not cut this backbone")
    if len(sites2) > 1:
        raise ValueError(f"{enzyme2} cuts multiple times at positions: {sites2}")

    cut1 = sites1[0]
    cut2 = sites2[0]

    # Ensure cut1 < cut2 for consistent handling
    if cut1 > cut2:
        cut1, cut2 = cut2, cut1
        enzyme1, enzyme2 = enzyme2, enzyme1

    # Check minimum distance
    if cut2 - cut1 < min_overlap:
        raise ValueError(
            f"Cut sites too close ({cut2 - cut1} bp). "
            f"Minimum required: {min_overlap} bp"
        )

    # Create linearized sequence
    # Keep: [0:cut1] + [cut2:end]
    # Remove: [cut1:cut2] (the insert dropout)
    linear_seq = sequence[:cut1] + sequence[cut2:]

    # Adjust features for new coordinates
    adjusted_features = []
    dropout_size = cut2 - cut1

    for feature in backbone_seqrecord.features:
        start = int(feature.location.start)
        end = int(feature.location.end)

        # Skip features in the dropout region
        if start >= cut1 and end <= cut2:
            continue

        # Adjust features after the dropout
        if start >= cut2:
            new_start = start - dropout_size
            new_end = end - dropout_size
        elif end <= cut1:
            new_start = start
            new_end = end
        else:
            # Feature spans the cut site - this is complex
            # For now, skip it (TODO: handle split features properly)
            # See issue #8
            continue

        new_feature = SeqFeature(
            FeatureLocation(new_start, new_end, strand=feature.location.strand),
            type=feature.type,
            qualifiers=feature.qualifiers.copy()
        )
        adjusted_features.append(new_feature)

    # Create new SeqRecord
    linear_record = SeqRecord(
        Seq(linear_seq),
        id=backbone_seqrecord.id + "_linear",
        name=backbone_seqrecord.name,
        description=f"{backbone_seqrecord.description} linearized with {enzyme1}/{enzyme2}",
        features=adjusted_features,
        annotations=backbone_seqrecord.annotations.copy()
    )

    # Update topology annotation
    linear_record.annotations["topology"] = "linear"

    # Add metadata about the cut
    linear_record.annotations["linearization"] = {
        "enzyme1": enzyme1,
        "enzyme2": enzyme2,
        "cut1_position": cut1,
        "cut2_position": cut2,
        "dropout_size": dropout_size
    }

    return linear_record


def get_overhang(enzyme: str) -> Tuple[str, str]:
    """
    Get the 5' and 3' overhangs produced by an enzyme.

    Args:
        enzyme: Enzyme name

    Returns:
        Tuple of (5' overhang, 3' overhang)
        Empty string if blunt end
    """
    if enzyme not in RESTRICTION_ENZYMES:
        raise ValueError(f"Unknown enzyme: {enzyme}")

    cut_pattern = RESTRICTION_ENZYMES[enzyme]["cut"]
    recognition = RESTRICTION_ENZYMES[enzyme]["recognition"]
    cut_pos = cut_pattern.index("/")

    # For most enzymes, the overhang is the uncut portion
    if cut_pos == 0:
        # Cuts at very beginning - 5' overhang
        return (recognition, "")
    elif cut_pos == len(recognition):
        # Cuts at very end - 3' overhang
        return ("", recognition)
    else:
        # Internal cut - creates sticky ends
        five_prime = recognition[cut_pos:]
        return (five_prime, "")

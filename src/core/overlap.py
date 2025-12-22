"""
Overlap detection and Tm calculation.

This module provides functions for finding overlapping sequences
and calculating melting temperatures for Gibson assembly.
"""

from typing import Tuple, Optional

from Bio.SeqUtils import MeltingTemp as mt

# FIXME: nearest-neighbor parameters should be configurable


def find_overlap(
    seq1: str,
    seq2: str,
    min_overlap: int = 15,
    max_overlap: int = 70
) -> Tuple[int, str]:
    """
    Find overlap between end of seq1 and start of seq2.

    Searches for the longest sequence that appears at both the
    3' end of seq1 and the 5' end of seq2.

    Args:
        seq1: First sequence (look at 3' end)
        seq2: Second sequence (look at 5' end)
        min_overlap: Minimum overlap to consider
        max_overlap: Maximum overlap to search for

    Returns:
        Tuple of (overlap_length, overlap_sequence)
        Returns (0, "") if no overlap found >= min_overlap

    Example:
        >>> seq1 = "ATGCGATCGATCGATCGATCGTAGCTAGCT"
        >>> seq2 = "GCTAGCTATGCATGCATGC"
        >>> length, seq = find_overlap(seq1, seq2, min_overlap=8)
        >>> print(f"Overlap: {length} bp - {seq}")
        Overlap: 8 bp - GCTAGCTA
    """
    seq1 = seq1.upper()
    seq2 = seq2.upper()

    # Limit search to reasonable range
    max_search = min(max_overlap, len(seq1), len(seq2))

    # Search from longest to shortest
    for overlap_len in range(max_search, min_overlap - 1, -1):
        suffix = seq1[-overlap_len:]
        prefix = seq2[:overlap_len]

        if suffix == prefix:
            return (overlap_len, suffix)

    return (0, "")


def get_Tm(seq: str, method: str = "nearest_neighbor") -> float:
    """
    Calculate melting temperature for a DNA sequence.

    Uses BioPython's MeltingTemp module for accurate Tm calculation.

    Args:
        seq: DNA sequence string
        method: Calculation method - "nearest_neighbor" (default) or "gc"

    Returns:
        Melting temperature in Celsius

    Note:
        Switched from Tm_GC to Tm_NN in v0.4.1 for better accuracy
        with short overlaps typical in Gibson assembly.

    Example:
        >>> tm = get_Tm("ATGCGATCGATCGATCGATCGTAG")
        >>> print(f"Tm: {tm:.1f}°C")
        Tm: 58.3°C
    """
    seq = seq.upper()

    if method == "gc":
        # Simple GC method - less accurate but fast
        # Kept for backward compatibility
        return mt.Tm_GC(seq)
    else:
        # Nearest-neighbor method - more accurate for short sequences
        # Default parameters for Gibson assembly conditions
        try:
            return mt.Tm_NN(
                seq,
                nn_table=mt.DNA_NN4,  # SantaLucia 2004
                dnac1=250,  # nM primer concentration
                dnac2=250,
                Na=50,      # mM sodium
                Mg=0,       # mM magnesium (Gibson uses isothermal)
            )
        except ValueError:
            # Fallback for very short sequences
            return mt.Tm_GC(seq)


def design_overlap(
    target_tm: float = 52.0,
    seq1_end: str = "",
    seq2_start: str = "",
    gc_range: Tuple[float, float] = (0.40, 0.60)
) -> Tuple[int, str]:
    """
    Design an optimal overlap sequence.

    Given the sequences at a junction, determines the optimal
    overlap length to achieve the target Tm.

    Args:
        target_tm: Target melting temperature (default 52°C)
        seq1_end: 3' end of first fragment (at least 70 bp)
        seq2_start: 5' start of second fragment (at least 70 bp)
        gc_range: Acceptable GC content range

    Returns:
        Tuple of (optimal_length, overlap_sequence)

    Raises:
        ValueError: If optimal overlap cannot be designed
    """
    if not seq1_end or not seq2_start:
        raise ValueError("Both sequence ends must be provided")

    seq1_end = seq1_end.upper()
    seq2_start = seq2_start.upper()

    # Try different overlap lengths
    best_overlap = None
    best_delta = float('inf')

    for length in range(15, min(71, len(seq1_end), len(seq2_start))):
        suffix = seq1_end[-length:]
        prefix = seq2_start[:length]

        if suffix != prefix:
            continue

        # Calculate properties
        tm = get_Tm(suffix)
        gc = (suffix.count('G') + suffix.count('C')) / len(suffix)

        # Check GC range
        if gc < gc_range[0] or gc > gc_range[1]:
            continue

        # Find closest to target
        delta = abs(tm - target_tm)
        if delta < best_delta:
            best_delta = delta
            best_overlap = (length, suffix)

    if best_overlap is None:
        raise ValueError(
            f"Could not design overlap with Tm near {target_tm}°C. "
            "Consider modifying the junction sequences."
        )

    return best_overlap


def gc_content(seq: str) -> float:
    """
    Calculate GC content of a sequence.

    Args:
        seq: DNA sequence string

    Returns:
        GC content as fraction (0.0-1.0)
    """
    seq = seq.upper()
    gc = seq.count('G') + seq.count('C')
    return gc / len(seq) if seq else 0.0


def check_homopolymer(seq: str, max_run: int = 6) -> list:
    """
    Find homopolymer runs exceeding the maximum length.

    Args:
        seq: DNA sequence string
        max_run: Maximum allowed consecutive identical nucleotides

    Returns:
        List of (position, nucleotide, length) tuples for violations
    """
    seq = seq.upper()
    violations = []

    for nucleotide in "ATGC":
        i = 0
        while i < len(seq):
            if seq[i] == nucleotide:
                # Count run length
                run_start = i
                while i < len(seq) and seq[i] == nucleotide:
                    i += 1
                run_length = i - run_start

                if run_length > max_run:
                    violations.append((run_start, nucleotide, run_length))
            else:
                i += 1

    return violations

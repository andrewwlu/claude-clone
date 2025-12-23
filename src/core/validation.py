"""
Sequence and assembly validation.

This module provides comprehensive validation for DNA sequences
and Gibson assembly designs.
"""

import re
from typing import List, Dict, Any, Tuple, Optional

from .overlap import get_Tm, gc_content, check_homopolymer, find_overlap


# TODO: add JSON schema validation for config files

class ValidationResult:
    """Container for validation results."""

    def __init__(self, passed: bool, errors: List[str], warnings: List[str]):
        self.passed = passed
        self.errors = errors
        self.warnings = warnings

    def __bool__(self):
        return self.passed

    def __repr__(self):
        status = "PASS" if self.passed else "FAIL"
        return f"ValidationResult({status}, errors={len(self.errors)}, warnings={len(self.warnings)})"


def validate_sequence(
    sequence: str,
    check_synthesis: bool = True,
    check_expression: bool = False,
    min_length: int = 50,
    max_length: int = 10000
) -> ValidationResult:
    """
    Validate a DNA sequence for cloning and synthesis.

    Args:
        sequence: DNA sequence string
        check_synthesis: Check synthesis constraints (default True)
        check_expression: Check expression requirements (default False)
        min_length: Minimum sequence length
        max_length: Maximum sequence length

    Returns:
        ValidationResult with pass/fail status and messages

    Example:
        >>> result = validate_sequence("ATGCGATCGATCG")
        >>> if result:
        ...     print("Sequence is valid")
        >>> else:
        ...     print("Errors:", result.errors)
    """
    errors = []
    warnings = []
    seq = sequence.upper()

    # Basic validation
    if not seq:
        errors.append("Empty sequence")
        return ValidationResult(False, errors, warnings)

    # Check valid characters
    invalid_chars = set(seq) - set("ATGC")
    if invalid_chars:
        errors.append(f"Invalid characters: {invalid_chars}")

    # Check length
    if len(seq) < min_length:
        errors.append(f"Sequence too short: {len(seq)} bp (min: {min_length})")
    if len(seq) > max_length:
        errors.append(f"Sequence too long: {len(seq)} bp (max: {max_length})")

    # GC content
    gc = gc_content(seq)
    if gc < 0.25:
        errors.append(f"GC content too low: {gc*100:.1f}% (min: 25%)")
    elif gc < 0.40:
        warnings.append(f"GC content low: {gc*100:.1f}% (prefer 40-60%)")
    elif gc > 0.75:
        errors.append(f"GC content too high: {gc*100:.1f}% (max: 75%)")
    elif gc > 0.60:
        warnings.append(f"GC content high: {gc*100:.1f}% (prefer 40-60%)")

    if check_synthesis:
        # Homopolymer runs
        homopolymers = check_homopolymer(seq, max_run=6)
        for pos, nuc, length in homopolymers:
            if length > 8:
                errors.append(f"Long homopolymer: {nuc}x{length} at position {pos}")
            else:
                warnings.append(f"Homopolymer run: {nuc}x{length} at position {pos}")

        # Check for problematic repeats
        # FIXME: This is a simplified check - should use more sophisticated algorithm
        for repeat_len in [10, 12, 15]:
            for i in range(len(seq) - repeat_len * 2):
                subseq = seq[i:i + repeat_len]
                if seq.count(subseq) > 1:
                    warnings.append(f"Repeated sequence ({repeat_len} bp) at position {i}")
                    break

    if check_expression:
        # Check start codon
        if not seq.startswith("ATG"):
            warnings.append("No start codon (ATG) at beginning")

        # Check stop codon
        if seq[-3:] not in ["TAA", "TAG", "TGA"]:
            warnings.append("No stop codon at end")

        # Check reading frame
        if len(seq) % 3 != 0:
            warnings.append(f"Length not divisible by 3: {len(seq)} bp")

        # Check for internal stops
        for i in range(0, len(seq) - 3, 3):
            codon = seq[i:i+3]
            if codon in ["TAA", "TAG", "TGA"] and i < len(seq) - 3:
                errors.append(f"Internal stop codon ({codon}) at position {i}")

    passed = len(errors) == 0
    return ValidationResult(passed, errors, warnings)


def validate_overlap(
    overlap_seq: str,
    min_size: int = 15,
    max_size: int = 70,
    min_tm: float = 46.0,
    max_tm: float = 60.0
) -> ValidationResult:
    """
    Validate an overlap sequence for Gibson assembly.

    Args:
        overlap_seq: Overlap sequence string
        min_size: Minimum overlap length
        max_size: Maximum overlap length
        min_tm: Minimum melting temperature
        max_tm: Maximum melting temperature

    Returns:
        ValidationResult with pass/fail status and messages
    """
    errors = []
    warnings = []
    seq = overlap_seq.upper()

    # Check length
    if len(seq) < min_size:
        errors.append(f"Overlap too short: {len(seq)} bp (min: {min_size})")
    elif len(seq) < 20:
        warnings.append(f"Overlap on short side: {len(seq)} bp (prefer 25+ bp)")

    if len(seq) > max_size:
        errors.append(f"Overlap too long: {len(seq)} bp (max: {max_size})")

    # Check Tm
    if len(seq) >= 10:  # Need minimum length for Tm calculation
        tm = get_Tm(seq)
        if tm < min_tm:
            errors.append(f"Overlap Tm too low: {tm:.1f}°C (min: {min_tm}°C)")
        elif tm < 48:
            warnings.append(f"Overlap Tm low: {tm:.1f}°C (prefer 48-56°C)")

        if tm > max_tm:
            errors.append(f"Overlap Tm too high: {tm:.1f}°C (max: {max_tm}°C)")
        elif tm > 56:
            warnings.append(f"Overlap Tm high: {tm:.1f}°C (prefer 48-56°C)")

    # Check GC content
    gc = gc_content(seq)
    if gc < 0.30:
        errors.append(f"Overlap GC too low: {gc*100:.1f}%")
    elif gc < 0.40:
        warnings.append(f"Overlap GC low: {gc*100:.1f}%")

    if gc > 0.70:
        errors.append(f"Overlap GC too high: {gc*100:.1f}%")
    elif gc > 0.60:
        warnings.append(f"Overlap GC high: {gc*100:.1f}%")

    # Check for homopolymers
    for nuc in "ATGC":
        if nuc * 5 in seq:
            warnings.append(f"Homopolymer in overlap: {nuc}x5+")

    # Check for self-complementarity (hairpin risk)
    if len(seq) >= 16:
        half = len(seq) // 2
        first_half = seq[:half]
        second_half = seq[-half:]
        # Simple complement check - not full hairpin analysis
        complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
        rev_comp = "".join(complement.get(n, n) for n in reversed(second_half))
        matching = sum(1 for a, b in zip(first_half, rev_comp) if a == b)
        if matching > half * 0.7:
            warnings.append("Potential hairpin structure in overlap")

    passed = len(errors) == 0
    return ValidationResult(passed, errors, warnings)


def validate_assembly(
    backbone_seq: str,
    fragment_seqs: List[str],
    circular: bool = True
) -> ValidationResult:
    """
    Validate a complete Gibson assembly design.

    Args:
        backbone_seq: Linearized backbone sequence
        fragment_seqs: List of fragment sequences in order
        circular: Whether final product should be circular

    Returns:
        ValidationResult with pass/fail status and messages
    """
    errors = []
    warnings = []

    all_seqs = [backbone_seq] + fragment_seqs
    if circular:
        all_seqs.append(backbone_seq)  # Close the circle

    # Validate each junction
    tms = []
    for i in range(len(all_seqs) - 1):
        seq1 = all_seqs[i]
        seq2 = all_seqs[i + 1]

        overlap_len, overlap_seq = find_overlap(seq1, seq2, min_overlap=10)

        if overlap_len < 15:
            errors.append(f"Junction {i+1}: Overlap too short ({overlap_len} bp)")
        elif overlap_len < 20:
            warnings.append(f"Junction {i+1}: Overlap on short side ({overlap_len} bp)")

        if overlap_len >= 10:
            tm = get_Tm(overlap_seq)
            tms.append(tm)

            if tm < 46:
                errors.append(f"Junction {i+1}: Tm too low ({tm:.1f}°C)")
            elif tm > 60:
                errors.append(f"Junction {i+1}: Tm too high ({tm:.1f}°C)")

    # Check Tm balance
    if len(tms) > 1:
        tm_delta = max(tms) - min(tms)
        if tm_delta > 15:
            errors.append(f"Tm imbalance: {tm_delta:.1f}°C difference between junctions")
        elif tm_delta > 10:
            warnings.append(f"Consider balancing Tm: {tm_delta:.1f}°C difference")

    # Check fragment count
    if len(fragment_seqs) > 6:
        warnings.append(f"High fragment count ({len(fragment_seqs)}): efficiency may be reduced")

    # Check total length
    total_length = len(backbone_seq) + sum(len(f) for f in fragment_seqs)
    # Subtract approximate overlap lengths
    total_length -= (len(fragment_seqs) + 1) * 25  # Rough estimate

    if total_length > 15000:
        warnings.append(f"Large assembly ({total_length} bp): consider verification strategy")

    passed = len(errors) == 0
    return ValidationResult(passed, errors, warnings)

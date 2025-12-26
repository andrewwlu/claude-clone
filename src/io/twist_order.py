"""
Twist Bioscience order generation.

This module handles creating order files for Twist gene fragment
synthesis, including validation against Twist constraints.
"""

import csv
import json
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Any, Optional, Tuple

from ..core.overlap import gc_content, check_homopolymer


# Twist constraints
TWIST_MIN_LENGTH = 300  # bp
TWIST_MAX_LENGTH = 1800  # bp
TWIST_PREFERRED_GC_MIN = 0.25
TWIST_PREFERRED_GC_MAX = 0.75
TWIST_HOMOPOLYMER_LIMIT = 6

# Pricing (approximate, subject to change)
# TODO: fetch current pricing from Twist API
TWIST_PRICE_PER_BP = 0.09  # Standard
TWIST_PRICE_PER_BP_EXPRESS = 0.15  # Express (5-day)


def validate_twist_constraints(
    sequence: str,
    fragment_id: str = ""
) -> Tuple[bool, List[str], List[str]]:
    """
    Validate a sequence against Twist synthesis constraints.

    Args:
        sequence: DNA sequence
        fragment_id: Fragment identifier for error messages

    Returns:
        Tuple of (is_valid, errors, warnings)

    Example:
        >>> valid, errors, warnings = validate_twist_constraints("ATGC...")
        >>> if not valid:
        ...     print("Cannot order:", errors)
    """
    errors = []
    warnings = []
    seq = sequence.upper()
    prefix = f"[{fragment_id}] " if fragment_id else ""

    # Check length
    if len(seq) < TWIST_MIN_LENGTH:
        errors.append(f"{prefix}Too short: {len(seq)} bp (min: {TWIST_MIN_LENGTH} bp)")
    elif len(seq) > TWIST_MAX_LENGTH:
        errors.append(f"{prefix}Too long: {len(seq)} bp (max: {TWIST_MAX_LENGTH} bp)")

    # Check GC content
    gc = gc_content(seq)
    if gc < TWIST_PREFERRED_GC_MIN:
        errors.append(f"{prefix}GC too low: {gc*100:.1f}% (min: {TWIST_PREFERRED_GC_MIN*100}%)")
    elif gc < 0.30:
        warnings.append(f"{prefix}GC low: {gc*100:.1f}% (may affect yield)")

    if gc > TWIST_PREFERRED_GC_MAX:
        errors.append(f"{prefix}GC too high: {gc*100:.1f}% (max: {TWIST_PREFERRED_GC_MAX*100}%)")
    elif gc > 0.70:
        warnings.append(f"{prefix}GC high: {gc*100:.1f}% (may affect yield)")

    # Check homopolymers
    homopolymers = check_homopolymer(seq, max_run=TWIST_HOMOPOLYMER_LIMIT)
    for pos, nuc, length in homopolymers:
        if length > 10:
            errors.append(f"{prefix}Homopolymer {nuc}x{length} at pos {pos}")
        else:
            warnings.append(f"{prefix}Homopolymer {nuc}x{length} at pos {pos} (may affect synthesis)")

    # Check for valid nucleotides
    invalid = set(seq) - set("ATGC")
    if invalid:
        errors.append(f"{prefix}Invalid nucleotides: {invalid}")

    return (len(errors) == 0, errors, warnings)


def calculate_cost(
    fragments: List[Dict[str, Any]],
    express: bool = False
) -> Dict[str, Any]:
    """
    Calculate estimated cost for a Twist order.

    Args:
        fragments: List of fragment dicts with 'sequence' key
        express: Use express pricing (5-day delivery)

    Returns:
        Cost breakdown dictionary
    """
    price_per_bp = TWIST_PRICE_PER_BP_EXPRESS if express else TWIST_PRICE_PER_BP

    total_bp = 0
    fragment_costs = []

    for frag in fragments:
        seq = frag.get("sequence", "")
        bp = len(seq)
        cost = bp * price_per_bp
        total_bp += bp
        fragment_costs.append({
            "id": frag.get("id", frag.get("name", "unknown")),
            "bp": bp,
            "cost": round(cost, 2)
        })

    return {
        "fragments": fragment_costs,
        "total_bp": total_bp,
        "price_per_bp": price_per_bp,
        "total_cost": round(total_bp * price_per_bp, 2),
        "pricing_type": "express" if express else "standard",
        "estimated_delivery": "5 business days" if express else "13-17 business days"
    }


def generate_twist_order(
    fragments: List[Dict[str, Any]],
    order_id: Optional[str] = None,
    output_dir: Optional[Path] = None,
    validate: bool = True
) -> Dict[str, Any]:
    """
    Generate Twist order CSV file.

    Args:
        fragments: List of fragment dicts with 'id', 'sequence', optional 'notes'
        order_id: Order identifier (default: generated from timestamp)
        output_dir: Output directory (default: data/orders/)
        validate: Validate sequences before generating (default: True)

    Returns:
        Order result with file path and validation results

    Example:
        >>> frags = [{"id": "FRG-001", "sequence": "ATGC...", "notes": "Gene X"}]
        >>> result = generate_twist_order(frags, "ORDER-2025-001")
        >>> print(f"Order saved to: {result['file_path']}")
    """
    if order_id is None:
        order_id = f"twist_order_{datetime.now().strftime('%Y%m%d_%H%M%S')}"

    if output_dir is None:
        output_dir = Path(__file__).parent.parent.parent / "data" / "orders"

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Validate if requested
    validation_results = []
    all_valid = True

    if validate:
        for frag in fragments:
            frag_id = frag.get("id", frag.get("name", "unknown"))
            seq = frag.get("sequence", "")
            valid, errors, warnings = validate_twist_constraints(seq, frag_id)
            validation_results.append({
                "id": frag_id,
                "valid": valid,
                "errors": errors,
                "warnings": warnings
            })
            if not valid:
                all_valid = False

    if not all_valid:
        return {
            "success": False,
            "error": "Validation failed for one or more fragments",
            "validation": validation_results
        }

    # Generate CSV
    csv_path = output_dir / f"{order_id}.csv"

    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Name", "Sequence", "Notes"])

        for frag in fragments:
            name = frag.get("id", frag.get("name", ""))
            sequence = frag.get("sequence", "")
            notes = frag.get("notes", frag.get("description", ""))
            writer.writerow([name, sequence, notes])

    # Calculate cost
    cost_info = calculate_cost(fragments)

    return {
        "success": True,
        "order_id": order_id,
        "file_path": str(csv_path),
        "fragment_count": len(fragments),
        "total_bp": cost_info["total_bp"],
        "estimated_cost": cost_info["total_cost"],
        "estimated_delivery": cost_info["estimated_delivery"],
        "validation": validation_results if validate else None
    }


def split_for_twist(
    sequence: str,
    fragment_id_prefix: str,
    overlap_size: int = 25,
    target_size: int = 1500
) -> List[Dict[str, Any]]:
    """
    Split a long sequence into Twist-compatible fragments.

    Args:
        sequence: Long DNA sequence to split
        fragment_id_prefix: Prefix for fragment IDs
        overlap_size: Size of overlaps between fragments
        target_size: Target fragment size (will be between 300-1800)

    Returns:
        List of fragment dicts ready for ordering

    Example:
        >>> long_gene = "ATGC..." * 1000  # 4000 bp
        >>> frags = split_for_twist(long_gene, "GENE1")
        >>> print(f"Split into {len(frags)} fragments")
    """
    seq = sequence.upper()
    fragments = []

    if len(seq) <= TWIST_MAX_LENGTH:
        # No splitting needed
        return [{
            "id": f"{fragment_id_prefix}_1",
            "sequence": seq,
            "notes": "Full sequence"
        }]

    # Calculate number of fragments needed
    effective_size = target_size - overlap_size
    num_fragments = (len(seq) + effective_size - 1) // effective_size

    for i in range(num_fragments):
        start = i * effective_size
        end = min(start + target_size, len(seq))

        # Adjust last fragment to ensure it's long enough
        if i == num_fragments - 1 and end - start < TWIST_MIN_LENGTH:
            start = max(0, len(seq) - TWIST_MIN_LENGTH)
            end = len(seq)

        frag_seq = seq[start:end]

        fragments.append({
            "id": f"{fragment_id_prefix}_{i+1}",
            "sequence": frag_seq,
            "notes": f"Fragment {i+1} of {num_fragments}, pos {start+1}-{end}"
        })

    return fragments

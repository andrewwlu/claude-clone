"""
Core Gibson assembly logic.

This module contains the core functions for:
- Backbone cutting with restriction enzymes
- Fragment assembly simulation
- Overlap detection and Tm calculation
- Sequence validation
"""

from .backbone import cut_backbone_with_enzymes
from .assembly import assemble_plasmid, rotate_assembly, legacy_assemble
from .overlap import find_overlap, get_Tm, design_overlap
from .validation import validate_sequence, validate_overlap, validate_assembly

__all__ = [
    "cut_backbone_with_enzymes",
    "assemble_plasmid",
    "rotate_assembly",
    "legacy_assemble",
    "find_overlap",
    "get_Tm",
    "design_overlap",
    "validate_sequence",
    "validate_overlap",
    "validate_assembly",
]

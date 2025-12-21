"""
Gibson Assembly Pipeline
========================

A Claude Code-powered pipeline for designing Gibson assembly reactions.

Author: Victoria Tobin (vtobin@caltech.edu)
Affiliation: California Institute of Technology
"""

__version__ = "0.4.1"
__author__ = "Victoria Tobin"
__email__ = "vtobin@caltech.edu"

from .core.assembly import assemble_plasmid, rotate_assembly
from .core.backbone import cut_backbone_with_enzymes
from .core.overlap import find_overlap, get_Tm
from .core.validation import validate_sequence, validate_assembly

__all__ = [
    "assemble_plasmid",
    "rotate_assembly",
    "cut_backbone_with_enzymes",
    "find_overlap",
    "get_Tm",
    "validate_sequence",
    "validate_assembly",
]

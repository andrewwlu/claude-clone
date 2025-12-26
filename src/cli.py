"""
Command-line interface for Gibson assembly pipeline.

This module provides a simple CLI for running Gibson assembly
operations outside of Claude Code context.
"""

import argparse
import sys
import json
from pathlib import Path

from . import __version__
from .core.assembly import assemble_plasmid
from .core.backbone import cut_backbone_with_enzymes
from .core.validation import validate_sequence, validate_assembly
from .io.genbank import read_seq_file, write_gb_file
from .io.sheets_sync import get_cache_metadata, check_cache_freshness
from .io.twist_order import generate_twist_order, validate_twist_constraints


def cmd_validate(args):
    """Validate a sequence file."""
    record = read_seq_file(args.file)
    result = validate_sequence(
        str(record.seq),
        check_synthesis=True,
        check_expression=args.check_expression
    )

    print(f"\nValidation: {record.id}")
    print(f"Length: {len(record.seq)} bp")
    print(f"Status: {'PASS' if result else 'FAIL'}")

    if result.errors:
        print("\nErrors:")
        for err in result.errors:
            print(f"  - {err}")

    if result.warnings:
        print("\nWarnings:")
        for warn in result.warnings:
            print(f"  - {warn}")

    return 0 if result else 1


def cmd_cut(args):
    """Cut a backbone with restriction enzymes."""
    record = read_seq_file(args.backbone)
    linear = cut_backbone_with_enzymes(
        record,
        enzyme1=args.enzyme1,
        enzyme2=args.enzyme2
    )

    print(f"\nCut {record.id} with {args.enzyme1}/{args.enzyme2}")
    print(f"Original: {len(record.seq)} bp (circular)")
    print(f"Linear: {len(linear.seq)} bp")

    if args.output:
        out_path = write_gb_file(linear, Path(args.output).parent, Path(args.output).name)
        print(f"Saved to: {out_path}")

    return 0


def cmd_cache_status(args):
    """Show cache status."""
    metadata = get_cache_metadata()
    freshness = check_cache_freshness(max_age_hours=24)

    print("\nCache Status:")
    print("-" * 50)

    for db_name, info in metadata.items():
        status = "fresh" if freshness.get(db_name) else "stale"
        synced = info.get("synced_at", "never")
        count = info.get("record_count", 0)
        print(f"{db_name:12} | {count:5} records | {synced} | {status}")

    return 0


def cmd_twist_validate(args):
    """Validate sequences for Twist ordering."""
    record = read_seq_file(args.file)
    valid, errors, warnings = validate_twist_constraints(str(record.seq), record.id)

    print(f"\nTwist Validation: {record.id}")
    print(f"Length: {len(record.seq)} bp")
    print(f"Status: {'PASS' if valid else 'FAIL'}")

    if errors:
        print("\nErrors:")
        for err in errors:
            print(f"  - {err}")

    if warnings:
        print("\nWarnings:")
        for warn in warnings:
            print(f"  - {warn}")

    return 0 if valid else 1


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Gibson Assembly Pipeline CLI",
        prog="gibson"
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}"
    )

    subparsers = parser.add_subparsers(dest="command", help="Commands")

    # validate command
    validate_parser = subparsers.add_parser("validate", help="Validate a sequence file")
    validate_parser.add_argument("file", help="Sequence file to validate")
    validate_parser.add_argument("--check-expression", action="store_true",
                                  help="Also check expression requirements")
    validate_parser.set_defaults(func=cmd_validate)

    # cut command
    cut_parser = subparsers.add_parser("cut", help="Cut backbone with enzymes")
    cut_parser.add_argument("backbone", help="Backbone sequence file")
    cut_parser.add_argument("--enzyme1", default="XhoI", help="First enzyme")
    cut_parser.add_argument("--enzyme2", default="HindIII", help="Second enzyme")
    cut_parser.add_argument("--output", "-o", help="Output file path")
    cut_parser.set_defaults(func=cmd_cut)

    # cache-status command
    cache_parser = subparsers.add_parser("cache-status", help="Show cache status")
    cache_parser.set_defaults(func=cmd_cache_status)

    # twist-validate command
    twist_parser = subparsers.add_parser("twist-validate",
                                          help="Validate for Twist ordering")
    twist_parser.add_argument("file", help="Sequence file to validate")
    twist_parser.set_defaults(func=cmd_twist_validate)

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        return 0

    try:
        return args.func(args)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())

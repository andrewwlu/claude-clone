#!/usr/bin/env python3
"""
PostToolUse hook: Validate GenBank files.

This hook runs after Write operations and validates any
GenBank files that were written.

Usage: python validate_genbank.py <file_path>

Exit codes:
- 0: Valid or not a GenBank file
- 2: Invalid GenBank file (blocks with message)
"""

import sys
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from src.io.genbank import validate_genbank


def main():
    if len(sys.argv) < 2:
        return 0

    file_path = Path(sys.argv[1])

    # Only validate .gb and .gbk files
    if file_path.suffix.lower() not in [".gb", ".gbk", ".genbank"]:
        return 0

    if not file_path.exists():
        return 0

    is_valid, errors, warnings = validate_genbank(file_path)

    if warnings:
        for warn in warnings:
            print(f"GenBank warning: {warn}", file=sys.stderr)

    if not is_valid:
        print(f"GenBank validation failed for {file_path.name}:", file=sys.stderr)
        for err in errors:
            print(f"  - {err}", file=sys.stderr)
        return 2  # Block with error message

    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
PreToolUse hook: Validate before writing files.

This hook runs before Write operations and performs
pre-write validation checks.

Usage: python pre_write_validator.py <file_path>

Exit codes:
- 0: OK to write
- 2: Block write with message
"""

import sys
import os
from pathlib import Path


# Files that should never be overwritten
PROTECTED_FILES = [
    ".env",
    ".env.local",
    "credentials.json",
    "service-account.json",
]

# Directories that should not have writes
PROTECTED_DIRS = [
    ".git",
    "node_modules",
    "__pycache__",
]


def main():
    if len(sys.argv) < 2:
        return 0

    file_path = Path(sys.argv[1])

    # Check protected files
    if file_path.name in PROTECTED_FILES:
        print(f"Cannot write to protected file: {file_path.name}", file=sys.stderr)
        print("This file may contain sensitive credentials.", file=sys.stderr)
        return 2

    # Check protected directories
    for protected in PROTECTED_DIRS:
        if protected in file_path.parts:
            print(f"Cannot write to protected directory: {protected}", file=sys.stderr)
            return 2

    # Check for secrets in content (if we had access to content)
    # This would require the hook to receive the content being written

    # Check file extension for certain validations
    if file_path.suffix.lower() in [".gb", ".gbk"]:
        # GenBank files will be validated post-write
        pass

    return 0


if __name__ == "__main__":
    sys.exit(main())

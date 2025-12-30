#!/usr/bin/env python3
"""
PostToolUse hook: Log file operations to audit trail.

This hook runs after Write and Edit operations and logs
the action to the audit trail.

Usage: python audit_log.py <tool_name> <file_path>

Exit codes:
- 0: Success (always succeed to not block operations)
"""

import sys
import os
import json
from datetime import datetime
from pathlib import Path


def main():
    if len(sys.argv) < 3:
        return 0

    tool_name = sys.argv[1]
    file_path = sys.argv[2]

    try:
        audit_file = Path(__file__).parent.parent.parent / "data" / "audit.jsonl"
        audit_file.parent.mkdir(parents=True, exist_ok=True)

        entry = {
            "timestamp": datetime.utcnow().isoformat() + "Z",
            "tool": tool_name,
            "file": file_path,
            "cwd": os.getcwd(),
        }

        with open(audit_file, "a") as f:
            f.write(json.dumps(entry) + "\n")

        return 0

    except Exception as e:
        # Never block on audit failures
        print(f"Audit log warning: {e}", file=sys.stderr)
        return 0


if __name__ == "__main__":
    sys.exit(main())

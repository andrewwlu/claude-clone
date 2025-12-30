#!/usr/bin/env python3
"""
SessionEnd hook: Archive session information.

This hook runs when a Claude Code session ends and archives
relevant session data for later analysis.

Exit codes:
- 0: Success
- Non-zero: Error (but don't block session end)
"""

import sys
import os
import json
from datetime import datetime
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))


def main():
    """Archive session information."""
    try:
        archive_dir = Path(__file__).parent.parent.parent / "data" / "sessions"
        archive_dir.mkdir(parents=True, exist_ok=True)

        # Create session archive entry
        session_info = {
            "ended_at": datetime.utcnow().isoformat() + "Z",
            "working_dir": os.getcwd(),
            "environment": {
                "python_version": sys.version,
                "platform": sys.platform,
            }
        }

        # Append to session log
        log_file = archive_dir / "session_log.jsonl"
        with open(log_file, "a") as f:
            f.write(json.dumps(session_info) + "\n")

        return 0

    except Exception as e:
        # Don't block session end on errors
        print(f"Session archive warning: {e}", file=sys.stderr)
        return 0


if __name__ == "__main__":
    sys.exit(main())

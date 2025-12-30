#!/usr/bin/env python3
"""
SessionStart hook: Check cache freshness.

This hook runs when a Claude Code session starts and checks
if the local data cache is fresh. If stale, it alerts the user.

Exit codes:
- 0: Success (cache is fresh or warning given)
- 2: Block action (should not block session start)
"""

import sys
import os
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from src.io.sheets_sync import check_cache_freshness, get_cache_metadata


def main():
    """Check cache freshness and report status."""
    try:
        freshness = check_cache_freshness(max_age_hours=24)
        metadata = get_cache_metadata()

        stale_dbs = [db for db, fresh in freshness.items() if not fresh]

        if stale_dbs:
            print(f"Gibson Pipeline: Cache stale for: {', '.join(stale_dbs)}")
            print("Consider running /sync-data to refresh.")

        # Check for missing caches
        for db, info in metadata.items():
            if info.get("status") == "not_synced":
                print(f"Gibson Pipeline: {db} database not synced yet.")

        # Always succeed - don't block session start
        return 0

    except Exception as e:
        # Don't block session start on errors
        print(f"Gibson Pipeline: Cache check warning: {e}", file=sys.stderr)
        return 0


if __name__ == "__main__":
    sys.exit(main())

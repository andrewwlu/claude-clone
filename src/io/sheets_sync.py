"""
Google Sheets synchronization via MCP.

This module handles syncing data between Google Sheets databases
and local cache files. It's designed to work with the Google Sheets
MCP server for Claude Code integration.

Note: Direct MCP calls are handled by Claude Code. These functions
support the /sync-data command workflow.
"""

import json
import os
from pathlib import Path
from datetime import datetime
from typing import Dict, Any, List, Optional


# Cache directory
CACHE_DIR = Path(__file__).parent.parent.parent / "data" / "cache"

# TODO: add retry logic for flaky network connections


def get_cache_path(database: str) -> Path:
    """Get the cache file path for a database."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    return CACHE_DIR / f"{database}.json"


def read_cached_data(database: str) -> Optional[Dict[str, Any]]:
    """
    Read cached data from local file.

    Args:
        database: Database name (fragments, plasmids, backbones)

    Returns:
        Cached data as dictionary, or None if not found

    Example:
        >>> fragments = read_cached_data("fragments")
        >>> if fragments:
        ...     print(f"Loaded {len(fragments['records'])} fragments")
    """
    cache_path = get_cache_path(database)

    if not cache_path.exists():
        return None

    with open(cache_path) as f:
        return json.load(f)


def update_cache(database: str, data: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Update local cache with new data.

    Args:
        database: Database name
        data: List of records to cache

    Returns:
        Cache metadata including timestamp and record count

    Example:
        >>> meta = update_cache("fragments", fragments_from_sheets)
        >>> print(f"Cached {meta['record_count']} records")
    """
    cache_path = get_cache_path(database)

    cache_content = {
        "database": database,
        "synced_at": datetime.utcnow().isoformat() + "Z",
        "record_count": len(data),
        "records": data
    }

    with open(cache_path, "w") as f:
        json.dump(cache_content, f, indent=2)

    return {
        "database": database,
        "synced_at": cache_content["synced_at"],
        "record_count": len(data),
        "cache_path": str(cache_path)
    }


def get_cache_metadata() -> Dict[str, Any]:
    """
    Get metadata about all cached databases.

    Returns:
        Dictionary with metadata for each cached database
    """
    metadata = {}

    for db_name in ["fragments", "plasmids", "backbones"]:
        cache_path = get_cache_path(db_name)
        if cache_path.exists():
            with open(cache_path) as f:
                data = json.load(f)
                metadata[db_name] = {
                    "synced_at": data.get("synced_at"),
                    "record_count": data.get("record_count", 0),
                    "cache_path": str(cache_path)
                }
        else:
            metadata[db_name] = {
                "synced_at": None,
                "record_count": 0,
                "cache_path": str(cache_path),
                "status": "not_synced"
            }

    return metadata


def check_cache_freshness(max_age_hours: int = 24) -> Dict[str, bool]:
    """
    Check if cached data is fresh enough.

    Args:
        max_age_hours: Maximum age in hours before data is stale

    Returns:
        Dictionary mapping database names to freshness status
    """
    freshness = {}
    now = datetime.utcnow()

    for db_name in ["fragments", "plasmids", "backbones"]:
        cache = read_cached_data(db_name)
        if cache and "synced_at" in cache:
            synced_at = datetime.fromisoformat(cache["synced_at"].replace("Z", "+00:00"))
            age_hours = (now - synced_at.replace(tzinfo=None)).total_seconds() / 3600
            freshness[db_name] = age_hours < max_age_hours
        else:
            freshness[db_name] = False

    return freshness


def sync_from_sheets(
    mcp_data: Dict[str, List[Dict[str, Any]]],
    databases: Optional[List[str]] = None
) -> Dict[str, Any]:
    """
    Sync data from MCP Google Sheets response to local cache.

    This function is called after Claude Code fetches data from
    Google Sheets via the MCP server.

    Args:
        mcp_data: Data from MCP organized by database name
        databases: List of databases to sync (default: all)

    Returns:
        Sync results with metadata for each database

    Example:
        # After MCP fetch:
        >>> result = sync_from_sheets({
        ...     "fragments": [{"id": "FRG-001", ...}],
        ...     "plasmids": [{"id": "P-001", ...}]
        ... })
    """
    if databases is None:
        databases = list(mcp_data.keys())

    results = {}
    for db_name in databases:
        if db_name in mcp_data:
            results[db_name] = update_cache(db_name, mcp_data[db_name])
        else:
            results[db_name] = {"error": f"No data provided for {db_name}"}

    return results


def lookup_fragment(fragment_id: str) -> Optional[Dict[str, Any]]:
    """
    Look up a fragment by ID from cache.

    Args:
        fragment_id: Fragment ID (e.g., "FRG-0526")

    Returns:
        Fragment record or None if not found
    """
    cache = read_cached_data("fragments")
    if not cache:
        return None

    for record in cache.get("records", []):
        if record.get("id") == fragment_id or record.get("fragment_id") == fragment_id:
            return record

    return None


def lookup_plasmid(plasmid_id: str) -> Optional[Dict[str, Any]]:
    """
    Look up a plasmid by ID from cache.

    Args:
        plasmid_id: Plasmid ID (e.g., "P-0654")

    Returns:
        Plasmid record or None if not found
    """
    cache = read_cached_data("plasmids")
    if not cache:
        return None

    for record in cache.get("records", []):
        if record.get("id") == plasmid_id or record.get("plasmid_id") == plasmid_id:
            return record

    return None


def lookup_backbone(backbone_id: str) -> Optional[Dict[str, Any]]:
    """
    Look up a backbone by ID from cache.

    Args:
        backbone_id: Backbone ID (e.g., "BB-0001")

    Returns:
        Backbone record or None if not found
    """
    cache = read_cached_data("backbones")
    if not cache:
        return None

    for record in cache.get("records", []):
        if record.get("id") == backbone_id or record.get("backbone_id") == backbone_id:
            return record

    return None

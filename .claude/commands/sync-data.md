---
description: Sync fragment, plasmid, and backbone databases from Google Sheets
allowed-tools: Read, Bash, mcp__google-sheets__*
---

# Sync Data from Google Sheets

Refresh the local cache from the Google Sheets databases.

## Process

1. Connect to Google Sheets MCP server
2. Fetch data from each database:
   - FragmentDB (1NkHj2kUESMA6JdYEkofW3iKv2vARPeFw0FQjYFzrfk4)
   - PlasmidDB (15IHOpvQ2-CQXBc9F_oQ5IXDpUMkRuY5aBtp4XXzOqTk)
   - BackboneDB (1H9yiM4tK3D1fCF9g2qcoWlmb0r8UCJSpalzVs1D8NXc)
3. Validate fetched data
4. Write to local cache files in `data/cache/`
5. Log sync operation to audit trail

## Expected Output

After syncing, report:
- Number of fragments synced
- Number of plasmids synced
- Number of backbones synced
- Any validation warnings
- Timestamp of last sync

## Cache Files

Write the following files:
- `data/cache/fragments.json`
- `data/cache/plasmids.json`
- `data/cache/backbones.json`
- `data/cache/sync_metadata.json`

## Error Handling

If MCP server is not available:
1. Check if environment variables are set:
   - GIBSON_SERVICE_ACCOUNT_PATH
   - GIBSON_DRIVE_FOLDER_ID
2. Report which variables are missing
3. Provide setup instructions

## Human-in-the-Loop

Before overwriting existing cache:
1. Show summary of changes (new/modified/deleted records)
2. Ask for confirmation to proceed
3. Only sync after approval

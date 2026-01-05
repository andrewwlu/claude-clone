# Gibson Assembly Automation Pipeline

**Author**: Andrew Lu (andrew.lu.chun@gmail.com)
**Affiliation**: California Institute of Technology
**Version**: 0.4.1

A natural-language interface for designing and simulating Gibson assembly reactions,
integrated with lab fragment/plasmid databases.

## Quick Start

- `/sync-data` - Refresh local cache from Google Sheets
- `/gibson-assemble P-0654` - Run Gibson assembly for a plasmid
- `/council-review P-0654` - Get 10-expert parallel review
- `/validate-sequences FRG-0526` - Validate a fragment sequence

## Key Constraints

### Gibson Assembly Parameters
| Parameter | Min | Max | Optimal |
|-----------|-----|-----|---------|
| Overlap size | 15 bp | 70 bp | 25 bp |
| Overlap Tm | 46°C | 60°C | 52°C |
| Max Tm delta | - | 15°C | - |

### Twist Order Constraints
| Parameter | Min | Max |
|-----------|-----|-----|
| Fragment size | 300 bp | 1800 bp |

## Human-in-the-Loop Protocol

**CRITICAL**: Always show assembly plans before execution:
1. Present proposed actions with all overlaps and Tm values
2. Wait for explicit user confirmation ("yes", "proceed")
3. Execute only after approval
4. Report results with verification steps

## Data Architecture

- **Local cache**: `data/cache/*.json` - Synced from Google Sheets
- **Output files**: `data/plasmids/*.gb`, `data/orders/*.csv`
- **Audit trail**: `data/audit.jsonl`

## Code Reference

- Core logic: `src/core/assembly.py`, `src/core/overlap.py`
- I/O: `src/io/genbank.py`, `src/io/sheets_sync.py`
- Hooks: `src/hooks/`

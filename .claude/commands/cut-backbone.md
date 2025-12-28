---
description: Simulate restriction digest of a backbone
allowed-tools: Read, Grep, Glob
argument-hint: <backbone_id> [enzyme1] [enzyme2]
---

# Cut Backbone: $ARGUMENTS

Simulate restriction enzyme digestion of the specified backbone.

## Input Parsing

Parse $ARGUMENTS:
- First argument: backbone ID (e.g., BB-0001)
- Second argument: enzyme 1 (default: XhoI)
- Third argument: enzyme 2 (default: HindIII)

## Process

1. **Retrieve Backbone Sequence**
   - Look up backbone in cache
   - Load sequence from data/backbones/ if available

2. **Find Cut Sites**
   - Search for enzyme 1 recognition sequence
   - Search for enzyme 2 recognition sequence
   - Report positions of all sites

3. **Simulate Digest**
   - Calculate resulting fragment sizes
   - Identify linearized backbone
   - Report overhangs generated

4. **Validate for Gibson**
   - Check that exactly one site per enzyme (or report alternatives)
   - Ensure cut sites are in appropriate locations
   - Verify insert region is correct

## Enzyme Reference

From config/gibson_config.json:
```
XhoI: C/TCGAG (5' overhang: TCGA)
HindIII: A/AGCTT (5' overhang: AGCT)
EcoRI: G/AATTC (5' overhang: AATT)
BamHI: G/GATCC (5' overhang: GATC)
```

## Output Format

```
## Backbone Digest: [backbone_id]

### Enzyme Sites
| Enzyme | Recognition | Positions | Count |
|--------|-------------|-----------|-------|
| XhoI   | CTCGAG      | 1234      | 1     |
| HindIII| AAGCTT      | 5678      | 1     |

### Digest Products
| Fragment | Size (bp) | Description |
|----------|-----------|-------------|
| Linear backbone | 4444 | Insert region removed |
| Insert dropout | 1234 | Will be replaced |

### Overhangs
- 5' overhang: TCGA (from XhoI)
- 3' overhang: AGCT (from HindIII)

### Gibson Compatibility
[Assessment of compatibility with Gibson assembly]
```

## Error Handling

- If backbone not found: suggest similar IDs
- If enzyme not recognized: show available enzymes
- If multiple cut sites: show all and ask which to use
- If no cut sites: report and suggest alternatives

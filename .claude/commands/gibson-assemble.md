---
description: Run Gibson assembly simulation for a plasmid
allowed-tools: Read, Write, Grep, Glob, Task, mcp__google-sheets__*
argument-hint: <plasmid_id>
---

# Gibson Assembly: $ARGUMENTS

Simulate Gibson assembly for the specified plasmid.

## Process

1. **Retrieve Plasmid Design**
   - Look up $ARGUMENTS in PlasmidDB cache
   - Get backbone ID and fragment IDs
   - Fetch all sequences

2. **Validate Components**
   - Check all fragments exist in cache
   - Verify backbone sequence available
   - Validate fragment sizes (300-1800 bp for Twist)

3. **Design Assembly**
   - Use assembly-planner agent to design overlaps
   - Calculate Tm for all junctions
   - Validate against constraints:
     - Overlap: 15-70 bp (optimal 25 bp)
     - Tm: 46-60°C (optimal 52°C)
     - Max Tm delta: 15°C

4. **Human Approval**
   - Present assembly plan with all details
   - Show junction table with Tm values
   - Wait for explicit approval ("yes", "proceed")

5. **Execute Assembly (after approval)**
   - Simulate assembly in silico
   - Cut backbone with specified enzymes
   - Join fragments with overlaps
   - Rotate to canonical start position

6. **Generate Output**
   - Create GenBank file in data/plasmids/
   - Update PlasmidDB status
   - Log operation to audit trail

## Output Format

```
## Assembly Plan: $ARGUMENTS

### Components
- Backbone: [ID] ([size] bp)
- Fragments: [list with sizes]

### Junction Analysis
| Junction | Overlap (bp) | Tm (°C) | Status |
|----------|--------------|---------|--------|
| BB-F1    | 25           | 52.3    | ✓      |
| F1-F2    | 25           | 51.8    | ✓      |
| F2-BB    | 25           | 53.1    | ✓      |

### Estimated Success: [%]

Proceed with assembly? [yes/no]
```

## Error Handling

- If plasmid not found: suggest similar IDs
- If fragments missing: list missing IDs and suggest /sync-data
- If constraints violated: show specific issues and suggest redesign

---
description: Get comprehensive 10-expert parallel review of an assembly design
allowed-tools: Read, Grep, Glob, Task
argument-hint: <plasmid_id>
---

# Assembly Council Review: $ARGUMENTS

Launch all 10 expert agents in parallel to review the assembly design.

## Expert Panel

Launch these agents IN PARALLEL using the Task tool:

1. **thermodynamic-expert** (haiku) - Overlap Tm analysis
2. **expression-expert** (haiku) - Protein expression potential
3. **primer-designer** (haiku) - Verification primers, colony PCR
4. **manufacturing-expert** (haiku) - Synthesis feasibility
5. **annotation-expert** (haiku) - GenBank compliance
6. **strategy-expert** (sonnet) - Alternative approaches
7. **cost-optimizer** (haiku) - Cost analysis
8. **timeline-expert** (haiku) - Project timeline
9. **host-compatibility-expert** (sonnet) - Expression host analysis
10. **troubleshooting-expert** (sonnet) - Failure diagnosis, fixes

## Execution

```
Use Task tool to spawn all 10 agents in parallel.
Each agent analyzes the same plasmid design independently.
Collect all results and synthesize into unified report.
```

## Report Synthesis

After all experts return, synthesize their reports:

```
## Council Review: $ARGUMENTS

### Summary
| Expert | Verdict | Key Concern |
|--------|---------|-------------|
| Thermodynamic | PASS/WARN/FAIL | [summary] |
| Expression | PASS/WARN/FAIL | [summary] |
| ... | ... | ... |

### Critical Issues (FAIL verdicts)
[List any FAIL verdicts with details]

### Warnings
[List all warnings from experts]

### Recommendations
[Consolidated recommendations]

### Overall Verdict: [APPROVED/NEEDS_REVISION/REJECTED]

### Required Actions
[If not APPROVED, list what must be addressed]
```

## Decision Rules

- **APPROVED**: All experts PASS or WARN (no FAIL)
- **NEEDS_REVISION**: 1-2 FAIL verdicts with clear fixes
- **REJECTED**: 3+ FAILs or fundamental design problems

## Human-in-the-Loop

After presenting council results:
1. Highlight any FAIL verdicts prominently
2. Ask user how to proceed
3. Offer to address specific concerns
4. Do not proceed with assembly until concerns resolved

---
description: Interactively design a new gene fragment with overlaps
allowed-tools: Read, Write, Grep, Glob
argument-hint: <gene_name> [backbone_id]
---

# Design Fragment: $ARGUMENTS

Interactive workflow for designing a new gene fragment for Gibson assembly.

## Input Parsing

- First argument: gene name or sequence source
- Second argument (optional): target backbone ID

## Workflow

### Step 1: Sequence Input
Ask user for:
- Gene sequence (paste or file path)
- Source organism
- Codon optimization target (if needed)

### Step 2: Fragment Boundaries
Analyze sequence and propose:
- Fragment count (based on size limits 300-1800 bp)
- Split points (at natural boundaries if possible)
- Overlap regions

### Step 3: Overlap Design
For each junction:
- Design 25 bp overlap
- Calculate Tm
- Optimize GC content
- Check for secondary structure

### Step 4: Validation
Run validation checks:
- Synthesis compatibility
- Assembly compatibility
- Expression readiness

### Step 5: Human Review
Present complete design:
- All fragments with sequences
- All overlaps with Tm
- Estimated costs
- Potential issues

### Step 6: Output (after approval)
- Create fragment entries for database
- Generate sequences with overlaps
- Save to data/fragments/

## Design Constraints

```
Fragment size: 300-1800 bp (Twist)
Overlap size: 25 bp (optimal)
Overlap Tm: 52°C (target)
Max fragments per assembly: 6
```

## Output Format

```
## Fragment Design: $ARGUMENTS

### Overview
- Gene: [name]
- Total length: [bp]
- Fragments: [count]
- Backbone: [id]

### Fragment Details

#### Fragment 1: [ID]
- Size: [bp]
- GC: [%]
- 5' overlap: [sequence] (Tm: [°C])
- 3' overlap: [sequence] (Tm: [°C])
- Sequence:
```
[sequence in FASTA format]
```

#### Fragment 2: [ID]
...

### Assembly Preview
[ASCII diagram of assembly]

### Estimated Costs
- Synthesis: $[amount]
- Total: $[amount]

Save fragment design? [yes/no]
```

## Codon Optimization

If codon optimization requested:
1. Ask for target organism (E. coli, Human, Yeast)
2. Apply organism-specific codon preferences
3. Maintain CAI > 0.8
4. Remove rare codon clusters
5. Report optimization changes

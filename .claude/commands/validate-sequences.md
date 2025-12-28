---
description: Validate sequences for synthesis and assembly compatibility
allowed-tools: Read, Grep, Glob
argument-hint: <sequence_id_or_file>
---

# Validate Sequences: $ARGUMENTS

Run comprehensive validation on the specified sequence(s).

## Input Handling

- If $ARGUMENTS is a fragment ID (e.g., FRG-0526): look up in cache
- If $ARGUMENTS is a file path: read sequence from file
- If $ARGUMENTS is "all": validate all cached fragments

## Validation Checks

### Basic Validation
- [ ] Valid nucleotides only (ATGC)
- [ ] No ambiguous bases (N, R, Y, etc.)
- [ ] Length within synthesis limits

### Synthesis Compatibility
- [ ] GC content 25-75%
- [ ] No homopolymer runs >6 bp
- [ ] No problematic repeats >8 bp
- [ ] Complexity score acceptable

### Assembly Compatibility
- [ ] Overlaps present and correct
- [ ] Overlap Tm in range (46-60°C)
- [ ] No conflicting restriction sites
- [ ] Reading frame maintained (if CDS)

### Expression Readiness (if CDS)
- [ ] Start codon present
- [ ] Stop codon present
- [ ] Length divisible by 3
- [ ] No internal stop codons

## Output Format

```
## Validation Report: $ARGUMENTS

### Summary
- Sequences validated: [count]
- Passed: [count]
- Warnings: [count]
- Failed: [count]

### Results

#### [Sequence ID 1]
- Length: [bp]
- GC content: [%]
- Status: [PASS/WARN/FAIL]
- Issues: [list]

#### [Sequence ID 2]
...

### Action Required
[List sequences needing attention]
```

## Batch Mode

When validating multiple sequences:
1. Process all sequences
2. Sort results by status (FAIL first, then WARN, then PASS)
3. Summarize common issues
4. Suggest batch fixes if patterns detected

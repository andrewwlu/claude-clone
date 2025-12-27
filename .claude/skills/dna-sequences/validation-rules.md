# Sequence Validation Rules

Comprehensive validation rules for DNA sequences in the Gibson assembly pipeline.

## Pre-Synthesis Validation

### Required Checks

```python
validation_rules = {
    "valid_characters": {
        "pattern": r"^[ATGC]+$",
        "severity": "ERROR",
        "message": "Sequence contains invalid characters"
    },
    "min_length": {
        "value": 50,
        "severity": "ERROR",
        "message": "Sequence too short for synthesis"
    },
    "max_length": {
        "value": 5000,
        "severity": "WARNING",
        "message": "Sequence may need splitting"
    },
    "gc_content": {
        "min": 25,
        "max": 75,
        "severity": "WARNING",
        "message": "Extreme GC content may cause synthesis issues"
    }
}
```

### Homopolymer Detection

```python
def find_homopolymers(seq, max_length=6):
    """Find homopolymer runs exceeding max_length."""
    issues = []
    for nucleotide in "ATGC":
        pattern = nucleotide * (max_length + 1)
        pos = 0
        while True:
            pos = seq.find(pattern, pos)
            if pos == -1:
                break
            # Find actual length
            end = pos
            while end < len(seq) and seq[end] == nucleotide:
                end += 1
            issues.append({
                "type": "homopolymer",
                "nucleotide": nucleotide,
                "position": pos,
                "length": end - pos
            })
            pos = end
    return issues
```

### Repeat Detection

```python
def find_repeats(seq, min_length=8):
    """Find direct repeats that may cause issues."""
    repeats = []
    for length in range(min_length, len(seq) // 2):
        seen = {}
        for i in range(len(seq) - length + 1):
            subseq = seq[i:i+length]
            if subseq in seen:
                repeats.append({
                    "sequence": subseq,
                    "positions": [seen[subseq], i],
                    "length": length
                })
            else:
                seen[subseq] = i
    return repeats
```

## Assembly Validation

### Overlap Validation

```python
def validate_overlap(overlap_seq):
    """Validate overlap sequence for Gibson assembly."""
    errors = []
    warnings = []

    # Length check
    if len(overlap_seq) < 15:
        errors.append("Overlap too short (<15 bp)")
    elif len(overlap_seq) > 70:
        errors.append("Overlap too long (>70 bp)")
    elif len(overlap_seq) < 20:
        warnings.append("Overlap on short side (consider 25+ bp)")

    # GC content
    gc = (overlap_seq.count('G') + overlap_seq.count('C')) / len(overlap_seq)
    if gc < 0.25:
        errors.append(f"GC content too low ({gc*100:.1f}%)")
    elif gc > 0.75:
        errors.append(f"GC content too high ({gc*100:.1f}%)")
    elif gc < 0.40 or gc > 0.60:
        warnings.append(f"GC content suboptimal ({gc*100:.1f}%)")

    # Homopolymers
    for nuc in "ATGC":
        if nuc * 5 in overlap_seq:
            warnings.append(f"Homopolymer run of {nuc}")

    return errors, warnings
```

### Junction Validation

```python
def validate_junction(frag1_end, frag2_start, expected_overlap):
    """Validate that fragments have matching overlaps."""
    # Check that 3' of frag1 matches expected overlap
    if frag1_end[-len(expected_overlap):] != expected_overlap:
        return False, "Fragment 1 end doesn't match expected overlap"

    # Check that 5' of frag2 matches expected overlap
    if frag2_start[:len(expected_overlap)] != expected_overlap:
        return False, "Fragment 2 start doesn't match expected overlap"

    return True, "Junction valid"
```

## Post-Assembly Validation

### Circular Sequence Check

```python
def validate_circular(seq, expected_length, tolerance=10):
    """Validate circularized assembly."""
    if abs(len(seq) - expected_length) > tolerance:
        return False, f"Length mismatch: {len(seq)} vs {expected_length}"
    return True, "Length correct"
```

### Feature Preservation

```python
def validate_features(seq, expected_features):
    """Verify all expected features are present."""
    missing = []
    for feature in expected_features:
        if feature["sequence"] not in seq:
            missing.append(feature["name"])
    if missing:
        return False, f"Missing features: {', '.join(missing)}"
    return True, "All features present"
```

---
name: dna-sequences
description: Validate and analyze DNA sequences, check for synthesis issues, codon optimization. Use when working with raw sequences or preparing for gene synthesis.
version: 0.4.1
allowed-tools: Read, Grep, Glob
---

# DNA Sequence Analysis

Expert knowledge for validating and analyzing DNA sequences.

## Sequence Validation

### Valid Nucleotides

Standard DNA alphabet:
- **A** - Adenine
- **T** - Thymine
- **G** - Guanine
- **C** - Cytosine

IUPAC ambiguity codes (avoid in synthesis):
- R (A/G), Y (C/T), M (A/C), K (G/T)
- S (G/C), W (A/T), B (C/G/T), D (A/G/T)
- H (A/C/T), V (A/C/G), N (any)

### GC Content

```python
def gc_content(seq):
    """Calculate GC content as percentage."""
    gc = seq.count('G') + seq.count('C')
    return (gc / len(seq)) * 100
```

| Range | Status | Notes |
|-------|--------|-------|
| 40-60% | Optimal | Best for synthesis and PCR |
| 30-40% | Acceptable | May need optimization |
| 60-70% | Acceptable | Watch for secondary structure |
| <30% | Problematic | AT-rich, low Tm |
| >70% | Problematic | GC-rich, synthesis issues |

## Synthesis Constraints

### Homopolymer Limits

| Nucleotide | Max Run | Notes |
|------------|---------|-------|
| A | 6 | Longer runs cause polymerase slippage |
| T | 6 | Especially problematic in promoters |
| G | 4 | G-quadruplex formation risk |
| C | 6 | Less problematic than G |

### Repeat Sequences

- **Direct repeats**: Max 8 bp identical sequences
- **Inverted repeats**: Max 8 bp (hairpin formation)
- **Tandem repeats**: Avoid >3 copies of any motif

### Problem Sequences

```python
problem_patterns = {
    "GGGG": "G-quadruplex risk",
    "TATA": "May act as promoter element",
    "AATAAA": "Polyadenylation signal",
    "GGTACC": "KpnI site - may interfere",
}
```

## Codon Optimization

### Codon Adaptation Index (CAI)

```python
def calculate_cai(seq, codon_table):
    """Calculate CAI for expression optimization."""
    codons = [seq[i:i+3] for i in range(0, len(seq)-2, 3)]
    weights = [codon_table.get(c, {}).get('weight', 0) for c in codons]
    return geometric_mean(weights)
```

### Common Host Preferences

| Host | Preferred | Avoid |
|------|-----------|-------|
| E. coli | GCG (Ala), CTG (Leu) | AGG (Arg), AGA (Arg) |
| Human | GCC (Ala), CTG (Leu) | TCG (Ser), CCG (Pro) |
| Yeast | GCT (Ala), TTG (Leu) | CGC (Arg), CCC (Pro) |

## Reading Frame Analysis

### Start Codons
- **ATG** - Standard start (Met)
- **GTG** - Alternative in prokaryotes (Val→Met)
- **TTG** - Rare alternative in prokaryotes

### Stop Codons
- **TAA** - Ochre (most efficient in E. coli)
- **TAG** - Amber
- **TGA** - Opal (least efficient in E. coli)

### Frame Validation

```python
def validate_frame(seq):
    """Check that CDS is valid."""
    if len(seq) % 3 != 0:
        return False, "Length not divisible by 3"
    if seq[:3] != "ATG":
        return False, "Missing start codon"
    if seq[-3:] not in ["TAA", "TAG", "TGA"]:
        return False, "Missing stop codon"
    return True, "Valid"
```

## See Also

- @validation-rules.md - Detailed validation rules

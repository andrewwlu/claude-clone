---
name: gibson-assembly
description: Design Gibson assembly reactions, calculate overlaps and Tm, validate junctions. Use when designing plasmid constructs or troubleshooting assemblies.
version: 0.4.1
allowed-tools: Read, Grep, Glob, Edit, Write
---

# Gibson Assembly Design

Expert knowledge for designing and executing Gibson assembly reactions.

## Core Principles

Gibson assembly uses a 3-enzyme cocktail:
1. **T5 exonuclease** - Creates 3' overhangs by chewing back 5' ends
2. **Phusion polymerase** - Fills in gaps
3. **Taq ligase** - Seals nicks

## Overlap Design Rules

### Size Requirements
| Parameter | Min | Max | Optimal |
|-----------|-----|-----|---------|
| Overlap | 15 bp | 70 bp | 25 bp |

### Melting Temperature
| Parameter | Min | Max | Optimal |
|-----------|-----|-----|---------|
| Overlap Tm | 46°C | 60°C | 52°C |
| Max delta | - | 15°C | <5°C |

### GC Content
- Optimal: 40-60%
- Avoid: <25% or >75%
- 3' end: Should end in G or C (GC clamp)

## Tm Calculation

Use nearest-neighbor method (not simple GC):

```python
from Bio.SeqUtils import MeltingTemp as mt

def get_Tm(seq):
    """Calculate Tm using nearest-neighbor method."""
    return mt.Tm_NN(seq, nn_table=mt.DNA_NN4)
```

**Note**: We switched from `Tm_GC` to `Tm_NN` in v0.4.1 for accuracy.

## Fragment Limits

- **Maximum fragments**: 6 (efficiency drops significantly beyond)
- **Optimal fragments**: 2-4
- **Minimum fragment size**: 200 bp (for handling)

## Common Failure Modes

1. **Low Tm overlaps** - Assembly fails, no colonies
2. **Unbalanced Tm** - Preferential assembly of some junctions
3. **Secondary structure** - Overlaps form hairpins
4. **Homopolymers** - T5 exonuclease stalls

## Protocol Reference

Standard Gibson Assembly (NEB):
1. Set up reaction on ice
2. Add fragments in equimolar ratio
3. Incubate at 50°C for 15-60 min
4. Transform 2 µL into competent cells

## See Also

- @overlap-rules.md - Detailed overlap design
- @restriction-enzymes.md - Enzyme cut sites

# Overlap Design Rules

Detailed rules for designing Gibson assembly overlaps.

## Optimal Overlap Sequence

### Sequence Composition

1. **GC Content**: 40-60%
   - Too low (<40%): Weak binding, low Tm
   - Too high (>60%): Secondary structure risk

2. **GC Clamp**: End overlaps with G or C
   - Stabilizes 3' end binding
   - Improves polymerase extension

3. **Avoid Homopolymers**
   - Max consecutive identical bases: 4
   - T5 exonuclease can stall on long runs

4. **Avoid Repeats**
   - No direct repeats >6 bp
   - No inverted repeats >8 bp (hairpin risk)

## Overlap Positioning

### At CDS Junctions

```
Fragment 1:  ...ATG-GAC-CAT-GGC-[overlap]
Fragment 2:            [overlap]-AAA-CGT-TAA...
```

- Maintain reading frame across junction
- Avoid splitting codons if possible
- Use silent mutations to optimize overlap

### At Non-Coding Junctions

- More flexibility in sequence choice
- Prioritize Tm optimization
- Can use synthetic bridging sequences

## Tm Balancing

All overlaps in an assembly should have similar Tm:

```
Junction 1 Tm: 52°C
Junction 2 Tm: 54°C  ← Ideal (within 5°C)
Junction 3 Tm: 51°C
```

vs.

```
Junction 1 Tm: 48°C
Junction 2 Tm: 58°C  ← Problematic (10°C delta)
Junction 3 Tm: 52°C
```

## Calculation Methods

### Nearest-Neighbor (Preferred)

```python
from Bio.SeqUtils import MeltingTemp as mt

tm = mt.Tm_NN(
    seq,
    nn_table=mt.DNA_NN4,
    dnac1=250,  # primer concentration nM
    dnac2=250,
    Na=50,      # sodium concentration mM
    Mg=0        # magnesium concentration mM
)
```

### Quick Estimate (4+2 Rule)

For rough estimates only:
```
Tm ≈ 4(G+C) + 2(A+T)
```

**Warning**: This is inaccurate for Gibson overlaps. Use nearest-neighbor.

## Troubleshooting

| Problem | Cause | Solution |
|---------|-------|----------|
| No colonies | Low Tm | Increase overlap length |
| Wrong assembly | Tm imbalance | Balance all Tm values |
| Low efficiency | Secondary structure | Redesign overlap sequence |
| Multiple products | Repeat sequences | Use unique overlaps |

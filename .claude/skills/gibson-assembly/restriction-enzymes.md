# Restriction Enzymes Reference

Common restriction enzymes used in Gibson assembly workflows.

## Standard Type II Enzymes

### For Backbone Linearization

| Enzyme | Recognition | Cut Pattern | Overhang |
|--------|-------------|-------------|----------|
| XhoI | CTCGAG | C/TCGAG | 5' TCGA |
| HindIII | AAGCTT | A/AGCTT | 5' AGCT |
| EcoRI | GAATTC | G/AATTC | 5' AATT |
| BamHI | GGATCC | G/GATCC | 5' GATC |
| NcoI | CCATGG | C/CATGG | 5' CATG |
| NdeI | CATATG | CA/TATG | 5' TA |
| XbaI | TCTAGA | T/CTAGA | 5' CTAG |
| SpeI | ACTAGT | A/CTAGT | 5' CTAG |
| PstI | CTGCAG | CTGCA/G | 3' TGCA |
| SalI | GTCGAC | G/TCGAC | 5' TCGA |

### Blunt Cutters

| Enzyme | Recognition | Cut Pattern |
|--------|-------------|-------------|
| EcoRV | GATATC | GAT/ATC |
| SmaI | CCCGGG | CCC/GGG |
| HpaI | GTTAAC | GTT/AAC |

## PCR Primer Sites (Custom)

For amplification-based cloning:

```python
pcr_sites = {
    "PCR1": {
        "recognition": "ACTTAAGCTTCGCCACCATG",
        "cut": "/ACTTAAGCTTCGCCACCATG",
        "note": "Kozak-HindIII fusion"
    },
    "PCR2": {
        "recognition": "GCATGGACGAGCTGTACAAG",
        "cut": "GCATGGACGAGCTGTACAAG/",
        "note": "C-terminal fusion linker"
    }
}
```

## Enzyme Compatibility

### Buffer Compatibility Matrix

Common double digests that work in single buffer:

| Enzyme 1 | Enzyme 2 | Buffer | Efficiency |
|----------|----------|--------|------------|
| XhoI | HindIII | CutSmart | 100%/100% |
| EcoRI | BamHI | CutSmart | 100%/100% |
| NcoI | XhoI | CutSmart | 100%/100% |
| EcoRI | HindIII | CutSmart | 100%/100% |

### Heat Inactivation

| Enzyme | Temperature | Time |
|--------|-------------|------|
| XhoI | 65°C | 20 min |
| HindIII | 80°C | 20 min |
| EcoRI | 65°C | 20 min |
| BamHI | 65°C | 20 min |

## Star Activity Warning

These enzymes show star activity (non-specific cutting) under certain conditions:

- **EcoRI** - High glycerol, low ionic strength
- **BamHI** - Extended incubation, high enzyme
- **HindIII** - High glycerol

**Prevention**: Use NEB CutSmart buffer, don't over-digest.

## Configuration Reference

See `config/gibson_config.json` for enzyme definitions used in this pipeline:

```json
{
  "restriction_enzymes": {
    "XhoI": {"recognition": "CTCGAG", "cut": "C/TCGAG"},
    "HindIII": {"recognition": "AAGCTT", "cut": "A/AGCTT"}
  }
}
```

---
description: Generate Twist Bioscience order CSV for gene fragments
allowed-tools: Read, Write, Grep, Glob
argument-hint: <plasmid_id_or_fragment_list>
---

# Generate Twist Order: $ARGUMENTS

Create a CSV file formatted for Twist Bioscience gene fragment ordering.

## Input Handling

- If $ARGUMENTS is a plasmid ID: get all fragments for that plasmid
- If $ARGUMENTS is comma-separated fragment IDs: use those fragments
- If $ARGUMENTS is "pending": get all fragments with status "not ordered"

## Twist Constraints

Validate all fragments against Twist limits:
- Minimum size: 300 bp
- Maximum size: 1800 bp
- GC content: 25-75% (preferred)
- No homopolymers >6 bp

## CSV Format

Twist Gene Fragments format:
```csv
Name,Sequence,Notes
FRG-0526,ATGCGATCG...,MAP2K1 fragment 1
FRG-0527,ATGCGATCG...,MAP2K1 fragment 2
```

## Process

1. **Retrieve Fragments**
   - Look up all requested fragments
   - Validate against Twist constraints

2. **Pre-Order Validation**
   - Check all fragments meet size requirements
   - Flag any complexity issues
   - Calculate estimated cost

3. **Human Approval**
   - Show order summary with costs
   - List any warnings
   - Wait for confirmation

4. **Generate Order (after approval)**
   - Create CSV in data/orders/
   - Update fragment status in database
   - Log order to audit trail

## Output Format

```
## Twist Order: [order_id]

### Order Summary
| Fragment | Size (bp) | GC% | Status | Est. Cost |
|----------|-----------|-----|--------|-----------|
| FRG-0526 | 1200      | 52% | ✓      | $108      |
| FRG-0527 | 1500      | 48% | ✓      | $135      |

**Total Fragments:** 2
**Total Base Pairs:** 2700 bp
**Estimated Cost:** $243 (at $0.09/bp)
**Estimated Delivery:** 13-17 business days

### Warnings
[Any issues to note]

Proceed with order generation? [yes/no]
```

## Output File

Write to: `data/orders/twist_order_[YYYYMMDD]_[plasmid_id].csv`

## Cost Estimation

- Standard pricing: $0.09/bp
- Express pricing: $0.15/bp (5-day)
- Minimum order: $0 (no minimum)

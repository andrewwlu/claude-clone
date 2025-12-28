---
name: genbank-files
description: Read, write, and manipulate GenBank files. Use when working with .gb/.gbk files, adding annotations, or exporting sequences.
version: 0.4.1
allowed-tools: Read, Write, Glob, Edit
---

# GenBank File Handling

Expert knowledge for working with GenBank format sequence files.

## File Structure

```
LOCUS       name        length bp    mol    topology    division  date
DEFINITION  description
ACCESSION   accession_number
VERSION     version
KEYWORDS    keywords
SOURCE      source_organism
  ORGANISM  organism_name
            classification
REFERENCE   reference_number
  AUTHORS   author_names
  TITLE     title
  JOURNAL   journal
FEATURES             Location/Qualifiers
     feature_key     location
                     /qualifier="value"
ORIGIN
        1 sequence...
//
```

## Common Feature Keys

### Structural Features
| Key | Description | Required Qualifiers |
|-----|-------------|---------------------|
| `source` | Entire sequence | /organism, /mol_type |
| `gene` | Gene region | /gene |
| `CDS` | Coding sequence | /gene, /product, /translation |
| `mRNA` | Messenger RNA | /gene, /product |
| `exon` | Exon region | /gene, /number |
| `intron` | Intron region | /gene, /number |

### Regulatory Features
| Key | Description | Common Qualifiers |
|-----|-------------|-------------------|
| `promoter` | Promoter region | /note |
| `terminator` | Terminator region | /note |
| `enhancer` | Enhancer region | /note |
| `regulatory` | Generic regulatory | /regulatory_class |

### Other Features
| Key | Description | Common Qualifiers |
|-----|-------------|-------------------|
| `rep_origin` | Replication origin | /note |
| `misc_feature` | Miscellaneous | /note |
| `primer_bind` | Primer binding site | /note |

## Location Syntax

### Simple Locations
```
100..200        # From 100 to 200 (inclusive)
<100..200       # 5' partial
100..>200       # 3' partial
complement(100..200)  # Minus strand
```

### Complex Locations
```
join(1..100,200..300)           # Split feature
complement(join(1..100,200..300))  # Minus strand split
order(1..100,200..300)          # Uncertain order
```

### Circular Features
```
join(4000..5000,1..500)         # Crosses origin
```

## Reading GenBank Files

```python
from Bio import SeqIO

def read_genbank(filepath):
    """Read GenBank file and return SeqRecord."""
    record = SeqIO.read(filepath, "genbank")
    return record

def extract_cds(record):
    """Extract all CDS features."""
    cds_features = []
    for feature in record.features:
        if feature.type == "CDS":
            cds_features.append({
                "gene": feature.qualifiers.get("gene", ["unknown"])[0],
                "product": feature.qualifiers.get("product", ["unknown"])[0],
                "location": str(feature.location),
                "translation": feature.qualifiers.get("translation", [""])[0]
            })
    return cds_features
```

## Writing GenBank Files

```python
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

def create_genbank(seq, features, metadata):
    """Create GenBank file from sequence and features."""
    record = SeqRecord(
        Seq(seq),
        id=metadata["id"],
        name=metadata["name"],
        description=metadata["description"],
        annotations={
            "molecule_type": "DNA",
            "topology": metadata.get("topology", "circular"),
            "organism": metadata.get("organism", "synthetic construct")
        }
    )

    # Add features
    for f in features:
        feature = SeqFeature(
            FeatureLocation(f["start"], f["end"], strand=f.get("strand", 1)),
            type=f["type"],
            qualifiers=f.get("qualifiers", {})
        )
        record.features.append(feature)

    return record

def write_genbank(record, filepath):
    """Write SeqRecord to GenBank file."""
    SeqIO.write(record, filepath, "genbank")
```

## Validation

```python
def validate_genbank(filepath):
    """Validate GenBank file format."""
    try:
        record = SeqIO.read(filepath, "genbank")

        # Check required elements
        errors = []
        if not record.id:
            errors.append("Missing LOCUS ID")
        if not record.seq:
            errors.append("Missing sequence")
        if not any(f.type == "source" for f in record.features):
            errors.append("Missing source feature")

        return len(errors) == 0, errors
    except Exception as e:
        return False, [str(e)]
```

## Output Directory

GenBank files should be written to: `data/plasmids/`

Naming convention: `[PLASMID_ID].gb` (e.g., `P-0654.gb`)

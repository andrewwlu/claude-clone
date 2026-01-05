# Gibson Assembly Pipeline

Design plasmids with natural language. Automates Gibson assembly design, validation, and ordering for gene synthesis.

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## Example

**You say:**
```
I need to clone human BRCA1 exons 11-15 into pET28a for bacterial expression.
The insert is 2.4kb. Design the assembly.
```

**Claude responds:**

```
Analyzing BRCA1 exons 11-15 (2,847 bp) for pET28a backbone...

Fragment design (Twist synthesis):
┌─────────────────────────────────────────────────────────────────────┐
│  pET28a          FRG-1           FRG-2           pET28a            │
│  (backbone) ───▶ (1,423 bp) ───▶ (1,449 bp) ───▶ (backbone)        │
│             25bp          25bp            25bp                      │
│            overlap       overlap         overlap                    │
│            Tm=54°C       Tm=52°C         Tm=55°C                   │
└─────────────────────────────────────────────────────────────────────┘

Launching Assembly Council (10 experts in parallel)...

┌────────────────────┬────────┬─────────────────────────────────────┐
│ Expert             │ Verdict│ Notes                               │
├────────────────────┼────────┼─────────────────────────────────────┤
│ Thermodynamic      │ PASS   │ All overlaps 52-55°C, ΔTm = 3°C     │
│ Expression         │ PASS   │ CAI 0.76, no rare codon clusters    │
│ Manufacturing      │ PASS   │ Both fragments within Twist limits  │
│ Primer Designer    │ PASS   │ 4 primers designed for verification │
│ Host Compatibility │ PASS   │ BL21(DE3) recommended               │
│ Troubleshooting    │ WARN   │ GC-rich region at 1.2kb, monitor    │
│ Strategy           │ PASS   │ 2-fragment approach optimal         │
│ Annotation         │ PASS   │ All features properly annotated     │
│ Cost Optimizer     │ PASS   │ $258 total ($0.09/bp × 2,872 bp)    │
│ Timeline           │ PASS   │ 15 days synthesis + 3 days cloning  │
└────────────────────┴────────┴─────────────────────────────────────┘

Council Verdict: APPROVED

Shall I generate the GenBank file and Twist order?
```

**After approval:**

```
Created:
  → data/plasmids/P-0892.gb          (annotated GenBank file)
  → data/orders/twist_2025-01-04.csv (ready to upload to Twist)

Verification primers:
  P0892_colF  5'-GCTTAGATCGATCGATCGAT-3'  Tm=58°C  ┐
  P0892_colR  5'-ATCGATCGATCGATCGATCG-3'  Tm=57°C  ┘ → 1,247 bp (spans junction 1)

  P0892_seq1  5'-TAATACGACTCACTATAGGG-3'  (T7 promoter, standard)
  P0892_seq2  5'-GATCGATCGATCGATCGATC-3'  Tm=58°C  (insert coverage)
  P0892_seq3  5'-CGATCGATCGATCGATCGAT-3'  Tm=57°C  (insert coverage)
```

---

## Quick Start

```bash
pip install -r requirements.txt
./scripts/setup_mcp.sh
claude
```

---

## Data Layer

The pipeline connects to three Google Sheets databases via MCP:

| Database | Contents | Example |
|----------|----------|---------|
| Fragments | Synthesized DNA fragments | FRG-0526: BRCA1 exon 11 |
| Plasmids | Completed constructs | P-0654: pET28a-BRCA1 |
| Backbones | Vector sequences | pET28a, pUC19, pBAD |

Data flows through three tiers:

```
Google Sheets (source of truth)
       ↓ MCP sync
Local cache (data/cache/*.json)
       ↓ Python modules
Assembly engine (src/core/)
```

The `/sync-data` command refreshes the local cache. Session hooks check cache freshness on startup and warn if data is stale.

---

## Architecture

```
Claude Code
├── Skills
│   ├── gibson-assembly      Overlap design, Tm calculation
│   ├── dna-sequences        Validation, codon optimization
│   └── genbank-files        File I/O, annotations
│
├── Commands
│   ├── /sync-data           Refresh from Google Sheets
│   ├── /gibson-assemble     Design and simulate assembly
│   ├── /council-review      10-expert parallel review
│   ├── /validate-sequences  Check synthesis compatibility
│   ├── /cut-backbone        Restriction digest simulation
│   ├── /generate-twist-order Create order CSV
│   └── /design-fragment     Interactive fragment design
│
├── Assembly Council (10 parallel agents)
│   ├── thermodynamic-expert    haiku     Overlap Tm, annealing
│   ├── expression-expert       haiku     Codon usage, stability
│   ├── primer-designer         haiku     Verification primers, colony PCR
│   ├── manufacturing-expert    haiku     Twist constraints
│   ├── annotation-expert       haiku     GenBank compliance
│   ├── strategy-expert         sonnet    Alternative approaches
│   ├── cost-optimizer          haiku     Pricing analysis
│   ├── timeline-expert         haiku     Project scheduling
│   ├── host-compatibility      sonnet    Expression hosts
│   └── troubleshooting-expert  sonnet    Diagnose failures, suggest fixes
│
├── Hooks
│   ├── SessionStart           Cache freshness check
│   ├── PostToolUse            GenBank validation, audit log
│   └── PreToolUse             Write protection
│
└── MCP
    └── google-sheets          Fragment, plasmid, backbone DBs
```

---

## Project Structure

```
.claude/
├── agents/assembly-council/   10 expert agents
├── commands/                  7 slash commands
├── skills/                    3 knowledge domains
└── settings.json              Hooks and permissions

src/
├── core/                      Assembly, backbone, overlap, validation
├── io/                        GenBank, Sheets sync, Twist orders
└── hooks/                     Session and validation scripts

tests/                         Pytest suite with fixtures
```

---

## Key Constraints

| Parameter | Min | Max | Optimal |
|-----------|-----|-----|---------|
| Overlap size | 15 bp | 70 bp | 25 bp |
| Overlap Tm | 46°C | 60°C | 52°C |
| Tm delta | — | 15°C | — |
| Fragment size (Twist) | 300 bp | 1800 bp | — |

---

## Requirements

- Python 3.10+
- BioPython 1.77
- Google Sheets MCP server (for database sync)

---

## Citation

```
Lu, A. (2025). Gibson Assembly Pipeline. California Institute of Technology.
```

---

## License

MIT

---

Built with [Claude Code](https://claude.ai/claude-code)


# Dashboard

The interactive dashboard provides a complete browser-based analysis environment. It does not require any coding — all parameters are configurable through the UI.

---

## Launching

```bash
# From the repository root
python insilico_pcr/webapp/run.py

# Or with uvicorn directly
uvicorn insilico_pcr.webapp.app:app --host 0.0.0.0 --port 8765 --reload
```

Open `http://localhost:8765` in any modern browser.

---

## Quick Demo

1. Click **"Load demo"** — runs 2 primer pairs against a 10 kbp synthetic genome
2. Explore all tabs — results populate immediately
3. Go to **"Live Parameters"** tab — adjust mismatch slider, click "Re-run"

---

## Analysis Panels

### Run Setup

| Field | Description |
|---|---|
| Pair name | Label for this primer pair |
| Forward primer | 5′→3′ sequence, uppercase ACGT |
| Reverse primer | 5′→3′ sequence, uppercase ACGT |
| Genome source | Paste sequence, or upload FASTA file |
| Max mismatches | 0–6; controls search sensitivity |
| Min/Max amplicon | Size filter in base pairs |
| [Na⁺] | Total sodium concentration in mM |
| [Mg²⁺] | Total magnesium concentration in mM (0 = Na⁺-only mode) |
| [dNTP] | dNTP concentration in mM (chelates Mg²⁺ 1:1) |
| Primer conc | Total primer strand concentration in nM |
| 3′ strict mode | Reject primers with 3′-terminal mismatches |
| Hairpin analysis | Compute hairpin ΔG for each primer |
| Primer dimer | Detect homo/hetero dimers via NN model |

### Amplicons & Alignment

Ranked amplicon table with:
- Position, length, GC%, composite score, mismatch badges
- Click any row to view the nucleotide alignment in the panel below
- Alignment viewer shows primer-to-template base pairing with mismatches highlighted in red

### Thermodynamics

- Tm comparison bar chart (all primer pairs side-by-side)
- Primer quality radar chart (Tm, ΔG, GC, length, clamp)
- Per-amplicon binding site details: Tm, ΔG, SW score, binding score

### Off-target Explorer

- Total hit count and specificity index
- Risk classification per hit (high/medium/low)
- Annotated reasons (unexpected position, poor Tm, etc.)
- Score distribution histogram

### Primer Quality

For each primer (forward and reverse):
- Tm, ΔG₃₇, GC fraction, length
- GC clamp status
- Maximum homopolymer run
- Low-complexity detection

### Genome Overview

- Amplicon positions plotted at genomic coordinates
- Amplicon length histogram
- GC content vs length scatter

### Live Parameters

Adjust parameters and re-run without leaving the page:
- Mismatches slider (0–6)
- [Na⁺] slider (10–200 mM)
- [Mg²⁺] slider (0–10 mM)
- Edit primer sequences directly

Results update all panels on re-run (no page reload).

### Export

- **JSON**: full structured output (machine-readable)
- **CSV**: amplicon table (spreadsheet-compatible)
- **HTML**: self-contained report with raw JSON

---

## API Endpoints

The dashboard backend exposes a REST API:

| Method | Endpoint | Description |
|---|---|---|
| GET | `/` | Dashboard HTML |
| POST | `/api/run` | Run full pipeline |
| GET | `/api/demo` | Load pre-computed demo results |
| GET | `/api/demo/genome` | Fetch demo genome string for live panel |
| POST | `/api/params` | Live re-run with new parameters |
| GET | `/api/benchmark` | Return cached benchmark metrics |
| POST | `/api/export/json` | Export results as JSON download |
| POST | `/api/export/csv` | Export amplicon table as CSV |

### POST /api/params — Live re-run

```json
{
  "primer_name":   "ACTB",
  "fwd_primer":    "GCACTGGTGGCATCGATCTA",
  "rev_primer":    "TAGCTAGCATGCTAGCTAGC",
  "genome_string": "ATCGATCG...",
  "mismatches":    3,
  "na_conc":       50.0,
  "mg_conc":       2.0,
  "min_size":      50,
  "max_size":      3000
}
```

Response:
```json
{
  "success": true,
  "results": { "primer_pairs": [...] },
  "summary": {
    "amplicons_found": 3,
    "offtargets": 1,
    "top_score": 35.05
  },
  "elapsed": 0.063
}
```

---

## Technology Stack

| Layer | Technology |
|---|---|
| Backend | FastAPI (Python) |
| Frontend | Vanilla JavaScript (no framework) |
| Plots | Plotly.js (served locally) |
| Fonts | Inter + JetBrains Mono (Google Fonts) |
| Theme | Custom scientific light theme (CSS variables) |
| Templates | Jinja2 |

---

## Screenshots

See `assets/screenshots/` for dashboard screenshots. To regenerate:

1. Run `python insilico_pcr/webapp/run.py`
2. Load the demo
3. Take screenshots at 1440×900 px
4. Save as `assets/screenshots/dashboard_overview.png`, `alignment_viewer.png`, etc.

Recommended screenshot names:

| File | Panel |
|---|---|
| `dashboard_overview.png` | Full dashboard, run panel visible |
| `alignment_viewer.png` | Amplicons tab, alignment panel open |
| `thermodynamics_panel.png` | Thermodynamics tab |
| `offtarget_explorer.png` | Off-target tab |
| `genome_track.png` | Genome overview tab |
| `live_parameters.png` | Live parameters tab with sliders |

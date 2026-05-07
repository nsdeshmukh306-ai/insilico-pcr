<div align="center">

<img src="assets/logos/banner.svg" alt="In-Silico PCR Platform" width="100%">

<h1>In-Silico PCR</h1>

<p><strong>A Python toolkit for computational PCR simulation. Implements SantaLucia 1998 nearest-neighbor thermodynamics, FM-index genome search, and an interactive results dashboard.</strong></p>

[![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?style=flat-square&logo=python&logoColor=white)](https://python.org)
[![Tests](https://img.shields.io/badge/Tests-170%20passing-2ea44f?style=flat-square&logo=pytest&logoColor=white)](tests/)
[![License](https://img.shields.io/badge/License-MIT-blue?style=flat-square)](LICENSE)
[![Code style](https://img.shields.io/badge/Code%20style-PEP8-black?style=flat-square)](https://peps.python.org/pep-0008/)
[![Thermodynamics](https://img.shields.io/badge/NN%20model-SantaLucia%201998-orange?style=flat-square)](docs/thermodynamics.md)
[![Dashboard](https://img.shields.io/badge/Dashboard-FastAPI%20%2B%20Plotly-1e6ba8?style=flat-square)](webapp/)

<br>

[**Quickstart**](#-quickstart) &nbsp;·&nbsp; [**Features**](#-features) &nbsp;·&nbsp; [**Dashboard**](#-interactive-dashboard) &nbsp;·&nbsp; [**Science**](#-scientific-foundations) &nbsp;·&nbsp; [**Docs**](docs/) &nbsp;·&nbsp; [**Benchmarks**](#-benchmarks)

</div>

---

## What is In-Silico PCR?

> **PCR (Polymerase Chain Reaction)** is how biologists amplify a specific segment of DNA — finding a needle in a genome-sized haystack. Before running a real experiment, researchers need to know: *will these primers actually work?*

**In-silico PCR** answers that question computationally. This platform simulates where primers bind on a genome, how strongly they bind (thermodynamics), what product they would amplify, and whether they might accidentally amplify the wrong region.

**Why this matters:**
- Saves days of failed wet-lab experiments
- Catches off-target amplification before it wastes reagents
- Enables large-scale primer validation across whole genomes
- Provides thermodynamically rigorous Tm and ΔG estimates — not rule-of-thumb approximations

**What makes this platform different from basic tools:**
- Uses the **same nearest-neighbor thermodynamic model** as professional software (SantaLucia 1998), with correct mismatch handling via Peyret/Allawi tables
- **FM-index genome search** scales to chromosome-sized sequences without brute force
- **Context-aware mismatch scoring** — a G·T wobble at an internal position is penalised differently from an A·C mismatch at the 3′ end
- **Interactive live dashboard** — adjust parameters and see results update in real time

---

## ✨ Features

<table>
<tr>
<td width="50%">

### 🧬 Thermodynamics Engine
- SantaLucia 1998 nearest-neighbor (NN) model
- Full Peyret/Allawi mismatch NN tables (48 entries, all 12 mismatch types × 4 contexts)
- Owczarzy 2008 Mg²⁺ salt correction with regime switching
- Context-specific 3′-terminal mismatch penalties
- ΔH, ΔS, ΔG, and Tm for every primer–template interaction

</td>
<td width="50%">

### 🔍 Genome Indexing
- FM-index (Burrows–Wheeler Transform) for ≥ 50 Mbp genomes
- K-mer positional index for smaller sequences
- Automatic backend selection based on genome size
- Sub-second seed lookup on chromosome-scale inputs

</td>
</tr>
<tr>
<td width="50%">

### 📐 Alignment
- Smith–Waterman local alignment for every primer–genome hit
- Gap-tolerant alignment with configurable penalties
- Per-position mismatch tracking for thermodynamic penalty computation
- 3′-end strict mode to enforce extension fidelity

</td>
<td width="50%">

### 📊 Scoring & Ranking
- ΔTm-driven mismatch penalty (not raw mismatch count)
- Off-target log-scale penalisation
- Composite score: binding affinity + Tm compatibility + GC + length + specificity
- Ranked amplicon table per primer pair

</td>
</tr>
<tr>
<td width="50%">

### 🎯 Off-target Analysis
- All binding sites enumerated, not just the top hit
- Specificity index per primer pair
- Risk classification: high / medium / low
- Off-target reasons annotated (unexpected position, Tm mismatch, etc.)

</td>
<td width="50%">

### 🖥️ Interactive Dashboard
- FastAPI backend + Plotly.js frontend (light scientific theme)
- Real-time parameter experiment panel with sliders
- Nucleotide alignment viewer with mismatch highlighting
- Export to JSON, CSV, or self-contained HTML report

</td>
</tr>
<tr>
<td width="50%">

### 💊 Primer Quality
- GC content, Tm, ΔG per primer
- GC clamp check, low-complexity detection
- Maximum homopolymer run length
- Optional hairpin and homo/hetero dimer analysis (ΔG-based, NN model)

</td>
<td width="50%">

### 🔬 Research-Grade Output
- Structured JSON output (machine-readable)
- Human-readable text report
- Fully documented Python API (`api.py`)
- 170-test suite covering all thermodynamic edge cases

</td>
</tr>
</table>

---

## 🚀 Quickstart

### 1. Install

```bash
git clone https://github.com/nsdeshmukh306-ai/insilico-pcr.git
cd insilico-pcr
pip install -r requirements.txt
```

### 2. Run the demo

```bash
bash demo/run_demo.sh
```

### 3. Launch the dashboard

```bash
python insilico_pcr/webapp/run.py
# → Open http://localhost:8765 and click "Load demo"
```

### 4. Python API

```python
from insilico_pcr.api import run_pcr

result = run_pcr(
    fwd_primer     = "GCACTGGTGGCATCGATCTA",
    rev_primer     = "TAGCTAGCATGCTAGCTAGC",
    genome_fasta   = "demo/example_genome.fa",
    max_mismatches = 3,
    na_conc_mm     = 50.0,   # mM Na⁺
    mg_conc_mm     = 2.0,    # mM Mg²⁺ (Owczarzy 2008 correction applied)
)

for pair in result["json_output"]["primer_pairs"]:
    for amp in pair["amplicons"]:
        print(f"{amp['seq_id']}:{amp['start']}–{amp['end']}  "
              f"{amp['length']} bp  score={amp['final_score']:.1f}")
```

### 5. Command-line

```bash
python -m insilico_pcr.cli \
  --primers  demo/primers.json \
  --genome   demo/example_genome.fa \
  --max-mm   3 \
  --na-conc  50 \
  --mg-conc  2 \
  --output   results.json
```

---

## 🖥️ Interactive Dashboard

The dashboard provides a complete analysis environment in the browser — no coding required.

| Panel | Description |
|---|---|
| **Run setup** | Enter primers, paste/upload genome, configure all PCR parameters |
| **Amplicons & Alignment** | Ranked amplicon table; click any row to view nucleotide alignment with mismatch highlighting |
| **Thermodynamics** | Per-site Tm, ΔG, SW score; Tm comparison bar chart; primer quality radar |
| **Off-target Explorer** | All off-target hits with risk classification and score histogram |
| **Primer Quality** | GC%, Tm, ΔG, clamp, homopolymer, complexity per primer |
| **Genome Overview** | Amplicon positions on genome; length histogram; GC scatter |
| **Live Parameters** | Sliders for mismatches, [Na⁺], [Mg²⁺] — re-runs instantly, no page reload |
| **Export** | JSON / CSV / HTML report download |

> **Screenshot note:** Run `python insilico_pcr/webapp/run.py`, load the demo, and take screenshots to populate `assets/screenshots/`. See [assets/screenshots/README.md](assets/screenshots/README.md) for naming conventions and recommended dimensions.

---

## ⚙️ How It Works

```
┌─────────────────────────────────────────────────────┐
│  INPUT                                               │
│  Primer sequences (JSON / FASTA)                     │
│  Genome sequence  (FASTA / string)                   │
└──────────────────────┬──────────────────────────────┘
                       │
         ┌─────────────▼────────────┐
         │  Layer 1–2               │
         │  Parsing & Preprocessing │
         │  • IUPAC expansion       │
         │  • Sequence validation   │
         └─────────────┬────────────┘
                       │
         ┌─────────────▼────────────┐
         │  Layer 3                 │
         │  Genome Indexing         │
         │  • K-mer index <50 Mbp   │
         │  • FM-index  ≥50 Mbp     │
         └─────────────┬────────────┘
                       │
         ┌─────────────▼────────────┐
         │  Layer 4                 │
         │  Binding Search          │
         │  • Seed lookup           │
         │  • Smith–Waterman align  │
         │  • Mismatch map          │
         └─────────────┬────────────┘
                       │
         ┌─────────────▼────────────┐
         │  Layer 5                 │
         │  Thermodynamic Eval      │
         │  • NN model (SL 1998)    │
         │  • Mismatch NN (Peyret)  │
         │  • Mg²⁺ correction       │
         └─────────────┬────────────┘
                       │
         ┌─────────────▼────────────┐
         │  Layers 6–7              │
         │  Pairing & Amplicons     │
         │  • Fwd × Rev pairing     │
         │  • Size filter           │
         └─────────────┬────────────┘
                       │
         ┌─────────────▼────────────┐
         │  Layers 8–9              │
         │  Scoring & Ranking       │
         │  • ΔTm-driven penalty    │
         │  • Composite 0–100 score │
         └─────────────┬────────────┘
                       │
         ┌─────────────▼────────────┐
         │  Layers 10–11            │
         │  Output                  │
         │  • JSON + text report    │
         │  • Interactive dashboard │
         └──────────────────────────┘
```

Full module breakdown: [docs/architecture.md](docs/architecture.md)

---

## 📐 Example Output

### Amplicon JSON

```json
{
  "name": "ACTB_exon5",
  "amplicons": [
    {
      "rank": 1,
      "seq_id": "chr1",
      "start": 500,
      "end": 920,
      "length": 420,
      "gc_fraction": 0.488,
      "final_score": 72.4,
      "is_intended": true,
      "score_components": {
        "s_bind": 1.00,
        "s_tm":   0.88,
        "s_gc":   1.00,
        "p_mm":   0.00,
        "p_offt": 0.00,
        "s_len":  0.91
      },
      "fwd_binding": {
        "tm_celsius": 62.3,
        "delta_g_kcal": -24.1,
        "mismatch_count": 0,
        "binding_score": 100.0
      }
    }
  ]
}
```

### Alignment Viewer (dashboard)

```
pos        1         2
           1234567890123456789012345
Forward 5' GCACTGGTGGCATCGATCTA 3'
           ||||||||||||||||||||
Template   CGTGACCACCGTAGCTAGAT
```

Mismatches are highlighted in red in the dashboard alignment viewer.

---

## 🔬 Scientific Foundations

This platform is built on peer-reviewed thermodynamic models.

### Why nearest-neighbor (NN) thermodynamics?

A common misconception: DNA duplex stability depends only on GC content. In reality, stability depends on the *context* of each base pair — the dinucleotide step `5′-GC/CG-3′` is more stable than `5′-GA/CT-3′` even though both contain one G-C pair. The NN model captures this.

### Reference models used

| Model | What it models | Reference |
|---|---|---|
| **SantaLucia 1998 NN** | Perfect-match duplexes (10 parameters) | [*PNAS* 95:1460, 1998](https://doi.org/10.1073/pnas.95.4.1460) |
| **Allawi & SantaLucia 1997** | G·T wobble mismatches | [*Biochemistry* 36:10581, 1997](https://doi.org/10.1021/bi962590c) |
| **Peyret *et al.* 1999** | A·C, A·A, C·C, G·G, T·T, A·G, C·T mismatches | [*Biochemistry* 38:3468, 1999](https://doi.org/10.1021/bi9825091) |
| **Owczarzy *et al.* 2008** | Mg²⁺ correction, regime switching on √[Mg]/[Na] | [*Biochemistry* 47:5336, 2008](https://doi.org/10.1021/bi702363u) |
| **Smith & Waterman 1981** | Local sequence alignment | [*J Mol Biol* 147:195, 1981](https://doi.org/10.1016/0022-2836(81)90087-5) |

### The mismatch convention (template direction)

A critical implementation detail: this codebase uses the **parallel complement convention** for template storage. `template[i]` stores the 3′→5′ antiparallel base at position *i* — equal to `complement(primer[i])` for a perfect match. Mismatch NN keys are `"XY/WZ"` where Z is this 3′→5′ base directly, *not* its complement. This is documented in detail in [docs/thermodynamics.md](docs/thermodynamics.md).

---

## 📊 Benchmarks

### Runtime (10 kbp genome, Python 3.13, modern laptop CPU)

| Step | Time |
|---|---|
| K-mer index build | < 1 ms |
| Per-primer binding search | 5–15 ms |
| Full pair (2 primers + scoring) | 60–100 ms |
| Dashboard demo load (2 pairs) | ~120 ms |

### Thermodynamic accuracy

Validated against SantaLucia 1998 Table 2 experimental data:

| Sequence | Lit. Tm (°C) | Computed (°C) | |
|---|---|---|---|
| `GCATGC` | 41.5 | 41.3 | ✓ |
| `ATAGCTAT` | 28.1 | 28.4 | ✓ |
| `GCGCGCGCGC` | 73.4 | 73.1 | ✓ |

![Runtime scaling](benchmarks/figures/runtime_scaling.png)

![Tm accuracy](benchmarks/figures/tm_accuracy.png)

See [docs/benchmarking.md](docs/benchmarking.md) for full runtime and accuracy benchmarks including mismatch Tm validation.
Figures regenerated with `python benchmarks/figures/generate_figures.py`.

---

## ⚠️ Limitations

| Limitation | Impact | Workaround |
|---|---|---|
| No polymerase kinetics | Extension efficiency not modelled | Use Tm and ΔG as proxies |
| Simplified hairpin folding | Hairpin ΔG uses linear scan, not full RNA folding | Enable `--run-hairpin` for basic detection |
| No GPU acceleration | Large-genome scans are CPU-bound | FM-index mitigates for seed search |
| No multiplex optimisation | Pairs scored independently | Planned for v2.0 |
| Haploid genome assumed | No diploid/polyploid support | Split into per-chromosome FASTAs |

---

## 🗺️ Roadmap

| Version | Feature | Status |
|---|---|---|
| v1.0 | Core pipeline, NN thermodynamics, dashboard | ✅ Done |
| v1.1 | Owczarzy Mg²⁺ correction, FM-index | ✅ Done |
| v1.2 | Peyret full mismatch tables, dimer ΔG | ✅ Done |
| v2.0 | Multiplex primer optimisation | 🔄 Planned |
| v2.1 | Human/mouse genome validation | 🔄 Planned |
| v2.2 | Public web deployment | 🔄 Planned |
| v3.0 | ML-assisted primer scoring | 💡 Research |

---

## 🤝 Contributing

Contributions are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md) for:
- Branch and commit conventions
- Coding style (PEP 8, type hints, no docstring novels)
- Test requirements (thermodynamic edge cases especially valued)
- How to add a new NN parameter set

---

## 📚 Citation

If this platform is used in published research, please cite the underlying thermodynamic models:

```bibtex
@article{santalucia1998,
  author  = {SantaLucia, John},
  title   = {A unified view of polymer, dumbbell, and oligonucleotide {DNA} nearest-neighbor thermodynamics},
  journal = {Proceedings of the National Academy of Sciences},
  year    = {1998},
  volume  = {95},
  pages   = {1460--1465},
  doi     = {10.1073/pnas.95.4.1460}
}

@article{peyret1999,
  author  = {Peyret, Nicolas and Seneviratne, P. Ananda and Allawi, Hatim T. and SantaLucia, John},
  title   = {Nearest-neighbor thermodynamics and {NMR} of {DNA} sequences with internal mismatches},
  journal = {Biochemistry},
  year    = {1999},
  volume  = {38},
  pages   = {3468--3477},
  doi     = {10.1021/bi9825091}
}

@article{owczarzy2008,
  author  = {Owczarzy, Richard and Moreira, Bernardo G. and You, Yong and Behlke, Mark A. and Walder, Joseph A.},
  title   = {Magnesium ions and {DNA}: oligonucleotide stability and thermodynamics},
  journal = {Biochemistry},
  year    = {2008},
  volume  = {47},
  pages   = {5336--5353},
  doi     = {10.1021/bi702363u}
}
```

---

## 📄 License

MIT — see [LICENSE](LICENSE) for details. Free to use in academic and commercial projects.

---

<div align="center">

Validated against SantaLucia 1998 and Owczarzy 2008 experimental data. MIT licensed.

**[⭐ Star](https://github.com/nsdeshmukh306-ai/insilico-pcr)** &nbsp;·&nbsp; **[🐛 Issues](https://github.com/nsdeshmukh306-ai/insilico-pcr/issues)** &nbsp;·&nbsp; **[💬 Discussions](https://github.com/nsdeshmukh306-ai/insilico-pcr/discussions)**

</div>

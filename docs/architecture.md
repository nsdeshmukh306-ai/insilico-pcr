# Architecture

The pipeline is structured as 12 discrete, independently testable layers. Each layer has a single responsibility and communicates through typed data structures.

---

## Layer Map

```
Layer  Module                  Responsibility
─────  ──────────────────────  ────────────────────────────────────────────
  1    input_handler.py        Parse FASTA/JSON, validate sequences
  2    preprocessor.py         IUPAC expansion, sanitisation, deduplication
  3    genome_index.py         K-mer positional index (< 50 Mbp)
  3b   genome_index_fm.py      FM-index / BWT index  (≥ 50 Mbp)
  4    binding_search.py       Seed lookup + Smith–Waterman alignment
  5    binding_eval.py         Per-site thermodynamic evaluation (Tm, ΔG)
  6    pairing_engine.py       Forward × reverse site pairing + orientation check
  7    amplicon.py             Amplicon sequence extraction, length filter
  8    thermodynamics.py       NN model, mismatch NN, salt correction
  9    scoring.py              Composite 0–100 score, ranking
 10    offtarget.py            Off-target classification, specificity index
 11    output_handler.py       JSON + text report generation
 12    advanced/               Hairpin detection, primer dimer ΔG
```

---

## Data Flow

```
PrimerInput (Layer 1)
    │
    ▼
ProcessedPrimer (Layer 2) — IUPAC variants, strand, GC, Tm
    │
    ▼
GenomeIndex (Layer 3/3b) — k-mer dict or FM-index
    │
    ▼
BindingSite[] (Layer 4) — per-hit: position, strand, SW score, mismatch_positions
    │
    ▼
EvaluatedSite[] (Layer 5) — Tm, ΔG, ΔH, ΔS added to each site
    │
    ▼
AmpliconHit[] (Layer 6) — fwd_site + rev_site pair, orientation confirmed
    │
    ▼
Amplicon[] (Layer 7) — sequence extracted, length filtered
    │
    ▼
ScoredAmplicon[] (Layer 9) — composite score computed
    │
    ▼
OfftargetSummary (Layer 10) — risk classification, specificity index
    │
    ▼
OutputJSON + TextReport (Layer 11)
```

---

## Key Data Classes

```python
# modules/amplicon.py
@dataclass
class Amplicon:
    seq_id: str
    start: int
    end: int
    length: int
    sequence: str
    gc_fraction: float
    hit: AmpliconHit         # contains fwd_site, rev_site
    fwd_binding_score: float
    rev_binding_score: float
    fwd_tm: float
    rev_tm: float

# modules/scoring.py
@dataclass
class ScoredAmplicon:
    amplicon: Amplicon
    s_bind: float            # 0–1
    s_tm: float              # 0–1
    s_gc: float              # 0–1
    p_mm: float              # 0–1 (penalty)
    p_offt: float            # 0–1 (penalty)
    s_len: float             # 0–1
    final_score: float       # 0–100
    is_intended: bool
```

---

## Genome Index Selection

```python
GENOME_SIZE_THRESHOLD = 50_000_000  # 50 Mbp

if genome_length >= GENOME_SIZE_THRESHOLD:
    index = build_fm_index(genome)      # BWT-based, O(n) space
else:
    index = build_kmer_index(genome)    # dict-based, O(n·k) space
```

The FM-index uses suffix arrays with a compressed BWT representation. Seed lookup is O(m log n) for a primer of length m on genome of length n. The k-mer index is a dict mapping each k-mer to a list of genome positions — O(1) lookup but O(n) memory.

---

## Thermodynamics Integration

`binding_eval.py` calls `thermodynamics.py` for every binding site. The call passes:
- `primer_seq`: the primer sequence (5′→3′)
- `template`: the aligned genome region in **parallel complement convention** (3′→5′ direction)
- `mismatch_positions`: list of 0-based positions where `template[i] ≠ complement(primer[i])`
- `three_prime_mismatch`: bool, triggers extra ΔH penalty

This keeps thermodynamics as a pure function with no genome I/O.

---

## Scoring Formula

```
Score = (w_bind·S_bind + w_tm·S_tm + w_gc·S_gc + w_len·S_len
         - w_mm·P_mm - w_offt·P_offt)
        ────────────────────────────────────────────────────── × 100
                    (w_bind + w_tm + w_gc + w_len)

Default weights: w_bind=0.30, w_tm=0.25, w_gc=0.10, w_len=0.05
                 w_mm=0.20,   w_offt=0.10
```

`P_mm` is ΔTm-driven, not count-based:
```
P_mm = min(1, (ΔTm_fwd + ΔTm_rev) / 20°C)
```

See [scoring.md](scoring.md) for full derivation.

---

## Testing Strategy

Each layer has a dedicated test module. Thermodynamic layers are tested against published literature values (SantaLucia 1998 Table 2, Peyret 1999 supplementary):

```
tests/
├── test_thermodynamics.py           # Perfect-match NN, Tm
├── test_mismatch_thermodynamics.py  # Peyret mismatch tables
├── test_salt_mg.py                  # Owczarzy 2008 Mg²⁺ correction
├── test_fm_index.py                 # FM-index correctness
├── test_dimer_dg.py                 # Primer dimer ΔG
├── test_binding_search.py           # SW alignment
├── test_genome_index.py             # K-mer index
├── test_pairing_scoring.py          # Pairing + scoring
├── test_integration.py              # End-to-end pipeline
└── test_advanced.py                 # Hairpin, dimer
```

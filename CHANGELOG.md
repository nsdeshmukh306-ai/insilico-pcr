# Changelog

All notable changes to this project are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [Unreleased]

### Planned
- Multiplex primer optimisation
- Human/mouse genome validation suite
- Public web deployment

---

## [1.2.0] — 2026-05-07

### Added
- **Full Peyret/Allawi mismatch NN table** — 48 entries covering all 12 mismatch types × 4 contexts (`_MM_NN_PARAMS` in `thermodynamics.py`)
- **Owczarzy 2008 Mg²⁺ salt correction** — regime switching on √[Mg²⁺]/[Na⁺] ratio; `calc_tm()` accepts `mg_conc` and `dntp_conc`
- **FM-index genome indexing** (`modules/genome_index_fm.py`) — auto-switches to BWT-based index for genomes ≥ 50 Mbp
- **ΔG-based primer dimer detection** (`modules/advanced/primer_dimer.py`) — NN model for duplex ΔG, strong 3′-end penalty
- **ΔTm-driven mismatch penalty** in scoring — replaces old count-based formula; `P_mm = min(1, ΔTm_combined / 20°C)`
- **Interactive FastAPI dashboard** (`webapp/`) — Plotly.js, light scientific theme, 10 analysis panels
- **Live parameter experiment panel** — re-runs pipeline with slider-adjusted params, no page reload
- **Expanded test suite** — 4 new test files: mismatch thermodynamics, Mg²⁺ salt, FM-index, dimer ΔG (170 tests total)

### Fixed
- **Critical bug in `calc_nn_thermodynamics`** — was passing `complement(t3)` to `_mm_nn_lookup` instead of `t3` directly (template stored in 3′→5′ parallel complement convention; Z = t3, not complement(t3))
- `test_ac_mismatch_ca_context` was calling `_mm_nn_lookup("CA", "G")` (A·G mismatch) instead of correct `_mm_nn_lookup("CA", "C")` (A·C mismatch)
- `test_gt_mismatch_less_destabilising_than_ac` used `'A'` in template for G·T wobble (gives G·A mismatch); corrected to `'T'`
- `test_multiple_mismatches_lower_tm_progressively` mutated position 0 (5′-terminal, not captured by NN loop); fixed to start at position 1

### Changed
- Scoring `P_mm` formula changed from count-based to ΔTm-driven (see `modules/scoring.py` for rationale comment)
- `/api/params` backend model renamed `LivePCRRequest` — removed stale `result_json` required field; frontend payload now matches backend schema exactly

---

## [1.1.0] — 2026-04-15

### Added
- Smith–Waterman alignment in `modules/binding_search.py`
- Off-target risk classification (high / medium / low) in `modules/offtarget.py`
- Specificity index per primer pair
- JSON output schema (`modules/output_handler.py`)
- CLI entry point (`cli.py`)
- Python API (`api.py`)

### Fixed
- Primer complement direction bug in pairing engine

---

## [1.0.0] — 2026-03-01

### Added
- Initial 12-layer pipeline architecture
- SantaLucia 1998 NN thermodynamics (perfect match)
- K-mer positional genome index
- Basic Tm calculation with Na⁺-only salt correction
- Amplicon extraction and size filtering
- Basic composite scoring
- Text report output
- Initial test suite (100 tests)

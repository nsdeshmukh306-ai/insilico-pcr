# Contributing to In-Silico PCR

This project is research software. Correctness and reproducibility come before feature velocity.

---

## Types of Contributions Welcome

| Type | Notes |
|---|---|
| **Bug fixes** | Especially thermodynamic calculation errors |
| **New test cases** | Edge cases for mismatch NN, salt correction |
| **Large-genome benchmarks** | Human / mouse genome validation results |
| **Documentation** | Clearer explanations, diagrams, examples |
| **Dashboard UI** | Clean, scientific aesthetic only |
| **Performance** | Profiling data required with PR |

---

## Getting Started

```bash
# Fork and clone
git clone https://github.com/nsdeshmukh306-ai/insilico-pcr.git
cd insilico-pcr

# Create a branch
git checkout -b fix/mismatch-tm-calculation

# Install in editable mode with dev deps
pip install -e .
pip install pytest pytest-cov ruff
```

---

## Code Style

- **PEP 8** enforced via `ruff`
- **Type hints** on all public function signatures
- **No docstring walls** — one short line is better than five vague sentences
- **No comments** explaining *what* code does — only *why* (invariants, workarounds, hidden constraints)
- Line length: 100 characters

Run linter before committing:

```bash
ruff check insilico_pcr/ tests/
```

---

## Test Requirements

All PRs must maintain the full test suite passing:

```bash
pytest tests/ -v
```

**For thermodynamic changes specifically:**
- Add a test in `tests/test_thermodynamics.py` or `tests/test_mismatch_thermodynamics.py`
- Cite the literature value you are validating against
- Include the expected ΔH, ΔS, and (where possible) Tm

**Template convention reminder:** `template[i]` stores the 3′→5′ antiparallel base at position *i* (the "parallel complement convention"). See `docs/thermodynamics.md` for details. This trips up contributors — read the docs before touching `modules/thermodynamics.py`.

---

## Branch Conventions

| Branch prefix | Use |
|---|---|
| `fix/` | Bug fixes |
| `feat/` | New features |
| `docs/` | Documentation only |
| `bench/` | Benchmarks |
| `refactor/` | No behaviour change |

---

## Pull Request Checklist

- [ ] Tests pass (`pytest tests/`)
- [ ] Linter clean (`ruff check`)
- [ ] New behaviour covered by a test
- [ ] CHANGELOG.md updated under `[Unreleased]`
- [ ] Literature references cited if thermodynamics is touched

---

## Adding a New NN Parameter Set

1. Add entries to `_MM_NN_PARAMS` in `modules/thermodynamics.py`
2. Follow the existing key format: `"XY/WZ"` (primer 5′→3′ / template 3′→5′)
3. Add symmetry fallback handling in `_mm_nn_lookup()`
4. Cite the exact table, row, and paper for every value
5. Add tests that validate at least 3 lookup values from the paper

---

## Reporting Bugs

For thermodynamic bugs, please include:
- Primer sequence
- Template sequence (and which direction it is stored)
- Expected Tm / ΔH / ΔS (and the reference paper)
- Observed output

For dashboard bugs, include the browser console log and the JSON payload sent.

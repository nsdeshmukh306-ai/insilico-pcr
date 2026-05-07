# Thermodynamics

This document explains the thermodynamic models implemented in `modules/thermodynamics.py`.

---

## Overview

The melting temperature (Tm) of a primer–template duplex is calculated using the **nearest-neighbor (NN) model**, which is the standard method used by Primer3, OligoCalc, and IDT OligoAnalyzer. This implementation uses:

1. **SantaLucia 1998** — unified NN parameters for perfect-match duplexes
2. **Peyret *et al.* 1999 + Allawi & SantaLucia 1997** — context-specific mismatch NN parameters
3. **Owczarzy *et al.* 2008** — Mg²⁺ salt correction with regime switching

---

## 1. Perfect-Match Duplexes (SantaLucia 1998)

### Why nearest-neighbor?

A naive approach estimates Tm from GC content only (e.g. the 4°C/GC + 2°C/AT rule). This is wrong for short oligonucleotides. The NN model accounts for the fact that base-pair stability depends on both the base and its immediate neighbours.

### The 10 NN parameters

There are 10 unique dinucleotide steps (under Watson-Crick symmetry). Each has experimentally measured ΔH° and ΔS° values from optical melting experiments:

| Sequence | ΔH° (kcal/mol) | ΔS° (cal/mol·K) |
|---|---|---|
| AA/TT | −7.9 | −22.2 |
| AT/TA | −7.2 | −20.4 |
| TA/AT | −7.2 | −21.3 |
| CA/GT | −8.5 | −22.7 |
| GT/CA | −8.4 | −22.4 |
| CT/GA | −7.8 | −21.0 |
| GA/CT | −8.2 | −22.2 |
| CG/GC | −10.6 | −27.2 |
| GC/CG | −9.8 | −24.4 |
| GG/CC | −8.0 | −19.9 |

*Source: SantaLucia, PNAS 95:1460, 1998, Table 2*

### Tm formula

```
        ΔH°_total
Tm = ─────────────────────── − 273.15   (°C)
      ΔS°_total + R·ln(C_T/4)
```

Where:
- `ΔH°_total` = sum of all NN dinucleotide ΔH values + terminal correction
- `ΔS°_total` = sum of all NN dinucleotide ΔS values + terminal correction
- `R` = 1.987 cal/mol·K (gas constant)
- `C_T` = total strand concentration (typically 250 nM for primer design)

Terminal corrections (initiation parameters): `+0.2 kcal/mol ΔH` and `−5.7 cal/mol·K ΔS` per terminal GC; `+2.2 kcal/mol ΔH` and `+6.9 cal/mol·K ΔS` per terminal AT.

### ΔG at 37°C

```python
delta_g = delta_h - (310.15 * delta_s / 1000)   # kcal/mol
```

---

## 2. Mismatch Nearest-Neighbor Model

When a primer has mismatches relative to the template, the matched NN parameters are replaced with mismatch-specific values from Peyret *et al.* 1999 (A·A, C·C, G·G, T·T, A·C, A·G, C·T) and Allawi & SantaLucia 1997 (G·T wobble).

### The 48-entry table

There are 12 mismatch types × 4 contexts (the flanking base on the 5′ side of the primer) = 48 unique entries in `_MM_NN_PARAMS`.

Key format: `"XY/WZ"`
- `X` = 5′ primer base of the step (flanking)
- `Y` = 3′ primer base of the step (the mismatch base)
- `W` = complement(X) = expected template base
- `Z` = actual template base (3′→5′ convention)

Example: `"AG/TT"` means primer `5′-AG-3′`, template `3′-TT-5′` → G·T wobble in A context.

### Template convention (critical)

This codebase uses the **parallel complement convention**: `template[i]` stores the *3′→5′ antiparallel base* at position *i*. For a perfect match, `template[i] = complement(primer[i])`. For a mismatch, `template[i]` is the actual base at that position in the 3′→5′ direction.

The mismatch NN lookup receives `t3 = template[i+1]` directly — **not** `complement(template[i+1])`. This is because the key format uses the 3′→5′ base as Z.

```python
# CORRECT — template[i] is already in 3′→5′ convention
nn_dh, nn_ds = _mm_nn_lookup(dinuc, t3)

# WRONG (previous bug) — double-complement corrupts the key
nn_dh, nn_ds = _mm_nn_lookup(dinuc, complement(t3))
```

The bug was fixed in v1.2.0. The fix was validated against Allawi & SantaLucia 1997 Table 1 and Peyret 1999 supplementary data.

### 3′-terminal mismatch penalty

A mismatch at the 3′ end of a primer is especially destabilising because it is directly adjacent to the extension site. An extra `+2.5 kcal/mol ΔH` penalty is applied when `three_prime_mismatch=True`.

---

## 3. Salt Correction

### Na⁺-only (Marmur–Schildkraut–Doty)

For Na⁺-only conditions, the Tm is adjusted by:

```
Tm_corrected = Tm_1M + (16.6 × log10([Na⁺]))
```

Where `Tm_1M` is the Tm at 1 M Na⁺ (computed from NN parameters).

### Mg²⁺ correction (Owczarzy 2008)

When Mg²⁺ is present, the Owczarzy 2008 model selects one of three regimes based on the ratio `√[Mg²⁺] / [Na⁺]`:

| Condition | Regime |
|---|---|
| `√[Mg] / [Na] < 0.22` | Na⁺-dominated — use standard Na correction |
| `0.22 ≤ √[Mg] / [Na] < 6.0` | Mixed — weighted combination |
| `√[Mg] / [Na] ≥ 6.0` | Mg²⁺-dominated — use pure Mg correction |

The Mg²⁺ correction coefficients are from Owczarzy et al. 2008 Eq. 16:

```
1/Tm_Mg = 1/Tm_Na + a + b·ln([Mg]) + (g + h·ln([Mg]))·ln²([Mg])
```

where `a, b, c, d, e, f, g, h` are empirically fitted constants.

**dNTP effect:** dNTPs chelate Mg²⁺ 1:1. Effective free Mg²⁺ is:
```
[Mg_free] = max(0, [Mg_total] - [dNTP])
```

---

## 4. Public API

```python
from insilico_pcr.modules.thermodynamics import (
    calc_nn_thermodynamics,   # (primer, template=None) → (ΔH, ΔS)
    calc_tm,                  # (primer, ...) → Tm °C
    calc_delta_g,             # (primer, ...) → ΔG kcal/mol
    gc_content,               # (seq) → float 0–1
    _mm_nn_lookup,            # (dinuc, t3) → (ΔH, ΔS)  [for testing]
    _nn_lookup,               # (dinuc) → (ΔH, ΔS)       [for testing]
)
```

### `calc_tm` full signature

```python
def calc_tm(
    primer_seq: str,
    template: str = None,          # enables mismatch NN; 3′→5′ convention
    mismatch_positions: list = None, # fallback if no template
    three_prime_mismatch: bool = False,
    na_conc: float = 0.05,         # M  (50 mM default)
    mg_conc: float = 0.0,          # M  (0 = Na⁺ only)
    dntp_conc: float = 0.0,        # M
    primer_conc: float = 250e-9,   # M  (250 nM default)
) -> float
```

---

## References

1. SantaLucia, J. (1998). A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. *PNAS*, 95(4), 1460–1465. https://doi.org/10.1073/pnas.95.4.1460

2. Allawi, H. T. & SantaLucia, J. (1997). Thermodynamics and NMR of internal G·T mismatches in DNA. *Biochemistry*, 36(34), 10581–10594. https://doi.org/10.1021/bi962590c

3. Peyret, N., Seneviratne, P. A., Allawi, H. T., & SantaLucia, J. (1999). Nearest-neighbor thermodynamics and NMR of DNA sequences with internal A·A, C·C, G·G, and T·T mismatches. *Biochemistry*, 38(12), 3468–3477. https://doi.org/10.1021/bi9825091

4. Owczarzy, R., Moreira, B. G., You, Y., Behlke, M. A., & Walder, J. A. (2008). Magnesium ions and DNA: oligonucleotide stability, thermodynamics, and implications for PCR performance. *Biochemistry*, 47(19), 5336–5353. https://doi.org/10.1021/bi702363u

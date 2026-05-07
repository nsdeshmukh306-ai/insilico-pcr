"""
Thermodynamics Module
=====================
Computes melting temperature (Tm) and approximate Gibbs free energy (ΔG)
for primer–template duplexes using nearest-neighbor (NN) thermodynamics.

Primary references
------------------
Matched NN parameters:
  SantaLucia J. Jr. (1998). A unified view of polymer, dumbbell, and
  oligonucleotide DNA nearest-neighbor thermodynamics.  PNAS 95:1460–1465.

Mismatch NN parameters (full context-dependent tables):
  G·T / T·G  — Allawi & SantaLucia (1997). Biochemistry 36:10581–10594.
  A·C / C·A  — Allawi & SantaLucia (1998). Biochemistry 37:2170–2179.
  C·T / T·C  — Allawi & SantaLucia (1998). Nucleic Acids Res 26:2694–2701.
  A·A / C·C / G·G / T·T / A·G / G·A
             — Peyret et al. (1999). Biochemistry 38:3468–3477.

Salt correction:
  Na⁺-only  — SantaLucia (1998): Tm_corr = Tm_1M + 16.6·log10([Na⁺])
  Mg²⁺/mixed — Owczarzy et al. (2008). Biochemistry 47:5336–5353.
               Regime selection by √[Mg²⁺]/[Na⁺] ratio;
               free [Mg²⁺] = max(0, [Mg²⁺] − [dNTP]) after chelation.

Assumptions / Limitations
--------------------------
  - DNA/DNA duplexes only (no RNA, no LNA, no modified bases).
  - Mismatch NN parameters cover 12 mismatch types × 4 5'-neighbour contexts.
    Each mismatch contributes its context-specific ΔH/ΔS to the NN step where
    it occupies the 3'-position; the flanking step (where it is at 5') uses
    the standard matched lookup (standard bioinformatics approximation).
  - A·G / G·A parameters from Peyret 1999; least precisely constrained —
    a context average fallback is used for unlisted entries.
  - Polymerase kinetics and processivity are NOT modelled.
"""

import math
from typing import Optional, Tuple

# ---------------------------------------------------------------------------
# SantaLucia 1998 Table 2: matched NN parameters for DNA/DNA
# Key: "XY/X'Y'" where XY is 5'→3' top (primer), X'Y' is 3'→5' bottom
# Value: (ΔH kcal/mol, ΔS cal/mol/K)
# ---------------------------------------------------------------------------
_NN_PARAMS: dict[str, Tuple[float, float]] = {
    "AA/TT": (-7.9, -22.2),
    "AT/TA": (-7.2, -20.4),
    "TA/AT": (-7.2, -21.3),
    "CA/GT": (-8.5, -22.7),
    "GT/CA": (-8.4, -22.4),
    "CT/GA": (-7.8, -21.0),
    "GA/CT": (-8.2, -22.2),
    "CG/GC": (-10.6, -27.2),
    "GC/CG": (-9.8, -24.4),
    "GG/CC": (-8.0, -19.9),
}

# Initiation parameters (SantaLucia 1998 Table 2)
_INIT_GC = (0.1,  -2.8)   # duplex terminating with G or C
_INIT_AT = (2.3,   4.1)   # duplex terminating with A or T

# Gas constant in cal/(mol·K)
_R = 1.987

# ---------------------------------------------------------------------------
# Peyret / Allawi & SantaLucia mismatch NN parameters
# ---------------------------------------------------------------------------
# Key: "XY/WZ" where
#   XY = 5'→3' primer dinucleotide
#   WZ = 3'→5' template dinucleotide at this step
#   X:W  always Watson-Crick  (5' base pair is matched)
#   Y:Z  is the mismatch      (Z ≠ Watson-Crick complement of Y)
# Value: (ΔH kcal/mol, ΔS cal/mol/K)
# ---------------------------------------------------------------------------
_MM_NN_PARAMS: dict[str, Tuple[float, float]] = {

    # --- G·T mismatches (primer G, template T) — Allawi & SantaLucia 1997 ---
    "AG/TT": ( 1.0,   0.9),
    "CG/GT": (-4.1, -11.7),
    "GG/CT": (-3.8, -12.2),
    "TG/AT": ( 5.2,  13.5),

    # --- T·G mismatches (primer T, template G) — Allawi & SantaLucia 1997 ---
    "AT/TG": (-2.5,  -8.3),
    "CT/GG": (-0.1,  -1.0),
    "GT/CG": (-1.5,  -2.7),
    "TT/AG": ( 3.9,   9.4),

    # --- A·C mismatches (primer A, template C) — Allawi & SantaLucia 1998a ---
    "AA/TC": ( 2.3,   4.6),
    "CA/GC": ( 0.6,  -0.6),
    "GA/CC": (-0.7,  -3.8),
    "TA/AC": ( 3.4,   8.0),

    # --- C·A mismatches (primer C, template A) — Allawi & SantaLucia 1998a ---
    "AC/TA": ( 5.3,  14.6),
    "CC/GA": ( 1.9,   3.7),
    "GC/CA": ( 5.2,  14.2),
    "TC/AA": ( 7.6,  20.2),

    # --- C·T mismatches (primer C, template T) — Allawi & SantaLucia 1998b ---
    "AC/TT": ( 0.7,   0.2),
    "CC/GT": (-0.8,  -4.5),
    "GC/CT": (-3.4, -10.8),
    "TC/AT": ( 1.2,   0.7),

    # --- T·C mismatches (primer T, template C) — Allawi & SantaLucia 1998b ---
    "AT/TC": (-1.2,  -6.2),
    "CT/GC": ( 2.3,   5.4),
    "GT/CC": ( 5.2,  13.5),
    "TT/AC": ( 1.0,   0.7),

    # --- A·A mismatches (primer A, template A) — Peyret et al. 1999 ---
    "AA/TA": ( 1.2,   1.7),
    "CA/GA": (-0.9,  -4.2),
    "GA/CA": (-2.9,  -9.8),
    "TA/AA": ( 4.7,  12.9),

    # --- C·C mismatches (primer C, template C) — Peyret et al. 1999 ---
    "AC/TC": ( 0.0,  -4.4),
    "CC/GC": (-1.5,  -7.2),
    "GC/CC": ( 3.6,  -0.1),
    "TC/AC": ( 6.1,  16.4),

    # --- G·G mismatches (primer G, template G) — Peyret et al. 1999 ---
    "AG/TG": (-3.1,  -9.5),
    "CG/GG": (-4.9, -15.3),
    "GG/CG": (-6.0, -15.8),
    "TG/AG": ( 1.6,   3.6),

    # --- T·T mismatches (primer T, template T) — Peyret et al. 1999 ---
    "AT/TT": (-2.7, -10.8),
    "CT/GT": (-5.0, -15.8),
    "GT/CT": (-2.2,  -8.4),
    "TT/AT": ( 0.2,  -1.5),

    # --- A·G mismatches (primer A, template G) — Peyret et al. 1999 ---
    "AA/TG": (-0.6,  -2.3),
    "CA/GG": (-3.1,  -9.5),
    "GA/CG": (-4.3, -10.7),
    "TA/AG": ( 4.7,  14.2),

    # --- G·A mismatches (primer G, template A) — Peyret et al. 1999 ---
    "AG/TA": (-0.7,  -2.3),
    "CG/GA": (-2.9,  -9.8),
    "GG/CA": (-0.6,  -0.7),
    "TG/AA": ( 1.6,   3.6),
}

# Owczarzy 2008 Mg²⁺-correction constants (Table 2, equation 16)
# Units: all in K⁻¹
_OWC_A = 3.92e-5
_OWC_B = -9.11e-6
_OWC_C = 6.26e-5
_OWC_D = 1.42e-5
_OWC_E = -4.82e-4
_OWC_F =  5.25e-4
_OWC_G =  8.31e-5


# ---------------------------------------------------------------------------
# Helper: complement / reverse-complement
# ---------------------------------------------------------------------------
_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)


def reverse_complement(seq: str) -> str:
    return complement(seq)[::-1]


# ---------------------------------------------------------------------------
# Matched NN lookup
# ---------------------------------------------------------------------------
def _nn_lookup(dinuc: str) -> Tuple[float, float]:
    """Return (ΔH, ΔS) for a 5'→3' dinucleotide from the matched NN table."""
    comp = complement(dinuc)          # bottom strand 3'→5'
    key  = f"{dinuc}/{comp}"
    if key in _NN_PARAMS:
        return _NN_PARAMS[key]
    sym = f"{comp[::-1]}/{dinuc[::-1]}"
    if sym in _NN_PARAMS:
        return _NN_PARAMS[sym]
    avg_dh = sum(v[0] for v in _NN_PARAMS.values()) / len(_NN_PARAMS)
    avg_ds = sum(v[1] for v in _NN_PARAMS.values()) / len(_NN_PARAMS)
    return (avg_dh, avg_ds)


# ---------------------------------------------------------------------------
# Mismatch NN lookup (Peyret / Allawi & SantaLucia)
# ---------------------------------------------------------------------------
def _mm_nn_lookup(primer_dinuc: str, template_3prime_base: str) -> Tuple[float, float]:
    """
    Return (ΔH, ΔS) for a dinucleotide step where the 3'-position is mismatched.

    primer_dinuc          : 2-char 5'→3' primer dinucleotide  (e.g. "AG")
    template_3prime_base  : template base opposite the 3' primer base (e.g. "T")
                            This is read in the 3'→5' direction of the template,
                            i.e. the raw template character opposite primer[i+1].
    """
    X = primer_dinuc[0]
    Y = primer_dinuc[1]
    W = complement(X)                 # 5' base is always matched
    Z = template_3prime_base.upper()
    key = f"{X}{Y}/{W}{Z}"
    if key in _MM_NN_PARAMS:
        return _MM_NN_PARAMS[key]

    # Symmetry: same mismatch read from the other strand
    sym = f"{complement(Z)}{complement(W)}/{complement(Y)}{complement(X)}"
    if sym in _MM_NN_PARAMS:
        return _MM_NN_PARAMS[sym]

    # Context-average fallback: average all entries with same mismatch type (Y·Z)
    matching = [v for k, v in _MM_NN_PARAMS.items()
                if len(k) == 5 and k[1] == Y and k[4] == Z]
    if matching:
        return (sum(v[0] for v in matching) / len(matching),
                sum(v[1] for v in matching) / len(matching))

    # Last resort: global average of all mismatch params
    all_vals = list(_MM_NN_PARAMS.values())
    return (sum(v[0] for v in all_vals) / len(all_vals),
            sum(v[1] for v in all_vals) / len(all_vals))


# ---------------------------------------------------------------------------
# Core NN thermodynamics
# ---------------------------------------------------------------------------
def calc_nn_thermodynamics(
    primer: str,
    template: Optional[str] = None,
    mismatch_positions: Optional[list] = None,
    three_prime_mismatch: bool = False,
) -> Tuple[float, float]:
    """
    Compute total ΔH (kcal/mol) and ΔS (cal/mol/K) by nearest-neighbor summing.

    When `template` is provided (aligned, same length as primer, no gap chars),
    the full Peyret/Allawi mismatch NN tables are used for every mismatched
    dinucleotide step — the most accurate mode.

    When only `mismatch_positions` is given (backward-compatible API), a
    context-averaged mismatch NN value is used for affected steps.

    When neither is given, the primer is evaluated against its perfect complement
    (pure matched NN — correct for Tm without mismatches).

    Parameters
    ----------
    primer : str
        5'→3' primer sequence (uppercase DNA, ACGT only).
    template : str, optional
        Aligned template sequence (5'→3'), same length as primer.
        Mismatches are detected by comparison with complement(primer[i]).
    mismatch_positions : list of int, optional
        0-based positions within the primer where mismatches occur.
        Only used when `template` is not provided.
    three_prime_mismatch : bool
        If True, adds an extra +2.5 kcal/mol initiation-like ΔH penalty
        for the 3'-terminal mismatch (blocks polymerase extension).

    Returns
    -------
    (dH_total, dS_total) : (float, float)
    """
    primer = primer.upper().replace("U", "T")
    if len(primer) < 2:
        return (0.0, 0.0)

    if template is not None:
        template = template.upper().replace("U", "T")
        if len(template) < len(primer):
            template = template + "N" * (len(primer) - len(template))
        template = template[: len(primer)]

    mismatch_set = set(mismatch_positions) if mismatch_positions else set()

    # Initiation penalties (both termini)
    dH, dS = 0.0, 0.0
    for terminal in (primer[0], primer[-1]):
        if terminal in ("G", "C"):
            dH += _INIT_GC[0];  dS += _INIT_GC[1]
        else:
            dH += _INIT_AT[0];  dS += _INIT_AT[1]

    # Nearest-neighbour summation
    for i in range(len(primer) - 1):
        dinuc = primer[i: i + 2]

        if template is not None:
            # Full Peyret/Allawi lookup
            p3 = primer[i + 1]
            t3 = template[i + 1]   # template base opposite primer[i+1], read 5'→3'
            # Template is 5'→3'; to pair with primer at position i+1, we need
            # complement(t3) == p3 for a match.  The 3'→5' template base here
            # is complement(t3), so match condition: complement(p3) == t3.
            if t3 == "N" or complement(p3) == t3:
                nn_dh, nn_ds = _nn_lookup(dinuc)
            else:
                # Mismatch at 3' position of this step.
                # template is stored in "parallel complement" convention:
                # template[i] = complement(primer[i]) for a perfect match,
                # i.e. each element is the 3'→5' antiparallel base at position i.
                # The NN table keys are "XY/WZ" where Z is the 3'→5' template
                # base at the 3' position of the step — so Z = t3 directly.
                nn_dh, nn_ds = _mm_nn_lookup(dinuc, t3)

        elif mismatch_set:
            # Backward-compatible: position-only mismatch info
            if (i + 1) in mismatch_set:
                # Mismatch at 3' position of this step; use context-averaged params
                X = dinuc[0]
                context_vals = [v for k, v in _MM_NN_PARAMS.items()
                                if k[0] == X]
                if context_vals:
                    nn_dh = sum(v[0] for v in context_vals) / len(context_vals)
                    nn_ds = sum(v[1] for v in context_vals) / len(context_vals)
                else:
                    nn_dh, nn_ds = _nn_lookup(dinuc)
            else:
                nn_dh, nn_ds = _nn_lookup(dinuc)
        else:
            # Perfect duplex
            nn_dh, nn_ds = _nn_lookup(dinuc)

        dH += nn_dh
        dS += nn_ds

    # Additional destabilisation for 3'-terminal mismatch
    # (blocks polymerase; disrupts terminal stacking beyond NN contribution)
    if three_prime_mismatch:
        dH += 2.5    # kcal/mol, positive → destabilising

    return (dH, dS)


# ---------------------------------------------------------------------------
# Salt correction  — Owczarzy et al. (2008)
# ---------------------------------------------------------------------------
def _owczarzy_salt_correction(
    tm_1m_celsius: float,
    fGC: float,
    n_bases: int,
    na_conc: float,
    mg_conc: float,
    dntp_conc: float,
) -> float:
    """
    Apply salt correction per Owczarzy et al. (2008) Biochemistry 47:5336.

    Free Mg²⁺ is reduced by dNTP chelation (1:1 stoichiometry).
    Regime is selected by the √[free Mg²⁺] / [Na⁺] ratio:
      ratio < 0.22  → Na⁺-only correction (SantaLucia 1998)
      0.22 ≤ ratio < 6.0 → mixed (linear interpolation between Na and Mg formulae)
      ratio ≥ 6.0   → Mg²⁺-only correction (Owczarzy eq. 16)

    Parameters
    ----------
    tm_1m_celsius : float   Tm at 1 M NaCl in °C
    fGC           : float   GC fraction (0–1)
    n_bases       : int     Primer length
    na_conc       : float   [Na⁺] in mol/L
    mg_conc       : float   [Mg²⁺] (total) in mol/L
    dntp_conc     : float   [dNTP] in mol/L (chelates Mg²⁺ 1:1)

    Returns
    -------
    Tm corrected (°C)
    """
    free_mg = max(0.0, mg_conc - dntp_conc)

    def _mg_formula(tm1m_c: float) -> float:
        """Owczarzy 2008 eq. 16 — Mg²⁺-only."""
        if free_mg <= 0:
            return tm1m_c
        ln_mg = math.log(free_mg)
        inv = (1.0 / (tm1m_c + 273.15)
               + _OWC_A + _OWC_B * ln_mg
               + fGC * (_OWC_C + _OWC_D * ln_mg)
               + (1.0 / (2.0 * (n_bases - 1)))
               * (_OWC_E + _OWC_F * ln_mg + _OWC_G * ln_mg ** 2))
        return (1.0 / inv) - 273.15

    def _na_formula(tm1m_c: float) -> float:
        """SantaLucia 1998 / Wetmur 1991 — Na⁺-only."""
        if na_conc <= 0:
            return tm1m_c
        return tm1m_c + 16.6 * math.log10(na_conc)

    if free_mg == 0.0:
        # No free Mg²⁺ — Na⁺-only
        return _na_formula(tm_1m_celsius)

    if na_conc <= 0.0:
        # No Na⁺ — Mg²⁺-only
        return _mg_formula(tm_1m_celsius)

    ratio = math.sqrt(free_mg) / na_conc

    if ratio < 0.22:
        return _na_formula(tm_1m_celsius)
    elif ratio >= 6.0:
        return _mg_formula(tm_1m_celsius)
    else:
        # Mixed regime: linear blend by ratio position in [0.22, 6.0)
        blend = (ratio - 0.22) / (6.0 - 0.22)
        tm_na = _na_formula(tm_1m_celsius)
        tm_mg = _mg_formula(tm_1m_celsius)
        return tm_na + blend * (tm_mg - tm_na)


# ---------------------------------------------------------------------------
# GC content
# ---------------------------------------------------------------------------
def gc_content(seq: str) -> float:
    """Return GC fraction (0.0–1.0) for a DNA sequence."""
    seq = seq.upper()
    if not seq:
        return 0.0
    return sum(1 for b in seq if b in ("G", "C")) / len(seq)


# ---------------------------------------------------------------------------
# Melting temperature
# ---------------------------------------------------------------------------
def calc_tm(
    primer: str,
    template: Optional[str] = None,
    mismatch_positions: Optional[list] = None,
    three_prime_mismatch: bool = False,
    primer_conc: float = 250e-9,   # 250 nM, typical PCR primer concentration
    na_conc: float    = 0.05,      # 50 mM Na⁺, typical PCR buffer
    mg_conc: float    = 0.0,       # Mg²⁺ (mol/L); 0 → Na⁺-only correction
    dntp_conc: float  = 0.0,       # dNTP (mol/L); chelates Mg²⁺ 1:1
) -> float:
    """
    Calculate melting temperature (°C) using nearest-neighbor thermodynamics
    with Owczarzy et al. (2008) salt correction.

    Parameters
    ----------
    primer : str
        5'→3' primer sequence.
    template : str, optional
        Aligned template (5'→3', same length as primer).
        Enables full Peyret/Allawi mismatch NN lookup.
    mismatch_positions : list of int, optional
        0-based mismatch positions (backward-compatible; used when template
        is not provided).
    three_prime_mismatch : bool
        3'-terminal mismatch flag (extra destabilisation).
    primer_conc : float
        Total strand concentration in mol/L (default 250 nM).
    na_conc : float
        Monovalent Na⁺ concentration in mol/L (default 50 mM).
    mg_conc : float
        Total Mg²⁺ concentration in mol/L (default 0 — Na-only mode).
    dntp_conc : float
        Total dNTP concentration in mol/L (default 0).
        Free Mg²⁺ = max(0, mg_conc − dntp_conc).

    Returns
    -------
    Tm : float  (°C)
    """
    dH, dS = calc_nn_thermodynamics(
        primer, template, mismatch_positions, three_prime_mismatch
    )
    if dH == 0.0:
        return 0.0

    c_t = primer_conc / 4.0          # non-self-complementary strand conc
    dH_cal = dH * 1000.0             # kcal/mol → cal/mol
    tm_1m = dH_cal / (dS + _R * math.log(c_t)) - 273.15   # °C at 1 M NaCl

    fGC = gc_content(primer)
    n   = len(primer)
    return round(_owczarzy_salt_correction(tm_1m, fGC, n, na_conc, mg_conc, dntp_conc), 2)


# ---------------------------------------------------------------------------
# Gibbs free energy
# ---------------------------------------------------------------------------
def calc_delta_g(
    primer: str,
    template: Optional[str] = None,
    mismatch_positions: Optional[list] = None,
    three_prime_mismatch: bool = False,
    temperature: float = 37.0,   # °C — physiological / PCR annealing temperature
) -> float:
    """
    Calculate ΔG (kcal/mol) at a given temperature.

    ΔG = ΔH − T·ΔS   (ΔH in kcal/mol, ΔS in cal/mol/K → divide by 1000)

    Parameters
    ----------
    primer : str
    template : str, optional
    mismatch_positions : list of int, optional
    three_prime_mismatch : bool
    temperature : float  (°C, default 37)

    Returns
    -------
    ΔG : float (kcal/mol)
    """
    dH, dS = calc_nn_thermodynamics(
        primer, template, mismatch_positions, three_prime_mismatch
    )
    T_k = temperature + 273.15
    return round(dH - T_k * (dS / 1000.0), 3)


# ---------------------------------------------------------------------------
# Simple Tm shortcut (Wallace rule — fast pre-filter, short oligos only)
# ---------------------------------------------------------------------------
def calc_tm_basic(seq: str) -> float:
    """Wallace rule: Tm = 2·(A+T) + 4·(G+C).  Valid only for <14 bp."""
    seq = seq.upper()
    at = sum(1 for b in seq if b in ("A", "T"))
    gc = sum(1 for b in seq if b in ("G", "C"))
    return float(2 * at + 4 * gc)


# ---------------------------------------------------------------------------
# Primer quality checks
# ---------------------------------------------------------------------------
def check_primer_quality(
    primer: str,
    min_tm: float = 50.0,
    max_tm: float = 72.0,
    min_gc: float = 0.40,
    max_gc: float = 0.65,
    na_conc: float = 0.05,
    mg_conc: float = 0.0,
    dntp_conc: float = 0.0,
) -> dict:
    """
    Return a dict of quality flags for a primer sequence.

    Checks performed:
      - Length (18–30 bp recommended)
      - GC content (40–65%)
      - Tm (50–72 °C, at specified salt conditions)
      - 3'-end GC clamp (last 5 bases should include 1–3 G/C)
      - Homopolymer runs (≥4 identical consecutive bases flagged)
    """
    primer = primer.upper()
    length = len(primer)
    gc     = gc_content(primer)
    tm     = calc_tm(primer, na_conc=na_conc, mg_conc=mg_conc, dntp_conc=dntp_conc)
    dg     = calc_delta_g(primer)

    tail     = primer[-5:] if len(primer) >= 5 else primer
    tail_gc  = sum(1 for b in tail if b in ("G", "C"))
    gc_clamp = 1 <= tail_gc <= 3

    max_run = run = 1
    for i in range(1, len(primer)):
        run = run + 1 if primer[i] == primer[i - 1] else 1
        if run > max_run:
            max_run = run

    return {
        "sequence":       primer,
        "length":         length,
        "gc_fraction":    round(gc, 3),
        "tm_celsius":     tm,
        "delta_g_37":     dg,
        "length_ok":      18 <= length <= 30,
        "gc_ok":          min_gc <= gc <= max_gc,
        "tm_ok":          min_tm <= tm <= max_tm,
        "gc_clamp_ok":    gc_clamp,
        "max_run":        max_run,
        "low_complexity": max_run >= 4,
    }

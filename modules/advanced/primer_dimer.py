"""
Primer-Dimer Detection — Thermodynamic ΔG Model
================================================
Detects homo- and hetero-dimer formation between primers by computing the
nearest-neighbour duplex ΔG for every possible alignment between primer A
and the reverse complement of primer B.

Method
------
For each alignment offset, the overlapping region forms a partial duplex.
The NN model (SantaLucia 1998 matched + Peyret/Allawi mismatch tables) is
used to compute:

  ΔG_duplex = ΔH − T · ΔS   at 37 °C (physiological / approximate annealing)

The alignment with the MOST NEGATIVE ΔG is reported as the primary dimer.

3'-end penalty
--------------
If the 3' end of either primer falls inside the overlapping region, an
additional ΔG contribution is applied:

  ΔG_3prime_adjusted = ΔG_duplex + Δ_3prime

where Δ_3prime = −PRIME3_EXTRA_DG kcal/mol (further stabilises the dimer,
making it more negative — worse for PCR specificity).  The adjusted value is
what drives the risk flag and is stored in dimer_score as abs(ΔG_adjusted).

Threshold
---------
A dimer is flagged when abs(ΔG_adjusted) ≥ DG_THRESHOLD_KCAL (default 3.5
kcal/mol, ≈ Tm ~ 25 °C for short partial duplexes, enough for extension risk).

This replaces the previous heuristic complementarity count (dimer_score ≥ 4).
The backward-compatible `dimer_score` field equals abs(ΔG_adjusted) so existing
code that checks `dimer_score >= 0` still works.
"""

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

# Thermodynamic engine
import sys, os as _os
sys.path.insert(0, _os.path.join(_os.path.dirname(__file__), "..", ".."))
from modules.thermodynamics import calc_nn_thermodynamics, complement as _compl

# Temperature for ΔG calculation (°C)
_TEMPERATURE = 37.0

# Extra ΔG stabilisation (kcal/mol, negative → more stable dimer) when the
# 3' end of a primer is in the duplex region (can be extended by polymerase).
PRIME3_EXTRA_DG: float = -1.5

# Flag threshold: abs(ΔG_adjusted) must exceed this to generate a warning.
DG_THRESHOLD_KCAL: float = 3.5

# Minimum overlap length to evaluate (too-short windows don't form stable duplexes).
MIN_OVERLAP: int = 4


_COMP_MAP = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}


def _rc(seq: str) -> str:
    return "".join(_COMP_MAP.get(b, "N") for b in reversed(seq.upper()))


def _aln_string(sub_a: str, sub_rc_b: str) -> str:
    """Return a 3-line ASCII alignment depiction."""
    top  = list(sub_a)
    bot  = list(sub_rc_b)
    conn = []
    for a, b in zip(sub_a, sub_rc_b):
        wc = _COMP_MAP.get(a) == b
        if wc:
            conn.append("|")
        elif (a == "G" and b == "T") or (a == "T" and b == "G"):
            conn.append(".")
        else:
            conn.append(" ")
    return (
        "5'-" + "".join(top)  + "-3'\n"
        "   " + "".join(conn) + "\n"
        "3'-" + "".join(bot)  + "-5'"
    )


def _duplex_dg(sub_a: str, sub_rc_b: str, temperature: float = _TEMPERATURE) -> float:
    """
    Compute ΔG (kcal/mol) of a partial duplex formed by sub_a (5'→3') pairing
    with sub_rc_b (5'→3', same reading direction).

    sub_rc_b is a portion of rc(primer_b), so each position is what the
    antiparallel primer_b presents as a template base to sub_a.
    """
    dH, dS = calc_nn_thermodynamics(
        primer   = sub_a,
        template = sub_rc_b,
    )
    T_k = temperature + 273.15
    return dH - T_k * (dS / 1000.0)


# ---------------------------------------------------------------------------
# Result dataclass
# ---------------------------------------------------------------------------
@dataclass
class DimerResult:
    primer_a:       str
    primer_b:       str
    dimer_type:     str     # "homo" or "hetero"
    dimer_score:    float   # abs(ΔG_adjusted) in kcal/mol; higher = worse risk
    delta_g:        float   # raw duplex ΔG (kcal/mol); more negative = more stable dimer
    delta_g_adj:    float   # ΔG with 3'-end penalty applied
    three_prime_a:  bool    # 3' end of primer_a inside duplex region
    three_prime_b:  bool    # 3' end of primer_b inside duplex region
    alignment:      str     # ASCII alignment depiction of best window
    warning:        str


# ---------------------------------------------------------------------------
# Core: best-ΔG alignment
# ---------------------------------------------------------------------------
def _best_dg_alignment(
    primer_a: str,
    rc_b:     str,
    temperature: float = _TEMPERATURE,
) -> Tuple[float, float, bool, bool, str]:
    """
    Slide rc_b across primer_a and find the alignment with the most negative
    ΔG (best dimer candidate).

    Returns
    -------
    (best_dg, best_dg_adj, three_prime_a, three_prime_b, alignment_str)
    """
    len_a = len(primer_a)
    len_b = len(rc_b)

    best_dg     = 0.0          # no dimer → ΔG = 0 (convention)
    best_dg_adj = 0.0
    best_3pa    = False
    best_3pb    = False
    best_aln    = ""

    for offset in range(-(len_b - 1), len_a):
        a_start = max(0, offset)
        a_end   = min(len_a, offset + len_b)
        b_start = max(0, -offset)
        b_end   = b_start + (a_end - a_start)

        if a_end - a_start < MIN_OVERLAP:
            continue

        sub_a    = primer_a[a_start: a_end]
        sub_rc_b = rc_b[b_start: b_end]

        # Quick pre-filter: at least 50% complementarity needed for ΔG calc
        wc_count = sum(
            1 for x, y in zip(sub_a, sub_rc_b) if _COMP_MAP.get(x) == y
        )
        if wc_count < MIN_OVERLAP // 2:
            continue

        dg = _duplex_dg(sub_a, sub_rc_b, temperature)

        # 3'-end involvement
        three_pa = (a_end >= len_a - 4)         # 3' end of A in window
        three_pb = (b_start <= 4)               # 3' end of B (5' of rc_b) in window

        dg_adj = dg + (PRIME3_EXTRA_DG if (three_pa or three_pb) else 0.0)

        if dg_adj < best_dg_adj:
            best_dg     = dg
            best_dg_adj = dg_adj
            best_3pa    = three_pa
            best_3pb    = three_pb
            best_aln    = _aln_string(sub_a, sub_rc_b)

    return best_dg, best_dg_adj, best_3pa, best_3pb, best_aln


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
def detect_primer_dimer(
    primer_a: str,
    primer_b: str,
    name_a: str = "FWD",
    name_b: str = "REV",
    dg_threshold: float = DG_THRESHOLD_KCAL,
    temperature: float  = _TEMPERATURE,
) -> DimerResult:
    """
    Check for dimer formation between two primers using NN ΔG.

    Parameters
    ----------
    primer_a, primer_b : str  5'→3' primer sequences (uppercase ACGT)
    name_a, name_b     : str  Labels for warning messages
    dg_threshold       : float  abs(ΔG_adj) must exceed this to flag (kcal/mol)
    temperature        : float  °C for ΔG calculation (default 37)

    Returns
    -------
    DimerResult
      - delta_g      : raw ΔG of best alignment (most negative = most stable)
      - delta_g_adj  : ΔG with 3'-end stabilisation penalty applied
      - dimer_score  : abs(delta_g_adj)   — 0 if no dimer, positive = dimer risk
      - warning      : non-empty string if risk ≥ threshold
    """
    primer_a = primer_a.upper()
    primer_b = primer_b.upper()
    dimer_type = "homo" if primer_a == primer_b else "hetero"

    rc_b = _rc(primer_b)
    dg, dg_adj, tp_a, tp_b, aln = _best_dg_alignment(primer_a, rc_b, temperature)

    score = abs(dg_adj)   # positive; 0 means no predicted dimer

    warning = ""
    if score >= dg_threshold:
        parts = []
        if tp_a:
            parts.append(f"3'-end of {name_a} involved")
        if tp_b:
            parts.append(f"3'-end of {name_b} involved")
        why = "; ".join(parts) if parts else "internal"
        warning = (
            f"{dimer_type.capitalize()} primer dimer risk "
            f"(ΔG={dg_adj:.2f} kcal/mol, {why}). "
            "May reduce PCR efficiency or generate artefact products."
        )

    return DimerResult(
        primer_a     = primer_a,
        primer_b     = primer_b,
        dimer_type   = dimer_type,
        dimer_score  = round(score, 3),
        delta_g      = round(dg, 3),
        delta_g_adj  = round(dg_adj, 3),
        three_prime_a = tp_a,
        three_prime_b = tp_b,
        alignment    = aln,
        warning      = warning,
    )


def check_all_dimers(
    forward: str,
    reverse: str,
    temperature: float = _TEMPERATURE,
) -> List[DimerResult]:
    """Check all three dimer combinations: FWD-FWD, REV-REV, FWD-REV."""
    return [
        detect_primer_dimer(forward, forward, "FWD", "FWD", temperature=temperature),
        detect_primer_dimer(reverse, reverse, "REV", "REV", temperature=temperature),
        detect_primer_dimer(forward, reverse, "FWD", "REV", temperature=temperature),
    ]

"""
Hairpin Prediction
===================
Detects potential hairpin (snap-back) structures within a primer sequence.

Method:
  A hairpin forms when the primer's 3' end folds back and base-pairs with
  an internal region of the same strand. We scan all possible stem–loop
  configurations:
    - Stem: 4–12 bp
    - Loop: 3–8 nt (thermodynamically required minimum loop = 3 nt)

Scoring:
  ΔG_hairpin ≈ ΔG_stem + ΔG_loop (approximate)
  Using Turner 2004 loop penalty table (approximated).

Limitation:
  Full RNA/DNA thermodynamic hairpin folding (e.g. mfold) is NOT used.
  This is a simplified stem–loop scan intended as a fast quality flag.
  Primers with ΔG_hairpin < -3 kcal/mol are flagged as problematic.

References:
  Turner 2004 approximation; SantaLucia & Hicks (2004) for stem ΔG.
"""

import math
from dataclasses import dataclass, field
from typing import List, Optional

# Approximate ΔG for loop closure (Turner 2004, interpolated, DNA)
# Key: loop length (nt), Value: ΔG penalty (kcal/mol, positive = destabilising)
_LOOP_DG = {
    3: 5.4,
    4: 4.9,
    5: 4.4,
    6: 4.1,
    7: 3.8,
    8: 3.6,
}
_LOOP_DG_DEFAULT = 3.4  # For loops > 8 nt

# Watson-Crick complement mapping
_WC = {"A": "T", "T": "A", "G": "C", "C": "G"}


def _wc_pair(a: str, b: str) -> bool:
    return _WC.get(a) == b


def _stem_dg(stem_seq: str) -> float:
    """
    Approximate ΔG for a DNA/DNA stem using simplified NN parameters.
    Uses the average NN ΔG of -1.5 kcal/mol per bp (rough approximation).
    More accurate: use thermodynamics.calc_delta_g on the stem.
    """
    from ..thermodynamics import calc_delta_g
    return calc_delta_g(stem_seq, temperature=37.0)


@dataclass
class HairpinResult:
    primer:        str
    has_hairpin:   bool
    delta_g:       float    # kcal/mol, most stable hairpin
    stem_start:    int      # 0-based position of stem start in primer
    stem_end:      int      # 0-based position of stem end
    loop_start:    int
    loop_end:      int
    stem_seq:      str
    loop_seq:      str
    warning:       str      # human-readable warning text


def detect_hairpin(
    primer: str,
    min_stem: int = 4,
    max_stem: int = 12,
    min_loop: int = 3,
    max_loop: int = 8,
    dg_threshold: float = -3.0,
) -> HairpinResult:
    """
    Scan a primer for hairpin structures.

    Parameters
    ----------
    primer : str  (uppercase ACGT)
    min_stem, max_stem : int  stem length range
    min_loop, max_loop : int  loop length range
    dg_threshold : float  (kcal/mol) flag if ΔG below this value

    Returns
    -------
    HairpinResult
    """
    primer = primer.upper()
    n = len(primer)
    best_dg = 0.0
    best_stem_start = 0
    best_stem_end   = 0
    best_loop_start = 0
    best_loop_end   = 0
    best_stem_seq   = ""
    best_loop_seq   = ""

    for stem_len in range(min_stem, max_stem + 1):
        for loop_len in range(min_loop, max_loop + 1):
            # Hairpin: stem1 (5' end) + loop + stem2 (3' end, complement of stem1)
            # Total length needed: stem_len + loop_len + stem_len
            total = 2 * stem_len + loop_len
            if total > n:
                continue

            # Scan all positions where this hairpin could form
            for i in range(n - total + 1):
                stem1 = primer[i : i + stem_len]
                loop  = primer[i + stem_len : i + stem_len + loop_len]
                stem2 = primer[i + stem_len + loop_len : i + total]

                # stem2 must be reverse-complement of stem1
                rc_stem1 = stem1[::-1]
                matches = sum(_wc_pair(rc_stem1[k], stem2[k]) for k in range(stem_len))
                if matches < stem_len - 1:   # allow one mismatch in stem
                    continue

                # Compute ΔG
                loop_dg = _LOOP_DG.get(loop_len, _LOOP_DG_DEFAULT)
                stem_dg = _stem_dg(stem1)   # use stem1 as the duplex

                dg = stem_dg + loop_dg

                if dg < best_dg:
                    best_dg         = dg
                    best_stem_start = i
                    best_stem_end   = i + stem_len
                    best_loop_start = i + stem_len
                    best_loop_end   = i + stem_len + loop_len
                    best_stem_seq   = stem1
                    best_loop_seq   = loop

    has_hairpin = best_dg < dg_threshold
    if has_hairpin:
        warning = (
            f"Potential hairpin detected (ΔG={best_dg:.2f} kcal/mol, "
            f"stem={best_stem_seq}, loop={best_loop_seq}). "
            "May reduce effective primer concentration."
        )
    else:
        warning = ""

    return HairpinResult(
        primer      = primer,
        has_hairpin = has_hairpin,
        delta_g     = round(best_dg, 3),
        stem_start  = best_stem_start,
        stem_end    = best_stem_end,
        loop_start  = best_loop_start,
        loop_end    = best_loop_end,
        stem_seq    = best_stem_seq,
        loop_seq    = best_loop_seq,
        warning     = warning,
    )

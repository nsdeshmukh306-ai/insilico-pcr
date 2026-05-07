"""
Scoring & Ranking Engine (Layer 9)
====================================
Computes a composite confidence score for each amplicon and ranks them.

Score formula
-------------
Score = w_bind * S_bind
      + w_tm   * S_tm
      + w_gc   * S_gc
      - w_mm   * P_mm
      - w_offt * P_offt
      + w_len  * S_len

Component definitions
~~~~~~~~~~~~~~~~~~~~~
S_bind  : mean normalised binding score of fwd+rev primers (0–1)
          = (fwd_binding_score + rev_binding_score) / 200

S_tm    : Tm compatibility score (0–1)
          Penalises Tm outside 50–72 °C and Tm delta >5 °C between primers
          S_tm = tm_range_factor * (1 - delta_tm_penalty)
          tm_range_factor: 1.0 if 55–68°C, declines to 0 outside 45–78°C
          delta_tm_penalty: min(1, |fwd_Tm - rev_Tm| / 20)

S_gc    : amplicon GC score (0–1)
          1.0 if 40–60% GC; declines outside

P_mm    : thermodynamic mismatch penalty (0–1) — replaces old count-based formula.

          Previous (removed): P_mm = (fwd_mm + rev_mm) / (2 * max_mismatches + 1)
          This was purely count-based and gave equal weight to all mismatches
          regardless of their thermodynamic impact.

          Current (ΔTm-driven):
            ΔTm_fwd = max(0, Tm_perfect_fwd − Tm_observed_fwd)
            ΔTm_rev = max(0, Tm_perfect_rev − Tm_observed_rev)
            P_mm = min(1, (ΔTm_fwd + ΔTm_rev) / DeltaTm_MAX)

          where DeltaTm_MAX = 2 * MISMATCH_TM_SCALE (default 20 °C) is the
          combined ΔTm at which the penalty saturates to 1.0.

          Rationale (Peyret/SantaLucia):
          - Each mismatch drops Tm by ~2–8 °C depending on type and context.
          - A G·T wobble at an internal position is far less destabilising than
            an A·C mismatch at the 3'-terminal position.
          - By using actual Tm drop (computed from Peyret NN tables via
            calc_tm(primer, template=aligned_template)), P_mm directly reflects
            the duplex stability loss rather than a raw count.
          - 3'-terminal mismatches are already penalised via three_prime_mm flag
            in calc_tm (extra +2.5 kcal/mol ΔH), so no separate additive boost
            is needed here.

          MISMATCH_TM_SCALE = 10 °C
          At ΔTm_combined ≥ 20 °C the penalty saturates at 1.0.
          A single 5 °C Tm drop gives P_mm = 0.25 (vs 0.33 for one mismatch
          in the old formula with max_mismatches=3).

P_offt  : off-target penalty (0–1)
          = log(1 + off_target_count) / log(1 + 100)

S_len   : amplicon length suitability (0–1) from amplicon.amplicon_length_score()

Default weights (sum to 1 for binding/scores, penalties are subtractive)
  w_bind = 0.30
  w_tm   = 0.25
  w_gc   = 0.10
  w_mm   = 0.20   (penalty weight)
  w_offt = 0.10   (penalty weight)
  w_len  = 0.05

Final score is clipped to [0, 100].

Design rationale
~~~~~~~~~~~~~~~~
- Binding score dominates (30%): reflects primer–template affinity.
- Tm compatibility (25%): mismatched Tm between primers reduces efficiency.
- Mismatch penalty (20%): thermodynamically calibrated — Peyret NN ΔTm shift.
- Off-target penalty (10%): penalises non-specific amplification.
- GC and length are supporting factors (15% combined).
"""

import math
import logging
from dataclasses import dataclass
from typing import List, Dict, Optional

from .amplicon import Amplicon, amplicon_length_score
from .thermodynamics import calc_tm, gc_content

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# ΔTm-based mismatch penalty calibration
# ---------------------------------------------------------------------------
# MISMATCH_TM_SCALE: a combined (fwd + rev) ΔTm of 2× this value saturates P_mm at 1.0.
# Empirically, a single 3-mismatch primer (max_mismatches=3) typically drops
# Tm by ~6–15 °C. Setting scale=10 °C means P_mm saturates at 20 °C combined drop.
MISMATCH_TM_SCALE: float = 10.0   # °C per primer; combined scale = 2 × 10 = 20 °C


def _perfect_tm(primer_seq: str, na_conc: float = 0.05, mg_conc: float = 0.0,
                dntp_conc: float = 0.0, primer_conc: float = 250e-9) -> float:
    """
    Tm of the primer against its perfect complement (no mismatches).
    Used as the baseline for computing ΔTm.
    """
    return calc_tm(primer_seq, na_conc=na_conc, mg_conc=mg_conc,
                   dntp_conc=dntp_conc, primer_conc=primer_conc)


def _delta_tm_penalty(
    primer_seq: str,
    aligned_template: Optional[str],
    mismatch_positions: list,
    three_prime_mm: bool,
    na_conc: float = 0.05,
    mg_conc: float = 0.0,
    dntp_conc: float = 0.0,
    primer_conc: float = 250e-9,
) -> float:
    """
    Compute ΔTm = Tm_perfect − Tm_observed for one primer.

    Uses the Peyret/Allawi mismatch NN tables (via calc_tm with template)
    when the aligned template is available; falls back to position-only
    mismatch_positions when not.

    Returns ΔTm in °C, clamped to [0, ∞).
    """
    tm_perfect = _perfect_tm(primer_seq, na_conc, mg_conc, dntp_conc, primer_conc)
    tm_observed = calc_tm(
        primer_seq,
        template             = aligned_template if aligned_template else None,
        mismatch_positions   = mismatch_positions if not aligned_template else None,
        three_prime_mismatch = three_prime_mm,
        na_conc              = na_conc,
        mg_conc              = mg_conc,
        dntp_conc            = dntp_conc,
        primer_conc          = primer_conc,
    )
    return max(0.0, tm_perfect - tm_observed)


# ---------------------------------------------------------------------------
# Default weights
# ---------------------------------------------------------------------------
DEFAULT_WEIGHTS = {
    "w_bind": 0.30,
    "w_tm":   0.25,
    "w_gc":   0.10,
    "w_mm":   0.20,
    "w_offt": 0.10,
    "w_len":  0.05,
}


@dataclass
class ScoredAmplicon:
    """An Amplicon annotated with individual score components and a final score."""
    amplicon:        Amplicon
    s_bind:          float    # binding score component (0–1)
    s_tm:            float    # Tm compatibility component (0–1)
    s_gc:            float    # GC content component (0–1)
    p_mm:            float    # mismatch penalty (0–1, higher = more mismatches)
    p_offt:          float    # off-target penalty (0–1)
    s_len:           float    # length suitability (0–1)
    final_score:     float    # composite score (0–100)
    off_target_count: int     # number of off-target amplicons for this pair
    is_intended:     bool     # True if this is the primary target amplicon


# ---------------------------------------------------------------------------
# Component calculators
# ---------------------------------------------------------------------------
def _tm_range_factor(tm: float) -> float:
    """Score Tm for being in the optimal range [55, 68] °C."""
    if 55.0 <= tm <= 68.0:
        return 1.0
    elif 45.0 <= tm < 55.0:
        return (tm - 45.0) / 10.0
    elif 68.0 < tm <= 78.0:
        return (78.0 - tm) / 10.0
    return 0.0


def _gc_score(gc: float) -> float:
    """Score GC content: peak 0.40–0.60, zero below 0.20 or above 0.80."""
    if 0.40 <= gc <= 0.60:
        return 1.0
    elif 0.20 <= gc < 0.40:
        return (gc - 0.20) / 0.20
    elif 0.60 < gc <= 0.80:
        return (0.80 - gc) / 0.20
    return 0.0


def score_amplicon(
    amp: Amplicon,
    off_target_count: int = 0,
    is_intended: bool = True,
    max_mismatches: int = 3,
    weights: Dict[str, float] = None,
    na_conc: float = 0.05,
    mg_conc: float = 0.0,
    dntp_conc: float = 0.0,
    primer_conc: float = 250e-9,
) -> ScoredAmplicon:
    """
    Compute all score components and final composite score for one amplicon.

    Parameters
    ----------
    amp : Amplicon
    off_target_count : int
        Number of OTHER amplicons produced by the same primer pair
        (not counting this one). Used to penalise non-specific amplification.
    is_intended : bool
        Mark this as the target amplicon (True) or an off-target (False).
    max_mismatches : int
        Pipeline max_mismatches setting (kept for API compatibility; P_mm is
        now driven by ΔTm, not raw count).
    weights : dict, optional
        Override default score weights.
    na_conc, mg_conc, dntp_conc, primer_conc : float
        Salt / concentration conditions forwarded to ΔTm calculation.
    """
    w = {**DEFAULT_WEIGHTS, **(weights or {})}

    fwd_site = amp.hit.fwd_site.site
    rev_site = amp.hit.rev_site.site

    # S_bind: average normalised binding score
    s_bind = (amp.fwd_binding_score + amp.rev_binding_score) / 200.0
    s_bind = max(0.0, min(1.0, s_bind))

    # S_tm: Tm compatibility
    fwd_tm_factor = _tm_range_factor(amp.fwd_tm)
    rev_tm_factor = _tm_range_factor(amp.rev_tm)
    avg_tm_factor = (fwd_tm_factor + rev_tm_factor) / 2.0
    delta_tm = abs(amp.fwd_tm - amp.rev_tm)
    delta_tm_penalty = min(1.0, delta_tm / 20.0)
    s_tm = avg_tm_factor * (1.0 - delta_tm_penalty)

    # S_gc: amplicon GC content
    s_gc = _gc_score(amp.gc_fraction)

    # P_mm: thermodynamic mismatch penalty (ΔTm-driven, Peyret NN tables)
    #
    # For each primer, compute ΔTm = Tm_perfect − Tm_observed.
    # Tm_observed uses the Peyret/Allawi mismatch NN tables when the aligned
    # template sequence is available (via calc_tm(template=...)).
    # Combined penalty saturates at 1.0 when ΔTm_total ≥ 2×MISMATCH_TM_SCALE.
    delta_tm_fwd = _delta_tm_penalty(
        primer_seq        = fwd_site.primer_seq,
        aligned_template  = fwd_site.aligned_template or None,
        mismatch_positions= fwd_site.mismatch_pos,
        three_prime_mm    = fwd_site.three_prime_mm,
        na_conc           = na_conc,
        mg_conc           = mg_conc,
        dntp_conc         = dntp_conc,
        primer_conc       = primer_conc,
    )
    delta_tm_rev = _delta_tm_penalty(
        primer_seq        = rev_site.primer_seq,
        aligned_template  = rev_site.aligned_template or None,
        mismatch_positions= rev_site.mismatch_pos,
        three_prime_mm    = rev_site.three_prime_mm,
        na_conc           = na_conc,
        mg_conc           = mg_conc,
        dntp_conc         = dntp_conc,
        primer_conc       = primer_conc,
    )
    combined_delta_tm = delta_tm_fwd + delta_tm_rev
    p_mm = min(1.0, combined_delta_tm / (2.0 * MISMATCH_TM_SCALE))

    # P_offt: off-target penalty
    p_offt = math.log1p(off_target_count) / math.log1p(100)
    p_offt = min(1.0, p_offt)

    # S_len: amplicon length suitability
    s_len = amplicon_length_score(amp.length)

    # Composite score (0–100)
    positive = (
        w["w_bind"] * s_bind +
        w["w_tm"]   * s_tm   +
        w["w_gc"]   * s_gc   +
        w["w_len"]  * s_len
    )
    penalties = (
        w["w_mm"]   * p_mm  +
        w["w_offt"] * p_offt
    )
    raw_score = (positive - penalties) / (
        w["w_bind"] + w["w_tm"] + w["w_gc"] + w["w_len"]
    )
    final_score = max(0.0, min(100.0, raw_score * 100.0))

    return ScoredAmplicon(
        amplicon         = amp,
        s_bind           = round(s_bind,  4),
        s_tm             = round(s_tm,    4),
        s_gc             = round(s_gc,    4),
        p_mm             = round(p_mm,    4),
        p_offt           = round(p_offt,  4),
        s_len            = round(s_len,   4),
        final_score      = round(final_score, 2),
        off_target_count = off_target_count,
        is_intended      = is_intended,
    )


def rank_amplicons(
    amplicons: List[Amplicon],
    intended_region: tuple = None,
    max_mismatches: int = 3,
    weights: Dict[str, float] = None,
) -> List[ScoredAmplicon]:
    """
    Score and rank all amplicons for a primer pair.

    Parameters
    ----------
    amplicons : list of Amplicon
    intended_region : (seq_id, start, end) tuple, optional
        If provided, the amplicon closest to this region is flagged as intended.
    max_mismatches : int
    weights : dict, optional

    Returns
    -------
    List of ScoredAmplicon sorted by final_score descending.
    """
    if not amplicons:
        return []

    total = len(amplicons)

    # Determine which amplicon is "intended" (best binding, no off-target context yet)
    # If intended_region given, use overlap; otherwise pick top binding score
    intended_idx = 0
    if intended_region is not None:
        i_seq, i_start, i_end = intended_region
        best_overlap = -1
        for idx, amp in enumerate(amplicons):
            if amp.seq_id != i_seq:
                continue
            overlap = max(0, min(amp.end, i_end) - max(amp.start, i_start))
            if overlap > best_overlap:
                best_overlap = overlap
                intended_idx = idx
    else:
        best_bs = -1
        for idx, amp in enumerate(amplicons):
            bs = amp.fwd_binding_score + amp.rev_binding_score
            if bs > best_bs:
                best_bs = bs
                intended_idx = idx

    scored = []
    for idx, amp in enumerate(amplicons):
        off_count = total - 1  # all others count as off-target for scoring purposes
        is_intended = (idx == intended_idx)
        sa = score_amplicon(
            amp,
            off_target_count = off_count if not is_intended else 0,
            is_intended      = is_intended,
            max_mismatches   = max_mismatches,
            weights          = weights,
        )
        scored.append(sa)

    scored.sort(key=lambda x: x.final_score, reverse=True)
    log.debug("Ranked %d amplicons; top score=%.1f", len(scored), scored[0].final_score if scored else 0)
    return scored

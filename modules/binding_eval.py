"""
Binding Evaluation Module (Layer 5)
=====================================
Computes full thermodynamic and structural properties for each BindingSite:
  - Mismatch count and positions (already in BindingSite)
  - 3'-end mismatch penalty weight
  - GC content of primer and of the binding window
  - Melting temperature (Tm) accounting for mismatches
  - ΔG at annealing temperature
  - Binding score (0–100, higher = better binding)

Binding score rationale:
  binding_score = 100 * (sw_score / perfect_sw_score)
               * (1 - 3prime_mm_penalty)
               * gc_factor
  where:
    perfect_sw_score = primer_len * SW_MATCH
    3prime_mm_penalty = 0.5 if 3'-end mismatch else 0
    gc_factor = 1 if 40%≤GC≤65%, otherwise scales down

This score is NOT a probability; it is a normalised heuristic for ranking.
"""

from dataclasses import dataclass
from typing import List

from .binding_search import BindingSite, SW_MATCH
from .thermodynamics import (
    calc_delta_g,
    calc_tm,
    gc_content,
)


@dataclass
class EvaluatedSite:
    """BindingSite enriched with thermodynamic properties."""
    site:          BindingSite
    tm:            float   # Melting temperature (°C), mismatch-adjusted
    delta_g:       float   # ΔG at 37 °C (kcal/mol)
    gc_primer:     float   # GC fraction of primer sequence
    gc_template:   float   # GC fraction of bound template window
    binding_score: float   # 0–100 composite binding quality score
    mm_weight:     float   # Weighted mismatch penalty (3' end = 3×)


def _mismatch_weight(mismatch_pos: List[int], primer_len: int) -> float:
    """
    Compute a weighted mismatch penalty.
    3'-end positions (last 5 bases) are penalised 3× more than internal.
    Position weights decay linearly from 1.0 (5' end) to 3.0 (3' end)
    over the last 5 bases.

    Returns normalised weight in [0, 1]: 0 = no mismatches, 1 = completely mismatched.
    """
    if not mismatch_pos or primer_len == 0:
        return 0.0

    total_weight = 0.0
    for p in mismatch_pos:
        distance_from_3prime = primer_len - 1 - p
        if distance_from_3prime < 5:
            w = 3.0 - 0.4 * distance_from_3prime   # 3.0 at 3' end, 1.4 at 5th from end
        else:
            w = 1.0
        total_weight += w

    # Normalise by maximum possible weight (all bases mismatched)
    max_weight = sum(
        3.0 - 0.4 * (primer_len - 1 - i) if (primer_len - 1 - i) < 5 else 1.0
        for i in range(primer_len)
    )
    return min(1.0, total_weight / max_weight)


def evaluate_site(
    site: BindingSite,
    na_conc: float = 0.05,
    mg_conc: float = 0.0,
    dntp_conc: float = 0.0,
    primer_conc: float = 250e-9,
    anneal_temp: float = 55.0,
) -> EvaluatedSite:
    """
    Compute thermodynamic and quality properties for a single BindingSite.

    Parameters
    ----------
    site : BindingSite
    na_conc : float    Na⁺ in mol/L (default 50 mM)
    mg_conc : float    Mg²⁺ total in mol/L (default 0 → Na-only)
    dntp_conc : float  dNTP total in mol/L (chelates Mg²⁺ 1:1)
    primer_conc : float   strand concentration in mol/L (default 250 nM)
    anneal_temp : float   °C (default 55 — typical PCR annealing temperature)
    """
    primer = site.primer_seq
    primer_len = len(primer)

    # Use aligned_template for full Peyret mismatch NN lookup when available
    aligned_tmpl = site.aligned_template if site.aligned_template else None

    # Thermodynamics (full Peyret mismatch NN when template available)
    tm = calc_tm(
        primer,
        template            = aligned_tmpl,
        mismatch_positions  = site.mismatch_pos,
        three_prime_mismatch= site.three_prime_mm,
        na_conc             = na_conc,
        mg_conc             = mg_conc,
        dntp_conc           = dntp_conc,
        primer_conc         = primer_conc,
    )
    dg = calc_delta_g(
        primer,
        template            = aligned_tmpl,
        mismatch_positions  = site.mismatch_pos,
        three_prime_mismatch= site.three_prime_mm,
        temperature         = anneal_temp,
    )

    gc_p = gc_content(primer)
    gc_t = gc_content(site.aligned_template) if site.aligned_template else gc_p

    # Mismatch weight
    mm_w = _mismatch_weight(site.mismatch_pos, primer_len)

    # SW score normalised to perfect score
    perfect = float(primer_len * SW_MATCH)
    sw_norm = max(0.0, site.sw_score / perfect) if perfect > 0 else 0.0

    # 3'-end penalty
    prime3_pen = 0.5 if site.three_prime_mm else 0.0

    # GC factor: peaks at 0.525 (52.5% GC), declines outside 40–65%
    gc_center = 0.525
    gc_width  = 0.125    # half-width of the "good" GC range
    gc_dist   = abs(gc_p - gc_center)
    gc_factor = max(0.0, 1.0 - max(0.0, gc_dist - gc_width) / gc_width)

    binding_score = 100.0 * sw_norm * (1.0 - prime3_pen) * (1.0 - mm_w * 0.5) * (0.5 + 0.5 * gc_factor)

    return EvaluatedSite(
        site          = site,
        tm            = tm,
        delta_g       = dg,
        gc_primer     = round(gc_p, 3),
        gc_template   = round(gc_t, 3),
        binding_score = round(binding_score, 2),
        mm_weight     = round(mm_w, 4),
    )


def evaluate_sites(
    sites: List[BindingSite],
    na_conc: float = 0.05,
    mg_conc: float = 0.0,
    dntp_conc: float = 0.0,
    primer_conc: float = 250e-9,
    anneal_temp: float = 55.0,
) -> List[EvaluatedSite]:
    """Evaluate a list of BindingSite objects."""
    return [evaluate_site(s, na_conc, mg_conc, dntp_conc, primer_conc, anneal_temp)
            for s in sites]

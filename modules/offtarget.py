"""
Off-Target Analysis (Layer 10)
================================
Identifies and characterises all off-target amplicons produced by a primer pair.

Definitions:
  - Primary / intended amplicon: the highest-scored amplicon, or the one
    the user designates as the target.
  - Off-target amplicons: ALL other amplicons from the same primer pair.

Off-target report includes:
  - Genomic location (seq_id, start, end, size)
  - Off-target score (0–100; higher = more likely to amplify)
  - Reason: too_many_mm | non_target_locus | unexpected_size
  - Mismatch summary for each primer in the off-target binding

Off-target likelihood model:
  Likelihood is computed from the ScoredAmplicon.final_score of each off-target.
  A high final_score off-target is a serious concern; low-score ones may amplify
  inefficiently but are still reported (the user decides the threshold).

NOTE: This module does NOT attempt to model competitive PCR inhibition or
      amplification efficiency across multiple templates. That would require
      a kinetic model that is beyond scope.
"""

import logging
from dataclasses import dataclass, field
from typing import List, Optional

from .scoring import ScoredAmplicon

log = logging.getLogger(__name__)


@dataclass
class OffTargetHit:
    """Summary record for a single off-target amplicon."""
    pair_name:    str
    seq_id:       str
    start:        int
    end:          int
    size:         int
    sequence:     str
    gc_fraction:  float
    fwd_mm:       int
    rev_mm:       int
    fwd_tm:       float
    rev_tm:       float
    offtarget_score: float  # 0–100; higher = more likely to amplify
    reasons:      List[str] = field(default_factory=list)


def _classify_offtarget(sa: ScoredAmplicon, intended: Optional[ScoredAmplicon]) -> List[str]:
    """Return a list of reason strings explaining why this is an off-target."""
    reasons = []
    amp = sa.amplicon

    if intended is not None:
        int_amp = intended.amplicon
        if amp.seq_id != int_amp.seq_id:
            reasons.append("non_target_locus")
        elif amp.start != int_amp.start or amp.end != int_amp.end:
            reasons.append("unexpected_genomic_position")

    if amp.fwd_mm + amp.rev_mm > 0:
        reasons.append(f"mismatch_binding(fwd={amp.fwd_mm},rev={amp.rev_mm})")

    if amp.length < 50 or amp.length > 2000:
        reasons.append(f"unexpected_size({amp.length}bp)")

    if amp.hit.fwd_site.site.three_prime_mm or amp.hit.rev_site.site.three_prime_mm:
        reasons.append("3prime_mismatch")

    if not reasons:
        reasons.append("alternate_locus")

    return reasons


def analyse_offtargets(
    scored_amplicons: List[ScoredAmplicon],
    max_offtargets: int = 50,
) -> List[OffTargetHit]:
    """
    Identify off-target amplicons from a ranked amplicon list.

    The first (highest-scored) amplicon marked as intended is the primary target.
    All other amplicons are off-targets.

    Parameters
    ----------
    scored_amplicons : list of ScoredAmplicon (sorted by final_score desc)
    max_offtargets : int
        Cap on how many off-target records to return.

    Returns
    -------
    List of OffTargetHit sorted by offtarget_score descending.
    """
    if not scored_amplicons:
        return []

    # Identify the intended (primary) amplicon
    intended = next((sa for sa in scored_amplicons if sa.is_intended), scored_amplicons[0])

    offtargets: List[OffTargetHit] = []
    for sa in scored_amplicons:
        if sa.is_intended:
            continue

        amp = sa.amplicon
        reasons = _classify_offtarget(sa, intended)

        ot = OffTargetHit(
            pair_name       = amp.pair_name,
            seq_id          = amp.seq_id,
            start           = amp.start,
            end             = amp.end,
            size            = amp.length,
            sequence        = amp.sequence[:80] + ("..." if len(amp.sequence) > 80 else ""),
            gc_fraction     = amp.gc_fraction,
            fwd_mm          = amp.fwd_mm,
            rev_mm          = amp.rev_mm,
            fwd_tm          = amp.fwd_tm,
            rev_tm          = amp.rev_tm,
            offtarget_score = sa.final_score,
            reasons         = reasons,
        )
        offtargets.append(ot)

        if len(offtargets) >= max_offtargets:
            log.warning("Off-target cap (%d) reached for pair '%s'.", max_offtargets, amp.pair_name)
            break

    offtargets.sort(key=lambda o: o.offtarget_score, reverse=True)
    return offtargets


@dataclass
class OffTargetSummary:
    """Aggregated off-target statistics for a primer pair."""
    pair_name:          str
    total_offtargets:   int
    high_risk:          int    # off-target score >= 70
    medium_risk:        int    # off-target score 40–69
    low_risk:           int    # off-target score < 40
    specificity_index:  float  # 0–100; 100 = perfectly specific
    hits:               List[OffTargetHit] = field(default_factory=list)


def summarise_offtargets(
    pair_name: str,
    offtargets: List[OffTargetHit],
) -> OffTargetSummary:
    """
    Aggregate off-target hits into a summary with a specificity index.

    Specificity index:
      SI = 100 * exp(-0.05 * total_offtargets - 0.1 * high_risk)
      Ranges 0–100; decreases rapidly with many high-risk off-targets.
    """
    import math
    total  = len(offtargets)
    high   = sum(1 for o in offtargets if o.offtarget_score >= 70)
    medium = sum(1 for o in offtargets if 40 <= o.offtarget_score < 70)
    low    = sum(1 for o in offtargets if o.offtarget_score < 40)

    si = 100.0 * math.exp(-0.05 * total - 0.1 * high)

    return OffTargetSummary(
        pair_name         = pair_name,
        total_offtargets  = total,
        high_risk         = high,
        medium_risk       = medium,
        low_risk          = low,
        specificity_index = round(si, 2),
        hits              = offtargets,
    )

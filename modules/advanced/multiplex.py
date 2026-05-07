"""
Multiplex PCR Simulation
=========================
Simulates multiplexed PCR where multiple primer pairs are combined in a
single reaction. The key concern is cross-pair interactions:
  - Inter-pair dimer formation
  - Competing amplicons from different pairs that overlap in size
  - Tm incompatibility across pairs

This module:
  1. Checks all cross-pair primer dimer combinations
  2. Flags pairs with overlapping amplicon sizes (difficult to resolve on gel)
  3. Reports a multiplex compatibility score

Limitation:
  No competitive inhibition modelling. Relative amplification efficiency
  in multiplex PCR depends on template abundance and many kinetic factors
  that are outside the scope of this pipeline.
"""

from dataclasses import dataclass, field
from typing import List, Tuple

from .primer_dimer import check_all_dimers, DimerResult
from ..scoring import ScoredAmplicon


@dataclass
class MultiplexConflict:
    pair_a:   str
    pair_b:   str
    conflict: str    # description of the conflict
    severity: str    # "high", "medium", "low"


@dataclass
class MultiplexReport:
    pair_names:      List[str]
    conflicts:       List[MultiplexConflict] = field(default_factory=list)
    dimer_warnings:  List[str]               = field(default_factory=list)
    compat_score:    float = 100.0           # 0–100; 100 = fully compatible
    recommendation:  str  = ""


def check_multiplex_compatibility(
    pairs: List[dict],   # list of {"name": str, "fwd": str, "rev": str, "amplicons": List[ScoredAmplicon]}
    size_overlap_threshold: int = 20,  # bp; amplicons within this range flagged
    tm_diff_threshold: float = 5.0,   # °C; Tm difference between pairs flagged
) -> MultiplexReport:
    """
    Check compatibility of multiple primer pairs for multiplex PCR.

    Parameters
    ----------
    pairs : list of dicts, each with keys:
        name     : str
        fwd      : str  forward primer sequence
        rev      : str  reverse primer sequence
        amplicons: list of ScoredAmplicon objects (can be empty)
    size_overlap_threshold : int
    tm_diff_threshold : float

    Returns
    -------
    MultiplexReport
    """
    pair_names = [p["name"] for p in pairs]
    conflicts: List[MultiplexConflict] = []
    dimer_warnings: List[str] = []

    # 1. Cross-pair dimer checks
    for i in range(len(pairs)):
        for j in range(i, len(pairs)):
            pa, pb = pairs[i], pairs[j]
            if i == j:
                # Same pair: fwd vs rev
                dimers = check_all_dimers(pa["fwd"], pa["rev"])
            else:
                # Cross pair: all 4 combinations
                dimers = [
                    check_all_dimers(pa["fwd"], pb["fwd"])[2],  # hetero
                    check_all_dimers(pa["fwd"], pb["rev"])[2],
                    check_all_dimers(pa["rev"], pb["fwd"])[2],
                    check_all_dimers(pa["rev"], pb["rev"])[2],
                ]

            for dr in (dimers if isinstance(dimers, list) else [dimers]):
                if dr.warning:
                    label = f"Pairs {pa['name']}×{pb['name']}" if i != j else f"Pair {pa['name']}"
                    dimer_warnings.append(f"{label}: {dr.warning}")
                    severity = "high" if (dr.three_prime_a or dr.three_prime_b) else "medium"
                    conflicts.append(MultiplexConflict(
                        pair_a   = pa["name"],
                        pair_b   = pb["name"],
                        conflict = f"Primer dimer (score={dr.dimer_score:.1f}, type={dr.dimer_type})",
                        severity = severity,
                    ))

    # 2. Amplicon size overlap
    all_amplicons: List[Tuple[str, int]] = []   # (pair_name, size)
    for p in pairs:
        for sa in p.get("amplicons", []):
            if sa.is_intended:
                all_amplicons.append((p["name"], sa.amplicon.length))

    for i in range(len(all_amplicons)):
        for j in range(i + 1, len(all_amplicons)):
            name_i, size_i = all_amplicons[i]
            name_j, size_j = all_amplicons[j]
            if abs(size_i - size_j) <= size_overlap_threshold:
                conflicts.append(MultiplexConflict(
                    pair_a   = name_i,
                    pair_b   = name_j,
                    conflict = (
                        f"Amplicon size overlap: {size_i} bp vs {size_j} bp "
                        f"(diff={abs(size_i-size_j)} bp ≤ threshold {size_overlap_threshold} bp). "
                        "Difficult to resolve on gel."
                    ),
                    severity = "medium",
                ))

    # 3. Tm compatibility across pairs
    tm_values: List[Tuple[str, float]] = []
    for p in pairs:
        for sa in p.get("amplicons", []):
            if sa.is_intended:
                avg_tm = (sa.amplicon.fwd_tm + sa.amplicon.rev_tm) / 2
                tm_values.append((p["name"], avg_tm))

    for i in range(len(tm_values)):
        for j in range(i + 1, len(tm_values)):
            name_i, tm_i = tm_values[i]
            name_j, tm_j = tm_values[j]
            if abs(tm_i - tm_j) > tm_diff_threshold:
                conflicts.append(MultiplexConflict(
                    pair_a   = name_i,
                    pair_b   = name_j,
                    conflict = (
                        f"Tm incompatibility: pair {name_i} avg Tm={tm_i:.1f}°C vs "
                        f"pair {name_j} avg Tm={tm_j:.1f}°C "
                        f"(diff={abs(tm_i-tm_j):.1f}°C > {tm_diff_threshold}°C)."
                    ),
                    severity = "low",
                ))

    # Compatibility score: starts at 100, penalised per conflict
    penalty_map = {"high": 20, "medium": 10, "low": 5}
    score = 100.0
    for c in conflicts:
        score -= penalty_map.get(c.severity, 5)
    score = max(0.0, score)

    high_count = sum(1 for c in conflicts if c.severity == "high")
    if high_count > 0:
        rec = f"⚠ {high_count} HIGH-severity conflict(s). Redesign primers before multiplexing."
    elif conflicts:
        rec = "Moderate conflicts detected. Test empirically with gradient annealing."
    else:
        rec = "Primer pairs appear compatible for multiplex PCR."

    return MultiplexReport(
        pair_names     = pair_names,
        conflicts      = conflicts,
        dimer_warnings = dimer_warnings,
        compat_score   = round(score, 1),
        recommendation = rec,
    )

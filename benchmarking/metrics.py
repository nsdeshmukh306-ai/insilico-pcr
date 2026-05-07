"""
Benchmarking Metrics
====================
Computes quantitative comparison metrics between this pipeline's output
and reference results from external tools (UCSC In-Silico PCR, Primer-BLAST)
or ground-truth coordinates.

Metrics
-------
amplicon_coord_deviation : int
    |predicted_start - reference_start| + |predicted_end - reference_end|
    In bp. Zero means exact coordinate match.

tm_deviation : float
    |predicted_Tm - reference_Tm| in °C.
    Reference Tm may be from the external tool or from a wet-lab measurement.

off_target_count_delta : int
    predicted_off_target_count - reference_off_target_count
    Positive = this tool predicts more off-targets than reference (conservative).
    Negative = this tool misses off-targets present in reference (liberal).

runtime_seconds : float
    Wall-clock seconds for the pipeline run.

sensitivity : float
    Fraction of reference amplicons recovered by this pipeline (0–1).
    sensitivity = TP / (TP + FN)

precision : float
    Fraction of predicted amplicons that match a reference (0–1).
    precision = TP / (TP + FP)
    A predicted amplicon is a true positive if its coordinates overlap a
    reference amplicon by ≥ MIN_OVERLAP_FRACTION of the reference length.

f1 : float
    Harmonic mean of sensitivity and precision.

Usage
-----
>>> from benchmarking.metrics import compute_metrics, BenchmarkResult
>>> result = compute_metrics(predicted_amplicons, reference_amplicons)
>>> print(result.summary_text())
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

# Minimum fraction of reference amplicon length that must be overlapped by a
# predicted amplicon for it to count as a true positive.
MIN_OVERLAP_FRACTION: float = 0.80


@dataclass
class AmpliconCoord:
    """Minimal amplicon representation for benchmarking (no pipeline deps)."""
    seq_id: str
    start:  int       # 0-based
    end:    int       # 0-based exclusive
    tm_fwd: float = 0.0
    tm_rev: float = 0.0
    name:   str   = ""

    @property
    def length(self) -> int:
        return max(0, self.end - self.start)


@dataclass
class PairBenchmarkResult:
    """Benchmark result for a single primer pair."""
    pair_name:               str
    predicted:               List[AmpliconCoord]
    reference:               List[AmpliconCoord]
    true_positives:          int   = 0
    false_positives:         int   = 0
    false_negatives:         int   = 0
    sensitivity:             float = 0.0
    precision:               float = 0.0
    f1:                      float = 0.0
    mean_coord_deviation_bp: float = 0.0   # avg over matched TPs
    mean_tm_deviation_c:     float = 0.0   # avg |ΔTm| over matched TPs
    off_target_delta:        int   = 0     # predicted OT - reference OT
    runtime_seconds:         float = 0.0
    notes:                   str   = ""


@dataclass
class BenchmarkReport:
    """Aggregated benchmark across all primer pairs in a run."""
    tool_version:      str = "insilico_pcr-2.0"
    reference_tool:    str = "unknown"
    genome:            str = ""
    pair_results:      List[PairBenchmarkResult] = field(default_factory=list)
    total_runtime_s:   float = 0.0

    # Aggregate metrics (computed by finalize())
    macro_sensitivity: float = 0.0
    macro_precision:   float = 0.0
    macro_f1:          float = 0.0
    mean_coord_dev:    float = 0.0
    mean_tm_dev:       float = 0.0

    def finalize(self) -> None:
        """Compute macro-averaged metrics from pair_results."""
        if not self.pair_results:
            return
        self.macro_sensitivity = sum(r.sensitivity for r in self.pair_results) / len(self.pair_results)
        self.macro_precision   = sum(r.precision   for r in self.pair_results) / len(self.pair_results)
        self.macro_f1          = sum(r.f1           for r in self.pair_results) / len(self.pair_results)
        coord_vals = [r.mean_coord_deviation_bp for r in self.pair_results if r.true_positives > 0]
        tm_vals    = [r.mean_tm_deviation_c      for r in self.pair_results if r.true_positives > 0]
        self.mean_coord_dev = sum(coord_vals) / len(coord_vals) if coord_vals else 0.0
        self.mean_tm_dev    = sum(tm_vals)    / len(tm_vals)    if tm_vals    else 0.0

    def to_dict(self) -> dict:
        return {
            "tool_version":      self.tool_version,
            "reference_tool":    self.reference_tool,
            "genome":            self.genome,
            "total_runtime_s":   round(self.total_runtime_s, 3),
            "macro_sensitivity": round(self.macro_sensitivity, 4),
            "macro_precision":   round(self.macro_precision,   4),
            "macro_f1":          round(self.macro_f1,           4),
            "mean_coord_dev_bp": round(self.mean_coord_dev,    2),
            "mean_tm_dev_c":     round(self.mean_tm_dev,       2),
            "pair_results": [
                {
                    "pair_name":     r.pair_name,
                    "true_positives":  r.true_positives,
                    "false_positives": r.false_positives,
                    "false_negatives": r.false_negatives,
                    "sensitivity":     round(r.sensitivity, 4),
                    "precision":       round(r.precision,   4),
                    "f1":              round(r.f1,           4),
                    "mean_coord_dev_bp": round(r.mean_coord_deviation_bp, 2),
                    "mean_tm_dev_c":     round(r.mean_tm_deviation_c,     2),
                    "off_target_delta":  r.off_target_delta,
                    "runtime_s":         round(r.runtime_seconds, 3),
                    "notes":             r.notes,
                }
                for r in self.pair_results
            ],
        }

    def to_text(self) -> str:
        lines = [
            "=" * 72,
            "IN-SILICO PCR BENCHMARK REPORT",
            "=" * 72,
            f"Tool version   : {self.tool_version}",
            f"Reference tool : {self.reference_tool}",
            f"Genome         : {self.genome}",
            f"Total runtime  : {self.total_runtime_s:.2f} s",
            "",
            "AGGREGATE METRICS",
            "-" * 40,
            f"  Macro sensitivity : {self.macro_sensitivity:.4f}",
            f"  Macro precision   : {self.macro_precision:.4f}",
            f"  Macro F1          : {self.macro_f1:.4f}",
            f"  Mean coord dev    : {self.mean_coord_dev:.1f} bp",
            f"  Mean Tm deviation : {self.mean_tm_dev:.2f} °C",
            "",
            "PER-PAIR RESULTS",
            "-" * 40,
        ]
        for r in self.pair_results:
            lines += [
                f"  Pair: {r.pair_name}",
                f"    TP={r.true_positives}  FP={r.false_positives}  FN={r.false_negatives}",
                f"    Sensitivity={r.sensitivity:.3f}  Precision={r.precision:.3f}  F1={r.f1:.3f}",
                f"    Coord dev={r.mean_coord_deviation_bp:.1f} bp  Tm dev={r.mean_tm_deviation_c:.2f} °C",
                f"    Off-target delta={r.off_target_delta:+d}  Runtime={r.runtime_seconds:.3f} s",
            ]
            if r.notes:
                lines.append(f"    Notes: {r.notes}")
            lines.append("")
        lines.append("=" * 72)
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Core metric computation
# ---------------------------------------------------------------------------
def _overlaps(pred: AmpliconCoord, ref: AmpliconCoord,
              min_frac: float = MIN_OVERLAP_FRACTION) -> bool:
    """True if pred overlaps ref by at least min_frac of ref.length."""
    if pred.seq_id != ref.seq_id:
        return False
    overlap = max(0, min(pred.end, ref.end) - max(pred.start, ref.start))
    if ref.length == 0:
        return False
    return (overlap / ref.length) >= min_frac


def compute_pair_metrics(
    pair_name:   str,
    predicted:   List[AmpliconCoord],
    reference:   List[AmpliconCoord],
    runtime_s:   float = 0.0,
    notes:       str   = "",
) -> PairBenchmarkResult:
    """
    Compare predicted amplicons against reference amplicons for one primer pair.

    A predicted amplicon is a true positive if it overlaps a reference amplicon
    by ≥ MIN_OVERLAP_FRACTION of the reference length (greedy 1-to-1 matching).

    Parameters
    ----------
    pair_name  : str
    predicted  : list of AmpliconCoord  (from this pipeline)
    reference  : list of AmpliconCoord  (from UCSC, Primer-BLAST, or wet-lab)
    runtime_s  : float  wall-clock seconds for this pair
    notes      : str    free-text annotation

    Returns
    -------
    PairBenchmarkResult
    """
    matched_ref   = set()
    true_pos_pairs: List[Tuple[AmpliconCoord, AmpliconCoord]] = []

    for pred in predicted:
        for i, ref in enumerate(reference):
            if i in matched_ref:
                continue
            if _overlaps(pred, ref):
                matched_ref.add(i)
                true_pos_pairs.append((pred, ref))
                break

    tp = len(true_pos_pairs)
    fp = len(predicted) - tp
    fn = len(reference)  - tp

    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    precision   = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    f1 = (2 * sensitivity * precision / (sensitivity + precision)
          if (sensitivity + precision) > 0 else 0.0)

    # Coordinate and Tm deviations (only over matched TPs)
    coord_devs, tm_devs = [], []
    for pred, ref in true_pos_pairs:
        coord_devs.append(abs(pred.start - ref.start) + abs(pred.end - ref.end))
        # Tm deviation: use average of fwd+rev Tm if available
        pred_tm = (pred.tm_fwd + pred.tm_rev) / 2.0
        ref_tm  = (ref.tm_fwd  + ref.tm_rev)  / 2.0
        if ref_tm > 0:
            tm_devs.append(abs(pred_tm - ref_tm))

    mean_coord = sum(coord_devs) / len(coord_devs) if coord_devs else 0.0
    mean_tm    = sum(tm_devs)    / len(tm_devs)    if tm_devs    else 0.0

    # Off-target: non-primary amplicons
    ref_ot   = max(0, len(reference) - 1)
    pred_ot  = max(0, len(predicted) - 1)

    return PairBenchmarkResult(
        pair_name               = pair_name,
        predicted               = predicted,
        reference               = reference,
        true_positives          = tp,
        false_positives         = fp,
        false_negatives         = fn,
        sensitivity             = round(sensitivity, 4),
        precision               = round(precision,   4),
        f1                      = round(f1,           4),
        mean_coord_deviation_bp = round(mean_coord,   2),
        mean_tm_deviation_c     = round(mean_tm,      2),
        off_target_delta        = pred_ot - ref_ot,
        runtime_seconds         = round(runtime_s,    3),
        notes                   = notes,
    )

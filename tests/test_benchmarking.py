"""
Tests for the benchmarking module (metrics computation).
Does NOT require BioPython or the full pipeline — tests metrics logic only.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest
from benchmarking.metrics import (
    AmpliconCoord,
    compute_pair_metrics,
    BenchmarkReport,
)


def _coord(seq_id, start, end, name=""):
    return AmpliconCoord(seq_id=seq_id, start=start, end=end, name=name)


class TestAmpliconCoord:
    def test_length(self):
        c = _coord("chr1", 100, 520)
        assert c.length == 420

    def test_zero_length(self):
        c = _coord("chr1", 100, 100)
        assert c.length == 0


class TestComputePairMetrics:
    def test_perfect_match(self):
        """Exactly matching prediction → TP=1, FP=0, FN=0, sensitivity=1, precision=1."""
        pred = [_coord("chr1", 500, 920)]
        ref  = [_coord("chr1", 500, 920)]
        r    = compute_pair_metrics("pair1", pred, ref)
        assert r.true_positives  == 1
        assert r.false_positives == 0
        assert r.false_negatives == 0
        assert r.sensitivity     == 1.0
        assert r.precision       == 1.0
        assert r.f1              == 1.0

    def test_no_predictions(self):
        """No predictions → TP=0, FP=0, FN=1, sensitivity=0, precision=0."""
        pred = []
        ref  = [_coord("chr1", 500, 920)]
        r    = compute_pair_metrics("pair1", pred, ref)
        assert r.true_positives  == 0
        assert r.false_positives == 0
        assert r.false_negatives == 1
        assert r.sensitivity     == 0.0
        assert r.precision       == 0.0
        assert r.f1              == 0.0

    def test_no_reference(self):
        """No reference → TP=0, FP=1, FN=0, sensitivity=0, precision=0."""
        pred = [_coord("chr1", 500, 920)]
        ref  = []
        r    = compute_pair_metrics("pair1", pred, ref)
        assert r.true_positives  == 0
        assert r.false_positives == 1
        assert r.false_negatives == 0
        assert r.sensitivity     == 0.0
        assert r.precision       == 0.0

    def test_near_match_counts_as_tp(self):
        """Prediction overlapping ≥80% of reference counts as TP."""
        ref  = [_coord("chr1", 500, 920)]     # length=420
        # Pred shifted by 5 bp: overlap=415/420 ≈ 98.8% → TP
        pred = [_coord("chr1", 505, 925)]
        r    = compute_pair_metrics("pair1", pred, ref)
        assert r.true_positives == 1
        assert r.false_positives == 0

    def test_small_overlap_counts_as_fp(self):
        """Prediction overlapping <80% of reference is FP."""
        ref  = [_coord("chr1", 500, 920)]     # length=420
        # Pred covers only 100 bp of the reference out of 420 → 24% overlap → FP
        pred = [_coord("chr1", 820, 950)]     # overlap with ref: 820-920 = 100 bp
        r    = compute_pair_metrics("pair1", pred, ref)
        assert r.true_positives  == 0
        assert r.false_positives == 1
        assert r.false_negatives == 1

    def test_wrong_chromosome_is_fp(self):
        """Prediction on wrong chromosome cannot match reference → FP."""
        ref  = [_coord("chr1", 500, 920)]
        pred = [_coord("chr2", 500, 920)]
        r    = compute_pair_metrics("pair1", pred, ref)
        assert r.true_positives == 0

    def test_coordinate_deviation_computed(self):
        """Coord deviation = |pred_start - ref_start| + |pred_end - ref_end|."""
        ref  = [_coord("chr1", 500, 920)]
        pred = [_coord("chr1", 503, 917)]   # start +3, end -3 → dev=6
        r    = compute_pair_metrics("pair1", pred, ref)
        assert r.true_positives == 1
        assert abs(r.mean_coord_deviation_bp - 6.0) < 0.1

    def test_f1_score_formula(self):
        """F1 = 2*P*S / (P+S)."""
        ref  = [_coord("chr1", 100, 500), _coord("chr2", 200, 600)]
        pred = [_coord("chr1", 100, 500)]   # misses chr2 → FN=1
        r    = compute_pair_metrics("p", pred, ref)
        expected_f1 = 2 * 1.0 * 0.5 / (1.0 + 0.5)
        assert abs(r.f1 - expected_f1) < 1e-4

    def test_multiple_tp(self):
        """Multiple matching amplicons all counted as TP."""
        ref  = [_coord("chr1", 100, 400), _coord("chr2", 200, 600)]
        pred = [_coord("chr1", 100, 400), _coord("chr2", 200, 600)]
        r    = compute_pair_metrics("p", pred, ref)
        assert r.true_positives  == 2
        assert r.false_positives == 0
        assert r.false_negatives == 0
        assert r.sensitivity     == 1.0
        assert r.precision       == 1.0

    def test_runtime_stored(self):
        r = compute_pair_metrics("p", [], [], runtime_s=1.234)
        assert abs(r.runtime_seconds - 1.234) < 0.001


class TestBenchmarkReport:
    def _make_report(self):
        r1 = compute_pair_metrics("p1",
            [_coord("chr1", 100, 500)],
            [_coord("chr1", 100, 500)],
        )
        r2 = compute_pair_metrics("p2",
            [_coord("chr2", 200, 600)],
            [_coord("chr2", 250, 650)],   # slightly off
        )
        report = BenchmarkReport(reference_tool="test", genome="test_genome")
        report.pair_results = [r1, r2]
        report.finalize()
        return report

    def test_finalize_computes_macro_metrics(self):
        report = self._make_report()
        assert 0 <= report.macro_sensitivity <= 1
        assert 0 <= report.macro_precision   <= 1
        assert 0 <= report.macro_f1          <= 1

    def test_to_dict_has_required_keys(self):
        report = self._make_report()
        d = report.to_dict()
        for key in ("tool_version", "reference_tool", "genome",
                    "macro_sensitivity", "macro_precision", "macro_f1",
                    "mean_coord_dev_bp", "mean_tm_dev_c", "pair_results"):
            assert key in d, f"Missing key: {key}"

    def test_to_text_contains_header(self):
        report = self._make_report()
        text = report.to_text()
        assert "BENCHMARK REPORT" in text
        assert "macro_sensitivity" in text.lower() or "Macro sensitivity" in text

    def test_pair_results_in_dict(self):
        report = self._make_report()
        d = report.to_dict()
        assert len(d["pair_results"]) == 2
        names = [p["pair_name"] for p in d["pair_results"]]
        assert "p1" in names and "p2" in names

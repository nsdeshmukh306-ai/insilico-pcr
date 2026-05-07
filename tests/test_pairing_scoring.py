"""Unit tests for pairing engine, amplicon extraction, and scoring."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest
from modules.binding_search import BindingSite
from modules.binding_eval import EvaluatedSite, evaluate_site
from modules.pairing_engine import pair_binding_sites
from modules.amplicon import extract_amplicon, amplicon_length_score
from modules.scoring import score_amplicon, DEFAULT_WEIGHTS


def make_evaluated_site(
    seq_id="chr1",
    start=0,
    end=20,
    strand="+",
    primer="GCTAGCTAGCTAGCTAGCTA",
    mm=0,
    mm_pos=None,
    sw_score=40.0,
):
    """Helper: build an EvaluatedSite for testing."""
    raw_site = BindingSite(
        seq_id           = seq_id,
        start            = start,
        end              = end,
        strand           = strand,
        primer_seq       = primer,
        aligned_template = primer,
        sw_score         = sw_score,
        mismatch_count   = mm,
        mismatch_pos     = mm_pos or [],
        gap_count        = 0,
        three_prime_mm   = False,
    )
    return evaluate_site(raw_site)


class TestPairingEngine:
    def test_valid_pair(self):
        fwd = make_evaluated_site(start=100, end=120, strand="+")
        rev = make_evaluated_site(start=500, end=520, strand="-")
        hits = pair_binding_sites("p1", [fwd], [rev], min_amplicon_size=50, max_amplicon_size=1000)
        assert len(hits) == 1
        assert hits[0].amplicon_size == 520 - 100

    def test_wrong_strand_rejected(self):
        # Both fwd — should not pair
        fwd1 = make_evaluated_site(start=100, end=120, strand="+")
        fwd2 = make_evaluated_site(start=500, end=520, strand="+")  # wrong strand
        hits = pair_binding_sites("p1", [fwd1], [fwd2], min_amplicon_size=50, max_amplicon_size=1000)
        assert len(hits) == 0

    def test_too_small_amplicon_rejected(self):
        fwd = make_evaluated_site(start=100, end=120, strand="+")
        rev = make_evaluated_site(start=130, end=150, strand="-")  # 50 bp total
        hits = pair_binding_sites("p1", [fwd], [rev], min_amplicon_size=200, max_amplicon_size=1000)
        assert len(hits) == 0

    def test_too_large_amplicon_rejected(self):
        fwd = make_evaluated_site(start=100, end=120, strand="+")
        rev = make_evaluated_site(start=5000, end=5020, strand="-")
        hits = pair_binding_sites("p1", [fwd], [rev], min_amplicon_size=50, max_amplicon_size=1000)
        assert len(hits) == 0

    def test_wrong_order_rejected(self):
        # Rev primer upstream of fwd — invalid
        fwd = make_evaluated_site(start=500, end=520, strand="+")
        rev = make_evaluated_site(start=100, end=120, strand="-")
        hits = pair_binding_sites("p1", [fwd], [rev], min_amplicon_size=50, max_amplicon_size=1000)
        assert len(hits) == 0

    def test_different_contigs_not_paired(self):
        fwd = make_evaluated_site(seq_id="chr1", start=100, end=120, strand="+")
        rev = make_evaluated_site(seq_id="chr2", start=500, end=520, strand="-")
        hits = pair_binding_sites("p1", [fwd], [rev])
        assert len(hits) == 0


class TestAmpliconLengthScore:
    def test_optimal_length(self):
        # ~200 bp should give a high score
        score = amplicon_length_score(200)
        assert score > 0.9

    def test_very_short(self):
        score = amplicon_length_score(30)
        assert score < 0.15

    def test_very_long(self):
        score = amplicon_length_score(5000)
        assert score < 0.2

    def test_zero(self):
        assert amplicon_length_score(0) == 0.0


class TestScoring:
    def _make_amplicon(self):
        from modules.pairing_engine import PrimerPairHit
        from modules.amplicon import Amplicon
        fwd = make_evaluated_site(start=100, end=120, strand="+", sw_score=40.0)
        rev = make_evaluated_site(start=500, end=520, strand="-", sw_score=40.0)
        hit = PrimerPairHit(
            pair_name      = "test",
            seq_id         = "chr1",
            fwd_site       = fwd,
            rev_site       = rev,
            amplicon_start = 100,
            amplicon_end   = 520,
            amplicon_size  = 420,
        )
        return Amplicon(
            pair_name         = "test",
            seq_id            = "chr1",
            start             = 100,
            end               = 520,
            length            = 420,
            sequence          = "A" * 420,
            gc_fraction       = 0.50,
            fwd_primer        = fwd.site.primer_seq,
            rev_primer        = rev.site.primer_seq,
            fwd_tm            = fwd.tm,
            rev_tm            = rev.tm,
            fwd_mm            = 0,
            rev_mm            = 0,
            fwd_binding_score = fwd.binding_score,
            rev_binding_score = rev.binding_score,
            hit               = hit,
        )

    def test_score_range(self):
        amp = self._make_amplicon()
        sa = score_amplicon(amp, off_target_count=0, is_intended=True, max_mismatches=3)
        assert 0 <= sa.final_score <= 100

    def test_more_mismatches_lower_score(self):
        from modules.amplicon import Amplicon
        from modules.pairing_engine import PrimerPairHit

        fwd0 = make_evaluated_site(start=100, end=120, strand="+", mm=0)
        fwd2 = make_evaluated_site(start=100, end=120, strand="+", mm=2,
                                    mm_pos=[3, 7], sw_score=36.0)
        rev  = make_evaluated_site(start=500, end=520, strand="-")

        def make_amp(fwd_es):
            hit = PrimerPairHit("t","chr1",fwd_es,rev,100,520,420)
            return Amplicon("t","chr1",100,520,420,"A"*420,0.5,
                            fwd_es.site.primer_seq, rev.site.primer_seq,
                            fwd_es.tm, rev.tm,
                            fwd_es.site.mismatch_count, rev.site.mismatch_count,
                            fwd_es.binding_score, rev.binding_score, hit)

        sa0 = score_amplicon(make_amp(fwd0), max_mismatches=3)
        sa2 = score_amplicon(make_amp(fwd2), max_mismatches=3)
        assert sa0.final_score > sa2.final_score

    def test_off_target_lowers_score(self):
        amp = self._make_amplicon()
        sa_0ot = score_amplicon(amp, off_target_count=0, is_intended=True)
        sa_5ot = score_amplicon(amp, off_target_count=5, is_intended=False)
        assert sa_0ot.final_score > sa_5ot.final_score

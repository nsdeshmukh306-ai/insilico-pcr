"""Unit tests for the thermodynamics module."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest
from modules.thermodynamics import (
    calc_tm,
    calc_delta_g,
    gc_content,
    complement,
    reverse_complement,
    calc_nn_thermodynamics,
    check_primer_quality,
    calc_tm_basic,
)


class TestComplement:
    def test_complement_simple(self):
        assert complement("ACGT") == "TGCA"

    def test_reverse_complement(self):
        assert reverse_complement("ACGT") == "ACGT"
        assert reverse_complement("AAAA") == "TTTT"
        assert reverse_complement("GCTA") == "TAGC"

    def test_reverse_complement_palindrome(self):
        # AATT is its own reverse complement
        assert reverse_complement("AATT") == "AATT"


class TestGCContent:
    def test_all_gc(self):
        assert gc_content("GCGCGC") == 1.0

    def test_all_at(self):
        assert gc_content("ATATAT") == 0.0

    def test_mixed(self):
        assert abs(gc_content("ACGT") - 0.5) < 1e-9

    def test_empty(self):
        assert gc_content("") == 0.0


class TestNNThermodynamics:
    def test_returns_tuple(self):
        dh, ds = calc_nn_thermodynamics("ATCGATCG")
        assert isinstance(dh, float)
        assert isinstance(ds, float)

    def test_dh_negative(self):
        # ΔH should be negative for stable duplexes
        dh, _ = calc_nn_thermodynamics("GCGCGCGCGC")
        assert dh < 0

    def test_mismatch_reduces_stability(self):
        dh_perfect, ds_perfect = calc_nn_thermodynamics("ATCGATCG")
        dh_mm, ds_mm = calc_nn_thermodynamics("ATCGATCG", mismatch_positions=[3])
        # Mismatch should reduce |ΔH| (less negative)
        assert dh_mm > dh_perfect


class TestTm:
    def test_tm_range(self):
        # A typical 20 bp primer should have Tm between 40–80 °C
        tm = calc_tm("GCTAGCTAGCTAGCTAGCTA")
        assert 40 < tm < 80, f"Tm={tm} out of expected range"

    def test_gc_rich_higher_tm(self):
        tm_gc = calc_tm("GCGCGCGCGCGCGCGCGCGC")  # all-GC
        tm_at = calc_tm("ATATATATATATATATATATAT")  # all-AT
        assert tm_gc > tm_at, "GC-rich primer should have higher Tm"

    def test_mismatch_lowers_tm(self):
        seq = "GCTAGCTAGCTAGCTAGCTA"
        tm_perfect = calc_tm(seq)
        tm_mm = calc_tm(seq, mismatch_positions=[5, 10])
        assert tm_mm < tm_perfect, "Mismatches should lower Tm"

    def test_three_prime_mm_extra_penalty(self):
        # The three_prime_mismatch flag applies an extra penalty multiplier
        # to the same terminal position — compare with/without the flag.
        seq = "GCTAGCTAGCTAGCTAGCTA"
        last = len(seq) - 1
        tm_no_flag   = calc_tm(seq, mismatch_positions=[last], three_prime_mismatch=False)
        tm_with_flag = calc_tm(seq, mismatch_positions=[last], three_prime_mismatch=True)
        assert tm_with_flag < tm_no_flag, "three_prime_mismatch=True should lower Tm further"

    def test_salt_correction_direction(self):
        """Lower [Na+] should give lower Tm (less stabilising ionic environment)."""
        tm_50mm = calc_tm("GCTAGCTAGCTAGCTAGCTA", na_conc=0.05)
        tm_200mm = calc_tm("GCTAGCTAGCTAGCTAGCTA", na_conc=0.20)
        assert tm_200mm > tm_50mm, "Higher [Na+] should raise Tm"

    def test_wallace_rule_sanity(self):
        # Wallace rule is valid only for short oligos; we just check direction
        seq = "ATCGATCG"
        tm_basic = calc_tm_basic(seq)
        assert tm_basic > 0


class TestDeltaG:
    def test_dg_negative(self):
        """ΔG at 37 °C should be negative for stable duplexes."""
        dg = calc_delta_g("GCTAGCTAGCTAGCTAGCTA")
        assert dg < 0, f"ΔG={dg} should be negative"

    def test_dg_less_negative_with_mismatches(self):
        seq = "GCTAGCTAGCTAGCTAGCTA"
        dg_perfect = calc_delta_g(seq)
        dg_mm = calc_delta_g(seq, mismatch_positions=[3, 7])
        assert dg_mm > dg_perfect, "Mismatches should make ΔG less negative"

    def test_higher_temp_less_stable(self):
        seq = "GCTAGCTAGCTAGCTAGCTA"
        dg_37 = calc_delta_g(seq, temperature=37.0)
        dg_72 = calc_delta_g(seq, temperature=72.0)
        assert dg_72 > dg_37, "At higher temperature ΔG should be less negative"


class TestPrimerQuality:
    def test_good_primer(self):
        q = check_primer_quality("GCTAGCTAGCTAGCTAGCTA")
        assert q["length_ok"]
        assert not q["low_complexity"]

    def test_short_primer_fails(self):
        q = check_primer_quality("ATCGATCG")
        assert not q["length_ok"]

    def test_homopolymer_flagged(self):
        q = check_primer_quality("AAAAAAAAAAAAAAAAAAAAATCG")
        assert q["low_complexity"]

    def test_gc_clamp(self):
        # Primer ending in GC should have gc_clamp_ok
        q = check_primer_quality("ATCGATCGATCGATCGATGC")
        assert q["gc_clamp_ok"]

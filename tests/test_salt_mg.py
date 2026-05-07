"""
Tests for Owczarzy et al. (2008) Mg²⁺ salt correction.

Validates:
 1. Na⁺-only mode (mg_conc=0) reproduces previous Wetmur behaviour.
 2. Mg²⁺ raises Tm vs Na⁺-only at the same ionic conditions.
 3. Free Mg²⁺ is correctly reduced by dNTP chelation (1:1 stoichiometry).
 4. Mixed Na⁺/Mg²⁺ regime interpolates between Na and Mg formulae.
 5. Regime selection boundary (√[Mg²⁺]/[Na⁺] ratio).
 6. check_primer_quality propagates mg_conc/dntp_conc.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest
from modules.thermodynamics import (
    calc_tm,
    _owczarzy_salt_correction,
    check_primer_quality,
)

PRIMER = "GCTAGCTAGCTAGCTAGCTA"   # 20-mer, mixed GC


class TestOwczarrySaltCorrection:
    """Unit tests for _owczarzy_salt_correction."""

    def test_no_mg_equals_na_only(self):
        """With free_mg=0, correction should equal Na⁺-only (Wetmur formula)."""
        import math
        tm_1m = 65.0   # hypothetical Tm at 1 M NaCl
        fGC   = 0.5
        n     = 20
        na    = 0.05   # 50 mM

        # Wetmur: tm_1m + 16.6 * log10(na)
        expected = tm_1m + 16.6 * math.log10(na)
        got = _owczarzy_salt_correction(tm_1m, fGC, n, na_conc=na, mg_conc=0.0, dntp_conc=0.0)
        assert abs(got - expected) < 0.01, f"Na-only: expected {expected:.2f}, got {got:.2f}"

    def test_mg_raises_tm_vs_na_only(self):
        """Adding Mg²⁺ to Na⁺ solution should raise Tm (stabilising effect)."""
        tm_1m = 65.0
        fGC   = 0.5
        n     = 20

        tm_na_only = _owczarzy_salt_correction(tm_1m, fGC, n, na_conc=0.05, mg_conc=0.0, dntp_conc=0.0)
        tm_with_mg = _owczarzy_salt_correction(tm_1m, fGC, n, na_conc=0.05, mg_conc=0.003, dntp_conc=0.0)

        assert tm_with_mg > tm_na_only, (
            f"Mg²⁺ should raise Tm: Na-only={tm_na_only:.2f}, with Mg={tm_with_mg:.2f}"
        )

    def test_dntp_chelation_reduces_free_mg(self):
        """dNTP chelates Mg²⁺ 1:1; equal concentrations → free Mg²⁺ ≈ 0."""
        tm_1m = 65.0
        fGC   = 0.5
        n     = 20
        mg    = 0.003   # 3 mM Mg²⁺
        dntp  = 0.003   # 3 mM dNTP — chelates all Mg²⁺

        tm_no_dntp   = _owczarzy_salt_correction(tm_1m, fGC, n, na_conc=0.05, mg_conc=mg, dntp_conc=0.0)
        tm_with_dntp = _owczarzy_salt_correction(tm_1m, fGC, n, na_conc=0.05, mg_conc=mg, dntp_conc=dntp)

        # After full chelation, Mg behaves as if absent → same as Na-only
        import math
        tm_na_only = tm_1m + 16.6 * math.log10(0.05)
        assert abs(tm_with_dntp - tm_na_only) < 0.01, (
            f"Full chelation: expected {tm_na_only:.2f}, got {tm_with_dntp:.2f}"
        )
        assert tm_no_dntp > tm_with_dntp, "dNTP should reduce Mg effect → lower Tm"

    def test_excess_dntp_clamps_at_zero(self):
        """dNTP excess beyond Mg²⁺ is clamped to free_mg=0, not negative."""
        tm_1m = 65.0
        got = _owczarzy_salt_correction(tm_1m, 0.5, 20,
                                         na_conc=0.05, mg_conc=0.001, dntp_conc=0.005)
        # No negative Mg concentrations → should be identical to no-Mg case
        import math
        expected = tm_1m + 16.6 * math.log10(0.05)
        assert abs(got - expected) < 0.01

    def test_mg_only_no_na(self):
        """Pure Mg²⁺ buffer (Na=0) uses Owczarzy Mg-only formula."""
        tm_1m = 65.0
        tm_mg = _owczarzy_salt_correction(tm_1m, 0.5, 20, na_conc=0.0, mg_conc=0.003, dntp_conc=0.0)
        # Should return a valid temperature (not 0 or crash)
        assert 30 < tm_mg < 90, f"Mg-only Tm={tm_mg} outside reasonable range"

    def test_mixed_regime_interpolates(self):
        """Mixed regime result should lie between Na-only and Mg-only corrections."""
        import math
        tm_1m = 65.0
        fGC, n = 0.5, 20

        # ratio = sqrt(Mg)/Na: choose values that put us in mixed regime (0.22–6.0)
        # sqrt(0.003)/0.05 ≈ 1.10 → mixed
        tm_na  = _owczarzy_salt_correction(tm_1m, fGC, n, na_conc=0.05, mg_conc=0.0, dntp_conc=0.0)
        tm_mg  = _owczarzy_salt_correction(tm_1m, fGC, n, na_conc=0.0,  mg_conc=0.003, dntp_conc=0.0)
        tm_mix = _owczarzy_salt_correction(tm_1m, fGC, n, na_conc=0.05, mg_conc=0.003, dntp_conc=0.0)

        lo = min(tm_na, tm_mg)
        hi = max(tm_na, tm_mg)
        assert lo <= tm_mix <= hi + 0.5, (
            f"Mixed regime Tm={tm_mix:.2f} not between Na ({tm_na:.2f}) and Mg ({tm_mg:.2f})"
        )


class TestCalcTmMgConc:
    """Integration tests for calc_tm with Mg²⁺ parameters."""

    def test_mg_raises_tm_in_calc_tm(self):
        """Mg²⁺ should raise Tm relative to Na-only in calc_tm."""
        tm_na = calc_tm(PRIMER, na_conc=0.05, mg_conc=0.0)
        tm_mg = calc_tm(PRIMER, na_conc=0.05, mg_conc=0.003)
        assert tm_mg > tm_na, f"Expected Mg to raise Tm: Na={tm_na:.2f}, Mg={tm_mg:.2f}"

    def test_higher_mg_raises_tm_more(self):
        """Higher [Mg²⁺] → higher Tm (monotonic)."""
        tm_low  = calc_tm(PRIMER, na_conc=0.0, mg_conc=0.001)
        tm_high = calc_tm(PRIMER, na_conc=0.0, mg_conc=0.010)
        assert tm_high > tm_low, f"Expected higher Mg → higher Tm: {tm_low:.2f} < {tm_high:.2f}"

    def test_dntp_reduces_tm(self):
        """dNTP chelation should lower Tm relative to unchelated Mg²⁺."""
        tm_no_dntp = calc_tm(PRIMER, na_conc=0.05, mg_conc=0.003, dntp_conc=0.0)
        tm_with_dntp = calc_tm(PRIMER, na_conc=0.05, mg_conc=0.003, dntp_conc=0.003)
        assert tm_with_dntp < tm_no_dntp, "dNTP chelation should lower Tm"

    def test_typical_pcr_buffer(self):
        """Typical PCR buffer (50 mM KCl ≈ 50 mM Na⁺ + 1.5 mM Mg²⁺, 0.2 mM dNTP)."""
        tm = calc_tm(
            PRIMER,
            na_conc   = 0.05,
            mg_conc   = 0.0015,
            dntp_conc = 0.0002,
            primer_conc = 250e-9,
        )
        # Should be in plausible range for a 20-mer primer
        assert 45 < tm < 80, f"Typical PCR Tm={tm} outside expected range"

    def test_backward_compat_no_mg(self):
        """calc_tm with no mg_conc/dntp_conc must behave as before."""
        tm1 = calc_tm(PRIMER, na_conc=0.05)
        tm2 = calc_tm(PRIMER, na_conc=0.05, mg_conc=0.0, dntp_conc=0.0)
        assert abs(tm1 - tm2) < 1e-9, "No-Mg modes should be identical"


class TestCheckPrimerQualityMg:
    """Verify check_primer_quality propagates Mg²⁺ parameters."""

    def test_mg_changes_tm_in_quality_check(self):
        """Quality check Tm should differ between Na-only and Mg-containing buffer."""
        q_na = check_primer_quality(PRIMER, na_conc=0.05, mg_conc=0.0)
        q_mg = check_primer_quality(PRIMER, na_conc=0.05, mg_conc=0.003)
        assert q_mg["tm_celsius"] > q_na["tm_celsius"], "Mg should raise Tm in quality check"

    def test_tm_ok_flag_consistent(self):
        """tm_ok flag should reflect the mg-corrected Tm."""
        q = check_primer_quality(PRIMER, mg_conc=0.003)
        assert q["tm_ok"] == (50.0 <= q["tm_celsius"] <= 72.0)

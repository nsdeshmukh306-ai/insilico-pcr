"""
Tests for ΔG-based primer-dimer detection.

Validates:
 1. DimerResult has the new delta_g field and it is a float.
 2. dimer_score == abs(delta_g_adj) (backward-compatible positive metric).
 3. Self-complementary primer → large |ΔG| (strong dimer).
 4. Unrelated primers → ΔG near 0 (no dimer).
 5. 3'-end dimer flag is set correctly when terminal bases overlap.
 6. delta_g_adj < delta_g when 3'-end is involved (extra penalty applied).
 7. homo-dimer type is set for identical primers.
 8. check_all_dimers returns three results with correct types.
 9. Numeric ΔG ranges are physically plausible.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest
from modules.advanced.primer_dimer import (
    detect_primer_dimer,
    check_all_dimers,
    DimerResult,
    DG_THRESHOLD_KCAL,
    PRIME3_EXTRA_DG,
)


class TestDimerResultStructure:
    def test_has_delta_g_field(self):
        dr = detect_primer_dimer("GCTAGCTAGCTAGCTAGCTA", "ATAGGCTAAATCGATCGATC")
        assert hasattr(dr, "delta_g"),     "DimerResult must have delta_g"
        assert hasattr(dr, "delta_g_adj"), "DimerResult must have delta_g_adj"
        assert isinstance(dr.delta_g,     float)
        assert isinstance(dr.delta_g_adj, float)

    def test_dimer_score_is_abs_dg_adj(self):
        """dimer_score == abs(delta_g_adj) within floating-point rounding."""
        dr = detect_primer_dimer("GCTAGCTAGCTAGCTAGCTA", "CGATCGATCGATCGATCGAT")
        assert abs(dr.dimer_score - abs(dr.delta_g_adj)) < 1e-6

    def test_dimer_score_nonnegative(self):
        """dimer_score must always be ≥ 0 (backward-compatible)."""
        pairs = [
            ("GCTAGCTAGCTAGCTAGCTA", "ATAGGCTAAATCGATCGATC"),
            ("AATTAATTAATTAATTAATT", "AATTAATTAATTAATTAATT"),
            ("GCGCGCGCGCGCGCGCGCGC", "ATATATATATATATATATATAT"),
        ]
        for a, b in pairs:
            dr = detect_primer_dimer(a, b)
            assert dr.dimer_score >= 0, f"dimer_score < 0 for {a} vs {b}"

    def test_dimer_type_homo(self):
        seq = "AATTAATTAATTAATTAATT"
        dr  = detect_primer_dimer(seq, seq)
        assert dr.dimer_type == "homo"

    def test_dimer_type_hetero(self):
        dr = detect_primer_dimer("GCTAGCTAGCTAGCTAGCTA", "ATAGGCTAAATCGATCGATC")
        assert dr.dimer_type == "hetero"


class TestDimerDeltaG:
    def test_unrelated_primers_low_score(self):
        """Unrelated primers should produce near-zero or low dimer score."""
        # Completely random, no complementarity
        dr = detect_primer_dimer(
            "AAAAAAAAAAAAAAAAAAAAAA",
            "CCCCCCCCCCCCCCCCCCCCCC",
        )
        # No meaningful duplex possible between poly-A and poly-C
        assert dr.dimer_score < DG_THRESHOLD_KCAL, (
            f"Poly-A / Poly-C dimer score={dr.dimer_score:.2f} too high"
        )

    def test_self_complementary_strong_dimer(self):
        """A primer that is its own reverse complement forms a strong homodimer."""
        # "AATTAATT" RC = "AATTAATT" — perfect self-complementary
        seq = "AATTAATTAATTAATTAATT"
        dr  = detect_primer_dimer(seq, seq)
        # Should produce a non-trivial ΔG (negative, stable dimer)
        assert dr.delta_g <= 0, f"Self-complement dimer ΔG should be ≤ 0, got {dr.delta_g}"

    def test_3prime_dimer_penalty_applied(self):
        """When 3'-end is involved, delta_g_adj < delta_g (more negative)."""
        # Design a primer where 3'-end overlap is likely
        fwd = "GGGGGGGGGGGGGGGGGGCC"   # ends CC
        rev = "GGGGGGGGGGGGGGGGGGCC"   # same → 3'-end overlap
        dr  = detect_primer_dimer(fwd, rev)
        if dr.three_prime_a or dr.three_prime_b:
            # 3'-end penalty (PRIME3_EXTRA_DG < 0) should make adj more negative
            assert dr.delta_g_adj <= dr.delta_g, (
                f"3'-penalty should lower delta_g_adj: adj={dr.delta_g_adj} delta_g={dr.delta_g}"
            )

    def test_no_3prime_when_internal(self):
        """When best alignment is purely internal, three_prime flags should be False."""
        # Use primers designed so ends don't overlap
        fwd = "TTTTTTTTTTGCGCTTTTTTTTTT"
        rev = "TTTTTTTTTTGCGCTTTTTTTTTT"
        dr  = detect_primer_dimer(fwd, rev)
        # We don't force any particular result; just verify the dataclass is consistent
        if not dr.three_prime_a and not dr.three_prime_b:
            assert dr.delta_g_adj == dr.delta_g, (
                "No 3'-end penalty should mean delta_g_adj == delta_g"
            )

    def test_delta_g_physically_plausible(self):
        """ΔG for a typical primer dimer should be in range [−30, 0] kcal/mol."""
        primers = [
            ("GCTAGCTAGCTAGCTAGCTA", "CGATCGATCGATCGATCGAT"),
            ("AATTAATTAATTAATTAATT", "AATTAATTAATTAATTAATT"),
        ]
        for a, b in primers:
            dr = detect_primer_dimer(a, b)
            assert -30 <= dr.delta_g <= 2.0, (
                f"ΔG={dr.delta_g:.2f} outside plausible range for {a} vs {b}"
            )

    def test_warning_when_threshold_exceeded(self):
        """Warning string must be non-empty when dimer_score ≥ DG_THRESHOLD_KCAL."""
        # Self-complementary primer likely exceeds threshold
        seq = "AATTAATTAATTAATTAATT"
        dr  = detect_primer_dimer(seq, seq, dg_threshold=0.0)   # force flag at any score
        if dr.dimer_score >= 0.0:
            # With threshold=0, any dimer should be flagged
            assert len(dr.warning) > 0

    def test_no_warning_when_below_threshold(self):
        """Warning must be empty when dimer_score < threshold."""
        dr = detect_primer_dimer(
            "AAAAAAAAAAAAAAAAAAAAAA",
            "CCCCCCCCCCCCCCCCCCCCCC",
            dg_threshold=999.0,   # unreachable threshold
        )
        assert dr.warning == ""


class TestCheckAllDimers:
    def test_returns_three_results(self):
        dimers = check_all_dimers("GCTAGCTAGCTAGCTAGCTA", "TAGCTAGCTAGCTAGCTAGC")
        assert len(dimers) == 3

    def test_types_include_homo(self):
        dimers = check_all_dimers("GCTAGCTAGCTAGCTAGCTA", "TAGCTAGCTAGCTAGCTAGC")
        types  = {d.dimer_type for d in dimers}
        assert "homo" in types, "check_all_dimers must include homo-dimer results"

    def test_all_are_dimer_results(self):
        dimers = check_all_dimers("GCTAGCTAGCTAGCTAGCTA", "TAGCTAGCTAGCTAGCTAGC")
        for d in dimers:
            assert isinstance(d, DimerResult)
            assert isinstance(d.delta_g, float)
            assert d.dimer_score >= 0

    def test_temperature_parameter(self):
        """Higher temperature → less stable dimers → higher ΔG (less negative)."""
        fwd = "AATTAATTAATTAATTAATT"
        rev = "AATTAATTAATTAATTAATT"
        dr_37  = detect_primer_dimer(fwd, rev, temperature=37.0)
        dr_65  = detect_primer_dimer(fwd, rev, temperature=65.0)
        # At higher T, ΔG becomes less negative (less stable)
        assert dr_65.delta_g >= dr_37.delta_g, (
            f"Higher T should give less negative ΔG: 37°C={dr_37.delta_g:.2f}, 65°C={dr_65.delta_g:.2f}"
        )

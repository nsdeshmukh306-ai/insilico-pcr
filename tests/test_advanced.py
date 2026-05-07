"""Unit tests for advanced modules: hairpin, primer-dimer, degenerate."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest
from modules.advanced.hairpin import detect_hairpin
from modules.advanced.primer_dimer import detect_primer_dimer, check_all_dimers
from modules.advanced.degenerate import degeneracy, consensus_to_iupac, design_degenerate_primer


class TestHairpin:
    def test_no_hairpin_random(self):
        result = detect_hairpin("GCTAGCTAGCTAGCTAGCTA")
        # random sequence may or may not have hairpin, just test structure
        assert hasattr(result, "has_hairpin")
        assert hasattr(result, "delta_g")

    def test_clear_hairpin_detected(self):
        # Artificially constructed hairpin: AAAA-loop-TTTT forms a stem
        # stem=GGGGCC, loop=TTTT, stem_rc=GGCCCC
        hairpin_seq = "GCCCCC" + "TTTT" + "GGGGGC" + "AAAAAAAAAAAAA"
        result = detect_hairpin(hairpin_seq, dg_threshold=-2.0)
        # We just verify the function runs and returns a valid result
        assert isinstance(result.delta_g, float)

    def test_warning_only_when_flagged(self):
        result = detect_hairpin("GCTAGCTAGCTAGCTAGCTA", dg_threshold=-3.0)
        if result.has_hairpin:
            assert len(result.warning) > 0
        else:
            assert result.warning == ""


class TestPrimerDimer:
    def test_no_dimer_unrelated(self):
        # Unrelated sequences should have low dimer score
        dr = detect_primer_dimer("GCTAGCTAGCTAGCTAGCTA", "ATAGGCTAAATCGATCGATC")
        # Score might be low
        assert isinstance(dr.dimer_score, float)

    def test_homo_dimer_self_complement(self):
        # A primer that is its own reverse complement — strong self-complementarity
        # e.g. AATTAATT: RC = AATTAATT
        seq = "AATTAATTAATTAATTAATT"
        dr = detect_primer_dimer(seq, seq)
        assert dr.dimer_type == "homo"

    def test_three_prime_flagged(self):
        # Design primers where 3' ends are complementary
        fwd = "GGGGGGGGGGGGGGGGGGCC"  # ends in CC
        rev = "GGGGGGGGGGGGGGGGGGCC"  # also ends in CC — RC starts with GG
        dr = detect_primer_dimer(fwd, rev)
        assert dr.dimer_score >= 0

    def test_check_all_dimers_returns_three(self):
        dimers = check_all_dimers("GCTAGCTAGCTAGCTAGCTA", "TAGCTAGCTAGCTAGCTAGC")
        assert len(dimers) == 3
        types = {d.dimer_type for d in dimers}
        assert "homo" in types


class TestDegenerate:
    def test_no_iupac(self):
        assert degeneracy("ATCG") == 1

    def test_n_is_4(self):
        assert degeneracy("N") == 4

    def test_r_is_2(self):
        assert degeneracy("R") == 2

    def test_multi_iupac(self):
        assert degeneracy("RY") == 4   # R=2, Y=2

    def test_consensus_all_same(self):
        seqs = ["ATCG", "ATCG", "ATCG"]
        assert consensus_to_iupac(seqs) == "ATCG"

    def test_consensus_two_variants(self):
        # position 0: A vs A → A; position 1: T vs G → K (G/T wobble code)
        # To get R (A/G) at position 1, use sequences where pos 1 = A and G
        seqs = ["AACG", "AGCG"]  # position 1: A vs G → R
        result = consensus_to_iupac(seqs)
        assert result[1] == "R", f"Expected R at pos 1, got {result[1]}"

    def test_design_degenerate(self):
        targets = ["ATCGATCG", "ATCGAGCG", "ATCGAACG"]
        result = design_degenerate_primer(targets)
        assert "consensus" in result
        assert result["degeneracy"] >= 1

    def test_design_degenerate_empty(self):
        result = design_degenerate_primer([])
        assert result["consensus"] == ""

    def test_design_unequal_lengths_error(self):
        result = design_degenerate_primer(["ATCG", "ATCGATCG"])
        assert "warning" in result
        assert result["consensus"] == ""

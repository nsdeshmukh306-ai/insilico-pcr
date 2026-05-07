"""Unit tests for binding_search module."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest
from modules.binding_search import smith_waterman_align, ungapped_align, find_binding_sites


class TestSmithWaterman:
    def test_identical_sequences(self):
        score, aq, at, ts, te = smith_waterman_align("ATCG", "ATCG")
        assert score > 0
        assert aq == "ATCG"

    def test_one_mismatch(self):
        score_perfect, _, _, _, _ = smith_waterman_align("ATCG", "ATCG")
        score_mm,      _, _, _, _ = smith_waterman_align("ATCG", "AGCG")  # 1 mismatch
        assert score_mm < score_perfect

    def test_query_in_longer_target(self):
        query  = "ATCG"
        target = "XXXXATCGXXXX"
        score, aq, at, ts, te = smith_waterman_align(query, target)
        assert aq == "ATCG"
        # ts and te should correspond to the ATCG region
        assert target[ts:te] == "ATCG"

    def test_zero_similarity(self):
        score, _, _, _, _ = smith_waterman_align("AAAA", "TTTT")
        # All mismatches — SW can return 0 (stops at zero)
        assert score == 0.0


class TestUngappedAlign:
    def test_perfect_match(self):
        score, mm_pos, gaps = ungapped_align("ATCG", "ATCG")
        assert mm_pos == []
        assert gaps == 0
        assert score > 0

    def test_single_mismatch(self):
        score, mm_pos, gaps = ungapped_align("ATCG", "AGCG")
        assert 1 in mm_pos
        assert len(mm_pos) == 1

    def test_all_mismatches(self):
        _, mm_pos, _ = ungapped_align("AAAA", "TTTT")
        assert len(mm_pos) == 4


class TestFindBindingSites:
    def test_exact_match(self):
        primer = "ATCGATCG"
        genome = "NNNN" + primer + "NNNN"
        sites = find_binding_sites(
            primer_seq         = primer,
            strand             = "+",
            seq_id             = "chr1",
            genome_seq         = genome,
            candidate_starts   = [4],
            primer_len         = len(primer),
            max_mismatches     = 0,
            three_prime_strict = True,
        )
        assert len(sites) == 1
        assert sites[0].mismatch_count == 0
        assert sites[0].start == 4

    def test_too_many_mismatches_filtered(self):
        primer = "ATCGATCG"
        genome = "NNNN" + "TTTTTTTT" + "NNNN"  # 8 mismatches
        sites = find_binding_sites(
            primer_seq       = primer,
            strand           = "+",
            seq_id           = "chr1",
            genome_seq       = genome,
            candidate_starts = [4],
            primer_len       = len(primer),
            max_mismatches   = 2,
            three_prime_strict = False,
        )
        assert len(sites) == 0

    def test_allowed_mismatches_pass(self):
        primer = "ATCGATCG"
        # 2 mismatches
        target = list(primer)
        target[2] = "A"
        target[5] = "T"
        genome = "NNNN" + "".join(target) + "NNNN"
        sites = find_binding_sites(
            primer_seq       = primer,
            strand           = "+",
            seq_id           = "chr1",
            genome_seq       = genome,
            candidate_starts = [4],
            primer_len       = len(primer),
            max_mismatches   = 2,
            three_prime_strict = False,
            use_sw           = False,  # use fast ungapped for determinism
        )
        assert len(sites) >= 1

    def test_three_prime_mismatch_blocked(self):
        primer = "ATCGATCG"
        target = list(primer)
        target[-1] = "T"  # 3'-terminal mismatch
        genome = "NNNN" + "".join(target) + "NNNN"
        sites = find_binding_sites(
            primer_seq       = primer,
            strand           = "+",
            seq_id           = "chr1",
            genome_seq       = genome,
            candidate_starts = [4],
            primer_len       = len(primer),
            max_mismatches   = 3,
            three_prime_strict = True,
            use_sw           = False,
        )
        assert len(sites) == 0

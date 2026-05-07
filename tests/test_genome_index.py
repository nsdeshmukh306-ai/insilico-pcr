"""Unit tests for genome_index module."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from modules.genome_index import (
    build_kmer_index,
    lookup_seeds,
    get_sequence_at,
    records_to_map,
    _rc,
)


def make_record(seq_id, seq):
    return SeqRecord(Seq(seq.upper()), id=seq_id, description="")


class TestRC:
    def test_simple(self):
        assert _rc("ATCG") == "CGAT"

    def test_palindrome(self):
        assert _rc("AATT") == "AATT"

    def test_all_gc(self):
        assert _rc("GGGG") == "CCCC"


class TestBuildKmerIndex:
    def test_basic_index(self):
        rec = make_record("chr1", "ATCGATCGATCGATCG")
        idx = build_kmer_index([rec], k=4)
        assert "ATCG" in idx

    def test_index_has_both_strands(self):
        rec = make_record("chr1", "ATCGATCG")
        idx = build_kmer_index([rec], k=4)
        # RC of ATCG is CGAT, should also be in index
        assert "CGAT" in idx

    def test_n_bases_excluded(self):
        rec = make_record("chr1", "ATCNNNNATCG")
        idx = build_kmer_index([rec], k=4)
        # k-mers spanning N should not be indexed
        for kmer in idx:
            assert "N" not in kmer

    def test_multiple_records(self):
        r1 = make_record("chr1", "ATCGATCG")
        r2 = make_record("chr2", "GCTAGCTA")
        idx = build_kmer_index([r1, r2], k=4)
        # Both records should be represented
        seqs_in_index = set(pos[0] for kmer_hits in idx.values() for pos in kmer_hits)
        assert "chr1" in seqs_in_index
        assert "chr2" in seqs_in_index


class TestLookupSeeds:
    def test_exact_hit(self):
        seq = "AAAAAATCGATCGTTTTTT"
        rec = make_record("chr1", seq)
        idx = build_kmer_index([rec], k=4)
        hits = lookup_seeds("ATCGATCG", idx, k=4, strand_filter="+")
        assert any(h[0] == "chr1" for h in hits)

    def test_no_hit(self):
        seq = "AAAAAAAAAAAAAAAAAAA"
        rec = make_record("chr1", seq)
        idx = build_kmer_index([rec], k=4)
        hits = lookup_seeds("GCGCGCGC", idx, k=4)
        # GCGC might not appear in all-A sequence
        assert all(h[0] != "chr1" or seq.find("GCGC") >= 0 for h in hits)


class TestGetSequenceAt:
    def test_forward(self):
        gmap = {"chr1": "ATCGATCG"}
        sub = get_sequence_at(gmap, "chr1", 0, 4, "+")
        assert sub == "ATCG"

    def test_reverse(self):
        gmap = {"chr1": "ATCGATCG"}
        sub = get_sequence_at(gmap, "chr1", 0, 4, "-")
        assert sub == _rc("ATCG")

    def test_out_of_bounds(self):
        gmap = {"chr1": "ATCG"}
        sub = get_sequence_at(gmap, "chr1", 10, 20, "+")
        assert sub == ""

    def test_missing_id(self):
        gmap = {"chr1": "ATCG"}
        sub = get_sequence_at(gmap, "chr2", 0, 4, "+")
        assert sub == ""

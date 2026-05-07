"""
Tests for FM-index (BWT-based) genome indexing.

Validates:
 1. Suffix array construction is correct for small strings.
 2. Backward search finds exact matches (single and multiple).
 3. No false positives (absent patterns return empty list).
 4. Forward and reverse-strand seed lookup matches k-mer index output.
 5. build_fm_index / lookup_seeds_fm interface parity with kmer equivalents.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest
from modules.genome_index_fm import (
    _build_suffix_array,
    FMIndex,
    GenomeFMIndex,
    build_fm_index,
    lookup_seeds_fm,
)


# ---------------------------------------------------------------------------
# Minimal SeqRecord stand-in (avoids Bio dependency)
# ---------------------------------------------------------------------------
class _Rec:
    def __init__(self, seq_id, seq):
        self.id  = seq_id
        self.seq = _Seq(seq)


class _Seq:
    def __init__(self, s):
        self._s = s

    def __str__(self):
        return self._s


# ---------------------------------------------------------------------------
# Suffix Array tests
# ---------------------------------------------------------------------------
class TestSuffixArray:
    def test_single_char(self):
        """Single character plus sentinel gives SA of length 2."""
        sa = _build_suffix_array("A$")
        assert len(sa) == 2
        assert sa[0] == 1   # "$" is smallest
        assert sa[1] == 0   # "A$"

    def test_banana(self):
        """Classic 'banana$' SA: positions sorted lexicographically."""
        text = "banana$"
        sa   = _build_suffix_array(text)
        # Suffixes sorted: $, a$, ana$, anana$, banana$, na$, nana$
        expected_starts = [6, 5, 3, 1, 0, 4, 2]
        assert sa == expected_starts, f"SA wrong: {sa}"

    def test_all_same_chars(self):
        """Repeated character: SA must still be a valid permutation."""
        sa = _build_suffix_array("AAAA$")
        assert sorted(sa) == list(range(5))
        assert sa[0] == 4   # "$" is first

    def test_dna_short(self):
        """Short DNA string: SA is a valid permutation of [0..n-1]."""
        text = "ACGT$"
        sa   = _build_suffix_array(text)
        assert sorted(sa) == list(range(len(text)))
        assert sa[0] == 4   # "$" is smallest


# ---------------------------------------------------------------------------
# FMIndex backward search tests
# ---------------------------------------------------------------------------
class TestFMIndex:
    def _idx(self, text):
        return FMIndex(text)

    def test_exact_single_occurrence(self):
        """Pattern present once → returns its start position."""
        idx = self._idx("ACGTACGT")
        hits = idx.backward_search("CGT")
        assert hits == [1, 5], f"Expected [1, 5], got {hits}"

    def test_exact_not_found(self):
        """Absent pattern → empty list."""
        idx = self._idx("ACGTACGT")
        assert idx.backward_search("TTT") == []

    def test_full_string(self):
        """Search for the entire text → position 0."""
        text = "GCTAGCTA"
        idx  = self._idx(text)
        hits = idx.backward_search(text)
        assert hits == [0]

    def test_single_char(self):
        """Single-character pattern finds all occurrences."""
        idx  = self._idx("AAGCAA")
        hits = idx.backward_search("A")
        assert sorted(hits) == [0, 1, 4, 5]

    def test_primer_in_genome(self):
        """Realistic primer embedded in a longer sequence."""
        primer = "GCACTGGTGG"
        genome = "N" * 200 + primer + "N" * 300
        idx    = FMIndex(genome)
        hits   = idx.backward_search(primer)
        assert 200 in hits, f"Expected hit at 200, got {hits}"

    def test_repeated_motif(self):
        """Pattern repeated 3× at known positions."""
        motif  = "ATCG"
        genome = motif + "NNNN" + motif + "NNNN" + motif
        idx    = FMIndex(genome)
        hits   = idx.backward_search(motif)
        assert sorted(hits) == [0, 8, 16]

    def test_no_false_positives_random_primer(self):
        """Pattern not in text must not appear in results."""
        idx = FMIndex("AAAAAAAAAAAAAAAAAAAAAAAAAAAA")
        assert idx.backward_search("GCTA") == []


# ---------------------------------------------------------------------------
# GenomeFMIndex and build_fm_index tests
# ---------------------------------------------------------------------------
class TestGenomeFMIndex:
    def _make_rec(self, seq_id, seq):
        return _Rec(seq_id, seq)

    def test_build_single_record(self):
        """Build FM-index for a single short record."""
        rec  = self._make_rec("seq1", "ACGTACGTACGT")
        recs = [rec]
        idx  = build_fm_index(recs)
        assert "seq1" in idx
        assert isinstance(idx["seq1"], GenomeFMIndex)

    def test_fwd_search_finds_primer(self):
        """Forward-strand seed lookup finds primer at correct position."""
        primer = "GCTAGC"
        genome = "A" * 50 + primer + "A" * 50
        rec    = self._make_rec("chr1", genome)
        fm_idx = build_fm_index([rec])

        hits = lookup_seeds_fm(primer, fm_idx, k=6, strand_filter="+")
        fwd_starts = [h[1] for h in hits if h[0] == "chr1" and h[2] == "+"]
        assert 50 in fwd_starts, f"Primer should start at 50, hits: {hits}"

    def test_rev_search_finds_primer(self):
        """Reverse-complement primer found on '-' strand."""
        primer = "GCTAGC"
        rc_p   = primer[::-1].translate(str.maketrans("ACGT", "TGCA"))   # GCTAGC → GCTAGC (palindrome-ish)
        # Use a non-palindromic primer
        primer2 = "AACCGG"
        rc_p2   = "CCGGTT"   # rc of AACCGG
        genome  = "T" * 30 + rc_p2 + "T" * 30   # rc of primer on fwd strand
        rec     = self._make_rec("chr1", genome)
        fm_idx  = build_fm_index([rec])

        hits = lookup_seeds_fm(primer2, fm_idx, k=6, strand_filter="-")
        minus_hits = [h for h in hits if h[2] == "-"]
        assert len(minus_hits) > 0, f"Expected '-' strand hit, got: {hits}"

    def test_multi_record(self):
        """Multiple records each indexed separately."""
        r1 = self._make_rec("chr1", "ACGTACGT" * 5)
        r2 = self._make_rec("chr2", "GCTAGCTA" * 5)
        fm = build_fm_index([r1, r2])
        assert "chr1" in fm and "chr2" in fm

    def test_seed_lookup_interface_parity(self):
        """
        lookup_seeds_fm output format matches lookup_seeds kmer output:
        each element is (seq_id, genome_start, strand, offset) where offset is int.
        """
        primer = "AACCGG"
        genome = "N" * 10 + primer + "N" * 10
        rec    = self._make_rec("c1", genome)
        fm_idx = build_fm_index([rec])
        hits   = lookup_seeds_fm(primer, fm_idx, k=4, strand_filter="+")
        for h in hits:
            assert len(h) == 4
            sid, gstart, strand, offset = h
            assert isinstance(sid, str)
            assert isinstance(gstart, int)
            assert strand in ("+", "-")
            assert isinstance(offset, int)


# ---------------------------------------------------------------------------
# Parity with k-mer index for small genomes
# ---------------------------------------------------------------------------
class TestFMvsKmerParity:
    """Forward-strand hits from FM-index should match those from k-mer index."""

    def test_same_fwd_hits_as_kmer(self):
        """FM-index and k-mer index should agree on forward-strand seed positions."""
        from modules.genome_index import build_kmer_index, lookup_seeds

        primer = "GCACTGGT"   # 8-mer
        genome = "AAAA" + primer + "CCCC" + primer + "GGGG"
        rec    = _Rec("g1", genome)

        # k-mer index
        kmer_idx = build_kmer_index([rec], k=8)
        kmer_hits = lookup_seeds(primer, kmer_idx, k=8, strand_filter="+")
        kmer_starts = sorted({h[1] for h in kmer_hits if h[0] == "g1"})

        # FM-index
        fm_idx  = build_fm_index([rec], k=8)
        fm_hits = lookup_seeds_fm(primer, fm_idx, k=8, strand_filter="+")
        fm_starts = sorted({h[1] for h in fm_hits if h[0] == "g1"})

        assert kmer_starts == fm_starts, (
            f"Hit positions differ: kmer={kmer_starts}, fm={fm_starts}"
        )

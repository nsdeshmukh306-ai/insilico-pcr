"""
FM-Index (Burrows-Wheeler Transform) Genome Index
==================================================
Provides sub-linear search for exact primer k-mer seeds on large genomes.

Architecture
------------
1. Suffix Array (SA) — built with prefix-doubling in O(n log² n) time.
2. Burrows-Wheeler Transform (BWT) — BWT[i] = text[SA[i]-1].
3. C array — cumulative character counts (first-column of the BWT matrix).
4. Sampled rank checkpoints — every CHKPT positions, storing per-character
   occurrence counts. Rank queries run in O(CHKPT) rather than O(n).
5. Backward search — O(m × CHKPT) per query of length m.

One FMIndex is built per strand per genome record.  Two indexes are held:
  fwd_index : built on the forward strand sequence
  rev_index : built on the reverse-complement sequence

Search returns start positions on the FORWARD strand (consistent with the
k-mer index interface).

Performance note
----------------
The suffix array is built in pure Python (prefix doubling, O(n log² n)).
For genomes ≤ 10 Mbp this is fast enough (seconds).  For 10–50 Mbp it may
take 1–5 minutes; for larger genomes consider installing a compiled SA library
and replacing `_build_suffix_array`.  The auto-switch in genome_index.py is
set at 50 Mbp.

Interface
---------
The module exposes `build_fm_index` and `lookup_seeds_fm`, which replicate
the signatures of the k-mer equivalents in genome_index.py so callers are
unaffected.
"""

import logging
from typing import Dict, List, Optional, Tuple

log = logging.getLogger(__name__)

SENTINEL    = "$"       # Appended to text; must be lex-smallest char
CHKPT       = 64        # Rank checkpoint interval

_COMP = str.maketrans("ACGTN", "TGCAN")


def _rc(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


# ---------------------------------------------------------------------------
# Suffix Array (prefix doubling)
# ---------------------------------------------------------------------------
def _build_suffix_array(text: str) -> List[int]:
    """
    Build suffix array by prefix doubling.  O(n log² n) time and O(n) space.

    Returns list of length n where SA[i] is the start of the i-th
    lexicographically smallest suffix of `text`.
    """
    n = len(text)
    if n == 0:
        return []

    # Initial ranks: ordinal of each character
    rank = [ord(c) for c in text]
    sa   = list(range(n))
    k    = 1

    while k < n:
        # Sort SA by (rank[i], rank[i+k] or -1 for out-of-bounds)
        rank_snapshot = rank[:]

        def _key(i: int, r=rank_snapshot, step=k, length=n):
            return (r[i], r[i + step] if i + step < length else -1)

        sa.sort(key=_key)

        # Update ranks from sorted order
        new_rank = [0] * n
        new_rank[sa[0]] = 0
        for j in range(1, n):
            p, c = sa[j - 1], sa[j]
            rp = (rank_snapshot[p], rank_snapshot[p + k] if p + k < n else -1)
            rc = (rank_snapshot[c], rank_snapshot[c + k] if c + k < n else -1)
            new_rank[c] = new_rank[p] + (0 if rc == rp else 1)
        rank = new_rank

        if rank[sa[-1]] == n - 1:
            break  # All suffixes have unique rank — done
        k *= 2

    return sa


# ---------------------------------------------------------------------------
# FM-Index class
# ---------------------------------------------------------------------------
class FMIndex:
    """
    FM-Index over a single DNA sequence string.

    Parameters
    ----------
    text : str
        Forward-strand genome sequence (uppercase ACGTN, no whitespace).
    """

    def __init__(self, text: str) -> None:
        self.text = text
        self.n    = len(text)
        self._build()

    def _build(self) -> None:
        n    = self.n
        full = self.text + SENTINEL   # SENTINEL is lex-smallest

        sa = _build_suffix_array(full)
        self.sa = sa

        # BWT: character immediately before SA[i] in the text (circular = sentinel wraps)
        bwt = [full[sa[i] - 1] if sa[i] > 0 else SENTINEL for i in range(len(sa))]
        self.bwt = bwt

        # Alphabet (sorted, SENTINEL first)
        alph = sorted(set(full))
        self.alph = alph

        # C array: cumulative counts (how many chars in full < c)
        char_count: Dict[str, int] = {c: 0 for c in alph}
        for ch in full:
            char_count[ch] += 1
        self.C: Dict[str, int] = {}
        cum = 0
        for c in alph:
            self.C[c] = cum
            cum += char_count[c]

        # Rank checkpoints: occ_ckpt[c][b] = count of c in bwt[0 : b*CHKPT]
        nbwt  = len(bwt)
        nblks = (nbwt + CHKPT - 1) // CHKPT + 1
        self.occ_ckpt: Dict[str, List[int]] = {c: [0] * nblks for c in alph}
        running: Dict[str, int] = {c: 0 for c in alph}
        for i, ch in enumerate(bwt):
            running[ch] += 1
            if (i + 1) % CHKPT == 0:
                blk = (i + 1) // CHKPT
                for c in alph:
                    self.occ_ckpt[c][blk] = running[c]
        # Final block
        last_blk = (nbwt + CHKPT - 1) // CHKPT
        for c in alph:
            self.occ_ckpt[c][last_blk] = running[c]

    def rank(self, c: str, i: int) -> int:
        """
        Number of occurrences of character c in BWT[0..i] (0-based, inclusive).
        Returns 0 if i < 0 or c not in alphabet.
        """
        if i < 0 or c not in self.occ_ckpt:
            return 0
        blk   = i // CHKPT
        count = self.occ_ckpt[c][blk]
        start = blk * CHKPT
        for j in range(start, i + 1):
            if self.bwt[j] == c:
                count += 1
        return count

    def backward_search(self, pattern: str) -> List[int]:
        """
        Find all starting positions of `pattern` in the indexed text.
        Returns a sorted list of 0-based start positions.
        """
        lo = 0
        hi = len(self.bwt) - 1

        for i in range(len(pattern) - 1, -1, -1):
            c = pattern[i]
            if c not in self.C:
                return []
            lo = self.C[c] + self.rank(c, lo - 1)
            hi = self.C[c] + self.rank(c, hi) - 1
            if lo > hi:
                return []

        return sorted(self.sa[lo: hi + 1])


# ---------------------------------------------------------------------------
# Per-record FM-index holder
# ---------------------------------------------------------------------------
class GenomeFMIndex:
    """Holds forward + reverse FM-indexes for a single genome record."""

    def __init__(self, seq_id: str, seq: str) -> None:
        self.seq_id = seq_id
        self.seq    = seq
        n = len(seq)
        log.info("Building FM-index for %s (n=%d) …", seq_id, n)
        self.fwd = FMIndex(seq)
        self.rev = FMIndex(_rc(seq))
        log.info("FM-index ready for %s.", seq_id)


# ---------------------------------------------------------------------------
# Build function (mirrors build_kmer_index interface)
# ---------------------------------------------------------------------------
def build_fm_index(records: list, k: int = 8) -> Dict[str, "GenomeFMIndex"]:
    """
    Build FM-indexes for all genome records.

    Parameters
    ----------
    records : list of SeqRecord
    k       : seed length (kept for API parity; not used internally)

    Returns
    -------
    dict mapping seq_id → GenomeFMIndex
    """
    result: Dict[str, GenomeFMIndex] = {}
    for rec in records:
        sid = rec.id
        seq = str(rec.seq).upper()
        result[sid] = GenomeFMIndex(sid, seq)
    return result


# ---------------------------------------------------------------------------
# Seed lookup (mirrors lookup_seeds interface)
# ---------------------------------------------------------------------------
def lookup_seeds_fm(
    primer: str,
    fm_indexes: Dict[str, GenomeFMIndex],
    k: int = 8,
    step: int = 1,
    strand_filter: str = "both",
) -> List[Tuple[str, int, str, int]]:
    """
    Find all genome positions where any k-mer from `primer` matches via FM-index.

    Returns
    -------
    List of (seq_id, genome_start_fwd, strand, primer_offset) —
    identical signature to `lookup_seeds` in genome_index.py.
    """
    primer = primer.upper()
    m      = len(primer)
    hits: List[Tuple[str, int, str, int]] = []
    seen  = set()

    for sid, gfm in fm_indexes.items():
        seq_len = gfm.fwd.n

        for offset in range(0, m - k + 1, step):
            kmer = primer[offset: offset + k]
            if "N" in kmer:
                continue

            # --- Forward strand ---
            if strand_filter in ("+", "both"):
                for gpos in gfm.fwd.backward_search(kmer):
                    primer_start = gpos - offset
                    key = (sid, primer_start, "+")
                    if key not in seen:
                        seen.add(key)
                        hits.append((sid, primer_start, "+", offset))

            # --- Reverse-complement strand ---
            if strand_filter in ("-", "both"):
                # RC index is built on rc(fwd_seq).
                # A hit at position p in rc(seq) corresponds to fwd-strand
                # position: seq_len - p - k
                for rc_pos in gfm.rev.backward_search(kmer):
                    fwd_pos = seq_len - rc_pos - k
                    primer_start = fwd_pos - (m - k - offset)
                    key = (sid, primer_start, "-")
                    if key not in seen:
                        seen.add(key)
                        hits.append((sid, primer_start, "-", offset))

    return hits

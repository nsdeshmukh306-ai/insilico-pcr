"""
Genome Indexing (Layer 3)
=========================
Builds a positional k-mer index of all genome sequences for efficient
seed-based primer binding site discovery.

Design
------
Index structure:
    kmer_index[kmer_string] = [(seq_id, start_pos, strand), ...]

  - 'strand' is '+' (forward) or '-' (reverse complement)
  - start_pos is 0-based, refers to the forward strand coordinate
  - For the reverse strand, start_pos is the position where the k-mer starts
    on the forward strand (counting from the 5' end of the forward strand)

Memory note:
  - At k=8: ~65,000 possible kmers × ~1 match per 350 bp → ~3M entries per
    100 Mbp genome. Manageable in RAM.
  - At k=12: fewer false positives, ~260,000 possible kmers. Recommended for
    genomes >50 Mbp. This module uses k=8 by default; change via PCRParams.

Index is built once per genome and can be reused for all primer pairs.
"""

import logging
from collections import defaultdict
from typing import Dict, List, Tuple

from Bio.SeqRecord import SeqRecord

log = logging.getLogger(__name__)

# Type alias: position tuple
PosEntry = Tuple[str, int, str]   # (seq_id, start_pos_fwd, strand)
KmerIndex = Dict[str, List[PosEntry]]


# ---------------------------------------------------------------------------
# Reverse complement (local, no import loop)
# ---------------------------------------------------------------------------
_COMP = str.maketrans("ACGTN", "TGCAN")


def _rc(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


# ---------------------------------------------------------------------------
# Index builder
# ---------------------------------------------------------------------------
def build_kmer_index(
    records: List[SeqRecord],
    k: int = 8,
) -> KmerIndex:
    """
    Build a positional k-mer index from a list of SeqRecords.

    Parameters
    ----------
    records : list of SeqRecord
        Preprocessed (uppercase, sanitised) genome/transcript sequences.
    k : int
        Seed length (default 8). Increase for specificity on large genomes.

    Returns
    -------
    KmerIndex : dict mapping kmer → list of (seq_id, start_pos, strand)
    """
    index: KmerIndex = defaultdict(list)
    total_kmers = 0

    for rec in records:
        seq_str = str(rec.seq).upper()
        seq_id  = rec.id
        n       = len(seq_str)

        # Forward strand
        for i in range(n - k + 1):
            kmer = seq_str[i : i + k]
            if "N" not in kmer:   # skip Ns to avoid spurious hits
                index[kmer].append((seq_id, i, "+"))
                total_kmers += 1

        # Reverse complement strand
        rc_seq = _rc(seq_str)
        for i in range(n - k + 1):
            kmer = rc_seq[i : i + k]
            if "N" not in kmer:
                # Convert rc position back to forward-strand coordinates
                fwd_pos = n - i - k
                index[kmer].append((seq_id, fwd_pos, "-"))
                total_kmers += 1

    log.info(
        "Built k=%d index: %d unique kmers, %d total entries from %d records.",
        k, len(index), total_kmers, len(records),
    )
    return dict(index)


# ---------------------------------------------------------------------------
# Seed lookup
# ---------------------------------------------------------------------------
def lookup_seeds(
    primer: str,
    index: KmerIndex,
    k: int = 8,
    step: int = 1,
    strand_filter: str = "both",   # "+", "-", or "both"
) -> List[Tuple[str, int, str, int]]:
    """
    Find all genome positions where any k-mer from `primer` matches the index.

    Returns
    -------
    List of (seq_id, genome_start_fwd, strand, primer_offset)
      - genome_start_fwd : 0-based position on the forward strand where the
                           k-mer starts.
      - primer_offset    : 0-based position within the primer where the seed starts.
    The caller extends from these seed positions using SW alignment.
    """
    hits: List[Tuple[str, int, str, int]] = []
    seen = set()   # deduplicate (same genome locus, different seeds can hit it)

    primer = primer.upper()
    n = len(primer)

    for offset in range(0, n - k + 1, step):
        kmer = primer[offset : offset + k]
        if "N" in kmer:
            continue

        for (seq_id, gpos, strand) in index.get(kmer, []):
            if strand_filter != "both" and strand != strand_filter:
                continue

            # Estimate where primer would START on the genome from this seed
            if strand == "+":
                primer_genome_start = gpos - offset
            else:
                # On minus strand: kmer at gpos means primer 5'→3' is on the
                # rc strand; seed is at primer_offset from the 5' end of primer
                # (the rc of the primer).
                primer_genome_start = gpos - (n - k - offset)

            dedup_key = (seq_id, primer_genome_start, strand)
            if dedup_key not in seen:
                seen.add(dedup_key)
                hits.append((seq_id, primer_genome_start, strand, offset))

    return hits


# ---------------------------------------------------------------------------
# Sequence retrieval helper
# ---------------------------------------------------------------------------
def get_sequence_at(
    records_map: Dict[str, str],
    seq_id: str,
    start: int,
    end: int,
    strand: str,
) -> str:
    """
    Retrieve a sub-sequence from a genome.

    Parameters
    ----------
    records_map : dict seq_id → sequence string (forward strand)
    start, end  : 0-based half-open interval on the forward strand
    strand      : '+' or '-'

    Returns the subsequence, reverse-complemented if strand == '-'.
    Returns '' if coordinates are out of bounds.
    """
    seq = records_map.get(seq_id, "")
    if not seq:
        return ""
    n = len(seq)
    s = max(0, start)
    e = min(n, end)
    if s >= e:
        return ""
    sub = seq[s:e]
    if strand == "-":
        return _rc(sub)
    return sub


def records_to_map(records: List[SeqRecord]) -> Dict[str, str]:
    """Convert a list of SeqRecords to {id: sequence_string} for fast lookup."""
    return {rec.id: str(rec.seq) for rec in records}


# ---------------------------------------------------------------------------
# Auto-switching index builder
# ---------------------------------------------------------------------------
LARGE_GENOME_THRESHOLD = 50_000_000   # 50 Mbp — switch to FM-index above this


def build_index(records: List[SeqRecord], k: int = 8):
    """
    Build a genome search index, automatically selecting between k-mer and FM-index.

    For total genome size ≤ LARGE_GENOME_THRESHOLD (50 Mbp):
        Returns a k-mer KmerIndex dict (fast construction, O(1) lookup).
    For total genome size > LARGE_GENOME_THRESHOLD:
        Returns a dict mapping seq_id → GenomeFMIndex (lower memory for large genomes).
        Callers must use lookup_seeds_unified() instead of lookup_seeds().

    Returns
    -------
    (index_obj, index_type)  where index_type is "kmer" or "fm".
    """
    total_bp = sum(len(str(r.seq)) for r in records)
    if total_bp > LARGE_GENOME_THRESHOLD:
        log.info(
            "Genome size %d bp > threshold %d bp — using FM-index.",
            total_bp, LARGE_GENOME_THRESHOLD,
        )
        from .genome_index_fm import build_fm_index
        return build_fm_index(records, k=k), "fm"
    else:
        log.info(
            "Genome size %d bp ≤ threshold %d bp — using k-mer index.",
            total_bp, LARGE_GENOME_THRESHOLD,
        )
        return build_kmer_index(records, k=k), "kmer"


def lookup_seeds_unified(
    primer: str,
    index,
    index_type: str,
    k: int = 8,
    step: int = 1,
    strand_filter: str = "both",
) -> List[Tuple[str, int, str, int]]:
    """
    Dispatch seed lookup to the appropriate backend based on index_type.

    Parameters
    ----------
    primer       : 5'→3' primer sequence
    index        : return value of build_index() — either KmerIndex or FM dict
    index_type   : "kmer" or "fm"
    k, step, strand_filter : passed through to the backend

    Returns
    -------
    List of (seq_id, genome_start_fwd, strand, primer_offset)
    """
    if index_type == "fm":
        from .genome_index_fm import lookup_seeds_fm
        return lookup_seeds_fm(primer, index, k=k, step=step, strand_filter=strand_filter)
    return lookup_seeds(primer, index, k=k, step=step, strand_filter=strand_filter)

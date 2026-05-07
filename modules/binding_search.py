"""
Primer Binding Search Engine (Layer 4)
=======================================
Finds all candidate primer binding sites on a genome using a two-phase strategy:

Phase 1 — Seed
  K-mer index lookup gives candidate positions (fast; O(1) per kmer lookup).

Phase 2 — Extend (Smith-Waterman local alignment)
  Smith–Waterman (SW) local alignment scores each candidate position.
  High gap-open penalty discourages bulge loops (biologically rare in PCR).
  The alignment window is primer_length ± max_mismatches to allow indels.

Why SW and not BWT/Burrows-Wheeler?
  SW is universally accepted, easy to implement correctly without libraries,
  and produces interpretable alignment scores + mismatch maps.
  For primers (≤60 bp), the O(m*n) cost per candidate window is negligible.

Assumptions:
  - Only DNA/DNA duplexes (no RNA).
  - Gaps (insertions/deletions vs template) are penalised heavily — they are
    biologically possible (bulge loops) but uncommon in PCR primer binding.
  - The alignment is semi-global at the primer ends: the full primer must be
    covered (global in query), but the template is local.

Strand convention (IMPORTANT):
  ┌─────────────────────────────────────────────────────────────┐
  │ Forward primer  → binds forward strand  (sense strand)      │
  │ Reverse primer  → binds reverse strand (antisense strand)   │
  │                                                             │
  │ To search:                                                  │
  │   fwd primer as-is       vs fwd strand genome               │
  │   rev primer as-is       vs fwd strand genome BUT we        │
  │     actually look for rc(rev_primer) on the fwd strand,     │
  │     because the rev primer's 3' end faces LEFT on the       │
  │     forward strand.  Equivalently, we search the rev primer │
  │     sequence against the reverse-complement of the genome.  │
  └─────────────────────────────────────────────────────────────┘
  This module uses the second interpretation: all searches are against
  the FORWARD strand; reverse primer is searched as its reverse complement.
"""

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Smith-Waterman scoring constants
# ---------------------------------------------------------------------------
SW_MATCH    =  2
SW_MISMATCH = -1
SW_GAP_OPEN = -5   # Affine gap: high penalty discourages bulge loops
SW_GAP_EXT  = -1


# ---------------------------------------------------------------------------
# Result dataclass
# ---------------------------------------------------------------------------
@dataclass
class BindingSite:
    """
    A candidate primer binding location on the genome.

    Coordinates are always on the FORWARD strand (0-based, half-open).
    strand indicates which template strand the primer binds:
      '+' → primer binds fwd strand (forward primer)
      '-' → primer binds rev strand (reverse primer's binding site)

    The 'primer_seq' stored is the primer 5'→3' sequence.
    The 'aligned_template' is the matching genomic window (fwd strand coords),
    reverse-complemented if strand=='-' so it reads 3'→5' relative to primer.
    """
    seq_id:           str
    start:            int    # 0-based start on fwd strand
    end:              int    # 0-based end (exclusive) on fwd strand
    strand:           str    # '+' or '-'
    primer_seq:       str    # 5'→3' primer sequence
    aligned_template: str    # genomic window matching primer, same orientation
    sw_score:         float  = 0.0
    mismatch_count:   int    = 0
    mismatch_pos:     List[int] = field(default_factory=list)  # 0-based in primer
    gap_count:        int    = 0
    three_prime_mm:   bool   = False   # mismatch at 3'-terminal base


# ---------------------------------------------------------------------------
# Core Smith-Waterman (semi-global: full query coverage)
# ---------------------------------------------------------------------------
def smith_waterman_align(
    query: str,
    target: str,
    match: int    = SW_MATCH,
    mismatch: int = SW_MISMATCH,
    gap_open: int = SW_GAP_OPEN,
    gap_ext: int  = SW_GAP_EXT,
) -> Tuple[float, str, str, int, int]:
    """
    Smith-Waterman local alignment with affine gap penalties.
    Uses numpy for the DP matrix.

    Returns
    -------
    (score, aligned_query, aligned_target, target_start, target_end)
      - score         : alignment score
      - aligned_query : query string with '-' for gaps
      - aligned_target: target string with '-' for gaps
      - target_start  : 0-based start in target (before alignment)
      - target_end    : 0-based end (exclusive) in target
    """
    m = len(query)
    n = len(target)

    # Score matrix H, gap matrices E (horizontal gap), F (vertical gap)
    H = np.zeros((m + 1, n + 1), dtype=np.float32)
    E = np.full((m + 1, n + 1), -np.inf, dtype=np.float32)
    F = np.full((m + 1, n + 1), -np.inf, dtype=np.float32)
    # Traceback
    TB = np.zeros((m + 1, n + 1), dtype=np.int8)
    # TB codes: 0=stop, 1=diag(match/mismatch), 2=up(gap in target), 3=left(gap in query)

    best_score = 0.0
    best_i, best_j = 0, 0

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            s = match if query[i - 1] == target[j - 1] else mismatch

            # Affine gap
            E[i][j] = max(
                E[i][j - 1] + gap_ext,
                H[i][j - 1] + gap_open + gap_ext,
            )
            F[i][j] = max(
                F[i - 1][j] + gap_ext,
                H[i - 1][j] + gap_open + gap_ext,
            )

            candidates = [
                0,                          # stop
                H[i - 1][j - 1] + s,       # diagonal
                F[i][j],                    # gap in target (insertion in query)
                E[i][j],                    # gap in query (deletion from query)
            ]
            best_val = max(candidates)
            H[i][j] = best_val
            TB[i][j] = np.argmax(candidates)

            if best_val > best_score:
                best_score = best_val
                best_i, best_j = i, j

    # Traceback
    aq, at = [], []
    i, j = best_i, best_j
    while i > 0 and j > 0 and H[i][j] > 0:
        tb = TB[i][j]
        if tb == 1:
            aq.append(query[i - 1])
            at.append(target[j - 1])
            i -= 1
            j -= 1
        elif tb == 2:
            aq.append(query[i - 1])
            at.append("-")
            i -= 1
        elif tb == 3:
            aq.append("-")
            at.append(target[j - 1])
            j -= 1
        else:
            break

    aq = "".join(reversed(aq))
    at = "".join(reversed(at))
    target_start = j
    target_end   = best_j

    return float(best_score), aq, at, target_start, target_end


# ---------------------------------------------------------------------------
# Fast ungapped alignment (used when max_mismatches < 3 and no gap needed)
# ---------------------------------------------------------------------------
def ungapped_align(query: str, target_window: str) -> Tuple[float, List[int], int]:
    """
    Align `query` directly against `target_window` (same length, no gaps).
    Returns (score, mismatch_positions, gap_count=0).
    Score = matches * SW_MATCH + mismatches * SW_MISMATCH.
    """
    if len(query) != len(target_window):
        # Pad/trim target to query length
        target_window = target_window[:len(query)].ljust(len(query), "N")

    mm_pos = []
    score = 0.0
    for i, (q, t) in enumerate(zip(query, target_window)):
        if q == t:
            score += SW_MATCH
        else:
            score += SW_MISMATCH
            mm_pos.append(i)
    return score, mm_pos, 0


# ---------------------------------------------------------------------------
# Main search function
# ---------------------------------------------------------------------------
def find_binding_sites(
    primer_seq: str,
    strand: str,
    seq_id: str,
    genome_seq: str,
    candidate_starts: List[int],
    primer_len: int,
    max_mismatches: int,
    three_prime_strict: bool = True,
    use_sw: bool = True,
) -> List[BindingSite]:
    """
    For each candidate start position, align the primer and return BindingSite
    objects that pass the mismatch filter.

    Parameters
    ----------
    primer_seq : str
        5'→3' primer sequence (already RC'd if searching on '-' strand).
    strand : str
        '+' or '-' — which strand the primer is designed to bind.
    seq_id : str
        Genome record identifier.
    genome_seq : str
        Full forward-strand sequence of the genome record.
    candidate_starts : list of int
        0-based start positions on the FORWARD strand where the primer is
        estimated to begin (from seed lookup).
    primer_len : int
        Length of primer in bp.
    max_mismatches : int
    three_prime_strict : bool
        If True, reject hits with a 3'-terminal mismatch.
    use_sw : bool
        Use Smith-Waterman (True) or fast ungapped alignment (False).

    Returns
    -------
    List of BindingSite objects (filtered, scored).
    """
    results: List[BindingSite] = []
    genome_len = len(genome_seq)

    # Padding around candidate for SW to handle small indels
    pad = max_mismatches + 2

    for cand_start in candidate_starts:
        # Extract template window (forward strand)
        win_start = max(0, cand_start - pad)
        win_end   = min(genome_len, cand_start + primer_len + pad)
        template_window = genome_seq[win_start:win_end]

        if not template_window:
            continue

        # If strand=='-', the primer binds the reverse-complement of the genome.
        # We've been passed the RC of the primer, so align it against the
        # reverse-complement of the template window.
        if strand == "-":
            from .genome_index import _rc
            align_target = _rc(template_window)
        else:
            align_target = template_window

        if use_sw and (max_mismatches > 0 or len(primer_seq) != len(align_target)):
            score, aq, at, t_start, t_end = smith_waterman_align(primer_seq, align_target)
            aligned_primer    = aq
            aligned_template  = at
            # Compute mismatches from alignment
            mm_pos = []
            gap_count = 0
            ppos = 0   # position in primer (0-based)
            for c_q, c_t in zip(aq, at):
                if c_q == "-":
                    gap_count += 1
                    continue
                if c_t == "-":
                    gap_count += 1
                    ppos += 1
                    continue
                if c_q != c_t:
                    mm_pos.append(ppos)
                ppos += 1

            # Actual genomic coordinates of the aligned region
            if strand == "+":
                g_start = win_start + t_start
                g_end   = win_start + t_end
            else:
                # rc alignment: t_start/t_end are in RC coords
                window_len = len(template_window)
                g_end   = win_end - t_start
                g_start = win_end - t_end
        else:
            # Fast ungapped
            # Align primer directly against same-length window
            direct_window = genome_seq[cand_start : cand_start + primer_len]
            if strand == "-":
                from .genome_index import _rc
                direct_window = _rc(direct_window)
            score, mm_pos, gap_count = ungapped_align(primer_seq, direct_window)
            aligned_template = direct_window
            g_start = cand_start
            g_end   = cand_start + primer_len

        # Coverage check: alignment must cover nearly the full primer
        aligned_query_len = len(aq) - aq.count("-") if use_sw else primer_len
        if aligned_query_len < primer_len - max_mismatches - 1:
            continue

        mismatch_count = len(mm_pos)
        if mismatch_count > max_mismatches:
            continue

        # 3'-end mismatch check
        three_prime_mm = (len(mm_pos) > 0 and mm_pos[-1] == len(primer_seq) - 1)
        if three_prime_strict and three_prime_mm:
            continue

        # Get the raw template subsequence for reporting
        raw_template_subseq = genome_seq[g_start:g_end] if g_end <= genome_len else ""

        site = BindingSite(
            seq_id           = seq_id,
            start            = g_start,
            end              = g_end,
            strand           = strand,
            primer_seq       = primer_seq,
            aligned_template = raw_template_subseq,
            sw_score         = float(score),
            mismatch_count   = mismatch_count,
            mismatch_pos     = mm_pos,
            gap_count        = gap_count,
            three_prime_mm   = three_prime_mm,
        )
        results.append(site)

    return results


# ---------------------------------------------------------------------------
# High-level search: one primer against one genome record
# ---------------------------------------------------------------------------
def search_primer_on_record(
    primer_seq: str,
    strand: str,
    seq_id: str,
    genome_seq: str,
    index: dict,
    k: int,
    max_mismatches: int,
    three_prime_strict: bool = True,
) -> List[BindingSite]:
    """
    Orchestrates seed lookup + SW extension for a single primer × genome record.

    Strand convention for k-mer index lookup:
      The index stores k-mers from the FORWARD strand under strand='+', and
      k-mers from rc(forward strand) under strand='-'.  When the rev primer
      (5'→3') is indexed on the '-' entry, its own k-mers appear verbatim in
      the '-' strand index at the corresponding fwd-strand coordinates.
      Therefore we always look up `primer_seq` directly in its own strand's
      index — no rc transformation is needed here.
    """
    from .genome_index import lookup_seeds

    # Search primer_seq in the matching strand of the index.
    raw_hits = lookup_seeds(primer_seq, index, k=k, strand_filter=strand)

    # Deduplicate candidate starts (multiple seeds → same locus)
    starts_set = set()
    for (sid, gpos, st, offset) in raw_hits:
        if sid == seq_id:
            starts_set.add(gpos)

    if not starts_set:
        return []

    return find_binding_sites(
        primer_seq       = primer_seq,
        strand           = strand,
        seq_id           = seq_id,
        genome_seq       = genome_seq,
        candidate_starts = list(starts_set),
        primer_len       = len(primer_seq),
        max_mismatches   = max_mismatches,
        three_prime_strict = three_prime_strict,
        use_sw           = True,
    )

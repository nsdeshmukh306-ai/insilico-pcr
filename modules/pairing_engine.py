"""
Pairing Engine (Layer 6)
=========================
Pairs forward and reverse primer binding sites to identify valid amplicons.

Pairing rules (biologically grounded):
  1. Forward primer must bind the FORWARD strand (strand='+').
  2. Reverse primer must bind the REVERSE strand (strand='-').
  3. They must be on the SAME chromosome/contig (same seq_id).
  4. The forward primer's start must be UPSTREAM (smaller coordinate) of
     the reverse primer's start.
  5. The amplicon size (distance between outer ends) must fall within
     [min_amplicon_size, max_amplicon_size].

Orientation diagram:
    fwd_start          rev_end
    5'--[fwd]---------- --[rev_rc]--3'  (fwd strand)
    3'------- ----------[rev]------5'  (rev strand)
         →→→→→→→→→→→→→→→
         amplicon region

The amplicon spans: fwd_start → rev_end (inclusive of both primers).
Amplicon length = rev_end - fwd_start.

Note: the reverse primer binding site on the forward strand occupies
  [rev_start, rev_end] where rev_start < rev_end.
  Since it binds the REVERSE strand, its 3' end points LEFT (towards fwd_start).
  So we need: rev_start >= fwd_end  (primers do not overlap)
  And: rev_end - fwd_start ∈ [min_size, max_size].
"""

import logging
from dataclasses import dataclass, field
from typing import List, Tuple

from .binding_eval import EvaluatedSite

log = logging.getLogger(__name__)


@dataclass
class PrimerPairHit:
    """
    A valid forward+reverse primer pair hit that can produce an amplicon.
    """
    pair_name:    str
    seq_id:       str
    fwd_site:     EvaluatedSite
    rev_site:     EvaluatedSite
    amplicon_start: int    # = fwd_site.site.start (0-based, fwd strand)
    amplicon_end:   int    # = rev_site.site.end   (0-based, fwd strand, exclusive)
    amplicon_size:  int    # = amplicon_end - amplicon_start

    @property
    def combined_binding_score(self) -> float:
        """Average binding score of both primers."""
        return (self.fwd_site.binding_score + self.rev_site.binding_score) / 2.0

    @property
    def tm_delta(self) -> float:
        """Absolute Tm difference between fwd and rev primers. Ideally <5°C."""
        return abs(self.fwd_site.tm - self.rev_site.tm)


def pair_binding_sites(
    pair_name: str,
    fwd_sites: List[EvaluatedSite],
    rev_sites: List[EvaluatedSite],
    min_amplicon_size: int = 50,
    max_amplicon_size: int = 3000,
) -> List[PrimerPairHit]:
    """
    Find all valid forward × reverse pairings.

    Parameters
    ----------
    pair_name : str
    fwd_sites : EvaluatedSite list (strand must be '+')
    rev_sites : EvaluatedSite list (strand must be '-')
    min_amplicon_size, max_amplicon_size : int

    Returns
    -------
    List of PrimerPairHit objects, sorted by amplicon_start.
    """
    hits: List[PrimerPairHit] = []

    # Group by seq_id for efficient pairing
    fwd_by_chr: dict[str, list] = {}
    for fs in fwd_sites:
        fwd_by_chr.setdefault(fs.site.seq_id, []).append(fs)

    rev_by_chr: dict[str, list] = {}
    for rs in rev_sites:
        rev_by_chr.setdefault(rs.site.seq_id, []).append(rs)

    for seq_id in set(fwd_by_chr) & set(rev_by_chr):
        fwds = sorted(fwd_by_chr[seq_id], key=lambda x: x.site.start)
        revs = sorted(rev_by_chr[seq_id], key=lambda x: x.site.start)

        for fs in fwds:
            f_start = fs.site.start
            f_end   = fs.site.end

            for rs in revs:
                r_start = rs.site.start
                r_end   = rs.site.end

                # Primer orientation check:
                # fwd binds '+', rev binds '-'.
                # The rev primer's 3' end is at r_start (leftmost on fwd strand).
                # We need the amplicon to be: fwd_start → rev_end.
                if rs.site.strand != "-":
                    continue
                if fs.site.strand != "+":
                    continue

                # Fwd primer must be upstream of rev primer
                if f_start >= r_start:
                    continue

                # Primers should not overlap
                if f_end > r_start:
                    continue

                amplicon_size = r_end - f_start
                if amplicon_size < min_amplicon_size or amplicon_size > max_amplicon_size:
                    continue

                hit = PrimerPairHit(
                    pair_name       = pair_name,
                    seq_id          = seq_id,
                    fwd_site        = fs,
                    rev_site        = rs,
                    amplicon_start  = f_start,
                    amplicon_end    = r_end,
                    amplicon_size   = amplicon_size,
                )
                hits.append(hit)

    hits.sort(key=lambda h: (h.seq_id, h.amplicon_start))
    log.debug("Pair '%s': found %d valid amplicon(s).", pair_name, len(hits))
    return hits

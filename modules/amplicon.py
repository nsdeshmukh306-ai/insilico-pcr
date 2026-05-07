"""
Amplicon Extraction (Layer 7)
==============================
Extracts the full amplicon sequence from the genome for each PrimerPairHit
and computes amplicon-level properties.

The amplicon sequence INCLUDES the primer sequences at both ends, reflecting
the standard PCR product definition (primers are incorporated into the product).
"""

from dataclasses import dataclass
from typing import Dict, List

from .pairing_engine import PrimerPairHit
from .thermodynamics import gc_content


@dataclass
class Amplicon:
    """
    A fully characterised PCR amplicon.
    """
    pair_name:    str
    seq_id:       str
    start:        int        # 0-based, fwd strand (inclusive of fwd primer)
    end:          int        # 0-based, fwd strand (exclusive; includes rev primer)
    length:       int
    sequence:     str        # Amplicon sequence on the forward strand
    gc_fraction:  float
    fwd_primer:   str        # Forward primer sequence
    rev_primer:   str        # Reverse primer sequence (as supplied, 5'→3')
    fwd_tm:       float      # Tm of fwd primer (°C)
    rev_tm:       float      # Tm of rev primer (°C)
    fwd_mm:       int        # Mismatches in fwd binding
    rev_mm:       int        # Mismatches in rev binding
    fwd_binding_score: float
    rev_binding_score: float
    hit:          PrimerPairHit   # Reference to the pairing hit


def extract_amplicon(
    hit: PrimerPairHit,
    genome_map: Dict[str, str],
) -> Amplicon:
    """
    Extract the amplicon sequence for a PrimerPairHit.

    Parameters
    ----------
    hit : PrimerPairHit
    genome_map : dict  {seq_id: forward_strand_sequence}

    Returns
    -------
    Amplicon
    """
    genome_seq = genome_map.get(hit.seq_id, "")
    seq = genome_seq[hit.amplicon_start : hit.amplicon_end]

    return Amplicon(
        pair_name         = hit.pair_name,
        seq_id            = hit.seq_id,
        start             = hit.amplicon_start,
        end               = hit.amplicon_end,
        length            = hit.amplicon_size,
        sequence          = seq,
        gc_fraction       = round(gc_content(seq), 3),
        fwd_primer        = hit.fwd_site.site.primer_seq,
        rev_primer        = hit.rev_site.site.primer_seq,
        fwd_tm            = hit.fwd_site.tm,
        rev_tm            = hit.rev_site.tm,
        fwd_mm            = hit.fwd_site.site.mismatch_count,
        rev_mm            = hit.rev_site.site.mismatch_count,
        fwd_binding_score = hit.fwd_site.binding_score,
        rev_binding_score = hit.rev_site.binding_score,
        hit               = hit,
    )


def extract_amplicons(
    hits: List[PrimerPairHit],
    genome_map: Dict[str, str],
) -> List[Amplicon]:
    """Extract amplicons for a list of PrimerPairHit objects."""
    return [extract_amplicon(h, genome_map) for h in hits]


def amplicon_length_score(length: int) -> float:
    """
    Score amplicon length suitability (0–1).
    Peak at 200 bp; graceful decline for longer products.
    Standard diagnostic PCR: 100–1000 bp.
    """
    if length <= 0:
        return 0.0
    # Gaussian-like peak at 200 bp, σ=400 bp (wide, to allow 50–3000 bp)
    import math
    mu    = 200.0
    sigma = 500.0
    score = math.exp(-0.5 * ((length - mu) / sigma) ** 2)
    # Penalise very short (<50 bp) and very long (>2000 bp) amplicons extra
    if length < 50:
        score *= 0.1
    elif length > 2000:
        score *= max(0.1, 1.0 - (length - 2000) / 1000.0)
    return round(score, 4)

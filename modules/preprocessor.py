"""
Preprocessor (Layer 2)
======================
Cleans and normalises genome sequences and primer sequences before indexing
and alignment. Also expands IUPAC ambiguity codes in primers.

Operations:
  1. Sequence sanitisation (remove non-DNA characters, log Ns)
  2. Primer Tm pre-computation (stored on PrimerPair objects)
  3. IUPAC expansion → list of unambiguous variant sequences per primer
  4. Reverse-complement generation for each primer variant

Limitations:
  - IUPAC expansion is combinatorial; very degenerate primers (many N/R/Y bases)
    can produce hundreds of variants. An expansion cap is enforced (default 64).
"""

import itertools
import logging
from typing import List, Tuple

from Bio.SeqRecord import SeqRecord

from .input_handler import IUPAC_BASES, PCRParams, PrimerPair
from .thermodynamics import calc_tm, reverse_complement

log = logging.getLogger(__name__)

MAX_EXPANSION = 64   # Maximum IUPAC variant sequences per primer position


# ---------------------------------------------------------------------------
# Sequence sanitisation
# ---------------------------------------------------------------------------
def sanitise_sequence(seq: str) -> Tuple[str, int]:
    """
    Return (cleaned_seq, n_count).
    Replaces non-ACGTN characters with 'N' and logs a warning.
    """
    cleaned = []
    replaced = 0
    for base in seq.upper():
        if base in ("A", "C", "G", "T", "N"):
            cleaned.append(base)
        else:
            cleaned.append("N")
            replaced += 1
    if replaced:
        log.warning("Replaced %d non-standard characters with N.", replaced)
    n_count = cleaned.count("N")
    return "".join(cleaned), n_count


def preprocess_genome_record(record: SeqRecord) -> SeqRecord:
    """
    Sanitise a genome SeqRecord in place. Returns the modified record.
    Logs statistics about ambiguous bases.
    """
    raw = str(record.seq)
    clean, n_count = sanitise_sequence(raw)
    frac_n = n_count / len(clean) if clean else 0.0
    if frac_n > 0.05:
        log.warning(
            "Record %s has %.1f%% ambiguous (N) bases — "
            "primer binding in these regions will be unreliable.",
            record.id, frac_n * 100,
        )
    record.seq = record.seq.__class__(clean)
    return record


# ---------------------------------------------------------------------------
# IUPAC expansion
# ---------------------------------------------------------------------------
def expand_iupac(seq: str, cap: int = MAX_EXPANSION) -> List[str]:
    """
    Expand a potentially IUPAC-coded primer into a list of unambiguous
    DNA sequences (all combinations of ambiguous bases).

    Parameters
    ----------
    seq : str
        Uppercase primer sequence, may contain IUPAC codes.
    cap : int
        Maximum number of expansions returned. If the true count exceeds cap,
        a random sample of `cap` combinations is returned with a warning.

    Returns
    -------
    List of strings (each an unambiguous ACGT sequence).
    """
    # Build list of possibilities for each position
    choices = [IUPAC_BASES.get(b, [b]) for b in seq.upper()]

    # Quick count
    total = 1
    for c in choices:
        total *= len(c)

    if total > cap:
        log.warning(
            "IUPAC expansion of '%s' yields %d variants; capping at %d.",
            seq, total, cap,
        )
        # Take the first `cap` combinations (deterministic)
        variants = []
        for combo in itertools.islice(itertools.product(*choices), cap):
            variants.append("".join(combo))
        return variants

    return ["".join(combo) for combo in itertools.product(*choices)]


# ---------------------------------------------------------------------------
# Primer preprocessing
# ---------------------------------------------------------------------------
def preprocess_primer_pair(
    pair: PrimerPair,
    params: PCRParams,
) -> Tuple[List[str], List[str]]:
    """
    Expand IUPAC codes and compute Tm on a PrimerPair.

    Returns
    -------
    (fwd_variants, rev_variants) : lists of unambiguous primer sequences
        fwd_variants : forward primer variants (5'→3')
        rev_variants : reverse primer variants (5'→3' of REVERSE primer as given)
    """
    fwd_variants = expand_iupac(pair.forward)
    rev_variants = expand_iupac(pair.reverse)

    # Compute Tm on the primary (first) variant and store on pair
    pair.fwd_tm = calc_tm(
        fwd_variants[0],
        na_conc=params.na_conc,
        primer_conc=params.primer_conc,
    )
    pair.rev_tm = calc_tm(
        rev_variants[0],
        na_conc=params.na_conc,
        primer_conc=params.primer_conc,
    )

    log.debug(
        "Pair '%s': fwd_Tm=%.1f°C  rev_Tm=%.1f°C  fwd_variants=%d  rev_variants=%d",
        pair.name, pair.fwd_tm, pair.rev_tm,
        len(fwd_variants), len(rev_variants),
    )
    return fwd_variants, rev_variants


def get_search_sequences(primer_seq: str) -> Tuple[str, str]:
    """
    For a given primer sequence, return:
      (as_given_5to3, reverse_complement)

    The 'as_given' form is used to search the forward strand.
    The reverse-complement form is used to search the reverse strand.
    """
    rc = reverse_complement(primer_seq)
    return primer_seq, rc

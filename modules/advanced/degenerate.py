"""
Degenerate Primer Support (IUPAC)
===================================
Utilities for handling degenerate (IUPAC) primers beyond simple expansion.

This module provides:
  1. Degeneracy statistics (number of variants, effective primer concentration)
  2. Consensus sequence generation from a set of primers
  3. Degenerate primer design helper: given a set of target sequences,
     compute the most parsimonious degenerate primer covering all variants

IUPAC codes:
  R=AG, Y=CT, S=GC, W=AT, K=GT, M=AC
  B=CGT, D=AGT, H=ACT, V=ACG, N=ACGT
"""

import itertools
from typing import List, Set

IUPAC_MAP = {
    frozenset("A"):    "A",
    frozenset("C"):    "C",
    frozenset("G"):    "G",
    frozenset("T"):    "T",
    frozenset("AG"):   "R",
    frozenset("CT"):   "Y",
    frozenset("GC"):   "S",
    frozenset("AT"):   "W",
    frozenset("GT"):   "K",
    frozenset("AC"):   "M",
    frozenset("CGT"):  "B",
    frozenset("AGT"):  "D",
    frozenset("ACT"):  "H",
    frozenset("ACG"):  "V",
    frozenset("ACGT"): "N",
}

IUPAC_EXPAND = {v: k for k, v in IUPAC_MAP.items()}


def degeneracy(seq: str) -> int:
    """Return total number of unambiguous sequences represented by a degenerate primer."""
    total = 1
    for b in seq.upper():
        bases = IUPAC_EXPAND.get(b, frozenset(b))
        total *= len(bases)
    return total


def effective_conc(primer_conc: float, degen: int) -> float:
    """
    Effective concentration of each variant when degeneracy is degen.
    Assumes equal representation of all variants.
    """
    return primer_conc / max(1, degen)


def consensus_to_iupac(sequences: List[str]) -> str:
    """
    Given a list of aligned sequences (same length), return the IUPAC
    degenerate consensus sequence that covers all variants at each position.
    Sequences must be the same length.
    """
    if not sequences:
        return ""
    length = len(sequences[0])
    result = []
    for i in range(length):
        bases = frozenset(s[i].upper() for s in sequences if i < len(s))
        # Remove Ns from the set before encoding
        clean_bases = bases - {"N"}
        if not clean_bases:
            clean_bases = frozenset("N")
        iupac_code = IUPAC_MAP.get(clean_bases, "N")
        result.append(iupac_code)
    return "".join(result)


def design_degenerate_primer(
    target_sequences: List[str],
    max_degeneracy: int = 64,
) -> dict:
    """
    Design the most parsimonious degenerate primer for a set of target sequences.
    All sequences must be the same length (pre-aligned).

    Parameters
    ----------
    target_sequences : list of str (same length, pre-aligned)
    max_degeneracy : int  cap on degeneracy (default 64)

    Returns
    -------
    dict with keys:
      consensus : str  IUPAC degenerate primer
      degeneracy: int
      variants  : list of str  (if degeneracy <= max_degeneracy)
      warning   : str
    """
    if not target_sequences:
        return {"consensus": "", "degeneracy": 0, "variants": [], "warning": "No sequences provided."}

    lengths = set(len(s) for s in target_sequences)
    if len(lengths) > 1:
        return {
            "consensus": "",
            "degeneracy": 0,
            "variants": [],
            "warning": f"Sequences have different lengths: {lengths}. Please align first.",
        }

    consensus = consensus_to_iupac(target_sequences)
    degen = degeneracy(consensus)

    warning = ""
    if degen > max_degeneracy:
        warning = (
            f"Degeneracy {degen} exceeds maximum {max_degeneracy}. "
            "Consider splitting into multiple specific primers."
        )

    from ..preprocessor import expand_iupac
    variants = expand_iupac(consensus, cap=max_degeneracy) if degen <= max_degeneracy else []

    return {
        "consensus":  consensus,
        "degeneracy": degen,
        "variants":   variants,
        "warning":    warning,
    }

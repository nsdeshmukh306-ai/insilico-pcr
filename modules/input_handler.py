"""
Input Handler (Layer 1)
=======================
Parses and validates all inputs to the in-silico PCR pipeline:
  - Primer sequences (forward, reverse)
  - Genome / transcript FASTA files
  - Run parameters (mismatch tolerance, amplicon size limits, etc.)

Supported input formats:
  - Primer sequences: raw string, FASTA, JSON
  - Genome: FASTA (gzipped or plain)

No external alignment tools assumed. Biopython SeqIO is used for FASTA parsing.
"""

import gzip
import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator, List, Optional, Union

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


# ---------------------------------------------------------------------------
# IUPAC ambiguity codes
# ---------------------------------------------------------------------------
IUPAC_BASES = {
    "A": ["A"],
    "C": ["C"],
    "G": ["G"],
    "T": ["T"],
    "R": ["A", "G"],
    "Y": ["C", "T"],
    "S": ["G", "C"],
    "W": ["A", "T"],
    "K": ["G", "T"],
    "M": ["A", "C"],
    "B": ["C", "G", "T"],
    "D": ["A", "G", "T"],
    "H": ["A", "C", "T"],
    "V": ["A", "C", "G"],
    "N": ["A", "C", "G", "T"],
}

_VALID_IUPAC = set(IUPAC_BASES.keys())
_VALID_STRICT = {"A", "C", "G", "T"}


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------
@dataclass
class PrimerPair:
    """A forward/reverse primer pair with metadata."""
    name:       str
    forward:    str          # 5'→3' sequence, uppercase
    reverse:    str          # 5'→3' sequence of the reverse primer, uppercase
    # Thermodynamics populated later
    fwd_tm:     float = 0.0
    rev_tm:     float = 0.0

    def __post_init__(self):
        self.forward = self.forward.upper().replace("U", "T")
        self.reverse = self.reverse.upper().replace("U", "T")


@dataclass
class PCRParams:
    """Run-time parameters controlling PCR simulation behaviour."""
    max_mismatches:     int   = 3      # Max allowed mismatches per primer
    min_amplicon_size:  int   = 50     # bp
    max_amplicon_size:  int   = 3000   # bp
    primer_conc:        float = 250e-9 # 250 nM
    na_conc:            float = 0.05   # 50 mM Na⁺  (mol/L)
    mg_conc:            float = 0.0    # Mg²⁺ total  (mol/L); 0 → Na-only correction
    dntp_conc:          float = 0.0    # dNTP total  (mol/L); chelates Mg²⁺ 1:1
    seed_length:        int   = 8      # k-mer seed length for indexing
    three_prime_strict: bool  = True   # Reject 3'-terminal mismatches by default
    max_off_target:     int   = 50     # Limit off-target hits reported
    n_threads:          int   = 1      # Parallel workers


# ---------------------------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------------------------
def _clean_seq(seq: str) -> str:
    """Strip whitespace, convert to uppercase."""
    return re.sub(r"\s+", "", seq).upper().replace("U", "T")


def validate_primer(seq: str, allow_iupac: bool = True) -> str:
    """
    Validate a primer sequence. Returns cleaned sequence or raises ValueError.
    """
    seq = _clean_seq(seq)
    valid = _VALID_IUPAC if allow_iupac else _VALID_STRICT
    bad = set(seq) - valid
    if bad:
        raise ValueError(f"Invalid bases in primer: {bad!r}  sequence: {seq}")
    if len(seq) < 10:
        raise ValueError(f"Primer too short ({len(seq)} bp): {seq}")
    if len(seq) > 60:
        raise ValueError(f"Primer unusually long ({len(seq)} bp); check input.")
    return seq


# ---------------------------------------------------------------------------
# Primer input parsers
# ---------------------------------------------------------------------------
def parse_primers_from_strings(
    fwd: str,
    rev: str,
    name: str = "primer_pair_1",
    allow_iupac: bool = True,
) -> PrimerPair:
    """Parse a single primer pair supplied as raw strings."""
    return PrimerPair(
        name    = name,
        forward = validate_primer(fwd, allow_iupac),
        reverse = validate_primer(rev, allow_iupac),
    )


def parse_primers_from_fasta(path: Union[str, Path]) -> List[PrimerPair]:
    """
    Parse primer pairs from a FASTA file.
    Expects alternating records: forward then reverse.
    Record IDs ending in '_F', '_fwd', '_forward' are treated as forward primers;
    '_R', '_rev', '_reverse' as reverse primers.
    If no suffix, records are treated as alternating F/R.
    """
    path = Path(path)
    records = list(SeqIO.parse(str(path), "fasta"))
    if len(records) % 2 != 0:
        raise ValueError(f"FASTA file {path} has odd number of sequences; need F/R pairs.")
    pairs = []
    for i in range(0, len(records), 2):
        fwd_rec = records[i]
        rev_rec = records[i + 1]
        pair = PrimerPair(
            name    = fwd_rec.id.rstrip("_FfwdForward").rstrip("_"),
            forward = validate_primer(str(fwd_rec.seq)),
            reverse = validate_primer(str(rev_rec.seq)),
        )
        pairs.append(pair)
    return pairs


def parse_primers_from_json(path: Union[str, Path]) -> List[PrimerPair]:
    """
    Parse primer pairs from a JSON file.

    Expected format:
    [
      {"name": "pair1", "forward": "ATCG...", "reverse": "CGTA..."},
      ...
    ]
    """
    path = Path(path)
    with open(path) as fh:
        data = json.load(fh)
    if isinstance(data, dict):
        data = [data]
    pairs = []
    for i, entry in enumerate(data):
        pair = PrimerPair(
            name    = entry.get("name", f"pair_{i+1}"),
            forward = validate_primer(entry["forward"]),
            reverse = validate_primer(entry["reverse"]),
        )
        pairs.append(pair)
    return pairs


# ---------------------------------------------------------------------------
# Genome / FASTA input
# ---------------------------------------------------------------------------
def iter_genome_records(
    fasta_path: Union[str, Path],
    max_records: Optional[int] = None,
) -> Iterator[SeqRecord]:
    """
    Lazily yield SeqRecord objects from a FASTA file (plain or gzipped).
    Sequences are converted to uppercase strings.

    Parameters
    ----------
    fasta_path : path-like
    max_records : int, optional
        Stop after this many records (useful for testing).
    """
    fasta_path = Path(fasta_path)
    opener = gzip.open if fasta_path.suffix in (".gz", ".gzip") else open
    with opener(str(fasta_path), "rt") as fh:
        for idx, rec in enumerate(SeqIO.parse(fh, "fasta")):
            if max_records is not None and idx >= max_records:
                break
            rec.seq = rec.seq.__class__(str(rec.seq).upper().replace("U", "T"))
            yield rec


def load_genome_records(
    fasta_path: Union[str, Path],
    max_records: Optional[int] = None,
) -> List[SeqRecord]:
    """Load all genome records into memory. Use iter_genome_records for large genomes."""
    return list(iter_genome_records(fasta_path, max_records))


# ---------------------------------------------------------------------------
# Parameter loading
# ---------------------------------------------------------------------------
def parse_params_from_json(path: Union[str, Path]) -> PCRParams:
    """Load PCRParams from a JSON config file."""
    path = Path(path)
    with open(path) as fh:
        d = json.load(fh)
    return PCRParams(**{k: v for k, v in d.items() if k in PCRParams.__dataclass_fields__})


def default_params() -> PCRParams:
    return PCRParams()

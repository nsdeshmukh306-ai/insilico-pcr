#!/usr/bin/env python3
"""
Generate a synthetic example genome FASTA for testing.

Contains:
  - 3 chromosomes (chr1: 5000 bp, chr2: 3000 bp, chr3: 2000 bp)
  - Known primer sites embedded at predictable positions:
      ACTB_exon5:
        FWD: GCACTGGTGGCATCGATCTA  at chr1:500–520 (fwd strand)
        REV: GCTAGCTAGCATGCTAGCTA  at chr1:900–920 (rev strand; RC = TAGCTAGCATGCTAGCTAGC)
        → amplicon chr1:500–920 = 420 bp
      GAPDH_qPCR:
        FWD: GTCTCCTCTGACTTCAACAGCG at chr2:200–222 (fwd strand)
        REV: TTGGCTACAGCAACAGGGTGG  at chr2:500–521 (rev strand; RC = ACCACCCTGTTGCTGTAGCCAA)
        → amplicon chr2:200–521 = 321 bp
"""

import random

random.seed(42)

def random_dna(n: int) -> str:
    return "".join(random.choice("ACGT") for _ in range(n))

def embed(genome: list, pos: int, seq: str) -> None:
    """Embed seq into genome (list of chars) at pos."""
    for i, b in enumerate(seq):
        if pos + i < len(genome):
            genome[pos + i] = b

def rc(seq: str) -> str:
    comp = {"A":"T","T":"A","G":"C","C":"G"}
    return "".join(comp[b] for b in reversed(seq.upper()))

# chr1: 5000 bp
chr1 = list(random_dna(5000))
fwd1 = "GCACTGGTGGCATCGATCTA"
rev1 = "TAGCTAGCATGCTAGCTAGC"
embed(chr1, 500, fwd1)
embed(chr1, 900, rc(rev1))   # RC of rev primer on fwd strand
# Off-target site with 2 mismatches (chr1:2000)
fwd1_mm = list(fwd1)
fwd1_mm[5]  = "G"  # mismatch
fwd1_mm[12] = "T"  # mismatch
embed(chr1, 2000, "".join(fwd1_mm))
embed(chr1, 2350, rc(rev1))

# chr2: 3000 bp
chr2 = list(random_dna(3000))
fwd2 = "GTCTCCTCTGACTTCAACAGCG"
rev2 = "ACCACCCTGTTGCTGTAGCCAA"
embed(chr2, 200, fwd2)
embed(chr2, 500, rc(rev2))

# chr3: 2000 bp (no primers here — tests specificity)
chr3 = list(random_dna(2000))

records = [
    ("chr1", "".join(chr1)),
    ("chr2", "".join(chr2)),
    ("chr3", "".join(chr3)),
]

with open("example_genome.fa", "w") as fh:
    for name, seq in records:
        fh.write(f">{name}\n")
        # Write in 80-char lines
        for i in range(0, len(seq), 80):
            fh.write(seq[i:i+80] + "\n")

print("Generated example_genome.fa")
print("  chr1: 5000 bp  (ACTB_exon5 amplicon at 500-920, 420 bp)")
print("  chr2: 3000 bp  (GAPDH_qPCR amplicon at 200-521, 321 bp)")
print("  chr3: 2000 bp  (no intended amplicons)")

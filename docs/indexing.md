# Genome Indexing

Two index backends are used, selected automatically based on genome size.

---

## Selection Logic

```python
THRESHOLD = 50_000_000  # 50 Mbp

if len(genome) >= THRESHOLD:
    index = build_fm_index(genome)     # BWT-based; O(n) space
else:
    index = build_kmer_index(genome)   # Hash-based; fast for small genomes
```

---

## K-mer Index (`genome_index.py`)

Used for genomes < 50 Mbp.

### How it works

A Python `dict` maps every k-mer (default k=8) to a list of positions where it occurs in the genome.

```
Genome:   ATCGATCGATCG...
k=4       ATCG → [0, 4, 8, ...]
          TCGA → [1, 5, 9, ...]
          CGAT → [2, 6, ...]
          ...
```

### Search

For a primer query, the k-mer at position 0 (seed) is looked up, returning candidate positions. Smith–Waterman alignment is then run at each candidate to score the full primer–template match.

### Complexity

| Operation | Time | Space |
|---|---|---|
| Build | O(n) | O(n·k) |
| Seed lookup | O(1) | — |
| Full alignment | O(m²) per candidate | — |

---

## FM-Index (`genome_index_fm.py`)

Used for genomes ≥ 50 Mbp (e.g. chromosome-scale).

### What is an FM-index?

The FM-index is based on the **Burrows–Wheeler Transform (BWT)** — a reversible string transformation that groups similar substrings together, enabling pattern search in O(m) time (m = query length) using only O(n) space for the index.

### Burrows–Wheeler Transform

The BWT rotates all cyclic shifts of the genome and sorts them. The last column of this sorted matrix is the BWT. It has the property that identical characters cluster together in runs, enabling compression and efficient backward search.

```
Genome BANANA$ →  sorted rotations →  BWT: ANNB$AA
```

### Backward search

Query `PRIMER` is searched right-to-left:
1. Start with the full suffix array range
2. For each character (rightmost first), narrow the range using the FM-index
3. Result: all positions in O(m log n) time

### In this implementation

```python
class FMIndex:
    def __init__(self, text: str): ...
    def search(self, pattern: str) -> List[int]: ...

class GenomeFMIndex:
    def __init__(self, genome: dict[str, str]): ...   # per-chromosome
    def lookup_seeds(self, primer: str, k: int) -> List[tuple]: ...
```

The `GenomeFMIndex` wraps per-chromosome `FMIndex` objects and handles multi-chromosome genomes transparently.

### Complexity

| Operation | Time | Space |
|---|---|---|
| Build | O(n log n) | O(n) |
| Pattern search | O(m log n) | — |

---

## Seed-and-Extend Strategy

Both indices use a **seed-and-extend** approach:

1. **Seed:** look up a short k-mer from the primer's 5′ end
2. **Filter:** discard candidates more than `max_mismatches` away
3. **Extend:** run Smith–Waterman alignment over a window `[pos - flank, pos + primer_len + flank]`
4. **Evaluate:** keep hits with SW score above threshold

The flank size accounts for indels. Default: `flank = max_mismatches + 2`.

---

## References

- Ferragina, P. & Manzini, G. (2000). Opportunistic data structures with applications. *FOCS 2000*. The original FM-index paper.
- Burrows, M. & Wheeler, D. (1994). A block-sorting lossless data compression algorithm. *Technical Report 124*, Digital Equipment Corporation.
- Smith, T. F. & Waterman, M. S. (1981). Identification of common molecular subsequences. *Journal of Molecular Biology*, 147(1), 195–197.

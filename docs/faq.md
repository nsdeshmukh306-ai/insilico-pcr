# FAQ

Common questions about the platform.

---

## General

**Q: What is in-silico PCR?**

PCR amplifies a specific DNA segment using two short primer sequences. In-silico PCR computationally simulates this process — finding where primers bind on a genome, how strongly, what size product results, and whether they might amplify unintended regions.

**Q: Who is this tool for?**

- Molecular biologists designing primers before running wet-lab experiments
- Computational biologists validating primer sets at scale
- Bioinformatics students learning thermodynamic models
- Anyone who wants to understand *why* a primer pair works or fails

**Q: How does it compare to Primer-BLAST or UCSC In-Silico PCR?**

| Feature | This tool | Primer-BLAST | UCSC In-Silico PCR |
|---|---|---|---|
| NN thermodynamics | Full SantaLucia 1998 | Yes | Simplified |
| Mismatch NN tables | Peyret 1999 full | Partial | No |
| Mg²⁺ correction | Owczarzy 2008 | Yes | No |
| Interactive dashboard | Yes (local) | Web only | Web only |
| Open source | Yes | No | No |
| Custom genome support | Yes (FASTA/string) | Limited | Yes |

---

## Thermodynamics

**Q: Why is the computed Tm different from another tool?**

Tm depends on: NN parameters used, salt correction formula, primer concentration, and whether mismatches are accounted for. This tool uses the SantaLucia 1998 unified NN set with Owczarzy 2008 Mg²⁺ correction. Other tools may use older or different parameter sets.

**Q: What does ΔG₃₇ mean?**

ΔG is the Gibbs free energy of duplex formation at 37°C. More negative = more stable. For PCR primers, typical values are −15 to −35 kcal/mol. Very stable duplexes (ΔG < −40) may form unintended hairpins or dimers.

**Q: Why does a G·T mismatch cause less Tm drop than an A·C mismatch?**

G·T forms a "wobble" base pair — a non-Watson-Crick pairing that still provides some stabilisation. A·C mismatches have no stabilising geometry and are more destabilising. This is modelled using the Allawi & SantaLucia 1997 (G·T) and Peyret 1999 (A·C) NN tables.

**Q: What is the "parallel complement convention" in the template?**

`template[i]` stores the 3′→5′ antiparallel base at position *i*. For a perfect match, `template[i] = complement(primer[i])`. This means: if `primer[i] = 'A'`, a perfect template has `template[i] = 'T'`. A mismatch at position *i* means `template[i] ≠ complement(primer[i])`. This is the convention expected by the NN lookup function. See `docs/thermodynamics.md` for details.

---

## Usage

**Q: Can I use my own genome?**

Yes. Provide a FASTA file via `--genome` on the CLI, `genome_fasta=` in the API, or paste/upload in the dashboard. Multi-chromosome FASTAs are supported.

**Q: What primer format is accepted?**

JSON: `[{"name": "pair_1", "forward": "ACGT...", "reverse": "ACGT..."}]`

FASTA: two sequences, first forward then reverse.

Single pair on CLI: `--fwd ACGT... --rev ACGT...`

**Q: What does "max mismatches" control?**

It sets the maximum number of primer–template mismatches allowed when reporting a binding site. Higher values find more (potentially off-target) hits but are slower. For typical primer validation, 3 is a good default. Set to 0 for exact-match only.

**Q: Why does the live re-run panel say "Genome sequence unavailable"?**

If you loaded results via the demo button, the genome string is fetched automatically. If you pasted a genome in the run panel, it is stored for live re-runs. If the error persists, go to the Run panel, paste the genome, and click "Run analysis" once.

---

## Scores

**Q: What does a score of 35 mean?**

Scores reflect the composite quality of a primer pair at a given binding site. A score of 35 is low — it typically indicates poor Tm compatibility, off-target amplification, or mismatch-driven instability. Click the amplicon row to see the score component breakdown.

**Q: Why are the scores low for the demo dataset?**

The demo uses a synthetic 10 kbp genome designed for testing, not realistic primer design. The example primers have poor Tm for this genome (the SantaLucia model gives unrealistic values for primers tested outside their design temperature range). In real use with appropriate primers, scores of 60–90 are common.

---

## Performance

**Q: How fast is it on a real genome?**

The FM-index build for a 50 Mbp chromosome takes ~2 seconds. Per-pair search takes ~500 ms. For human-genome scale (~3 Gbp), full FM-index build takes ~2 minutes; this is built once per genome and can be cached. Full human-genome validation is on the roadmap for v2.1.

**Q: Can I run multiple primer pairs?**

Yes — the JSON primer input format accepts any number of pairs, and each is analysed independently. The dashboard shows all pairs as selectable tabs.

---

## Development

**Q: How do I add a new NN mismatch parameter?**

1. Add the `"XY/WZ"` key-value pair to `_MM_NN_PARAMS` in `modules/thermodynamics.py`
2. Cite the paper and table number in a comment
3. Add a test in `tests/test_mismatch_thermodynamics.py` validating the lookup value

**Q: How do I run the tests?**

```bash
pytest tests/ -v
```

All 170 tests should pass. Thermodynamic tests include values from the original literature.

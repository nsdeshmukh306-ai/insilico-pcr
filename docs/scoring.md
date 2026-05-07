# Scoring

Each amplicon receives a composite score from 0 to 100. Higher scores indicate primer pairs more likely to work reliably in a real PCR experiment.

---

## Score Formula

```
Score = (w_bind·S_bind + w_tm·S_tm + w_gc·S_gc + w_len·S_len
         − w_mm·P_mm − w_offt·P_offt)
        ──────────────────────────────────────────────────────── × 100
                      (w_bind + w_tm + w_gc + w_len)

Score is clipped to [0, 100].
```

---

## Component Definitions

### S_bind — Binding Score (weight 0.30)

```
S_bind = (fwd_binding_score + rev_binding_score) / 200
```

`binding_score` is the normalised Smith–Waterman score (0–100) from the alignment module. A perfect match gives `binding_score = 100`.

### S_tm — Tm Compatibility (weight 0.25)

```
avg_tm_factor = mean of tm_range_factor(fwd_Tm), tm_range_factor(rev_Tm)
delta_tm_penalty = min(1, |fwd_Tm − rev_Tm| / 20°C)
S_tm = avg_tm_factor × (1 − delta_tm_penalty)
```

`tm_range_factor`:
- 1.0 if Tm ∈ [55, 68] °C (optimal range)
- Linear decline: 0 below 45 °C, 0 above 78 °C

Two primers with Tm difference > 10 °C are heavily penalised.

### S_gc — GC Content (weight 0.10)

```
S_gc = 1.0    if GC ∈ [0.40, 0.60]
S_gc = linear decline to 0 below 0.20 or above 0.80
```

### S_len — Amplicon Length (weight 0.05)

```
S_len = amplicon_length_score(length)
```

Peaks for amplicons in the range 100–1000 bp (typical PCR product). Very short (< 50 bp) or very long (> 3000 bp) amplicons score lower.

### P_mm — Mismatch Penalty (weight 0.20)

This is the most thermodynamically rigorous component.

**Previous (removed) formula:**
```
P_mm = (fwd_mm + rev_mm) / (2 × max_mismatches + 1)
```
This was purely count-based — a G·T wobble counted the same as an A·C mismatch.

**Current formula (v1.2.0+):**
```
ΔTm_fwd = max(0, Tm_perfect_fwd − Tm_observed_fwd)
ΔTm_rev = max(0, Tm_perfect_rev − Tm_observed_rev)
P_mm = min(1, (ΔTm_fwd + ΔTm_rev) / 20°C)
```

`Tm_observed` is calculated using Peyret mismatch NN tables via `calc_tm(primer, template=aligned_template)`. The aligned template uses the 3′→5′ parallel complement convention so the NN table keys resolve correctly.

**Why ΔTm is better:**
- A single G·T wobble at an internal position typically drops Tm by ~2–3 °C → P_mm ≈ 0.1–0.15
- An A·C mismatch at the 3′ terminal position drops Tm by ~8–12 °C → P_mm ≈ 0.4–0.6
- The old formula gave both cases the same P_mm = 0.33 (for max_mismatches=3)

### P_offt — Off-target Penalty (weight 0.10)

```
P_offt = log(1 + off_target_count) / log(1 + 100)
```

Log scale: 1 off-target gives P_offt ≈ 0.15; 100 off-targets gives P_offt = 1.0.

---

## Default Weights

| Component | Symbol | Default Weight | Direction |
|---|---|---|---|
| Binding score | S_bind | 0.30 | Positive |
| Tm compatibility | S_tm | 0.25 | Positive |
| GC content | S_gc | 0.10 | Positive |
| Amplicon length | S_len | 0.05 | Positive |
| Mismatch penalty | P_mm | 0.20 | Negative |
| Off-target penalty | P_offt | 0.10 | Negative |

Weights can be overridden by passing a `weights` dict to `score_amplicon()` or `rank_amplicons()`.

---

## Score Interpretation

| Score | Interpretation |
|---|---|
| 80–100 | Excellent — high confidence, good specificity |
| 60–79 | Good — likely to work, minor concerns |
| 40–59 | Fair — review Tm match and off-targets |
| 20–39 | Poor — significant issues |
| 0–19 | Not recommended |

---

## Ranking

`rank_amplicons(amplicons, intended_region=None)` scores all amplicons for a primer pair, identifies the "intended" amplicon (by genomic overlap if a region is provided, otherwise by best binding score), and returns the list sorted by `final_score` descending.

The intended amplicon has `off_target_count = 0`; all others count as off-targets for P_offt computation.

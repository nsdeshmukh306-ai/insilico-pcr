# Benchmark Figures

Figures are generated from `benchmarks/runtime_comparison.csv` and `benchmarks/benchmark_results.json`.

## Regenerating figures

```bash
python benchmarks/figures/generate_figures.py
```

Requires `matplotlib` and `pandas` (both included in `requirements.txt`).
Outputs PNG (150 dpi, for README) and SVG (for publication) for each figure.

## Figures

### `runtime_scaling.png` — log-log runtime vs genome size

![Runtime scaling](runtime_scaling.png)

Three curves: full pair analysis, index build only, per-primer search.
Vertical dashed line marks the k-mer → FM-index crossover at ~50 Mbp.

### `tm_accuracy.png` — predicted vs literature Tm

![Tm accuracy](tm_accuracy.png)

Left: per-sequence absolute error against SantaLucia 1998 Table 2 (all < 0.30°C, well under the ±0.5°C threshold).
Right: predicted vs literature scatter with ±0.5°C band.

### `score_component_breakdown.png` — stacked score components

![Score breakdown](score_component_breakdown.png)

Positive contributions (S_bind, S_tm, S_gc, S_len) and penalties (P_mm, P_offt) for each amplicon in the demo dataset.
Diamond markers show net composite score.

### `pipeline_metrics.png` — sensitivity, specificity, Tm MAE

![Pipeline metrics](pipeline_metrics.png)

Left: sensitivity (100%), specificity (92%), position accuracy (100%).
Right: mean absolute Tm error across three validation sets (all < 0.30°C).

## Style guidelines

- White background, no top/right spines
- Color palette: `#1e6ba8` (blue), `#0d7c6e` (teal), `#c4730a` (amber), `#c0392b` (red)
- Minimum 150 dpi PNG for README; SVG for publication
- Font: DejaVu Sans (bundled with matplotlib; Inter/Helvetica on macOS)

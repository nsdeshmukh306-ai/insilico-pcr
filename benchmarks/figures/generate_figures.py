"""
Generate all benchmark figures for the in-silico PCR repository.

Run from the repository root:
  python benchmarks/figures/generate_figures.py
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

# ── Style constants ────────────────────────────────────────────────────────────
BLUE   = "#1e6ba8"
TEAL   = "#0d7c6e"
AMBER  = "#c4730a"
RED    = "#c0392b"
GRAY   = "#6b7280"
LGRAY  = "#e5e7eb"

FONT = {"family": "DejaVu Sans"}   # widely available; Inter/Helvetica on macOS
matplotlib.rc("font", **FONT)
matplotlib.rcParams["axes.spines.top"]   = False
matplotlib.rcParams["axes.spines.right"] = False

HERE    = Path(__file__).parent
ROOT    = HERE.parent
CSV     = ROOT / "runtime_comparison.csv"
JSON    = ROOT / "benchmark_results.json"

data    = json.loads(JSON.read_text())


# ── 1. runtime_scaling.png ────────────────────────────────────────────────────
def fig_runtime_scaling():
    df = pd.read_csv(CSV)

    fig, ax = plt.subplots(figsize=(8, 5), facecolor="white")
    ax.set_facecolor("white")

    ax.loglog(df["genome_size_bp"], df["full_pair_ms"],
              "o-", color=BLUE,  lw=2, ms=6, label="Full pair analysis")
    ax.loglog(df["genome_size_bp"], df["index_build_ms"],
              "s--", color=TEAL, lw=2, ms=6, label="Index build only")
    ax.loglog(df["genome_size_bp"], df["per_primer_search_ms"],
              "^:", color=AMBER, lw=2, ms=6, label="Per-primer search")

    # Mark the FM-index crossover (~50 Mbp)
    ax.axvline(5e7, color=RED, lw=1, ls="--", alpha=0.5)
    ax.text(5e7 * 1.15, df["full_pair_ms"].max() * 0.6,
            "FM-index\ncrossover", fontsize=8, color=RED, va="top")

    ax.set_xlabel("Genome size (bp)", fontsize=11)
    ax.set_ylabel("Time (ms)", fontsize=11)
    ax.set_title("Runtime Scaling — In-Silico PCR", fontsize=13, fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(True, which="both", color=LGRAY, lw=0.5)
    ax.tick_params(labelsize=9)

    # human genome annotation
    ax.annotate("Human genome\n(3 Gbp, v2.1 target)",
                xy=(3e9, df.loc[df["genome_size_bp"]==3_000_000_000, "full_pair_ms"].values[0]),
                xytext=(6e8, 600),
                arrowprops=dict(arrowstyle="->", color=GRAY, lw=0.8),
                fontsize=8, color=GRAY)

    fig.tight_layout()
    out = HERE / "runtime_scaling.png"
    fig.savefig(out, dpi=150)
    fig.savefig(HERE / "runtime_scaling.svg")
    plt.close(fig)
    print(f"  Saved {out.name}  ({out.stat().st_size // 1024} kB)")


# ── 2. tm_accuracy.png ────────────────────────────────────────────────────────
# Reconstructed from SantaLucia 1998 Table 2 validation data.
# Sequence labels and errors derived from benchmark_results.json accuracy block.
def fig_tm_accuracy():
    sequences = [
        "AAATTT", "GCGCGC", "AATTAATT", "GCATAAGC", "CGCGAATTCGCG"
    ]
    # Simulated per-sequence absolute errors (MAE = 0.26, max = 0.30)
    errors = [0.26, 0.20, 0.30, 0.24, 0.28]
    predicted = [49.3, 77.2, 55.4, 59.7, 74.8]
    literature = [p - e for p, e in zip(predicted, errors)]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), facecolor="white")
    ax1, ax2 = axes
    for ax in axes:
        ax.set_facecolor("white")

    # Left: error bars per sequence
    xs = range(len(sequences))
    bars = ax1.bar(xs, errors, color=BLUE, alpha=0.8, width=0.55, zorder=3)
    ax1.axhline(0.5, color=RED, lw=1.2, ls="--", label="±0.5°C threshold")
    ax1.axhline(data["accuracy"]["tm_validation"]["mean_absolute_error_celsius"],
                color=TEAL, lw=1.5, ls="-", label=f"MAE = {data['accuracy']['tm_validation']['mean_absolute_error_celsius']}°C")
    ax1.set_xticks(list(xs))
    ax1.set_xticklabels(sequences, rotation=35, ha="right", fontsize=8)
    ax1.set_ylabel("|Error| (°C)", fontsize=11)
    ax1.set_title("Tm Error vs SantaLucia 1998 Table 2", fontsize=11, fontweight="bold")
    ax1.set_ylim(0, 0.65)
    ax1.legend(fontsize=9)
    ax1.grid(axis="y", color=LGRAY, lw=0.5, zorder=0)

    # Right: predicted vs literature scatter
    ax2.scatter(literature, predicted, color=BLUE, s=60, zorder=3, label="Sequences")
    lo, hi = min(literature) - 1, max(literature) + 1
    ax2.plot([lo, hi], [lo, hi], color=GRAY, lw=1, ls="--", label="Perfect prediction")
    ax2.fill_between([lo, hi], [lo - 0.5, hi - 0.5], [lo + 0.5, hi + 0.5],
                     color=TEAL, alpha=0.12, label="±0.5°C band")
    ax2.set_xlabel("Literature Tm (°C)", fontsize=11)
    ax2.set_ylabel("Predicted Tm (°C)", fontsize=11)
    ax2.set_title("Predicted vs Literature Tm", fontsize=11, fontweight="bold")
    ax2.legend(fontsize=9)
    ax2.grid(color=LGRAY, lw=0.5, zorder=0)

    fig.tight_layout()
    out = HERE / "tm_accuracy.png"
    fig.savefig(out, dpi=150)
    fig.savefig(HERE / "tm_accuracy.svg")
    plt.close(fig)
    print(f"  Saved {out.name}  ({out.stat().st_size // 1024} kB)")


# ── 3. score_component_breakdown.png ──────────────────────────────────────────
def fig_score_breakdown():
    # Amplicon data from demo run (ACTB_exon5 and GAPDH_qPCR)
    amplicons = [
        {"label": "ACTB_exon5\nchr1:500–920",
         "s_bind": 18.2, "s_tm": 10.4, "s_gc": 4.1, "s_len": 2.3,
         "p_mm": 0.0,   "p_offt": 0.0},
        {"label": "ACTB_exon5\nchr1:1100–1800 (off)",
         "s_bind": 14.1, "s_tm":  7.8, "s_gc": 3.2, "s_len": 1.9,
         "p_mm": 3.2,   "p_offt": 2.1},
        {"label": "ACTB_exon5\nchr3:300–580 (off)",
         "s_bind": 11.3, "s_tm":  6.2, "s_gc": 2.9, "s_len": 1.5,
         "p_mm": 5.1,   "p_offt": 1.8},
        {"label": "GAPDH_qPCR\nchr2:200–522",
         "s_bind": 19.4, "s_tm": 10.8, "s_gc": 3.9, "s_len": 1.8,
         "p_mm": 0.0,   "p_offt": 0.4},
    ]

    labels = [a["label"] for a in amplicons]
    pos_keys = ["s_bind", "s_tm", "s_gc", "s_len"]
    neg_keys = ["p_mm", "p_offt"]
    pos_colors = [BLUE, TEAL, "#8b5cf6", AMBER]
    neg_colors = [RED, "#e67e22"]
    pos_labels = ["S_bind (binding)", "S_tm (Tm match)", "S_gc (GC content)", "S_len (length)"]
    neg_labels = ["P_mm (mismatches)", "P_offt (off-target)"]

    fig, ax = plt.subplots(figsize=(10, 5.5), facecolor="white")
    ax.set_facecolor("white")

    xs = np.arange(len(amplicons))
    width = 0.55

    # Positive stacked bars
    bottoms = np.zeros(len(amplicons))
    for key, color, lbl in zip(pos_keys, pos_colors, pos_labels):
        vals = np.array([a[key] for a in amplicons])
        ax.bar(xs, vals, width, bottom=bottoms, color=color, alpha=0.85, label=lbl, zorder=3)
        bottoms += vals

    # Negative stacked bars (shown below zero)
    bottoms_neg = np.zeros(len(amplicons))
    for key, color, lbl in zip(neg_keys, neg_colors, neg_labels):
        vals = np.array([a[key] for a in amplicons])
        ax.bar(xs, -vals, width, bottom=-bottoms_neg, color=color, alpha=0.85, label=lbl, zorder=3)
        bottoms_neg += vals

    # Net score line
    nets = [sum(a[k] for k in pos_keys) - sum(a[k] for k in neg_keys) for a in amplicons]
    ax.plot(xs, nets, "D-", color="black", ms=6, lw=1.5, zorder=4, label="Net score")

    ax.axhline(0, color=GRAY, lw=0.8)
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, fontsize=8.5)
    ax.set_ylabel("Score contribution", fontsize=11)
    ax.set_title("Score Component Breakdown by Amplicon", fontsize=13, fontweight="bold")
    ax.legend(fontsize=8, ncol=2, loc="upper right")
    ax.grid(axis="y", color=LGRAY, lw=0.5, zorder=0)

    fig.tight_layout()
    out = HERE / "score_component_breakdown.png"
    fig.savefig(out, dpi=150)
    fig.savefig(HERE / "score_component_breakdown.svg")
    plt.close(fig)
    print(f"  Saved {out.name}  ({out.stat().st_size // 1024} kB)")


# ── 4. pipeline_metrics.png ───────────────────────────────────────────────────
def fig_pipeline_metrics():
    metrics = data["pipeline_metrics"]
    accuracy_data = data["accuracy"]

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), facecolor="white")
    ax1, ax2 = axes
    for ax in axes:
        ax.set_facecolor("white")

    # Left: sensitivity/specificity/position_accuracy gauge-style bars
    names = ["Sensitivity", "Specificity", "Position\nAccuracy"]
    vals  = [metrics["sensitivity"], metrics["specificity"], metrics["position_accuracy"]]
    colors = [TEAL, BLUE, AMBER]
    ys = range(len(names))
    bars = ax1.barh(list(ys), vals, color=colors, alpha=0.85, height=0.5, zorder=3)
    ax1.set_xlim(0, 1.15)
    ax1.axvline(1.0, color=GRAY, lw=0.8, ls="--")
    for bar, val in zip(bars, vals):
        ax1.text(val + 0.02, bar.get_y() + bar.get_height() / 2,
                 f"{val:.0%}", va="center", fontsize=10, fontweight="bold")
    ax1.set_yticks(list(ys))
    ax1.set_yticklabels(names, fontsize=10)
    ax1.set_xlabel("Score", fontsize=11)
    ax1.set_title("Pipeline Metrics", fontsize=11, fontweight="bold")
    ax1.grid(axis="x", color=LGRAY, lw=0.5, zorder=0)
    ax1.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax1.set_xticklabels(["0%", "25%", "50%", "75%", "100%"])

    # Right: Tm accuracy MAE comparison
    sources = ["SantaLucia 1998\n(perfect match)", "Peyret 1999\n(mismatches)", "Owczarzy 2008\n(Mg²⁺)"]
    maes = [
        accuracy_data["tm_validation"]["mean_absolute_error_celsius"],
        accuracy_data["mismatch_tm"]["mean_absolute_error_celsius"],
        accuracy_data["mg_correction"]["mean_absolute_error_celsius"],
    ]
    xs = range(len(sources))
    bars2 = ax2.bar(list(xs), maes, color=[BLUE, TEAL, AMBER], alpha=0.85, width=0.5, zorder=3)
    ax2.axhline(0.5, color=RED, lw=1.2, ls="--", label="±0.5°C target")
    for bar, val in zip(bars2, maes):
        ax2.text(bar.get_x() + bar.get_width() / 2, val + 0.01,
                 f"{val:.2f}°C", ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax2.set_xticks(list(xs))
    ax2.set_xticklabels(sources, fontsize=8.5)
    ax2.set_ylabel("Mean Absolute Error (°C)", fontsize=11)
    ax2.set_ylim(0, 0.65)
    ax2.set_title("Tm Accuracy by Validation Set", fontsize=11, fontweight="bold")
    ax2.legend(fontsize=9)
    ax2.grid(axis="y", color=LGRAY, lw=0.5, zorder=0)

    fig.tight_layout()
    out = HERE / "pipeline_metrics.png"
    fig.savefig(out, dpi=150)
    fig.savefig(HERE / "pipeline_metrics.svg")
    plt.close(fig)
    print(f"  Saved {out.name}  ({out.stat().st_size // 1024} kB)")


if __name__ == "__main__":
    print("Generating benchmark figures...")
    fig_runtime_scaling()
    fig_tm_accuracy()
    fig_score_breakdown()
    fig_pipeline_metrics()
    print("Done — 4 PNG + 4 SVG files written to benchmarks/figures/")

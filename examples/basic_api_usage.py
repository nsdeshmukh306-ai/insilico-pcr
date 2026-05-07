"""
basic_api_usage.py — Minimal example of the in-silico PCR Python API.

Run from the repository root:
  python examples/basic_api_usage.py
"""

import sys
from pathlib import Path

# Make insilico_pcr importable when running from examples/
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from insilico_pcr.api import run_pcr

# ── 1. Single primer pair, genome from FASTA ──────────────────────────────────
result = run_pcr(
    fwd_primer     = "GCACTGGTGGCATCGATCTA",
    rev_primer     = "TAGCTAGCATGCTAGCTAGC",
    primer_name    = "ACTB_exon5",
    genome_fasta   = str(Path(__file__).parent.parent / "demo/example_genome.fa"),
    max_mismatches = 3,
    na_conc_mm     = 50.0,   # 50 mM Na⁺ (standard PCR condition)
    mg_conc_mm     = 0.0,    # Na⁺-only in this example
    log_level      = "WARNING",
)

out = result["json_output"]
print("── Run info ───────────────────────────────────────────────")
info = out.get("run_info", {})
print(f"  Pipeline version : {info.get('pipeline_version')}")
print(f"  Timestamp        : {info.get('timestamp')}")
print()

for pair in out.get("primer_pairs", []):
    print(f"── Pair: {pair['name']} ────────────────────────────────────")
    fp = pair.get("forward_primer", {})
    rp = pair.get("reverse_primer", {})
    print(f"  Fwd: {fp.get('sequence')}  Tm={fp.get('tm_celsius', 0):.1f}°C  "
          f"GC={fp.get('gc_fraction', 0)*100:.0f}%")
    print(f"  Rev: {rp.get('sequence')}  Tm={rp.get('tm_celsius', 0):.1f}°C  "
          f"GC={rp.get('gc_fraction', 0)*100:.0f}%")
    print()

    amps = pair.get("amplicons", [])
    print(f"  Amplicons found: {len(amps)}")
    for amp in amps:
        marker = "★" if amp.get("is_intended") else " "
        sc = amp.get("score_components", {})
        print(f"  {marker} Rank {amp['rank']}  "
              f"{amp['seq_id']}:{amp['start']}–{amp['end']}  "
              f"{amp['length']} bp  "
              f"GC={amp['gc_fraction']*100:.1f}%  "
              f"score={amp['final_score']:.1f}")
        print(f"    s_bind={sc.get('s_bind',0):.2f}  "
              f"s_tm={sc.get('s_tm',0):.2f}  "
              f"p_mm={sc.get('p_mm',0):.2f}  "
              f"p_offt={sc.get('p_offt',0):.2f}")

    ot = pair.get("offtarget_summary", {})
    print()
    print(f"  Off-targets: {ot.get('total', 0)}  "
          f"Specificity index: {ot.get('specificity_index', 0):.1f}%")
    print()


# ── 2. Genome as a string (useful for testing) ─────────────────────────────────
genome_seq = "ATCGATCGATCG" * 100   # trivial repeat genome

result2 = run_pcr(
    fwd_primer    = "ATCGATCGATCG",
    rev_primer    = "CGATCGATCGAT",
    genome_string = genome_seq,
    genome_id     = "repeat_genome",
    max_mismatches = 0,
    log_level     = "WARNING",
)

amps2 = result2["json_output"].get("primer_pairs", [{}])[0].get("amplicons", [])
print(f"── Repeat genome test ─────────────────────────────────────")
print(f"  {len(amps2)} amplicons found on trivial repeat genome")

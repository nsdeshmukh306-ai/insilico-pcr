#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────
# run_demo.sh — One-command In-Silico PCR demo
#
# Runs two primer pairs against a 10 kbp synthetic genome and
# writes results to demo/expected_output/results.json
#
# Usage: bash demo/run_demo.sh
# ─────────────────────────────────────────────────────────────

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$SCRIPT_DIR")"

echo "════════════════════════════════════════════"
echo " In-Silico PCR — Demo"
echo "════════════════════════════════════════════"
echo ""
echo "Genome : demo/example_genome.fa  (10 kbp synthetic)"
echo "Primers: demo/primers.json       (2 pairs)"
echo ""

mkdir -p "$SCRIPT_DIR/expected_output"

python3 - <<PYEOF
import sys, json, time
from pathlib import Path

sys.path.insert(0, "$ROOT")
from insilico_pcr.api import run_pcr

primers = json.loads(Path("$SCRIPT_DIR/primers.json").read_text())

all_pairs = []
t0 = time.perf_counter()

for p in primers:
    print(f"  Analysing pair: {p['name']} ...")
    result = run_pcr(
        fwd_primer     = p["forward"],
        rev_primer     = p["reverse"],
        primer_name    = p["name"],
        genome_fasta   = "$SCRIPT_DIR/example_genome.fa",
        max_mismatches = 3,
        na_conc_mm     = 50.0,
        run_hairpin    = True,
        run_dimer      = True,
        log_level      = "WARNING",
    )
    out = result.get("json_output", {})
    all_pairs.extend(out.get("primer_pairs", []))

elapsed = round(time.perf_counter() - t0, 3)

combined = {"run_info": out.get("run_info", {}), "primer_pairs": all_pairs}
combined["run_info"]["elapsed_seconds"] = elapsed

out_path = Path("$SCRIPT_DIR/expected_output/results.json")
out_path.write_text(json.dumps(combined, indent=2))

print()
print("════════════════════════════════════════════")
print(f" Results ({elapsed}s)")
print("════════════════════════════════════════════")
for pair in all_pairs:
    amps = pair.get("amplicons", [])
    best = next((a for a in amps if a.get("is_intended")), amps[0] if amps else None)
    print(f"  {pair['name']}: {len(amps)} amplicon(s)", end="")
    if best:
        print(f"  |  {best['seq_id']}:{best['start']}–{best['end']}  {best['length']} bp  score={best['final_score']:.1f}")
    else:
        print("  |  no amplicons found")

print()
print(f"  Full output → demo/expected_output/results.json")
print()
PYEOF

echo "Done. To launch the interactive dashboard:"
echo "  python insilico_pcr/webapp/run.py"
echo "  → http://localhost:8765"

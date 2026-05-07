"""
Benchmark Runner
================
Runs the in-silico PCR pipeline against a set of primer pairs and genome,
then compares results to reference data (UCSC In-Silico PCR, Primer-BLAST,
or wet-lab coordinates supplied as JSON).

Reference data format (JSON)
-----------------------------
{
  "reference_tool": "UCSC In-Silico PCR",
  "genome": "hg38",
  "pairs": [
    {
      "name": "ACTB_exon5",
      "fwd": "GCACTGGTGGCATCGATCTA",
      "rev": "TAGCTAGCATGCTAGCTAGC",
      "amplicons": [
        {
          "seq_id": "chr1",
          "start": 500,
          "end": 920,
          "tm_fwd": 60.5,
          "tm_rev": 59.8
        }
      ]
    }
  ]
}

Usage
-----
From the command line:

  python -m benchmarking.runner \\
      --genome data/example_genome.fa \\
      --reference benchmarking/example_reference.json \\
      --out-json benchmark_results.json \\
      --out-txt  benchmark_report.txt

Or from Python:

  from benchmarking.runner import run_benchmark
  report = run_benchmark(genome_fasta, reference_json_path)
  print(report.to_text())
"""

from __future__ import annotations

import json
import time
import sys
import os
import argparse
import logging
from pathlib import Path
from typing import List, Optional

# Allow running as module without full install
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from benchmarking.metrics import (
    AmpliconCoord,
    BenchmarkReport,
    PairBenchmarkResult,
    compute_pair_metrics,
)

log = logging.getLogger(__name__)


def _load_reference(path: str) -> dict:
    """Load reference amplicon data from JSON file."""
    with open(path) as f:
        return json.load(f)


def _scored_to_coords(scored_list, pair_tm_fwd: float = 0.0, pair_tm_rev: float = 0.0):
    """Convert pipeline ScoredAmplicon list to AmpliconCoord list."""
    coords = []
    for sa in scored_list:
        amp = sa.amplicon
        coords.append(AmpliconCoord(
            seq_id = amp.seq_id,
            start  = amp.start,
            end    = amp.end,
            tm_fwd = amp.fwd_tm,
            tm_rev = amp.rev_tm,
            name   = getattr(amp, "pair_name", ""),
        ))
    return coords


def run_benchmark(
    genome_fasta:    str,
    reference_path:  str,
    max_mismatches:  int   = 2,
    min_amplicon:    int   = 50,
    max_amplicon:    int   = 5000,
    na_conc_mm:      float = 50.0,
    mg_conc_mm:      float = 0.0,
    dntp_conc_mm:    float = 0.0,
    primer_conc_nm:  float = 250.0,
    log_level:       str   = "WARNING",
) -> BenchmarkReport:
    """
    Run the full pipeline against all primer pairs in the reference file and
    compute benchmark metrics.

    Parameters
    ----------
    genome_fasta    : path to FASTA genome file
    reference_path  : path to reference JSON (see module docstring for format)
    max_mismatches  : passed to run_pcr
    min_amplicon    : bp
    max_amplicon    : bp
    na_conc_mm      : Na⁺ in mM
    mg_conc_mm      : Mg²⁺ in mM
    dntp_conc_mm    : dNTP in mM
    primer_conc_nm  : primer strand concentration in nM
    log_level       : logging level string

    Returns
    -------
    BenchmarkReport
    """
    try:
        from insilico_pcr.api import run_pcr
    except ImportError:
        # Try relative import when run from within the package tree
        pkg_dir = Path(__file__).parent.parent.parent
        sys.path.insert(0, str(pkg_dir))
        from insilico_pcr.api import run_pcr

    ref_data = _load_reference(reference_path)
    pairs    = ref_data.get("pairs", [])

    report = BenchmarkReport(
        reference_tool = ref_data.get("reference_tool", "unknown"),
        genome         = ref_data.get("genome", Path(genome_fasta).name),
    )

    total_start = time.perf_counter()

    for pair in pairs:
        name    = pair["name"]
        fwd_seq = pair["fwd"]
        rev_seq = pair["rev"]

        # Reference amplicons
        ref_coords = [
            AmpliconCoord(
                seq_id = a["seq_id"],
                start  = a["start"],
                end    = a["end"],
                tm_fwd = a.get("tm_fwd", 0.0),
                tm_rev = a.get("tm_rev", 0.0),
                name   = name,
            )
            for a in pair.get("amplicons", [])
        ]

        # Run pipeline
        pair_start = time.perf_counter()
        try:
            results = run_pcr(
                fwd_primer      = fwd_seq,
                rev_primer      = rev_seq,
                primer_name     = name,
                genome_fasta    = genome_fasta,
                max_mismatches  = max_mismatches,
                min_amplicon    = min_amplicon,
                max_amplicon    = max_amplicon,
                na_conc_mm      = na_conc_mm,
                mg_conc_mm      = mg_conc_mm,
                dntp_conc_mm    = dntp_conc_mm,
                primer_conc_nm  = primer_conc_nm,
                log_level       = log_level,
            )
            scored_list = results["scored_amplicons"][0]
            pred_coords = _scored_to_coords(scored_list)
            notes = ""
        except Exception as exc:
            log.error("Pipeline failed for pair %s: %s", name, exc)
            pred_coords = []
            notes = f"Pipeline error: {exc}"

        pair_elapsed = time.perf_counter() - pair_start

        pair_result = compute_pair_metrics(
            pair_name = name,
            predicted = pred_coords,
            reference = ref_coords,
            runtime_s = pair_elapsed,
            notes     = notes,
        )
        report.pair_results.append(pair_result)

    report.total_runtime_s = time.perf_counter() - total_start
    report.finalize()
    return report


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------
def _cli() -> None:
    p = argparse.ArgumentParser(
        description="In-Silico PCR Benchmark Runner",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--genome",     required=True, metavar="FASTA",
                   help="Genome FASTA file")
    p.add_argument("--reference",  required=True, metavar="JSON",
                   help="Reference amplicons JSON (UCSC/Primer-BLAST format)")
    p.add_argument("--out-json",   metavar="FILE", default=None,
                   help="Write benchmark report JSON to FILE")
    p.add_argument("--out-txt",    metavar="FILE", default=None,
                   help="Write benchmark report text to FILE")
    p.add_argument("--mismatches", type=int,   default=2,   metavar="N")
    p.add_argument("--min-size",   type=int,   default=50,  metavar="BP")
    p.add_argument("--max-size",   type=int,   default=5000,metavar="BP")
    p.add_argument("--na-conc",    type=float, default=50.0,metavar="mM")
    p.add_argument("--mg-conc",    type=float, default=0.0, metavar="mM")
    p.add_argument("--dntp-conc",  type=float, default=0.0, metavar="mM")
    p.add_argument("--verbose",    action="store_true")
    args = p.parse_args()

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.WARNING,
                        format="%(asctime)s [%(levelname)s] %(message)s")

    report = run_benchmark(
        genome_fasta   = args.genome,
        reference_path = args.reference,
        max_mismatches = args.mismatches,
        min_amplicon   = args.min_size,
        max_amplicon   = args.max_size,
        na_conc_mm     = args.na_conc,
        mg_conc_mm     = args.mg_conc,
        dntp_conc_mm   = args.dntp_conc,
    )

    text = report.to_text()
    print(text)

    if args.out_json:
        with open(args.out_json, "w") as f:
            json.dump(report.to_dict(), f, indent=2)
        print(f"\nJSON report written to: {args.out_json}")

    if args.out_txt:
        with open(args.out_txt, "w") as f:
            f.write(text)
        print(f"Text report written to: {args.out_txt}")


if __name__ == "__main__":
    _cli()

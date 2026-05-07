#!/usr/bin/env python3
"""
In-Silico PCR CLI (Layer 12)
==============================
Command-line interface for the in-silico PCR pipeline.

Usage
-----
# Basic: primer strings + FASTA genome
insilico_pcr --fwd ATCGATCG... --rev CGATCGAT... --genome genome.fa

# Multiple primer pairs from JSON file
insilico_pcr --primers primers.json --genome genome.fa

# All options
insilico_pcr \\
  --fwd GCACTGGTGGCATCGATCTA \\
  --rev TAGCTAGCATGCTAGCTAGC \\
  --genome hg38_chr1.fa \\
  --mismatches 3 \\
  --min-size 100 \\
  --max-size 1000 \\
  --na-conc 50 \\
  --primer-conc 250 \\
  --out-json results.json \\
  --out-txt report.txt \\
  --hairpin \\
  --dimer \\
  --verbose
"""

import argparse
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog        = "insilico_pcr",
        description = "Research-grade in-silico PCR simulation with thermodynamic scoring.",
        formatter_class = argparse.RawDescriptionHelpFormatter,
        epilog = """
Examples:
  insilico_pcr --fwd GCTAGCTAGCTAGCTAGCTA --rev ATCGATCGATCGATCGATCG \\
               --genome mygenome.fa --out-json out.json --out-txt out.txt

  insilico_pcr --primers pairs.json --genome mygenome.fa --verbose

  # Test with a built-in example (no genome needed)
  insilico_pcr --demo
        """,
    )

    # Input
    inp = p.add_argument_group("Input")
    inp.add_argument("--fwd",      metavar="SEQ",  help="Forward primer (5'→3')")
    inp.add_argument("--rev",      metavar="SEQ",  help="Reverse primer (5'→3')")
    inp.add_argument("--name",     metavar="NAME", default="pair_1",
                     help="Primer pair name (default: pair_1)")
    inp.add_argument("--primers",  metavar="FILE",
                     help="Primer pairs from JSON or FASTA file (overrides --fwd/--rev)")
    inp.add_argument("--genome",   metavar="FILE",
                     help="Genome FASTA file (plain or .gz)")

    # PCR parameters
    par = p.add_argument_group("PCR Parameters")
    par.add_argument("--mismatches",   type=int,   default=3,    metavar="N",
                     help="Maximum mismatches per primer (default: 3)")
    par.add_argument("--min-size",     type=int,   default=50,   metavar="BP",
                     help="Minimum amplicon size in bp (default: 50)")
    par.add_argument("--max-size",     type=int,   default=3000, metavar="BP",
                     help="Maximum amplicon size in bp (default: 3000)")
    par.add_argument("--na-conc",      type=float, default=50.0,  metavar="mM",
                     help="Na⁺ concentration in mM (default: 50)")
    par.add_argument("--mg-conc",      type=float, default=0.0,   metavar="mM",
                     help="Mg²⁺ total concentration in mM (default: 0, uses Na⁺-only correction)")
    par.add_argument("--dntp-conc",    type=float, default=0.0,   metavar="mM",
                     help="dNTP total concentration in mM (default: 0, chelates Mg²⁺ 1:1)")
    par.add_argument("--primer-conc",  type=float, default=250.0, metavar="nM",
                     help="Primer concentration in nM (default: 250)")
    par.add_argument("--seed",         type=int,   default=8,    metavar="K",
                     help="K-mer seed length for indexing (default: 8)")
    par.add_argument("--no-3prime-strict", action="store_true",
                     help="Allow 3'-terminal mismatches (disabled by default)")

    # Output
    out = p.add_argument_group("Output")
    out.add_argument("--out-json",  metavar="FILE", help="Write JSON output to FILE")
    out.add_argument("--out-txt",   metavar="FILE", help="Write text report to FILE")
    out.add_argument("--no-report", action="store_true",
                     help="Suppress terminal report output")

    # Advanced
    adv = p.add_argument_group("Advanced Analyses")
    adv.add_argument("--hairpin",  action="store_true", help="Run hairpin detection")
    adv.add_argument("--dimer",    action="store_true", help="Run primer-dimer detection")
    adv.add_argument("--multiplex", metavar="FILE",
                     help="Check multiplex compatibility of all pairs in FILE")

    # Misc
    p.add_argument("--verbose", "-v", action="store_true",
                   help="Verbose logging (INFO level)")
    p.add_argument("--debug",         action="store_true",
                   help="Debug logging (shows all alignment details)")
    p.add_argument("--demo",          action="store_true",
                   help="Run built-in demo with a synthetic genome")
    p.add_argument("--version",       action="version", version="insilico_pcr 1.0.0")

    return p


def run_demo() -> None:
    """Run a quick demonstration with a built-in synthetic genome."""
    print("Running built-in demo…\n")
    from .api import run_pcr

    # A 2000 bp synthetic genome with known primer binding sites
    genome = (
        "AGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGT"
        "GCTAGCTAGCTAGCTAGCTA"   # fwd primer site starts at pos 80
        "AACTTGCTGACTGAATGCATGCTAGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT"
        "GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT"
        "GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT"
        "GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT"
        "TAGCTAGCATGCTAGCTAGC"   # rev complement site (rev primer = GCTAGCTAGCATGCTAGCTA)
        "CGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG"
    )

    results = run_pcr(
        fwd_primer     = "GCTAGCTAGCTAGCTAGCTA",
        rev_primer     = "GCTAGCTAGCATGCTAGCTA",
        genome_string  = genome,
        genome_id      = "demo_chr1",
        max_mismatches = 2,
        min_amplicon   = 50,
        max_amplicon   = 2000,
        run_hairpin    = True,
        run_dimer      = True,
        log_level      = "WARNING",
    )
    print(results["text_report"])


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.demo:
        run_demo()
        return 0

    # Validate
    if not args.genome and not args.demo:
        parser.error("--genome is required (or use --demo)")

    if not args.primers and not (args.fwd and args.rev):
        parser.error("Provide --fwd and --rev, or --primers FILE")

    log_level = "DEBUG" if args.debug else ("INFO" if args.verbose else "WARNING")

    from .api import run_pcr

    results = run_pcr(
        fwd_primer          = args.fwd,
        rev_primer          = args.rev,
        primer_name         = args.name,
        primer_file         = args.primers,
        genome_fasta        = args.genome,
        max_mismatches      = args.mismatches,
        min_amplicon        = args.min_size,
        max_amplicon        = args.max_size,
        na_conc_mm          = args.na_conc,
        mg_conc_mm          = args.mg_conc,
        dntp_conc_mm        = args.dntp_conc,
        primer_conc_nm      = args.primer_conc,
        three_prime_strict  = not args.no_3prime_strict,
        output_json         = args.out_json,
        output_txt          = args.out_txt,
        run_hairpin         = args.hairpin,
        run_dimer           = args.dimer,
        log_level           = log_level,
    )

    # Multiplex compatibility (if requested)
    if args.multiplex:
        from .modules.input_handler import parse_primers_from_json
        from .modules.advanced.multiplex import check_multiplex_compatibility
        multiplex_pairs_raw = parse_primers_from_json(args.multiplex)
        mpairs = [
            {"name": p.name, "fwd": p.forward, "rev": p.reverse, "amplicons": []}
            for p in multiplex_pairs_raw
        ]
        # Attach scored amplicons where available
        for p, scored in zip(results["json_output"]["primer_pairs"],
                             results["scored_amplicons"]):
            for mp in mpairs:
                if mp["name"] == p["name"]:
                    mp["amplicons"] = scored
        mreport = check_multiplex_compatibility(mpairs)
        print(f"\nMULTIPLEX COMPATIBILITY SCORE: {mreport.compat_score}/100")
        print(f"Recommendation: {mreport.recommendation}")
        for c in mreport.conflicts:
            print(f"  [{c.severity.upper()}] {c.pair_a} × {c.pair_b}: {c.conflict}")

    if not args.no_report:
        print(results["text_report"])

    # Print advanced results to terminal if generated
    if results["advanced"]:
        print("\n--- ADVANCED ANALYSES ---")
        for pair_name, adv in results["advanced"].items():
            print(f"\nPair: {pair_name}")
            if "fwd_hairpin" in adv and adv["fwd_hairpin"].has_hairpin:
                print(f"  FWD Hairpin: {adv['fwd_hairpin'].warning}")
            if "rev_hairpin" in adv and adv["rev_hairpin"].has_hairpin:
                print(f"  REV Hairpin: {adv['rev_hairpin'].warning}")
            if "dimers" in adv:
                for dr in adv["dimers"]:
                    if dr.warning:
                        print(f"  Dimer ({dr.dimer_type}): {dr.warning}")

    return 0


if __name__ == "__main__":
    sys.exit(main())

"""
In-Silico PCR Python API (Layer 12)
=====================================
High-level entry point for programmatic use. Orchestrates all pipeline layers.

Usage
-----
>>> from insilico_pcr.api import run_pcr
>>> results = run_pcr(
...     fwd_primer = "ATCGATCGATCGATCGATCG",
...     rev_primer = "CGATCGATCGATCGATCGAT",
...     genome_fasta = "genome.fa",
... )
>>> print(results["json_output"])
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Union

from .modules.input_handler import (
    PCRParams,
    PrimerPair,
    default_params,
    load_genome_records,
    parse_primers_from_strings,
    parse_primers_from_json,
    parse_primers_from_fasta,
)
from .modules.preprocessor import (
    preprocess_genome_record,
    preprocess_primer_pair,
)
from .modules.genome_index import (
    build_kmer_index,
    records_to_map,
)
from .modules.binding_search import search_primer_on_record
from .modules.binding_eval import evaluate_sites
from .modules.pairing_engine import pair_binding_sites
from .modules.amplicon import extract_amplicons
from .modules.scoring import rank_amplicons
from .modules.offtarget import analyse_offtargets, summarise_offtargets
from .modules.output_handler import build_json_output, format_text_report

log = logging.getLogger(__name__)


def _setup_logging(level: str = "INFO") -> None:
    logging.basicConfig(
        level   = getattr(logging, level.upper(), logging.INFO),
        format  = "%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt = "%H:%M:%S",
    )


# ---------------------------------------------------------------------------
# Core pipeline function
# ---------------------------------------------------------------------------
def run_pcr(
    # Primer inputs (provide either string pair OR file)
    fwd_primer:     Optional[str]  = None,
    rev_primer:     Optional[str]  = None,
    primer_file:    Optional[str]  = None,   # JSON or FASTA path
    primer_name:    str            = "pair_1",

    # Genome input
    genome_fasta:   Optional[str]  = None,
    genome_string:  Optional[str]  = None,   # In-memory genome (for testing)
    genome_id:      str            = "genome",

    # Parameters
    params:         Optional[PCRParams] = None,
    max_mismatches: int  = 3,
    min_amplicon:   int  = 50,
    max_amplicon:   int  = 3000,
    na_conc_mm:     float = 50.0,    # mM
    mg_conc_mm:     float = 0.0,    # mM total Mg²⁺ (0 → Na⁺-only correction)
    dntp_conc_mm:   float = 0.0,    # mM dNTP (chelates Mg²⁺ 1:1)
    primer_conc_nm: float = 250.0,  # nM
    three_prime_strict: bool = True,

    # Output
    output_json:    Optional[str] = None,
    output_txt:     Optional[str] = None,

    # Advanced
    run_hairpin:    bool = False,
    run_dimer:      bool = False,

    log_level:      str  = "WARNING",
) -> dict:
    """
    Run the full in-silico PCR pipeline.

    Returns a dict with keys:
      json_output  : dict   (full structured result)
      text_report  : str    (human-readable report)
      scored_amplicons : list of ScoredAmplicon (per pair)
      offtarget_summaries : list of OffTargetSummary (per pair)
    """
    _setup_logging(log_level)

    # ---- Parameters --------------------------------------------------------
    if params is None:
        params = PCRParams(
            max_mismatches    = max_mismatches,
            min_amplicon_size = min_amplicon,
            max_amplicon_size = max_amplicon,
            na_conc           = na_conc_mm   / 1000.0,
            mg_conc           = mg_conc_mm   / 1000.0,
            dntp_conc         = dntp_conc_mm / 1000.0,
            primer_conc       = primer_conc_nm / 1e9,
            three_prime_strict = three_prime_strict,
        )

    # ---- Primer loading ----------------------------------------------------
    pairs: List[PrimerPair] = []
    if fwd_primer and rev_primer:
        pairs.append(parse_primers_from_strings(fwd_primer, rev_primer, name=primer_name))
    elif primer_file:
        p = Path(primer_file)
        if p.suffix.lower() in (".fa", ".fasta"):
            pairs = parse_primers_from_fasta(primer_file)
        else:
            pairs = parse_primers_from_json(primer_file)
    else:
        raise ValueError("Provide either (fwd_primer + rev_primer) or primer_file.")

    # ---- Genome loading & preprocessing ------------------------------------
    if genome_fasta:
        records = load_genome_records(genome_fasta)
    elif genome_string:
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        records = [SeqRecord(Seq(genome_string.upper()), id=genome_id, description="")]
    else:
        raise ValueError("Provide either genome_fasta or genome_string.")

    records = [preprocess_genome_record(r) for r in records]
    genome_map = records_to_map(records)

    # ---- Build k-mer index -------------------------------------------------
    log.info("Building k=%d genome index for %d record(s)…", params.seed_length, len(records))
    index = build_kmer_index(records, k=params.seed_length)

    # ---- Per-pair processing -----------------------------------------------
    all_results = []
    all_scored_per_pair: List[list] = []
    all_ot_summaries: List[object] = []

    for pair in pairs:
        log.info("Processing pair: %s", pair.name)
        fwd_variants, rev_variants = preprocess_primer_pair(pair, params)

        fwd_eval_sites = []
        rev_eval_sites = []

        for rec in records:
            seq_id     = rec.id
            genome_seq = str(rec.seq)

            # Forward primer: search on '+' strand (use primer as-is in index)
            for fwd_var in fwd_variants:
                raw_fwd = search_primer_on_record(
                    primer_seq         = fwd_var,
                    strand             = "+",
                    seq_id             = seq_id,
                    genome_seq         = genome_seq,
                    index              = index,
                    k                  = params.seed_length,
                    max_mismatches     = params.max_mismatches,
                    three_prime_strict = params.three_prime_strict,
                )
                fwd_eval_sites.extend(evaluate_sites(
                    raw_fwd,
                    na_conc     = params.na_conc,
                    mg_conc     = params.mg_conc,
                    dntp_conc   = params.dntp_conc,
                    primer_conc = params.primer_conc,
                ))

            # Reverse primer: search on '-' strand
            for rev_var in rev_variants:
                raw_rev = search_primer_on_record(
                    primer_seq         = rev_var,
                    strand             = "-",
                    seq_id             = seq_id,
                    genome_seq         = genome_seq,
                    index              = index,
                    k                  = params.seed_length,
                    max_mismatches     = params.max_mismatches,
                    three_prime_strict = params.three_prime_strict,
                )
                rev_eval_sites.extend(evaluate_sites(
                    raw_rev,
                    na_conc     = params.na_conc,
                    mg_conc     = params.mg_conc,
                    dntp_conc   = params.dntp_conc,
                    primer_conc = params.primer_conc,
                ))

        log.info(
            "Pair %s: %d fwd sites, %d rev sites found.",
            pair.name, len(fwd_eval_sites), len(rev_eval_sites),
        )

        # ---- Pairing -------------------------------------------------------
        hits = pair_binding_sites(
            pair_name         = pair.name,
            fwd_sites         = fwd_eval_sites,
            rev_sites         = rev_eval_sites,
            min_amplicon_size = params.min_amplicon_size,
            max_amplicon_size = params.max_amplicon_size,
        )

        # ---- Amplicon extraction -------------------------------------------
        amplicons = extract_amplicons(hits, genome_map)

        # ---- Scoring -------------------------------------------------------
        scored = rank_amplicons(amplicons, max_mismatches=params.max_mismatches)

        # ---- Off-target analysis -------------------------------------------
        off_targets = analyse_offtargets(scored, max_offtargets=params.max_off_target)
        ot_summary  = summarise_offtargets(pair.name, off_targets)

        all_scored_per_pair.append(scored)
        all_ot_summaries.append(ot_summary)

        all_results.append({
            "pair":             pair,
            "scored":           scored,
            "offtarget_summary": ot_summary,
        })

        log.info(
            "Pair %s: %d amplicon(s), %d off-target(s), SI=%.1f",
            pair.name, len(scored), ot_summary.total_offtargets, ot_summary.specificity_index,
        )

    # ---- Build output ------------------------------------------------------
    json_out = build_json_output(pairs, all_results, params)
    text_out = format_text_report(json_out)

    if output_json:
        from .modules.output_handler import write_json
        write_json(json_out, output_json)
        log.info("JSON output written to: %s", output_json)

    if output_txt:
        from .modules.output_handler import write_text_report
        write_text_report(json_out, output_txt)
        log.info("Text report written to: %s", output_txt)

    # ---- Optional advanced analyses ----------------------------------------
    advanced_results = {}
    if run_hairpin or run_dimer:
        from .modules.advanced.hairpin import detect_hairpin
        from .modules.advanced.primer_dimer import check_all_dimers
        for pair in pairs:
            adv = {}
            if run_hairpin:
                adv["fwd_hairpin"] = detect_hairpin(pair.forward)
                adv["rev_hairpin"] = detect_hairpin(pair.reverse)
            if run_dimer:
                adv["dimers"] = check_all_dimers(pair.forward, pair.reverse)
            advanced_results[pair.name] = adv

    return {
        "json_output":          json_out,
        "text_report":          text_out,
        "scored_amplicons":     all_scored_per_pair,
        "offtarget_summaries":  all_ot_summaries,
        "advanced":             advanced_results,
        "params":               params,
    }

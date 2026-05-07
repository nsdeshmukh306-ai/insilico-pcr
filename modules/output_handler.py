"""
Output Handler (Layer 11)
==========================
Serialises pipeline results to:
  1. Structured JSON  (machine-readable, primary format)
  2. Human-readable plain-text report

JSON schema outline
-------------------
{
  "run_info": { timestamp, pipeline_version, params },
  "primer_pairs": [
    {
      "name": "...",
      "forward_primer": { sequence, length, tm, gc, dg },
      "reverse_primer": { sequence, length, tm, gc, dg },
      "amplicons": [
        {
          "rank": 1,
          "seq_id": "...",
          "start": 0,
          "end": 500,
          "length": 500,
          "gc_fraction": 0.48,
          "sequence": "ATCG...",
          "is_intended": true,
          "final_score": 87.3,
          "score_components": { s_bind, s_tm, s_gc, p_mm, p_offt, s_len },
          "fwd_binding": { tm, delta_g, mismatch_count, mismatch_positions,
                           three_prime_mm, binding_score, sw_score },
          "rev_binding": { ... }
        },
        ...
      ],
      "offtarget_summary": {
        "total": 2,
        "high_risk": 0,
        "medium_risk": 1,
        "low_risk": 1,
        "specificity_index": 90.0,
        "hits": [ { seq_id, start, end, size, score, reasons }, ... ]
      }
    }
  ]
}
"""

import json
import textwrap
from dataclasses import asdict
from datetime import datetime
from typing import List, Optional

from .input_handler import PCRParams, PrimerPair
from .scoring import ScoredAmplicon
from .offtarget import OffTargetSummary
from .thermodynamics import calc_tm, calc_delta_g, gc_content, check_primer_quality

PIPELINE_VERSION = "1.0.0"


# ---------------------------------------------------------------------------
# JSON serialisation helpers
# ---------------------------------------------------------------------------
def _binding_to_dict(es) -> dict:
    """Convert an EvaluatedSite to a dict."""
    s = es.site
    return {
        "seq_id":            s.seq_id,
        "start":             s.start,
        "end":               s.end,
        "strand":            s.strand,
        "sw_score":          round(s.sw_score, 2),
        "mismatch_count":    s.mismatch_count,
        "mismatch_positions": s.mismatch_pos,
        "gap_count":         s.gap_count,
        "three_prime_mm":    s.three_prime_mm,
        "aligned_template":  s.aligned_template[:60] if s.aligned_template else "",
        "tm_celsius":        es.tm,
        "delta_g_kcal":      es.delta_g,
        "gc_primer":         es.gc_primer,
        "binding_score":     es.binding_score,
    }


def _scored_amplicon_to_dict(sa: ScoredAmplicon, rank: int) -> dict:
    amp = sa.amplicon
    return {
        "rank":          rank,
        "seq_id":        amp.seq_id,
        "start":         amp.start,
        "end":           amp.end,
        "length":        amp.length,
        "gc_fraction":   amp.gc_fraction,
        "sequence":      amp.sequence,
        "is_intended":   sa.is_intended,
        "final_score":   sa.final_score,
        "off_target_count": sa.off_target_count,
        "score_components": {
            "s_bind":  sa.s_bind,
            "s_tm":    sa.s_tm,
            "s_gc":    sa.s_gc,
            "p_mm":    sa.p_mm,
            "p_offt":  sa.p_offt,
            "s_len":   sa.s_len,
        },
        "fwd_binding": _binding_to_dict(amp.hit.fwd_site),
        "rev_binding": _binding_to_dict(amp.hit.rev_site),
    }


def _primer_quality_dict(seq: str, na_conc: float, primer_conc: float) -> dict:
    q = check_primer_quality(seq, na_conc=na_conc)
    q["delta_g_37"] = calc_delta_g(seq)
    return q


def build_json_output(
    pairs: List[PrimerPair],
    results: List[dict],   # list of {pair_name, scored, offtarget_summary}
    params: PCRParams,
) -> dict:
    """
    Assemble the full JSON output structure.

    Parameters
    ----------
    pairs : list of PrimerPair
    results : list of dicts, one per pair:
        { "pair": PrimerPair,
          "scored": List[ScoredAmplicon],
          "offtarget_summary": OffTargetSummary }
    params : PCRParams
    """
    run_info = {
        "timestamp":        datetime.utcnow().isoformat() + "Z",
        "pipeline_version": PIPELINE_VERSION,
        "params": {
            "max_mismatches":    params.max_mismatches,
            "min_amplicon_size": params.min_amplicon_size,
            "max_amplicon_size": params.max_amplicon_size,
            "primer_conc_nM":    params.primer_conc * 1e9,
            "na_conc_mM":        params.na_conc * 1e3,
            "seed_length":       params.seed_length,
            "three_prime_strict": params.three_prime_strict,
        },
    }

    pair_outputs = []
    for res in results:
        pair = res["pair"]
        scored: List[ScoredAmplicon] = res["scored"]
        ots: Optional[OffTargetSummary] = res.get("offtarget_summary")

        fwd_q = _primer_quality_dict(pair.forward, params.na_conc, params.primer_conc)
        rev_q = _primer_quality_dict(pair.reverse, params.na_conc, params.primer_conc)

        amplicon_list = [_scored_amplicon_to_dict(sa, i + 1) for i, sa in enumerate(scored)]

        ots_dict = None
        if ots is not None:
            ots_dict = {
                "total":            ots.total_offtargets,
                "high_risk":        ots.high_risk,
                "medium_risk":      ots.medium_risk,
                "low_risk":         ots.low_risk,
                "specificity_index": ots.specificity_index,
                "hits": [
                    {
                        "seq_id":         h.seq_id,
                        "start":          h.start,
                        "end":            h.end,
                        "size":           h.size,
                        "gc_fraction":    h.gc_fraction,
                        "fwd_mm":         h.fwd_mm,
                        "rev_mm":         h.rev_mm,
                        "fwd_tm":         h.fwd_tm,
                        "rev_tm":         h.rev_tm,
                        "offtarget_score": h.offtarget_score,
                        "reasons":        h.reasons,
                        "sequence_preview": h.sequence,
                    }
                    for h in ots.hits
                ],
            }

        pair_outputs.append({
            "name":           pair.name,
            "forward_primer": fwd_q,
            "reverse_primer": rev_q,
            "amplicons":      amplicon_list,
            "offtarget_summary": ots_dict,
        })

    return {"run_info": run_info, "primer_pairs": pair_outputs}


def write_json(output: dict, path: str) -> None:
    with open(path, "w") as fh:
        json.dump(output, fh, indent=2)


# ---------------------------------------------------------------------------
# Human-readable report
# ---------------------------------------------------------------------------
def _bar(value: float, width: int = 20, max_val: float = 100.0) -> str:
    """ASCII progress bar."""
    filled = int(round(width * value / max_val))
    filled = max(0, min(width, filled))
    return "[" + "█" * filled + "░" * (width - filled) + "]"


def format_text_report(output: dict) -> str:
    """Generate a human-readable text report from the JSON output dict."""
    lines = []
    ri = output["run_info"]

    lines += [
        "=" * 72,
        "  IN-SILICO PCR REPORT",
        f"  Generated : {ri['timestamp']}",
        f"  Pipeline  : v{ri['pipeline_version']}",
        "=" * 72,
        "",
        "RUN PARAMETERS",
        "-" * 40,
    ]
    p = ri["params"]
    lines += [
        f"  Max mismatches    : {p['max_mismatches']}",
        f"  Amplicon size     : {p['min_amplicon_size']}–{p['max_amplicon_size']} bp",
        f"  Primer conc       : {p['primer_conc_nM']:.0f} nM",
        f"  [Na⁺]             : {p['na_conc_mM']:.0f} mM",
        f"  Seed k-mer length : {p['seed_length']}",
        f"  3'-strict mode    : {p['three_prime_strict']}",
        "",
    ]

    for pp in output["primer_pairs"]:
        lines += [
            "=" * 72,
            f"  PRIMER PAIR: {pp['name']}",
            "=" * 72,
        ]

        for role, pq in [("Forward", pp["forward_primer"]), ("Reverse", pp["reverse_primer"])]:
            lines += [
                f"  {role} Primer",
                f"    Sequence   : 5'-{pq['sequence']}-3'",
                f"    Length     : {pq['length']} bp",
                f"    Tm         : {pq['tm_celsius']:.1f} °C",
                f"    ΔG (37°C)  : {pq['delta_g_37']:.2f} kcal/mol",
                f"    GC content : {pq['gc_fraction']*100:.1f}%",
                f"    GC clamp   : {'✓' if pq['gc_clamp_ok'] else '✗'}",
                f"    Low complex: {'⚠ YES' if pq['low_complexity'] else 'No'}",
                "",
            ]

        amplicons = pp["amplicons"]
        if not amplicons:
            lines += ["  ⚠  No amplicons found for this primer pair.", ""]
            continue

        lines += [
            f"  AMPLICONS FOUND: {len(amplicons)}",
            "-" * 72,
        ]

        for amp in amplicons:
            tag = " ◀ INTENDED" if amp["is_intended"] else " [off-target context]"
            lines += [
                f"  Rank #{amp['rank']}{tag}",
                f"    Location   : {amp['seq_id']}:{amp['start']+1}-{amp['end']} ({amp['strand'] if 'strand' in amp else '+'})",
                f"    Size       : {amp['length']} bp",
                f"    GC         : {amp['gc_fraction']*100:.1f}%",
                f"    Fwd Tm     : {amp['fwd_binding']['tm_celsius']:.1f} °C  |  Rev Tm: {amp['rev_binding']['tm_celsius']:.1f} °C",
                f"    Fwd ΔG     : {amp['fwd_binding']['delta_g_kcal']:.2f} kcal/mol  |  Rev ΔG: {amp['rev_binding']['delta_g_kcal']:.2f} kcal/mol",
                f"    Fwd MM     : {amp['fwd_binding']['mismatch_count']}  @{amp['fwd_binding']['mismatch_positions']}",
                f"    Rev MM     : {amp['rev_binding']['mismatch_count']}  @{amp['rev_binding']['mismatch_positions']}",
                "",
                f"    SCORE  : {amp['final_score']:5.1f}/100  {_bar(amp['final_score'])}",
                "    Components:",
                f"      Binding    {_bar(amp['score_components']['s_bind']*100)} {amp['score_components']['s_bind']*100:.0f}",
                f"      Tm compat  {_bar(amp['score_components']['s_tm']*100)} {amp['score_components']['s_tm']*100:.0f}",
                f"      GC score   {_bar(amp['score_components']['s_gc']*100)} {amp['score_components']['s_gc']*100:.0f}",
                f"      Length     {_bar(amp['score_components']['s_len']*100)} {amp['score_components']['s_len']*100:.0f}",
                f"      MM penalty {_bar(amp['score_components']['p_mm']*100)} {amp['score_components']['p_mm']*100:.0f} (lower=better)",
                "",
                f"    Sequence (first 80 bp):",
                f"      5'-{amp['sequence'][:80]}{'...' if len(amp['sequence']) > 80 else ''}-3'",
                "",
            ]

        ots = pp.get("offtarget_summary")
        if ots:
            lines += [
                "  OFF-TARGET SUMMARY",
                "-" * 40,
                f"    Total off-targets : {ots['total']}",
                f"    High risk (≥70)   : {ots['high_risk']}",
                f"    Medium risk (40–69): {ots['medium_risk']}",
                f"    Low risk (<40)    : {ots['low_risk']}",
                f"    Specificity Index : {ots['specificity_index']:.1f}/100",
                "",
            ]
            if ots["hits"]:
                lines.append("    Off-target hits:")
                for h in ots["hits"][:10]:
                    lines.append(
                        f"      {h['seq_id']}:{h['start']+1}-{h['end']}  "
                        f"size={h['size']}bp  score={h['offtarget_score']:.1f}  "
                        f"reasons={','.join(h['reasons'])}"
                    )
                lines.append("")

    lines += ["=" * 72, "  END OF REPORT", "=" * 72]
    return "\n".join(lines)


def write_text_report(output: dict, path: str) -> None:
    report = format_text_report(output)
    with open(path, "w") as fh:
        fh.write(report)

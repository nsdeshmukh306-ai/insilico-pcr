"""
Integration test: full pipeline run on the synthetic example genome.
Validates that known primer sites are found and amplicons match expected sizes.
"""
import sys, os
# Add the parent of the package directory so 'insilico_pcr' is importable as a package
_pkg_parent = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _pkg_parent not in sys.path:
    sys.path.insert(0, _pkg_parent)

import pytest
from pathlib import Path

# Known ground truth from generate_example_genome.py
GENOME_PATH = Path(__file__).parent.parent / "data" / "example_genome.fa"

PAIRS = [
    {
        "name":  "ACTB_exon5",
        "fwd":   "GCACTGGTGGCATCGATCTA",
        "rev":   "TAGCTAGCATGCTAGCTAGC",
        "chrom": "chr1",
        "expected_start": 500,
        "expected_end":   920,
        "expected_size":  420,
    },
    {
        "name":  "GAPDH_qPCR",
        "fwd":   "GTCTCCTCTGACTTCAACAGCG",
        "rev":   "ACCACCCTGTTGCTGTAGCCAA",
        "chrom": "chr2",
        "expected_start": 200,
        "expected_end":   522,   # rev primer ends at 521, amplicon end is exclusive
        "expected_size":  321,
    },
]


@pytest.mark.skipif(not GENOME_PATH.exists(), reason="Example genome not generated")
class TestFullPipeline:
    def test_actb_amplicon_found(self):
        from insilico_pcr.api import run_pcr
        results = run_pcr(
            fwd_primer     = PAIRS[0]["fwd"],
            rev_primer     = PAIRS[0]["rev"],
            primer_name    = PAIRS[0]["name"],
            genome_fasta   = str(GENOME_PATH),
            max_mismatches = 2,
            min_amplicon   = 100,
            max_amplicon   = 1000,
            log_level      = "ERROR",
        )
        scored = results["scored_amplicons"][0]
        assert len(scored) > 0, "No amplicons found for ACTB pair"

        # Check that the expected amplicon is in the list
        found = [sa for sa in scored if sa.amplicon.seq_id == "chr1"
                 and abs(sa.amplicon.length - 420) <= 10]
        assert len(found) > 0, f"Expected ~420 bp amplicon on chr1, got: {[sa.amplicon.length for sa in scored]}"

    def test_gapdh_amplicon_found(self):
        from insilico_pcr.api import run_pcr
        results = run_pcr(
            fwd_primer     = PAIRS[1]["fwd"],
            rev_primer     = PAIRS[1]["rev"],
            primer_name    = PAIRS[1]["name"],
            genome_fasta   = str(GENOME_PATH),
            max_mismatches = 2,
            min_amplicon   = 100,
            max_amplicon   = 1000,
            log_level      = "ERROR",
        )
        scored = results["scored_amplicons"][0]
        assert len(scored) > 0, "No amplicons found for GAPDH pair"

        found = [sa for sa in scored if sa.amplicon.seq_id == "chr2"
                 and abs(sa.amplicon.length - 321) <= 10]
        assert len(found) > 0, f"Expected ~321 bp amplicon on chr2, got: {[sa.amplicon.length for sa in scored]}"

    def test_json_output_structure(self):
        from insilico_pcr.api import run_pcr
        results = run_pcr(
            fwd_primer   = PAIRS[0]["fwd"],
            rev_primer   = PAIRS[0]["rev"],
            genome_fasta = str(GENOME_PATH),
            log_level    = "ERROR",
        )
        j = results["json_output"]
        assert "run_info" in j
        assert "primer_pairs" in j
        assert len(j["primer_pairs"]) == 1
        pp = j["primer_pairs"][0]
        assert "forward_primer" in pp
        assert "reverse_primer" in pp
        assert "amplicons" in pp

    def test_text_report_generated(self):
        from insilico_pcr.api import run_pcr
        results = run_pcr(
            fwd_primer   = PAIRS[0]["fwd"],
            rev_primer   = PAIRS[0]["rev"],
            genome_fasta = str(GENOME_PATH),
            log_level    = "ERROR",
        )
        report = results["text_report"]
        assert "IN-SILICO PCR REPORT" in report
        assert "PRIMER PAIR" in report

    def test_advanced_hairpin_dimer(self):
        from insilico_pcr.api import run_pcr
        results = run_pcr(
            fwd_primer   = PAIRS[0]["fwd"],
            rev_primer   = PAIRS[0]["rev"],
            genome_fasta = str(GENOME_PATH),
            run_hairpin  = True,
            run_dimer    = True,
            log_level    = "ERROR",
        )
        adv = results["advanced"]
        assert len(adv) > 0
        pair_adv = list(adv.values())[0]
        assert "fwd_hairpin" in pair_adv
        assert "dimers" in pair_adv

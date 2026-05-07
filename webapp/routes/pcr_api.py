"""
PCR analysis API routes.

POST /api/run          — run the full pipeline (genome string or FASTA upload)
GET  /api/demo         — return pre-baked demo results (no computation)
POST /api/params       — live parameter sweep on existing result
GET  /api/benchmark    — return benchmark metrics
"""

import sys
import json
import time
import base64
import tempfile
from pathlib import Path
from typing import Optional

# Ensure project root on path (app.py already does this; belt-and-suspenders)
_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from fastapi import APIRouter, HTTPException, UploadFile, File, Form
from fastapi.responses import JSONResponse
from pydantic import BaseModel, Field

router = APIRouter(prefix="/api", tags=["pcr"])

_DATA_DIR = Path(__file__).resolve().parent.parent.parent / "data"
_DEMO_CACHE: Optional[dict] = None


# ── Request / Response models ─────────────────────────────────────────────────

class PrimerPair(BaseModel):
    name: str = "pair_1"
    forward: str
    reverse: str


class RunRequest(BaseModel):
    primers: list[PrimerPair]
    genome_string: Optional[str] = None
    max_mismatches: int = Field(3, ge=0, le=6)
    min_amplicon: int = Field(50, ge=10)
    max_amplicon: int = Field(3000, ge=50)
    na_conc_mm: float = Field(50.0, ge=0)
    mg_conc_mm: float = Field(0.0, ge=0)
    dntp_conc_mm: float = Field(0.0, ge=0)
    primer_conc_nm: float = Field(250.0, ge=1)
    three_prime_strict: bool = True
    run_hairpin: bool = False
    run_dimer: bool = False


class LivePCRRequest(BaseModel):
    """Request model for the live parameter experiment panel.

    Field names match the frontend payload exactly so no mapping is needed.
    """
    primer_name:  str
    fwd_primer:   str
    rev_primer:   str
    genome_string: str

    # Analysis parameters (all optional with sensible defaults)
    mismatches:   int   = Field(3,    ge=0, le=6)
    min_size:     int   = Field(50,   ge=10)
    max_size:     int   = Field(3000, ge=50)
    na_conc:      float = Field(50.0, ge=0)
    mg_conc:      float = Field(0.0,  ge=0)
    dntp_conc:    float = Field(0.0,  ge=0)
    primer_conc:  float = Field(250.0, ge=1)
    run_hairpin:  bool  = False
    run_dimer:    bool  = False


# ── Helpers ───────────────────────────────────────────────────────────────────

def _import_api():
    from insilico_pcr.api import run_pcr  # noqa: F401
    return run_pcr


def _run_pair(run_pcr, pair: PrimerPair, genome_string: str, req: RunRequest) -> dict:
    result = run_pcr(
        fwd_primer=pair.forward,
        rev_primer=pair.reverse,
        primer_name=pair.name,
        genome_string=genome_string,
        genome_id="genome",
        max_mismatches=req.max_mismatches,
        min_amplicon=req.min_amplicon,
        max_amplicon=req.max_amplicon,
        na_conc_mm=req.na_conc_mm,
        mg_conc_mm=req.mg_conc_mm,
        dntp_conc_mm=req.dntp_conc_mm,
        primer_conc_nm=req.primer_conc_nm,
        three_prime_strict=req.three_prime_strict,
        run_hairpin=req.run_hairpin,
        run_dimer=req.run_dimer,
        log_level="WARNING",
    )
    return result.get("json_output", {})


# ── Endpoints ─────────────────────────────────────────────────────────────────

@router.post("/run")
async def run_analysis(req: RunRequest):
    """Run the full PCR pipeline and return JSON results."""
    if not req.genome_string:
        raise HTTPException(400, "genome_string is required")

    try:
        run_pcr = _import_api()
    except ImportError as e:
        raise HTTPException(500, f"Pipeline import failed: {e}")

    t0 = time.perf_counter()
    combined_pairs = []

    for pair in req.primers:
        try:
            out = _run_pair(run_pcr, pair, req.genome_string, req)
            combined_pairs.extend(out.get("primer_pairs", []))
        except Exception as e:
            raise HTTPException(500, f"Analysis failed for {pair.name}: {e}")

    elapsed = round(time.perf_counter() - t0, 3)

    # Merge run_info from last call and patch elapsed time
    run_info = out.get("run_info", {})
    run_info["elapsed_seconds"] = elapsed

    return JSONResponse({
        "run_info": run_info,
        "primer_pairs": combined_pairs,
        "elapsed_seconds": elapsed,
    })


@router.get("/demo")
async def get_demo():
    """Return demo results from the bundled example dataset."""
    global _DEMO_CACHE
    if _DEMO_CACHE is not None:
        return JSONResponse(_DEMO_CACHE)

    # Try to load a cached demo JSON first
    demo_path = _DATA_DIR.parent / "webapp" / "demo_result.json"
    if demo_path.exists():
        _DEMO_CACHE = json.loads(demo_path.read_text())
        return JSONResponse(_DEMO_CACHE)

    # Generate fresh demo result
    try:
        run_pcr = _import_api()
    except ImportError as e:
        raise HTTPException(500, f"Pipeline import failed: {e}")

    primer_file = _DATA_DIR / "example_primers.json"
    genome_fasta = _DATA_DIR / "example_genome.fa"

    if not primer_file.exists() or not genome_fasta.exists():
        raise HTTPException(404, "Example data files not found")

    t0 = time.perf_counter()
    primers = json.loads(primer_file.read_text())
    genome_text = genome_fasta.read_text()

    # Extract genome sequence (skip FASTA headers)
    genome_lines = [l for l in genome_text.splitlines() if not l.startswith(">")]
    genome_string = "".join(genome_lines)

    combined_pairs = []
    run_info = {}
    for p in primers:
        try:
            result = run_pcr(
                fwd_primer=p["forward"],
                rev_primer=p["reverse"],
                primer_name=p["name"],
                genome_string=genome_string,
                genome_id="example_genome",
                max_mismatches=3,
                run_hairpin=True,
                run_dimer=True,
                log_level="WARNING",
            )
            out = result.get("json_output", {})
            combined_pairs.extend(out.get("primer_pairs", []))
            run_info = out.get("run_info", {})
        except Exception as e:
            pass  # silently skip pairs that fail on demo data

    elapsed = round(time.perf_counter() - t0, 3)
    run_info["elapsed_seconds"] = elapsed

    _DEMO_CACHE = {
        "run_info": run_info,
        "primer_pairs": combined_pairs,
        "elapsed_seconds": elapsed,
    }

    # Persist for next startup
    demo_path.parent.mkdir(parents=True, exist_ok=True)
    demo_path.write_text(json.dumps(_DEMO_CACHE, indent=2))

    return JSONResponse(_DEMO_CACHE)


@router.get("/demo/genome")
async def get_demo_genome():
    """Return the example genome string so the live panel can re-run without pasting.

    The genome is small (~3 kbp example), so returning it inline is fine.
    """
    genome_fasta = _DATA_DIR / "example_genome.fa"
    if not genome_fasta.exists():
        raise HTTPException(404, "Example genome not found")
    lines = genome_fasta.read_text().splitlines()
    genome_string = "".join(l for l in lines if not l.startswith(">"))
    return JSONResponse({"genome_string": genome_string, "length": len(genome_string)})


@router.post("/params")
async def live_rerun(req: LivePCRRequest):
    """Live parameter experiment: re-run a single primer pair with new settings.

    Returns a structured {success, results, summary} envelope so the frontend
    can update plots without a page reload.  Never exposes raw tracebacks.
    """
    import logging, time
    log = logging.getLogger(__name__)
    log.info("live_rerun: received request for pair=%r", req.primer_name)

    try:
        run_pcr = _import_api()
    except ImportError as e:
        log.error("live_rerun: pipeline import failed: %s", e)
        return JSONResponse({"success": False, "error": "Pipeline unavailable — check server logs."}, status_code=500)

    # Basic primer validation (catches obvious front-end errors before the pipeline sees them)
    _VALID_BASES = set("ACGTacgtNn")
    for label, seq in [("forward primer", req.fwd_primer), ("reverse primer", req.rev_primer)]:
        seq = seq.strip()
        if not seq:
            return JSONResponse({"success": False, "error": f"Empty {label}."}, status_code=422)
        if not all(b in _VALID_BASES for b in seq):
            bad = next(b for b in seq if b not in _VALID_BASES)
            return JSONResponse({"success": False, "error": f"Invalid base '{bad}' in {label}."}, status_code=422)
        if len(seq) < 10:
            return JSONResponse({"success": False, "error": f"{label.capitalize()} is too short (min 10 nt)."}, status_code=422)

    t0 = time.perf_counter()
    log.info("live_rerun: starting pipeline run for %r", req.primer_name)

    try:
        result = run_pcr(
            fwd_primer     = req.fwd_primer.strip().upper(),
            rev_primer     = req.rev_primer.strip().upper(),
            primer_name    = req.primer_name,
            genome_string  = req.genome_string,
            genome_id      = "live",
            max_mismatches = req.mismatches,
            min_amplicon   = req.min_size,
            max_amplicon   = req.max_size,
            na_conc_mm     = req.na_conc,
            mg_conc_mm     = req.mg_conc,
            dntp_conc_mm   = req.dntp_conc,
            primer_conc_nm = req.primer_conc,
            run_hairpin    = req.run_hairpin,
            run_dimer      = req.run_dimer,
            log_level      = "WARNING",
        )
    except Exception as e:
        log.error("live_rerun: pipeline error: %s", e, exc_info=True)
        # Return a clean message — never a traceback
        msg = str(e)
        if len(msg) > 200:
            msg = msg[:200] + "…"
        return JSONResponse({"success": False, "error": msg}, status_code=500)

    elapsed = round(time.perf_counter() - t0, 3)
    out = result.get("json_output", {})
    pairs = out.get("primer_pairs", [])
    amps  = pairs[0].get("amplicons", []) if pairs else []
    best  = next((a for a in amps if a.get("is_intended")), amps[0] if amps else None)

    log.info("live_rerun: done in %.3fs — %d amplicon(s)", elapsed, len(amps))

    return JSONResponse({
        "success":  True,
        "results":  out,
        "elapsed":  elapsed,
        "summary": {
            "amplicons_found": len(amps),
            "offtargets":      out.get("primer_pairs", [{}])[0].get("offtarget_summary", {}).get("total", 0) if pairs else 0,
            "top_score":       round(best["final_score"], 2) if best else None,
            "top_length":      best["length"] if best else None,
        },
    })


@router.get("/benchmark")
async def get_benchmark():
    """Return any cached benchmark metrics from the benchmarking module."""
    bench_dir = _DATA_DIR.parent / "benchmarking"
    report_file = bench_dir / "benchmark_report.json"
    if report_file.exists():
        return JSONResponse(json.loads(report_file.read_text()))
    # Return a placeholder structure if no real benchmarks run
    return JSONResponse({
        "status": "no_benchmark_run",
        "message": "Run benchmarking/runner.py to generate metrics.",
        "metrics": {}
    })

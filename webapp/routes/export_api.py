"""
Export routes: JSON, CSV, and HTML report download.
"""
import io
import csv
import json
from fastapi import APIRouter
from fastapi.responses import StreamingResponse, JSONResponse
from pydantic import BaseModel

router = APIRouter(prefix="/api/export", tags=["export"])


class ExportRequest(BaseModel):
    result: dict
    format: str = "json"   # json | csv | html


@router.post("/json")
async def export_json(req: ExportRequest):
    content = json.dumps(req.result, indent=2)
    return StreamingResponse(
        io.BytesIO(content.encode()),
        media_type="application/json",
        headers={"Content-Disposition": 'attachment; filename="pcr_results.json"'},
    )


@router.post("/csv")
async def export_csv(req: ExportRequest):
    buf = io.StringIO()
    w = csv.writer(buf)
    w.writerow([
        "pair_name", "rank", "seq_id", "start", "end", "length",
        "gc_fraction", "final_score", "fwd_tm", "rev_tm",
        "fwd_mismatches", "rev_mismatches", "is_intended",
    ])
    for pair in req.result.get("primer_pairs", []):
        name = pair.get("name", "")
        for amp in pair.get("amplicons", []):
            fb = amp.get("fwd_binding", {})
            rb = amp.get("rev_binding", {})
            w.writerow([
                name,
                amp.get("rank", ""),
                amp.get("seq_id", ""),
                amp.get("start", ""),
                amp.get("end", ""),
                amp.get("length", ""),
                amp.get("gc_fraction", ""),
                amp.get("final_score", ""),
                fb.get("tm_celsius", ""),
                rb.get("tm_celsius", ""),
                fb.get("mismatch_count", ""),
                rb.get("mismatch_count", ""),
                amp.get("is_intended", ""),
            ])
    buf.seek(0)
    return StreamingResponse(
        io.BytesIO(buf.getvalue().encode()),
        media_type="text/csv",
        headers={"Content-Disposition": 'attachment; filename="pcr_results.csv"'},
    )

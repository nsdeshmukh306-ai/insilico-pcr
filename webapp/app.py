"""
In-Silico PCR Dashboard — FastAPI application entry point.
"""

import sys
import os
from pathlib import Path

# Make insilico_pcr importable when running directly from webapp/
_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent   # /home/workshop/niraj
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from fastapi.middleware.cors import CORSMiddleware

from webapp.routes import pcr_api, export_api

app = FastAPI(
    title="In-Silico PCR Dashboard",
    description="Interactive analysis dashboard for the publication-grade in-silico PCR pipeline",
    version="1.0.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

_STATIC = Path(__file__).parent / "static"
_TEMPLATES = Path(__file__).parent / "templates"

app.mount("/static", StaticFiles(directory=str(_STATIC)), name="static")
templates = Jinja2Templates(directory=str(_TEMPLATES))

app.include_router(pcr_api.router)
app.include_router(export_api.router)


# Root — serve the dashboard SPA
from fastapi import Request
from fastapi.responses import HTMLResponse


@app.get("/", response_class=HTMLResponse)
async def index(request: Request):
    return templates.TemplateResponse(request, "index.html")

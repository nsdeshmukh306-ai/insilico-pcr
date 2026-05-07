"""
Convenience launcher: python3 webapp/run.py
"""
import sys
from pathlib import Path

# Ensure the project root (parent of insilico_pcr/) is on sys.path
ROOT = Path(__file__).resolve().parent.parent.parent   # /home/workshop/niraj
sys.path.insert(0, str(ROOT))

import uvicorn

if __name__ == "__main__":
    uvicorn.run(
        "insilico_pcr.webapp.app:app",
        host="0.0.0.0",
        port=8765,
        reload=True,
        log_level="info",
    )

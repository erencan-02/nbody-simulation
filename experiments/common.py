"""Shared helpers for the evaluation experiments.

Every experiment writes its figures to ``thesis/figures/generated`` and its
numerical results to ``thesis/data`` so that the thesis is always built from
the latest results (tables are read directly by Typst).
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from nbody.integrators import create_integrator  # noqa: E402
from nbody.viz.style import apply_style, integrator_style  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "thesis" / "figures" / "generated"
DATA_DIR = ROOT / "thesis" / "data"


def setup() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    # Note: the thesis fonts are CFF-flavoured OpenType files, which matplotlib
    # cannot convert to SVG paths reliably; figures use matplotlib's built-in
    # DejaVu Sans (text) and Computer Modern (math) instead.
    apply_style()


def save_fig(fig, name: str, fmt: str = "svg") -> Path:
    path = FIG_DIR / f"{name}.{fmt}"
    fig.savefig(path, dpi=200 if fmt == "png" else None)
    plt.close(fig)
    return path


def write_csv(name: str, header: list[str], rows: list[list]) -> Path:
    path = DATA_DIR / f"{name}.csv"
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(header)
        for row in rows:
            w.writerow([_fmt(v) for v in row])
    return path


def sanitize(obj):
    """Make results JSON-compliant: numpy scalars -> Python, nan/inf -> None."""
    if isinstance(obj, dict):
        return {str(k): sanitize(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [sanitize(v) for v in obj]
    if hasattr(obj, "item") and not isinstance(obj, (str, bytes)):
        obj = obj.item()
    if isinstance(obj, float) and not math.isfinite(obj):
        return None
    return obj


def write_json(name: str, obj) -> Path:
    path = DATA_DIR / f"{name}.json"
    path.write_text(json.dumps(sanitize(obj), indent=2), encoding="utf-8")
    return path


def _fmt(v):
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, float):
        return f"{v:.6g}"
    return v


def label(name: str) -> str:
    return create_integrator(name).label


def style(name: str, **kw) -> dict:
    s = integrator_style(name)
    s.update(kw)
    return s


def log(msg: str) -> None:
    print(f"    {msg}", flush=True)

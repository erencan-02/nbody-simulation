"""Run all evaluation experiments and regenerate the thesis figures and tables.

Usage (from the repository root)::

    python -m experiments.run_all            # everything (~8 min)
    python -m experiments.run_all e02 e05    # a subset
"""

from __future__ import annotations

import argparse
import importlib
import json
import time

from .common import DATA_DIR, sanitize

EXPERIMENTS = {
    "e01": ("e01_harmonic", "Harmonic oscillator: phase space and linear stability"),
    "e02": ("e02_convergence", "Convergence order and work-precision (Kepler)"),
    "e03": ("e03_long_term", "Long-term behaviour at equal cost (Kepler)"),
    "e04": ("e04_legacy", "Revisiting the original course experiment"),
    "e05": ("e05_figure_eight", "Figure-eight choreography vs. DOP853 reference"),
    "e06": ("e06_pythagorean", "Pythagorean problem: adaptivity and chaos"),
    "e07": ("e07_solar_system", "Solar System over 1000 years"),
    "e08": ("e08_performance", "Force-kernel performance"),
    "e09": ("e09_viewer", "Viewer screenshots with trajectory overlays"),
}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("only", nargs="*", help="experiment ids, e.g. e02 e05")
    args = parser.parse_args(argv)
    selected = args.only or list(EXPERIMENTS)
    unknown = [e for e in selected if e not in EXPERIMENTS]
    if unknown:
        parser.error(f"unknown experiments {unknown}; choose from {list(EXPERIMENTS)}")

    summary_path = DATA_DIR / "summary.json"
    summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}
    t_total = time.perf_counter()
    for key in selected:
        module_name, title = EXPERIMENTS[key]
        print(f"[{key}] {title}", flush=True)
        t0 = time.perf_counter()
        module = importlib.import_module(f"experiments.{module_name}")
        summary[key] = module.run()
        print(f"    done in {time.perf_counter() - t0:.1f} s", flush=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(sanitize(summary), indent=2), encoding="utf-8")
    print(f"all done in {time.perf_counter() - t_total:.1f} s; results in {DATA_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""E8 - Performance of the force evaluation.

Wall-clock time of one force evaluation versus N for
* the pair loop of the original implementation (``acceleration_loop``),
* the dense broadcast kernel and
* the pair-list kernel (used automatically for small N).
Timings are machine dependent; the ratios and the O(N^2) scaling are not.
"""

from __future__ import annotations

import platform
import time

import matplotlib.pyplot as plt
import numpy as np

from nbody.core import GravityModel, acceleration_loop

from .common import log, save_fig, setup, write_csv, write_json


def _time(fn, min_time=0.2, max_reps=20000):
    fn()  # warm-up
    reps, total = 0, 0.0
    t0 = time.perf_counter()
    while total < min_time and reps < max_reps:
        fn()
        reps += 1
        total = time.perf_counter() - t0
    return total / reps


def run() -> dict:
    setup()
    rng = np.random.default_rng(0)
    ns = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]
    rows = []
    res = {"loop": {}, "dense": {}, "pairs": {}}
    for n in ns:
        x = rng.normal(size=(n, 2))
        m = rng.uniform(0.5, 1.5, n)
        dense = GravityModel(m, kernel="dense")
        t_dense = _time(lambda: dense.acceleration(x))
        t_pairs = float("nan")
        if n <= 128:
            pairs = GravityModel(m, kernel="pairs")
            t_pairs = _time(lambda: pairs.acceleration(x))
        t_loop = _time(lambda: acceleration_loop(x, m), max_reps=200) if n <= 512 else float("nan")
        res["loop"][n], res["dense"][n], res["pairs"][n] = t_loop, t_dense, t_pairs
        rows.append([n, t_loop, t_dense, t_pairs, t_loop / min(t_dense, np.nan_to_num(t_pairs, nan=np.inf))])
        log(f"N={n:5d}: loop {t_loop * 1e6:10.1f} us, dense {t_dense * 1e6:9.1f} us, "
            f"pairs {t_pairs * 1e6:8.1f} us")
    write_csv("performance", ["N", "t_loop_s", "t_dense_s", "t_pairs_s", "speedup_vs_loop"], rows)

    fig, ax = plt.subplots(figsize=(5.8, 2.6))
    for key, lab, col, mk in [("loop", "pair loop (original)", "#E6001A", "o"),
                              ("dense", "dense broadcast kernel", "#243572", "s"),
                              ("pairs", "pair-list kernel", "#0083CC", "^")]:
        nn = np.array([n for n in ns if np.isfinite(res[key][n])])
        tt = np.array([res[key][n] for n in nn])
        ax.loglog(nn, tt * 1e6, marker=mk, color=col, lw=1.0, ms=3.5, label=lab)
    nref = np.array([64, 2048])
    t0 = res["dense"][64] * 1e6
    ax.loglog(nref, t0 * (nref / 64) ** 2, color="black", lw=0.6, ls="--", label=r"$\propto N^2$")
    ax.set_xlabel("number of bodies $N$")
    ax.set_ylabel(r"time per force evaluation [$\mu$s]")
    ax.legend(fontsize=6.5)
    fig.tight_layout()
    save_fig(fig, "performance")
    out = {
        "machine": f"{platform.processor() or platform.machine()}, Python {platform.python_version()}, "
        f"NumPy {np.__version__}",
        "speedup_n512": res["loop"][512] / res["dense"][512],
        "speedup_n2": res["loop"][2] / res["pairs"][2],
        "speedup_n128": res["loop"][128] / min(res["dense"][128], res["pairs"][128]),
        "t_dense_n2048_ms": res["dense"][2048] * 1e3,
    }
    write_json("performance", out)
    return out


if __name__ == "__main__":
    run()

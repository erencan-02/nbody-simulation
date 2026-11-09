"""E2 - Convergence order and work-precision on an eccentric Kepler orbit.

The maximum global phase-space error over two orbits (sampled ~256 times)
is compared with the exact solution for a sequence of step sizes.  Taking the
maximum over the interval avoids spurious superconvergence at special
instants (e.g. the O(h) errors of symplectic Euler cancel exactly at every
periapsis passage).  The empirical order is the slope of
a least-squares fit in log-log space restricted to the asymptotic regime.
The cost of each run is measured in force evaluations.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from nbody import Simulation, load_scenario, simulate
from nbody.analysis import estimate_order, phase_space_error
from nbody.analysis.kepler import kepler_orbit, orbital_period
from nbody.integrators import DormandPrince54, available_integrators, create_integrator

from .common import label, log, save_fig, setup, style, write_csv

E, Q = 0.6, 1e-3
# Fit window (min, max error) defining the asymptotic regime.  Explicit Euler
# only becomes asymptotic for errors below O(1) on this eccentric orbit.
FIT_WINDOW = (1e-11, 1e-2)
FIT_WINDOW_EULER = (1e-11, 1.0)
COST_TARGET = 1e-6


def run() -> dict:
    setup()
    sc = load_scenario("kepler", e=E, mass_ratio=Q)
    P = orbital_period(1.0, 1.0, Q)
    T = 2.0 * P

    def error(traj):
        x_ex, v_ex = kepler_orbit(traj.t, 1.0, Q, 1.0, E)
        return float(np.max(phase_space_error(traj.x, traj.v, x_ex, v_ex)))

    fixed = [n for n in available_integrators() if not create_integrator(n).adaptive]
    rows, orders = [], {}
    results = {}
    for name in fixed:
        integ = create_integrator(name)
        max_exp = {1: 17, 2: 14}.get(integ.order, 12)
        ns = 2 ** np.arange(4, max_exp + 1)  # steps per orbit
        hs, errs, evals = [], [], []
        for n in ns:
            model = sc.model()
            traj = simulate(model, sc.initial_state(), integ, T, P / n,
                            record_every=max(1, int(2 * n) // 256))
            hs.append(P / n)
            errs.append(error(traj))
            evals.append(model.n_evals)
            rows.append([name, int(n), P / n, model.n_evals, errs[-1]])
        fit = estimate_order(hs, errs, *(FIT_WINDOW_EULER if name == "euler" else FIT_WINDOW))
        orders[name] = fit
        results[name] = (np.array(hs), np.array(errs), np.array(evals))
        log(f"{name:<17s} measured order {fit.order:5.2f} (theory {integ.order}, {fit.n_points} pts)")

    # adaptive Dormand-Prince for a range of tolerances
    dp = {"rtol": [], "err": [], "evals": [], "steps": []}
    for rtol in np.logspace(-3, -13, 11):
        model = sc.model()
        sim = Simulation(model, sc.initial_state(), DormandPrince54(rtol=rtol, atol=rtol * 1e-3))
        traj = sim.run(T, 1e-2, record_every=1)
        dp["rtol"].append(rtol)
        dp["err"].append(error(traj))
        dp["evals"].append(model.n_evals)
        dp["steps"].append(traj.meta["n_steps"])
        rows.append(["dopri5", traj.meta["n_steps"], float("nan"), model.n_evals, dp["err"][-1]])
    write_csv("convergence", ["integrator", "steps_per_orbit", "dt", "force_evals", "error"], rows)

    # ---- figure ------------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(5.8, 3.0), sharey=True)
    for name, (hs, errs, evals) in results.items():
        ax1.loglog(hs, errs, ms=2.5, lw=0.9, label=label(name), **style(name), marker=_marker(name))
        ax2.loglog(evals, errs, ms=2.5, lw=0.9, **style(name), marker=_marker(name))
    ax2.loglog(dp["evals"], dp["err"], ms=3, lw=0.9, label=label("dopri5"), **style("dopri5"), marker="x")
    # reference slopes
    for p, (h0, e0) in {1: (2e-2, 3e-1), 2: (2e-2, 1.5e-3), 4: (8e-2, 1e-6), 6: (1.2e-1, 1e-9)}.items():
        h = np.array([h0 / 4, h0])
        ax1.loglog(h, e0 * (h / h0) ** p, color="black", lw=0.6)
        ax1.text(h[0] * 0.75, e0 * (0.25) ** p, f"$h^{p}$", fontsize=6.5, ha="right", va="center")
    ax1.set_xlabel("step size $h$")
    ax1.set_ylabel(r"max. phase-space error, $t\leq 2P$")
    ax1.set_title("(a) convergence")
    ax2.set_xlabel("force evaluations")
    ax2.set_title("(b) work-precision")
    ax1.set_ylim(1e-14, 1e2)
    handles, labels = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    fig.legend(handles + h2, labels + l2, loc="lower center", ncol=5, fontsize=6.5,
               bbox_to_anchor=(0.5, -0.06))
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    save_fig(fig, "convergence")

    table = []
    for name in fixed:
        integ = create_integrator(name)
        table.append([name, integ.label, integ.order, round(orders[name].order, 2),
                      orders[name].n_points, integ.evals_per_step,
                      integ.symplectic, integ.time_reversible])
    write_csv("orders", ["integrator", "label", "order_theory", "order_measured", "fit_points",
                         "evals_per_step", "symplectic", "time_reversible"], table)
    # cost to reach a given phase-space error (interpolated in log-log space)
    cost = {}
    for name, (hs, errs, evals) in list(results.items()) + [
        ("dopri5", (None, np.array(dp["err"]), np.array(dp["evals"])))
    ]:
        cost[name] = _cost_for_error(errs, evals, COST_TARGET)
    write_csv("cost_to_accuracy", ["integrator", "label", "force_evals"],
              [[n, label(n), c] for n, c in cost.items()])
    return {"orders": {k: v.order for k, v in orders.items()}, "cost_target": COST_TARGET, "cost": cost}


def _marker(name):
    from nbody.viz.style import INTEGRATOR_STYLE
    return INTEGRATOR_STYLE[name]["marker"]


def _cost_for_error(errs, evals, target):
    """Log-log interpolation of the evaluations needed for ``err = target``."""
    errs, evals = np.asarray(errs), np.asarray(evals, dtype=float)
    order = np.argsort(evals)
    errs, evals = errs[order], evals[order]
    for k in range(len(errs) - 1):
        if errs[k] >= target >= errs[k + 1]:
            f = (np.log(target) - np.log(errs[k])) / (np.log(errs[k + 1]) - np.log(errs[k]))
            return float(np.exp(np.log(evals[k]) + f * (np.log(evals[k + 1]) - np.log(evals[k]))))
    return float("nan")


if __name__ == "__main__":
    run()

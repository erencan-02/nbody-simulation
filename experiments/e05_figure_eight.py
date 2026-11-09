"""E5 - Figure-eight choreography: accuracy against an independent reference.

The reference solution is computed with SciPy's DOP853 (an independent,
eighth-order implementation) at rtol = atol = 1e-13.  All fixed-step
methods use the same budget of force evaluations per period.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

from nbody import Simulation, load_scenario
from nbody.analysis import conservation_report, phase_space_error
from nbody.core import acceleration
from nbody.integrators import DormandPrince54, create_integrator

from .common import label, log, save_fig, setup, style, write_csv

PERIOD = 6.32591398
N_PERIODS = 50
EVALS_PER_PERIOD = 1260  # divisible by 1, 2, 3, 4, 7, 9
METHODS = ["symplectic_euler", "leapfrog", "rk4", "yoshida4", "yoshida6a"]


def reference(sc, t_end, rtol=1e-13):
    """Dense-output DOP853 solution; returns a function t -> (x, v)."""
    m, G = sc.masses, sc.G
    n, d = sc.positions.shape

    def rhs(_t, y):
        x = y[: n * d].reshape(n, d)
        v = y[n * d :]
        return np.concatenate([v, acceleration(x, m, G).ravel()])

    y0 = np.concatenate([sc.positions.ravel(), sc.velocities.ravel()])
    sol = solve_ivp(rhs, (0.0, t_end), y0, method="DOP853", dense_output=True,
                    rtol=rtol, atol=rtol)
    assert sol.success

    def evaluate(t):
        y = sol.sol(np.asarray(t)).T
        return y[:, : n * d].reshape(-1, n, d), y[:, n * d :].reshape(-1, n, d)

    return evaluate


def run() -> dict:
    setup()
    sc = load_scenario("figure_eight")
    T = N_PERIODS * PERIOD
    ref = reference(sc, T)
    # self-check of the reference: compare with a run at a 10x looser tolerance
    t_chk = np.linspace(0, T, 500)
    ref_err = float(np.max(phase_space_error(*ref(t_chk), *reference(sc, T, 1e-12)(t_chk))))
    log(f"reference self-check (rtol 1e-13 vs 1e-12): {ref_err:.1e}")
    runs = {}
    for name in METHODS + ["dopri5"]:
        if name == "dopri5":
            integ, dt, rec = DormandPrince54(rtol=1e-9, atol=1e-12), PERIOD / 100, 5
        else:
            integ = create_integrator(name)
            dt = PERIOD * integ.evals_per_step / EVALS_PER_PERIOD
            rec = max(1, int(round(T / dt)) // 3000)
        traj = Simulation(sc.model(), sc.initial_state(), integ).run(T, dt, record_every=rec)
        x_ref, v_ref = ref(traj.t)
        err = phase_space_error(traj.x, traj.v, x_ref, v_ref)
        runs[name] = (traj, err)
        log(f"{name:<17s} final error {err[-1]:.2e}, evals/period "
            f"{traj.meta['n_force_evals'] / N_PERIODS:.0f}")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(5.8, 2.5), gridspec_kw={"width_ratios": [1.1, 1]})
    traj = runs["yoshida6a"][0]
    one = traj.t <= PERIOD
    for i in range(3):
        ax1.plot(traj.x[one, i, 0], traj.x[one, i, 1], color=sc.colors[i], lw=1.2,
                 alpha=0.85, label=sc.names[i])
        ax1.plot(*traj.x[0, i], "o", color=sc.colors[i], ms=4)
    ax1.set_aspect("equal")
    ax1.set_xlabel("$x$")
    ax1.set_ylabel("$y$")
    ax1.set_title("(a) figure-eight orbit (one period)")
    for name, (tr, err) in runs.items():
        ax2.loglog(tr.t[1:] / PERIOD, err[1:], lw=0.8, **style(name), label=label(name))
    ax2.set_xlabel("time [periods]")
    ax2.set_ylabel("phase-space error vs. DOP853")
    ax2.set_title(f"(b) error, {EVALS_PER_PERIOD} evaluations / period")
    h, l = ax2.get_legend_handles_labels()
    ax1.legend(h, l, fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.32), ncol=2)
    ax2.set_ylim(1e-12, 1e1)
    fig.tight_layout()
    save_fig(fig, "figure_eight")

    rows = []
    for name, (tr, err) in runs.items():
        rep = conservation_report(tr)
        rows.append([name, label(name), tr.meta["n_force_evals"], rep.energy_max,
                     rep.angular_momentum_max, float(err[-1])])
    write_csv("figure_eight", ["integrator", "label", "force_evals", "energy_max",
                               "angmom_max", "error_final"], rows)
    out = {r[0]: r[-1] for r in rows}
    out["reference_check"] = ref_err
    return out


if __name__ == "__main__":
    run()

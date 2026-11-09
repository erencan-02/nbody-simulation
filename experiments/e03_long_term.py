"""E3 - Long-term behaviour at equal computational cost.

All methods integrate an eccentric Kepler orbit for 200 periods with the
same budget of force evaluations per orbit, so differences reflect the
methods and not the amount of work.  Recorded: energy error, angular
momentum error and the global error with respect to the exact solution.
A separate forward-backward run measures time reversibility.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from nbody import Simulation, load_scenario
from nbody.analysis import conservation_report, phase_space_error, reversibility_error
from nbody.analysis.kepler import kepler_orbit, orbital_period
from nbody.integrators import DormandPrince54, create_integrator

from .common import label, log, save_fig, setup, style, write_csv

E, Q = 0.5, 1e-3
N_ORBITS = 200
EVALS_PER_ORBIT = 1008  # divisible by 1, 2, 3, 4, 7 and 9 evaluations per step
METHODS = ["euler", "symplectic_euler", "heun", "leapfrog", "rk4", "yoshida4", "yoshida6a"]


def run() -> dict:
    setup()
    sc = load_scenario("kepler", e=E, mass_ratio=Q)
    P = orbital_period(1.0, 1.0, Q)
    T = N_ORBITS * P
    runs = {}
    for name in METHODS + ["dopri5"]:
        if name == "dopri5":
            integ = DormandPrince54(rtol=1e-10, atol=1e-13)  # similar cost
            dt, rec = P / 100, 8
        else:
            integ = create_integrator(name)
            dt = P * integ.evals_per_step / EVALS_PER_ORBIT
            rec = max(1, int(round(T / dt)) // 5000)
        sim = Simulation(sc.model(), sc.initial_state(), integ)
        with np.errstate(all="ignore"):
            traj = sim.run(T, dt, record_every=rec)
        x_ex, v_ex = kepler_orbit(traj.t, 1.0, Q, 1.0, E)
        err = phase_space_error(traj.x, traj.v, x_ex, v_ex)
        runs[name] = (traj, err)
        log(f"{name:<17s} {traj.meta['n_force_evals'] / N_ORBITS:7.1f} evals/orbit, "
            f"final error {err[-1]:.2e}, {traj.meta['wall_time']:.1f} s")

    # ---- figure ------------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(5.8, 3.0))
    for name, (traj, err) in runs.items():
        tp = traj.t[1:] / P
        dE = np.abs(traj.energy() - traj.energy()[0]) / abs(traj.energy()[0])
        ax1.loglog(tp, np.maximum(dE[1:], 1e-16), lw=0.7, **style(name), label=label(name))
        ax2.loglog(tp, err[1:], lw=0.7, **style(name))
    t_ref = np.array([3.0, 300.0])
    ax2.loglog(t_ref, 2e-4 * t_ref, color="black", lw=0.6)
    ax2.loglog(t_ref, 2e-7 * t_ref**2, color="black", lw=0.6)
    ax2.text(300, 2e-4 * 300 * 1.6, r"$\propto t$", fontsize=7, ha="center")
    ax2.text(300, 2e-7 * 9e4 * 0.3, r"$\propto t^2$", fontsize=7, ha="left")
    ax1.set_ylim(1e-12, 1e1)
    ax2.set_ylim(1e-7, 1e2)
    ax1.set_xlabel("time [orbits]")
    ax2.set_xlabel("time [orbits]")
    ax1.set_ylabel(r"$|\Delta E/E_0|$")
    ax2.set_ylabel("phase-space error")
    ax1.set_title(f"(a) energy error, {EVALS_PER_ORBIT} evaluations / orbit")
    ax2.set_title("(b) global error vs. exact solution")
    h, l = ax1.get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=4, fontsize=6.5, bbox_to_anchor=(0.5, -0.05))
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    save_fig(fig, "long_term_kepler")

    # ---- table -----------------------------------------------------------------------
    rows = []
    out = {}
    for name, (traj, err) in runs.items():
        with np.errstate(all="ignore"):
            rep = conservation_report(traj)
        integ = create_integrator(name)
        rev = (
            float("nan") if integ.adaptive
            else reversibility_error(sc.model(), sc.initial_state(), integ, P / 200, 1000)
        )
        # growth exponent of the global error over the second half of the run
        half = traj.t > T / 10
        good = half & np.isfinite(err) & (err > 0) & (err < 1.0)
        slope = (
            float(np.polyfit(np.log(traj.t[good]), np.log(err[good]), 1)[0])
            if good.sum() > 10 and err[-1] < 2.0 else float("nan")  # nan: saturated
        )
        rows.append([name, label(name), traj.meta["n_force_evals"], rep.energy_max,
                     rep.energy_drift_rate * P, rep.angular_momentum_max, float(err[-1]),
                     slope, rev])
        out[name] = dict(energy_max=rep.energy_max, error_final=float(err[-1]),
                         error_growth=slope, reversibility=rev,
                         angular_momentum_max=rep.angular_momentum_max)
    write_csv(
        "long_term_kepler",
        ["integrator", "label", "force_evals", "energy_max", "energy_drift_per_orbit",
         "angmom_max", "error_final", "error_growth_exponent", "reversibility_error"],
        rows,
    )
    return out


if __name__ == "__main__":
    run()

"""E7 - The Solar System over 1000 years.

Leapfrog (h = 1 day), Yoshida-4 and RK4 with the same force-evaluation
budget.  Physical validation: simulated mean orbital periods versus Kepler's
third law with the two-body correction, P = 2 pi sqrt(a^3 / (G (M_sun + m))).
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from nbody import Simulation, load_scenario
from nbody.analysis import conservation_report, mean_motion_period, relative_energy_error
from nbody.analysis.kepler import orbital_period
from nbody.integrators import create_integrator
from nbody.scenarios import planetary_elements

from .common import label, log, save_fig, setup, style, write_csv

DAY = 1.0 / 365.25
T_END = 1000.0
METHODS = ["leapfrog", "yoshida4", "rk4"]


def run() -> dict:
    setup()
    sc = load_scenario("solar_system")
    runs = {}
    for name in METHODS:
        integ = create_integrator(name)
        dt = DAY * integ.evals_per_step  # equal cost: one evaluation per day
        rec = max(1, int(round(5 * DAY / dt)))  # record roughly every 5 days
        tr = Simulation(sc.model(), sc.initial_state(), integ).run(T_END, dt, record_every=rec)
        runs[name] = tr
        rep = conservation_report(tr)
        log(f"{name:<9s} dt={dt * 365.25:.0f} d, max|dE/E0|={rep.energy_max:.2e}, "
            f"{tr.meta['wall_time']:.0f} s")

    # ---- figure: orbits and energy -------------------------------------------------
    tr = runs["leapfrog"]
    fig, axs = plt.subplots(1, 3, figsize=(5.8, 2.3), gridspec_kw={"width_ratios": [1, 1, 1.35]})
    inner = tr.t <= 4.0
    outer = tr.t <= 170.0
    for i in range(1, sc.n_bodies):
        rel = tr.x[:, i] - tr.x[:, 0]
        if i <= 4:
            axs[0].plot(rel[inner, 0], rel[inner, 1], lw=0.7, color=sc.colors[i], label=sc.names[i])
        axs[1].plot(rel[outer, 0], rel[outer, 1], lw=0.5, color=sc.colors[i], label=sc.names[i])
    for ax, lim, title in [(axs[0], 1.8, "(a) inner planets"), (axs[1], 33, "(b) all planets")]:
        ax.plot(0, 0, "o", color=sc.colors[0], ms=3)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_title(title)
        ax.set_xlabel("$x$ [AU]")
        ax.tick_params(labelsize=6)
    axs[0].set_ylabel("$y$ [AU]")
    axs[0].legend(fontsize=5, loc="upper right", handlelength=1)
    ax = axs[2]
    for name, tr_m in runs.items():
        dE = np.abs(relative_energy_error(tr_m.energy()))
        # envelope: maximum of |dE| in 400 time windows (the raw signal oscillates
        # with Mercury's period and would fill the panel)
        edges = np.linspace(0, dE.size, 401).astype(int)
        t_env = np.array([tr_m.t[a:b].mean() for a, b in zip(edges[:-1], edges[1:])])
        e_env = np.array([dE[a:b].max() for a, b in zip(edges[:-1], edges[1:])])
        ax.semilogy(t_env, np.maximum(e_env, 1e-16), lw=1.0, **style(name),
                    label=f"{label(name)}, $h$={create_integrator(name).evals_per_step} d")
    ax.set_xlabel("time [yr]")
    ax.set_ylabel(r"$|\Delta E/E_0|$")
    ax.set_title(r"(c) energy error (envelope)")
    ax.set_ylim(1e-9, 1e0)
    ax.legend(fontsize=5.5, loc="lower right", frameon=True, framealpha=0.9, edgecolor="none")
    fig.tight_layout()
    save_fig(fig, "solar_system")

    # ---- periods ---------------------------------------------------------------------
    # Reference for the terrestrial planets: Yoshida 6A with h = 1 day over 100 yr
    # (Mercury: ~88 steps/orbit at order 6).  The giant planets are resolved by
    # >4000 steps per orbit in the 1000-yr leapfrog run, which serves as their
    # reference (their periods need the long time span).
    T_REF = 100.0
    ref = Simulation(sc.model(), sc.initial_state(), "yoshida6a").run(T_REF, DAY, record_every=2)
    elements = planetary_elements()
    rows, phase_rows = [], []
    for k, pname in enumerate(sc.names[1:], start=1):
        a = elements[pname][0]
        P_kep = orbital_period(a, 1.0, sc.masses[k], G=sc.G)
        terrestrial = a < 2.0
        src = ref if terrestrial else runs["leapfrog"]
        P_sim = mean_motion_period(src.t, src.x[:, k] - src.x[:, 0])
        rows.append([pname, a, P_kep, P_sim, (P_sim - P_kep) / P_kep,
                     "Yoshida 6A, 1 d, 100 yr" if terrestrial else "leapfrog, 1 d, 1000 yr"])
        if terrestrial:
            errs = []
            for name in METHODS:
                tr_m = runs[name]
                w = tr_m.t <= T_REF
                P_m = mean_motion_period(tr_m.t[w], tr_m.x[w, k] - tr_m.x[w, 0])
                errs.append((P_m - P_sim) / P_sim)
            phase_rows.append([pname, P_sim] + errs)
    write_csv("solar_periods", ["planet", "a_AU", "P_kepler_yr", "P_simulated_yr",
                                "rel_deviation", "source"], rows)
    write_csv("solar_phase_error", ["planet", "P_reference_yr"]
              + [f"rel_period_error_{n}" for n in METHODS], phase_rows)
    log("period error vs reference (Mercury): " + ", ".join(
        f"{n}={e:+.1e}" for n, e in zip(METHODS, phase_rows[0][2:])))

    table = []
    for name, tr_m in runs.items():
        rep = conservation_report(tr_m)
        table.append([name, label(name), create_integrator(name).evals_per_step,
                      tr_m.meta["n_force_evals"], rep.energy_max, rep.energy_drift_rate,
                      rep.angular_momentum_max, tr_m.meta["wall_time"]])
    write_csv("solar_conservation", ["integrator", "label", "dt_days", "force_evals", "energy_max",
                                     "energy_drift_per_yr", "angmom_max", "wall_time_s"], table)
    return {r[0]: r[-1] for r in rows}


if __name__ == "__main__":
    run()

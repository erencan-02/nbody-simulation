"""E4 - Revisiting the experiment of the original course project.

The original report evaluated the energy quotient q(t) = |E_t / E_0| of a
two-body system (G = 20, softening 5, dt = 0.1, 5000 leapfrog steps).  Its
force used Plummer softening but its energy diagnostic used the *unsoftened*
potential -G m_i m_j / r.  This experiment quantifies how much of the
reported energy fluctuation was an artefact of that inconsistency.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from nbody import load_scenario, simulate
from nbody.core import potential_energy

from .common import ROOT, label, log, save_fig, setup, style, write_json


def run() -> dict:
    setup()
    sc = load_scenario(str(ROOT / "scenarios" / "two_body_legacy.json"))
    n_steps, dt = 5000, 0.1
    traj = simulate(sc.model(), sc.initial_state(), "leapfrog", n_steps * dt, dt)
    K = traj.kinetic_energy()
    E_inconsistent = K + potential_energy(traj.x, sc.masses, sc.G, 0.0)
    E_consistent = K + traj.potential_energy()
    q_inc = np.abs(E_inconsistent / E_inconsistent[0])
    q_con = np.abs(E_consistent / E_consistent[0])
    r = np.linalg.norm(traj.x[:, 1] - traj.x[:, 0], axis=-1)
    steps = np.arange(traj.n_frames)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(5.8, 2.6))
    ax1.plot(steps, q_inc, color="#E6001A", lw=0.9, label="unsoftened $V$ (original)")
    ax1.plot(steps, q_con, color="#243572", lw=0.9, label="softened $V$ (consistent)")
    ax1.set_xlabel("time step")
    ax1.set_ylabel(r"$q(t)=|E_t/E_0|$")
    ax1.set_title("(a) energy quotient, leapfrog $h=0.1$")
    ax1.legend(loc="lower center", fontsize=6.5, frameon=True, framealpha=0.95, edgecolor="none")
    axr = ax1.twinx()
    axr.plot(steps, r, color="#B5B5B5", lw=0.5, zorder=0)
    axr.set_ylabel("separation $r$", color="#898989")
    axr.grid(False)
    axr.spines["right"].set_visible(True)

    res = {}
    for name, h in [("leapfrog", 0.1), ("leapfrog", 0.05), ("yoshida4", 0.1), ("rk4", 0.1)]:
        tr = simulate(sc.model(), sc.initial_state(), name, n_steps * 0.1, h,
                      record_every=max(1, int(round(0.1 / h))))
        E = tr.energy()
        dE = np.abs(E - E[0]) / abs(E[0])
        ax2.semilogy(tr.t, np.maximum(dE, 1e-16), lw=0.8,
                     **style(name, ls="-" if h == 0.1 else ":"), label=f"{label(name)}, $h={h}$")
        res[f"{name}_h{h}"] = float(dE.max())
    ax2.set_ylim(1e-13, 1e-2)
    ax2.set_xlabel("time $t$")
    ax2.set_ylabel(r"$|\Delta E/E_0|$ (consistent)")
    ax2.set_title("(b) consistent diagnostic")
    ax2.legend(fontsize=6, loc="lower right", ncol=1, frameon=True, framealpha=0.95, edgecolor="none")
    fig.tight_layout()
    save_fig(fig, "legacy_two_body")

    out = {
        "max_dev_inconsistent": float(np.max(np.abs(q_inc - 1))),
        "max_dev_consistent": float(np.max(np.abs(q_con - 1))),
        "artefact_factor": float(np.max(np.abs(q_inc - 1)) / np.max(np.abs(q_con - 1))),
        "r_min": float(r.min()),
        "r_max": float(r.max()),
        "softening": sc.softening,
        "max_energy_error": res,
    }
    write_json("legacy_two_body", out)
    log(f"max |q-1|: inconsistent {out['max_dev_inconsistent']:.2e}, consistent "
        f"{out['max_dev_consistent']:.2e} (factor {out['artefact_factor']:.0f})")
    return out


if __name__ == "__main__":
    run()

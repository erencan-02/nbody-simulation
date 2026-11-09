"""E1 - Linear test problem: harmonic oscillator.

(a) phase portraits and long-time energy error at a large step size;
(b) spectral radius of the one-step map M(h*omega) (linear stability).
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from nbody import HarmonicModel, State, simulate
from nbody.integrators import create_integrator

from .common import label, log, save_fig, setup, style, write_csv

METHODS = ["euler", "symplectic_euler", "leapfrog", "rk4", "yoshida4"]


def one_step_matrix(name: str, h: float) -> np.ndarray:
    """The integrators are linear on x'' = -x: build the 2x2 matrix of one step."""
    integ = create_integrator(name)
    M = np.zeros((2, 2))
    for col, (x0, v0) in enumerate([(1.0, 0.0), (0.0, 1.0)]):
        st = State([[x0]], [[v0]])
        integ.step(st, h, HarmonicModel([1.0]))
        M[:, col] = [st.x[0, 0], st.v[0, 0]]
    return M


def run() -> dict:
    setup()
    h = 0.2
    model = HarmonicModel([1.0])
    state = State([[1.0]], [[0.0]])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(5.8, 2.7))
    phi = np.linspace(0, 2 * np.pi, 400)
    ax1.plot(np.cos(phi), np.sin(phi), color="black", lw=0.6, ls=":", label="exact")
    for name in METHODS[:4]:
        traj = simulate(model.copy(), state, name, 2 * 2 * np.pi, h)
        ax1.plot(traj.x[:, 0, 0], traj.v[:, 0, 0], lw=0.9, **style(name), label=label(name))
    ax1.set_xlabel("$x$")
    ax1.set_ylabel("$v$")
    ax1.set_aspect("equal", adjustable="datalim")
    ax1.set_title(r"(a) phase space, $h\omega=0.2$, two periods")
    ax1.legend(loc="lower left", fontsize=6, frameon=True, framealpha=0.92, edgecolor="none")

    rows = []
    for name in METHODS:
        traj = simulate(model.copy(), state, name, 1000.0, h, record_every=5)
        E = traj.energy()
        dE = np.abs(E - E[0]) / abs(E[0])
        ax2.loglog(traj.t[1:], np.maximum(dE[1:], 1e-17), lw=0.9, **style(name), label=label(name))
        rows.append([name, label(name), float(dE.max()), float(dE[-1])])
    ax2.set_ylim(1e-9, 1e3)
    ax2.set_xlabel(r"time $t\,\omega$")
    ax2.set_ylabel(r"$|E(t)-E_0|/E_0$")
    ax2.set_title("(b) energy error")
    fig.tight_layout()
    save_fig(fig, "harmonic_phase_energy")
    write_csv("harmonic_energy", ["integrator", "label", "max_rel_energy_error", "final_rel_energy_error"], rows)

    # stability: spectral radius of M(h) for h*omega in (0, 3.2)
    hs = np.linspace(1e-3, 3.2, 1200)
    fig, ax = plt.subplots(figsize=(5.8, 2.3))
    limits = {}
    for name in METHODS:
        rho = np.array([max(abs(np.linalg.eigvals(one_step_matrix(name, hh)))) for hh in hs])
        ax.plot(hs, rho, lw=1.1, **style(name), label=label(name))
        unstable = hs[rho > 1 + 1e-9]
        limits[name] = float(unstable[0]) if unstable.size else float("inf")
    ax.axhline(1.0, color="black", lw=0.5)
    ax.set_ylim(0.0, 2.0)
    ax.set_xlim(0, 3.2)
    ax.set_xlabel(r"$h\omega$")
    ax.set_ylabel(r"spectral radius $\rho(M(h\omega))$")
    ax.legend(ncol=3, loc="upper left", fontsize=6.5)
    fig.tight_layout()
    save_fig(fig, "harmonic_stability")
    write_csv(
        "harmonic_stability",
        ["integrator", "label", "stability_limit"],
        [[n, label(n), limits[n]] for n in METHODS],
    )
    log("stability limits h*omega: " + ", ".join(f"{n}={v:.3f}" for n, v in limits.items()))
    return {"stability_limits": limits}


if __name__ == "__main__":
    run()

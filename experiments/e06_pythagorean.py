"""E6 - Pythagorean three-body problem: close encounters and chaos.

* Adaptive Dormand-Prince (rtol 1e-12) versus fixed-step leapfrog and
  Yoshida-4 with the *same* number of force evaluations.
* Sensitivity to initial conditions: a second adaptive run with the position
  of one body perturbed by 1e-8; the growth of the phase-space separation
  gives a finite-time Lyapunov exponent.
* Tolerance study: final state and qualitative outcome (which pair forms the
  final binary) as a function of the tolerance, compared with an independent
  SciPy DOP853 reference at rtol = atol = 3e-14.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from nbody import Simulation, load_scenario
from nbody.analysis import finite_time_lyapunov, phase_space_error, relative_energy_error
from nbody.integrators import DormandPrince54, create_integrator
from nbody.simulation import Trajectory

from .common import label, log, save_fig, setup, style, write_csv, write_json

T_END = 70.0
SAMPLE = 0.2  # sampling interval (also an upper bound for the adaptive step)
RTOL = 1e-12
DELTA0 = 1e-8
TOLERANCES = [1e-8, 1e-9, 1e-10, 1e-11]  # (+ RTOL from the nominal run)
PAIRS = ["m=3 + m=4", "m=3 + m=5", "m=4 + m=5"]


def binary_pair(x: np.ndarray) -> tuple[str, float]:
    """Closest pair of the final configuration (the binary) and its separation."""
    d = [np.linalg.norm(x[i] - x[j]) for i, j in [(0, 1), (0, 2), (1, 2)]]
    k = int(np.argmin(d))
    return PAIRS[k], float(d[k])


def scipy_reference(sc, t_end: float) -> tuple[np.ndarray, np.ndarray]:
    from scipy.integrate import solve_ivp

    from nbody.core import acceleration

    n, d = sc.positions.shape

    def rhs(_t, y):
        return np.concatenate([y[n * d:], acceleration(y[: n * d].reshape(n, d), sc.masses).ravel()])

    y0 = np.concatenate([sc.positions.ravel(), sc.velocities.ravel()])
    sol = solve_ivp(rhs, (0.0, t_end), y0, method="DOP853", rtol=3e-14, atol=3e-14)
    assert sol.success
    y = sol.y[:, -1]
    return y[: n * d].reshape(n, d), y[n * d:].reshape(n, d)


def run_segments(sim: Simulation, grid: np.ndarray, record_all: bool):
    """Advance an adaptive simulation through ``grid``; return samples (+ full path)."""
    xs, vs, full, steps = [sim.state.x.copy()], [sim.state.v.copy()], [], []
    for t in grid[1:]:
        tr = sim.run(t, 1e-3, record_every=1 if record_all else 10**9)
        xs.append(tr.x[-1])
        vs.append(tr.v[-1])
        if record_all:
            full.append(tr if not full else _drop_first(tr))
            steps.append(tr.step_sizes)
    return np.array(xs), np.array(vs), full, (np.concatenate(steps) if steps else None)


def _drop_first(tr: Trajectory) -> Trajectory:
    return Trajectory(tr.t[1:], tr.x[1:], tr.v[1:], tr.masses, tr.G, tr.softening)


def run() -> dict:
    setup()
    sc = load_scenario("pythagorean")
    grid = np.round(np.arange(0.0, T_END + SAMPLE / 2, SAMPLE), 10)

    # nominal adaptive run (full resolution + samples on the grid)
    model = sc.model()
    sim = Simulation(model, sc.initial_state(), DormandPrince54(rtol=RTOL, atol=RTOL * 1e-3))
    xs, vs, parts, steps = run_segments(sim, grid, record_all=True)
    t_full = np.concatenate([p.t for p in parts])
    x_full = np.concatenate([p.x for p in parts])
    budget = model.n_evals
    log(f"dopri5: {steps.size} steps, {budget} force evaluations, "
        f"h in [{steps.min():.1e}, {steps.max():.1e}]")

    # perturbed run
    st = sc.initial_state()
    st.x[0, 0] += DELTA0
    sim_p = Simulation(sc.model(), st, DormandPrince54(rtol=RTOL, atol=RTOL * 1e-3))
    xp, vp, _, _ = run_segments(sim_p, grid, record_all=False)
    sep = phase_space_error(xs, vs, xp, vp)  # separation in phase space
    sat = np.argmax(sep > 1e-2) if np.any(sep > 1e-2) else len(sep)
    lam = finite_time_lyapunov(grid[1:sat], sep[1:sat])
    log(f"finite-time Lyapunov exponent ~ {lam:.2f}, separation reaches 1e-2 at t={grid[sat - 1]:.1f}")

    # energy along the sampled adaptive solution
    m = sc.masses
    from nbody.analysis.diagnostics import total_energy

    dE_dp = relative_energy_error(total_energy(xs, vs, m, sc.G))

    # fixed-step methods with the same evaluation budget
    fixed = {}
    for name in ["leapfrog", "yoshida4"]:
        integ = create_integrator(name)
        dt = T_END * integ.evals_per_step / budget
        with np.errstate(all="ignore"):
            tr = Simulation(sc.model(), sc.initial_state(), integ).run(
                T_END, dt, record_every=max(1, int(SAMPLE / dt)))
            fixed[name] = (tr, relative_energy_error(tr.energy()))
        log(f"{name}: dt={dt:.2e}, final |dE/E0| = {abs(fixed[name][1][-1]):.2e}")

    # minimum pairwise distance along the full-resolution path
    d01 = np.linalg.norm(x_full[:, 0] - x_full[:, 1], axis=-1)
    d02 = np.linalg.norm(x_full[:, 0] - x_full[:, 2], axis=-1)
    d12 = np.linalg.norm(x_full[:, 1] - x_full[:, 2], axis=-1)
    dmin = np.minimum(np.minimum(d01, d02), d12)

    # ---- figure -------------------------------------------------------------------
    fig, axs = plt.subplots(2, 2, figsize=(5.8, 4.9))
    ax = axs[0, 0]
    for i in range(3):
        ax.plot(x_full[:, i, 0], x_full[:, i, 1], lw=0.4, color=sc.colors[i], label=sc.names[i])
        ax.plot(*sc.positions[i], "o", ms=3.5, color=sc.colors[i])
    ax.set_xlim(-6, 6)
    ax.set_ylim(-6, 6)
    ax.set_aspect("equal")
    ax.set_title(r"(a) trajectories, $0\leq t\leq 70$")
    ax.legend(fontsize=6, loc="lower left")

    ax = axs[0, 1]
    t_steps = t_full[1:]
    ax.semilogy(t_steps, steps, lw=0.6, color="#535353", label="step size $h$")
    ax.semilogy(t_full, dmin, lw=0.6, color="#E6001A", label=r"min. distance $r_{\min}$")
    ax.set_xlabel("time $t$")
    ax.set_title("(b) adaptive step size")
    ax.legend(fontsize=6, loc="lower right")

    ax = axs[1, 0]
    ax.semilogy(grid, np.maximum(np.abs(dE_dp), 1e-17), lw=0.7, **style("dopri5"),
                label=f"{label('dopri5')} (adaptive)")
    for name, (tr, dE) in fixed.items():
        ax.semilogy(tr.t, np.maximum(np.abs(dE), 1e-17), lw=0.7, **style(name),
                    label=f"{label(name)} (fixed $h$)")
    ax.set_ylim(1e-14, 1e2)
    ax.set_xlabel("time $t$")
    ax.set_ylabel(r"$|\Delta E/E_0|$")
    ax.set_title("(c) energy error at equal cost")
    ax.legend(fontsize=6, loc="lower right")

    ax = axs[1, 1]
    ax.semilogy(grid[1:], sep[1:], lw=0.8, color="#243572", label=r"$\delta(t)$")
    ax.semilogy(grid[1:sat], sep[1] * np.exp(lam * (grid[1:sat] - grid[1])), lw=0.6,
                ls="--", color="black", label=rf"$\propto e^{{{lam:.2f}\,t}}$")
    ax.set_xlabel("time $t$")
    ax.set_title(r"(d) separation, $\delta(0)=10^{-8}$")
    ax.legend(fontsize=6, loc="lower right")
    fig.tight_layout()
    save_fig(fig, "pythagorean")

    # ---- tolerance study: is the qualitative outcome converged? -------------------
    x_ref, v_ref = scipy_reference(sc, T_END)
    ref_pair, ref_sep = binary_pair(x_ref)
    log(f"reference (DOP853, 3e-14): binary {ref_pair}")
    tol_rows = []
    for rtol in TOLERANCES + [RTOL]:
        if rtol == RTOL:
            xf, vf, evals = xs[-1], vs[-1], budget
        else:
            mdl = sc.model()
            tr = Simulation(mdl, sc.initial_state(), DormandPrince54(rtol=rtol, atol=rtol * 1e-3)).run(
                T_END, 1e-3, record_every=10**9)
            xf, vf, evals = tr.x[-1], tr.v[-1], mdl.n_evals
        pair, sep_b = binary_pair(xf)
        err = float(phase_space_error(xf, vf, x_ref, v_ref))
        tol_rows.append([rtol, evals, pair, sep_b, err, pair == ref_pair])
        log(f"rtol={rtol:.0e}: {evals} evals, binary {pair}, final error vs reference {err:.2e}")
    write_csv("pythagorean_tolerance",
              ["rtol", "force_evals", "binary", "binary_separation", "final_error_vs_reference",
               "outcome_correct"], tol_rows)

    out = {
        "rtol": RTOL,
        "delta0": DELTA0,
        "reference_binary": ref_pair,
        "dopri_steps": int(steps.size),
        "dopri_evals": int(budget),
        "dopri_rejected": int(sim.n_rejected),
        "h_min": float(steps.min()),
        "h_max": float(steps.max()),
        "r_min": float(dmin.min()),
        "t_r_min": float(t_full[np.argmin(dmin)]),
        "dopri_energy_max": float(np.max(np.abs(dE_dp))),
        "lyapunov": float(lam),
        "t_separation_1e-2": float(grid[sat - 1]),
        "fixed": {
            name: {"dt": float(T_END * create_integrator(name).evals_per_step / budget),
                   "energy_max": float(np.nanmax(np.abs(dE)))}
            for name, (tr, dE) in fixed.items()
        },
    }
    write_json("pythagorean", out)
    return out


if __name__ == "__main__":
    run()

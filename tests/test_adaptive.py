import numpy as np
import pytest

from nbody import GravityModel, Simulation, State
from nbody.analysis import position_error
from nbody.analysis.kepler import kepler_initial_conditions, kepler_orbit
from nbody.integrators import DormandPrince54

E = 0.9  # very eccentric: step sizes must vary strongly along the orbit
P = 2 * np.pi


def run(rtol):
    x0, v0 = kepler_initial_conditions(1.0, 1e-3, a=1.0, e=E)
    model = GravityModel([1.0, 1e-3])
    sim = Simulation(model, State(x0, v0), DormandPrince54(rtol=rtol, atol=rtol * 1e-3))
    traj = sim.run(P * np.sqrt(1 / 1.001), 1e-2, record_every=1)
    return traj


def test_adaptive_run_hits_final_time_exactly():
    traj = run(1e-8)
    assert traj.t[-1] == pytest.approx(P * np.sqrt(1 / 1.001), abs=0, rel=0)
    assert traj.step_sizes is not None and traj.step_sizes.size == traj.meta["n_steps"]
    assert np.isclose(traj.step_sizes.sum(), traj.t[-1] - traj.t[0], rtol=1e-12)


def test_step_size_adapts_to_periapsis():
    traj = run(1e-8)
    h = traj.step_sizes
    assert h.max() / h.min() > 20  # small steps near periapsis, large at apoapsis


def test_error_decreases_with_tolerance():
    errs = []
    for rtol in (1e-6, 1e-8, 1e-10):
        traj = run(rtol)
        x_ex, _ = kepler_orbit(traj.t[-1:], 1.0, 1e-3, 1.0, E)
        errs.append(float(position_error(traj.x[-1], x_ex[0])))
    assert errs[0] > errs[1] > errs[2]
    assert errs[2] < 1e-6


def test_fsal_evaluation_count():
    traj = run(1e-8)
    n_acc, n_rej = traj.meta["n_steps"], traj.meta["n_rejected"]
    assert traj.meta["n_force_evals"] == 1 + 6 * (n_acc + n_rej)

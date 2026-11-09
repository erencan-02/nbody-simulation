import numpy as np
import pytest

from nbody import Simulation, load_scenario
from nbody.analysis import (
    estimate_order,
    finite_time_lyapunov,
    mean_motion_period,
    relative_energy_error,
)
from nbody.analysis.metrics import linear_drift_rate
from nbody.viz import TrailBuffer, TrajectoryPredictor


def test_relative_energy_error():
    np.testing.assert_allclose(relative_energy_error([-2.0, -2.2, -1.8]), [0.0, -0.1, 0.1])


def test_estimate_order_on_synthetic_data_ignores_plateau():
    h = np.logspace(-4, -1, 10)
    err = 3.0 * h**4
    err = np.maximum(err, 1e-13)  # artificial round-off plateau
    fit = estimate_order(h, err, err_min=1e-12)
    assert fit.order == pytest.approx(4.0, abs=1e-10)
    assert fit.n_points < h.size


def test_linear_drift_and_lyapunov_fit():
    t = np.linspace(0, 10, 101)
    assert linear_drift_rate(t, 2.5 * t + 1.0) == pytest.approx(2.5)
    assert finite_time_lyapunov(t, 1e-9 * np.exp(0.7 * t)) == pytest.approx(0.7)


def test_mean_motion_period_of_circular_orbit():
    t = np.linspace(0, 30, 3001)
    rel = np.stack([np.cos(2 * np.pi * t / 2.5), np.sin(2 * np.pi * t / 2.5)], axis=-1)
    assert mean_motion_period(t, rel) == pytest.approx(2.5, rel=1e-12)


def test_trail_buffer_ring_order():
    buf = TrailBuffer(n_bodies=1, dim=2, capacity=4)
    for k in range(6):
        buf.push(np.array([[k, -k]], dtype=float))
    out = buf.ordered()
    assert len(buf) == 4
    np.testing.assert_array_equal(out[:, 0, 0], [2, 3, 4, 5])
    buf.shift(np.array([1.0, 0.0]))
    np.testing.assert_array_equal(buf.ordered()[:, 0, 0], [3, 4, 5, 6])


def test_predictor_does_not_modify_running_simulation():
    sc = load_scenario("circumbinary")
    sim = Simulation(sc.model(), sc.initial_state(), "leapfrog")
    for _ in range(10):
        sim.step(0.01)
    before = sim.state.copy()
    evals = sim.model.n_evals
    pred = TrajectoryPredictor(horizon=1.0, dt=0.01, n_points=50).predict(sim)
    np.testing.assert_array_equal(sim.state.x, before.x)
    assert sim.state.t == before.t and sim.model.n_evals == evals
    # the forecast must coincide with actually integrating forward
    ref = sim.clone().run(sim.state.t + 1.0, 0.01, record_every=2)
    np.testing.assert_allclose(pred[-1], ref.x[-1], rtol=1e-12)
    np.testing.assert_allclose(pred[0], before.x)

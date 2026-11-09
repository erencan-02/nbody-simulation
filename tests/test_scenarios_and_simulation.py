from pathlib import Path

import numpy as np
import pytest

from nbody import Scenario, Simulation, Trajectory, available_scenarios, load_scenario, simulate
from nbody.analysis import conservation_report, diagnostics as diag, mean_motion_period
from nbody.analysis.kepler import orbital_period
from nbody.scenarios import planetary_elements

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("name", available_scenarios())
def test_builtin_scenarios_are_valid(name):
    sc = load_scenario(name)
    assert sc.n_bodies == len(sc.names)
    assert np.all(np.isfinite(sc.positions)) and np.all(sc.masses >= 0)
    sim = Simulation(sc.model(), sc.initial_state(), "leapfrog")
    sim.step(sc.recommended["dt"])


@pytest.mark.parametrize("name", ["solar_system", "circumbinary", "random_cluster", "kepler"])
def test_scenarios_are_in_the_barycentric_frame(name):
    sc = load_scenario(name)
    np.testing.assert_allclose(diag.linear_momentum(sc.velocities, sc.masses), 0, atol=1e-14)
    np.testing.assert_allclose(diag.center_of_mass(sc.positions, sc.masses), 0, atol=1e-14)


def test_random_cluster_is_virialised_and_reproducible():
    a, b = load_scenario("random_cluster", n=30, seed=7), load_scenario("random_cluster", n=30, seed=7)
    np.testing.assert_array_equal(a.positions, b.positions)
    Q = diag.virial_ratio(a.positions, a.velocities, a.masses, a.G, a.softening)
    assert Q == pytest.approx(0.5, rel=1e-12)


def test_json_roundtrip(tmp_path):
    sc = load_scenario("figure_eight")
    path = tmp_path / "fig8.json"
    sc.save_json(path)
    back = load_scenario(str(path))
    np.testing.assert_array_equal(back.positions, sc.positions)
    assert back.names == sc.names and back.recommended == sc.recommended


@pytest.mark.parametrize("fname", ["two_body_legacy.json", "three_body_legacy.json"])
def test_shipped_json_scenarios_load(fname):
    sc = Scenario.load_json(ROOT / "scenarios" / fname)
    assert sc.G == 20.0 and sc.softening == 5.0


def test_original_course_format_is_supported():
    legacy = {"G": 20.0, "bodies": [{"mass": 1000, "position": [0, 0], "velocity": [0, 0]},
                                    {"mass": 10, "position": [50, 0], "velocity": [0, -25]}]}
    sc = Scenario.from_dict(legacy)
    assert sc.n_bodies == 2 and sc.softening == 0.0 and sc.colors is None


def test_fixed_step_run_records_and_ends_exactly():
    sc = load_scenario("kepler")
    traj = simulate(sc.model(), sc.initial_state(), "leapfrog", 1.0, 0.03, record_every=4)
    n = int(np.ceil(1.0 / 0.03))
    assert traj.meta["n_steps"] == n
    assert traj.t[-1] == 1.0
    assert traj.n_frames == n // 4 + 1 + (1 if n % 4 else 0)
    np.testing.assert_allclose(np.diff(traj.t)[:-1], 4 * (1.0 / n))


def test_simulation_continues_from_current_time():
    sc = load_scenario("kepler")
    sim = Simulation(sc.model(), sc.initial_state(), "yoshida4")
    sim.run(1.0, 0.01)
    second = sim.run(2.0, 0.01)
    full = simulate(sc.model(), sc.initial_state(), "yoshida4", 2.0, 0.01)
    assert second.t[0] == 1.0
    np.testing.assert_allclose(second.x[-1], full.x[-1], rtol=1e-12, atol=1e-14)


def test_trajectory_save_load(tmp_path):
    sc = load_scenario("pythagorean")
    traj = simulate(sc.model(), sc.initial_state(), "dopri5", 1.0, 1e-3)
    path = tmp_path / "run.npz"
    traj.save(path)
    back = Trajectory.load(path)
    np.testing.assert_array_equal(back.x, traj.x)
    np.testing.assert_array_equal(back.step_sizes, traj.step_sizes)
    assert back.meta == traj.meta and back.G == traj.G


def test_solar_system_periods_follow_keplers_third_law():
    sc = load_scenario("solar_system", planets=("Mercury", "Venus", "Earth"))
    traj = simulate(sc.model(), sc.initial_state(), "yoshida4", 10.0, 1 / 365.25, record_every=2)
    for k, name in enumerate(sc.names[1:], start=1):
        a = planetary_elements()[name][0]
        P_sim = mean_motion_period(traj.t, traj.x[:, k] - traj.x[:, 0])
        P_kep = orbital_period(a, 1.0, sc.masses[k], G=sc.G)
        assert P_sim == pytest.approx(P_kep, rel=2e-3), name
    # 4th-order method, ~88 steps per Mercury orbit: bounded error of O(1e-6)
    assert conservation_report(traj).energy_max < 5e-6


def test_trajectory_energy_uses_the_generating_model():
    from nbody import HarmonicModel, State

    traj = simulate(HarmonicModel([2.0], omega=3.0), State([[1.0]], [[0.0]]), "leapfrog", 1.0, 0.01)
    E = traj.energy()
    assert E[0] == pytest.approx(0.5 * 2.0 * 9.0)  # (1/2) m omega^2 x^2
    assert np.max(np.abs(E - E[0])) / E[0] < 1e-3

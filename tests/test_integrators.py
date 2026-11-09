"""Check every integrator against the theory: order, symplecticity, symmetry, cost."""

import numpy as np
import pytest

from nbody import GravityModel, HarmonicModel, State, simulate
from nbody.analysis import conservation_report, estimate_order, phase_space_error, reversibility_error
from nbody.analysis.kepler import kepler_initial_conditions, kepler_orbit
from nbody.integrators import available_integrators, create_integrator, triple_jump_weights

FIXED = [n for n in available_integrators() if not create_integrator(n).adaptive]
SYMPLECTIC = [n for n in FIXED if create_integrator(n).symplectic]
SYMMETRIC = [n for n in FIXED if create_integrator(n).time_reversible]

M1, M2, E = 1.0, 0.5, 0.3
PERIOD = 2 * np.pi * np.sqrt(1.0 / (M1 + M2))
T_END = 1.37 * PERIOD  # deliberately not a multiple of the period


def kepler_setup():
    x0, v0 = kepler_initial_conditions(M1, M2, a=1.0, e=E)
    return GravityModel([M1, M2]), State(x0, v0)


def global_error(name, n_steps):
    model, state = kepler_setup()
    traj = simulate(model, state, name, T_END, T_END / n_steps, record_every=n_steps)
    x_exact, v_exact = kepler_orbit(np.array([T_END]), M1, M2, 1.0, E)
    return float(phase_space_error(traj.x[-1], traj.v[-1], x_exact[0], v_exact[0]))


@pytest.mark.parametrize("name", FIXED)
def test_convergence_order(name):
    integ = create_integrator(name)
    p = integ.order
    # choose resolutions in the asymptotic regime but above round-off
    base = {1: 2000, 2: 200, 4: 100, 6: 60}[p]
    ns = np.array([base, 2 * base, 4 * base])
    errs = np.array([global_error(name, n) for n in ns])
    fit = estimate_order(T_END / ns, errs)
    assert abs(fit.order - p) < 0.35, f"{name}: measured order {fit.order:.2f}, expected {p}"


@pytest.mark.parametrize("name", FIXED)
def test_force_evaluation_count(name):
    integ = create_integrator(name)
    model, state = kepler_setup()
    n = 17
    simulate(model, state, integ, 1.0, 1.0 / n, record_every=n)
    # FSAL methods need one extra evaluation for the very first step
    fsal_extra = 1 if name in SYMPLECTIC and integ.time_reversible else 0
    assert model.n_evals == n * integ.evals_per_step + fsal_extra


@pytest.mark.parametrize("name", SYMMETRIC)
def test_symmetric_methods_are_time_reversible(name):
    model, state = kepler_setup()
    err = reversibility_error(model, state, create_integrator(name), 0.05, 200)
    assert err < 1e-11


@pytest.mark.parametrize("name", ["euler", "symplectic_euler", "heun", "rk4"])
def test_non_symmetric_methods_are_not_reversible(name):
    model, state = kepler_setup()
    assert reversibility_error(model, state, create_integrator(name), 0.05, 200) > 1e-8


@pytest.mark.parametrize("name", SYMPLECTIC)
def test_symplectic_methods_have_bounded_energy_error(name):
    """Energy error of a symplectic method oscillates without secular drift."""
    model, state = kepler_setup()
    dt = T_END / 150
    traj = simulate(model, state, name, 100 * T_END, dt, record_every=5)
    dE = np.abs(traj.energy() - traj.energy()[0])
    first, last = dE[: dE.size // 10].max(), dE[-dE.size // 10 :].max()
    assert last < 3.0 * first + 1e-13


def test_explicit_euler_energy_grows_secularly():
    model, state = kepler_setup()
    traj = simulate(model, state, "euler", 20 * T_END, T_END / 2000, record_every=10)
    rep = conservation_report(traj)
    assert rep.energy_drift_rate > 0 and rep.energy_final > 1e-3


@pytest.mark.parametrize("name", FIXED)
def test_linear_momentum_is_conserved_by_all_methods(name):
    rng = np.random.default_rng(0)
    model = GravityModel(rng.uniform(0.5, 1.5, 5), softening=0.1)
    state = State(rng.normal(size=(5, 2)), rng.normal(size=(5, 2)) * 0.3)
    traj = simulate(model, state, name, 2.0, 0.01, record_every=50)
    assert conservation_report(traj).momentum_max < 1e-12


@pytest.mark.parametrize("name", SYMPLECTIC)
def test_symplectic_methods_conserve_angular_momentum(name):
    model, state = kepler_setup()
    traj = simulate(model, state, name, 5 * T_END, T_END / 100, record_every=10)
    assert conservation_report(traj).angular_momentum_max < 1e-12


def test_triple_jump_weights_conditions():
    for p in (2, 4, 6):
        g1, g2, g3 = triple_jump_weights(p)
        assert g1 == g3
        assert abs(g1 + g2 + g3 - 1.0) < 1e-14
        assert abs(g1 ** (p + 1) + g2 ** (p + 1) + g3 ** (p + 1)) < 1e-13


def test_symplectic_euler_is_area_preserving_on_harmonic_oscillator():
    """The one-step map of a symplectic method has a Jacobian with det = 1."""
    h = 0.3
    for name in SYMPLECTIC + ["euler", "rk4"]:
        integ = create_integrator(name)
        jac = np.zeros((2, 2))
        for col, (x0, v0) in enumerate([(1.0, 0.0), (0.0, 1.0)]):
            st = State([[x0]], [[v0]])
            integ.step(st, h, HarmonicModel([1.0]))
            jac[:, col] = [st.x[0, 0], st.v[0, 0]]  # the map is linear
        det = np.linalg.det(jac)
        if integ.symplectic:
            assert abs(det - 1.0) < 1e-13, name
        else:
            assert abs(det - 1.0) > 1e-7, name  # RK4: det = 1 - h^6/72 + O(h^8)


def test_registry_aliases_and_errors():
    assert create_integrator("verlet").name == "leapfrog"
    assert create_integrator("RK45").name == "dopri5"
    with pytest.raises(KeyError):
        create_integrator("does-not-exist")

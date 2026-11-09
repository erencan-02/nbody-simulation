import numpy as np
import pytest

from nbody.analysis import diagnostics as diag
from nbody.analysis.kepler import kepler_orbit, orbital_period, solve_kepler


@pytest.mark.parametrize("e", [0.0, 0.3, 0.9, 0.99])
def test_kepler_equation_is_solved(e):
    M = np.linspace(-20, 20, 1001)
    E = solve_kepler(M, e)
    np.testing.assert_allclose(E - e * np.sin(E), M, atol=1e-12)


@pytest.mark.parametrize("e", [0.0, 0.5, 0.9])
def test_exact_orbit_conserves_invariants_and_is_periodic(e):
    m1, m2 = 1.0, 0.3
    P = orbital_period(1.0, m1, m2)
    t = np.linspace(0, 2 * P, 777)
    x, v = kepler_orbit(t, m1, m2, a=1.0, e=e)
    m = np.array([m1, m2])
    E = diag.total_energy(x, v, m)
    L = diag.angular_momentum(x, v, m)
    np.testing.assert_allclose(E, -m1 * m2 / 2.0, rtol=1e-12)  # E = -G m1 m2 / (2a)
    np.testing.assert_allclose(L, L[0], rtol=1e-12)
    np.testing.assert_allclose(diag.linear_momentum(v, m), 0.0, atol=1e-14)
    np.testing.assert_allclose(x[-1], x[0], atol=1e-10)


def test_exact_orbit_satisfies_equations_of_motion():
    m1, m2, e = 1.0, 0.2, 0.6
    t = np.array([0.37])
    h = 1e-4
    xm, _ = kepler_orbit(t - h, m1, m2, e=e)
    x0, _ = kepler_orbit(t, m1, m2, e=e)
    xp, _ = kepler_orbit(t + h, m1, m2, e=e)
    acc_fd = (xp - 2 * x0 + xm)[0] / h**2
    from nbody.core import acceleration

    np.testing.assert_allclose(acc_fd, acceleration(x0[0], np.array([m1, m2])), rtol=1e-6)

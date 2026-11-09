import numpy as np
import pytest

from nbody.core import GravityModel, acceleration, acceleration_loop, potential_energy
from nbody.core import forces


def random_config(n=7, d=2, seed=0):
    rng = np.random.default_rng(seed)
    return rng.normal(size=(n, d)), rng.uniform(0.1, 2.0, size=n)


def test_two_body_acceleration_matches_newton():
    x = np.array([[0.0, 0.0], [2.0, 0.0]])
    m = np.array([3.0, 5.0])
    a = acceleration(x, m, G=1.5)
    # |a_1| = G m_2 / r^2 towards body 2, |a_2| = G m_1 / r^2 towards body 1
    np.testing.assert_allclose(a[0], [1.5 * 5.0 / 4.0, 0.0])
    np.testing.assert_allclose(a[1], [-1.5 * 3.0 / 4.0, 0.0])


@pytest.mark.parametrize("d", [2, 3])
@pytest.mark.parametrize("eps", [0.0, 0.1])
def test_vectorised_kernel_equals_pair_loop(d, eps):
    x, m = random_config(n=9, d=d)
    np.testing.assert_allclose(
        acceleration(x, m, 2.0, eps), acceleration_loop(x, m, 2.0, eps), rtol=1e-12, atol=1e-14
    )


@pytest.mark.parametrize("n", [2, 3, 10, 30])
def test_pair_kernel_equals_dense_kernel(n):
    x, m = random_config(n=n, d=3, seed=n)
    dense = GravityModel(m, G=0.7, softening=0.02, kernel="dense")
    pairs = GravityModel(m, G=0.7, softening=0.02, kernel="pairs")
    np.testing.assert_allclose(pairs.acceleration(x), dense.acceleration(x), rtol=1e-12, atol=1e-15)
    assert GravityModel(m).copy().kernel == "auto"


def test_blocked_kernel_equals_unblocked(monkeypatch):
    x, m = random_config(n=40, d=3, seed=4)
    ref = acceleration(x, m, 1.0, 0.05)
    monkeypatch.setattr(forces, "_MAX_PAIR_ELEMENTS", 3 * 40 * 7)  # blocks of 7 rows
    np.testing.assert_allclose(acceleration(x, m, 1.0, 0.05), ref, rtol=1e-13)


def test_newtons_third_law_total_force_vanishes():
    x, m = random_config(n=25, d=3, seed=2)
    a = acceleration(x, m, 1.0, 0.01)
    total = (m[:, None] * a).sum(axis=0)
    scale = np.abs(m[:, None] * a).sum()
    assert np.all(np.abs(total) < 1e-14 * scale)


@pytest.mark.parametrize("eps", [0.0, 0.3])
def test_force_is_negative_gradient_of_potential(eps):
    """Softened force and softened potential must be consistent: m_i a_i = -dV/dx_i."""
    x, m = random_config(n=5, d=2, seed=3)
    a = acceleration(x, m, 1.0, eps)
    h = 1e-6
    grad = np.zeros_like(x)
    for i in range(x.shape[0]):
        for k in range(x.shape[1]):
            xp, xm = x.copy(), x.copy()
            xp[i, k] += h
            xm[i, k] -= h
            grad[i, k] = (potential_energy(xp, m, 1.0, eps) - potential_energy(xm, m, 1.0, eps)) / (2 * h)
    np.testing.assert_allclose(m[:, None] * a, -grad, rtol=1e-7, atol=1e-8)


def test_potential_energy_batched_and_chunked(monkeypatch):
    rng = np.random.default_rng(5)
    xs = rng.normal(size=(13, 6, 2))
    m = rng.uniform(1, 2, size=6)
    single = np.array([potential_energy(x, m, 1.0, 0.1) for x in xs])
    monkeypatch.setattr(forces, "_MAX_PAIR_ELEMENTS", 60)  # force several chunks
    np.testing.assert_allclose(potential_energy(xs, m, 1.0, 0.1), single, rtol=1e-14)


def test_model_counts_evaluations_and_validates():
    x, m = random_config(n=3)
    model = GravityModel(m)
    for _ in range(4):
        model.acceleration(x)
    assert model.n_evals == 4
    assert model.copy().n_evals == 0
    with pytest.raises(ValueError):
        GravityModel([1.0, -1.0])
    with pytest.raises(ValueError):
        GravityModel([1.0], softening=-0.1)

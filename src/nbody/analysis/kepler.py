r"""Exact solution of the two-body (Kepler) problem.

The relative coordinate :math:`\mathbf r = \mathbf r_2 - \mathbf r_1` obeys
:math:`\ddot{\mathbf r} = -\mu\,\mathbf r/|\mathbf r|^3` with
:math:`\mu = G(m_1+m_2)`.  For a bound orbit with semi-major axis :math:`a`
and eccentricity :math:`e` starting at periapsis, the position at time
:math:`t` follows from Kepler's equation

.. math::  M = n t = E - e\sin E, \qquad n = \sqrt{\mu/a^3},

as :math:`x = a(\cos E - e)`, :math:`y = a\sqrt{1-e^2}\sin E`.  The bodies
themselves move about the (fixed) centre of mass:
:math:`\mathbf r_1 = -\tfrac{m_2}{M}\mathbf r`, :math:`\mathbf r_2 = \tfrac{m_1}{M}\mathbf r`.
"""

from __future__ import annotations

import numpy as np


def solve_kepler(M: np.ndarray, e: float, tol: float = 1e-15, max_iter: int = 50) -> np.ndarray:
    """Solve ``E - e sin E = M`` for the eccentric anomaly with Newton's method."""
    if not 0.0 <= e < 1.0:
        raise ValueError("only elliptic orbits (0 <= e < 1) are supported")
    M = np.asarray(M, dtype=np.float64)
    # reduce to [-pi, pi) for a robust starting value, keep the full turns
    turns = np.floor((M + np.pi) / (2.0 * np.pi))
    Mr = M - 2.0 * np.pi * turns
    E = np.where(e > 0.8, np.pi * np.sign(Mr) + 0.0 * Mr, Mr + e * np.sin(Mr))
    for _ in range(max_iter):
        f = E - e * np.sin(E) - Mr
        dE = f / (1.0 - e * np.cos(E))
        E = E - dE
        if np.all(np.abs(dE) < tol):
            break
    return E + 2.0 * np.pi * turns


def orbital_period(a: float, m1: float, m2: float, G: float = 1.0) -> float:
    """Kepler's third law, :math:`P = 2\\pi\\sqrt{a^3/(G(m_1+m_2))}`."""
    return 2.0 * np.pi * np.sqrt(a**3 / (G * (m1 + m2)))


def kepler_initial_conditions(
    m1: float, m2: float, a: float = 1.0, e: float = 0.0, G: float = 1.0
) -> tuple[np.ndarray, np.ndarray]:
    """Positions and velocities ``(2, 2)`` at periapsis in the centre-of-mass frame."""
    x, v = kepler_orbit(np.array([0.0]), m1, m2, a, e, G)
    return x[0], v[0]


def kepler_orbit(
    t: np.ndarray, m1: float, m2: float, a: float = 1.0, e: float = 0.0, G: float = 1.0
) -> tuple[np.ndarray, np.ndarray]:
    """Exact positions and velocities of both bodies, shapes ``(T, 2, 2)``.

    The orbit lies in the xy-plane, periapsis is on the positive x-axis at
    ``t = 0`` and the motion is counter-clockwise.
    """
    t = np.atleast_1d(np.asarray(t, dtype=np.float64))
    mu = G * (m1 + m2)
    n = np.sqrt(mu / a**3)
    E = solve_kepler(n * t, e)
    cosE, sinE = np.cos(E), np.sin(E)
    b = a * np.sqrt(1.0 - e * e)
    Edot = n / (1.0 - e * cosE)
    rel_x = np.stack([a * (cosE - e), b * sinE], axis=-1)
    rel_v = np.stack([-a * sinE * Edot, b * cosE * Edot], axis=-1)
    M = m1 + m2
    x = np.stack([-(m2 / M) * rel_x, (m1 / M) * rel_x], axis=1)
    v = np.stack([-(m2 / M) * rel_v, (m1 / M) * rel_v], axis=1)
    return x, v

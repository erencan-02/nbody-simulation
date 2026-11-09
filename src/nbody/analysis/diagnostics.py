r"""Physical diagnostics: the conserved quantities of the N-body problem.

For an isolated system with pairwise central forces the following quantities
are exact invariants of the flow (see the theory chapter for derivations):

* total energy :math:`E = T + V`,
* total linear momentum :math:`\mathbf P = \sum_i m_i \mathbf v_i`,
* total angular momentum :math:`\mathbf L = \sum_i m_i\, \mathbf r_i\times\mathbf v_i`,
* the centre of mass moves uniformly, :math:`\mathbf R(t) = \mathbf R_0 + t\,\mathbf P/M`.

All functions accept arrays with arbitrary leading (time) axes, i.e. shapes
``(..., N, d)``, so they work for single states and whole trajectories alike.
"""

from __future__ import annotations

import numpy as np

from ..core import forces


def kinetic_energy(v: np.ndarray, masses: np.ndarray) -> np.ndarray:
    return 0.5 * np.einsum("...ni,...ni,n->...", v, v, masses)


def potential_energy(
    x: np.ndarray, masses: np.ndarray, G: float = 1.0, softening: float = 0.0
) -> np.ndarray:
    return np.asarray(forces.potential_energy(x, masses, G, softening))


def total_energy(
    x: np.ndarray, v: np.ndarray, masses: np.ndarray, G: float = 1.0, softening: float = 0.0
) -> np.ndarray:
    return kinetic_energy(v, masses) + potential_energy(x, masses, G, softening)


def linear_momentum(v: np.ndarray, masses: np.ndarray) -> np.ndarray:
    """Total momentum, shape ``(..., d)``."""
    return np.einsum("n,...ni->...i", masses, v)


def angular_momentum(x: np.ndarray, v: np.ndarray, masses: np.ndarray) -> np.ndarray:
    """Total angular momentum about the origin.

    In two dimensions only the scalar :math:`L_z` exists and an array of shape
    ``(...)`` is returned; in three dimensions the result has shape ``(..., 3)``.
    """
    d = x.shape[-1]
    if d == 2:
        lz = x[..., 0] * v[..., 1] - x[..., 1] * v[..., 0]
        return np.einsum("n,...n->...", masses, lz)
    if d == 3:
        return np.einsum("n,...ni->...i", masses, np.cross(x, v))
    raise ValueError("angular momentum is only defined for d = 2 or d = 3")


def center_of_mass(x: np.ndarray, masses: np.ndarray) -> np.ndarray:
    return np.einsum("n,...ni->...i", masses, x) / masses.sum()


def center_of_mass_velocity(v: np.ndarray, masses: np.ndarray) -> np.ndarray:
    return np.einsum("n,...ni->...i", masses, v) / masses.sum()


def virial_ratio(
    x: np.ndarray, v: np.ndarray, masses: np.ndarray, G: float = 1.0, softening: float = 0.0
) -> np.ndarray:
    """Virial ratio :math:`Q = T/|V|`; a relaxed self-gravitating system has Q = 1/2."""
    return kinetic_energy(v, masses) / np.abs(potential_energy(x, masses, G, softening))

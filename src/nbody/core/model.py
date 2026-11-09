r"""Force models.

A *force model* defines a separable Hamiltonian

.. math::  H(\mathbf q, \mathbf p) = T(\mathbf p) + V(\mathbf q),
           \qquad T = \sum_i \frac{|\mathbf p_i|^2}{2 m_i},

through its potential energy :math:`V` and the corresponding accelerations
:math:`\mathbf a_i = -\nabla_i V / m_i`.  Integrators only talk to this
interface, so they are completely independent of the force law: the same
leapfrog code integrates planets and harmonic oscillators.

Every call of :meth:`ForceModel.acceleration` increments ``n_evals``.  Force
evaluations dominate the cost of every method studied here, which makes this
counter the natural, hardware-independent cost measure for the
work-precision comparisons.
"""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np

from . import forces
from .state import State


class ForceModel(ABC):
    """Abstract separable Hamiltonian system ``H = T(v) + V(x)``."""

    def __init__(self, masses: np.ndarray) -> None:
        m = np.array(masses, dtype=np.float64).ravel()
        if m.size == 0:
            raise ValueError("a model needs at least one body")
        if np.any(m < 0) or not np.all(np.isfinite(m)):
            raise ValueError("masses must be finite and non-negative")
        self.masses = m
        self.n_evals = 0

    # -- interface -------------------------------------------------------------
    @abstractmethod
    def _acceleration(self, x: np.ndarray) -> np.ndarray:
        """Accelerations at positions ``x`` (shape ``(N, d)``)."""

    @abstractmethod
    def potential_energy(self, x: np.ndarray) -> np.ndarray | float:
        """Potential energy; ``x`` may have leading batch axes ``(..., N, d)``."""

    @abstractmethod
    def copy(self) -> "ForceModel":
        """Independent copy with a fresh evaluation counter."""

    # -- shared functionality ----------------------------------------------------
    def acceleration(self, x: np.ndarray) -> np.ndarray:
        self.n_evals += 1
        return self._acceleration(x)

    @property
    def n_bodies(self) -> int:
        return self.masses.size

    @property
    def total_mass(self) -> float:
        return float(self.masses.sum())

    def kinetic_energy(self, v: np.ndarray) -> np.ndarray | float:
        """Kinetic energy; ``v`` may have leading batch axes ``(..., N, d)``."""
        e = 0.5 * np.einsum("...ni,...ni,n->...", v, v, self.masses)
        return float(e) if np.ndim(e) == 0 else e

    def energy(self, state: State) -> float:
        return float(self.kinetic_energy(state.v) + self.potential_energy(state.x))

    def reset_counter(self) -> None:
        self.n_evals = 0


class GravityModel(ForceModel):
    """Newtonian gravity with optional Plummer softening.

    Parameters
    ----------
    kernel:
        ``"dense"`` (broadcast over the full ``N x N`` interaction matrix),
        ``"pairs"`` (pair list with precomputed scatter matrix, fastest for
        small ``N``) or ``"auto"`` (pairs for ``N <= PAIR_KERNEL_MAX_N``).
    """

    PAIR_KERNEL_MAX_N = 24

    def __init__(
        self,
        masses: np.ndarray,
        G: float = 1.0,
        softening: float = 0.0,
        kernel: str = "auto",
    ) -> None:
        super().__init__(masses)
        if G <= 0:
            raise ValueError("G must be positive")
        if softening < 0:
            raise ValueError("softening must be non-negative")
        if kernel not in ("auto", "dense", "pairs"):
            raise ValueError("kernel must be 'auto', 'dense' or 'pairs'")
        self.G = float(G)
        self.softening = float(softening)
        self.kernel = kernel
        use_pairs = kernel == "pairs" or (
            kernel == "auto" and 2 <= self.n_bodies <= self.PAIR_KERNEL_MAX_N
        )
        self._pairs = forces.PairKernel(self.masses) if use_pairs else None

    def _acceleration(self, x: np.ndarray) -> np.ndarray:
        if self._pairs is not None:
            return self._pairs(x, self.G, self.softening)
        return forces.acceleration(x, self.masses, self.G, self.softening)

    def potential_energy(self, x: np.ndarray) -> np.ndarray | float:
        return forces.potential_energy(x, self.masses, self.G, self.softening)

    def copy(self) -> "GravityModel":
        return GravityModel(self.masses, self.G, self.softening, self.kernel)

    def __repr__(self) -> str:
        return (
            f"GravityModel(N={self.n_bodies}, G={self.G:g}, "
            f"softening={self.softening:g})"
        )


class HarmonicModel(ForceModel):
    r"""Independent isotropic harmonic oscillators, :math:`\ddot{\mathbf x} = -\omega^2 \mathbf x`.

    Not part of the N-body problem itself; it is the standard linear test
    problem of numerical analysis and is used to illustrate stability and
    phase-space behaviour of the integrators.
    """

    def __init__(self, masses: np.ndarray, omega: float = 1.0) -> None:
        super().__init__(masses)
        self.omega = float(omega)

    def _acceleration(self, x: np.ndarray) -> np.ndarray:
        return -(self.omega**2) * x

    def potential_energy(self, x: np.ndarray) -> np.ndarray | float:
        e = 0.5 * self.omega**2 * np.einsum("...ni,...ni,n->...", x, x, self.masses)
        return float(e) if np.ndim(e) == 0 else e

    def copy(self) -> "HarmonicModel":
        return HarmonicModel(self.masses, self.omega)

r"""Explicit Runge-Kutta methods defined by a Butcher tableau.

The second-order system :math:`\dot x = v,\ \dot v = a(x)` is integrated as a
first-order ODE for :math:`y = (x, v)`.  Stage :math:`i` reads

.. math::

    X_i = x_n + h\sum_{j<i} a_{ij} V_j, \qquad
    V_i = v_n + h\sum_{j<i} a_{ij} A_j, \qquad A_i = a(X_i),

and the update is :math:`x_{n+1} = x_n + h\sum_i b_i V_i`,
:math:`v_{n+1} = v_n + h\sum_i b_i A_i`.  Adding a new explicit method only
requires a new tableau.
"""

from __future__ import annotations

import numpy as np

from ..core.model import ForceModel
from ..core.state import State
from .base import Integrator, cached_acceleration


class ExplicitRungeKutta(Integrator):
    """Generic explicit Runge-Kutta method for ``x' = v, v' = a(x)``."""

    A: np.ndarray
    b: np.ndarray
    c: np.ndarray

    def __init_subclass__(cls, **kwargs) -> None:
        super().__init_subclass__(**kwargs)
        if hasattr(cls, "A"):
            cls.A = np.asarray(cls.A, dtype=np.float64)
            cls.b = np.asarray(cls.b, dtype=np.float64)
            cls.c = np.asarray(cls.c, dtype=np.float64)
            s = cls.b.size
            if cls.A.shape != (s, s) or cls.c.size != s:
                raise ValueError(f"inconsistent Butcher tableau in {cls.__name__}")
            if np.any(np.triu(cls.A) != 0):
                raise ValueError(f"{cls.__name__}: tableau is not explicit")
            if not np.allclose(cls.A.sum(axis=1), cls.c):
                raise ValueError(f"{cls.__name__}: row sums of A must equal c")
            cls.evals_per_step = s

    def stages(self, state: State, dt: float, model: ForceModel):
        """Return the stage derivatives ``(V_i, A_i)`` as two lists."""
        x0, v0 = state.x, state.v
        A = self.A
        kx: list[np.ndarray] = []
        kv: list[np.ndarray] = []
        for i in range(self.b.size):
            if i == 0:
                xi, vi = x0, v0
                ai = cached_acceleration(state, model)
            else:
                xi = x0.copy()
                vi = v0.copy()
                for j in range(i):
                    if A[i, j] != 0.0:
                        xi += (dt * A[i, j]) * kx[j]
                        vi += (dt * A[i, j]) * kv[j]
                ai = model.acceleration(xi)
            kx.append(vi)
            kv.append(ai)
        return kx, kv

    def step(self, state: State, dt: float, model: ForceModel) -> None:
        kx, kv = self.stages(state, dt, model)
        x_new = state.x.copy()
        v_new = state.v.copy()
        for bi, kxi, kvi in zip(self.b, kx, kv):
            if bi != 0.0:
                x_new += (dt * bi) * kxi
                v_new += (dt * bi) * kvi
        state.x = x_new
        state.v = v_new
        state.a = None
        state.t += dt


class Heun(ExplicitRungeKutta):
    """Heun's method / explicit trapezoidal rule (RK2, non-symplectic)."""

    name = "heun"
    label = "Heun (RK2)"
    order = 2
    symplectic = False
    time_reversible = False
    A = [[0.0, 0.0], [1.0, 0.0]]
    b = [0.5, 0.5]
    c = [0.0, 1.0]


class RK4(ExplicitRungeKutta):
    """Classical fourth-order Runge-Kutta method (non-symplectic)."""

    name = "rk4"
    label = "Runge-Kutta 4"
    order = 4
    symplectic = False
    time_reversible = False
    A = [
        [0.0, 0.0, 0.0, 0.0],
        [0.5, 0.0, 0.0, 0.0],
        [0.0, 0.5, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
    ]
    b = [1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0]
    c = [0.0, 0.5, 0.5, 1.0]

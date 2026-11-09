r"""First-order methods: explicit Euler and symplectic Euler.

Both methods use exactly one force evaluation per step and are of global
order one.  They differ only in *which* velocity updates the position, but
this tiny difference decides whether the discrete flow is symplectic:

* explicit Euler:   :math:`x_{n+1} = x_n + h v_n,\; v_{n+1} = v_n + h a(x_n)`
* symplectic Euler: :math:`v_{n+1} = v_n + h a(x_n),\; x_{n+1} = x_n + h v_{n+1}`
"""

from __future__ import annotations

from ..core.model import ForceModel
from ..core.state import State
from .base import Integrator, cached_acceleration


class ExplicitEuler(Integrator):
    """Forward Euler method (non-symplectic, order 1)."""

    name = "euler"
    label = "Explicit Euler"
    order = 1
    symplectic = False
    time_reversible = False
    evals_per_step = 1

    def step(self, state: State, dt: float, model: ForceModel) -> None:
        a = cached_acceleration(state, model)  # a(x_n)
        state.x += dt * state.v                # uses v_n
        state.v += dt * a                      # uses a(x_n)
        state.a = None
        state.t += dt


class SymplecticEuler(Integrator):
    """Symplectic (semi-implicit) Euler, kick-then-drift variant (order 1)."""

    name = "symplectic_euler"
    label = "Symplectic Euler"
    order = 1
    symplectic = True
    time_reversible = False
    evals_per_step = 1

    def step(self, state: State, dt: float, model: ForceModel) -> None:
        a = cached_acceleration(state, model)  # a(x_n)
        state.v += dt * a                      # v_{n+1}
        state.x += dt * state.v                # uses v_{n+1}
        state.a = None
        state.t += dt

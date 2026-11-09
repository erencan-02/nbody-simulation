r"""Störmer-Verlet / leapfrog in kick-drift-kick form (velocity Verlet).

.. math::

    v_{n+1/2} &= v_n + \tfrac h2\, a(x_n) \\
    x_{n+1}   &= x_n + h\, v_{n+1/2} \\
    v_{n+1}   &= v_{n+1/2} + \tfrac h2\, a(x_{n+1})

The scheme is the Strang splitting of :math:`H = T + V`; it is symplectic,
time-reversible and of order two.  The acceleration at the end of a step is
the acceleration at the start of the next one, so it is cached in the state
("first same as last") and every step costs a single force evaluation.  The
original course implementation evaluated the forces twice per step.
"""

from __future__ import annotations

from ..core.model import ForceModel
from ..core.state import State
from .base import Integrator, cached_acceleration


class VelocityVerlet(Integrator):
    """Leapfrog / velocity Verlet (symplectic, symmetric, order 2)."""

    name = "leapfrog"
    label = "Leapfrog (Verlet)"
    order = 2
    symplectic = True
    time_reversible = True
    evals_per_step = 1

    def step(self, state: State, dt: float, model: ForceModel) -> None:
        half = 0.5 * dt
        state.v += half * cached_acceleration(state, model)  # kick
        state.x += dt * state.v                              # drift
        a_new = model.acceleration(state.x)
        state.v += half * a_new                              # kick
        state.a = a_new                                      # FSAL cache
        state.t += dt

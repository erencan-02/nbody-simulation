r"""Higher-order symplectic integrators by composition of leapfrog steps.

If :math:`\Phi_h` is a symmetric method of even order :math:`p`, the
*triple jump*

.. math::  \Psi_h = \Phi_{\gamma_1 h}\circ\Phi_{\gamma_2 h}\circ\Phi_{\gamma_1 h},
   \qquad \gamma_1 = \frac{1}{2-2^{1/(p+1)}},\quad
          \gamma_2 = -\frac{2^{1/(p+1)}}{2-2^{1/(p+1)}},

is symmetric and of order :math:`p+2` (Creutz & Gocksch 1989, Yoshida 1990,
Suzuki 1990).  The conditions follow from
:math:`2\gamma_1+\gamma_2 = 1` (consistency) and
:math:`2\gamma_1^{p+1}+\gamma_2^{p+1} = 0` (cancellation of the leading error
term).  Since compositions of symplectic maps are symplectic, the result
inherits symplecticity from leapfrog.

Applying the construction to leapfrog yields Yoshida's fourth-order method;
applying it again yields a sixth-order method with :math:`3^2 = 9` stages.
Yoshida's optimised sixth-order "solution A" needs only seven stages.
"""

from __future__ import annotations

from collections.abc import Sequence

from ..core.model import ForceModel
from ..core.state import State
from .base import Integrator
from .verlet import VelocityVerlet


def triple_jump_weights(base_order: int) -> tuple[float, float, float]:
    """Triple-jump weights raising a symmetric method of order ``p`` to ``p+2``."""
    if base_order % 2 != 0 or base_order < 2:
        raise ValueError("the base method must be symmetric with even order >= 2")
    k = 2.0 ** (1.0 / (base_order + 1))
    g1 = 1.0 / (2.0 - k)
    g2 = -k * g1
    return (g1, g2, g1)


class Composition(Integrator):
    """Composition ``Phi_{w_s h} o ... o Phi_{w_1 h}`` of a base method."""

    def __init__(
        self,
        base: Integrator,
        weights: Sequence[float],
        *,
        order: int,
        name: str = "composition",
        label: str = "Composition",
    ) -> None:
        self.base = base
        self.weights = tuple(float(w) for w in weights)
        if abs(sum(self.weights) - 1.0) > 1e-12:
            raise ValueError("composition weights must sum to one (consistency)")
        self.order = order
        self.name = name
        self.label = label
        self.symplectic = base.symplectic
        palindromic = all(
            abs(a - b) < 1e-15 for a, b in zip(self.weights, reversed(self.weights))
        )
        self.time_reversible = base.time_reversible and palindromic
        # With FSAL caching each base step costs base.evals_per_step evaluations.
        self.evals_per_step = len(self.weights) * base.evals_per_step

    def step(self, state: State, dt: float, model: ForceModel) -> None:
        t0 = state.t
        for w in self.weights:
            self.base.step(state, w * dt, model)
        state.t = t0 + dt  # avoid accumulating round-off from the sub-steps

    def __repr__(self) -> str:
        return f"{type(self).__name__}(base={self.base!r}, stages={len(self.weights)})"


class Yoshida4(Composition):
    """Fourth-order symplectic method: triple jump of leapfrog (3 stages)."""

    def __init__(self) -> None:
        super().__init__(
            VelocityVerlet(),
            triple_jump_weights(2),
            order=4,
            name="yoshida4",
            label="Yoshida 4",
        )


class Yoshida6(Composition):
    """Sixth-order symplectic method: triple jump of :class:`Yoshida4` (9 stages)."""

    def __init__(self) -> None:
        y4 = Yoshida4()
        g = triple_jump_weights(4)
        weights = [gi * wj for gi in g for wj in y4.weights]
        super().__init__(
            VelocityVerlet(),
            weights,
            order=6,
            name="yoshida6",
            label="Yoshida 6 (triple jump)",
        )


# Yoshida (1990), Phys. Lett. A 150, 262, Table 1, solution A.
_YOSHIDA6A_W1 = -1.17767998417887
_YOSHIDA6A_W2 = 0.235573213359357
_YOSHIDA6A_W3 = 0.784513610477560
_YOSHIDA6A_W0 = 1.0 - 2.0 * (_YOSHIDA6A_W1 + _YOSHIDA6A_W2 + _YOSHIDA6A_W3)


class Yoshida6A(Composition):
    """Yoshida's optimised sixth-order method ("solution A", 7 stages)."""

    def __init__(self) -> None:
        w0, w1, w2, w3 = _YOSHIDA6A_W0, _YOSHIDA6A_W1, _YOSHIDA6A_W2, _YOSHIDA6A_W3
        super().__init__(
            VelocityVerlet(),
            (w3, w2, w1, w0, w1, w2, w3),
            order=6,
            name="yoshida6a",
            label="Yoshida 6 (solution A)",
        )

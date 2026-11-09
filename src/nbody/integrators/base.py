"""Common interface of all time integrators."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np

from ..core.model import ForceModel
from ..core.state import State


class Integrator(ABC):
    """A one-step method advancing a :class:`~nbody.core.state.State` in place.

    Subclasses declare their theoretical properties as attributes; they are
    used for documentation, for the evaluation tables and in the tests that
    check the implementation against the theory.

    Attributes
    ----------
    name:
        Registry key (e.g. ``"leapfrog"``).
    label:
        Human-readable name used in plots.
    order:
        Global order of accuracy ``p`` (global error ``O(dt^p)``).
    symplectic:
        Whether the one-step map is a symplectic transformation.
    time_reversible:
        Whether the method is symmetric, i.e. ``Phi_{-dt} o Phi_{dt} = id``.
    evals_per_step:
        Force evaluations per step in steady state (after the first step).
    """

    name: str = "integrator"
    label: str = "Integrator"
    order: int = 0
    symplectic: bool = False
    time_reversible: bool = False
    evals_per_step: int = 1
    adaptive: bool = False

    @abstractmethod
    def step(self, state: State, dt: float, model: ForceModel) -> None:
        """Advance ``state`` by one step of size ``dt`` (in place)."""

    def properties(self) -> dict:
        return {
            "name": self.name,
            "label": self.label,
            "order": self.order,
            "symplectic": self.symplectic,
            "time_reversible": self.time_reversible,
            "evals_per_step": self.evals_per_step,
            "adaptive": self.adaptive,
        }

    def __repr__(self) -> str:
        return f"{type(self).__name__}()"


@dataclass
class AdaptiveStepResult:
    """Outcome of one accepted step of an adaptive integrator."""

    dt_taken: float
    dt_next: float
    n_rejected: int


class AdaptiveIntegrator(Integrator):
    """Integrator with embedded error estimation and step-size control.

    :meth:`step` may take a *smaller* step than requested (after rejecting
    trial steps); the actually used step and a proposal for the next step are
    returned.
    """

    adaptive = True

    @abstractmethod
    def step(self, state: State, dt: float, model: ForceModel) -> AdaptiveStepResult:
        ...


def cached_acceleration(state: State, model: ForceModel) -> np.ndarray:
    """Acceleration at ``state.x``, re-using the FSAL cache when available."""
    if state.a is None:
        state.a = model.acceleration(state.x)
    return state.a

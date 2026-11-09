"""Phase-space state of an N-body system.

The state is stored as a *structure of arrays*: one ``(N, d)`` array for all
positions and one for all velocities.  This layout lets every numerical kernel
(force evaluation, integrator updates, diagnostics) operate on whole arrays at
once instead of looping over Python objects.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class State:
    """Positions, velocities and time of an N-body system.

    Attributes
    ----------
    x:
        Positions, shape ``(N, d)``.
    v:
        Velocities, shape ``(N, d)``.
    t:
        Current time.
    a:
        Optional cached acceleration evaluated at ``x`` (shape ``(N, d)``).
        Integrators whose last stage ends at the new position (velocity Verlet,
        Dormand-Prince) store it here so that the next step can reuse it
        ("first same as last").  Any code that modifies ``x`` must reset it
        to ``None``.
    """

    x: np.ndarray
    v: np.ndarray
    t: float = 0.0
    a: np.ndarray | None = field(default=None, repr=False)

    def __post_init__(self) -> None:
        self.x = np.array(self.x, dtype=np.float64)
        self.v = np.array(self.v, dtype=np.float64)
        if self.x.ndim != 2:
            raise ValueError(f"positions must have shape (N, d), got {self.x.shape}")
        if self.x.shape != self.v.shape:
            raise ValueError(
                f"positions {self.x.shape} and velocities {self.v.shape} differ in shape"
            )
        self.t = float(self.t)
        if self.a is not None:
            self.a = np.array(self.a, dtype=np.float64)
            if self.a.shape != self.x.shape:
                raise ValueError("cached acceleration has the wrong shape")

    @property
    def n_bodies(self) -> int:
        return self.x.shape[0]

    @property
    def dim(self) -> int:
        return self.x.shape[1]

    def copy(self) -> "State":
        """Deep copy (all arrays are copied by ``__post_init__``)."""
        return State(self.x, self.v, self.t, self.a)

    def invalidate_cache(self) -> None:
        """Forget the cached acceleration (call after modifying ``x`` by hand)."""
        self.a = None

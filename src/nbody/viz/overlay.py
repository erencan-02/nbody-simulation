"""Back-end independent logic of the trajectory overlays.

Two overlays are provided:

* :class:`TrailBuffer` - a fixed-size ring buffer with the recent *past*
  positions of every body, drawn as fading trails;
* :class:`TrajectoryPredictor` - a "ghost" integration of a copy of the
  current system that *simulates the future* trajectories of all bodies over
  a configurable horizon without touching the running simulation.

Both are pure NumPy and can be unit-tested without a display.
"""

from __future__ import annotations

import math

import numpy as np

from ..simulation import Simulation


class TrailBuffer:
    """Ring buffer holding the last ``capacity`` positions of ``n`` bodies."""

    def __init__(self, n_bodies: int, dim: int, capacity: int = 600) -> None:
        if capacity < 2:
            raise ValueError("capacity must be at least 2")
        self.capacity = int(capacity)
        self._data = np.zeros((self.capacity, n_bodies, dim))
        self._head = 0   # index of the next write
        self._size = 0

    def __len__(self) -> int:
        return self._size

    def push(self, x: np.ndarray) -> None:
        self._data[self._head] = x
        self._head = (self._head + 1) % self.capacity
        self._size = min(self._size + 1, self.capacity)

    def clear(self) -> None:
        self._head = 0
        self._size = 0

    def ordered(self) -> np.ndarray:
        """Stored positions from oldest to newest, shape ``(len, N, d)``."""
        if self._size < self.capacity:
            return self._data[: self._size].copy()
        return np.concatenate([self._data[self._head :], self._data[: self._head]])

    def shift(self, offset: np.ndarray) -> None:
        """Translate all stored points (used when the camera frame is re-centred)."""
        self._data[: self._size if self._size < self.capacity else None] += offset


class TrajectoryPredictor:
    """Forecast future trajectories by integrating a clone of the simulation.

    Parameters
    ----------
    horizon:
        Length of the forecast in simulation time units.
    dt:
        Step size of the ghost integration (may be coarser than the live one).
    n_points:
        Approximate number of samples returned per body.
    """

    def __init__(self, horizon: float, dt: float, n_points: int = 400) -> None:
        if horizon <= 0 or dt <= 0:
            raise ValueError("horizon and dt must be positive")
        self.horizon = float(horizon)
        self.dt = float(dt)
        self.n_points = int(n_points)
        self.last_n_evals = 0

    def predict(self, sim: Simulation) -> np.ndarray:
        """Predicted positions, shape ``(S, N, d)``; the input is not modified."""
        ghost = sim.clone()
        n_steps = max(1, math.ceil(self.horizon / self.dt))
        record_every = max(1, n_steps // max(1, self.n_points))
        traj = ghost.run(ghost.state.t + self.horizon, self.dt, record_every)
        self.last_n_evals = ghost.model.n_evals
        return traj.x

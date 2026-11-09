"""Scenario definition and JSON (de)serialisation.

JSON format (all keys except ``bodies`` are optional)::

    {
      "name": "two-body",
      "description": "...",
      "G": 1.0,
      "softening": 0.0,
      "units": {"length": "AU", "time": "yr", "mass": "M_sun"},
      "recommended": {"integrator": "leapfrog", "dt": 0.01, "t_end": 10.0},
      "bodies": [
        {"name": "Sun", "mass": 1.0, "position": [0, 0], "velocity": [0, 0],
         "color": "#FDCA00"}
      ]
    }

The format is backwards compatible with the scenario files of the original
course project, which only contained ``G`` and ``bodies``.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np

from ..analysis import diagnostics as diag
from ..core.model import GravityModel
from ..core.state import State


@dataclass
class Scenario:
    """Initial conditions plus the physical parameters of a simulation."""

    name: str
    masses: np.ndarray
    positions: np.ndarray
    velocities: np.ndarray
    G: float = 1.0
    softening: float = 0.0
    names: list[str] | None = None
    colors: list[str] | None = None
    description: str = ""
    units: dict = field(default_factory=dict)
    recommended: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.masses = np.array(self.masses, dtype=np.float64).ravel()
        self.positions = np.array(self.positions, dtype=np.float64)
        self.velocities = np.array(self.velocities, dtype=np.float64)
        n = self.masses.size
        if self.positions.ndim != 2 or self.positions.shape[0] != n:
            raise ValueError(f"positions must have shape ({n}, d)")
        if self.velocities.shape != self.positions.shape:
            raise ValueError("velocities must have the same shape as positions")
        if self.names is None:
            self.names = [f"body {i}" for i in range(n)]
        if len(self.names) != n:
            raise ValueError("one name per body is required")
        if self.colors is not None and len(self.colors) != n:
            raise ValueError("one colour per body is required")

    # -- derived objects ---------------------------------------------------------
    @property
    def n_bodies(self) -> int:
        return self.masses.size

    @property
    def dim(self) -> int:
        return self.positions.shape[1]

    def model(self) -> GravityModel:
        return GravityModel(self.masses, self.G, self.softening)

    def initial_state(self) -> State:
        return State(self.positions, self.velocities, 0.0)

    def barycentric(self) -> "Scenario":
        """Same scenario in the centre-of-mass frame (zero total momentum)."""
        R = diag.center_of_mass(self.positions, self.masses)
        V = diag.center_of_mass_velocity(self.velocities, self.masses)
        return replace(self, positions=self.positions - R, velocities=self.velocities - V)

    # -- serialisation -------------------------------------------------------------
    def to_dict(self) -> dict:
        bodies = []
        for i in range(self.n_bodies):
            b = {
                "name": self.names[i],
                "mass": float(self.masses[i]),
                "position": self.positions[i].tolist(),
                "velocity": self.velocities[i].tolist(),
            }
            if self.colors is not None:
                b["color"] = self.colors[i]
            bodies.append(b)
        return {
            "name": self.name,
            "description": self.description,
            "G": self.G,
            "softening": self.softening,
            "units": dict(self.units),
            "recommended": dict(self.recommended),
            "bodies": bodies,
        }

    @classmethod
    def from_dict(cls, data: dict, name: str = "custom") -> "Scenario":
        bodies = data.get("bodies")
        if not bodies:
            raise ValueError("a scenario needs a non-empty 'bodies' list")
        colors = [b.get("color") for b in bodies]
        return cls(
            name=data.get("name", name),
            masses=[b["mass"] for b in bodies],
            positions=[b["position"] for b in bodies],
            velocities=[b["velocity"] for b in bodies],
            G=float(data.get("G", 1.0)),
            softening=float(data.get("softening", 0.0)),
            names=[b.get("name", f"body {i}") for i, b in enumerate(bodies)],
            colors=colors if all(c is not None for c in colors) else None,
            description=data.get("description", ""),
            units=data.get("units", {}),
            recommended=data.get("recommended", {}),
        )

    def save_json(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")

    @classmethod
    def load_json(cls, path: str | Path) -> "Scenario":
        path = Path(path)
        data = json.loads(path.read_text(encoding="utf-8"))
        return cls.from_dict(data, name=path.stem)

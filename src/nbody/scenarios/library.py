"""Built-in scenarios.

Every builtin is a function returning a :class:`Scenario`; keyword arguments
parameterise it (e.g. ``load_scenario("kepler", e=0.9)``).  Scenarios are
returned in the centre-of-mass frame unless stated otherwise.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np

from ..analysis import diagnostics as diag
from ..analysis.kepler import kepler_initial_conditions, orbital_period
from .scenario import Scenario

_BUILTINS: dict[str, Callable[..., Scenario]] = {}

# TU Darmstadt corporate design colours (used for bodies in plots/viewer)
_PALETTE = ["#FDCA00", "#0083CC", "#E6001A", "#009D81", "#A60084",
            "#EC6500", "#99C000", "#721085", "#5D85C3", "#B5B5B5"]


def builtin(name: str):
    """Decorator registering a scenario factory under ``name``."""

    def deco(fn: Callable[..., Scenario]) -> Callable[..., Scenario]:
        _BUILTINS[name] = fn
        return fn

    return deco


def available_scenarios() -> list[str]:
    return list(_BUILTINS)


def scenario_description(name: str) -> str:
    doc = _BUILTINS[name].__doc__ or ""
    return doc.strip().splitlines()[0] if doc.strip() else ""


def load_scenario(name_or_path: str | Path, **params) -> Scenario:
    """Load a builtin scenario by name or a JSON scenario file by path."""
    key = str(name_or_path)
    if key in _BUILTINS:
        return _BUILTINS[key](**params)
    path = Path(key)
    if path.suffix.lower() == ".json" or path.exists():
        if params:
            raise TypeError("JSON scenarios do not take parameters")
        return Scenario.load_json(path)
    raise KeyError(
        f"unknown scenario '{key}'; builtins: {', '.join(available_scenarios())}"
    )


# ---------------------------------------------------------------------------
@builtin("kepler")
def kepler(e: float = 0.6, mass_ratio: float = 1e-3, a: float = 1.0) -> Scenario:
    """Two-body Kepler orbit (exact solution known), G = M = a = 1, P = 2 pi."""
    m1, m2 = 1.0, float(mass_ratio)
    x, v = kepler_initial_conditions(m1, m2, a=a, e=e, G=1.0)
    P = orbital_period(a, m1, m2)
    return Scenario(
        name="kepler",
        masses=[m1, m2],
        positions=x,
        velocities=v,
        G=1.0,
        names=["star", "planet"],
        colors=[_PALETTE[0], _PALETTE[1]],
        description=f"Kepler orbit with a={a}, e={e}, m2/m1={mass_ratio}; "
        "starts at periapsis.",
        recommended={"integrator": "leapfrog", "dt": P / 500, "t_end": 10 * P},
    )


@builtin("figure_eight")
def figure_eight() -> Scenario:
    """Chenciner-Montgomery figure-eight choreography of three equal masses."""
    x1 = np.array([0.97000436, -0.24308753])
    v3 = np.array([-0.93240737, -0.86473146])
    positions = np.array([x1, -x1, [0.0, 0.0]])
    velocities = np.array([-0.5 * v3, -0.5 * v3, v3])
    return Scenario(
        name="figure_eight",
        masses=[1.0, 1.0, 1.0],
        positions=positions,
        velocities=velocities,
        G=1.0,
        names=["A", "B", "C"],
        colors=_PALETTE[1:4],
        description="Periodic three-body orbit (Chenciner & Montgomery 2000), "
        "period T = 6.32591398.",
        recommended={"integrator": "yoshida4", "dt": 6.32591398 / 400, "t_end": 6.32591398 * 10},
    )


@builtin("pythagorean")
def pythagorean() -> Scenario:
    """Burrau's Pythagorean three-body problem (masses 3, 4, 5 at rest)."""
    return Scenario(
        name="pythagorean",
        masses=[3.0, 4.0, 5.0],
        positions=[[1.0, 3.0], [-2.0, -1.0], [1.0, -1.0]],
        velocities=np.zeros((3, 2)),
        G=1.0,
        names=["m=3", "m=4", "m=5"],
        colors=[_PALETTE[1], _PALETTE[2], _PALETTE[3]],
        description="Chaotic three-body problem with very close encounters "
        "(Szebehely & Peters 1967); ends in a binary plus an escaper.",
        recommended={"integrator": "dopri5", "dt": 1e-3, "t_end": 70.0},
    )


# J2000 mean orbital elements (Standish, JPL "Approximate positions of the
# planets", Table 1): semi-major axis [AU], eccentricity, longitude of
# perihelion [deg]; masses as inverse ratios Sun/planet (IAU 2009/2012).
_PLANETS = {
    #           a [AU]       e            varpi [deg]     M_sun / m
    "Mercury": (0.38709927, 0.20563593, 77.45779628, 6023600.0),
    "Venus": (0.72333566, 0.00677672, 131.60246718, 408523.71),
    "Earth": (1.00000261, 0.01671123, 102.93768193, 328900.56),  # Earth+Moon
    "Mars": (1.52371034, 0.09339410, -23.94362959, 3098703.59),
    "Jupiter": (5.20288700, 0.04838624, 14.72847983, 1047.348644),
    "Saturn": (9.53667594, 0.05386179, 92.59887831, 3497.9018),
    "Uranus": (19.18916464, 0.04725744, 170.95427630, 22902.98),
    "Neptune": (30.06992276, 0.00859048, 44.96476227, 19412.26),
}
_PLANET_COLORS = {
    "Sun": "#FDCA00", "Mercury": "#B5B5B5", "Venus": "#F5A300", "Earth": "#0083CC",
    "Mars": "#E6001A", "Jupiter": "#EC6500", "Saturn": "#D7AC00", "Uranus": "#50B695",
    "Neptune": "#5D85C3",
}


def planetary_elements() -> dict[str, tuple[float, float, float, float]]:
    """Orbital elements used by :func:`solar_system` (a, e, varpi, mass ratio)."""
    return dict(_PLANETS)


@builtin("solar_system")
def solar_system(planets: tuple[str, ...] | None = None) -> Scenario:
    """Sun and eight planets in AU / yr / M_sun (coplanar, all at perihelion)."""
    names = list(planets) if planets is not None else list(_PLANETS)
    G = 4.0 * np.pi**2  # AU^3 / (M_sun yr^2)
    masses = [1.0]
    pos = [[0.0, 0.0]]
    vel = [[0.0, 0.0]]
    for pname in names:
        a, e, varpi, inv_mass = _PLANETS[pname]
        m = 1.0 / inv_mass
        w = np.deg2rad(varpi)
        r_p = a * (1.0 - e)
        v_p = np.sqrt(G * (1.0 + m) * (1.0 + e) / (a * (1.0 - e)))  # vis-viva
        u = np.array([np.cos(w), np.sin(w)])
        masses.append(m)
        pos.append(r_p * u)
        vel.append(v_p * np.array([-u[1], u[0]]))
    sc = Scenario(
        name="solar_system" if planets is None else "solar_system_subset",
        masses=masses,
        positions=pos,
        velocities=vel,
        G=G,
        names=["Sun"] + names,
        colors=[_PLANET_COLORS[n] for n in ["Sun"] + names],
        description="Simplified coplanar Solar System from J2000 mean elements; "
        "every planet starts at its perihelion.",
        units={"length": "AU", "time": "yr", "mass": "M_sun"},
        recommended={"integrator": "leapfrog", "dt": 1.0 / 365.25, "t_end": 200.0},
    )
    return sc.barycentric()


@builtin("inner_solar_system")
def inner_solar_system() -> Scenario:
    """Sun, Mercury, Venus, Earth, Mars and Jupiter (fast-moving demo)."""
    sc = solar_system(("Mercury", "Venus", "Earth", "Mars", "Jupiter"))
    sc.name = "inner_solar_system"
    sc.recommended = {"integrator": "leapfrog", "dt": 1.0 / 365.25, "t_end": 30.0}
    return sc


@builtin("circumbinary")
def circumbinary(q: float = 0.5, r_planet: float = 3.5) -> Scenario:
    """Circular binary (separation 1) with a light circumbinary planet."""
    m1, m2 = 1.0 / (1.0 + q), q / (1.0 + q)
    xb, vb = kepler_initial_conditions(m1, m2, a=1.0, e=0.0)
    vp = np.sqrt(1.0 / r_planet)
    masses = [m1, m2, 1e-6]
    positions = np.vstack([xb, [[r_planet, 0.0]]])
    velocities = np.vstack([vb, [[0.0, vp]]])
    sc = Scenario(
        name="circumbinary",
        masses=masses,
        positions=positions,
        velocities=velocities,
        G=1.0,
        names=["star A", "star B", "planet"],
        colors=[_PALETTE[0], _PALETTE[5], _PALETTE[1]],
        description="Planet on an initially circular orbit around a binary star.",
        recommended={"integrator": "yoshida4", "dt": 2 * np.pi / 400, "t_end": 200.0},
    )
    return sc.barycentric()


@builtin("random_cluster")
def random_cluster(
    n: int = 50, seed: int = 1, softening: float = 0.05, virial_ratio: float = 0.5, dim: int = 2
) -> Scenario:
    """Gaussian star cluster of equal masses in virial equilibrium (Q = 1/2)."""
    rng = np.random.default_rng(seed)
    masses = np.full(n, 1.0 / n)
    x = rng.normal(size=(n, dim))
    v = rng.normal(size=(n, dim))
    x -= diag.center_of_mass(x, masses)
    v -= diag.center_of_mass_velocity(v, masses)
    T = float(diag.kinetic_energy(v, masses))
    V = float(diag.potential_energy(x, masses, 1.0, softening))
    v *= np.sqrt(virial_ratio * abs(V) / T)  # scale to T = Q |V|
    return Scenario(
        name="random_cluster",
        masses=masses,
        positions=x,
        velocities=v,
        G=1.0,
        softening=softening,
        names=[f"star {i}" for i in range(n)],
        colors=[_PALETTE[i % len(_PALETTE)] for i in range(n)],
        description=f"{n} equal-mass stars, Gaussian positions, virial ratio {virial_ratio}.",
        recommended={"integrator": "leapfrog", "dt": 1e-3, "t_end": 10.0},
    )

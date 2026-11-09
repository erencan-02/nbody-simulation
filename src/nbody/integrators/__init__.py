"""Time integrators, ordered by increasing sophistication.

=====================  =====  ==========  ==========  ===========
method                 order  symplectic  reversible  evals/step
=====================  =====  ==========  ==========  ===========
explicit Euler         1      no          no          1
symplectic Euler       1      yes         no          1
Heun (RK2)             2      no          no          2
leapfrog (Verlet)      2      yes         yes         1
Runge-Kutta 4          4      no          no          4
Yoshida 4              4      yes         yes         3
Yoshida 6 (A)          6      yes         yes         7
Yoshida 6 (triple)     6      yes         yes         9
Dormand-Prince 5(4)    5      no          no          6 (adaptive)
=====================  =====  ==========  ==========  ===========

New methods are added by subclassing :class:`Integrator` (or one of the
generic building blocks :class:`ExplicitRungeKutta` and :class:`Composition`)
and calling :func:`register`.
"""

from __future__ import annotations

from collections.abc import Callable

from .adaptive import DormandPrince54
from .base import AdaptiveIntegrator, AdaptiveStepResult, Integrator
from .composition import Composition, Yoshida4, Yoshida6, Yoshida6A, triple_jump_weights
from .euler import ExplicitEuler, SymplecticEuler
from .runge_kutta import RK4, ExplicitRungeKutta, Heun
from .verlet import VelocityVerlet

_REGISTRY: dict[str, Callable[..., Integrator]] = {}
_ALIASES: dict[str, str] = {}


def register(name: str, factory: Callable[..., Integrator], aliases: tuple[str, ...] = ()) -> None:
    """Make an integrator available under ``name`` (and optional aliases)."""
    key = name.lower()
    if key in _REGISTRY:
        raise KeyError(f"integrator '{name}' is already registered")
    _REGISTRY[key] = factory
    for alias in aliases:
        _ALIASES[alias.lower()] = key


def create_integrator(name: str, **kwargs) -> Integrator:
    """Instantiate a registered integrator by name or alias."""
    key = name.lower().replace("-", "_")
    key = _ALIASES.get(key, key)
    try:
        factory = _REGISTRY[key]
    except KeyError:
        raise KeyError(
            f"unknown integrator '{name}'; available: {', '.join(available_integrators())}"
        ) from None
    return factory(**kwargs)


def available_integrators() -> list[str]:
    """Registered integrator names in order of increasing sophistication."""
    return list(_REGISTRY)


register("euler", ExplicitEuler, aliases=("explicit_euler", "forward_euler"))
register("symplectic_euler", SymplecticEuler, aliases=("semi_implicit_euler",))
register("heun", Heun, aliases=("rk2",))
register("leapfrog", VelocityVerlet, aliases=("verlet", "velocity_verlet", "stormer_verlet"))
register("rk4", RK4, aliases=("runge_kutta", "runge_kutta4"))
register("yoshida4", Yoshida4, aliases=("forest_ruth",))
register("yoshida6a", Yoshida6A)
register("yoshida6", Yoshida6)
register("dopri5", DormandPrince54, aliases=("rk45", "dormand_prince"))

__all__ = [
    "Integrator",
    "AdaptiveIntegrator",
    "AdaptiveStepResult",
    "ExplicitEuler",
    "SymplecticEuler",
    "Heun",
    "VelocityVerlet",
    "RK4",
    "ExplicitRungeKutta",
    "Composition",
    "Yoshida4",
    "Yoshida6",
    "Yoshida6A",
    "DormandPrince54",
    "triple_jump_weights",
    "register",
    "create_integrator",
    "available_integrators",
]

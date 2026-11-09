"""nbody - structure-preserving integration of the gravitational N-body problem.

Quick start::

    from nbody import load_scenario, simulate
    from nbody.analysis import conservation_report

    sc = load_scenario("kepler", e=0.6)
    traj = simulate(sc.model(), sc.initial_state(), "yoshida4", t_end=100.0, dt=0.01)
    print(conservation_report(traj))
"""

from .core import ForceModel, GravityModel, HarmonicModel, State
from .integrators import available_integrators, create_integrator, register
from .scenarios import Scenario, available_scenarios, load_scenario
from .simulation import Simulation, Trajectory, simulate

__version__ = "1.0.0"

__all__ = [
    "State",
    "ForceModel",
    "GravityModel",
    "HarmonicModel",
    "Scenario",
    "Simulation",
    "Trajectory",
    "simulate",
    "load_scenario",
    "available_scenarios",
    "available_integrators",
    "create_integrator",
    "register",
]

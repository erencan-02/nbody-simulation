"""Visualisation: plot style, static plots and the interactive viewer.

The pygame viewer is imported lazily (``from nbody.viz.viewer import Viewer``)
so that the rest of the package works without a display or pygame.
"""

from .overlay import TrailBuffer, TrajectoryPredictor
from .style import INTEGRATOR_STYLE, TUDA, apply_style, integrator_style

__all__ = [
    "TrailBuffer",
    "TrajectoryPredictor",
    "INTEGRATOR_STYLE",
    "TUDA",
    "apply_style",
    "integrator_style",
]

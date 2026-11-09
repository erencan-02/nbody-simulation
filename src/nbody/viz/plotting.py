"""Static matplotlib plots of trajectories and diagnostics."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from ..analysis.metrics import relative_energy_error
from .style import integrator_style


def plot_orbits(
    ax,
    traj,
    names=None,
    colors=None,
    lw: float = 0.8,
    markers: bool = True,
    max_points: int = 20000,
    alpha: float = 1.0,
    label_bodies: bool = True,
):
    """Draw the orbits of all bodies of a trajectory (xy projection)."""
    stride = max(1, traj.n_frames // max_points)
    x = traj.x[::stride]
    for i in range(traj.n_bodies):
        kw = {}
        if colors is not None:
            kw["color"] = colors[i]
        label = names[i] if (names is not None and label_bodies) else None
        (line,) = ax.plot(x[:, i, 0], x[:, i, 1], lw=lw, alpha=alpha, label=label, **kw)
        if markers:
            ax.plot(traj.x[-1, i, 0], traj.x[-1, i, 1], "o", ms=3, color=line.get_color())
    ax.set_aspect("equal", adjustable="datalim")
    return ax


def plot_energy_error(ax, traj, label: str | None = None, absolute: bool = True, **kw):
    """Plot the relative energy error of a trajectory versus time."""
    dE = relative_energy_error(traj.energy())
    name = traj.meta.get("integrator", "")
    style = integrator_style(name) if name else {}
    style.update(kw)
    y = np.abs(dE) if absolute else dE
    ax.plot(traj.t, y, label=label or traj.meta.get("integrator_label", name), **style)
    if absolute:
        ax.set_yscale("log")
    ax.set_xlabel("time")
    ax.set_ylabel(r"$|\Delta E / E_0|$" if absolute else r"$\Delta E / E_0$")
    return ax


def quicklook(traj, names=None, colors=None, title: str = ""):
    """Orbits and energy error side by side (used by ``nbody simulate --plot``)."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
    plot_orbits(ax1, traj, names, colors)
    ax1.set_title("trajectories")
    if names is not None and len(names) <= 12:
        ax1.legend(loc="upper right")
    plot_energy_error(ax2, traj)
    ax2.set_title("relative energy error")
    fig.suptitle(title)
    fig.tight_layout()
    return fig

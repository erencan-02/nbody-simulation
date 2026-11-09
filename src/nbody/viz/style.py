"""Matplotlib style shared by all figures (matches the thesis layout)."""

from __future__ import annotations

import matplotlib as mpl

# TU Darmstadt corporate design colours
TUDA = {
    "darkblue": "#243572",  # 1d - thesis accent colour
    "blue": "#005AA9",      # 1b
    "lightblue": "#0083CC",  # 2b
    "cyan": "#009CDA",      # 2a
    "green": "#009D81",     # 3b
    "lime": "#99C000",      # 4b
    "yellow": "#FDCA00",    # 6b
    "orange": "#EC6500",    # 8b
    "amber": "#F5A300",     # 7b
    "red": "#E6001A",       # 9b
    "magenta": "#A60084",   # 10b
    "purple": "#721085",    # 11b
    "grey": "#535353",      # 0d
    "lightgrey": "#B5B5B5",  # 0b
}

# One fixed colour and line style per integrator.  Symplectic methods are drawn
# with solid lines, non-symplectic ones dashed, so the structural property can
# be read off every figure.
INTEGRATOR_STYLE = {
    "euler": dict(color=TUDA["red"], ls="--", marker="o"),
    "symplectic_euler": dict(color=TUDA["orange"], ls="-", marker="o"),
    "heun": dict(color=TUDA["amber"], ls="--", marker="s"),
    "leapfrog": dict(color=TUDA["lightblue"], ls="-", marker="s"),
    "rk4": dict(color=TUDA["green"], ls="--", marker="^"),
    "yoshida4": dict(color=TUDA["darkblue"], ls="-", marker="^"),
    "yoshida6a": dict(color=TUDA["magenta"], ls="-", marker="D"),
    "yoshida6": dict(color=TUDA["purple"], ls="-", marker="v"),
    "dopri5": dict(color=TUDA["grey"], ls=":", marker="x"),
}

TEXT_WIDTH_IN = 147.0 / 25.4  # A4 minus 2 x 31.5 mm margins of the template


def integrator_style(name: str, with_marker: bool = False) -> dict:
    style = dict(INTEGRATOR_STYLE.get(name, dict(color="black", ls="-", marker=".")))
    if not with_marker:
        style.pop("marker")
    return style


def apply_style() -> None:
    mpl.rcParams.update(
        {
            "figure.figsize": (TEXT_WIDTH_IN, TEXT_WIDTH_IN * 0.55),
            "figure.dpi": 150,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.02,
            "font.size": 8.5,
            "axes.titlesize": 9,
            "axes.labelsize": 8.5,
            "legend.fontsize": 7.5,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "axes.grid": True,
            "grid.alpha": 0.3,
            "grid.linewidth": 0.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.prop_cycle": mpl.cycler(
                color=[TUDA[k] for k in ("darkblue", "orange", "green", "red",
                                         "lightblue", "purple", "amber", "grey")]
            ),
            "lines.linewidth": 1.2,
            "lines.markersize": 3.5,
            "legend.frameon": False,
            "mathtext.fontset": "cm",
            "svg.fonttype": "path",
        }
    )

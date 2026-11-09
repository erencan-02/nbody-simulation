"""E9 - Screenshots of the interactive viewer (rendered headless).

Shows the trajectory overlays: fading trails (past) and the dotted forecast
obtained by integrating a ghost copy of the system (future).
"""

from __future__ import annotations

import os

from nbody import load_scenario

from .common import FIG_DIR, log, setup


def run() -> dict:
    setup()
    os.environ["SDL_VIDEODRIVER"] = "dummy"
    os.environ.setdefault("PYGAME_HIDE_SUPPORT_PROMPT", "1")
    from nbody.viz.viewer import Viewer

    shots = {
        "viewer_inner_solar_system": dict(
            scenario=load_scenario("inner_solar_system"), integrator="leapfrog",
            dt=1 / 365.25, steps_per_frame=3, frames=420, horizon=1.0, trail=100,
            zoom=2.4),
        "viewer_circumbinary": dict(
            scenario=load_scenario("circumbinary"), integrator="yoshida4",
            dt=0.01, steps_per_frame=4, frames=500, horizon=12.0, trail=110, zoom=1.0),
    }
    out = {}
    for name, cfg in shots.items():
        v = Viewer(cfg["scenario"], integrator=cfg["integrator"], dt=cfg["dt"],
                   steps_per_frame=cfg["steps_per_frame"], size=(1000, 760),
                   trail_length=cfg["trail"], prediction_horizon=cfg["horizon"])
        v.camera.scale *= cfg["zoom"]
        path = FIG_DIR / f"{name}.png"
        v.run(max_frames=cfg["frames"], screenshot=path)
        out[name] = str(path.name)
        log(f"saved {path.name}")
    return out


if __name__ == "__main__":
    run()

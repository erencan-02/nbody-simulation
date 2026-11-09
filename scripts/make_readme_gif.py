"""Render the animated README showcases (headless) with the real-time viewer.

    python scripts/make_readme_gif.py                 # all showcases
    python scripts/make_readme_gif.py sun_earth_moon  # a single one

Each showcase advances the simulation by a fixed amount of simulated time per
frame, chosen so that the motion is clearly visible: e.g. Mercury completes an
orbit in about one second of animation, the stylised Moon in about 0.8 s.
Trails (past) and the dotted forecast overlay (future) are switched on.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ["SDL_VIDEODRIVER"] = "dummy"
os.environ.setdefault("PYGAME_HIDE_SUPPORT_PROMPT", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import pygame  # noqa: E402
from PIL import Image  # noqa: E402

from nbody import load_scenario  # noqa: E402
from nbody.viz.viewer import Viewer  # noqa: E402

SIZE = (720, 540)
FPS = 25

SHOWCASES = {
    # inner Solar System, leapfrog with h = 1 day, 4 days per frame
    "showcase": dict(scenario="inner_solar_system", integrator="leapfrog", dt=1 / 365.25,
                     steps_per_frame=4, trail=70, horizon=0.8, zoom=2.3, seconds=9, warmup=30),
    # star - planet - moon, Yoshida 4, ~6 s per planetary orbit, ~0.8 s per lunar orbit
    "sun_earth_moon": dict(scenario="sun_earth_moon", integrator="yoshida4", dt=None,
                           steps_per_frame=10, trail=110, horizon=1.2, zoom=1.0, seconds=10,
                           warmup=60, inset=dict(body=1, zoom=4.0, size=(250, 190), trail=22,
                                                 label="around the Earth, x4")),
}


def draw_inset(viewer: Viewer, screen, font, cfg: dict) -> None:
    """Magnified view centred on one body, pasted into the lower right corner."""
    cam = viewer.camera
    saved = (cam.center.copy(), cam.scale, viewer.follow_com, viewer.show_trails,
             viewer.show_prediction)
    viewer.follow_com = viewer.show_trails = viewer.show_prediction = False
    centre = viewer.sim.state.x[cfg["body"], :2].copy()
    cam.center = centre
    cam.scale = saved[1] * cfg["zoom"]
    full = pygame.Surface(SIZE)
    viewer.render(full, pygame, font)
    # trails of the other bodies relative to the centred body (the Moon's orbit
    # around the Earth instead of its inertial epicycles)
    pts = viewer.trails.ordered()[-cfg.get("trail", 10**9):]
    if len(pts) >= 2:
        for i in range(pts.shape[1]):
            if i == cfg["body"]:
                continue
            rel = pts[:, i, :2] - pts[:, cfg["body"], :2] + centre
            scr = cam.to_screen(rel)
            n = len(scr)
            for c in range(6):
                a, b = c * n // 6, min(n, (c + 1) * n // 6 + 1)
                if b - a >= 2:
                    f = 0.15 + 0.75 * (c + 1) / 6
                    col = tuple(int(bg + (cc - bg) * f) for bg, cc in zip((11, 16, 38), viewer.colors[i]))
                    pygame.draw.aalines(full, col, False, scr[a:b].tolist())
        for i in range(pts.shape[1]):  # redraw bodies on top of the trails
            p = cam.to_screen(viewer.sim.state.x[i, :2])
            pygame.draw.circle(full, viewer.colors[i], (int(p[0]), int(p[1])), int(viewer.radii[i]))
    cam.center, cam.scale, viewer.follow_com, viewer.show_trails, viewer.show_prediction = saved
    w, h = cfg["size"]
    crop = full.subsurface(pygame.Rect(SIZE[0] // 2 - w // 2, SIZE[1] // 2 - h // 2, w, h))
    x0, y0 = SIZE[0] - w - 12, SIZE[1] - h - 12
    screen.blit(crop, (x0, y0))
    pygame.draw.rect(screen, (93, 133, 195), pygame.Rect(x0 - 1, y0 - 1, w + 2, h + 2), 1)
    screen.blit(font.render(cfg["label"], True, (180, 190, 215)), (x0 + 6, y0 + 4))


def render(name: str, cfg: dict, out: Path) -> None:
    sc = load_scenario(cfg["scenario"])
    viewer = Viewer(sc, integrator=cfg["integrator"], dt=cfg["dt"],
                    steps_per_frame=cfg["steps_per_frame"], size=SIZE,
                    trail_length=cfg["trail"], prediction_horizon=cfg["horizon"])
    viewer.camera.scale *= cfg["zoom"]
    pygame.init()
    screen = pygame.display.set_mode(SIZE)
    font = pygame.font.SysFont("consolas,dejavusansmono,monospace", 13)
    frames = []
    for k in range(cfg["warmup"] + FPS * cfg["seconds"]):
        viewer.advance()
        viewer.render(screen, pygame, font)
        if "inset" in cfg:
            draw_inset(viewer, screen, font, cfg["inset"])
        viewer.frame += 1
        if k >= cfg["warmup"]:
            frames.append(Image.frombytes("RGB", SIZE, pygame.image.tobytes(screen, "RGB")))
    pygame.quit()

    # one shared palette (the scene uses few colours) keeps the GIF small and flicker-free
    palette = frames[len(frames) // 2].quantize(colors=128, method=Image.Quantize.MEDIANCUT)
    gif = [f.quantize(palette=palette, dither=Image.Dither.NONE) for f in frames]
    gif[0].save(out, save_all=True, append_images=gif[1:], duration=int(1000 / FPS),
                loop=0, optimize=True, disposal=1)
    print(f"wrote {out} ({len(gif)} frames, {out.stat().st_size / 1e6:.1f} MB)")


def main(names: list[str]) -> None:
    out_dir = ROOT / "docs"
    out_dir.mkdir(exist_ok=True)
    for name in names or list(SHOWCASES):
        render(name, SHOWCASES[name], out_dir / f"{name}.gif")


if __name__ == "__main__":
    main(sys.argv[1:])

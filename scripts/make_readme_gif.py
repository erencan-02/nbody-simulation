"""Render the animated README showcase (headless) with the real-time viewer.

    python scripts/make_readme_gif.py [output.gif]

The inner Solar System is integrated with leapfrog (h = 1 day); every frame
advances the simulation by `DAYS_PER_FRAME` days, so Mercury completes an
orbit in about one second of animation and the Earth in about four.
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
SECONDS = 9
DAYS_PER_FRAME = 4
WARMUP_FRAMES = 30  # let trails build up before recording


def main(out: Path) -> None:
    viewer = Viewer(
        load_scenario("inner_solar_system"),
        integrator="leapfrog",
        dt=1 / 365.25,
        steps_per_frame=DAYS_PER_FRAME,
        size=SIZE,
        trail_length=70,
        prediction_horizon=0.8,
    )
    viewer.camera.scale *= 2.3
    pygame.init()
    screen = pygame.display.set_mode(SIZE)
    font = pygame.font.SysFont("consolas,dejavusansmono,monospace", 13)
    frames = []
    for k in range(WARMUP_FRAMES + FPS * SECONDS):
        viewer.advance()
        viewer.render(screen, pygame, font)
        viewer.frame += 1
        if k >= WARMUP_FRAMES:
            raw = pygame.image.tobytes(screen, "RGB")
            frames.append(Image.frombytes("RGB", SIZE, raw))
    pygame.quit()

    # one shared palette (the scene uses few colours) keeps the GIF small and flicker-free
    palette = frames[len(frames) // 2].quantize(colors=128, method=Image.Quantize.MEDIANCUT)
    gif = [f.quantize(palette=palette, dither=Image.Dither.NONE) for f in frames]
    gif[0].save(out, save_all=True, append_images=gif[1:], duration=int(1000 / FPS),
                loop=0, optimize=True, disposal=1)
    print(f"wrote {out} ({len(gif)} frames, {out.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "docs" / "showcase.gif")

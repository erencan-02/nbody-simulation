"""Interactive real-time viewer (pygame) with trajectory overlays.

Overlays
--------
* **trails** - fading polylines of the recent past positions;
* **prediction** - dotted forecast of the future trajectories obtained by
  integrating a clone of the running system (:class:`TrajectoryPredictor`);
* **vectors** - velocity arrows;
* **HUD** - integrator, step size, time and live conservation errors.

Controls
--------
``SPACE`` pause, ``T`` trails, ``P`` prediction, ``V`` velocity vectors,
``C`` follow centre of mass, ``F`` fit view, ``+``/``-``/wheel zoom,
arrow keys pan, ``.``/``,`` faster/slower, ``1``-``9`` switch integrator,
``R`` reset, ``S`` screenshot, ``H`` help, ``ESC``/``Q`` quit.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..analysis import diagnostics as diag
from ..integrators import available_integrators, create_integrator
from ..scenarios import Scenario
from ..simulation import Simulation
from .overlay import TrailBuffer, TrajectoryPredictor

BACKGROUND = (11, 16, 38)
GRID = (32, 42, 78)
HUD_TEXT = (220, 226, 240)
HUD_DIM = (140, 150, 180)
ACCENT = (93, 133, 195)


def _hex_to_rgb(h: str) -> tuple[int, int, int]:
    h = h.lstrip("#")
    return tuple(int(h[i : i + 2], 16) for i in (0, 2, 4))  # type: ignore[return-value]


def _blend(c: tuple[int, int, int], f: float) -> tuple[int, int, int]:
    """Blend colour ``c`` with the background (f = 1: pure colour)."""
    return tuple(int(BACKGROUND[k] + (c[k] - BACKGROUND[k]) * f) for k in range(3))  # type: ignore


@dataclass
class Camera:
    """World <-> screen transform (y axis pointing up)."""

    width: int
    height: int
    center: np.ndarray
    scale: float  # pixels per world length unit

    def to_screen(self, p: np.ndarray) -> np.ndarray:
        p = np.asarray(p)[..., :2]
        sx = self.width / 2 + (p[..., 0] - self.center[0]) * self.scale
        sy = self.height / 2 - (p[..., 1] - self.center[1]) * self.scale
        return np.stack([sx, sy], axis=-1)

    def fit(self, x: np.ndarray, center: np.ndarray | None = None, margin: float = 1.15) -> None:
        """Zoom so that all points ``x`` are visible around ``center``."""
        x = np.asarray(x)[..., :2].reshape(-1, 2)
        if center is None:
            center = 0.5 * (x.min(axis=0) + x.max(axis=0))
        self.center = np.array(center, dtype=float)[:2]
        extent = max(float(np.max(np.linalg.norm(x - self.center, axis=-1))), 1e-12) * margin
        self.scale = min(self.width, self.height) / 2.0 / extent


class Viewer:
    """Real-time visualisation of a :class:`Scenario`."""

    def __init__(
        self,
        scenario: Scenario,
        integrator: str = "leapfrog",
        dt: float | None = None,
        steps_per_frame: int = 1,
        size: tuple[int, int] = (1100, 850),
        trail_length: int = 400,
        prediction_horizon: float | None = None,
        fps: int = 60,
    ) -> None:
        try:
            import pygame  # noqa: F401
        except ImportError as exc:  # pragma: no cover - depends on the environment
            raise ImportError("the viewer needs pygame: pip install pygame") from exc

        self.scenario = scenario
        self.integrator_name = integrator
        self.dt = float(dt if dt is not None else scenario.recommended.get("dt", 1e-2))
        self.steps_per_frame = max(1, int(steps_per_frame))
        self.size = size
        self.fps = fps
        self.trail_length = trail_length
        horizon = prediction_horizon
        if horizon is None:
            t_end = scenario.recommended.get("t_end", 200 * self.dt)
            horizon = min(t_end, 300 * self.dt * self.steps_per_frame)
        self.predictor = TrajectoryPredictor(horizon, self.dt * 2, n_points=300)

        self.show_trails = True
        self.show_prediction = True
        self.show_vectors = False
        self.show_help = False
        self.follow_com = True
        self.paused = False
        self._screenshot_requested = False
        self.colors = [
            _hex_to_rgb(c) for c in (scenario.colors or ["#FFFFFF"] * scenario.n_bodies)
        ]
        m = scenario.masses
        rel = (m / m.max()) ** (1.0 / 3.0)
        self.radii = np.clip(2.0 + 9.0 * rel, 2.0, 11.0)
        self.reset()

    # -- state management --------------------------------------------------------
    def reset(self) -> None:
        self.sim = Simulation(
            self.scenario.model(), self.scenario.initial_state(), self.integrator_name
        )
        self.E0 = self.sim.energy()
        st = self.sim.state
        self.L0 = diag.angular_momentum(st.x, st.v, self.scenario.masses)
        self.trails = TrailBuffer(self.scenario.n_bodies, self.scenario.dim, self.trail_length)
        self.camera = Camera(self.size[0], self.size[1], np.zeros(2), 1.0)
        self.camera.fit(st.x, diag.center_of_mass(st.x, self.scenario.masses))
        self.prediction: np.ndarray | None = None
        self.frame = 0

    def set_integrator(self, name: str) -> None:
        self.integrator_name = name
        self.sim.set_integrator(create_integrator(name))
        self.prediction = None

    def advance(self) -> None:
        for _ in range(self.steps_per_frame):
            self.sim.step(self.dt)
        self.trails.push(self.sim.state.x)
        if self.show_prediction and (self.prediction is None or self.frame % 10 == 0):
            self.prediction = self.predictor.predict(self.sim)

    # -- drawing ---------------------------------------------------------------------
    def _draw_grid(self, screen, pygame) -> None:
        cam = self.camera
        span = self.size[0] / cam.scale
        step = 10.0 ** np.floor(np.log10(span / 4))
        x0 = cam.center[0] - span / 2
        y_span = self.size[1] / cam.scale
        for k in range(int(np.floor(x0 / step)), int(np.ceil((x0 + span) / step)) + 1):
            sx = cam.to_screen(np.array([k * step, 0.0]))[0]
            pygame.draw.line(screen, GRID, (sx, 0), (sx, self.size[1]), 1)
        y0 = cam.center[1] - y_span / 2
        for k in range(int(np.floor(y0 / step)), int(np.ceil((y0 + y_span) / step)) + 1):
            sy = cam.to_screen(np.array([0.0, k * step]))[1]
            pygame.draw.line(screen, GRID, (0, sy), (self.size[0], sy), 1)

    def _draw_trails(self, screen, pygame) -> None:
        pts = self.trails.ordered()
        if len(pts) < 2:
            return
        scr = self.camera.to_screen(pts)  # (L, N, 2)
        n_chunks = 6
        bounds = np.linspace(0, len(pts) - 1, n_chunks + 1).astype(int)
        for i in range(self.scenario.n_bodies):
            for c in range(n_chunks):
                a, b = bounds[c], bounds[c + 1] + 1
                if b - a < 2:
                    continue
                col = _blend(self.colors[i], 0.15 + 0.75 * (c + 1) / n_chunks)
                pygame.draw.aalines(screen, col, False, scr[a:b, i].tolist())

    def _draw_prediction(self, screen, pygame) -> None:
        if self.prediction is None or len(self.prediction) < 2:
            return
        scr = self.camera.to_screen(self.prediction)
        dash = 3
        for i in range(self.scenario.n_bodies):
            col = _blend(self.colors[i], 0.55)
            path = scr[:, i]
            for s in range(0, len(path) - 1, 2 * dash):
                seg = path[s : s + dash + 1]
                if len(seg) >= 2:
                    pygame.draw.aalines(screen, col, False, seg.tolist())
            end = path[-1]
            pygame.draw.circle(screen, col, (int(end[0]), int(end[1])), 4, 1)

    def _draw_bodies(self, screen, pygame, font) -> None:
        st = self.sim.state
        scr = self.camera.to_screen(st.x)
        for i in range(self.scenario.n_bodies):
            p = (int(scr[i, 0]), int(scr[i, 1]))
            r = int(self.radii[i])
            if self.show_vectors:
                vend = self.camera.to_screen(st.x[i] + st.v[i] * self.predictor.horizon * 0.05)
                pygame.draw.line(screen, (230, 0, 26), p, (int(vend[0]), int(vend[1])), 2)
            pygame.draw.circle(screen, self.colors[i], p, r)
            if self.scenario.n_bodies <= 12:
                label = font.render(self.scenario.names[i], True, HUD_DIM)
                screen.blit(label, (p[0] + r + 3, p[1] - r - 3))

    def _draw_hud(self, screen, font) -> None:
        st = self.sim.state
        integ = self.sim.integrator
        E = self.sim.energy()
        dE = (E - self.E0) / abs(self.E0)
        L = diag.angular_momentum(st.x, st.v, self.scenario.masses)
        nL = float(np.linalg.norm(self.L0))
        dL = float(np.linalg.norm(L - self.L0)) / nL if nL > 0 else float(np.linalg.norm(L))
        tunit = self.scenario.units.get("time", "")
        tu = f" {tunit}" if tunit else ""
        lines = [
            (f"{self.scenario.name}   N = {self.scenario.n_bodies}", HUD_TEXT),
            (f"integrator  {integ.label}  (order {integ.order}, "
             f"{'symplectic' if integ.symplectic else 'non-symplectic'})", HUD_TEXT),
            (f"dt = {self.dt:.3g}{tu}   x{self.steps_per_frame} steps/frame", HUD_DIM),
            (f"t  = {st.t:.4g}{tu}   force evals = {self.sim.model.n_evals}", HUD_DIM),
            (f"dE/E0 = {dE:+.3e}    |dL|/|L0| = {dL:.2e}", HUD_TEXT),
            (f"trails [T] {'on' if self.show_trails else 'off'}   "
             f"prediction [P] {'on' if self.show_prediction else 'off'} "
             f"({self.predictor.horizon:.3g}{tu})   help [H]", HUD_DIM),
        ]
        if self.paused:
            lines.append(("PAUSED", ACCENT))
        y = 10
        for text, col in lines:
            screen.blit(font.render(text, True, col), (12, y))
            y += 18
        if self.show_help:
            names = available_integrators()
            help_lines = [
                "SPACE pause   T trails   P prediction   V velocities",
                "C follow centre of mass   F fit view   R reset   S screenshot",
                "+/- or wheel zoom   arrows pan   . / , faster / slower",
            ] + [f"{k + 1}: {n}" for k, n in enumerate(names[:9])]
            y = self.size[1] - 18 * len(help_lines) - 10
            for text in help_lines:
                screen.blit(font.render(text, True, HUD_DIM), (12, y))
                y += 18

    def render(self, screen, pygame, font) -> None:
        if self.follow_com:
            self.camera.center = diag.center_of_mass(self.sim.state.x, self.scenario.masses)[:2]
        screen.fill(BACKGROUND)
        self._draw_grid(screen, pygame)
        if self.show_trails:
            self._draw_trails(screen, pygame)
        if self.show_prediction:
            self._draw_prediction(screen, pygame)
        self._draw_bodies(screen, pygame, font)
        self._draw_hud(screen, font)

    # -- event loop ----------------------------------------------------------------------
    def _handle_key(self, key, pygame) -> bool:
        names = available_integrators()
        if key in (pygame.K_ESCAPE, pygame.K_q):
            return False
        if key == pygame.K_SPACE:
            self.paused = not self.paused
        elif key == pygame.K_t:
            self.show_trails = not self.show_trails
        elif key == pygame.K_p:
            self.show_prediction = not self.show_prediction
            self.prediction = None
        elif key == pygame.K_v:
            self.show_vectors = not self.show_vectors
        elif key == pygame.K_c:
            self.follow_com = not self.follow_com
        elif key == pygame.K_h:
            self.show_help = not self.show_help
        elif key == pygame.K_f:
            st = self.sim.state
            self.camera.fit(st.x, diag.center_of_mass(st.x, self.scenario.masses))
        elif key == pygame.K_r:
            self.reset()
        elif key == pygame.K_s:
            self._screenshot_requested = True
        elif key in (pygame.K_PLUS, pygame.K_EQUALS, pygame.K_KP_PLUS):
            self.camera.scale *= 1.25
        elif key in (pygame.K_MINUS, pygame.K_KP_MINUS):
            self.camera.scale /= 1.25
        elif key == pygame.K_PERIOD:
            self.steps_per_frame = min(self.steps_per_frame * 2, 4096)
        elif key == pygame.K_COMMA:
            self.steps_per_frame = max(self.steps_per_frame // 2, 1)
        elif pygame.K_1 <= key <= pygame.K_9:
            idx = key - pygame.K_1
            if idx < len(names):
                self.set_integrator(names[idx])
        return True

    def run(self, max_frames: int | None = None, screenshot: str | Path | None = None) -> None:
        """Open the window and run the event loop.

        With ``max_frames`` the loop stops after that many frames; together
        with ``screenshot`` this renders an image without user interaction
        (set ``SDL_VIDEODRIVER=dummy`` for a fully headless run).
        """
        import pygame

        os.environ.setdefault("PYGAME_HIDE_SUPPORT_PROMPT", "1")
        pygame.init()
        try:
            screen = pygame.display.set_mode(self.size)
            pygame.display.set_caption(f"N-body viewer - {self.scenario.name}")
            font = pygame.font.SysFont("consolas,dejavusansmono,monospace", 14)
            clock = pygame.time.Clock()
            running = True
            while running:
                for event in pygame.event.get():
                    if event.type == pygame.QUIT:
                        running = False
                    elif event.type == pygame.KEYDOWN:
                        running = self._handle_key(event.key, pygame)
                    elif event.type == pygame.MOUSEWHEEL:
                        self.camera.scale *= 1.15 ** event.y
                keys = pygame.key.get_pressed()
                pan = 8.0 / self.camera.scale
                if keys[pygame.K_LEFT]:
                    self.camera.center[0] -= pan
                    self.follow_com = False
                if keys[pygame.K_RIGHT]:
                    self.camera.center[0] += pan
                    self.follow_com = False
                if keys[pygame.K_UP]:
                    self.camera.center[1] += pan
                    self.follow_com = False
                if keys[pygame.K_DOWN]:
                    self.camera.center[1] -= pan
                    self.follow_com = False
                if not self.paused:
                    self.advance()
                self.render(screen, pygame, font)
                pygame.display.flip()
                self.frame += 1
                if self._screenshot_requested:
                    self._screenshot_requested = False
                    pygame.image.save(screen, f"nbody_{self.scenario.name}_{self.frame}.png")
                if max_frames is not None and self.frame >= max_frames:
                    running = False
                if max_frames is None:
                    clock.tick(self.fps)
            if screenshot is not None:
                pygame.image.save(screen, str(screenshot))
        finally:
            pygame.quit()

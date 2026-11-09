"""Simulation driver and trajectory container."""

from __future__ import annotations

import json
import math
import time
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .analysis import diagnostics as diag
from .core.model import ForceModel, GravityModel
from .core.state import State
from .integrators import Integrator, create_integrator


@dataclass
class Trajectory:
    """Recorded time series of a simulation.

    Attributes
    ----------
    t:
        Sample times, shape ``(T,)``.
    x, v:
        Positions and velocities, shape ``(T, N, d)``.
    masses:
        Body masses, shape ``(N,)``.
    G, softening:
        Parameters of the gravity model that produced the data.
    step_sizes:
        All accepted step sizes (only recorded for adaptive integrators).
    meta:
        Free-form, JSON-serialisable metadata (integrator, cost, ...).
    """

    t: np.ndarray
    x: np.ndarray
    v: np.ndarray
    masses: np.ndarray
    G: float = 1.0
    softening: float = 0.0
    step_sizes: np.ndarray | None = None
    meta: dict = field(default_factory=dict)
    model: ForceModel | None = field(default=None, repr=False, compare=False)

    # -- shape helpers ---------------------------------------------------------
    @property
    def n_frames(self) -> int:
        return self.t.size

    @property
    def n_bodies(self) -> int:
        return self.masses.size

    @property
    def dim(self) -> int:
        return self.x.shape[-1]

    # -- diagnostics -------------------------------------------------------------
    def kinetic_energy(self) -> np.ndarray:
        return diag.kinetic_energy(self.v, self.masses)

    def potential_energy(self) -> np.ndarray:
        """Potential energy of the force model that produced the trajectory.

        Falls back to (softened) gravity if no model is attached, e.g. for
        trajectories loaded from disk.
        """
        if self.model is not None:
            return np.asarray(self.model.potential_energy(self.x))
        return diag.potential_energy(self.x, self.masses, self.G, self.softening)

    def energy(self) -> np.ndarray:
        return self.kinetic_energy() + self.potential_energy()

    def linear_momentum(self) -> np.ndarray:
        return diag.linear_momentum(self.v, self.masses)

    def angular_momentum(self) -> np.ndarray:
        return diag.angular_momentum(self.x, self.v, self.masses)

    def center_of_mass(self) -> np.ndarray:
        return diag.center_of_mass(self.x, self.masses)

    def state(self, k: int = -1) -> State:
        return State(self.x[k], self.v[k], float(self.t[k]))

    # -- persistence -------------------------------------------------------------
    def save(self, path: str | Path) -> None:
        """Store as compressed ``.npz`` (no pickling, safe to share)."""
        arrays = dict(t=self.t, x=self.x, v=self.v, masses=self.masses)
        if self.step_sizes is not None:
            arrays["step_sizes"] = self.step_sizes
        np.savez_compressed(
            path,
            **arrays,
            G=np.float64(self.G),
            softening=np.float64(self.softening),
            meta=np.array(json.dumps(self.meta)),
        )

    @classmethod
    def load(cls, path: str | Path) -> "Trajectory":
        with np.load(path, allow_pickle=False) as f:
            return cls(
                t=f["t"],
                x=f["x"],
                v=f["v"],
                masses=f["masses"],
                G=float(f["G"]),
                softening=float(f["softening"]),
                step_sizes=f["step_sizes"] if "step_sizes" in f.files else None,
                meta=json.loads(str(f["meta"])),
            )


class Simulation:
    """Couples a force model, a state and an integrator.

    The object can be advanced step by step (used by the interactive viewer)
    or run to a final time while recording a :class:`Trajectory`.
    """

    def __init__(
        self, model: ForceModel, state: State, integrator: Integrator | str = "leapfrog"
    ) -> None:
        if state.n_bodies != model.n_bodies:
            raise ValueError("state and model have a different number of bodies")
        self.model = model
        self.state = state.copy()
        self.integrator = (
            create_integrator(integrator) if isinstance(integrator, str) else integrator
        )
        self.n_steps = 0
        self.n_rejected = 0
        self._dt_next: float | None = None

    # -- stepping ----------------------------------------------------------------
    def step(self, dt: float) -> float:
        """Perform one step; returns the step size that was actually taken."""
        if self.integrator.adaptive:
            # ``dt`` acts as an upper bound; the controller's proposal is used if smaller.
            if self._dt_next is not None:
                dt = min(dt, self._dt_next)
            res = self.integrator.step(self.state, dt, self.model)
            self._dt_next = res.dt_next
            self.n_rejected += res.n_rejected
            taken = res.dt_taken
        else:
            self.integrator.step(self.state, dt, self.model)
            taken = dt
        self.n_steps += 1
        return taken

    def clone(self) -> "Simulation":
        """Independent copy (model counters are reset, the integrator is shared).

        Integrators are stateless, so sharing them is safe.
        """
        sim = Simulation(self.model.copy(), self.state, self.integrator)
        sim._dt_next = self._dt_next
        return sim

    def set_integrator(self, integrator: Integrator | str) -> None:
        self.integrator = (
            create_integrator(integrator) if isinstance(integrator, str) else integrator
        )
        self._dt_next = None

    def energy(self) -> float:
        return self.model.energy(self.state)

    # -- batch runs ----------------------------------------------------------------
    def run(
        self,
        t_end: float,
        dt: float,
        record_every: int = 1,
        callback: Callable[[State], None] | None = None,
    ) -> Trajectory:
        """Integrate from the current time to ``t_end``.

        Fixed-step methods take ``n = ceil((t_end - t0) / dt)`` steps of size
        ``(t_end - t0) / n <= dt`` so that ``t_end`` is hit exactly.  Adaptive
        methods use ``dt`` as the initial step and clip the last step.
        Every ``record_every``-th step and always the final state are stored.
        """
        if dt <= 0:
            raise ValueError("dt must be positive")
        if record_every < 1:
            raise ValueError("record_every must be >= 1")
        t0 = self.state.t
        if t_end <= t0:
            raise ValueError("t_end must be larger than the current time")
        evals0 = self.model.n_evals
        wall0 = time.perf_counter()
        if self.integrator.adaptive:
            traj = self._run_adaptive(t_end, dt, record_every, callback)
        else:
            traj = self._run_fixed(t_end, dt, record_every, callback)
        traj.meta.update(
            integrator=self.integrator.name,
            integrator_label=self.integrator.label,
            n_force_evals=int(self.model.n_evals - evals0),
            wall_time=time.perf_counter() - wall0,
        )
        return traj

    def _new_trajectory(self, t, x, v, step_sizes=None, **meta) -> Trajectory:
        G = getattr(self.model, "G", 1.0)
        eps = getattr(self.model, "softening", 0.0)
        if not isinstance(self.model, GravityModel):
            meta["model"] = type(self.model).__name__
        return Trajectory(
            t, x, v, self.model.masses.copy(), G, eps, step_sizes, meta, model=self.model.copy()
        )

    def _run_fixed(self, t_end, dt, record_every, callback) -> Trajectory:
        st = self.state
        t0 = st.t
        n_steps = max(1, math.ceil((t_end - t0) / dt - 1e-9))
        h = (t_end - t0) / n_steps
        n_rec = n_steps // record_every + 1 + (1 if n_steps % record_every else 0)
        ts = np.empty(n_rec)
        xs = np.empty((n_rec,) + st.x.shape)
        vs = np.empty((n_rec,) + st.v.shape)
        ts[0], xs[0], vs[0] = t0, st.x, st.v
        r = 1
        step = self.integrator.step
        model = self.model
        for k in range(1, n_steps + 1):
            step(st, h, model)
            st.t = t0 + k * h  # exact time, no accumulated round-off
            if callback is not None:
                callback(st)
            if k % record_every == 0 or k == n_steps:
                ts[r], xs[r], vs[r] = st.t, st.x, st.v
                r += 1
        self.n_steps += n_steps
        assert r == n_rec
        return self._new_trajectory(ts, xs, vs, dt=h, n_steps=n_steps)

    def _run_adaptive(self, t_end, dt, record_every, callback) -> Trajectory:
        st = self.state
        ts, xs, vs = [st.t], [st.x.copy()], [st.v.copy()]
        steps: list[float] = []
        dt_next = dt if self._dt_next is None else self._dt_next
        rejected0 = self.n_rejected
        k = 0
        tol = 1e-12 * max(1.0, abs(t_end))
        while t_end - st.t > tol:
            h = min(dt_next, t_end - st.t)
            res = self.integrator.step(st, h, self.model)
            self.n_rejected += res.n_rejected
            dt_next = res.dt_next
            steps.append(res.dt_taken)
            k += 1
            if t_end - st.t <= tol:
                st.t = t_end
            if callback is not None:
                callback(st)
            if k % record_every == 0 or st.t == t_end:
                ts.append(st.t)
                xs.append(st.x.copy())
                vs.append(st.v.copy())
        self._dt_next = dt_next
        self.n_steps += k
        return self._new_trajectory(
            np.array(ts),
            np.array(xs),
            np.array(vs),
            step_sizes=np.array(steps),
            n_steps=k,
            n_rejected=int(self.n_rejected - rejected0),
            rtol=getattr(self.integrator, "rtol", None),
            atol=getattr(self.integrator, "atol", None),
        )


def simulate(
    model: ForceModel,
    state: State,
    integrator: Integrator | str,
    t_end: float,
    dt: float,
    record_every: int = 1,
) -> Trajectory:
    """Convenience wrapper: run a fresh :class:`Simulation` and return the trajectory."""
    return Simulation(model, state, integrator).run(t_end, dt, record_every)

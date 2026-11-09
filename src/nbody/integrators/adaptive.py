r"""Adaptive Dormand-Prince 5(4) method with embedded error estimation.

Each step computes a fifth-order solution :math:`y_{n+1}` and an embedded
fourth-order solution :math:`\hat y_{n+1}` from the same seven stages.  Their
difference estimates the local error, which is measured in the scaled norm

.. math::  \mathrm{err} = \sqrt{\frac1n\sum_k
   \left(\frac{y_{n+1,k}-\hat y_{n+1,k}}
   {\mathrm{atol}+\mathrm{rtol}\max(|y_{n,k}|,|y_{n+1,k}|)}\right)^2}.

A step is accepted if :math:`\mathrm{err}\le 1`; in both cases the next step
size is :math:`h_\text{new} = h\,\min(f_\max,\max(f_\min, s\,\mathrm{err}^{-1/5}))`
(Hairer, Nørsett & Wanner 1993, Sec. II.4).  The last stage is evaluated at
the new solution (FSAL), so an accepted step costs six force evaluations.

Adaptivity resolves close encounters automatically, but a variable step size
destroys the symplectic structure: the method is neither symplectic nor
time-reversible.
"""

from __future__ import annotations

import numpy as np

from ..core.model import ForceModel
from ..core.state import State
from .base import AdaptiveIntegrator, AdaptiveStepResult, cached_acceleration

_C = np.array([0.0, 1 / 5, 3 / 10, 4 / 5, 8 / 9, 1.0, 1.0])
_A = np.array(
    [
        [0, 0, 0, 0, 0, 0, 0],
        [1 / 5, 0, 0, 0, 0, 0, 0],
        [3 / 40, 9 / 40, 0, 0, 0, 0, 0],
        [44 / 45, -56 / 15, 32 / 9, 0, 0, 0, 0],
        [19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729, 0, 0, 0],
        [9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656, 0, 0],
        [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84, 0],
    ],
    dtype=np.float64,
)
_B = _A[6].copy()  # fifth-order weights (= last row: FSAL)
_B_HAT = np.array(
    [5179 / 57600, 0, 7571 / 16695, 393 / 640, -92097 / 339200, 187 / 2100, 1 / 40]
)
_E = _B - _B_HAT


class DormandPrince54(AdaptiveIntegrator):
    """Embedded Runge-Kutta pair DOPRI5(4) with standard step-size control."""

    name = "dopri5"
    label = "Dormand-Prince 5(4)"
    order = 5
    symplectic = False
    time_reversible = False
    evals_per_step = 6

    def __init__(
        self,
        rtol: float = 1e-9,
        atol: float = 1e-12,
        safety: float = 0.9,
        fac_min: float = 0.2,
        fac_max: float = 5.0,
        dt_min: float = 1e-14,
        max_rejections: int = 50,
    ) -> None:
        if rtol <= 0 or atol < 0:
            raise ValueError("tolerances must be positive")
        self.rtol = float(rtol)
        self.atol = float(atol)
        self.safety = float(safety)
        self.fac_min = float(fac_min)
        self.fac_max = float(fac_max)
        self.dt_min = float(dt_min)
        self.max_rejections = int(max_rejections)

    def _attempt(self, state: State, dt: float, model: ForceModel):
        x0, v0 = state.x, state.v
        kx: list[np.ndarray] = [v0]
        kv: list[np.ndarray] = [cached_acceleration(state, model)]
        for i in range(1, 7):
            xi = x0.copy()
            vi = v0.copy()
            for j in range(i):
                aij = _A[i, j]
                if aij != 0.0:
                    xi += (dt * aij) * kx[j]
                    vi += (dt * aij) * kv[j]
            kx.append(vi)
            kv.append(model.acceleration(xi))
            if i == 6:
                x_new, v_new, a_new = xi, vi, kv[-1]
        ex = np.zeros_like(x0)
        ev = np.zeros_like(v0)
        for ei, kxi, kvi in zip(_E, kx, kv):
            if ei != 0.0:
                ex += (dt * ei) * kxi
                ev += (dt * ei) * kvi
        sx = self.atol + self.rtol * np.maximum(np.abs(x0), np.abs(x_new))
        sv = self.atol + self.rtol * np.maximum(np.abs(v0), np.abs(v_new))
        err = np.sqrt(
            (np.sum((ex / sx) ** 2) + np.sum((ev / sv) ** 2)) / (ex.size + ev.size)
        )
        return x_new, v_new, a_new, float(err)

    def step(self, state: State, dt: float, model: ForceModel) -> AdaptiveStepResult:
        if dt <= 0:
            raise ValueError("adaptive integration requires a positive step size")
        n_rejected = 0
        while True:
            x_new, v_new, a_new, err = self._attempt(state, dt, model)
            if not np.isfinite(err):
                err = np.inf
            if err <= 1.0:
                if err == 0.0:
                    factor = self.fac_max
                else:
                    factor = min(self.fac_max, max(self.fac_min, self.safety * err**-0.2))
                if n_rejected > 0:  # do not grow directly after a rejection
                    factor = min(factor, 1.0)
                state.x, state.v, state.a = x_new, v_new, a_new
                state.t += dt
                return AdaptiveStepResult(dt, dt * factor, n_rejected)
            n_rejected += 1
            if n_rejected > self.max_rejections or dt < self.dt_min:
                raise RuntimeError(
                    f"step size control failed at t={state.t:.6g} (dt={dt:.3g})"
                )
            factor = self.fac_min if not np.isfinite(err) else max(
                self.fac_min, self.safety * err**-0.2
            )
            dt *= factor

    def __repr__(self) -> str:
        return f"DormandPrince54(rtol={self.rtol:g}, atol={self.atol:g})"

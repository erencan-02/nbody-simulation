r"""Evaluation metrics for numerical N-body integrations.

The metrics fall into four groups:

1. **Conservation errors** - relative energy error
   :math:`\delta E(t) = (E(t)-E_0)/|E_0|`, its maximum, RMS and secular drift
   rate; momentum and angular-momentum errors.
2. **Accuracy against a reference** - global position / phase-space error
   against an exact or high-precision solution, and the empirical order of
   convergence obtained from a least-squares fit in log-log space.
3. **Structural properties** - time-reversibility error
   :math:`\|\Phi_{-h}^n\circ\Phi_h^n(y_0)-y_0\|`.
4. **Physical observables** - orbital periods (mean motion) and the
   finite-time Lyapunov exponent of chaotic configurations.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

from ..core.model import ForceModel
from ..core.state import State
from ..integrators import Integrator
from . import diagnostics as diag


# ---------------------------------------------------------------------------
# 1. Conservation
# ---------------------------------------------------------------------------
def relative_energy_error(E: np.ndarray) -> np.ndarray:
    r""":math:`\delta E(t) = (E(t) - E_0)/|E_0|`."""
    E = np.asarray(E, dtype=np.float64)
    return (E - E[0]) / abs(E[0])


def linear_drift_rate(t: np.ndarray, y: np.ndarray) -> float:
    """Slope of the least-squares line through ``y(t)`` (secular drift per unit time)."""
    t = np.asarray(t, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if t.size < 2:
        return 0.0
    tc = t - t.mean()
    return float(np.dot(tc, y - y.mean()) / np.dot(tc, tc))


@dataclass
class ConservationReport:
    energy_max: float          # max_t |dE/E0|
    energy_final: float        # dE/E0 at t_end
    energy_rms: float          # RMS of dE/E0
    energy_drift_rate: float   # slope of dE/E0 per unit time
    momentum_max: float        # max_t |P - P0| / sum_i m_i |v_i(0)|
    angular_momentum_max: float  # max_t |L - L0| / |L0|
    com_drift_max: float       # max_t |R(t) - R0 - t V0| / length scale

    def as_dict(self) -> dict:
        return asdict(self)


def conservation_report(traj) -> ConservationReport:
    """Compute all conservation metrics of a :class:`~nbody.simulation.Trajectory`."""
    t, x, v, m = traj.t, traj.x, traj.v, traj.masses
    dE = relative_energy_error(traj.energy())
    P = diag.linear_momentum(v, m)
    p_scale = float(np.sum(m * np.linalg.norm(v[0], axis=-1)))
    p_scale = p_scale if p_scale > 0 else 1.0
    L = diag.angular_momentum(x, v, m)
    L0 = L[0]
    L_norm0 = float(np.linalg.norm(L0))
    dL = np.linalg.norm(np.atleast_1d(L - L0).reshape(L.shape[0], -1), axis=-1)
    L_err = dL / L_norm0 if L_norm0 > 0 else dL
    R = diag.center_of_mass(x, m)
    V0 = diag.center_of_mass_velocity(v[0], m)
    length = float(np.max(np.linalg.norm(x[0] - R[0], axis=-1))) or 1.0
    com_dev = np.linalg.norm(R - R[0] - (t - t[0])[:, None] * V0, axis=-1) / length
    return ConservationReport(
        energy_max=float(np.max(np.abs(dE))),
        energy_final=float(dE[-1]),
        energy_rms=float(np.sqrt(np.mean(dE**2))),
        energy_drift_rate=linear_drift_rate(t, dE),
        momentum_max=float(np.max(np.linalg.norm(P - P[0], axis=-1)) / p_scale),
        angular_momentum_max=float(np.max(L_err)),
        com_drift_max=float(np.max(com_dev)),
    )


# ---------------------------------------------------------------------------
# 2. Accuracy and convergence
# ---------------------------------------------------------------------------
def position_error(x: np.ndarray, x_ref: np.ndarray) -> np.ndarray:
    """Largest Euclidean position deviation over all bodies, shape ``(...)``."""
    return np.max(np.linalg.norm(np.asarray(x) - np.asarray(x_ref), axis=-1), axis=-1)


def phase_space_error(x, v, x_ref, v_ref, length: float = 1.0, velocity: float = 1.0):
    """Euclidean error in scaled phase space (positions / length, velocities / velocity)."""
    dx = (np.asarray(x) - x_ref) / length
    dv = (np.asarray(v) - v_ref) / velocity
    return np.sqrt(np.sum(dx**2, axis=(-2, -1)) + np.sum(dv**2, axis=(-2, -1)))


@dataclass
class OrderFit:
    order: float
    log_constant: float
    n_points: int

    def as_dict(self) -> dict:
        return asdict(self)


def estimate_order(
    h: np.ndarray, err: np.ndarray, err_min: float = 0.0, err_max: float = np.inf
) -> OrderFit:
    r"""Empirical convergence order from :math:`\log\,\mathrm{err} = p\log h + c`.

    Only points with ``err_min < err < err_max`` enter the fit, which removes
    the round-off plateau (small errors) and the pre-asymptotic regime (large
    errors).  At least two points are required.
    """
    h = np.asarray(h, dtype=np.float64)
    err = np.asarray(err, dtype=np.float64)
    mask = np.isfinite(err) & (err > err_min) & (err < err_max) & (err > 0)
    if mask.sum() < 2:
        return OrderFit(float("nan"), float("nan"), int(mask.sum()))
    p, c = np.polyfit(np.log(h[mask]), np.log(err[mask]), 1)
    return OrderFit(float(p), float(c), int(mask.sum()))


# ---------------------------------------------------------------------------
# 3. Structure: time reversibility
# ---------------------------------------------------------------------------
def reversibility_error(
    model: ForceModel, state: State, integrator: Integrator, dt: float, n_steps: int
) -> float:
    r"""Error of integrating ``n_steps`` forward with ``dt`` and back with ``-dt``.

    For a symmetric method :math:`\Phi_{-h}=\Phi_h^{-1}` and the result is at
    round-off level; for non-symmetric methods it is of the size of the
    global error.  Returns the maximum position deviation from the start.
    """
    if integrator.adaptive:
        raise ValueError("reversibility is only defined for fixed step sizes")
    st = state.copy()
    st.invalidate_cache()
    for _ in range(n_steps):
        integrator.step(st, dt, model)
    for _ in range(n_steps):
        integrator.step(st, -dt, model)
    return float(position_error(st.x, state.x))


# ---------------------------------------------------------------------------
# 4. Physical observables
# ---------------------------------------------------------------------------
def mean_motion_period(t: np.ndarray, rel: np.ndarray) -> float:
    """Orbital period from a linear fit of the unwrapped polar angle of ``rel(t)``.

    ``rel`` is the (T, 2) relative position of a satellite with respect to its
    primary.  The fit averages over the anomalistic variation of the angular
    velocity, so the result is the mean (sidereal) period.
    """
    phi = np.unwrap(np.arctan2(rel[:, 1], rel[:, 0]))
    rate = linear_drift_rate(t, phi)
    return float(2.0 * np.pi / abs(rate))


def finite_time_lyapunov(t: np.ndarray, separation: np.ndarray, t_min: float = 0.0) -> float:
    r"""Largest finite-time Lyapunov exponent :math:`\lambda` from
    :math:`\ln\delta(t)\approx\lambda t + c`, fitted for ``t >= t_min``."""
    t = np.asarray(t, dtype=np.float64)
    s = np.asarray(separation, dtype=np.float64)
    mask = (t >= t_min) & (s > 0) & np.isfinite(s)
    return linear_drift_rate(t[mask], np.log(s[mask]))

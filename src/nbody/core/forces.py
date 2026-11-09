r"""Vectorised gravitational kernels.

All functions implement Newtonian gravity with Plummer softening

.. math::

    \Phi_{ij} = -\frac{G m_i m_j}{\sqrt{|\mathbf r_j-\mathbf r_i|^2+\varepsilon^2}},
    \qquad
    \mathbf a_i = G\sum_{j\neq i} m_j
        \frac{\mathbf r_j-\mathbf r_i}{\left(|\mathbf r_j-\mathbf r_i|^2+\varepsilon^2\right)^{3/2}} .

The acceleration is exactly the negative gradient of the softened potential
divided by the mass, so energy diagnostics computed with
:func:`potential_energy` are consistent with the dynamics (for
``softening = 0`` both reduce to the classical Newtonian expressions).

The direct summation costs :math:`\mathcal O(N^2)` operations and is evaluated
with NumPy broadcasting: the separation tensor ``dx[i, j] = x[j] - x[i]`` of
shape ``(N, N, d)`` is formed once and contracted with the weight matrix
``w[i, j] = m_j / (r_ij^2 + eps^2)^{3/2}``.
"""

from __future__ import annotations

import numpy as np

# Bodies x time-frames above which the potential is evaluated in chunks to
# bound the memory footprint of the pair-separation tensor.
_MAX_PAIR_ELEMENTS = 4_000_000


def pairwise_separations(x: np.ndarray) -> np.ndarray:
    """Return ``dx`` with ``dx[i, j] = x[j] - x[i]`` (shape ``(N, N, d)``)."""
    return x[np.newaxis, :, :] - x[:, np.newaxis, :]


def acceleration(
    x: np.ndarray, masses: np.ndarray, G: float = 1.0, softening: float = 0.0
) -> np.ndarray:
    """Gravitational acceleration of every body, shape ``(N, d)``.

    Parameters
    ----------
    x:
        Positions, shape ``(N, d)``.
    masses:
        Masses, shape ``(N,)``.  Zero masses (test particles) are allowed.
    G:
        Gravitational constant in the chosen unit system.
    softening:
        Plummer softening length :math:`\\varepsilon \\geq 0`.

    For large ``N`` the rows of the interaction matrix are processed in blocks
    so that the temporary separation tensor never exceeds a few ten MB.
    """
    n, d = x.shape
    block = max(1, _MAX_PAIR_ELEMENTS // max(1, n * d))
    if block >= n:
        return G * _acceleration_rows(x, x, masses, 0, softening)
    out = np.empty_like(x)
    for start in range(0, n, block):
        out[start : start + block] = _acceleration_rows(
            x[start : start + block], x, masses, start, softening
        )
    return G * out


def _acceleration_rows(
    xi: np.ndarray, x: np.ndarray, masses: np.ndarray, offset: int, softening: float
) -> np.ndarray:
    """Accelerations (without G) of the bodies ``offset .. offset+len(xi)``."""
    dx = x[np.newaxis, :, :] - xi[:, np.newaxis, :]  # (B, N, d)
    r2 = np.einsum("ijk,ijk->ij", dx, dx)
    if softening > 0.0:
        r2 += softening * softening
    # Self-interaction entries (j == offset + i) have dx = 0: set r2 = 1 there
    # to avoid 0/0 and remove them explicitly through the weights.
    rows = np.arange(xi.shape[0])
    r2[rows, rows + offset] = 1.0
    w = masses[np.newaxis, :] / (r2 * np.sqrt(r2))
    w[rows, rows + offset] = 0.0
    return np.einsum("ij,ijk->ik", w, dx)


class PairKernel:
    r"""Pair-list kernel for small ``N`` (lower call overhead than the dense one).

    Each unordered pair :math:`p=(i,j)`, :math:`i<j`, is evaluated once, and the
    pair forces are scattered to the bodies with a precomputed sparse-in-spirit
    matrix :math:`S \in \mathbb R^{N\times P}` with :math:`S_{i p} = m_j` and
    :math:`S_{j p} = -m_i`, so that :math:`\mathbf a = G\,S\,\mathbf f` with
    :math:`\mathbf f_p = (\mathbf r_j-\mathbf r_i)/(r_{ij}^2+\varepsilon^2)^{3/2}`.
    Newton's third law is built in.  The matrix has :math:`N\cdot N(N-1)/2`
    entries, so this kernel is only used for small systems.
    """

    def __init__(self, masses: np.ndarray) -> None:
        n = masses.size
        self.i, self.j = np.triu_indices(n, k=1)
        p = np.arange(self.i.size)
        self.scatter = np.zeros((n, self.i.size))
        self.scatter[self.i, p] = masses[self.j]
        self.scatter[self.j, p] = -masses[self.i]

    def __call__(self, x: np.ndarray, G: float, softening: float) -> np.ndarray:
        d = x[self.j] - x[self.i]
        r2 = np.einsum("pk,pk->p", d, d)
        if softening > 0.0:
            r2 += softening * softening
        d *= (1.0 / (r2 * np.sqrt(r2)))[:, np.newaxis]
        return G * (self.scatter @ d)


def potential_energy(
    x: np.ndarray, masses: np.ndarray, G: float = 1.0, softening: float = 0.0
) -> np.ndarray | float:
    """Total (softened) gravitational potential energy.

    ``x`` may carry leading batch axes, i.e. have shape ``(..., N, d)``; the
    result then has shape ``(...)``.  Each unordered pair is counted once.
    """
    x = np.asarray(x, dtype=np.float64)
    n = x.shape[-2]
    if n < 2:
        return np.zeros(x.shape[:-2]) if x.ndim > 2 else 0.0
    i, j = np.triu_indices(n, k=1)
    mm = masses[i] * masses[j]
    eps2 = softening * softening

    def _pot(xx: np.ndarray) -> np.ndarray:
        d = xx[..., j, :] - xx[..., i, :]
        r = np.sqrt(np.einsum("...k,...k->...", d, d) + eps2)
        return -G * np.sum(mm / r, axis=-1)

    if x.ndim == 2:
        return float(_pot(x))

    # Batched evaluation (e.g. a whole trajectory), chunked over the leading axis.
    batch_shape = x.shape[:-2]
    flat = x.reshape((-1,) + x.shape[-2:])
    chunk = max(1, _MAX_PAIR_ELEMENTS // max(1, len(i) * x.shape[-1]))
    out = np.empty(flat.shape[0])
    for start in range(0, flat.shape[0], chunk):
        out[start : start + chunk] = _pot(flat[start : start + chunk])
    return out.reshape(batch_shape)


def acceleration_loop(
    x: np.ndarray, masses: np.ndarray, G: float = 1.0, softening: float = 0.0
) -> np.ndarray:
    """Reference implementation with an explicit pair loop.

    It mirrors the algorithm of the original course project (Newton's third
    law is used to visit each pair once) and is kept only for validation and
    for the performance comparison in the evaluation chapter.
    """
    n, d = x.shape
    acc = np.zeros((n, d))
    eps2 = softening * softening
    for i in range(n - 1):
        for j in range(i + 1, n):
            r_vec = x[j] - x[i]
            r2 = float(r_vec @ r_vec) + eps2
            f = G / (r2 * np.sqrt(r2)) * r_vec
            acc[i] += masses[j] * f
            acc[j] -= masses[i] * f
    return acc

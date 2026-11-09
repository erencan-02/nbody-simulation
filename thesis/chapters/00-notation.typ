#heading(level: 1, numbering: none, outlined: false)[Notation]

#set par(justify: false)
#let sym-table(..rows) = table(
  columns: (5.2em, 1fr),
  stroke: none,
  inset: (x: 4pt, y: 2.2pt),
  ..rows,
)

#text(font: "Roboto", weight: "bold")[Physical quantities]
#sym-table(
  $N$, [number of bodies],
  $d$, [spatial dimension ($d = 2$ in all experiments, the code supports $d = 3$)],
  $m_i, M$, [mass of body $i$, total mass $M = sum_i m_i$],
  $bold(r)_i, bold(v)_i, bold(a)_i$, [position, velocity and acceleration of body $i$],
  $bold(p)_i$, [momentum $bold(p)_i = m_i bold(v)_i$],
  $bold(r)_(i j)$, [separation vector $bold(r)_(i j) = bold(r)_j - bold(r)_i$, $r_(i j) = |bold(r)_(i j)|$],
  $G$, [gravitational constant (unit dependent, e.g. $G = 1$ or $G = 4 pi^2 "AU"^3 M_dot.o^(-1) "yr"^(-2)$)],
  $epsilon$, [Plummer softening length],
  $T, V, H$, [kinetic energy, potential energy, Hamiltonian $H = T + V$],
  $E$, [total energy (value of $H$ along a solution)],
  $bold(P), bold(L), bold(R)$, [total linear momentum, angular momentum, centre of mass],
  $a, e$, [semi-major axis and eccentricity of a Kepler orbit],
  $P_"orb"$, [orbital period; $n = 2 pi \/ P_"orb"$ mean motion],
  $mu$, [gravitational parameter $mu = G(m_1 + m_2)$ of a two-body problem],
)

#text(font: "Roboto", weight: "bold")[Numerical quantities]
#sym-table(
  $bold(y)$, [phase-space vector $bold(y) = (bold(q), bold(p))$ or $(bold(x), bold(v))$],
  $h$, [time step (also $Delta t$ in the code: `dt`)],
  $t_n, bold(y)_n$, [discrete time $t_n = n h$ and numerical approximation of $bold(y)(t_n)$],
  $phi_t$, [exact flow of the differential equation],
  $Phi_h$, [one-step map of a numerical method],
  $p$, [order of accuracy (global error $cal(O)(h^p)$)],
  $J$, [canonical symplectic matrix $J = mat(0, I; -I, 0)$],
  $tilde(H)$, [modified Hamiltonian of a symplectic method],
  $delta E(t)$, [relative energy error $(E(t) - E_0) \/ |E_0|$],
  $gamma_k$, [weights of a composition method],
  $lambda$, [(finite-time) Lyapunov exponent],
)

#import "../lib/utils.typ": *
#import "@preview/cetz:0.4.2"

= Physical Background <ch:physics>

This chapter collects the physics that the simulation implements and that the evaluation relies on. We start from Newton's law of gravitation (@sec:newton), reformulate it in Hamiltonian form (@sec:hamilton), derive the conservation laws used as accuracy diagnostics (@sec:conservation), solve the two-body problem exactly to obtain a reference solution (@sec:kepler), discuss what changes for $N >= 3$ (@sec:chaos), and finally treat softening (@sec:softening) and units (@sec:units).

== Newtonian gravitation and equations of motion <sec:newton>

Consider $N$ point masses $m_i > 0$ at positions $bold(r)_i in RR^d$, $i = 1, dots, N$. By Newton's law of gravitation, body $j$ exerts on body $i$ the force

$ bold(F)_(i j) = G frac(m_i m_j, r_(i j)^3) bold(r)_(i j), quad bold(r)_(i j) = bold(r)_j - bold(r)_i, quad r_(i j) = |bold(r)_(i j)|, $ <eq:newton-force>

directed from $i$ towards $j$ (@fig:two-bodies). The force is antisymmetric, $bold(F)_(j i) = -bold(F)_(i j)$ (Newton's third law), and *central*, i.e. parallel to the line connecting the bodies. Gravitational forces superpose linearly, so Newton's second law $m_i dot.double(bold(r))_i = sum_(j != i) bold(F)_(i j)$ yields the equations of motion

$ dot.double(bold(r))_i = bold(a)_i (bold(r)_1, dots, bold(r)_N) = G sum_(j != i) m_j frac(bold(r)_j - bold(r)_i, |bold(r)_j - bold(r)_i|^3), quad i = 1, dots, N. $ <eq:eom>

These are $d N$ coupled, nonlinear second-order ordinary differential equations (ODEs), or $2 d N$ first-order ones for positions and velocities. The right-hand side is singular at collisions $r_(i j) -> 0$, and evaluating it costs $cal(O)(N^2)$ operations because every pair of bodies interacts.

#figure(
  cetz.canvas(length: 1cm, {
    import cetz.draw: *
    let blue = rgb("#243572")
    let red = rgb("#E6001A")
    let grey = rgb("#898989")
    let O = (0, 0)
    let ri = (1.2, 1.6)
    let rj = (5.6, 2.4)
    // origin and position vectors
    circle(O, radius: 0.04, fill: black)
    content((-0.25, -0.2), $O$)
    line(O, (ri.at(0) - 0.12, ri.at(1) - 0.16), mark: (end: ">", fill: grey), stroke: grey + 0.7pt)
    line(O, (rj.at(0) - 0.4, rj.at(1) - 0.17), mark: (end: ">", fill: grey), stroke: grey + 0.7pt)
    content((0.35, 1.0), text(fill: grey, $bold(r)_i$))
    content((3.2, 0.95), text(fill: grey, $bold(r)_j$))
    // separation vector
    line((ri.at(0) + 0.05, ri.at(1) + 0.6), (rj.at(0) - 0.1, rj.at(1) + 0.7),
      mark: (end: ">", fill: blue), stroke: (paint: blue, thickness: 0.6pt, dash: "dashed"))
    content((3.4, 3.45), text(fill: blue, $bold(r)_(i j) = bold(r)_j - bold(r)_i$))
    // forces
    line((ri.at(0) + 0.24, ri.at(1) + 0.04), (ri.at(0) + 1.5, ri.at(1) + 0.27),
      mark: (end: ">", fill: red), stroke: red + 1.1pt)
    content((2.25, 1.42), text(fill: red, $bold(F)_(i j)$))
    line((rj.at(0) - 0.42, rj.at(1) - 0.07), (rj.at(0) - 1.7, rj.at(1) - 0.3),
      mark: (end: ">", fill: red), stroke: red + 1.1pt)
    content((4.75, 1.6), text(fill: red, $bold(F)_(j i) = -bold(F)_(i j)$))
    // bodies
    circle(ri, radius: 0.24, fill: blue.lighten(75%), stroke: blue + 0.8pt)
    content(ri, text(size: 7pt, $m_i$))
    circle(rj, radius: 0.4, fill: blue.lighten(75%), stroke: blue + 0.8pt)
    content(rj, text(size: 8pt, $m_j$))
  }),
  caption: [Two bodies exert equal and opposite central forces on each other. The force on body $i$ points along the separation vector $bold(r)_(i j)$.],
) <two-bodies>

== Hamiltonian formulation <sec:hamilton>

The geometric properties that distinguish good integrators from bad ones are formulated most naturally in the Hamiltonian picture @goldstein2002classical @hairer2006geometric. Collect the positions in $bold(q) = (bold(r)_1, dots, bold(r)_N) in RR^(d N)$ and the momenta $bold(p)_i = m_i dot(bold(r))_i$ in $bold(p) in RR^(d N)$. The Hamiltonian is the total energy

$ H(bold(q), bold(p)) = underbrace(sum_(i=1)^N frac(|bold(p)_i|^2, 2 m_i), T(bold(p))) underbrace(- sum_(i < j) frac(G m_i m_j, r_(i j)), V(bold(q))), $ <eq:hamiltonian>

and the equations of motion are Hamilton's equations

$ dot(bold(q)) = frac(partial H, partial bold(p)) = M^(-1) bold(p), quad dot(bold(p)) = -frac(partial H, partial bold(q)) = -nabla V(bold(q)), $ <eq:hamilton>

with the diagonal mass matrix $M$. One checks that $-partial V \/ partial bold(r)_i = sum_(j != i) bold(F)_(i j)$, so @eq:hamilton is equivalent to @eq:eom. With the phase-space vector $bold(y) = (bold(q), bold(p)) in RR^(2 d N)$ and the canonical matrix

$ J = mat(0, I; -I, 0), quad J^T = J^(-1) = -J, $

Hamilton's equations take the compact form $dot(bold(y)) = J nabla H(bold(y))$.

Two structural features of @eq:hamiltonian are exploited by the integrators of @ch:numerics:

- *Separability.* $H = T(bold(p)) + V(bold(q))$. The sub-systems with Hamiltonian $T$ alone ("drift": $dot(bold(q)) = M^(-1) bold(p)$, $dot(bold(p)) = 0$) and $V$ alone ("kick": $dot(bold(q)) = 0$, $dot(bold(p)) = -nabla V$) can be solved *exactly*. Splitting methods are built from these exact sub-flows.
- *Symplecticity.* Let $phi_t$ denote the exact flow, i.e. the map that sends an initial value $bold(y)_0$ to the solution $bold(y)(t)$.

#definition[symplectic map][
  A differentiable map $Psi: RR^(2 m) -> RR^(2 m)$ is _symplectic_ if its Jacobian satisfies
  $ Psi'(bold(y))^T J Psi'(bold(y)) = J quad "for all" bold(y). $ <eq:symplectic>
]

*Theorem (Poincaré).* The flow $phi_t$ of a Hamiltonian system is symplectic for every $t$.

_Proof._ The Jacobian $Psi(t) = partial phi_t (bold(y)_0) \/ partial bold(y)_0$ solves the variational equation $dot(Psi) = J nabla^2 H(bold(y)(t)) Psi$ with $Psi(0) = I$. Using $J^T J = I$, $J J = -I$ and the symmetry of the Hessian $nabla^2 H$,
$ frac(dif, dif t) (Psi^T J Psi) = Psi^T nabla^2 H J^T J Psi + Psi^T J J nabla^2 H Psi = Psi^T nabla^2 H Psi - Psi^T nabla^2 H Psi = 0, $
so $Psi^T J Psi = Psi(0)^T J Psi(0) = J$. #h(1fr) $square$

Taking determinants of @eq:symplectic gives $det Psi' = 1$: symplectic maps preserve phase-space volume (Liouville's theorem). For $m = 1$ symplecticity is *equivalent* to area preservation; in higher dimensions it is a stronger property: the sum of oriented areas of projections onto the $(q_k, p_k)$ planes is preserved @hairer2006geometric. Symplecticity is the property a numerical method must inherit to show the long-time behaviour of Hamiltonian dynamics (@sec:backward-error).

== Conservation laws <sec:conservation>

The exact solution of @eq:eom has several invariants. They follow from symmetries of the Hamiltonian via Noether's theorem @noether1918invariante: invariance under time translation gives conservation of energy, invariance under spatial translation conservation of momentum, and invariance under rotation conservation of angular momentum. Since these invariants are known exactly, their numerical violation is the most important accuracy diagnostic of an $N$-body code.

*Energy.* Along any solution,
$ frac(dif H, dif t) = nabla H(bold(y))^T dot(bold(y)) = nabla H^T J nabla H = 0, $
because $J$ is skew-symmetric. Hence $E = T + V$ is constant.

*Linear momentum.* The total momentum $bold(P) = sum_i m_i bold(v)_i$ changes by the sum of all forces,
$ dot(bold(P)) = sum_i sum_(j != i) bold(F)_(i j) = sum_(i < j) (bold(F)_(i j) + bold(F)_(j i)) = bold(0) $
by Newton's third law. Consequently the centre of mass $bold(R) = M^(-1) sum_i m_i bold(r)_i$ moves uniformly, $bold(R)(t) = bold(R)_0 + t bold(P) \/ M$. All scenarios in this work are transformed to the barycentric frame $bold(R) = bold(P) = bold(0)$.

*Angular momentum.* With $bold(L) = sum_i bold(r)_i times m_i bold(v)_i$,
$ dot(bold(L)) = sum_i underbrace(bold(v)_i times m_i bold(v)_i, = bold(0)) + sum_i bold(r)_i times sum_(j != i) bold(F)_(i j) = sum_(i < j) (bold(r)_i - bold(r)_j) times bold(F)_(i j) = bold(0), $
because $bold(F)_(i j)$ is parallel to $bold(r)_j - bold(r)_i$. In two dimensions only the component $L_z = sum_i m_i (x_i v_(y,i) - y_i v_(x,i))$ exists.

*Virial theorem.* The moment of inertia $I = sum_i m_i |bold(r)_i|^2$ satisfies $dot(I) = 2 sum_i m_i bold(r)_i dot bold(v)_i$ and
$ frac(1, 2) dot.double(I) = sum_i m_i |bold(v)_i|^2 + sum_i bold(r)_i dot bold(F)_i = 2 T + sum_(i < j) (bold(r)_i - bold(r)_j) dot bold(F)_(i j) = 2 T + V, $
where we used $(bold(r)_i - bold(r)_j) dot bold(F)_(i j) = -G m_i m_j \/ r_(i j)$. For a bound system $dot(I)$ stays bounded, so the time average of $dot.double(I)$ vanishes and $2 chevron.l T chevron.r + chevron.l V chevron.r = 0$. A system in virial equilibrium therefore has the *virial ratio* $Q = T \/ |V| = 1\/2$. We use this to set up star clusters in equilibrium (@sec:scenarios).

== The two-body problem <sec:kepler>

For $N = 2$ the problem can be solved exactly. The solution is our most important reference for measuring global errors.

*Reduction.* Introduce the centre of mass $bold(R)$ and the relative coordinate $bold(r) = bold(r)_2 - bold(r)_1$. Subtracting the two equations of motion gives
$ dot.double(bold(r)) = -frac(mu, r^3) bold(r), quad mu = G (m_1 + m_2), $ <eq:kepler-ode>
and the bodies move as $bold(r)_1 = bold(R) - (m_2 \/ M) bold(r)$, $bold(r)_2 = bold(R) + (m_1 \/ M) bold(r)$.

*Planarity and Kepler's second law.* The specific angular momentum $bold(h) = bold(r) times dot(bold(r))$ is constant because the force is central. Hence the motion takes place in the plane perpendicular to $bold(h)$, and the area swept out per unit time, $|bold(h)| \/ 2$, is constant.

*Orbit equation.* The Laplace--Runge--Lenz vector $bold(A) = dot(bold(r)) times bold(h) - mu bold(r) \/ r$ is conserved, since
$ frac(dif, dif t) (dot(bold(r)) times bold(h)) = -frac(mu, r^3) bold(r) times (bold(r) times dot(bold(r))) = mu (frac(dot(bold(r)), r) - frac((bold(r) dot dot(bold(r))) bold(r), r^3)) = mu frac(dif, dif t) frac(bold(r), r). $
Let $f$ (the _true anomaly_) be the angle between $bold(r)$ and $bold(A) = mu bold(e)$. Using $bold(r) dot (dot(bold(r)) times bold(h)) = (bold(r) times dot(bold(r))) dot bold(h) = h^2$,
$ bold(r) dot bold(A) = h^2 - mu r = mu r e cos f quad ==> quad r(f) = frac(h^2 \/ mu, 1 + e cos f), $ <eq:conic>
a conic section with eccentricity $e = |bold(e)|$ and focus at the origin (Kepler's first law). For $0 <= e < 1$ the orbit is an ellipse with semi-major axis $a = h^2 \/ (mu (1 - e^2))$; its closest point to the focus (_periapsis_) lies in the direction of $bold(e)$ at distance $r_p = a (1 - e)$.

*Energy, vis-viva and Kepler's third law.* At periapsis $v_p = h \/ r_p$, so the specific orbital energy is
$ epsilon = frac(v_p^2, 2) - frac(mu, r_p) = frac(mu (1 + e), 2 a (1 - e)) - frac(mu, a (1 - e)) = -frac(mu, 2 a). $ <eq:kepler-energy>
Since $epsilon$ is conserved, we obtain the vis-viva equation $v^2 = mu (2\/r - 1\/a)$ for all points of the orbit. Dividing the area $pi a b$ of the ellipse, $b = a sqrt(1 - e^2)$, by the areal velocity $h\/2$ gives Kepler's third law
$ P_"orb" = frac(2 pi a b, h) = 2 pi sqrt(frac(a^3, mu)) = 2 pi sqrt(frac(a^3, G (m_1 + m_2))). $ <eq:kepler-third>

*Kepler's equation.* To obtain the position as a function of time, parametrise the ellipse by the _eccentric anomaly_ $E$: $x = a (cos E - e)$, $y = b sin E$, with periapsis on the $x$-axis. The angular momentum becomes
$ h = x dot(y) - y dot(x) = a b dot(E) (1 - e cos E), $
and with $h = sqrt(mu a (1 - e^2)) = n a b$, where $n = sqrt(mu \/ a^3)$ is the mean motion, we get $dot(E) (1 - e cos E) = n$. Integrating from periapsis at $t = 0$ yields Kepler's equation
$ E - e sin E = n t =: cal(M). $ <eq:kepler-equation>

#result[Exact reference solution][
  Given $t$, solve @eq:kepler-equation for $E$ with Newton's method, $E_(k+1) = E_k - (E_k - e sin E_k - cal(M)) \/ (1 - e cos E_k)$, which converges quadratically from $E_0 = cal(M) + e sin cal(M)$ (or $E_0 = pi$ for $e > 0.8$). Then
  $ bold(r)(t) = vec(a (cos E - e), b sin E), quad dot(bold(r))(t) = frac(n, 1 - e cos E) vec(-a sin E, b cos E). $
  This closed-form solution is implemented in `nbody.analysis.kepler` and is used to measure global errors and convergence orders (@sec:eval-convergence).
]

== Three and more bodies: non-integrability and chaos <sec:chaos>

For $N >= 3$ the conserved quantities of @sec:conservation are not enough to reduce the problem to quadratures. Only special solutions are known in closed form, for example Lagrange's equilateral configurations or the _figure-eight_ choreography of three equal masses, found numerically by Moore and proven to exist by Chenciner and Montgomery @chenciner2000remarkable. Generic configurations are chaotic. A famous example is Burrau's _Pythagorean problem_ @burrau1913numerische, in which three bodies with masses $3$, $4$, $5$ start at rest at the vertices of a right triangle with side lengths $3$, $4$, $5$. After a sequence of very close encounters the system breaks up into a binary and an escaping body @szebehely1967complete.

*Sensitivity to initial conditions.* Small perturbations $delta bold(y)$ of a trajectory evolve according to the linearised equation $delta dot(bold(y)) = J nabla^2 H(bold(y)(t)) delta bold(y)$. In chaotic regions they grow exponentially on average,
$ |delta bold(y)(t)| approx |delta bold(y)(0)| e^(lambda t), quad lambda = lim_(t -> oo) frac(1, t) ln frac(|delta bold(y)(t)|, |delta bold(y)(0)|), $ <eq:lyapunov>
where $lambda > 0$ is the largest Lyapunov exponent. Numerical errors act like such perturbations. If a solution must be known to a tolerance $Delta$ and the initial error is $delta_0$, the *predictability horizon* is
$ t_"pred" approx frac(1, lambda) ln frac(Delta, delta_0). $ <eq:predictability>
Each additional digit of accuracy extends $t_"pred"$ only by $ln(10) \/ lambda$. Beyond $t_"pred"$, comparing individual trajectories point by point is meaningless. A long chaotic integration can only be judged through invariants (energy, angular momentum) and statistical properties. This is why conservation errors are the primary metric of @ch:evaluation, and pointwise errors are only used for integrable or stable problems.

== Gravitational softening <sec:softening>

The singularity of @eq:newton-force at $r_(i j) = 0$ causes arbitrarily large accelerations during close encounters. For _collisionless_ systems such as star clusters or galaxies, where the bodies represent a smooth mass distribution, one replaces the point-mass potential by the Plummer potential @plummer1911problem @aarseth2003gravitational:
$ V_epsilon (bold(q)) = - sum_(i < j) frac(G m_i m_j, sqrt(r_(i j)^2 + epsilon^2)). $ <eq:plummer>
Differentiating with $partial r_(i j)^2 \/ partial bold(r)_i = -2 bold(r)_(i j)$ gives the softened acceleration
$ bold(a)_i = -frac(1, m_i) frac(partial V_epsilon, partial bold(r)_i) = G sum_(j != i) m_j frac(bold(r)_(i j), (r_(i j)^2 + epsilon^2)^(3\/2)), $ <eq:softened-acc>
which is bounded by $cal(O)(G m \/ epsilon^2)$ and equals the Newtonian force for $r_(i j) >> epsilon$. The choice of $epsilon$ trades force bias at small separations against discreteness noise @das2021optimal. For planetary systems and few-body problems we use $epsilon = 0$.

#callout(title: [Consistency of force and energy])[
  The softened system is again Hamiltonian, with $H_epsilon = T + V_epsilon$. All results of @sec:hamilton and @sec:conservation remain valid, but the conserved energy is $E_epsilon = T + V_epsilon$ and *not* $T + V$. Mixing the softened force @eq:softened-acc with the Newtonian potential in the energy diagnostic produces an apparent energy error even for the exact solution, because for $r >> epsilon$
  $ V - V_epsilon approx -sum_(i < j) frac(G m_i m_j epsilon^2, 2 r_(i j)^3) $
  varies along the orbit. The original course project used exactly this inconsistent diagnostic. @sec:eval-legacy quantifies the effect.
]

== Units and scaling <sec:units>

Newton's equations are invariant under the scaling $bold(r) -> alpha bold(r)$, $t -> beta t$, $m -> gamma m$ if $alpha^3 = gamma beta^2$. Results obtained in dimensionless units therefore transfer to any physical scale. We use two unit systems:

- *$N$-body units* with $G = 1$ and characteristic mass and length of order one, for test problems. The Kepler orbit with $a = 1$, $m_1 + m_2 approx 1$ then has period $P_"orb" approx 2 pi$.
- *Astronomical units* (AU, year, solar mass). Applying Kepler's third law @eq:kepler-third to the Earth ($a = 1 "AU"$, $P_"orb" = 1 "yr"$, $m_plus.o << M_dot.o$) gives
  $ G = 4 pi^2 "AU"^3 M_dot.o^(-1) "yr"^(-2). $

Choosing units in which the relevant quantities are of order one also keeps the floating-point arithmetic well conditioned. The original project had to rescale SI values by $1\/2.5 dot 10^11$ for plotting.

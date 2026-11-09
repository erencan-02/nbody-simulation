#import "../lib/utils.typ": *

= Numerical Integration of Hamiltonian Systems <ch:numerics>

This chapter develops the integrators studied in this thesis in order of increasing sophistication. For each method we derive its order of accuracy, its stability on the linear test problem and whether it preserves the symplectic structure. These are the three properties that the experiments of @ch:evaluation test. The presentation follows Hairer, Lubich and Wanner @hairer2006geometric and Hairer, Nørsett and Wanner @hairer1993solving.

== One-step methods, order and cost <sec:one-step>

We write the equations of motion as a first-order system $dot(bold(y)) = bold(f)(bold(y))$ with $bold(y) = (bold(x), bold(v))$ and $bold(f)(bold(y)) = (bold(v), bold(a)(bold(x)))$. A _one-step method_ is a map $Phi_h$ that advances the numerical solution by one time step $h$, $bold(y)_(n+1) = Phi_h (bold(y)_n)$. It approximates the exact flow $phi_h$.

#definition[order][
  The _local error_ of $Phi_h$ is $bold(delta)(bold(y), h) = Phi_h (bold(y)) - phi_h (bold(y))$. The method has _order_ $p$ if $bold(delta) = cal(O)(h^(p+1))$ for all smooth $bold(f)$.
]

*From local to global error.* Let $bold(e)_n = bold(y)_n - bold(y)(t_n)$ and assume $Phi_h$ is Lipschitz with constant $1 + h L$. Then $|bold(e)_(n+1)| <= (1 + h L)|bold(e)_n| + C h^(p+1)$. Summing this geometric recursion ("Lady Windermere's fan" @hairer1993solving) gives
$ |bold(e)_n| <= C h^(p+1) sum_(k=0)^(n-1) (1 + h L)^k <= frac(C h^p, L) (e^(L t_n) - 1). $ <eq:global-error>
A method of order $p$ therefore has global error $cal(O)(h^p)$ on bounded intervals. The constant may grow exponentially with $t$, and this bound says nothing about which errors accumulate. Answering that question requires the geometric analysis of @sec:backward-error.

*Cost.* In all methods considered here the evaluation of the accelerations $bold(a)(bold(x))$, which costs $cal(O)(N^2)$, dominates the cost. We therefore measure cost in _force evaluations_ and compare methods at equal cost rather than at equal step size. A method with $s$ evaluations per step must take steps that are $s$ times larger to be as cheap as a method with one evaluation per step.

*Linear test problem.* The harmonic oscillator $dot.double(x) = -omega^2 x$ is the model problem for oscillatory motion. Every method below is linear on it, $bold(y)_(n+1) = M(h omega) bold(y)_n$ with $bold(y) = (x, v)$. The method is *stable* if the powers $M^n$ stay bounded. For the $2 times 2$ matrices considered here, which have distinct eigenvalues away from the stability boundary, this is equivalent to a spectral radius $rho(M) <= 1$. The exact flow is a rotation in the $(omega x, v)$ plane with $|det| = 1$ and eigenvalues $e^(plus.minus i omega h)$ on the unit circle. Its energy $E = (v^2 + omega^2 x^2)\/2$ is conserved.

== Explicit Euler method

The simplest method follows from the Taylor expansion $bold(y)(t + h) = bold(y)(t) + h bold(f)(bold(y)(t)) + cal(O)(h^2)$:
$ bold(x)_(n+1) = bold(x)_n + h bold(v)_n, quad bold(v)_(n+1) = bold(v)_n + h bold(a)(bold(x)_n). $ <eq:euler>
Its local error is $frac(h^2, 2) dot.double(bold(y)) + cal(O)(h^3)$, so the method has order $p = 1$ and needs one force evaluation per step.

On the harmonic oscillator,
$ M_"E" = mat(1, h; -h omega^2, 1), quad det M_"E" = 1 + h^2 omega^2, quad rho(M_"E") = sqrt(1 + h^2 omega^2) > 1. $
The method is *unstable for every step size*. The energy grows by the factor $1 + h^2 omega^2$ in every step, so the numerical orbit spirals outwards (@fig:harmonic). For the Kepler problem the same mechanism lets the energy increase steadily until the planet escapes.

== Symplectic Euler method

Exchanging the order of the two updates in @eq:euler, so that the position is advanced with the _new_ velocity, gives the symplectic (semi-implicit) Euler method
$ bold(v)_(n+1) = bold(v)_n + h bold(a)(bold(x)_n), quad bold(x)_(n+1) = bold(x)_n + h bold(v)_(n+1). $ <eq:symplectic-euler>
It is still explicit, still of order one and still needs one force evaluation per step, yet its long-time behaviour is completely different.

*Symplecticity.* The first update of @eq:symplectic-euler is the exact time-$h$ flow $phi_h^V$ of the sub-Hamiltonian $V(bold(q))$ (a "kick": positions fixed, momenta change by $-h nabla V$). The second is the exact flow $phi_h^T$ of $T(bold(p))$ (a "drift"). Hence
$ Phi_h^"SE" = phi_h^T compose phi_h^V $
is a composition of two exact Hamiltonian flows. Both are symplectic by Poincaré's theorem (@sec:hamilton), and since $(A B)^T J (A B) = B^T (A^T J A) B = B^T J B = J$, so is their composition.

*Harmonic oscillator.* Here
$ M_"SE" = mat(1 - h^2 omega^2, h; -h omega^2, 1), quad det M_"SE" = 1, quad tr M_"SE" = 2 - h^2 omega^2. $
For $h omega < 2$ the eigenvalues are complex conjugates of modulus one, so the method is stable and neither damps nor amplifies. A direct calculation shows that the quadratic form
$ tilde(E)_"SE" = frac(1, 2) (v^2 + omega^2 x^2 - h omega^2 x v) $ <eq:se-invariant>
is *exactly* conserved: with $x' = (1 - h^2 omega^2) x + h v$ and $v' = v - h omega^2 x$ all terms of order $h$ and higher cancel in $v'^2 + omega^2 x'^2 - h omega^2 x' v'$. The true energy $E = tilde(E)_"SE" + frac(h omega^2, 2) x v$ therefore oscillates with amplitude $cal(O)(h)$ but never drifts. The numerical orbit is a slightly tilted ellipse rather than an outward spiral (@fig:harmonic).

#callout(title: [Remark: symplectic Euler is conjugate to leapfrog])[
  Symplectic Euler and the leapfrog method of the next section generate the same sequence of positions. Symplectic Euler's velocity $bold(v)_n$ approximates the leapfrog half-step velocity $bold(v)_(n-1\/2)$, which differs from $bold(v)(t_n)$ by $frac(h, 2) bold(a)(bold(x)_n)$. Its global error is therefore $cal(O)(h)$ in general. At special instants this $cal(O)(h)$ term can cancel, and the method then _appears_ to be of second order. For a Kepler orbit started at periapsis this happens after every full period, because the offset is perpendicular to the velocity there. For this reason @sec:eval-convergence measures the *maximum* error over the whole time interval.
]

== Störmer--Verlet / leapfrog <sec:leapfrog>

*Derivation by splitting.* The exact sub-flows $phi^T$ and $phi^V$ can be combined symmetrically (Strang splitting @strang1968construction):
$ Phi_h^"LF" = phi_(h\/2)^V compose phi_h^T compose phi_(h\/2)^V, $
which written out is the _velocity Verlet_ or _kick--drift--kick leapfrog_ method @verlet1967computer:
$ bold(v)_(n+1\/2) &= bold(v)_n + frac(h, 2) bold(a)(bold(x)_n), \
  bold(x)_(n+1) &= bold(x)_n + h bold(v)_(n+1\/2), \
  bold(v)_(n+1) &= bold(v)_(n+1\/2) + frac(h, 2) bold(a)(bold(x)_(n+1)). $ <eq:leapfrog>
Eliminating the velocities gives Störmer's two-step form $bold(x)_(n+1) - 2 bold(x)_n + bold(x)_(n-1) = h^2 bold(a)(bold(x)_n)$, the central-difference discretisation of $dot.double(bold(x)) = bold(a)$.

*Cost.* The acceleration $bold(a)(bold(x)_(n+1))$ computed at the end of a step is the one needed at the beginning of the next ("first same as last", FSAL). If it is cached, every step after the first costs *one* force evaluation, the same as the Euler methods. The original course implementation recomputed it and used two evaluations per step.

*Properties.*
- _Symplectic_, as a composition of exact Hamiltonian flows.
- _Symmetric (time-reversible)_: $Phi_(-h) compose Phi_h = "id"$. In $phi_(-h\/2)^V phi_(-h)^T phi_(-h\/2)^V phi_(h\/2)^V phi_h^T phi_(h\/2)^V$ the inner kicks cancel, then the drifts, then the outer kicks. Integrating forward and then backward with $-h$ returns exactly to the initial state, up to round-off.
- _Order two._ The method is consistent, and a symmetric method has even order: if $Phi_h = phi_h + C h^(p+1) + cal(O)(h^(p+2))$, then its adjoint $Phi_(-h)^(-1)$ has leading error $(-1)^p C h^(p+1)$ @hairer2006geometric. Symmetry, $Phi_h = Phi_(-h)^(-1)$, forces $C = (-1)^p C$, which for $C != 0$ requires $p$ even.

*Harmonic oscillator.* Multiplying the three linear maps gives
$ M_"LF" = mat(1 - frac(h^2 omega^2, 2), h; -h omega^2 (1 - frac(h^2 omega^2, 4)), 1 - frac(h^2 omega^2, 2)), quad det M_"LF" = 1, quad tr M_"LF" = 2 - h^2 omega^2. $
The eigenvalues are $e^(plus.minus i theta)$ with $cos theta = 1 - h^2 omega^2 \/ 2$, so the method is stable if and only if $h omega < 2$. Writing $theta = tilde(omega) h$ defines the *numerical frequency*
$ tilde(omega) = frac(2, h) arcsin frac(h omega, 2) = omega (1 + frac(h^2 omega^2, 24) + cal(O)(h^4)). $ <eq:numerical-frequency>
The amplitude is exact, but the phase advances slightly too fast, and the phase error grows *linearly* in time. As for symplectic Euler, a modified energy is conserved exactly:
$ tilde(E)_"LF" = frac(1, 2) (v^2 + omega^2 (1 - frac(h^2 omega^2, 4)) x^2). $ <eq:lf-invariant>
This can be checked by direct substitution: with $c = 1 - h^2 omega^2 \/ 2$ and $k = 1 - h^2 omega^2 \/ 4$ one has $c^2 + k h^2 omega^2 = 1$, which is exactly the condition that $v'^2 + k omega^2 x'^2 = v^2 + k omega^2 x^2$.

== Backward error analysis and long-time behaviour <sec:backward-error>

The modified energies @eq:se-invariant and @eq:lf-invariant are special cases of a general theorem that explains the excellent long-time behaviour of symplectic integrators.

*Modified equations.* For a method of order $p$ one can construct a _modified differential equation_
$ dot(tilde(bold(y))) = bold(f)(tilde(bold(y))) + h^p bold(f)_(p+1)(tilde(bold(y))) + h^(p+1) bold(f)_(p+2)(tilde(bold(y))) + dots $
whose exact solution reproduces the numerical solution, $tilde(bold(y))(n h) = bold(y)_n$, as a formal series. This is backward error analysis: the numerical solution is the exact solution of a nearby problem.

#result[Theorem (modified Hamiltonian, @hairer2006geometric[Ch. IX], @benettin1994hamiltonian)][
  If a method of order $p$ is symplectic and applied to $dot(bold(y)) = J nabla H(bold(y))$, its modified equation is again Hamiltonian, $dot(tilde(bold(y))) = J nabla tilde(H)(tilde(bold(y)))$, with
  $ tilde(H) = H + h^p H_(p+1) + h^(p+1) H_(p+2) + dots $
  After optimal truncation of the series, the numerical solution conserves $tilde(H)$ up to an error that is exponentially small in $1\/h$. Consequently
  $ |H(bold(y)_n) - H(bold(y)_0)| = cal(O)(h^p) quad "for exponentially long times" n h <= e^(c \/ h). $
]

A non-symplectic method has a modified equation that is not Hamiltonian. Generically the energy then drifts linearly, $H(bold(y)_n) - H(bold(y)_0) approx c h^p t$. Depending on the sign of $c$, the orbit gains energy (explicit Euler) or loses it (RK4 on the harmonic oscillator, @sec:runge-kutta).

*Error growth for integrable systems.* The Kepler problem shows how this difference affects the global error @calvo1993development @hairer2006geometric[Ch. X]. The mean motion depends on the energy, $n(E) = (-2E)^(3\/2) \/ mu$, which follows from @eq:kepler-energy and @eq:kepler-third. The error in the orbital phase $phi$ accumulates as
$ delta phi(t) = integral_0^t delta n dif s, quad delta n approx n'(E) delta E(s) + delta n_"num", $
where $delta n_"num" = cal(O)(h^p)$ is the frequency error of the method, as in @eq:numerical-frequency.

- For a *symplectic* method $delta E = cal(O)(h^p)$ stays bounded, so $delta phi approx c_1 h^p t$: the global error grows *linearly* in time.
- For a *non-symplectic* method $delta E approx c h^p t$, so $delta phi approx frac(1, 2) n'(E) c h^p t^2$: the global error grows *quadratically* in time.

Since the along-track position error is $a delta phi$, after long times the symplectic method wins regardless of the error constants. @sec:eval-long-term confirms both growth laws.

== Runge--Kutta methods <sec:runge-kutta>

An explicit $s$-stage Runge--Kutta method is defined by its Butcher tableau $(a_(i j))_(j < i)$, $b_i$, $c_i = sum_j a_(i j)$ @butcher2016numerical:
$ bold(k)_i = bold(f)(bold(y)_n + h sum_(j < i) a_(i j) bold(k)_j), quad i = 1, dots, s, quad bold(y)_(n+1) = bold(y)_n + h sum_(i=1)^s b_i bold(k)_i. $ <eq:rk>
Expanding the exact and numerical solutions in Taylor series and comparing coefficients gives the order conditions. For order four these are the eight equations
$ sum b_i = 1, quad sum b_i c_i = 1/2, quad sum b_i c_i^2 = 1/3, quad sum b_i a_(i j) c_j = 1/6, \
  sum b_i c_i^3 = 1/4, quad sum b_i c_i a_(i j) c_j = 1/8, quad sum b_i a_(i j) c_j^2 = 1/12, quad sum b_i a_(i j) a_(j k) c_k = 1/24. $
Two classical solutions are studied here:

#let tableau(n, ..cells) = table(
  columns: n,
  stroke: (x, y) => (if x == 0 { (right: 0.5pt) } else { (:) }) + (if y == n - 1 { (top: 0.5pt) } else { (:) }),
  inset: (x: 5pt, y: 3pt),
  align: center,
  ..cells,
)

#grid(
  columns: (1fr, 1fr),
  align: center + horizon,
  row-gutter: 6pt,
  [Heun's method (RK2), $p = 2$], [Classical RK4, $p = 4$],
  tableau(3,
    $0$, [], [],
    $1$, $1$, [],
    [], $1/2$, $1/2$),
  tableau(5,
    $0$, [], [], [], [],
    $1/2$, $1/2$, [], [], [],
    $1/2$, $0$, $1/2$, [], [],
    $1$, $0$, $0$, $1$, [],
    [], $1/6$, $1/3$, $1/3$, $1/6$),
)

For $dot(bold(x)) = bold(v)$, $dot(bold(v)) = bold(a)(bold(x))$ every stage needs one force evaluation, so Heun costs two and RK4 four evaluations per step.

*Stability function.* Applied to $dot(y) = lambda y$, @eq:rk gives $y_(n+1) = R(h lambda) y_n$ with $R(z) = 1 + z + z^2\/2 + dots + z^p \/ p!$ for the methods above. The harmonic oscillator has $lambda = plus.minus i omega$, so stability on the imaginary axis is what matters. With $z = i y$, $y = h omega$:
$ |R_"Euler" (i y)|^2 &= 1 + y^2, \
  |R_"Heun" (i y)|^2 &= (1 - y^2 \/ 2)^2 + y^2 = 1 + y^4 \/ 4, \
  |R_"RK4" (i y)|^2 &= (1 - y^2\/2 + y^4\/24)^2 + (y - y^3\/6)^2 = 1 - frac(y^6, 72) + frac(y^8, 576). $
Explicit Euler and Heun amplify for every $h > 0$, so their energy grows. RK4 is stable for $y^8 \/ 576 <= y^6 \/ 72$, i.e. $h omega <= 2 sqrt(2) approx 2.83$, but *dissipates* energy by about $(h omega)^6 \/ 72$ per step. RK4 has the largest stability interval of all methods considered here, but none of them is energy-preserving in the long run.

*No explicit Runge--Kutta method is symplectic.* A Runge--Kutta method is symplectic if and only if $b_i a_(i j) + b_j a_(j i) = b_i b_j$ for all $i, j$ @hairer2006geometric[Thm. VI.4.3]. For $i = j$ this requires $2 b_i a_(i i) = b_i^2$. An explicit method has $a_(i i) = 0$, so all $b_i = 0$, which contradicts $sum b_i = 1$. Symplectic methods must therefore be built differently: by splitting (@sec:leapfrog) or by composition (@sec:composition).

== Composition methods <sec:composition>

High-order symplectic methods are obtained by composing a basic symplectic method with suitably chosen step sizes @yoshida1990construction @suzuki1990fractal @creutz1989higher. Let $Phi_h$ be symmetric and of even order $p$, with local error $Phi_h = phi_h + h^(p+1) bold(C) + cal(O)(h^(p+2))$, and consider
$ Psi_h = Phi_(gamma_s h) compose dots.c compose Phi_(gamma_1 h). $
Since flows compose additively in time and the leading error terms add up at lowest order,
$ Psi_h = phi_((gamma_1 + dots + gamma_s) h) + (gamma_1^(p+1) + dots + gamma_s^(p+1)) h^(p+1) bold(C) + cal(O)(h^(p+2)). $
$Psi_h$ thus has order at least $p + 1$ if
$ sum_k gamma_k = 1 quad "(consistency)", quad sum_k gamma_k^(p+1) = 0 quad "(cancellation)". $ <eq:composition-conditions>
If the weights are also palindromic ($gamma_k = gamma_(s+1-k)$), $Psi_h$ is symmetric, so its order is even and hence at least $p + 2$.

*Triple jump.* With $s = 3$ and $gamma_1 = gamma_3$, @eq:composition-conditions becomes $2 gamma_1 + gamma_2 = 1$ and $2 gamma_1^(p+1) + gamma_2^(p+1) = 0$. Since $p + 1$ is odd, the second condition gives $gamma_2 = -2^(1\/(p+1)) gamma_1$, and therefore
$ gamma_1 = gamma_3 = frac(1, 2 - 2^(1\/(p+1))), quad gamma_2 = -frac(2^(1\/(p+1)), 2 - 2^(1\/(p+1))). $ <eq:triple-jump>

- Applied to leapfrog ($p = 2$), @eq:triple-jump gives $gamma_1 approx 1.3512$, $gamma_2 approx -1.7024$. This is *Yoshida's fourth-order method* @yoshida1990construction (equivalently Forest and Ruth @forest1990fourth). With FSAL it costs three force evaluations per step.
- Applied once more to Yoshida-4 ($p = 4$), the recursion yields a sixth-order method with $3^2 = 9$ leapfrog stages (9 evaluations).
- Yoshida also solved the order conditions for seven symmetric stages directly. This *sixth-order "solution A"* needs only 7 evaluations per step (coefficients in @app:coefficients).

All composition methods are symplectic (compositions of symplectic maps) and time-reversible (palindromic). Their price is the *negative* intermediate step $gamma_2 h$, which makes the error constants larger and the stability interval smaller. On the harmonic oscillator Yoshida-4 is stable only for $h omega < 1.575$, about $0.53$ per force evaluation, compared with $2$ for leapfrog (@sec:eval-harmonic). Increasing the order is therefore only worthwhile when high accuracy is required.

== Adaptive step-size control <sec:adaptive>

Fixed step sizes are inefficient when the dynamics have very different time scales, e.g. a highly eccentric orbit or a close encounter. *Embedded Runge--Kutta pairs* estimate the local error at no extra cost. The Dormand--Prince pair DOPRI5(4) @dormand1980family uses seven stages to compute a fifth-order solution $bold(y)_(n+1)$ and an embedded fourth-order solution $hat(bold(y))_(n+1)$. Their difference estimates the local error of the fourth-order solution, which is $cal(O)(h^5)$.

*Error norm.* Each component is scaled with mixed absolute and relative tolerances,
$ "err" = sqrt(frac(1, n) sum_(k=1)^n (frac(y_(n+1,k) - hat(y)_(n+1,k), "atol" + "rtol" dot max(|y_(n,k)|, |y_(n+1,k)|)))^2). $

*Step-size controller.* A step is accepted if $"err" <= 1$. Since $"err" approx C h^5$, the step size that would give $"err" = 1$ is $h dot "err"^(-1\/5)$. With a safety factor $s = 0.9$ and limits on the change,
$ h_"new" = h dot min(f_max, max(f_min, s dot "err"^(-1\/5))), quad f_min = 0.2, quad f_max = 5. $ <eq:controller>
After a rejected step the step is repeated with $h_"new"$, and no increase is allowed directly after a rejection @hairer1993solving. The seventh stage is evaluated at the new solution (FSAL), so an accepted step costs six force evaluations.

*Adaptivity versus geometry.* Adaptivity is essential for close encounters (@sec:eval-pythagorean), but it is not free. Even a symplectic method loses its long-time properties if the step size is changed according to the current state, because the modified Hamiltonian depends on $h$ @calvo1993development. DOPRI5 is in addition neither symplectic nor symmetric. Structure-preserving adaptivity, via time transformations or regularisation @hairer2006geometric[Ch. VIII] @aarseth2003gravitational, is discussed in @ch:conclusion.

== Summary of the methods

#let stab = load-csv("../data/harmonic_stability.csv")
#let stab-of(n) = {
  let r = row-of(stab, "integrator", n)
  if r == none { [--] } else {
    let v = num(r.stability_limit)
    if v < 0.01 { $0$ } else { fix(v, digits: 2) }
  }
}

@tbl:methods summarises the theoretical properties of all integrators. The stability limits were obtained numerically from the spectral radius of $M(h omega)$ (@sec:eval-harmonic) and agree with the analytic values $0$, $2$, $2$ and $2 sqrt(2)$ derived above.

#figure(
  data-table(
    columns: (1fr, auto, auto, auto, auto, auto),
    align: (left, center, center, center, center, center),
    header: ([Method], [Order $p$], [Symplectic], [Symmetric], [Evals/step], [Stability $h omega$]),
    [Explicit Euler], [1], [no], [no], [1], stab-of("euler"),
    [Symplectic Euler], [1], [yes], [no], [1], stab-of("symplectic_euler"),
    [Heun (RK2)], [2], [no], [no], [2], [$0$],
    [Leapfrog / Verlet], [2], [yes], [yes], [1], stab-of("leapfrog"),
    [Runge--Kutta 4], [4], [no], [no], [4], stab-of("rk4"),
    [Yoshida 4], [4], [yes], [yes], [3], stab-of("yoshida4"),
    [Yoshida 6 (solution A)], [6], [yes], [yes], [7], [--],
    [Yoshida 6 (triple jump)], [6], [yes], [yes], [9], [--],
    [Dormand--Prince 5(4)], [5], [no], [no], [6 (adaptive)], [--],
  ),
  caption: [Integrators studied in this thesis, in order of increasing sophistication. "Stability" is the largest $h omega$ for which the method is stable on the harmonic oscillator ($0$: unstable for every $h > 0$).],
) <methods>

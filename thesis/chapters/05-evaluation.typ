#import "../lib/utils.typ": *

#let orders = load-csv("../data/orders.csv")
#let cost = load-csv("../data/cost_to_accuracy.csv")
#let longterm = load-csv("../data/long_term_kepler.csv")
#let fig8 = load-csv("../data/figure_eight.csv")
#let solar-cons = load-csv("../data/solar_conservation.csv")
#let solar-periods = load-csv("../data/solar_periods.csv")
#let solar-phase = load-csv("../data/solar_phase_error.csv")
#let perf = load-csv("../data/performance.csv")
#let stab = load-csv("../data/harmonic_stability.csv")
#let legacy = json("../data/legacy_two_body.json")
#let pyth = json("../data/pythagorean.json")
#let perfj = json("../data/performance.json")
#let label-of(rows, n) = row-of(rows, "integrator", n).label

= Evaluation <ch:evaluation>

This chapter evaluates the integrators and the implementation on problems of increasing difficulty: a linear oscillator, an integrable two-body problem, a periodic three-body orbit, a chaotic three-body problem with close encounters, and the Solar System. It ends with a performance analysis. @sec:metrics defines the metrics; @sec:eval-harmonic to @sec:eval-performance present the experiments e01--e08, each implemented in `experiments/` and reproducible with one command.

== Metrics and methodology <sec:metrics>

*Conservation metrics.* With $E(t)$ computed from the *same* (softened) Hamiltonian that generates the dynamics:
- relative energy error $delta E(t) = (E(t) - E_0) \/ |E_0|$, its maximum $max_t |delta E|$ and its RMS;
- secular energy drift: the slope of the least-squares line through $delta E(t)$. Unlike the maximum, it separates bounded oscillation from systematic growth;
- momentum error $max_t |bold(P)(t) - bold(P)_0| \/ sum_i m_i |bold(v)_i (0)|$ (normalised because $bold(P)_0 = bold(0)$) and angular-momentum error $max_t |bold(L)(t) - bold(L)_0| \/ |bold(L)_0|$.

*Accuracy metrics.* Against a reference solution $(bold(x)^"ref", bold(v)^"ref")$:
- the phase-space error $e(t) = (|bold(x) - bold(x)^"ref"|^2 + |bold(v) - bold(v)^"ref"|^2)^(1\/2)$, summed over all bodies, in units where the orbital length and velocity scales are one;
- the empirical order $hat(p)$: the slope of a least-squares fit of $log e$ against $log h$, restricted to the asymptotic regime ($10^(-11) < e < 10^(-2)$) so that neither the round-off plateau nor the pre-asymptotic regime biases it;
- the error growth exponent $beta$ in $e(t) prop t^beta$, fitted over the last 90 % of a long run.

*Structural metrics.* The reversibility error $max_i |bold(x)_i - bold(x)_i (0)|$ after $n$ steps forward with $h$ and $n$ steps backward with $-h$.

*Physical observables.* Mean orbital periods, from a linear fit of the unwrapped polar angle (mean motion); finite-time Lyapunov exponents, from a fit of $ln delta(t)$.

*Cost.* The number of force evaluations, which is independent of hardware. Wherever methods are compared, they get the *same* evaluation budget, i.e. steps $s$ times larger for a method with $s$ evaluations per step.

*References.* The Kepler problem has the exact solution of @sec:kepler. For the figure-eight, an independent eighth-order integrator (SciPy's DOP853 @virtanen2020scipy) at $"rtol" = "atol" = 10^(-13)$ serves as reference; comparing it with a run at $10^(-12)$ shows that its own error is below $#sci(json("../data/summary.json").e05.reference_check)$.

== Linear stability and phase space <sec:eval-harmonic>

#figure(
  image("../figures/generated/harmonic_phase_energy.svg", width: 100%),
  caption: [Harmonic oscillator with $h omega = 0.2$. (a) Phase-space trajectories over two periods. (b) Relative energy error up to $t omega = 1000$.],
) <harmonic>

@fig:harmonic shows the behaviour predicted in @ch:numerics. Explicit Euler spirals outwards: its energy grows by the factor $1 + h^2 omega^2$ per step and exceeds $10^3$ after $t omega approx 100$. Symplectic Euler produces a slightly tilted ellipse, and its energy error oscillates with amplitude $approx h omega \/ 2 = 0.1$, exactly as the conserved modified energy @eq:se-invariant predicts. Leapfrog's error amplitude is $approx (h omega)^2 \/ 4 = 0.01$, as predicted by @eq:lf-invariant, and Yoshida-4's is smaller still. None of them drifts. RK4, in contrast, *loses* energy steadily, at the rate $(h omega)^6 \/ 72 approx 9 dot 10^(-7)$ per step derived in @sec:runge-kutta. Over the longest run this accumulates to $4 dot 10^(-3)$, more than Yoshida-4's bounded error although both are fourth order.

#figure(
  image("../figures/generated/harmonic_stability.svg", width: 100%),
  caption: [Spectral radius of the one-step matrix $M(h omega)$. Symplectic Euler and leapfrog have identical eigenvalues ($tr M = 2 - h^2 omega^2$ for both), so their curves coincide.],
) <stability>

@fig:stability shows the spectral radius of the one-step matrices, computed numerically from the implementation itself. The measured stability limits are #fix(row-of(stab, "integrator", "symplectic_euler").stability_limit, digits: 3) (symplectic Euler, leapfrog; theory: 2, the deviation is the grid resolution), #fix(row-of(stab, "integrator", "rk4").stability_limit, digits: 3) (RK4; theory $2 sqrt(2) = 2.828$) and #fix(row-of(stab, "integrator", "yoshida4").stability_limit, digits: 3) (Yoshida-4); explicit Euler is unstable for every $h > 0$. Inside their stability intervals the symplectic methods have $rho = 1$ exactly (no damping, no amplification), whereas RK4 has $rho < 1$.

== Convergence order and work-precision <sec:eval-convergence>

*Setup.* Kepler orbit with $e = 0.6$, mass ratio $10^(-3)$, $a = 1$, integrated over two periods with $2^4$ to $2^(17)$ steps per period. The error is the maximum phase-space error over $approx 256$ equally spaced samples, measured against the exact solution.

#figure(
  image("../figures/generated/convergence.svg", width: 100%),
  caption: [Kepler orbit ($e = 0.6$). (a) Maximum phase-space error versus step size; black lines indicate slopes $h^p$. (b) The same errors versus the number of force evaluations, including the adaptive Dormand--Prince method for tolerances $10^(-3)$ to $10^(-13)$.],
) <convergence>

#figure(
  data-table(
    columns: (auto, auto, auto, auto, auto),
    align: (left, center, center, center, right),
    header: ([Method], [$p$ (theory)], [$hat(p)$ (measured)], [fit points], [evals for $e = 10^(-6)$]),
    ..orders.map(r => (
      r.label, r.order_theory, fix(r.order_measured, digits: 2), r.fit_points,
      { let c = row-of(cost, "integrator", r.integrator); if c == none or num(c.force_evals) == none [not reached] else { int-fmt(c.force_evals) } },
    )).flatten(),
    [Dormand--Prince 5(4)], [5], [--], [--], int-fmt(row-of(cost, "integrator", "dopri5").force_evals),
  ),
  caption: [Theoretical and measured convergence orders on the Kepler orbit with $e = 0.6$, and the number of force evaluations needed for a maximum phase-space error of $10^(-6)$ over two periods (interpolated; "not reached": requires more than the largest resolution tested).],
) <orders>

*Orders.* All measured orders agree with the theory to within $0.3$ (@tbl:orders). Explicit Euler reaches its asymptotic regime only for errors below $cal(O)(1)$ and $2^(15)$ steps per period. Its order was therefore fitted in the window $e < 1$. This shows how slowly a first-order method converges on an eccentric orbit. RK4 shows a slightly higher apparent order ($approx 4.2$) because of pre-asymptotic contributions. As explained in @sec:leapfrog, symplectic Euler converges with order one only if the error is measured over the whole interval; at the final time alone it would appear second order.

*Work-precision.* At equal cost the high-order methods are orders of magnitude more accurate (@fig:convergence b). Two further observations matter for the choice of method:

+ *Over short times, symplecticity is not an advantage.* RK4 needs fewer evaluations than Yoshida-4 for $e = 10^(-6)$, because Yoshida-4 has a larger error constant (negative substep). The benefit of symplecticity lies in the long-time behaviour, studied next.
+ *Adaptivity pays off on eccentric orbits.* DOPRI5 concentrates its steps near periapsis and reaches $10^(-6)$ with the fewest evaluations of all methods, about a factor of #calc.round(num(row-of(cost, "integrator", "rk4").force_evals) / num(row-of(cost, "integrator", "dopri5").force_evals), digits: 1) fewer than RK4. Among the fixed-step methods Yoshida-6A is cheapest; it beats the 9-stage triple-jump sixth-order method by about a factor of two at equal accuracy.

== Long-term behaviour at equal cost <sec:eval-long-term>

*Setup.* Kepler orbit with $e = 0.5$ over 200 periods. Every fixed-step method gets *1008 force evaluations per period*, i.e. 1008 leapfrog steps, 252 RK4 steps, 336 Yoshida-4 steps, and so on. DOPRI5 uses $"rtol" = 10^(-10)$, which costs a similar #calc.round(num(row-of(longterm, "integrator", "dopri5").force_evals) / 200) evaluations per period.

#figure(
  image("../figures/generated/long_term_kepler.svg", width: 100%),
  caption: [Kepler orbit ($e = 0.5$) over 200 periods at equal cost. (a) Relative energy error. (b) Phase-space error with respect to the exact solution; black lines indicate linear and quadratic growth.],
) <long-term>

#figure(
  data-table(
    columns: (auto, auto, auto, auto, auto, auto, auto),
    align: (left, right, right, right, right, center, right),
    header: ([Method], [$max|delta E|$], [drift / orbit], [$max|delta L|$], [$e(T)$], [growth $beta$], [reversib.]),
    ..longterm.map(r => (
      r.label, sci(r.energy_max), sci(r.energy_drift_per_orbit), sci(r.angmom_max),
      sci(r.error_final), fix(r.error_growth_exponent, digits: 2), sci(r.reversibility_error),
    )).flatten(),
  ),
  caption: [Long-term metrics for the experiment of @fig:long-term. "Drift / orbit": slope of $delta E(t)$ per period. "Growth $beta$": exponent of $e(t) prop t^beta$ (-- where the error saturates at the orbit size). "Reversib.": position error after 1000 steps forward and backward with $h = P_"orb"\/200$.],
) <long-term-table>

@fig:long-term and @tbl:long-term-table show the central result of this thesis.

*Energy.* All symplectic methods keep the energy error *bounded*: the error oscillates with every orbit but does not grow, and the fitted drift is at the level of noise. Explicit Euler and Heun gain energy steadily; Euler's orbit is destroyed completely. RK4 and DOPRI5 drift linearly, $delta E prop t$, as predicted by backward error analysis.

*Global error.* The measured growth exponents are $beta approx 1.00$ for symplectic Euler, leapfrog, Yoshida-4 and Yoshida-6A, and $beta approx 2$ for RK4 and DOPRI5. This is exactly the linear-versus-quadratic law derived in @sec:backward-error. DOPRI5 is still the most accurate method after 200 orbits, thanks to its higher order and adaptive steps. Because its error grows quadratically, Yoshida-6A would overtake it after roughly $(e_"Y6A"(T) \/ e_"DP"(T)) dot T approx 10 dot 200 = 2000$ orbits.

*Invariants and reversibility.* Symplectic methods conserve the angular momentum to round-off ($10^(-14)$). This is a property of all symplectic partitioned methods for rotationally invariant problems, since they conserve quadratic invariants of this form @hairer2006geometric. RK4 and DOPRI5 do not. The symmetric methods return to the initial state to $10^(-13)$ when integrated backwards; the non-symmetric ones do not.

== Revisiting the original experiment <sec:eval-legacy>

The original course project evaluated a two-body system ($G = 20$, $m = 1000$ and $10$, softening $epsilon = 5$, $h = 0.1$, 5000 leapfrog steps) by the energy quotient $q(t) = |E_t \/ E_0|$, and observed periodic dips that it attributed to "large numerical errors when the satellite is strongly accelerated". The original energy diagnostic, however, used the *Newtonian* potential, while the dynamics used the *softened* force (@sec:softening).

#figure(
  image("../figures/generated/legacy_two_body.svg", width: 100%),
  caption: [Original two-body experiment. (a) Energy quotient with the inconsistent diagnostic of the original project (red) and with the consistent softened energy (blue); grey: separation of the bodies. (b) Consistent relative energy error for several methods and step sizes.],
) <legacy>

With the consistent diagnostic the maximum deviation $max|q - 1|$ drops from $#sci(legacy.max_dev_inconsistent)$ to $#sci(legacy.max_dev_consistent)$, a factor of #calc.round(legacy.artefact_factor) (@fig:legacy a). About 97 % of the reported "energy error" was therefore an artefact of the diagnostic, not of the integrator. The remaining true error behaves as theory predicts: halving the step size reduces it by a factor of $#calc.round(legacy.max_energy_error.at("leapfrog_h0.1") / legacy.max_energy_error.at("leapfrog_h0.05"), digits: 2) approx 2^2$, and the fourth-order methods reduce it to $10^(-7)$--$10^(-6)$ (@fig:legacy b). This finding motivated the strict rule in the new code that energy diagnostics are always computed from the force model that generates the dynamics.

== Periodic three-body orbit: the figure-eight <sec:eval-figure-eight>

*Setup.* Figure-eight initial conditions @chenciner2000remarkable, 50 periods, 1260 force evaluations per period for the fixed-step methods; DOPRI5 with $"rtol" = 10^(-9)$. The reference is the independent DOP853 solution.

#figure(
  image("../figures/generated/figure_eight.svg", width: 100%),
  caption: [Figure-eight choreography. (a) The three bodies chase each other along the same curve. (b) Phase-space error with respect to the DOP853 reference.],
) <figure-eight>

#figure(
  data-table(
    columns: (auto, auto, auto, auto, auto),
    align: (left, right, right, right, right),
    header: ([Method], [force evals], [$max|delta E|$], [$max|delta L|$], [$e(t = 50 P)$]),
    ..fig8.map(r => (r.label, int-fmt(r.force_evals), sci(r.energy_max), sci(r.angmom_max), sci(r.error_final))).flatten(),
  ),
  caption: [Figure-eight over 50 periods at (nearly) equal cost.],
) <figure-eight-table>

The figure-eight is a stable periodic orbit, so errors grow polynomially, and higher order translates directly into higher accuracy (@tbl:figure-eight-table). At equal cost, Yoshida-6A is about #calc.round(num(row-of(fig8, "integrator", "leapfrog").error_final) / num(row-of(fig8, "integrator", "yoshida6a").error_final)) times more accurate than leapfrog and more accurate than DOPRI5, which used 60 % more evaluations. Again, only the symplectic methods conserve the angular momentum to machine precision.

== Chaos and close encounters: the Pythagorean problem <sec:eval-pythagorean>

#let pyth-tol = load-csv("../data/pythagorean_tolerance.csv")

*Setup.* Burrau's problem @burrau1913numerische (masses 3, 4, 5 at rest), $0 <= t <= 70$. DOPRI5 with $"rtol" = 10^(-12)$ is compared with leapfrog and Yoshida-4 using the *same* number of force evaluations. A second DOPRI5 run with the position of one body perturbed by $10^(-8)$ measures the sensitivity to initial conditions. A tolerance study compares the final states for $"rtol" = 10^(-8), dots, 10^(-12)$ with an independent SciPy DOP853 solution at $"rtol" = "atol" = 3 dot 10^(-14)$.

#figure(
  image("../figures/generated/pythagorean.svg", width: 100%),
  caption: [Pythagorean three-body problem (DOPRI5, $"rtol" = 10^(-12)$). (a) Trajectories: after a sequence of close encounters the lightest body ($m = 3$) escapes and the two heavier bodies form a binary. (b) Adaptive step size and minimum pairwise distance. (c) Energy error at equal cost. (d) Phase-space separation of two runs whose initial conditions differ by $10^(-8)$.],
) <pythagorean>

*Adaptivity.* The closest approach, $r_min = #sci(pyth.r_min)$ at $t = #calc.round(pyth.t_r_min, digits: 2)$, forces DOPRI5 to reduce its step size to $#sci(pyth.h_min)$. Elsewhere it takes steps up to $#sci(pyth.h_max)$, so the step size varies over more than six orders of magnitude (@fig:pythagorean b). The step size follows the minimum distance closely; dimensionally, $h prop r_min^(3\/2)$ is the local dynamical time. With #int-fmt(pyth.dopri_evals) force evaluations (#int-fmt(pyth.dopri_steps) accepted and #int-fmt(pyth.dopri_rejected) rejected steps) the energy is conserved to $#sci(pyth.dopri_energy_max)$.

*Fixed steps fail.* With the same budget, leapfrog ($h = #sci(pyth.fixed.leapfrog.dt)$) and Yoshida-4 ($h = #sci(pyth.fixed.yoshida4.dt)$) cannot resolve the very first close encounter. The energy error jumps to $|delta E| = #sci(pyth.fixed.leapfrog.energy_max)$ and $#sci(pyth.fixed.yoshida4.energy_max)$, respectively, and the trajectories become meaningless (@fig:pythagorean c). Symplecticity does not help when the step size exceeds the time scale of the dynamics. This is precisely the regime in which the backward error expansion of @sec:backward-error breaks down.

*Chaos.* The separation of the two adaptive runs grows exponentially, with large jumps at close encounters (@fig:pythagorean d). A fit gives the finite-time Lyapunov exponent $lambda approx #calc.round(pyth.lyapunov, digits: 2)$. The separation reaches $10^(-2)$ at $t approx #calc.round(pyth.at("t_separation_1e-2"))$. Equation @eq:predictability predicts $t_"pred" approx ln(10^(-2) \/ 10^(-8)) \/ lambda approx #calc.round(calc.ln(1e6) / pyth.lyapunov)$, the right order of magnitude. It overestimates the horizon because the growth is concentrated in a few encounters rather than spread uniformly, as a single exponent assumes.

#figure(
  data-table(
    columns: (auto, auto, auto, auto, auto),
    align: (center, right, left, center, right),
    header: ([rtol], [force evals], [final binary], [outcome], [final error vs. reference]),
    ..pyth-tol.map(r => (sci(r.rtol, digits: 0), int-fmt(r.force_evals), r.binary,
      if r.outcome_correct == "yes" { [yes] } else { [no] }, sci(r.final_error_vs_reference))).flatten(),
  ),
  caption: [Pythagorean problem: final state at $t = 70$ as a function of the DOPRI5 tolerance. Reference: SciPy DOP853 at $3 dot 10^(-14)$, which yields the binary #pyth.reference_binary. "Outcome": the same pair forms the final binary as in the reference.],
) <pythagorean-tol>

*Which body escapes?* @tbl:pythagorean-tol makes the consequence of chaos explicit. With $"rtol" >= 10^(-10)$ the *qualitative* outcome, i.e. which body is ejected, is essentially random: $10^(-9)$ happens to agree with the reference, while $10^(-8)$ and $10^(-10)$ eject the wrong body. Only from $10^(-11)$ on do all runs agree with the reference and with the literature @szebehely1967complete: the bodies with $m = 4$ and $m = 5$ form the binary and $m = 3$ escapes. Even then the final positions still differ from the reference by $cal(O)(1)$, because small differences in the escape velocity accumulate linearly over the long escape. No integrator can deliver pointwise accuracy beyond the predictability horizon; only conserved quantities and robust qualitative features remain meaningful measures of accuracy there.

== The Solar System over 1000 years <sec:eval-solar>

*Setup.* Sun and eight planets (simplified coplanar model, @sec:scenarios), 1000 years, with a budget of one force evaluation per day: leapfrog with $h = 1$ d, Yoshida-4 with $h = 3$ d, RK4 with $h = 4$ d.

#figure(
  image("../figures/generated/solar_system.svg", width: 100%),
  caption: [Solar System. (a) Inner planets over four years and (b) all planets over 170 years (heliocentric, leapfrog). (c) Relative energy error over 1000 years at one force evaluation per day.],
) <solar>

#figure(
  data-table(
    columns: (auto, auto, auto, auto, auto, auto),
    align: (left, center, right, right, right, right),
    header: ([Method], [$h$ [d]], [$max|delta E|$], [drift / yr], [$max|delta L|$], [wall time]),
    ..solar-cons.map(r => (r.label, r.dt_days, sci(r.energy_max), sci(r.energy_drift_per_yr), sci(r.angmom_max), [#calc.round(num(r.wall_time_s)) s])).flatten(),
  ),
  caption: [Conservation over 1000 years of Solar System evolution (#int-fmt(row-of(solar-cons, "integrator", "leapfrog").force_evals) force evaluations per run).],
) <solar-conservation>

*Conservation.* Leapfrog and Yoshida-4 conserve the energy to $3 dot 10^(-6)$ over 1000 years without measurable drift, and the angular momentum to $10^(-13)$ (@tbl:solar-conservation). RK4 at the same cost takes only 22 steps per Mercury orbit. Its numerical damping lets Mercury spiral inwards over a few decades until a close passage to the Sun destroys the orbit, and the energy error jumps to $10^(-2)$ (@fig:solar c). The details of this catastrophe depend sensitively on the step size: in a 100-year run with a marginally different step, Mercury was ejected from the Solar System.

#figure(
  data-table(
    columns: (auto, auto, auto, auto, auto),
    align: (left, right, right, right, right),
    header: ([Planet], [$a$ [AU]], [$P_"Kepler"$ [yr]], [$P_"sim"$ [yr]], [rel. deviation]),
    ..solar-periods.map(r => (r.planet, fix(r.a_AU, digits: 3), fix(r.P_kepler_yr, digits: 4), fix(r.P_simulated_yr, digits: 4), sci(r.rel_deviation))).flatten(),
  ),
  caption: [Simulated mean orbital periods compared with Kepler's third law @eq:kepler-third (including the planet's own mass). Terrestrial planets: Yoshida-6A, $h$ = 1 d, 100 years; giant planets: leapfrog, $h$ = 1 d, 1000 years (more than 4000 steps per orbit).],
) <solar-periods>

*Physical validation.* The simulated periods of the terrestrial planets agree with Kepler's third law to better than $10^(-4)$ (@tbl:solar-periods). The giant planets deviate by up to 1 %. This is physics rather than numerical error, as the convergence of the terrestrial periods shows. The two-body law ignores the mutual perturbations of the planets; these are strongest for Saturn, owing to the near 5:2 resonance with Jupiter ("great inequality"), and for Uranus and Neptune, owing to their near 2:1 commensurability @murray1999solar. In addition, the mass enclosed by the outer orbits includes the inner giant planets ($approx 1.3 dot 10^(-3) M_dot.o$). Finally, the artificial initial condition, with every planet at perihelion simultaneously, does not correspond to the planets' mean orbits.

#figure(
  data-table(
    columns: (auto, auto, auto, auto, auto),
    align: (left, right, right, right, right),
    header: ([Planet], [$P_"ref"$ [yr]], [leapfrog (1 d)], [Yoshida 4 (3 d)], [RK4 (4 d)]),
    ..solar-phase.map(r => (r.planet, fix(r.P_reference_yr, digits: 4), sci(r.rel_period_error_leapfrog), sci(r.rel_period_error_yoshida4), sci(r.rel_period_error_rk4))).flatten(),
  ),
  caption: [Relative period error of the equal-cost runs with respect to the converged reference (first 100 years).],
) <solar-phase>

*Frequency errors.* Energy conservation does not imply phase accuracy. @tbl:solar-phase shows that Mercury's period is off by $#sci(row-of(solar-phase, "planet", "Mercury").rel_period_error_leapfrog)$ with leapfrog and $#sci(row-of(solar-phase, "planet", "Mercury").rel_period_error_yoshida4)$ with Yoshida-4, although both conserve the energy to $10^(-6)$. This is the frequency error $delta n_"num" = cal(O)(h^p)$ of @eq:numerical-frequency. It is invisible in the energy but makes the phase error grow linearly. It decreases rapidly for the outer planets, which are resolved with many more steps per orbit. For RK4 Mercury's period is wrong by #calc.round(100 * num(row-of(solar-phase, "planet", "Mercury").rel_period_error_rk4)) %, because the orbit has decayed.

== Performance <sec:eval-performance>

#figure(
  image("../figures/generated/performance.svg", width: 100%),
  caption: [Wall-clock time per force evaluation versus $N$ (#perfj.machine).],
) <performance>

@fig:performance compares the pair loop of the original implementation with the two vectorised kernels. All three scale as $cal(O)(N^2)$ for large $N$, but the vectorised kernels are faster by a factor of $approx #calc.round(perfj.speedup_n128)$ at $N = 128$ and $approx #calc.round(perfj.speedup_n512)$ at $N = 512$. For $N = 2$ the plain loop is actually faster (one pair; NumPy call overhead dominates). The pair-list kernel is the fastest choice up to $N approx 30$ and the dense kernel above, which is why `GravityModel` switches between them automatically at $N = 24$. Together with the halved number of force evaluations of the cached leapfrog, a leapfrog step is faster than in the original implementation by a factor of $approx #calc.round(2 * perfj.speedup_n128)$ at $N = 128$. With $#calc.round(perfj.t_dense_n2048_ms) "ms"$ per evaluation at $N = 2048$, direct summation in Python reaches its practical limit, which motivates the tree and multipole methods discussed in @ch:conclusion.

== Discussion

The experiments confirm the theory of @ch:numerics quantitatively, and they also show its limits:

- The *order* of a method determines its accuracy per unit cost on smooth problems (@sec:eval-convergence, @sec:eval-figure-eight). Sixth-order composition methods are the most efficient fixed-step methods once the accuracy requirement exceeds about $10^(-6)$.
- *Symplecticity* does not reduce the error over short times, but it changes how the error accumulates: bounded energy error, linear instead of quadratic phase-error growth, and exact conservation of angular momentum (@sec:eval-long-term, @sec:eval-solar). For long integrations this outweighs larger error constants.
- *Adaptivity* is indispensable when time scales vary strongly (@sec:eval-pythagorean). Combining it with symplecticity requires dedicated techniques (@ch:conclusion).
- *Correct diagnostics* are as important as correct integrators (@sec:eval-legacy). A metric that is not conserved by the exact dynamics cannot measure the quality of a numerical method.

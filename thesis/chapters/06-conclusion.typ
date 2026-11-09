#import "../lib/utils.typ": *

= Conclusion and Outlook <ch:conclusion>

== Summary

This thesis revisited the gravitational $N$-body problem of a computational physics course project and developed it into a systematic study of structure-preserving time integration.

On the *theoretical* side we derived the Hamiltonian structure of Newtonian gravity, its conservation laws and the exact solution of the Kepler problem. We then built a hierarchy of integrators of increasing sophistication, from explicit Euler through symplectic Euler, leapfrog, Heun and RK4 to Yoshida's composition methods of order four and six and an adaptive Dormand--Prince pair. For each method we derived the order, the linear stability limit and the presence or absence of symplecticity and time symmetry. Backward error analysis, in particular the modified Hamiltonian of symplectic methods, explains why the energy error of symplectic methods stays bounded. Combined with the energy dependence of the Kepler frequency, it predicts linear global-error growth for symplectic and quadratic growth for non-symplectic methods.

On the *implementation* side the original project was rebuilt as a modular package. It uses structure-of-arrays storage, vectorised force kernels, a common integrator interface with generic Runge--Kutta and composition building blocks, a registry, consistent energy diagnostics, a reproducible experiment pipeline, and an interactive viewer whose overlays show past trails and simulated future trajectories. A suite of 103 automated tests checks the implementation against the theory.

The *evaluation* confirmed the theory quantitatively:

- The measured convergence orders of all eight fixed-step methods agree with the theoretical ones to within $0.3$.
- At equal cost, every symplectic method kept the energy error bounded over 200 eccentric orbits and showed linear error growth (measured exponent $1.00$). RK4 and DOPRI5 showed quadratic growth (exponents near $2$).
- The original project's energy fluctuation was mostly an artefact of an inconsistent energy diagnostic, not a weakness of leapfrog.
- High order pays off for smooth problems: Yoshida-6A was the most accurate fixed-step method at equal cost on the Kepler and figure-eight problems.
- Close encounters, as in the Pythagorean problem, require adaptive step sizes. Fixed-step methods failed completely at equal cost there. Because the problem is chaotic, even the qualitative outcome (which body escapes) is only reproduced reliably for tolerances of $10^(-11)$ and below.
- A 1000-year Solar System integration conserved the energy to $3 dot 10^(-6)$ at one force evaluation per day, while RK4 at the same cost destroyed Mercury's orbit. The simulated periods of the terrestrial planets agree with Kepler's third law to better than $10^(-4)$.

== Recommendations

The experiments suggest a simple guideline for choosing an integrator:

#figure(
  data-table(
    columns: (auto, 1fr),
    align: (left, left),
    header: ([Situation], [Recommended method]),
    [quick, robust default; real-time visualisation], [leapfrog (1 evaluation/step, symplectic, reversible)],
    [long integrations of smooth systems, high accuracy], [Yoshida 6A (or Yoshida 4 at moderate accuracy)],
    [short integrations, moderate accuracy], [RK4 or DOPRI5 are competitive],
    [close encounters, widely varying time scales], [adaptive DOPRI5 (or a regularised method, see below)],
    [chaotic systems], [any method that conserves the invariants; judge results statistically],
    [teaching / never], [explicit Euler],
  ),
  caption: [Guideline for choosing an integrator, derived from the experiments.],
)

== Limitations and outlook

The study deliberately covers general-purpose methods. Several directions would extend it in a natural way:

- *Planetary systems.* Wisdom--Holman mixed-variable symplectic maps @wisdom1991symplectic split the Hamiltonian into Keplerian motion around the star, which is integrated exactly with the solver of @sec:kepler, and small planet--planet interactions. This permits steps of a sizeable fraction of the shortest orbital period. The `Composition` and `ForceModel` abstractions could host such a splitting directly.
- *Structure-preserving adaptivity.* Time transformations (Sundman, Poincaré) and the regularisations of Kustaanheimo--Stiefel or Burdet--Heggie make close encounters regular and keep the symplectic structure @hairer2006geometric @aarseth2003gravitational. IAS15 @rein2015ias15 shows that an adaptive 15th-order method can reach machine precision with errors that grow only as a random walk.
- *Large $N$.* The direct $cal(O)(N^2)$ summation limits the package to a few thousand bodies. Tree codes @barnes1986hierarchical and the fast multipole method @greengard1987fast reduce the cost to $cal(O)(N log N)$ or $cal(O)(N)$. Just-in-time compilation or GPU kernels would remove the Python overhead that dominates for small $N$.
- *Physics.* Three-dimensional initial conditions with inclinations from full ephemerides would allow a quantitative comparison with observed planetary motion. A first post-Newtonian correction would reproduce the relativistic perihelion advance of Mercury of 43 arcseconds per century.
- *Visualisation.* The forecast overlay could show an ensemble of perturbed forecasts, making the loss of predictability in chaotic systems (@eq:predictability) directly visible.

#import "../lib/utils.typ": *

= Introduction <ch:introduction>

== Motivation

The gravitational $N$-body problem asks for the motion of $N$ point masses that attract each other according to Newton's law of gravitation. It is one of the oldest problems of mathematical physics and is still central to modern astrophysics: planetary systems, star clusters, galaxies and the large-scale structure of the universe are all modelled as gravitational $N$-body systems @aarseth2003gravitational.

For $N = 2$ the problem is integrable and is solved by Kepler's laws. For $N >= 3$ no general closed-form solution exists. Poincaré showed that the three-body problem has, in general, no analytic first integrals beyond the classical ones and that its dynamics can be chaotic @poincare1890probleme. Sundman's convergent series solution @sundman1913memoire exists in principle but converges far too slowly to be of practical use. Numerical integration is therefore the method of choice.

Integrating the equations of motion numerically for a short time is easy. Integrating them faithfully over many orbital periods is not. Every numerical method introduces a small error in each step, and *how* these errors accumulate matters more than how large they are. A method can have a tiny local error and still let the total energy drift steadily, so that a planet slowly spirals into its star. Another method with a larger local error can keep the energy error bounded forever. The difference comes from the *geometric structure* of Hamiltonian dynamics, and in particular from the symplecticity of the flow @hairer2006geometric @sanzserna1994numerical. Understanding this structure and choosing algorithms that respect it is the subject of _geometric numerical integration_ and the central theme of this thesis.

== Origin and scope of this work

This thesis is a revised and substantially extended version of a course project carried out in the course _Computational Physics_ (winter term 2024/25) together with Chahid El Aaouadi and Simon Pfeiffer @can2025computerphysik. The original project implemented a two-dimensional $N$-body simulation with a leapfrog integrator and a real-time visualisation, and it assessed the quality of the simulation through the ratio $|E_t \/ E_0|$ of the current to the initial total energy.

The present work keeps the physical problem but rebuilds the project from the ground up with three goals:

+ *Theory.* Give a self-contained derivation of the physics (Hamiltonian formulation, conservation laws, exact two-body solution, softening, chaos) and of the numerical methods, ordered by increasing sophistication, including their order, stability and structure-preserving properties.
+ *Implementation.* Provide a clean, modular and efficient implementation in which new integrators, force models and scenarios can be added without touching existing code, and in which every theoretical claim is backed by an automated test.
+ *Evaluation.* Quantify the behaviour of all methods with well-defined metrics (energy, momentum and angular-momentum errors, global errors against exact or high-precision references, empirical convergence orders, time-reversibility, physical observables, computational cost) on a set of benchmark problems of increasing difficulty.

== Contributions

Compared with the original course project, the main contributions are:

- a hierarchy of nine integrators, from explicit Euler to sixth-order symplectic composition methods and an adaptive Dormand--Prince pair, all behind a common interface and a registry;
- derivations of the order conditions, stability limits, modified energies and error-growth laws that explain the observed behaviour;
- vectorised force kernels that are up to two orders of magnitude faster than the original pair loop, and a leapfrog that reuses the last force evaluation and so needs half as many evaluations per step;
- a consistent energy diagnostic for softened gravity. The re-evaluation in @sec:eval-legacy shows that most of the energy fluctuation reported in the original project was an artefact of an inconsistent potential;
- an evaluation framework with nine reproducible experiments whose figures and tables are generated automatically and read directly by this document;
- an interactive viewer with two trajectory overlays: fading trails of the past motion, and forecasts of the future motion computed by integrating a "ghost" copy of the system.

== Outline

@ch:physics introduces the physics of the gravitational $N$-body problem, from Newton's law to the Hamiltonian formulation, the conservation laws and the exact solution of the Kepler problem. @ch:numerics develops the numerical methods step by step and analyses their accuracy, stability and geometric properties. @ch:implementation describes the software architecture and the most important implementation decisions. @ch:evaluation defines the evaluation metrics and presents the experiments. @ch:conclusion summarises the results and discusses possible extensions. The appendices contain a short user guide and the coefficients of all methods.

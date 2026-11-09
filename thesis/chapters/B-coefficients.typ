#import "../lib/utils.typ": *

= Method Coefficients <app:coefficients>

// Note: the appendix tables are not numbered figures; the template's
// per-chapter numbering would label them with the chapter index instead of "B".
#let app-caption(body) = align(center, block(above: 1.2em, below: 0.6em, text(font: "Roboto", size: 0.95em, body)))

== Composition weights

All composition methods apply the leapfrog step @eq:leapfrog with step sizes $gamma_k h$, $k = 1, dots, s$. The weights are palindromic and sum to one.

#app-caption[Weights of the symplectic composition methods.]
#data-table(
  columns: (auto, auto, 1fr),
  align: (left, center, left),
  header: ([Method], [$s$], [Weights $gamma_1, dots, gamma_s$]),
  [Yoshida 4 (triple jump of leapfrog)], [3],
    [$gamma_1 = gamma_3 = 1\/(2 - 2^(1\/3)) = 1.351207191959658$ \ $gamma_2 = -2^(1\/3)\/(2 - 2^(1\/3)) = -1.702414383919315$],
  [Yoshida 6 (triple jump of Yoshida 4)], [9],
    [$gamma_(3(i-1)+j) = beta_i alpha_j$, $alpha$ = Yoshida-4 weights, \ $beta_1 = beta_3 = 1\/(2 - 2^(1\/5)) = 1.174671758089363$, \ $beta_2 = -2^(1\/5)\/(2 - 2^(1\/5)) = -1.349343516178726$],
  [Yoshida 6, solution A @yoshida1990construction], [7],
    [$(w_3, w_2, w_1, w_0, w_1, w_2, w_3)$ with \ $w_1 = -1.17767998417887$, $w_2 = 0.235573213359357$, \ $w_3 = 0.784513610477560$, $w_0 = 1 - 2(w_1 + w_2 + w_3)$],
)

== Dormand--Prince 5(4)

The Butcher tableau of DOPRI5(4) @dormand1980family. The row $b$ gives the fifth-order solution, which is propagated; it equals the last row of $A$ (FSAL). The row $hat(b)$ gives the embedded fourth-order solution used for the error estimate.

#app-caption[Butcher tableau of the Dormand--Prince 5(4) pair.]
#{
  set text(size: 8.6pt)
  let f(a, b) = $#a\/#b$
  align(center, table(
    columns: 8,
    stroke: (x, y) => (if x == 0 { (right: 0.5pt) } else { (:) }) + (if y == 7 { (top: 0.5pt) } else { (:) }),
    inset: (x: 4.5pt, y: 4pt),
    align: right + horizon,
    $0$, [], [], [], [], [], [], [],
    f(1, 5), f(1, 5), [], [], [], [], [], [],
    f(3, 10), f(3, 40), f(9, 40), [], [], [], [], [],
    f(4, 5), f(44, 45), f(-56, 15), f(32, 9), [], [], [], [],
    f(8, 9), f(19372, 6561), f(-25360, 2187), f(64448, 6561), f(-212, 729), [], [], [],
    $1$, f(9017, 3168), f(-355, 33), f(46732, 5247), f(49, 176), f(-5103, 18656), [], [],
    $1$, f(35, 384), $0$, f(500, 1113), f(125, 192), f(-2187, 6784), f(11, 84), [],
    $b$, f(35, 384), $0$, f(500, 1113), f(125, 192), f(-2187, 6784), f(11, 84), $0$,
    $hat(b)$, f(5179, 57600), $0$, f(7571, 16695), f(393, 640), f(-92097, 339200), f(187, 2100), f(1, 40),
  ))
}

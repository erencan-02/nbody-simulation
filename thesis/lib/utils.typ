// Shared helpers for the thesis: colours, number formatting and data tables.

// Accent colour of the thesis (TU Darmstadt corporate design colour 1d).
#let accent = rgb("#243572")
#let accent-light = rgb("#5D85C3")
#let accent-bg = rgb("#243572").lighten(92%)

// ---------------------------------------------------------------------------
// Number formatting
// ---------------------------------------------------------------------------

// Parse a CSV cell (string) or pass numbers through.  Returns none for nan.
#let num(x) = {
  if type(x) == str {
    let s = lower(x.trim())
    if s == "nan" or s == "" or s == "inf" { none } else { float(x) }
  } else { x }
}

// Scientific notation, e.g. sci(0.000123) -> 1.2 x 10^-4.
#let sci(x, digits: 1) = {
  let v = num(x)
  if v == none { return [--] }
  if v == 0 { return $0$ }
  let e = calc.floor(calc.log(calc.abs(v), base: 10))
  let m = calc.round(v / calc.pow(10.0, e), digits: digits)
  if calc.abs(m) >= 10 {
    m = m / 10
    e = e + 1
  }
  let ms = str(m)
  if digits > 0 and not ms.contains(".") { ms = ms + "." + "0" * digits }
  if e == 0 { $#ms$ } else { $#ms dot 10^(#e)$ }
}

// Fixed-point with a given number of decimals (keeps trailing zeros).
#let fix(x, digits: 2) = {
  let v = num(x)
  if v == none { return [--] }
  let r = calc.round(v, digits: digits)
  let s = str(r)
  if digits > 0 {
    if not s.contains(".") { s = s + "." }
    let decimals = s.split(".").at(1).len()
    s = s + "0" * (digits - decimals)
  }
  $#s$
}

// Integer with thin-space thousands separators.
#let int-fmt(x) = {
  let v = num(x)
  if v == none { return [--] }
  let s = str(calc.round(v))
  let out = ""
  let n = s.len()
  for (i, c) in s.clusters().enumerate() {
    if i > 0 and calc.rem(n - i, 3) == 0 { out += sym.space.thin }
    out += c
  }
  out
}

#let yesno(x) = if x == "yes" [yes] else [no]

// ---------------------------------------------------------------------------
// Data access
// ---------------------------------------------------------------------------

// Load a CSV produced by the experiments as a list of dictionaries.
#let load-csv(path) = {
  let rows = csv(path)
  let header = rows.at(0)
  rows.slice(1).map(r => {
    let d = (:)
    for (i, k) in header.enumerate() { d.insert(k, r.at(i)) }
    d
  })
}

// Find the row with column `key` equal to `value`.
#let row-of(rows, key, value) = rows.find(r => r.at(key) == value)

// ---------------------------------------------------------------------------
// Layout helpers
// ---------------------------------------------------------------------------

// Booktabs-like table styling.
#let data-table(columns: auto, align: auto, header: (), ..rows) = {
  let ncols = header.len()
  table(
    columns: columns,
    align: align,
    stroke: none,
    inset: (x: 5pt, y: 3.2pt),
    table.hline(stroke: 0.9pt),
    table.header(..header.map(h => text(weight: "bold", h))),
    table.hline(stroke: 0.5pt),
    ..rows,
    table.hline(stroke: 0.9pt),
  )
}

// Highlighted box for definitions, key results and remarks.
#let callout(title: none, body) = block(
  width: 100%,
  fill: accent-bg,
  stroke: (left: 2.2pt + accent),
  inset: (x: 10pt, y: 8pt),
  radius: (right: 2pt),
  breakable: true,
)[
  #if title != none [#text(font: "Roboto", weight: "bold", fill: accent, size: 0.95em, title) \ ]
  #body
]

#let definition(title, body) = callout(title: [Definition (#title)], body)
#let result(title, body) = callout(title: title, body)

// Small-caps-like label for code identifiers.
#let code(s) = raw(s)

// Structure-Preserving Numerical Integration of the Gravitational N-Body Problem
//
// Build (from the repository root):
//   python scripts/fetch_fonts.py
//   typst compile thesis/main.typ --font-path thesis/assets/fonts
//
// Template: TU Darmstadt thesis template for Typst
//   https://typst.app/universe/package/athena-tu-darmstadt-thesis

#import "@preview/athena-tu-darmstadt-thesis:0.1.2": *
#import "lib/utils.typ": accent

// The template only knows "master" and "bachelor" theses; this document is a
// course-project thesis, so the type line on the title page is adapted.
#show "Master thesis by": "Computational physics project by"
#show "Date of submission:": "Revised version:"

#show: tudapub.with(
  title: [Structure-Preserving Integration of the Gravitational N-Body Problem],
  title_german: [From explicit Euler to high-order symplectic composition methods],
  title_height: 5.2em,
  author: "Eren Can",
  thesis_type: "master",
  // TU Darmstadt corporate design colour 1d (dark blue) instead of 9c (red).
  accentcolor: "1d",
  date_of_submission: datetime(year: 2026, month: 11, day: 9),
  location: "Darmstadt",
  reviewer_names: (),
  language: "en",
  // Stand-in emblem; replace by image("assets/logos/tuda_logo.svg") if the
  // official logo is available (see README).
  logo_tuda: image("assets/logos/project-emblem.svg"),
  logo_sub_content_text: [
    Course project \
    Computational Physics \
    Winter term 2024/25 \
    \
    Revised edition 2026
  ],
  abstract: include "chapters/00-abstract.typ",
  bib: bibliography("references.bib"),
  margin: tud_page_margin_big,
  show_pages: (
    title_page: true,
    outline_table_of_contents: true,
    thesis_statement_pursuant: false,
  ),
  additional_pages_after_outline_table_of_contents: include "chapters/00-notation.typ",
  reduce_heading_space_when_first_on_page: false,
  figure_numbering_per_chapter: true,
)

// --- document-wide styling ---------------------------------------------------
#show link: set text(fill: accent)
#show cite: set text(fill: accent)
#show ref: it => {
  if it.element == none { it } else { text(fill: accent, it) }
}
#set math.equation(numbering: "(1)")
#show raw.where(block: true): it => block(
  width: 100%,
  breakable: false,
  fill: rgb("#f4f6fa"),
  stroke: (left: 1.5pt + accent),
  inset: (x: 8pt, y: 6pt),
  text(size: 8.2pt, it),
)
#show raw.where(block: false): set text(size: 0.92em)
#set figure(gap: 0.8em)
#show table: set text(size: 8.6pt)
#show figure.where(kind: table): set figure.caption(position: top)
#set list(indent: 0.6em)
#set enum(indent: 0.6em)

#include "chapters/01-introduction.typ"
#include "chapters/02-physics.typ"
#include "chapters/03-numerical-methods.typ"
#include "chapters/04-implementation.typ"
#include "chapters/05-evaluation.typ"
#include "chapters/06-conclusion.typ"

// --- appendix ------------------------------------------------------------------
#counter(heading).update(0)
#set heading(numbering: "A.1")
#include "chapters/A-usage.typ"
#include "chapters/B-coefficients.typ"

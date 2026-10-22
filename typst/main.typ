// Typst version of main.tex; see README.md for building and required fonts.

#import "template.typ": *

// Equivalent of \ifenpc: `typst compile --input enpc=true ...`
#let enpc = sys.inputs.at("enpc", default: "false") == "true"

#show: book.with(title: "Lecture notes in Numerical Analysis", author: "Urbain Vaes")

#if enpc {
  panic("The ENPC title page and syllabus have not been converted yet")
}

#include "frontmatter/title_nyu.typ"

#show: frontmatter
#part("license", output: false)
#include "frontmatter/license.typ"
#part("syllabus_nyu", output: false)
#include "frontmatter/syllabus_nyu.typ"
#part("acknowledgements", output: false)
#include "frontmatter/acknowledgements.typ"
#part("notations_enpc")
#include "frontmatter/notations_enpc.typ"
#part("contents", output: false)
#table-of-contents()

#show: mainmatter
#part("introduction")
#include "mainmatter/introduction.typ"

#part("roundoff")
#include "mainmatter/roundoff.typ"

#part("interpolation_approximation")
#include "mainmatter/interpolation_approximation.typ"

#part("numerical_integration")
#include "mainmatter/numerical_integration.typ"

#part("linear_equations")
#include "mainmatter/linear_equations.typ"

#part("nonlinear_equations")
#include "mainmatter/nonlinear_equations.typ"

#part("eigenproblems")
#include "mainmatter/eigenproblems.typ"

#part("odes")
#include "mainmatter/odes.typ"

#show: appendix
#part("matrices")
#include "mainmatter/matrices.typ"

#part("julia")
#include "mainmatter/julia.typ"

#part("chebyshev")
#include "mainmatter/chebyshev.typ"

#references()

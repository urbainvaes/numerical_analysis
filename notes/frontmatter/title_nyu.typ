#import "../template.typ": *

// \maketitle of the book class, with the title, author and date of title_nyu.tex.
// Vertical spacing calibrated against the LaTeX output.
#page(header: none, footer: none, numbering: none, {
  set par(first-line-indent: 0pt, justify: false)
  set align(center)
  v(60pt - 2cm + 15pt)
  image("/figures/nyu-logo.jpg", width: 30%)
  v(1.5cm + 3pt)
  // Bold Latin Modern exists only up to the 12pt design size
  text(font: "Latin Modern Roman 12", size: LARGE, weight: "bold")[MATH-UA 9252: Numerical Analysis]
  v(1cm)
  image("/figures/chebychev_cover.pdf", width: 70%)
  v(3em)
  set text(font: serif-at(large), size: large)
  [Urbain #textsc[Vaes] \ #tt[urbain.vaes\@nyu.edu]]
  v(1.5em + 1cm - 5pt)
  [#textsc[NYU Paris], Fall term 2026]
  v(2cm - 5pt)
  v(1fr)
  set align(left)
  [*Weekly schedule*:]
  list(
    [Lectures on Monday and Wednesday from 15:30 to 16:45 in room 406;],
    [Recitation on Monday and Wednesday from 16:50 to 17:35 in room 406;],
    [Office hour on Monday from 14:30 to 15:30 in room 606.],
  )
  v(25pt)
})

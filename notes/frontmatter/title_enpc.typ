#import "../template.typ": *

// \maketitle of the book class, with the title, author and date of title_enpc.tex.
// Spacing follows title_nyu.typ.
#page(header: none, footer: none, numbering: none, {
  set par(first-line-indent: 0pt, justify: false)
  set align(center)
  v(60pt - 2cm + 15pt)
  // \includegraphics[width=.25\textwidth,trim=50 50 50 50,clip]: the image is
  // 707.25bp × 900bp (943 × 1200 px at 96 dpi), with 50bp cropped on each side
  layout(size => {
    let w = 0.25 * size.width
    let scale = w / 607.25
    box(width: w, height: 800 * scale, clip: true,
      move(dx: -50 * scale, dy: -50 * scale,
        image("/figures/logo_enpc.jpg", width: 707.25 * scale)))
  })
  v(1cm)
  // Bold Latin Modern exists only up to the 12pt design size
  text(font: "Latin Modern Roman 12", size: LARGE, weight: "bold")[Practical scientific computing]
  v(.3cm)
  image("/figures/chebychev_cover.pdf", width: 70%)
  v(3em)
  set text(font: serif-at(large), size: large)
  [Urbain #textsc[Vaes] \ #tt[urbain.vaes\@inria.fr]]
  v(1.5em + .4cm - 5pt)
  [#textsc[École des Ponts], first year]
  v(.5cm)
  v(1fr)
  set text(font: serif, size: normalsize)
  rect(inset: 1em, width: 70% + 2em, {
    [*Schedule*: Tuesday from 8:30 to 11:15 on the following dates:]
    v(.5cm)
    set align(left)
    grid(columns: (20%, 20%, 20%), column-gutter: 0pt,
      list[8 Apr][15 Apr], list[29 Apr][6 May], list[13 May][20 May])
    v(.5cm)
    set align(center)
    [*Exam* on Tuesday 27 May]
  })
  v(25pt)
})

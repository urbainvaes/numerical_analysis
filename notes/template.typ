// Typst reproduction of preamble.tex and macros.tex.

// Icons (\faGears, \faLaptop, \faSearchPlus); needs the Font Awesome 7 fonts
// (Arch package otf-font-awesome)
#import "@preview/fontawesome:0.6.0": fa-icon

// Colors (preamble.tex)
#let lightblue = rgb("#f1f8fb")
#let lightgreen = rgb("#f1fbf4")
#let lightindigo = rgb("#f3f1fb")
#let darkblue = rgb("#315a88")
#let darkgrey = rgb(120, 120, 120)
#let lightred = rgb(100%, 70%, 71%)

// Fonts: the Latin Modern fonts used by lualatex, installed as system fonts
// (Arch packages otf-latin-modern and otf-latinmodern-math); see README.md
#let serif = "Latin Modern Roman"
#let heading-font = "Latin Modern Sans"
#let mono = "Latin Modern Mono"
#let mono-scale = 1

// LaTeX sizes for the 11pt class option
#let footnotesize = 9pt
#let small = 10pt
#let normalsize = 10.95pt
#let large = 12pt
#let Large = 14.4pt
#let LARGE = 17.28pt
#let huge = 20.74pt
#let Huge = 24.88pt

// \onehalfspacing at 11pt: \baselineskip = 1.213 × 13.6pt
#let baselineskip = 1.213 * 13.6pt
#let leading = baselineskip - normalsize
#let parindent = 17pt
#let ex = 4.72pt

// Optical size of Latin Modern Roman, as selected by the LaTeX font definitions
#let serif-at(size) = if size >= 15.5pt { "Latin Modern Roman 17" } else if size >= 11pt { "Latin Modern Roman 12" } else { serif }

// Which part of the book we are in: "front", "main" or "appendix"
#let matter = state("matter", "front")

// ---------------------------------------------------------------------------
// Inline helpers

#let tt(body) = text(font: mono, size: mono-scale * 1em, body)
// Latin Modern has no OpenType small caps, but a separate Caps font
#let textsc(body) = text(font: "Latin Modern Roman Caps", body)
// Indentation of a paragraph that follows an environment (a blank line after
// \end{...} in LaTeX); Typst only indents paragraphs that follow a paragraph
#let indent = h(parindent)
// \verb: monospace text without syntax highlighting
#let verb(code) = raw(code, lang: "txt")
// \moreinfo: marks sections that are for information only
#let moreinfo = super(typographic: false, text(fill: lightred, fa-icon("magnifying-glass-plus", solid: true)))

// \paragraph{...}: run-in heading
#let paragraph(title, body) = block(above: 3.25 * ex + leading, below: leading,
  [#text(font: heading-font, weight: "bold", title)#h(1em)#body])

// ---------------------------------------------------------------------------
// Math macros (macros.tex)

#let vect(x) = math.bold(x)
#let mat(x) = math.upright(math.sans(x))
#let real = math.upright(math.bold("R"))
#let complex = math.upright(math.bold("C"))
#let nat = math.upright(math.bold("N"))
#let integer = math.upright(math.bold("Z"))
#let floating = math.upright(math.bold("F"))
#let poly = math.upright(math.bold("P"))
#let expect = math.upright(math.bold("E"))
#let bigo = math.cal("O")
#let placeholder = math.class("normal", text(fill: luma(66%), $bullet$))
// \abs and \norm (\DeclarePairedDelimiter, unstarred): delimiters of fixed size
#let abs(x) = $class("opening", |) #x class("closing", |)$
#let norm(x) = $class("opening", ||) #x class("closing", ||)$
// \left ... \right: lr(...) scales the delimiters to the content. Unlike in TeX,
// Typst also scales matching delimiters without lr; this is disabled below.
#let lr = math.lr.with(size: 100%)
#let ip(..args) = $chevron.l #args.pos().join($,$) chevron.r$
// Machine arithmetic operations \madd, \msub, \mtimes, \mdiv
#let madd = math.class("binary", math.hat(sym.plus))
#let msub = math.class("binary", math.hat(sym.minus))
#let mtimes = math.class("binary", math.hat(sym.ast))
#let mdiv = math.class("binary", math.hat(sym.slash))
#let fl = math.op("fl")
// \dots and \cdots: dots separated by thin spaces, as in TeX
#let dots = math.class("normal", $thin . thin . thin . thin$)
#let cdots = math.class("normal", $dot.op thin dot.op thin dot.op$)
#let ee = math.upright("e")
// \backslash (set difference): an ordinary symbol, without space around it
#let backslash = math.class("normal", sym.backslash)
// \d, \D and \i: upright d (differentials), D and i
#let dd = math.upright("d")
#let DD = math.upright("D")
#let ii = math.upright("i")
#let vecnorm = norm
#let var = math.upright(math.bold("V"))
#let proba = math.upright(math.bold("P"))
#let laplacian = sym.triangle.stroked.t
// \DeclareMathOperator
#let diag = math.op("diag")
#let cond = math.op("cond")
#let Span = math.op("Span")
#let arccosh = math.op("arccosh")
#let sign = math.op("sign")
#let col = math.op("col")
#let spectrum = math.op("spectrum")
#let trace = math.op("tr", limits: true)
#let argmin = math.op("arg" + h(0.17em) + "min", limits: true)
#let argmax = math.op("arg" + h(0.17em) + "max", limits: true)
// Matrices with parentheses (pmatrix). Note that `mat` is \mat (sans serif) here.
#let pmat = math.mat.with(delim: "(")
#let bmat = math.mat.with(delim: "[")
// \bigl( ... \bigr) and friends: big((x)), or Big(\{) ... Big(\}) across lines
#let big(body) = math.lr(body, size: 1.4em)
#let Big(body) = math.lr(body, size: 1.85em)
#let bigg(body) = math.lr(body, size: 2.3em)
// \colon: punctuation, without space before
#let colon = $class("punctuation", :) #h(1em / 3)$

// Number of the chapter containing a location: "1", or "A" in the appendix
#let chapter-number(loc) = {
  let style = if matter.at(loc) == "appendix" { "A" } else { "1" }
  numbering(style, counter(heading).at(loc).first())
}

// Figure kinds numbered per chapter, e.g. "Figure 1.2" or "Exercise 1.4"
#let theorem-kinds = ("theorem", "definition", "notation", "example", "remark", "exercise", "task")
#let per-chapter-kinds = (image, table, "theorem", "definition", "notation", "example", "remark", "exercise")

// ---------------------------------------------------------------------------
// Headings

// Running header: "Chapter 1. Title", or the title of an unnumbered chapter.
// Pages on which a chapter starts have no header (LaTeX "plain" page style).
#let running-header() = context {
  let here-page = here().page()
  let starts-here = query(heading.where(level: 1)).any(h => h.location().page() == here-page)
  if starts-here { return }
  let previous = query(selector(heading.where(level: 1)).before(here()))
  if previous.len() == 0 { return }
  let h = previous.last()
  let title = if h.numbering != none {
    [#h.supplement #numbering(h.numbering, ..counter(heading).at(h.location()).slice(0, 1)).~~#h.body]
  } else { h.body }
  set text(fill: darkgrey, style: "italic")
  align(right, title)
}

#let chapter-heading(it) = {
  pagebreak(weak: true)
  set par(first-line-indent: 0pt, justify: false, leading: 0.5em)
  set text(font: heading-font, weight: "bold")
  // Spacing calibrated against the LaTeX output (titlesec: 50pt before, 40pt after)
  block(above: 0pt, below: 48.5pt, {
    counter(math.equation).update(0)
    for kind in per-chapter-kinds { counter(figure.where(kind: kind)).update(0) }
    if it.numbering != none {
      v(72pt)
      text(size: huge, [#it.supplement #counter(heading).display(it.numbering)])
      v(30pt)
    } else {
      v(74pt)
    }
    text(size: Huge, it.body)
  })
}

#let section-heading(it) = {
  set par(first-line-indent: 0pt, justify: false)
  set text(font: heading-font, weight: "bold", size: Large, fill: darkblue)
  let extra = 1.213 * 18pt - baselineskip
  block(above: 3.5 * ex + leading - extra, below: 2.3 * ex + leading, inset: (top: extra), sticky: true, {
    if it.numbering != none {
      counter(heading).display(it.numbering)
      h(1em)
    }
    it.body
  })
}

#let subsubsection-heading(it) = {
  set par(first-line-indent: 0pt, justify: false)
  set text(font: heading-font, weight: "bold", size: normalsize)
  block(above: 3.25 * ex + leading, below: 1.5 * ex + leading, sticky: true, it.body)
}

#let subsection-heading(it) = {
  set par(first-line-indent: 0pt, justify: false)
  set text(font: heading-font, weight: "bold", size: large)
  block(above: 3.25 * ex + leading, below: 1.5 * ex + leading, sticky: true, {
    if it.numbering != none {
      counter(heading).display(it.numbering)
      h(1em)
    }
    it.body
  })
}

// \chapter* and \section*: unnumbered and absent from the table of contents
#let chapter-star(title) = heading(level: 1, numbering: none, outlined: false, title)
#let section-star(title) = heading(level: 2, numbering: none, outlined: false, title)
#let subsection-star(title) = heading(level: 3, numbering: none, outlined: false, title)

// ---------------------------------------------------------------------------
// Table of contents (book class, tocdepth = 1, linktocpage)

#let toc-entry(it) = {
  let loc = it.element.location()
  let page-number = link(loc, text(fill: darkblue, it.page()))
  if it.level == 1 {
    block(above: 1em + leading, below: leading, strong({
      if it.element.numbering != none { box(width: 1.5em, it.prefix()) }
      it.body()
      h(1fr)
      page-number
    }))
  } else {
    let dots = box(width: 1fr, repeat(gap: 0.45em, justify: false)[.])
    block(spacing: leading, {
      h(1.5em)
      box(width: 2.3em, it.prefix())
      it.body()
      h(0.5em)
      dots
      box(width: 2.5em, align(right, page-number))
    })
  }
}

#let table-of-contents() = {
  chapter-star[Contents]
  v(1em)
  outline(title: none, depth: 2)
}

// \minitoc: sections and subsections of the current chapter
#let minitoc-entry(it) = {
  let loc = it.element.location()
  let dots = box(width: 1fr, repeat(gap: 0.45em, justify: false)[.])
  block(spacing: 7.4pt, {
    if it.level == 2 {
      strong[#box(width: 2.3em, it.prefix())#it.body()#h(0.5em)#dots#box(width: 2.2em, align(right, link(loc, text(fill: darkblue, it.page()))))]
    } else {
      h(2.3em)
      [#box(width: 3.2em, it.prefix())#it.body()#h(0.5em)#dots#box(width: 2.2em, align(right, link(loc, text(fill: darkblue, it.page()))))]
    }
  })
}

#let minitoc() = context {
  let start = here()
  let next = query(selector(heading.where(level: 1)).after(start))
  let target = selector(heading.where(level: 2)).or(heading.where(level: 3)).after(start)
  if next.len() > 0 { target = target.before(next.first().location(), inclusive: false) }
  // Spacing calibrated against the LaTeX output
  block(above: 0pt, below: 36.5pt, inset: (top: 17pt), pad(left: 3.8em, right: 2em, {
    set text(size: 10pt)
    show outline.entry: minitoc-entry
    outline(title: none, target: target)
  }))
}

// ---------------------------------------------------------------------------
// Per-chapter PDFs: `part("name")` before a chapter starts it on a new page and
// records that page; split.py then extracts build/typst/<name>.pdf from the book.

#let part(name, output: true) = {
  pagebreak(weak: true)
  context [#metadata((name: name, page: here().page(), output: output)) <part>]
}

// ---------------------------------------------------------------------------
// Parts of the book

#let frontmatter(body) = {
  matter.update("front")
  set heading(numbering: none)
  set page(numbering: "i")
  counter(page).update(1)
  body
}

#let mainmatter(body) = {
  pagebreak(weak: true)
  matter.update("main")
  set heading(numbering: "1.1")
  set page(numbering: "1")
  counter(page).update(1)
  counter(heading).update(0)
  body
}

#let appendix(body) = {
  pagebreak(weak: true)
  matter.update("appendix")
  set heading(numbering: "A.1")
  show heading.where(level: 1): set heading(supplement: [Appendix])
  counter(heading).update(0)
  body
}

// Equation numbers of the form (chapter.number)
#let equation-numbering(n) = context [(#chapter-number(here()).#n)]

// Display equations are unnumbered, like \[ ... \]; numbered ones, like the
// equation environment, are written #numbered($ ... $) <label>. With eq-last,
// the number is on the last line, as in an align environment with \notag on
// the other lines.
#let numbered(eq) = math.equation(block: true, numbering: equation-numbering, eq.body)
#let eq-last(eq) = math.equation(block: true, numbering: equation-numbering, number-align: end + bottom, eq.body)

// References: \eqref-style "(1.2)" for equations, "Exercise 1.4" for figures and
// environments, colored links otherwise. Numbers are computed at the location of
// the referenced element, not of the reference.
#let reference(it) = {
  let el = it.element
  if el == none { return it }
  let loc = el.location()
  if el.func() == math.equation {
    let n = counter(math.equation).at(loc).first()
    link(loc, text(fill: darkblue, style: "normal")[(#chapter-number(loc).#n)])
  } else if el.func() == figure {
    let n = el.counter.at(loc).first()
    let number = if el.kind in per-chapter-kinds [#chapter-number(loc).#n] else [#n]
    let supplement = if it.supplement == auto { el.supplement } else { it.supplement }
    link(loc, text(fill: darkblue)[#supplement #number])
  } else {
    text(fill: darkblue, it)
  }
}

// ---------------------------------------------------------------------------
// Theorem-like environments (amsthm styles in mdframed boxes)

// Floats are single-spaced (setspace)
#let float-leading = 13.6pt - normalsize

// Space around environments (\topsep), and around and inside the mdframed boxes;
// calibrated against the LaTeX output
#let env-spacing = leading + 9pt
#let box-above = leading + 6pt
#let box-below = leading
#let box-inset = (x: 0.5em, top: 0.3em + 1.5pt, bottom: 0.6em)

#let env-box(kind, name, style, fill, icon, title, body) = {
  let number = context {
    let n = counter(figure.where(kind: kind)).get().first()
    if kind in per-chapter-kinds [#chapter-number(here()).#n] else [#n]
  }
  let note = if title != none [ (#title)]
  let head = if style == "remark" {
    [#emph(name) #number#note#emph[.]]
  } else {
    [#if icon != none [#icon~]#strong[#name #number]#note#strong[.]]
  }
  let content = {
    set par(first-line-indent: (amount: parindent, all: false))
    set align(left)
    [#text(style: "normal", head)#h(5pt)]
    if style == "plain" { text(style: "italic", body) } else { body }
  }
  if fill == none {
    block(width: 100%, above: env-spacing, below: env-spacing, breakable: true, content)
  } else {
    block(width: 100%, above: box-above, below: box-below, breakable: true,
      fill: fill, inset: box-inset, content)
  }
}

#let environment(kind, name, style: "plain", fill: none, icon: none) = (title: none, body) => figure(
  kind: kind, supplement: name, numbering: "1", outlined: false,
  env-box(kind, name, style, fill, icon, title, body),
)

#let theorem = environment("theorem", [Theorem], fill: lightblue)
#let proposition = environment("theorem", [Proposition], fill: lightblue)
#let lemma = environment("theorem", [Lemma], fill: lightblue)
#let corollary = environment("theorem", [Corollary], fill: lightblue)
#let definition = environment("definition", [Definition], style: "definition", fill: lightblue)
#let notation = environment("notation", [Notation], style: "definition", fill: lightblue)
#let example = environment("example", [Example], style: "remark", fill: lightindigo)
#let remark = environment("remark", [Remark], style: "remark", fill: lightindigo)
#let exercise = environment("exercise", [Exercise], icon: fa-icon("gears", solid: true))
#let compexercise = environment("exercise", [Exercise], icon: fa-icon("laptop", solid: true))
// Tasks are numbered throughout the book, not per chapter
#let task = environment("task", [Task], icon: fa-icon("laptop", solid: true))

// \begin{proof}[Sketch of the proof] -> #proof(title: [Sketch of the proof])[...]
#let proof(title: [Proof], body) = block(width: 100%, above: leading + 6pt, below: leading + 6pt, breakable: true,
  [#emph[#title.]#h(0.5em)#body#h(1fr)#h(1em)$square.stroked$])

#let solution(body) = block(width: 100%, above: box-above, below: box-below,
  breakable: true, fill: lightgreen, inset: box-inset, {
    // \small
    set text(size: 10pt)
    set par(leading: 1.213 * 12pt - 10pt, spacing: 1.213 * 12pt - 10pt)
    [#emph[Solution.]#h(5pt)#body#h(1fr)#h(1em)$triangle.stroked.t$]
  })

// ---------------------------------------------------------------------------
// Algorithms (algorithm with the "ruled" float style, and algpseudocode)
//
// Each line is a pair (indentation level, content), or a triple with a label
// that can be referenced with #lineref(<label>).

#let algorithm(caption: none, ..lines) = figure(kind: "algorithm", supplement: [Algorithm],
  numbering: "1", outlined: false, caption: caption, {
    let cells = ()
    for (i, line) in lines.pos().enumerate() {
      let (level, body, ..rest) = line
      let number = [#text(size: 9pt)[#(i + 1):]]
      if rest.len() > 0 { number = [#number#metadata((algorithm-line: i + 1))#label(rest.first())] }
      cells.push(align(right, number))
      cells.push(pad(left: level * 1.5em, body))
    }
    grid(columns: (1.2em, 1fr), column-gutter: 0.5em, row-gutter: float-leading, ..cells)
  })

#let lineref(lbl) = context {
  let line = query(lbl).first()
  link(line.location(), text(fill: darkblue)[line #line.value.algorithm-line])
}

#let algorithm-figure(it) = block(width: 100%, above: env-spacing, below: env-spacing, {
  set align(left)
  set par(first-line-indent: 0pt)
  set block(spacing: 0pt)
  block(width: 100%, stroke: (top: 0.8pt, bottom: 0.4pt), inset: (y: 3pt),
    [*#it.supplement #it.counter.display()* #it.caption.body])
  block(width: 100%, stroke: (bottom: 0.4pt), inset: (y: 3pt), it.body)
})

// ---------------------------------------------------------------------------
// Figures and tables

// "Figure 1.2: ...", centered if it fits on one line (caption package defaults)
#let figure-caption(it) = context {
  let n = it.counter.get().first()
  let number = if it.kind in per-chapter-kinds [#chapter-number(here()).#n] else [#n]
  let caption = [#it.supplement #number#it.separator#it.body]
  layout(size => {
    if measure(caption).width <= size.width { align(center, caption) } else { align(left, caption) }
  })
}

// Julia listing with a caption (\begin{code} ... \end{code}), numbered throughout
// the book: #listing(caption: [...])[```julia ... ```] <label>
#let listing(caption: none, body) = figure(kind: raw, supplement: [Listing], caption: caption, body)

// LaTeX table in the style |c|c|...| with \hline between all rows. As with
// \arraystretch = 1, rows are one line high: struts of 0.7 and 0.3 times the line height.
#let tabular(stretch: 1, ..args) = context {
  let height = stretch * (par.leading + normalsize)
  set text(top-edge: 0.7 * height, bottom-edge: -0.3 * height)
  table(stroke: 0.4pt, align: center + horizon, inset: (x: 6pt, y: 0pt), ..args)
}

// ---------------------------------------------------------------------------
// Main template

#let book(title: none, author: none, body) = {
  set document(title: title, author: author)
  set page(
    paper: "a4",
    margin: 1in,
    header-ascent: 25pt,
    footer-descent: 30pt - normalsize,
    header: running-header(),
    footer: context {
      set align(center)
      counter(page).display()
    },
  )
  set text(font: serif, size: normalsize, lang: "en", top-edge: 1em, bottom-edge: "baseline")
  set par(justify: true, leading: leading, spacing: leading, first-line-indent: parindent)

  // Lists: \leftmargini = 2.5em, \labelsep = .5em, \itemsep + \parsep = 9pt, \topsep = 9pt
  // Second level: \leftmarginii = 2.2em, label \textendash
  set list(indent: 0pt, body-indent: 0.5em, spacing: leading + 9pt,
    marker: (box(width: 2em, align(right)[•]), box(width: 1.7em, align(right)[–])))
  show list: set block(above: leading + 9pt, below: leading + 9pt)
  // Paragraphs within items: \listparindent = 0pt and \parsep = 4.5pt
  show list: set par(first-line-indent: 0pt, spacing: leading + 4.5pt)
  show enum: set par(first-line-indent: 0pt, spacing: leading + 4.5pt)
  set enum(indent: 0pt, body-indent: 0.5em, spacing: leading + 9pt, numbering: n => box(width: 2em, align(right)[#n.]))
  show enum: set block(above: leading + 9pt, below: leading + 9pt)

  // Math
  show math.equation: set text(font: "Latin Modern Math")
  show math.equation: set math.lr(size: 1em)
  set math.mat(column-gap: 10pt)
  // \abovedisplayskip and \belowdisplayskip, calibrated against the LaTeX output
  show math.equation.where(block: true): set block(above: 16.3pt, below: 16.3pt)
  // TeX places sub- and superscripts of a closing parenthesis as for a letter,
  // whereas Typst shifts them by the depth and height of the parenthesis
  show math.attach: it => {
    if it.base.func() != math.lr or "size" in it.base.fields() or (it.b == none and it.t == none) { return it }
    let fields = it.fields()
    let base = fields.remove("base")
    base + math.attach(sym.zws, ..fields)
  }
  set math.equation(supplement: none)
  // \leq and \geq are \leqslant and \geqslant (macros.tex). In the LaTeX version these
  // come from the AMS symbol font, which is closer to New Computer Modern Math than to
  // the bolder glyphs of Latin Modern Math.
  show "≤": text(font: "New Computer Modern Math", "⩽")
  show "≥": text(font: "New Computer Modern Math", "⩾")

  // Raw text (Typst scales it by 0.8 by default): Julia code, colored as with minted
  set raw(lang: "julia", theme: "julia.tmTheme")
  show raw: set text(font: mono, size: mono-scale / 0.8 * 1em, style: "normal")
  show raw.where(block: true): set par(leading: baselineskip - mono-scale * normalsize, justify: false)
  show raw.where(block: true): set block(above: env-spacing, below: env-spacing)

  // Figures, tables and environments
  set figure(gap: 10pt)
  set figure.caption(separator: [: ])
  show figure.caption: figure-caption
  show figure.where(kind: image): set block(above: env-spacing, below: env-spacing)
  show figure.where(kind: image): set par(leading: float-leading, spacing: float-leading)
  show figure.where(kind: table): set par(leading: float-leading, spacing: float-leading)
  show figure.where(kind: "algorithm"): set par(leading: float-leading, spacing: float-leading)
  show figure.where(kind: raw): set block(above: env-spacing, below: env-spacing)
  show figure.where(kind: raw): set par(leading: float-leading, spacing: float-leading)
  show figure.where(kind: raw): it => {
    // \footnotesize code, left-aligned, caption below
    show raw.where(block: true): set text(size: footnotesize)
    show raw.where(block: true): set align(left)
    it
  }
  show figure.where(kind: table): set block(above: env-spacing, below: env-spacing)
  set place(clearance: 20pt)
  // Environments are figures, so that they are numbered and can be referenced
  show figure: it => if it.kind in theorem-kinds { it.body } else { it }
  show figure.where(kind: "algorithm"): algorithm-figure

  // Links and references (hyperref colorlinks)
  show link: set text(fill: darkblue)
  show ref: reference
  show smallcaps: it => textsc(it.body)
  // In citations, only the numbers are colored, as with biblatex and hyperref
  show cite: it => {
    if it.form != "normal" or it.style != auto { return it }
    let number = text(fill: darkblue, cite(it.key, style: "csl/number.csl"))
    [\[#number#if it.supplement != none [, #it.supplement]\]]
  }
  show heading.where(level: 1): set heading(supplement: [Chapter])
  show heading.where(level: 2): set heading(supplement: [Section])

  show heading.where(level: 1): chapter-heading
  show heading.where(level: 2): section-heading
  show heading.where(level: 3): subsection-heading
  show heading.where(level: 4): subsubsection-heading
  show heading.where(level: 4): set heading(numbering: none, outlined: false)

  show outline.entry: toc-entry

  body
}

// Bibliography: biblatex style trad-abbrv, approximated with a CSL style
#let references() = {
  part("references")
  show bibliography: set heading(outlined: false)
  // \bibfont = \small (10pt), and DOI/URL lines in \scriptsize (8pt)
  set text(size: 10pt)
  set par(first-line-indent: 0pt, leading: 1.213 * 12pt - 10pt, spacing: 1.213 * 12pt - 10pt + 4.5pt)
  show link: set text(font: mono, size: mono-scale * 8pt)
  show regex("^(doi|url)$"): set text(size: 8pt)
  // The line before a DOI or URL is not justified, as with \newline
  show "\u{2028}": linebreak()
  // Labels right-aligned, with the label width and separation of biblatex
  show regex("^\\[[0-9]+\\]$"): it => box(width: 1.6em, align(right, it)) + h(0.5em)
  // Adapt references.bib to Hayagriva: drop empty author fields, which would otherwise
  // be sorted first, and protect "and" in publishers, which biblatex treats as a list
  let bib = read("/references.bib")
  bib = bib.replace(regex("(?i)\bauthor\s*=\s*\{\s*\},?"), "")
  bib = bib.replace(regex("(?i)(\bpublisher\s*=\s*)\{([^{}]*)\}"), m => m.captures.at(0) + "{{" + m.captures.at(1) + "}}")
  bibliography(bytes(bib), title: [Bibliography], style: "csl/trad-abbrv.csl", full: true)
}

// \fullcite: same style, without DOI, URL and edition
#let fullcite(key) = cite(key, form: "full", style: "csl/fullcite.csl")

# Converting LaTeX chapters to Typst

Conventions used in `mainmatter/*.typ`. `mainmatter/roundoff.typ` is a complete
example; all helpers below are defined in `template.typ`, imported with
`#import "../template.typ": *`.

## Structure

| LaTeX                               | Typst                                          |
|-------------------------------------|------------------------------------------------|
| `\chapter{T}\label{cha:x}\minitoc`  | `= T <cha:x>` then `#minitoc()`                |
| `\section{T}\label{sec:x}`          | `== T <sec:x>`                                 |
| `\subsection{T}`                    | `=== T`                                        |
| `\subsubsection*{T}`                | `==== T` (never numbered)                      |
| `\section*{T}`                      | `#section-star[T]`                             |
| `\subsection*{T}`                   | `#subsection-star[T]`                          |
| `\paragraph{T} text...`             | `#paragraph[T][text...]` (body may contain displays) |
| `\moreinfo`                         | `#moreinfo`                                    |
| `\cref{x}`, `\Cref{x}`, `\eqref{x}` | `@x` (labels keep their LaTeX names; a label with spaces: `#ref(label("cha:a b"))`) |
| `\cite{k}`, `\cite[Ch. 3]{k}`       | `@k`, `@k[Ch. 3]`                              |
| `\cite{a,b}`                        | `@a @b` (each is rendered as its own `[n]`)    |
| `\fullcite{k}`                      | `#fullcite(<k>)`                               |
| `\url{u}`, `\href{u}{t}`            | `#link("u")`, `#link("u")[t]`                  |
| `\footnote{t}`                      | `#footnote[t]`                                 |

## Text

- `\emph{x}` → `_x_`, `\textbf{x}` → `*x*`, `\texttt{x}` → `#tt[x]`,
  `\textsc{x}` → `#textsc[x]`, `` ``x'' `` → `"x"`, `--` stays `--`,
  `~` stays `~`, `\%` → `%`, `i.e.\ ` → `i.e.~`.
- Escape Typst markup characters in text: `#`, `$`, `*`, `_`, `@`, `<`, `` ` ``,
  `\`; a line must not start with `- `, `+ `, `= ` or `1. ` unless it is a list
  or heading.
- `\julia{code}` and `\mintinline{julia}{code}` → `` `code` `` (Julia is the
  default language of raw text); `\verb?x?` → `#verb("x")`.
- `minted` blocks → ```` ```julia ```` blocks (flush left, as Typst dedents).
- `\begin{code}...\end{code}` with a caption → `#listing(caption: [...])[```julia ... ```] <label>`.
- Paragraph indentation: Typst only indents a paragraph that follows another
  paragraph. Where LaTeX has a blank line after `\end{env}` (theorem-like
  environments, proofs, solutions, algorithms, code) followed by text, start that
  paragraph with `#indent `. The same after a list (`\end{itemize}`) followed by a
  blank line. Not after headings or displays.
- `\small`, `\footnotesize` → `#text(size: small)[...]`, `#text(size: footnotesize)[...]`.
- `\vspace`, `\hspace` → `#v(...)`, `#h(...)`, only where they matter visually.

## Environments

`#theorem`, `#proposition`, `#lemma`, `#corollary`, `#definition`,
`#notation`, `#example`, `#remark`, `#exercise`, `#compexercise`, `#task`,
`#proof`, `#solution`:

```typst
#proposition(title: [Optional note])[
  Body.
] <prop:label>
```

The label goes right after the closing bracket.
`\begin{proof}[Sketch of the proof]` → `#proof(title: [Sketch of the proof])[...]`.

Lists: `itemize` → `- item`, `enumerate` → `+ item`; continuation lines of an
item are indented by two spaces; a blank line between items is allowed.

## Math

- Inline `$...$`/`\(...\)` → `$...$` with no space after the opening `$`.
- `\[...\]`, `equation*`, `align*` → `$ ... $` with spaces inside: unnumbered.
- `equation` → `#numbered($ ... $) <eq:label>`.
- `align` with a single number on the last line (other lines `\notag`/`\nonumber`)
  → `#eq-last($ ... $) <eq:label>`. An `align` with numbers on several lines
  is split into several numbered equations, or numbered once with `#eq-last`;
  report such cases.
- Alignment: `&` and line breaks `\`, at the top level of the equation.
- `\left( ... \right)` → `lr(( ... ))`, `\left\lvert ... \right\rvert` → `lr(|...|)`.
  Plain parentheses are NOT auto-scaled (as in TeX).
- `\bigl( ... \bigr)` → `big((...))`; `\Bigl\{` ... `\Bigr\}` across lines →
  `Big(\{)` ... `Big(\})`; also `bigg`.
- Spacing: `\,` → `thin`, `\quad` → `quad`, `\qquad` → `wide`, `\!` → `#h(-0.17em)`.
- Do not put a space between a symbol and its parenthesized argument when the
  LaTeX has none visible: `f(x)`, `fl(x)`, `floating(p, E)`; TeX ignores spaces,
  Typst does not always.
- Exception: a single-letter subscript or superscript followed by parentheses
  absorbs them: `T_n(x)` is T with subscript "n(x)". Write `T_n (x)`.
  (`fl_16(x)` and `x_(n+1)(y)` are fine.)
- A `~` in LaTeX math is a space: `#sym.space.nobreak` (or just a space), not `~`.
- Inside `pmat(...)` and `cases(...)`, a comma separates cells: write `\,` for a
  literal comma.
- Macros (from `macros.tex`): `\vect x` → `vect(x)`, `\mat A` → `mat(A)`,
  `\real` `\complex` `\nat` `\integer` `\floating` `\poly` `\expect` `\var` `\proba`
  → `real` `complex` `nat` `integer` `floating` `poly` `expect` `var` `proba`,
  `\bigo` → `bigo`, `\abs{x}` → `abs(x)`, `\norm{x}` → `norm(x)`,
  `\vecnorm{x}` → `vecnorm(x)`, `\ip{x,y}` → `ip(x, y)`, `\placeholder` →
  `placeholder`, `\madd` `\msub` `\mtimes` `\mdiv` → `madd` `msub` `mtimes` `mdiv`,
  `\d` → `dd`, `\D` → `DD`, `\i` → `ii`, `\e` → `ee`, `\t` (transpose) → `T`,
  `\laplacian` → `laplacian`, `{\rm fl}` → `fl`.
- Operators: `diag`, `cond`, `Span`, `arccosh`, `sign`, `col`, `spectrum`,
  `trace`, `argmin`, `argmax`; built-in: `sin`, `cos`, `exp`, `log`, `ln`, `det`,
  `dim`, `sup`, `inf`, `lim`, `limsup`, `liminf`, `max`, `min`, `Re`, `arccos`, ...
  Others: `op("name")`.
- `\begin{pmatrix} a & b \\ c & d \end{pmatrix}` → `pmat(a, b; c, d)`
  (NOT `mat`, which is `\mat`); `bmatrix` → `bmat(...)`; `cases` → `cases(...)`.
  Dots in matrices: `sym.dots.h`, `sym.dots.v`, `sym.dots.down` (`dots` is
  redefined by the template).
- Greek: `\varepsilon` → `epsilon`, `\epsilon` → `epsilon.alt`, `\varphi` →
  `phi`, `\phi` → `phi.alt`, `\ell` → `ell`.
- `\dots`, `\ldots`, `\dotsc`, `\hdots` → `dots`; `\cdots`, `\dotsb` → `cdots`.
- `\mathrm{x}`/`\text{x}` → `"x"` (or `upright(x)`); `\mathbf` → `bold(upright(..))`;
  `\mathcal{U}` → `cal(U)`; `\boldsymbol` → `bold`; `\overline` → `overline`;
  `\bar x` → `macron(x)`; `\widehat x` → `hat(x)`; `\widetilde x` → `tilde(x)`;
  `\underbrace{a}_{b}` → `underbrace(a, b)`; `\xrightarrow{t}` → `stretch(->)^(t)`.
- `\leq`, `\geq` → `<=`, `>=` (rendered slanted, as in the LaTeX version).
- `\textcolor{red}{x}` → `#text(fill: red, $x$)`.

## Figures, tables, algorithms

```typst
#figure(placement: top, image("/figures/f.pdf", width: 80%), caption: [...]) <fig:x>
#figure(placement: none, kind: table, tabular(columns: 3, [a], [b], [c], ...), caption: [...]) <table:x>
#algorithm(caption: [...], (0, [$i <- 1$]), (0, [*while* $x != 0$ *do*]), (1, [...], "line:label"), ...) <algo:x>
```

- Placement: use `top` for figures without placement or that LaTeX floats to the
  top of a page, `none` for those that LaTeX keeps in place (check the LaTeX PDF).
- `tabular` draws `|c|c|` tables with `\hline` between all rows; a table without
  frame: `table(stroke: none, ...)`.
- Algorithm lines: `(indentation level, content)` or with a label
  `(level, content, "line:x")`, referenced with `#lineref(<line:x>)`.
  Keywords in bold: `*for* ... *do*`, `*end for*`, `*return*`, `*if* ... *then*`.
- TikZ: redraw with `@preview/cetz:0.4.2`, `@preview/fletcher:0.5.8` (commutative
  diagrams, `tikzcd`) or `@preview/cetz-plot:0.1.3` (pgfplots).

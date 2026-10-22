# Conversion progress

## Done
All chapters (frontmatter, introduction, chapters 1-8, appendices A-C,
bibliography) and the ENPC title page and syllabus are converted and included
in `main.typ`; `typst compile --input enpc=true` builds the ENPC version.
Chapter 8 (optimization), commented out in `main.tex`, is included; it has no
LaTeX build to compare against. As in `main.tex`, the ENPC syllabus is
converted but not included. The
book compiles without errors or warnings. For every chapter, the numbering of
equations, theorem-like environments, exercises, figures, tables, algorithms
and tasks matches the LaTeX PDF exactly. Chapter 1 was compared page by page
with the LaTeX PDF; the other chapters were checked by their conversion agents
and by spot checks. The Typst book has 242 pages with chapter 8 (LaTeX, without it: 227).

## Remaining
- Known differences:
  - inline matrices are smaller than in LaTeX;
  - `\ge` (upright in LaTeX) is rendered slanted by the global `≥ → ⩾` rule;
  - `stretch(->)` does not lengthen under long text;
  - code blocks are flush left, whereas minted keeps the source indentation;
  - references use the name of each environment ("Proposition 1.1"), whereas
    cleveref prints "Theorem 1.1" for all environments sharing its counter.
- A few chapters define local helpers at the top of their file (e.g.
  `octagon-graph`, `algorithm-unnumbered`, `numbered-code`); some are
  duplicated and could move to `template.typ`.

Conventions: see CONVERTING.md.

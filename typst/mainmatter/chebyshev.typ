#import "../template.typ": *
#import "@preview/cetz:0.4.2"
#import "@preview/cetz-plot:0.1.3": plot

#heading(level: 1)[Chebyshev polynomials]#label("cha:Chebyshev polynomials")

The Chebyshev polynomials $(T_n)_(n in nat)$ are given on $[-1, 1]$ by the formula
#numbered($ forall x in [-1, 1], wide
  T_n (x) = cos big((n arccos(x))). $) <eq:chebyshev>
Although this formula makes sense only if $x in [-1, 1]$,
the polynomials are defined for all $x in real$.
Equivalently,
the Chebyshev polynomials can be defined from the equation
#numbered($ forall x in [1, infinity), wide
  T_n (x) = cosh big((n arccosh(x))), $) <eq:chebyshev_cosh>
where $cosh(theta) = 1/2 lr((ee^theta + ee^(-theta)))$
and $arccosh colon [1, infinity) -> [0, infinity)$ is the inverse function of $cosh$.
The first few Chebyshev polynomials are illustrated in~@fig:chebyshev.
It is immediate to show the following properties from~#ref(<eq:chebyshev>):
- The roots of $T_n$ are given by
  $ z_k = cos lr((pi/(2n) + (k pi)/n)), wide k = 0, dots, n-1. $
  These are illustrated in~@fig:chebyshev_roots.

- The polynomial $T_n$ takes the value 1 or -1 when evaluated at
  #numbered($ x_k = cos lr(((k pi)/n)), wide k = 0, dots, n. $) <eq:chebyshev_maxima>
  More precisely, it holds that $T_n (x_k) = (-1)^k$.

// pgfplots cycle list "exotic"
#figure(
  placement: top,
  cetz.canvas({
    import cetz.draw: *
    set-style(legend: (fill: white, stroke: 0.4pt, padding: 0.1, spacing: 0.3,
      item: (preview: (width: 0.6, margin: 0.15), spacing: 0.06)),
      axes: (stroke: 0.4pt, tick: (stroke: 0.4pt)))
    let styles = (
      (paint: rgb(0, 128, 128)),
      (paint: rgb(255, 128, 0)),
      (paint: rgb(67, 132, 161)),
      (paint: rgb(255, 76, 76)),
      (paint: rgb(153, 204, 0)),
      (paint: rgb(255, 0, 0), dash: (3pt, 2pt)),
      (paint: rgb(165, 160, 66), dash: (3pt, 2pt)),
    )
    plot.plot(
      size: (11.1, 4.8),
      axis-style: "scientific",
      x-min: -1.2, x-max: 1.2, x-tick-step: 0.2,
      y-min: -1.2, y-max: 1.2, y-tick-step: 1,
      x-label: none, y-label: none,
      legend: "north-east",
      {
        for n in range(7) {
          plot.add(
            domain: (-1, 1), samples: 500,
            x => calc.cos(n * calc.acos(calc.clamp(x, -1, 1))),
            style: (stroke: (thickness: 2pt, ..styles.at(n))),
            label: $T_#n (x)$,
          )
        }
      },
    )
  }),
  caption: [Illustration of the first few Chebyshev polynomials over the interval $[-1, 1]$.],
) <fig:chebyshev>

#exercise[
  Show that~@eq:chebyshev defines a polynomial of degree~$n$,
  and find its expression in the usual polynomial notation.
] <exercise:chebyshev>
#solution[
  The key idea is to rewrite the cosine function in terms of the complex exponential:
  $ cos(n theta) = 1/2 lr((ee^(ii n theta) + ee^(-ii n theta)))
    = 1/2 Big((big((cos(theta) + ii sin(theta)))^n + big((cos(theta) - ii sin(theta)))^n)). $
  By expanding the powers on the right-hand side,
  we obtain
  $
    big((cos(theta) + ii sin(theta)))^n
    &= sum_(j=0)^n binom(n, j) cos(theta)^(n-j) thin ii^j thin sin(theta)^j \
    big((cos(theta) - ii sin(theta)))^n
    &= sum_(j=0)^n binom(n, j) cos(theta)^(n-j) thin (-ii)^j thin sin(theta)^j.
  $
  The terms corresponding to odd values of~$j$ cancel out in the expression of $cos(n theta)$,
  and so we obtain the following expression for $cos(n theta)$ in terms of $cos(theta)$ and $sin(theta)$:
  $
    cos(n theta)
    &= sum_(j=0)^(floor.l n\/2 floor.r) binom(n, 2j) cos(theta)^(n-2j) ii^(2j) sin(theta)^(2j) \
    &= sum_(j=0)^(floor.l n\/2 floor.r) (-1)^j binom(n, 2j) cos(theta)^(n-2j) big((1 - cos(theta)^2))^j.
  $
  Therefore,
  we conclude that
  #numbered($ T_n (x) = sum_(j=0)^(floor.l n\/2 floor.r) binom(n, 2j) x^(n-2j) big((x^2 - 1))^j. $) <eq:expression_chebyshev>
]

#exercise[
  Show that the same polynomials are obtained from~@eq:chebyshev_cosh.
] <exercise:chebyshev_cosh>
#solution[
  Notice that
  $
    cosh(n xi)
    &= 1/2 lr((ee^(n xi) + ee^(-n xi))) \
    &= 1/2 Big((big((cosh(xi) + sinh(xi)))^n + big((cosh(xi) - sinh(xi)))^n)).
  $
  Using the binomial formula, we obtain
  $
    cosh(n xi)
    &= 1/2 sum_(j=0)^n binom(n, j) lr((cosh(xi)^(n-j) sinh(xi)^j + cosh(xi)^(n-j) (-1)^j sinh(xi)^j)) \
    &= 1/2 sum_(j=0)^n binom(n, j) cosh(xi)^(n-j) lr((sinh(xi)^j + (-1)^j sinh(xi)^j)).
  $
  The contributions of the odd values of~$j$ cancel out,
  and so we obtain
  $
    cosh(n xi)
    &= sum_(j=0)^(floor.l n \/ 2 floor.r) binom(n, 2j) cosh(xi)^(n-2j) sinh(xi)^(2j).
  $
  Since $cosh(xi)^2 - sinh(xi)^2 = 1$,
  we deduce that
  $ cosh(n xi)
    = sum_(j=0)^(floor.l n \/ 2 floor.r) binom(n, j) cosh(xi)^(n-2j) (cosh(xi)^2 - 1)^j, $
  which after the substitution of $xi = arccosh(x)$ leads to~@eq:expression_chebyshev.
]

#exercise(title: [Yet another expression for the Chebyshev polynomials])[
  Show that $T_n (x)$ may be defined from the formula
  #numbered($ T_n (x) =
    1/2 Big((x+sqrt(x^2-1)))^n + 1/2 Big((x-sqrt(x^2-1)))^n wide " for "#sym.space.nobreak |x| >= 1. $) <eq:chebyshev_polynomials>
] <exercise:yet_another_expression_cheby>
#solution[
  We showed in the solution of~@exercise:chebyshev_cosh that
  $ cosh(n xi) = 1/2 Big((big((cosh(xi) + sinh(xi)))^n + big((cosh(xi) - sinh(xi)))^n)). $
  Letting $xi = arccosh(x)$ in this equation and using that $cosh(xi)^2 - sinh(xi)^2 = 1$,
  we obtain
  $ T_n (x) = 1/2 Big((big((x + sqrt(x^2 - 1)))^n + big((x - sqrt(x^2 - 1)))^n)), $
  which is the required formula.
]

#exercise(title: [Recursion relation])[
  Show that the Chebyshev polynomials satisfy the relation
  #numbered($ forall n in {1, 2, dots},
    wide
    T_(n+1) = 2 x T_n - T_(n-1). $) <eq:recursion_chebyshev>
] <exercise:recursion_chebyshev>
#solution[
  It is sufficient to show the identity for $x in [-1, 1]$,
  where the formula~@eq:chebyshev applies.
  Using well-known trigonometric identities,
  we have
  $
    cos big(((n+1) theta)) &= cos(n theta) cos(theta) - sin(n theta) sin(theta) \
    cos big(((n-1) theta)) &= cos(n theta) cos(theta) + sin(n theta) sin(theta).
  $
  Adding both equations and rearranging, we obtain
  $ cos big(((n+1) theta)) = 2 cos(n theta) cos(theta) - cos big(((n-1) theta)). $
  Therefore, using this equation with $theta = arccos(x)$,
  we obtain the statement.
]

#remark[
  The recursion relation in~@exercise:recursion_chebyshev
  can be employed to show by recursion that $T_n (x)$ is indeed a polynomial of degree~$n$.
]

#exercise[
  Since $T_n colon real -> real$ is a polynomial,
  it may be written in the standard form
  $ T_n (x) = alpha^((n))_n x^n + dots + alpha^((n))_1 x + alpha^((n))_0. $
  Prove that~$alpha^((n))_n = 2^((n-1))$ provided that $n >= 1$.
] <exercise:chebyshev_leading_coefficient>
#solution[
  From the definition~@eq:chebyshev,
  the Chebyshev polynomials of degrees 0 and 1 are given by $T_0(x) = 1$ and $T_1(x) = x$.
  The statement then follows by recursion, using~@exercise:recursion_chebyshev.
]

#exercise[
  Let $xi in real backslash (-1, 1)$.
  Show that, among all the polynomials in $poly(n)$ that are bounded from above by 1 in absolute value uniformly over the interval $(-1, 1)$,
  the Chebyshev polynomial~$T_n$ achieves the largest absolute value when evaluated at~$xi$.
] <exercise:auxiliary_chebyshev>
#solution[
  Reasoning by contradiction,
  we assume that there exists $p in poly(n)$ that satisfies
  $ sup_(x in (-1, 1)) |p(x)| <= 1
    wide " and " wide |p(xi)| > |T_n (xi)|. $
  Let $q(x) = p(x) T_n (xi) \/ p(xi)$.
  Then by construction $q(xi) = T_n (xi)$
  and
  $ sup_(x in (-1, 1)) |q(x)| < 1. $
  Consequently, denoting by $x_k$ the points defined in~@eq:chebyshev_maxima,
  we have that
  $ forall k in {0, dots, n}, wide
    (-1)^k (T_n - q)(x_k) > 0. $
  In other words, the polynomial $T_n - q$ takes positive values at ${x_0, x_2, x_4, dots}$ and negative values at~${x_1, x_3, x_5, dots}$.
  Consequently, by the intermediate value theorem,
  $T_n - q$ possesses~$n$ distinct roots in the open interval $(-1, 1)$.
  Since, in addition, $(T_n - q)(xi) = 0$,
  we deduce that~$T_n - q$ has $n+1$ distinct roots,
  which is a contradiction given that $T_n - q$ is a nonzero polynomial of degree at most~$n$.
]

#figure(
  placement: top,
  cetz.canvas(length: 6.5cm, {
    import cetz.draw: *
    let n = 8
    // Axes (pgfplots "axis lines = middle")
    line((-1.1, 0), (1.1, 0), stroke: 0.4pt, mark: (end: "stealth", fill: black, scale: 0.5))
    line((0, -0.1), (0, 1.1), stroke: 0.4pt, mark: (end: "stealth", fill: black, scale: 0.5))
    line((0, -0.015), (0, 0.015), stroke: 0.4pt)
    line((-0.015, 1), (0.015, 1), stroke: 0.4pt)
    content((0, 1), anchor: "east", padding: 0.333em, $1$)
    // Half circle
    arc((1, 0), start: 0deg, stop: 180deg, radius: 1, stroke: (paint: red, thickness: 0.8pt))
    for k in range(n) {
      let theta = calc.pi / (2 * n) + k * calc.pi / n
      let (x, y) = (calc.cos(theta), calc.sin(theta))
      line((x, 0), (x, y), stroke: (paint: gray, thickness: 0.8pt, dash: "dashed"))
      circle((x, y), radius: 2pt, fill: black, stroke: none)
      circle((x, 0), radius: 2pt, fill: black, stroke: none)
      content((x, 0), anchor: "north", padding: 0.333em, $z_#k$)
    }
  }),
  caption: [Roots of the Chebyshev polynomial $T_8$.],
) <fig:chebyshev_roots>

#exercise[
  Assume that $0 < lambda_1 < lambda_2$.
  Prove that for any polynomial $p in poly(n)$ that satisfies $p(0) = 1$,
  it holds that
  $ sup_(lambda in (lambda_1, lambda_2)) abs(p(lambda)) >= 1/(T_n (xi)),
    wide xi := (lambda_2 + lambda_1)/(lambda_2 - lambda_1), $
  with equality for
  #numbered($ p_*(lambda) = frac(T_n lr(((lambda_1 + lambda_2 - 2 lambda)/(lambda_2 - lambda_1))), T_n lr(((lambda_1 + lambda_2)/(lambda_2 - lambda_1)))). $) <eq:optimal_polynomial>
] <exercise:conjugate_gradient_chebyshev>
#solution[
  Assume that $p in poly(n)$ is such that $p(0) = 1$,
  and let $q in poly(n)$ be given by
  $ q(mu) = p lr(((lambda_1 + lambda_2 - (lambda_2 - lambda_1) mu)/2))
    wide <=> wide
    p(lambda) = q lr(((lambda_1 + lambda_2 - 2 lambda)/(lambda_2 - lambda_1))). $
  Since $xi > 1$, it holds from~@eq:chebyshev_polynomials that $T_n (xi) > 0$ and
  it follows from~@exercise:auxiliary_chebyshev that
  $ p(0) = q(xi)
    <= T_n (xi) sup_(mu in (-1, 1)) |q(mu)|
    = T_n (xi) sup_(lambda in (lambda_1, lambda_2)) |p(lambda)|, $
  with equality when $q prop T_n$,
  i.e.~when
  $ p(lambda) prop T_n lr(((lambda_1 + lambda_2 - 2 lambda)/(lambda_2 - lambda_1))). $
  The expression~@eq:optimal_polynomial then follows from the fact that $p_*(0) = 1$.
]

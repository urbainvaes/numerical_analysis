#import "../template.typ": *
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

// minted with fontsize=\footnotesize, without caption
#let footnotesize-code(body) = text(size: footnotesize, body)

// Plus sign with reduced spacing, for an equation that TeX fits on the line
// by shrinking the spaces around binary operators (Typst does not shrink them)
#let tplus = $#h(0.12em)#math.class("normal", sym.plus)#h(0.12em)$

// Romberg scheme (tikzcd with options cramped, sep=small); `errors` is the
// last row, with the scaling of the error for each column
#let romberg-scheme(errors) = diagram(
  spacing: (12pt, 8pt),
  node-inset: 4pt,
  edge-stroke: 0.45pt,
  {
    let rows = (
      ($J(h)$,),
      ($J(h\/2)$, $J_1(h\/2)$),
      ($J(h\/4)$, $J_1(h\/4)$, $J_2(h\/4)$),
      ($J(h\/8)$, $J_1(h\/8)$, $J_2(h\/8)$, $J_3(h\/8)$),
    )
    for (i, row) in rows.enumerate() {
      for (j, entry) in row.enumerate() {
        node((j, i), entry)
        if j + 1 < row.len() { edge((j, i), (j + 1, i), "->") }
        if i + 1 < rows.len() { edge((j, i), (j + 1, i + 1), "->") }
      }
    }
    for (j, e) in errors.enumerate() { node((j, 4), e) }
  },
)

= Numerical integration <cha:quadrature>
#minitoc()

#section-star[Introduction]
Integrals are ubiquitous in science and mathematics.
In this chapter,
we are concerned with the problem of calculating numerically integrals of the form
#numbered($ I = integral_Omega u(vect(x)) thin dd vect(x), $) <eq:integral>
Perhaps somewhat surprisingly,
the numerical calculation of such integrals when $n >> 1$ is still a very active area of research today.
In this chapter,
however,
we will focus for simplicity on the one-dimensional setting where $Omega = [a, b] subset real$.
We assume throughout this chapter that the function~$u$ is Riemann-integrable.
Then, by definition,
$
  I = lim_(h -> 0)
  sum_(i=0)^(n-1) u(t_i) (z_(i+1) - z_i),
$
where $a = z_0 < cdots < z_n = b$ is a partition of the interval~$[a, b]$ such that the maximum spacing between successive~$x$ values is equal to~$h$,
and with $t_i in [x_i, x_(i+1)]$ for all $i in {0, dots, n-1}$.

All the numerical integration formulas that we present in this chapter are based on a deterministic approximation of the form
#numbered($ hat(I) = sum_(i=0)^n w_i u(x_i), $) <eq:deterministic_integration>
where~$x_0 < dots < x_n$ are the _integration points_ and $w_0, dots, w_n$ are the _integration weights_.
In many cases, integration formulas contain a small parameter that can be refined to improve the accuracy of the approximation.
In methods based on equidistant interpolation nodes,
for example, this parameter encodes the distance between nodes and is typically denoted by~$h$.
We shall often use the notation~$hat(I)_h$ to emphasize the dependence of the approximation on~$h$.
The difference~$E_h = I - hat(I)_h$ is called the _integration error_ or _discretization error_.
The _degree of precision_,
defined hereafter,
is an important measure of the quality of quadrature rule.
#definition[
  The _degree of precision_ of an integration method is the smallest integer number~$d$ such that
  the integration error is zero for all $u in poly(d)$,
  i.e. for all the polynomials of degree less than or equal to~$d$.
]

#indent We observe that,
without loss of generality,
we can assume that the integration interval is equal to~$[-1, 1]$.
Indeed, using the change of variable
#eq-last($
  zeta colon & [-1, 1] -> [a, b]; \
  & y |-> (b+a)/2 + ((b-a))/2 y,
$) <eq:change_of_variable_integration>
we have
#numbered($
  integral_a^b u(x) thin dd x
  = integral_(-1)^1 u big((zeta(y))) thin zeta'(y) thin dd y
  = (b-a)/2 integral_(-1)^1 u compose zeta (y) thin dd y,
$) <eq:change_of_variable_equivalence>
and the right-hand side is the integral of $u compose zeta$ over the interval~$[-1, 1]$.

== The closed Newton--Cotes method <sec:newton_cotes>

Given a set of equidistant points $-1 = x_0 < cdots < x_n = 1$,
a natural method for approximating the integral~@eq:integral of a function~$u colon [-1, 1] -> real$
is to first construct the interpolating polynomial~$hat(u)$ at the nodes,
and then calculate the exact integral of this polynomial.
By construction, this method is exact for polynomials of degree up to~$n$,
and so the degree of precision is equal to _at least_~$n$.
Let $phi_0, dots, phi_n$ denote the Lagrange polynomials associated with the integration nodes.
Then we have
$
  I approx integral_(-1)^1 hat(u)(x) thin dd x
  = integral_(-1)^1 sum_(i=0)^n u(x_i) phi_i (x) thin dd x
  = sum_(i=0)^n u(x_i) underbrace(integral_(-1)^1 phi_i (x) thin dd x, w_i).
$
The weights are independent of the function~$u$,
and so they can be calculated a priori.
The class of integration methods obtained using this approach are known as _Newton--Cotes methods_.
We present a few particular cases:
- $n = 1$, $d = 1$ (trapezoidal rule):
  #numbered($
    integral_(-1)^1 u(x) thin d x
    approx u(-1) + u(1).
  $) <eq:trapezoidal_rule>

- $n = 2$, $d = 3$ (Simpson's rule):
  #numbered($
    integral_(-1)^1 u(x) thin d x
    approx 1/3 u(-1) + 4/3 u(0) + 1/3 u(1).
  $) <eq:simpsons>

- $n = 3$, $d = 3$ (Simpson's $3/8$ rule):
  $
    integral_(-1)^1 u(x) thin d x
    approx 1/4 u(-1) + 3/4 u(-1\/3) + 3/4 u(1\/3) + 1/4 u(1).
  $

- $n = 4$, $d = 5$ (Boole's rule):
  $
    integral_(-1)^1 u(x) thin d x
    approx 7/45 u(-1) + 32/45 u lr((-1/2)) + 12/45 u lr((0)) + 32/45 u lr((1/2)) + 7/45 u(1).
  $

#remark[
  Note that,
  although it is based on a quadratic polynomial interpolation,
  Simpson's rule~@eq:simpsons has a degree of precision equal to~$3$.
  This is because any integration rule with nodes and weights symmetric around $x=0$ is exact for odd functions,
  in particular $x^3$.
  Likewise, the degree of precision of Boole's rule is equal to 5.
]

#indent In principle,
this approach could be employed in order to construct integration rules of arbitrary high degree of precision.
In practice, however, the weights become more and more imbalanced as the number of interpolation points increases,
with some of them becoming negative.
As a result, roundoff errors become increasingly detrimental to accuracy.
In addition, in cases where the interpolating polynomial does not converge to~$u$,
for example if~$u$ is Runge's function,
the approximate integral may not converge to the correct value in the limit as~$n -> infinity$,
even in exact arithmetic!

The integration rules presented in this section,
which are based on equidistant nodes that include the endpoints of the integration interval,
are called _closed Newton--Cotes_ methods.
A similar approach can be employed in order to construct to integration rules based on equidistant nodes that do not include the endpoints;
these are called _open Newton--Cotes_ methods.

== Composite methods with equidistant nodes <sec:composite_methods>
A natural alternative to the approach presented in~@sec:newton_cotes
is to construct an integration rule using piecewise polynomial interpolation,
which we studied in~@sub:piecewise_interpolation.
After partitioning the integration interval in a number of subintervals,
the integral can be approximated by using one of the rules presented in~@sec:newton_cotes within each subinterval.

#paragraph[Composite trapezoidal rule.][
  Let us illustrate the composite approach with an example.
  To this end,
  we introduce a partition $a = x_0 < cdots < x_n = b$ of the interval $[a, b]$ and
  assume that the nodes are equidistant with $x_(i+1) - x_i = h$.
  Using~@eq:change_of_variable_equivalence with $a = x_i$ and $b = x_(i+1)$,
  we first generalize~@eq:trapezoidal_rule to an interval~$[x_i, x_(i+1)]$ as follows:
  $
    integral_(x_i)^(x_(i+1)) u(x) thin dd x
    = h/2 integral_(-1)^1 u compose zeta (y) thin dd y
    approx h/2 big((u compose zeta(-1) + u compose zeta(1)))
    = h/2 big((u(x_i) + u(x_(i+1)))),
  $
  where $approx$ in this equation indicates approximation using the trapezoidal rule.
  Applying this approximation to each subinterval of the partition,
  we obtain the composite trapezoidal rule:
  #eq-last($
    integral_a^b u(x) thin dd x
    &= sum_(i=0)^(n-1) integral_(x_i)^(x_(i+1)) u(x) thin dd x
    approx
    h/2 sum_(i=0)^(n-1) big((u(x_i) + u(x_(i+1)))) \
    &= h/2 big((u(x_0) + 2 u(x_1) + 2 u(x_2) + cdots + 2 u(x_(n-2)) + 2 u(x_(n-1)) + u(x_n))).
  $) <eq:composite_trapezoidal_rule>
  Like the trapezoidal rule~@eq:trapezoidal_rule,
  the composite trapezoidal rule~@eq:composite_trapezoidal_rule has a degree of precision equal to 1.
  However,
  the integration error of the method depends on the parameter~$h$,
  which represents the width of each subinterval:
  for very small~$h$,
  equation~@eq:composite_trapezoidal_rule is expected to provide a good approximation of the integral.
  An error estimate can be obtained directly from the formula in~@theorem:interpolation_error for the interpolation error,
  provided that we assume that~$u in C^2[a, b]$.

  #theorem(title: [Integration error for the composite trapezoidal rule])[
    Let $hat(I)_h$ denote the approximate integral calculated using~@eq:composite_trapezoidal_rule.
    Then
    #numbered($
      big(|I - hat(I)_h|)
      <= (b-a)/12 C_2 h^2,
      wide
      C_2 := sup_(xi in [a, b]) lr(|u''(xi)|).
    $) <eq:error_bound_trapezoidal>
  ] <theorem:integration_error_trapeze>

  #proof[
    Denoting by~$hat(u)_h$ the piecewise linear interpolation of~$u$,
    we have
    $
      integral_(x_i)^(x_(i+1)) u(x) - hat(u)(x) thin dd x
      = 1/2 integral_(x_i)^(x_(i+1)) u'' big((xi(x))) (x - x_i) (x - x_(i+1)) thin dd x.
    $
    Since $(x - x_i) (x - x_(i+1))$ is nonpositive over the interval $[x_i, x_(i+1)]$,
    we deduce that
    $
      lr(|integral_(x_i)^(x_(i+1)) u(x) - hat(u)(x) thin dd x|)
      <= 1/2 lr((sup_(xi in [a, b]) lr(|u''(xi)|))) integral_(x_i)^(x_(i+1)) (x - x_i) (x_(i+1) - x) thin dd x
      = C_2 h^3/12.
    $
    Summing the contributions of all the intervals,
    we obtain
    $
      big(|I - hat(I)_h|)
      <= sum_(i=0)^(n-1) lr(|integral_(x_i)^(x_(i+1)) u(x) - hat(u)(x) thin dd x|)
      <= n times C_2 h^3/12 = (b-a)/12 C_2 h^2,
    $
    which concludes the proof.
  ]

  #indent The integration error therefore scales as~$bigo(h^2)$.
  (Strictly speaking, we have shown only that the integration error admits an upper bound that scales at $bigo(h^2)$,
  but it turns out that the dependence on~$h$ of this bound is optimal).

  #example(title: [Integration error for the composite trapezoidal rule])[
    Let us approximate the integral of $x^2$ over $[0, 1]$
    via the trapezoidal rule.
    In other words, let us calculate the right-hand side of
    $
      integral_0^1 x^2 thin dd x
      approx
      h/2 lr((x_0^2 + 2 x_1^2 + dots + 2 x_(n - 1)^2 + x_n^2)),
    $
    which we shall denote $hat(I)_h$ hereafter.
    Using $x_i = i h$, and $x_n = 1$, we obtain
    $
      hat(I)_h = h/2 lr((2 h^2 sum_(i = 1)^(n - 1) i^2 + 1))
      = h/2 lr((h^2/3 (n - 1) n (2n - 1) + 1))
      = h/2 lr((h^2/3 (2 n^3 - 3 n^2 + n) + 1)).
    $
    Now, observing that $n = 1 \/ h$ yields
    $
      hat(I)_h
      = h/2 lr((h^2/3 (2 h^(-3) - 3 h^(-2) + h^(-1)) + 1))
      = h/2 lr((2/3 h^(-1) - 1 + h/3 + 1))
      = 1/3 + h^2/6,
    $
    which is exactly the error shown by~@theorem:integration_error_trapeze.
  ]
]

#paragraph[Composite Simpson rule.][
  The composite Simpson rule is derived in~@exercise:composite_simpson.
  Given an odd number $n+1$ of equidistant points $a = x_0 < x_1 < cdots < x_n = b$,
  this rule is given by
  // The LaTeX source has the typo u_(x_4) instead of u(x_4); tplus and the
  // trailing space: see above, to avoid overlapping the equation number
  #numbered($
    hat(I)_h = h/3 Big((u(x_0) tplus 4 u(x_1) tplus 2 u(x_2) tplus 4 u(x_3) tplus 2 u(x_4) tplus cdots tplus 2 u(x_(n-2)) tplus 4 u(x_(n-1)) tplus u(x_n))). #h(2.2em)
  $) <eq:composite_simpson>
  This approximation is obtained by integrating the piecewise quadratic interpolant
  over a partition of the integration interval into $n\/2$ subintervals of equal width.
  Obtaining an optimal error estimate,
  in terms of the dependence on~$h$,
  for this integration formula is slightly more involved.

  #theorem(title: [Integration error for the composite Simpson rule])[
    Let $hat(I)_h$ denote the approximate integral calculated using~@eq:composite_simpson.
    Then
    #numbered($
      big(|I - hat(I)_h|)
      <= (b-a) (C_4 h^4)/180,
      wide
      C_4 := sup_(xi in [a, b]) lr(|u^((4))(xi)|).
    $) <eq:error_bound_simpson>
  ] <theorem:integration_error_simpson>
  #proof[
    For a given subinterval~$[x_(2i), x_(2i+2)]$,
    let us denote by~$hat(u)_2$ the quadratic interpolating polynomial at $x_(2i), x_(2i+1), x_(2i+2)$,
    and by $hat(u)_3(placeholder; alpha)$ the cubic interpolating polynomial relative to the nodes $x_(2i), x_(2i+1), x_(2i+2), alpha$,
    for some $alpha in [x_(2i), x_(2i+1)]$ that does not coincide with the integration nodes.
    We have
    #numbered($
      integral_(x_(2i))^(x_(2i+2)) u(x) - hat(u)_2(x) thin dd x
      = integral_(x_(2i))^(x_(2i+2)) u(x) - hat(u)_3(x; alpha) thin dd x + integral_(x_(2i))^(x_(2i+2)) hat(u)_3(x; alpha) - hat(u)_2(x) thin dd x.
    $) <eq:composite_simpson_error>
    The second term on the right-hand side is zero,
    because the integrand is a cubic polynomial with zeros at $x_(2i)$, $x_(2i+1)$ and $x_(2i+2)$,
    and because
    $
      integral_(x_(2i))^(x_(2i+2)) (x - x_(2i)) (x - x_(2i+1)) (x - x_(2i+2)) = 0.
    $
    // (The integral is zero because the integrand is odd around the vertical axis $x = x_(2i+1)$.)
    By~@theorem:interpolation_error,
    the first term in~@eq:composite_simpson_error is bounded from above as follows:
    $
      lr(|integral_(x_(2i))^(x_(2i+2)) u(x) - hat(u)_3(x; alpha) thin dd x|)
      &<= integral_(x_(2i))^(x_(2i+2)) lr(|(u^((4)) big((xi(x))))/24 (x-x_(2i))(x- x_(2i+1)) (x-x_(2i+2)) (x - alpha)|) thin dd x \
      &<= C_4/24 integral_(x_(2i))^(x_(2i+2)) big(|(x-x_(2i))(x- x_(2i+1)) (x-x_(2i+2)) (x - alpha)|) thin dd x.
    $
    This inequality is valid for any $alpha in A := [x_(2i), x_(2i+2)] backslash {x_(2i), x_(2i+1), x_(2i+2)}$.
    Denoting by $h(alpha)$ the integral on the right-hand side,
    we observe that
    $
      lim_(alpha -> x_(2i+1)) h(alpha) =
      integral_(x_(2i))^(x_(2i+2)) (x-x_(2i))(x- x_(2i+1))^2 (x_(2i+2) - x) thin dd x = 4/15 h^5.
    $
    Therefore, we conclude that
    $
      lr(|integral_(x_(2i))^(x_(2i+2)) u(x) - hat(u)_2(x) thin dd x|)
      <= inf_(alpha in A) lr(|integral_(x_(2i))^(x_(2i+2)) u(x) - hat(u)_3(x; alpha) thin dd x|) <= C_4/90 h^5.
    $
    Summing the contributions of all the subintervals,
    we finally obtain
    #numbered($
      abs(I - hat(I)_h) <= n/2 times (C_4 h^5)/90 = (b-a) (C_4 h^4)/180,
    $) <eq:scaling_error>
    which concludes the proof.
  ]

  #remark[
    The cancellation of the second term in~@eq:composite_simpson_error also follows from the fact that
    the degree of precision of the Simpson rule~@eq:simpsons is equal to~3,
    and so
    $
      integral_(x_(2i))^(x_(2i+2)) hat(u)_3(x) - hat(u)_2(x) thin dd x = 1/3 (hat(u)_3 - hat(u)_2)(x_(2i)) + 4/3 (hat(u)_3 - hat(u)_2)(x_(2i+1)) + 1/3 (hat(u)_3 - hat(u)_2)(x_(2i+2)) = 0,
    $
    where we used the short-hand notation $hat(u)_3(x) = hat(u)_3(x; alpha)$.
  ]
]

#paragraph[General composite quadrature rules.][
  In view of~@eq:change_of_variable_equivalence,
  any quadrature rule with $N + 1$ integration points of the form
  #numbered($
    integral_(-1)^1 u(x) thin dd x
    approx sum_(j=0)^N w_j u(x_j)
  $) <eq:one-interval_rule>
  admits a composite version for an arbitrary integration interval $[a, b]$ obtained by applying the rule locally in $M$ equally-sized subintervals:
  #numbered($
    integral_a^b u(x) thin dd x
    approx h/2 sum_(i=0)^(M - 1) sum_(j=0)^N w_j u lr((a + i h + h/2 + (x_j h)/2)),
    wide h := (b - a)/M.
  $) <eq:composite_rule>
  Clearly, the single-interval rule~@eq:one-interval_rule and the composite rule~@eq:composite_rule have the same degree of precision.
  To conclude this section,
  we prove a relation between this degree of precision of a rule and
  the convergence rate with respect to~$h$ of the composite rule~@eq:composite_rule.

  #theorem[
    Suppose that the degree of precision of the single-interval rule~@eq:one-interval_rule is equal to~$d$.
    Then for any $u in C^((d+1))[a, b]$,
    there is $c > 0$ such that
    #numbered($
      forall M >= 1, wide
      Big(|I[u] - hat(I)_h [u]|)
      <= c h^(d+1),
      wide c := (C_(d+1) (b-a))/((d+1)!) lr((1 + 1/2 sum_(i=0)^J |w_j|)).
    $) <eq:theorem_order_from_precision>
    Here $I[u]$ and $hat(I)_h [u]$ denote respectively the left-hand side and the right-hand side of~@eq:composite_rule,
    and $C_(d+1)$ is defined in~@eq:local_error.
  ] <theorem:relation_precision_order>
  #proof[
    Fix $i in {0, dots M - 1}$,
    and let $p_i$ denote the Taylor expansion of degree~$d$ of the function~$u$ around $z_i := a + i h$:
    #numbered($
      p_i (x) = u(z_i) + u'(z_i) (x - z_i) + dots + (u^((d))(z_i))/(d!) (x - z_i)^d.
    $) <eq:taylor_for_integration>
    By the mean-value form of the remainder,
    it is simple to prove that
    #numbered($
      forall x in [z_i, z_(i+1)], wide
      big(|u(x) - p_i (x)|) = C_(d+1)/((d+1)!) h^(d+1), wide C_(d+1) = sup_(x in [a, b]) Big(|u^((d+1))(x)|).
    $) <eq:local_error>
    Let us define the local contributions to the total integral as follows:
    $
      I^((i))[u] := integral_(z_i)^(z_(i+1)) u(x) thin dd x,
      wide
      hat(I)_h^((i))[u] = h/2 sum_(i=0)^N w_i u lr((z_i + h/2 + (x_i h)/2)),
    $
    so that $I[u] = sum_(i=1)^(M - 1) I^((i))[u]$ and $hat(I)_h [u] = sum_(i=1)^(M - 1) hat(I)_h^((i))[u]$.
    Since $p$ is a polynomial of degree at most~$d$,
    it follows that
    $
      lr(|I^((i))[u] - hat(I)^((i))_h [u]|)
      &= Big(|I^((i))[u] - I^((i))[p] + hat(I)_h^((i))[p] - hat(I)_h^((i))[u]|) \
      &<=
      Big(|I^((i))[u] - I^((i))[p]|)
      + Big(|hat(I)_h^((i))[u] - hat(I)_h^((i))[p]|)
      =
      Big(|I^((i))[u - p]|)
      + Big(|hat(I)_h^((i))[u - p]|),
    $
    where we used the triangle inequality and linearity in the second line.
    Thus, using~@eq:local_error we deduce that
    $
      lr(|I^((i))[u] - hat(I)^((i))_h [u]|)
      &<= C_(d+1)/((d+1)!) h^(d+1) integral_(z_i)^(z_(i+1)) 1 thin dd x + h/2 C_(d+1)/((d+1)!) h^(d+1) sum_(i=0)^J |w_j| \
      &= C_(d+1)/((d+1)!) h^(d+2) lr((1 + 1/2 sum_(i=0)^J |w_j|)).
    $
    Summing the contributions of the local errors,
    we obtain
    $
      lr(|I[u] - hat(I)_h [u]|)
      <= sum_(i=0)^(M-1) lr(|I^((i))[u] - hat(I)_h^((i))[u]|)
      <= M C_(d+1)/((d+1)!) h^(d+2) lr((1 + 1/2 sum_(i=0)^J |w_j|)),
    $
    which leads to the result since $M h = b - a$.
  ]

  #remark[
    A few comments are in order.
    - The constant $c$ in~@eq:theorem_order_from_precision is not sharp in general,
      but the power of~$h$ is optimal.

    - This result may be viewed as a generalization of~@theorem:integration_error_trapeze[Theorems] and~#ref(<theorem:integration_error_simpson>, supplement: []),
      albeit with a worse constant prefactor.

    - Instead of the polynomial~$p_i$ in~@eq:taylor_for_integration,
      we could have used any polynomial approximation of~$u$ such that~@eq:local_error is satisfied,
      with possibly a different prefactor but the same power of~$h$ on the right-hand side.
      A natural choice, for example,
      would have been to define~$p_i$ by interpolation through (possibly a subset or superset of) the local integration points.

    - @theorem:relation_precision_order further motivates why the degree of precision of an integration rule is an interesting metric.
  ]
] <par:General_composite_quadrature_rules>

#paragraph[Estimating the error a posteriori.][
  In practice,
  it is useful to be able to estimate the integration error so that,
  if the error is deemed too large,
  a better approximation of the integral can be calculated by using a smaller value for the step size~$h$.
  Calculating the exact error $I - hat(I)_h$ is impossible in general,
  because this would require to know the exact value of the integral,
  but it is possible to calculate a rough approximation of the error based on two numerical approximations of the integral,
  as we illustrate formally hereafter for the composite Simpson rule.

  Suppose that $hat(I)_(2h)$ and $hat(I)_h$ are two approximations of the integral,
  calculated using the composite Simpson rule with step size $2h$ and $h$, respectively.
  If we assume that the error proportionally to~$bigo(h^4)$ as~@eq:scaling_error suggests,
  then it holds approximately that
  #numbered($
    I - hat(I)_h approx 1/2^4 (I - hat(I)_(2h)).
  $) <eq:relation_errors>
  This implies that
  $
    I - hat(I)_(2h) = (I - hat(I)_h) + (I_h - hat(I)_(2h)) approx 1/16 (I - hat(I)_(2h)) + (hat(I)_h - hat(I)_(2h)).
  $
  Rearranging this equation gives an approximation of the error for $hat(I)_(2h)$:
  $
    I - hat(I)_(2h) approx 16/15 (hat(I)_h - hat(I)_(2h)).
  $
  Using~@eq:relation_errors,
  we can then derive an error estimate for $hat(I)_h$:
  #numbered($
    abs(I - hat(I)_h) approx 1/15 abs(hat(I)_h - hat(I)_(2h)).
  $) <eq:integration_error_approximation>
  The right-hand side can be calculated numerically,
  because it does not depend on the exact value of the integral.
  In practice,
  the two sides of~@eq:integration_error_approximation are often very close for small~$h$.
  In the code example below,
  we approximate the integral
  #numbered($
    I = integral_0^(pi/2) cos(x) thin dd x = 1
  $) <eq:integral_example>
  for different step sizes
  and compare the exact error with the approximate error obtained using~@eq:integration_error_approximation.
  The results obtained are summarized in~@tab:discretization_error,
  which shows a good match between the two quantities.
  #[
    // Caption above the table, as in the LaTeX source
    #set figure.caption(position: top)
    #figure(
      placement: none,
      kind: table,
      // \arraystretch = 1.5
      tabular(
        columns: 3, stretch: 1.5,
        [$h$], [Exact error $abs(I - hat(I)_h)$], [Approximate error $1/15 abs(hat(I)_h - hat(I)_(2h))$],
        [$2^(-4)$], [$5.166847063531321 times 10^(-7)$], [$5.185892840930961 times 10^(-7)$],
        [$2^(-5)$], [$3.226500089326123 times 10^(-8)$], [$3.229464703065806 times 10^(-8)$],
        [$2^(-6)$], [$2.0161285974040766 times 10^(-9)$], [$2.016591486390477 times 10^(-9)$],
        [$2^(-7)$], [$1.2600120946615334 times 10^(-10)$], [$1.260084925291949 times 10^(-10)$],
      ),
      caption: [Comparison between the exact integration error and the approximate integration error calculated using~@eq:integration_error_approximation.],
    ) <tab:discretization_error>
  ]

  #footnotesize-code[
```julia
# Composite Simpson's rule
function composite_simpson(u, a, b, n)
    # Integration nodes
    x = LinRange(a, b, n + 1)
    # Evaluation of u at the nodes
    ux = u.(x)
    # Step size
    h = x[2] - x[1]
    # Approximation of the integral
    return (h/3) * sum([ux[1]; ux[end]; 4ux[2:2:end-1]; 2ux[3:2:end-2]])
end
# Function to integrate
u(x) = cos(x)
# Integration bounds
a, b = 0, π/2
# Exact integral
I = 1.0
# Number of subintervals
ns = [8; 16; 32; 64; 128]
# Approximate integrals
Î = composite_simpson.(u, a, b, ns)
# Calculate exact and approximate errors
for i in 2:length(ns)
    println("Exact error: $(I - Î[i]), ",
            "Approx error: $((Î[i] - Î[i-1])/15)")
end
```
  ]
]

== Richardson extrapolation and Romberg's method
In the previous section,
we showed how the integration error could be approximated based on two approximations of the integral with different step sizes.
The aim of this section is to show that,
by cleverly combining two approximations $hat(I)_h$ and $hat(I)_(2h)$ of an integral,
an approximation even better than $hat(I)_h$ can be constructed.

This approach is based on _Richardson's extrapolation_,
which is a general method for accelerating the convergence of sequences,
with applications beyond numerical integration.
The idea is the following:
assume that $J(h)$ is an approximation with step size~$h$ of some unknown quantity~$J_* = lim_(h -> 0) J(h)$,
and that we have access to evaluations of $J$ at $h, h\/2, h\/4, h\/8 dots$.
If $J$ extends to a smooth function over $[0, H]$,
then by Taylor expansion it holds that
$
  J(eta) = J(0) + J'(0) eta + J''(0) eta^2/2 + J^((3))(0) eta^3/3! + cdots + J^((k))(0) eta^k/k! + bigo(eta^(k+1)).
$

#paragraph[Elimination of the linear error term.][
  Let us assume that $J'(0) != 0$,
  so that the leading order term after the constant $J(0)$ scales as~$eta$.
  Then we have
  $
    J(h) &= J(0) + J'(0) h + bigo(h^2) \
    J(h\/2) &= J(0) + J'(0) h/2 + bigo(h^2).
  $
  We now ask the following question:
  can we combine linearly $J(h)$ and $J(h\/2)$ in order to construct an approximation $J_1(h\/2)$ of~$J(0)$ with an error scaling as $bigo(h^2)$?
  Employing the ansatz
  $J_1(h\/2) = alpha J(h) + beta J(h\/2)$,
  we calculate
  #numbered($
    J_1(h\/2) = (alpha + beta) J(0) + J'(0) h lr((alpha + (1 - alpha)/2)) + bigo(h^2).
  $) <eq:expansion_richardson>
  Since we want this expression to approximate $J(0)$ for small~$h$,
  we need to impose that $alpha + beta = 1$.
  Then, in order for the term multiplying~$h$ to cancel out,
  we require that
  $
    alpha + (1 - alpha)/2 = 0
    wide <=> wide
    alpha = - 1.
  $
  This yields the formula
  #numbered($
    J_1(h\/2) = 2J(h\/2) - J(h).
  $) <eq:richardson>
  Notice that, in the case where $J$ is a linear function,
  $J_1(h\/2)$ is exactly equal to~$J(0)$.
  This reveals a geometric interpretation of~@eq:richardson:
  the approximation $J_1(h\/2)$ is simply the $y$ intercept of the straight line passing through the points $big((h\/2, J(h\/2)))$ and $big((h, J(h)))$.
]

#paragraph[Elimination of the quadratic error term.][
  If we had tracked the coefficient of $h^2$ in the previous paragraph,
  we would have obtained instead of~@eq:expansion_richardson the following equation:
  $
    J_1(h\/2) = J(0) - J^((2))(0) h^2/4 + bigo(h^3).
  $
  Provided that we have access also to $J(h\/4)$,
  we can also calculate
  $
    J_1(h\/4) = 2J(h\/4) - J(h\/2) = J(0) - J^((3))(0) h^2/16 + bigo(h^3).
  $
  At this point,
  it is natural to wonder whether we can combine $J_1(h\/2)$ and $J_1(h\/4)$ in order to produce an even better approximation of $J(0)$.
  Applying the same reasoning as in the previous section leads us to introduce
  $
    J_2(h\/4) = (4 J_1(h\/4) - J_1(h\/2))/(4 - 1) = J(0) + bigo(h^3).
  $
  This is an exact approximation of $J(0)$ if $J$ is a quadratic polynomial,
  indicating that $J_2(h\/4)$ is simply the $y$ intercept of the quadratic polynomial interpolating the function~$J$ through the three points
  $big((h\/4, J(h\/4)))$, $big((h\/2, J(h\/2)))$ and $big((h, J(h)))$.
]

#paragraph[Elimination of higher order terms.][
  The procedure above can be repeated in order to eliminate terms of higher and higher orders.
  The following schematic illustrates,
  for example, the calculation of an approximation $J_3(h\/8) = J(0) + bigo(h^4)$.
  #block(above: 16.3pt, below: 16.3pt, width: 100%, align(center,
    romberg-scheme(($bigo(h)$, $bigo(h^2)$, $bigo(h^3)$, $bigo(h^4).$))))
  Here, the last row indicates the scaling of the error with respect to the parameter $h$ in the limit as~$h -> 0$.
  The linear combination in order to calculate $J_i (h\/2^i)$ is always of the form
  $
    J_i (h\/2^i) = (2^i J_(i-1) (h\/2^i) - J_(i-1)(h\/2^(i-1)))/(2^i - 1), wide J_0 = J.
  $
  In practice we calculate the values taken by $J, J_1, J_2, dots$ at specific values of~$h$,
  but these are in fact functions of~$h$.
  In~@fig:richardson,
  we plot these functions when~$J(h) = 1 + sin(h)$.
  It appears clearly from the figure that,
  for sufficiently small~$h$,
  $J_3(h)$ provides the most precise approximation of $J(0) = 1$.
  Constructing the functions in Julia can be achieved in just a few lines of code.
  #footnotesize-code[
```julia
J(h) = 1 + sin(h)
J_1(h) = 2J(h) - J(2h)
J_2(h) = (4J_1(h) - J_1(2h))/3
J_3(h) = (8J_2(h) - J_2(2h))/7
```
  ]
  #figure(
    placement: top,
    image("/figures/richardson.pdf", width: 80%),
    caption: [Illustration of the functions $J_1$, $J_2$ and $J_3$ constructed by Richardson extrapolation.],
  ) <fig:richardson>
]

#paragraph[Generalization.][
  Sometimes,
  it is known a priori that the Taylor development of the function~$J$ around zero contains only even powers of~$h$.
  In this case, the Richardson extrapolation procedure can be slightly modified in order to produce approximations with errors scaling as~$bigo(h^4)$,
  then~$bigo(h^6)$, then~$bigo(h^8)$, etc.
  This procedure is illustrated below:
  #block(above: 16.3pt, below: 16.3pt, width: 100%, align(center,
    romberg-scheme(($bigo(h^2)$, $bigo(h^4)$, $bigo(h^6)$, $bigo(h^8).$))))
  This time, the linear combinations required for populating this table are given by
  #numbered($
    J_i (h\/2^i) = (2^(2i) J_(i-1) (h\/2^i) - J_(i-1)(h\/2^(i-1)))/(2^(2i) - 1).
  $) <eq:generalized_richardson>
]

#paragraph[Application to integration: Romberg's method][
  Romberg's integration method consists of applying Richardson's extrapolation to the function
  $
    J(h) = hat(I)_h = u(x_0) + 2 u(x_1) + 2 u(x_2) + cdots + 2 u(x_(n-1)) + 2 u(x_n), wide h in lr({(b-a)/n: n in nat}).
  $
  where $a = x_0 < x_1 < cdots < x_n = b$ are equidistant nodes.
  The right-hand side of this equation is simply the composite trapezoidal rule with step size~$h$.
  It is possible to show that~$J(h)$ may be expanded as follows:
  #numbered($
    forall k in nat, wide
    J(h) = I + alpha_1 h^2 + alpha_2 h^4 + cdots + alpha_k h^(2k) + bigo(h^(2k+2)).
  $) <eq:expansion_error_composite_trapezium>
  This is the content of the following result.
  #lemma[
    Let $J(h)$ denote the approximation of the integral $I[u]$ by the composite trapezium rule,
    and assume that $u in C^infinity [a, b]$.
    Then $J(h)$ may be expanded as in @eq:expansion_error_composite_trapezium.
  ]
  #proof[
    Using the same notation as in~@theorem:relation_precision_order,
    we introduce
    $
      I^((i))[u] = integral_(x_i)^(x_(i+1)) u(x) thin dd x,
      wide
      hat(I)_h^((i))[u] = h/2 Big((u(x_i) + u(x_(i+1)))),
      wide i = 0, dots, n-1.
    $
    Fix $i in {0, dots, n-1}$ and let $p_i$ denote the Taylor expansion of~$u$ around $x_(i + 1/2) := x_i + h/2$,
    of degree $2k + 1$:
    $
      p_i (x) = u big((x_(i + 1/2))) + u' big((x_(i + 1/2))) Big((x - x_(i + 1/2))) + cdots + (u^((2k+1)) big((x_(i + 1/2))))/((2k+1)!) Big((x - x_(i + 1/2)))^(2k+1).
    $
    From the mean-value form of the remainder,
    it is clear that
    $
      I^((i))[u] = I^((i))[p_i] + bigo(h^(2k+3)),
      wide
      hat(I)_h^((i))[u] = hat(I)_h^((i))[p_i] + bigo(h^(2k+3)).
    $
    Substituting~$p_i$ and noting that odd powers cancel when integrating,
    we have that
    $
      I^((i))[p_i] &= u big((x_(i + 1/2))) h + omega_2 h^3 + omega_4 h^5 + cdots + omega_(2k) h^(2k+1), \
      hat(I)_h^((i))[p_i] &= u big((x_(i + 1/2))) h + eta_2 h^3 + eta_4 h^5 + cdots + eta_(2k) h^(2k+1),
    $
    for appropriate coefficients.
    Thus we obtain that
    $
      I[u] - hat(I)_h [u]
      &= sum_(i=0)^n I^((i))[p_i] - I_h^((i))[p_i] \
      &= (omega_2 - eta_2) h^2 + (omega_2 - eta_2) h^4 + dots + (omega_(2k) - eta_(2k)) h^(2k) + bigo(h^(2k+2)),
    $
    which concludes the proof.
  ]

  #indent Richardson's extrapolation~@eq:generalized_richardson can therefore be employed in order to compute approximations of the integral with increasing accuracy.
  The convergence of Romberg's method for calculating the integral~@eq:integral_example is illustrated in~@fig:romberg.
  #figure(
    placement: none,
    image("/figures/romberg.pdf", width: 95%),
    caption: [Convergence of Romberg's method. The straight lines correspond to the monomial functions $f(h) = C_i h^i$,
      with $i = 2, 4, 6, 8$ and for appropriate constants $C_i$.
      We observe a good agreement between the observed and theoretical convergence rates.],
  ) <fig:romberg>
]

== Methods with non-equidistant nodes
The Newton--Cotes method relies on equidistant integration nodes,
and the only degrees of freedom are the integration weights.
If the nodes are not fixed,
then additional degrees of freedom are available,
and these can be leveraged in order to construct a better integration formula.
The total number of degrees of freedom for a general integration rule of the form~@eq:deterministic_integration is $2n + 2$ which,
in principle, should enable to construct an integration rule with a degree of precision equal to~$2n+1$.

A necessary condition for an integration rule of the form~@eq:deterministic_integration to have a degree of precision
equal to $2 n +1$ is that it integrates exactly all the monomials of degree $0$ to~$2 n+1$.
This condition is also sufficient because,
assuming that it is satisfied,
we have by linearity of the functionals~$I$ and $hat(I)$ that
$
  hat(I) lr((alpha_0 + alpha_1 x + cdots + alpha_(2n+1) x^(2n+1)))
  &= alpha_0 hat(I)(1) + alpha_1 hat(I)(x) + cdots + alpha_(2n+1) hat(I) lr((x^(2n+1))) \
  &= alpha_0 I(1) + alpha_1 I(x) + cdots + alpha_(2n+1) I lr((x^(2n+1))) \
  &= I lr((alpha_0 + alpha_1 x + cdots + alpha_(2n+1) x^(2n+1))),
$
Here $I(u)$ and $hat(I)(u)$ denote respectively the exact integral of~$u$ and its approximate integral using~@eq:deterministic_integration.
In order to find the nodes and weights of the integration rule,
we can therefore solve the following nonlinear system of $2n+2$ equations with $2n+2$ unknowns:
#numbered($
  sum_(i=0)^n w_i x_i^d &= integral_(-1)^1 x^d thin dd x, wide d = 0, dots, 2n+1.
$) <eq:system_gauss_legendre>
The quadrature rule thus obtained is called the _Gauss--Legendre quadrature_.
#example[
  Let us derive the Gauss--Legendre quadrature with $n + 1 = 2$ nodes.
  The system of equations that we need to solve in this case is the following:
  $
    w_0 + w_1 = 2,
    wide w_0 x_0 + w_1 x_1 = 0,
    wide w_0 x_0^2 + w_1 x_1^2 = 2/3,
    wide w_0 x_0^3 + w_1 x_1^3 = 0.
  $
  The solution to these equations is given by
  $
    - x_0 = x_1 = sqrt(3)/3, wide w_0 = w_1 = 1.
  $
]

#paragraph[Connection with orthogonal polynomials.][
  Let $(L_n)_(n in nat)$ denote the Legendre polynomials,
  i.e. the orthogonal polynomials with for the inner product
  $
    ip(f, g) = integral_(-1)^1 f(x) g(x) thin dd x.
  $
  The nodes and weights of the Gauss--Legendre quadrature rules can be obtained constructively from Legendre polynomials,
  as shown in~@sub:orthogonal_integration.
  We shall now demonstrate this connection in much more direct manner.
  Specifically, we prove that the integration nodes are given by the roots of a Legendre polynomial.
  #theorem[
    For every $n in nat$,
    there exists a unique solution to the system of equations~@eq:system_gauss_legendre.
    The nodes $(x_i)_(i in {0, dots, n})$ are the roots of $L_(n+1)$
    and the weights are given by
    #numbered($
      w_i = integral_(-1)^1 ell_i (x) thin dd x, wide i in {1, dots, n},
    $) <eq:weights_gauss>
    where $ell_i$ is the Lagrange polynomial
    $
      ((x - x_0) cdots (x - x_(i-1))(x- x_(i+1)) cdots (x - x_n))
      /((x_i - x_0) cdots (x_i - x_(i-1))(x_i- x_(i+1)) cdots (x_i - x_n)).
    $
    In addition, the weights are all positive.
  ] <theorem:gauss_legendre_connection_polynomials>
  #proof[
    We show first that the solution to~@eq:system_gauss_legendre exists,
    then that it is unique,
    and finally that the weights are positive.
    #paragraph[Existence.][
      We begin by showing that,
      if $(x_i)_(i in {0, dots n})$ are the roots of~$L_(n+1)$ and the weights are defined from~@eq:weights_gauss,
      then the equations~@eq:system_gauss_legendre are satisfied.
      To this end,
      it is sufficient to show that for all $p in poly(2n + 1)$,
      #numbered($
        integral_(-1)^1 p(x) thin dd x = sum_(i=0)^n w_i p(x_i).
      $) <eq:gauss_legendre_equation>
      Take $p in poly(2n + 1)$ and let $q in poly(n)$ and $r in poly(n)$ be the polynomials such that
      $
        p(x) = q(x) L_(n+1)(x) + r(x).
      $
      The quotient $q$ and the remainder $r$ can be obtained by Euclidean division of $p$ by $L_(n+1)$.
      Since~$L_(n+1)$ is orthogonal to any polynomial in~$poly(n)$,
      in particular $q$,
      and the nodes $(x_i)_(i in {0, dots n})$ are the roots of~$L_(n+1)$,
      it holds that
      #numbered($
        integral_(-1)^1 p(x) thin dd x
        - sum_(i=0)^n w_i p(x_i)
        = integral_(-1)^1 r(x) thin dd x
        - sum_(i=0)^n w_i r(x_i).
      $) <eq:intermediate_gauss>
      Given that $r in poly(n)$,
      the remainder $r$ must coincides with its polynomial interpolation at the points $x_0, dots x_m$.
      Therefore,
      $
        r = r(x_0) ell_0 + cdots + r(x_n) ell_n,
      $
      and so
      $
        integral_(-1)^1 r(x) thin dd x
        = sum_(i=0)^n r(x_i) integral_(-1)^1 ell_i (x) thin dd x
        = sum_(i=0)^n r(x_i) w_i,
      $
      where we used~@eq:weights_gauss in the last equality.
      Consequently, the right-hand side of~@eq:intermediate_gauss is zero,
      and since $p$ was arbitrary this implies that~@eq:gauss_legendre_equation is satisfied.
    ]

    #paragraph[Uniqueness.][
      Next,
      we show that the nodes are necessarily the roots of $L_(n+1)$.
      To this end,
      assume that $x_0, dots, x_n$ and weights $w_0, dots, w_n$ are such that the equations~@eq:system_gauss_legendre are satisfied,
      and let
      $
        q(x) = (x - x_0) dots (x - x_n).
      $
      Our goal is to show that $q(x)$ coincides with $L_(n+1)$ up to a constant factor.
      In order to show this,
      it is sufficient to prove that $q(x)$ is orthogonal to $x^d$ for all $d = 0, dots, n$,
      because the only polynomial in~$poly(n+1)$ that satisfies these orthogonality relations is the Legendre polynomial~$L_(n+1)$,
      or a multiple thereof.
      For the values of~$d$ considered, the polynomial~$q(x) x^d$ belongs to $poly(2n+1)$.
      Given that the integration rule with nodes $x_0, dots, x_n$ and weights $w_0, dots, w_n$ is exact for all the elements of $poly(2n+1)$ by assumption,
      we deduce that
      $
        forall d in {0, dots, n},
        wide integral_(-1)^1 q(x) x^d thin dd x = sum_(i=0)^n w_i q(x_i) x_i^d = 0.
      $
      Finally, we show that the weights are necessarily given by~@eq:weights_gauss.
      Since the integration rule must be exact for any Lagrange polynomial $ell_j$,
      we have
      $
        integral_(-1)^1 ell_j (x) thin dd x = sum_(i=0)^n w_i ell_j (x_i) = w_j,
      $
      which concludes the proof of uniqueness.
    ]
    #paragraph[Positivity of the weights.][
      Since the integration rule is exact for all the polynomials in~$poly(2n+1)$ and $ell_j (x)^2 in poly(2n+1)$,
      we deduce that
      $
        integral_(-1)^(-1) big(|ell_j (x)|)^2 thin dd x
        = sum_(i=0)^n w_i big(|ell_j (x_i)|)^2 = w_j.
      $
      The left-hand side is positive,
      and so $w_j$ must also be positive.
    ]
  ]

  #indent Since the integration weights are all positive,
  the Gauss--Legendre quadrature rules are less susceptible to roundoff errors than the Newton--Cotes methods.
  In addition, we have the following result.

  #theorem[
    Assume that $u in C([-1, 1])$,
    and let $hat(I)_n (u)$ denote the approximation of~$I(u)$ using the Gauss--Legendre quadrature.
    Then $hat(I)_n (u) -> I(u)$ in the limit as $n -> infinity$.
  ]
  #proof[
    The positivity of the weights is crucial for the proof.
    By the Weierstrass approximation theorem,
    for all $epsilon > 0$ there exists polynomial $p$ such that
    $
      E := max_(x in [-1, 1]) Big(|p(x) - u(x)|) <= epsilon.
    $
    Since the degree of precision of the Gauss--Legendre quadrature with $n+1$ integration nodes is~$2n+1$,
    there exists $N$ sufficiently large such that $hat(I)_n (p) = I(p)$ for all $n >= N$.
    Thus
    $
      forall n >= N,
      wide
      Big(|hat(I)_n (u) - I(u)|)
      &<=
      Big(|hat(I)_n (u) - hat(I)_n (p)|)
      + Big(|hat(I)_n (p) - I(p)|)
      + Big(|I(p) - I(u)|) \
      &=
      Big(|hat(I)_n (u) - hat(I)_n (p)|)
      + Big(|I(p) - I(u)|) \
      &<= sum_(i=0)^n |w_i| E + integral_(-1)^1 E thin dd x
      = 2 epsilon + 2 epsilon = 4 epsilon.
    $
    Indeed, all the weights are positive
    and the Gauss--Legendre quadrature rule is exact for the constant function~$u(x) = 1$,
    which implies that the weights add up to 2.
    Since $epsilon > 0$ was arbitrary,
    this concludes the proof.
  ]
]

#paragraph[Computation of the integration nodes.][
  We proved that the integration nodes are given by the roots of the Legendre polynomials.
  In practice, calculating these roots can be achieved by calculating the eigenvalues of a tridiagonal matrix, see~@sub:orthogonal_integration.
  This is known as the Golub--Welsch algorithm.
]

#paragraph[Generalization to higher dimensions.][
  Gauss--Legendre integration is ubiquitous in numerical methods for partial differential equations,
  in particular the _finite element method_.
  Its generalization to higher dimensions is immediate:
  for a function $u colon [-1, 1] times [-1, 1] -> real$,
  we have
  $
    integral_0^1 integral_0^1 u(x, y) thin dd y dd x approx sum_(i=0)^n sum_(j=0)^n w_i w_j u(x_i, y_i).
  $
  The degree of precision of this integration rule is the same as
  that of the corresponding one-dimensional rule,
  and this approach can be generalized to any dimension~$d$.
  The associated computational cost,
  however, scales as $n^d$,
  and so it is not a good idea to use a deterministic method of this type in the high dimensional setting.
  The explosion of the computational cost as the dimension increases is known as~_the curse of dimensionality_.
]

== Introduction to probabilistic integration methods
So far in this chapter,
we covered only _deterministic_ integration formulas.
Much of the research around the calculation of high-dimensional integrals today is concerned with~_probabilistic_ integration methods using probabilistic approaches.
These methods are based on the connection between integrals and expectations.
To illustrate the simplest probabilistic approach,
called the _Monte Carlo method_,
we consider the problem of approximating the integral
$
  I = integral_0^1 u(x) thin dd x.
$
This integral may be expressed as the expectation $expect big([u(X)])$,
where~$expect$ is the expectation operator and~$X tilde.op cal(U)(0, 1)$ is a uniformly distributed random variable over the interval~$[0, 1]$.
Therefore, in practice,
the integral $I$ may be approximated by generating a large number of _independent_ samples $X_1, X_2, dots$ drawn from the distribution~$cal(U)(0, 1)$
and averaging $f(X_i)$ over all of these samples:
$
  hat(I)_N = 1/N sum_(n=1)^N u(X_n).
$
The quantity $hat(I)_N$,
where the subscript~$N$ denotes the number of samples employed,
is itself a random variable;
it is called an _estimator_ of the exact integral~$I$.
In Julia, the calculation of the integral using the Monte Carlo method can be achieved with the following code.
```julia
N = 1000
u(x) = x^2
X = rand(N)
Î = (1/N) * sum(u.(X))
```
Since the expectation operator is linear,
we calculate that
$
  expect [hat(I)_N] = 1/N sum_(n=1)^N expect big([u(X_n)]) = 1/N sum_(n=1)^N integral_0^1 u(x) thin dd x = 1/N sum_(n=1)^N I = I.
$
The estimator $hat(I)_N$ is therefore said to be _unbiased_,
because the bias $expect [hat(I)_N] - I$ is zero.
Now assume that the function $u$ is square integrable over the interval $[0, 1]$ and let
$
  sigma^2 := var big([u(X)]) := expect bigg([Big((u(X) - expect big([u(X)])))^2]) = integral_0^1 |u(x) - I|^2 thin dd x,
$
where $var [X]$ denotes the variance of the random variable $X$.
We can calculate explicitly the variance of the estimator $hat(I)_N$:
$
  var [hat(I)_N]
  = expect lr([lr((1/N sum_(n=1)^N u(X_n) - I))^2])
  = expect lr([1/N^2 sum_(n=1)^N sum_(m=1)^N big((u(X_n) - I)) big((u(X_m) - I))]).
$
Since the samples $X_1, X_2, dots$ are independent,
it holds that
$
  expect Big([big((u(X_n) - I)) big((u(X_m) - I))]) = delta_(n m) sigma^2,
$
where $delta_(m n)$ is the Kronecker delta.
Consequently, we deduce that
$
  var [hat(I)_N]
  &= 1/N^2 sum_(n=1)^N sum_(m=1)^N expect Big([big((u(X_n) - I)) big((u(X_m) - I))])
  = 1/N^2 sum_(n=1)^N sum_(m=1)^N delta_(m n) sigma^2 = sigma^2/N.
$
Therefore, the variance of $hat(I)_N$ decreases as $1\/N$ when the number of samples $N$ increases,
indicating that the estimator becomes increasingly accurate.
To state this more precisely,
we will use Chebyshev's inequality.
#theorem(title: [Chebyshev's inequality])[
  Let $Z$ be a random variable with mean $m$ and variance $s^2$.
  Then for any real $k > 0$,
  $
    proba big([|Z - m| >= k s]) <= 1/k^2.
  $
]
Let $epsilon > 0$.
Employing Chebyshev's inequality with $k = epsilon \/ sqrt(var [hat(I)_N])$,
we obtain
#numbered($
  proba Big([|hat(I)_N - I| >= epsilon])
  <= (var [hat(I)_N])/epsilon^2
  = sigma^2/(N epsilon^2).
$) <eq:chebyshev_monte_carlo>
Consequently, for any $epsilon > 0$,
the probability that the integration error $big(|hat(I)_N - I|)$ is greater than or equal to $epsilon$ tends to zero in the limit as $N -> infinity$.
In probabilistic jargon, we say that $hat(I)_N$ converges _in probability_ to~$I$.
Equation~@eq:chebyshev_monte_carlo can also be employed in order to construct a _confidence interval_ for the exact integral,
as we demonstrate in the following paragraph.

#paragraph[Construction of a $(1 - alpha)$ confidence interval.][
  By definition,
  a $(1 - alpha)$ confidence interval for~$I$ is
  an interval $[z_1, z_2]$,
  the endpoints of which being random variables,
  such that the probability that $I in [z_1, z_2]$ is greater than or equal to~$1 - alpha$.
  To construct such an interval,
  we begin by finding~$epsilon$ such that the right-hand side of~@eq:chebyshev_monte_carlo is equal to $alpha$,
  which gives $epsilon = sqrt(sigma^2 \/ N alpha)$.
  For this value of~$epsilon$,
  we have
  $
    proba Big([|hat(I)_N - I| >= epsilon]) <= alpha.
  $
  Since $|hat(I)_N - I| >= epsilon$ if and only if $I in.not (hat(I)_N - epsilon, hat(I)_N + epsilon)$,
  we conclude that
  $
    proba Big([I in lr((hat(I)_N - epsilon, hat(I)_N + epsilon))])
    >= 1 - alpha.
  $
  We have thus shown that
  $
    lr((hat(I)_N - sqrt(sigma^2/(N alpha)), hat(I)_N + sqrt(sigma^2/(N alpha))))
  $
  is a $(1 - alpha)$ confidence interval for $I$.
]

// \paragraph{Central limit theorem}
== Exercises
#exercise[
  Derive the Simpson's integration rule~@eq:simpsons.
]

#exercise[
  Derive the composite Simpson integration rule~@eq:composite_simpson.
] <exercise:composite_simpson>

#exercise[
  Consider the integration rule
  $
    integral_0^1 u(x) thin dd x approx w_1 u(0) + w_2 u(1) + w_3 u'(0).
  $
  Find $w_1$, $w_2$ and $w_3$ so that this integration rule has the highest possible degree of precision.
]

#exercise[
  Consider the integration rule
  $
    integral_(-1)^1 u(x) thin dd x approx w_1 u(x_1) + w_2 u'(x_1).
  $
  Find $w_1$, $w_2$ and $x_1$ so that this integration rule has the highest possible degree of precision.
]

#exercise[
  What is the degree of precision of the following quadrature rule?
  $
    integral_(-1)^1 u(x) thin dd x approx 2/3 lr((2 u lr((-1/2)) - u(0) + 2 u lr((1/2)))).
  $
]

#exercise[
  The Gauss--Hermite quadrature rule with~$n+1$ nodes is an approximation of the form
  $
    integral_(-infinity)^infinity u(x) thin ee^(- x^2/2) thin dd x approx sum_(i=0)^n w_i u(x_i),
  $
  such that the rule is exact for all polynomials of degree less than or equal to $2n+1$.
  Find the Gauss--Hermite rule with two nodes.
]

#exercise[
  Use Romberg's method to construct an integration rule with an error term scaling as~$bigo(h^4)$.
  Is there a link between the method you obtained and another integration rule seen in class?
]

#exercise(title: [Improving the error bound for the composite trapezoidal rule])[
  The notation used in this exercise is the same as in @sec:composite_methods.
  In particular, $hat(I)_h$ denotes the approximate integral obtained by using the composite trapezoidal rule~@eq:composite_trapezoidal_rule,
  and $hat(u)_h$ is the corresponding piecewise linear interpolant.

  A version of the mean value theorem states that,
  if $g colon [a, b] -> real$ is a non-negative integrable function
  and $f colon [a, b] -> real$ is continuous,
  then there exists $xi in (a, b)$ such that
  #numbered($
    integral_a^b f(x) g(x) thin dd x = f(xi) integral_a^b g(x) thin dd x.
  $) <eq:mean_value>
  - Using~@eq:mean_value,
    show that, for all $i in {0, dots, n-1}$,
    there exists $xi_i in (x_i, x_(i+1))$ such that
    $
      integral_(x_i)^(x_(i+1)) u(x) - hat(u)_h (x) thin dd x =
      - u''(xi_i) h^3/12.
    $

  - Prove, by using the intermediate value theorem,
    that if $f colon [a, b] -> real$ is a continuous function,
    then for any set $xi_0, dots, xi_(n-1)$ of points within the interval $(a, b)$,
    there exists $c in (a, b)$ such that
    $
      1/n sum_(i=0)^(n-1) f(xi_i) = f(c).
    $

  - Combining the previous items,
    conclude that there exists $xi in (a, b)$ such that
    $
      I - hat(I)_h = - u''(xi) (b-a) h^2/12,
    $
    which is a more precise expression of the error than that obtained in~@eq:error_bound_trapezoidal.
]
#remark[
  One may convince oneself of~@eq:mean_value by rewriting this equation as
  $
    (integral_a^b f(x) g(x) thin dd x)/(integral_a^b g(x) thin dd x) = f(xi).
  $
  The left-hand side is the average of $f(x)$ with respect to the probability measure
  with density given by
  $
    x |->
    g(x)/(integral_a^b g(x) thin dd x).
  $
]

#exercise(title: [From the final exam of Spring 2022])[
  Construct an integration rule of the form
  $
    integral_(-1)^1 u(x) thin dd x approx w_1 u lr((-1/2)) + w_2 u(0) + w_3 u lr((1/2))
  $
  with a degree of precision equal to at least 2.
  What is the degree of precision of the rule constructed?
]
#solution[
  The Lagrange polynomials associated with $-1\/2$, $0$ and $1\/2$ are respectively
  $
    p_1(x) &= 2x lr((x - 1/2)), \
    p_2(x) &= - 4 lr((x + 1/2)) lr((x - 1/2)), \
    p_3(x) &= 2 lr((x + 1/2)) x.
  $
  We deduce that
  $
    w_1 = integral_(-1)^1 p_1(x) = 4/3, wide
    w_2 = integral_(-1)^1 p_2(x) = - 2/3, wide
    w_3 = integral_(-1)^1 p_3(x) = 4/3.
  $
  By construction,
  the degree of precision is at least 2.
  However, the integration rule is exact also when $u(x) = x^3$.
  Since the rule is not exact for $u(x) = x^4$,
  we conclude that the degree of precision is 3.
]

#exercise(title: [From the final exam of Spring 2022])[
  The Gauss--Laguerre quadrature rule with~$n$ nodes is an approximation of the form
  $
    integral_0^infinity u(x) thin ee^(-x) thin dd x approx sum_(i=1)^n w_i u(x_i),
  $
  such that the rule is exact when $u$ is a polynomial of degree less than or equal to $2n-1$.
  - Find the Gauss--Laguerre rule with one node ($n = 1$).

  - Find the Gauss--Laguerre quadrature rule with two nodes ($n = 2$).
    You may find it useful to first calculate the Laguerre polynomial of degree 2.
]
#solution[
  Below are the derivations of the Gauss--Laguerre rules with 1 and 2 nodes.
  #paragraph[Gauss--Laguerre rule with 1 node.][
    We are looking for $w_1$ and $x_1$ such that
    $
      forall (a,b) in real^2, wide
      integral_0^infinity (a + b x) thin ee^(-x) thin dd x = w_1 (a + b x_1).
    $
    The left-hand side is equal to
    $
      a integral_0^infinity ee^(-x) thin dd x + b integral_0^infinity x ee^(-x) thin dd x
      = a + b integral_0^infinity x ee^(-x) thin dd x.
    $
    Using integration by parts,
    we find the value of the remaining integral on the right-hand side:
    $
      integral_0^infinity x ee^(-x)
      &= integral_0^infinity - (x ee^(-x))' + ee^(-x) thin dd x \
      &= - (x ee^(-x)) Big(|)_(x = infinity) + (x ee^(-x)) Big(|)_(x = 0) + integral_0^infinity ee^(-x) thin dd x \
      &= 0 + 0 + 1.
    $
    (To be fully rigorous, we would need to write the first term on the second line as a limit.)
    Therefore, we obtain
    $
      a + b = w_1(a + b x_1),
    $
    which implies that $w_1 = x_1 = 1$.
  ]

  #paragraph[Gauss--Laguerre rule with 2 nodes.][
    The integration nodes are given by the roots of the Laguerre polynomials,
    which are the orthogonal polynomials for the inner product
    $
      ip(f, g) :=
      integral_0^infinity f(x) g(x) thin ee^(-x) thin dd x.
    $
    The first polynomial is $ell_0(x) = 1$.
    It is simple to check that the only linear monic polynomial orthogonal to $ell_0$ is given by $ell_1(x) = x - 1$.
    Next, by integration by parts we calculate that
    $
      integral_0^infinity x^2 thin ee^(-x) thin dd x
      = integral_0^infinity - (x^2 ee^(-x))' + 2 x ee^(-x) thin dd x = 2.
    $
    and, similarly,
    $
      integral_0^infinity x^3 thin ee^(-x) thin dd x
      = integral_0^infinity - (x^3 ee^(-x))' + 3 x^2 ee^(-x) thin dd x = 6.
    $
    Consider the ansatz $ell_2(x) = x^2 + a ell_1(x) + b$.
    In order for $ell_2$ to be orthogonal to $ell_0$ and~$ell_1$,
    it is necessary that
    $
      0 &= integral_0^infinity ell_2(x) thin ell_0(x) thin ee^(-x) thin dd x = 2 + b, \
      0 &= integral_0^infinity ell_2(x) thin ell_1(x) thin ee^(-x) thin dd x
      = 4 + a integral_0^infinity ell_1(x) ell_1(x) thin space d x = 4 + a.
    $
    Therefore, we conclude that $a = -4$ and $b=-2$,
    which gives
    $
      ell_2(x) = x^2 - 4 x + 2.
    $
    The roots are given by $2 plus.minus sqrt(2)$,
    so we have $x_1 = 2 - sqrt(2)$ and $x_2 = 2 + sqrt(2)$.
    It remains to find the weights.
    To this end, we need only two additional equations;
    it is sufficient to require that, for any $(a,b) in real^2$,
    $
      a + b =
      integral_0^infinity (a + b x) thin ee^(-x) thin dd x
      &= w_1 (a + b x_1) + w_2 (a + b x_2) \
      &= a (w_1 + w_2) + 2 b (w_1 + w_2) + sqrt(2) b (w_2 - w_1).
    $
    Letting $a = 1$ and $b = 0$,
    we obtain $w_1 + w_2 = 1$.
    Then, letting $a = 0$ and $b = 1$,
    we deduce
    $
      1 = 2 + sqrt(2) (w_2 - w_1) wide <=> w_2 - w_1 = - sqrt(2)/2.
    $
    Therefore
    $
      w_1 = (2 + sqrt(2))/4, wide w_2 = (2 - sqrt(2))/4,
    $
    which concludes the exercise.
  ]
]

#compexercise(title: [Calculating the volume of hyperballs])[
  Let $B_d$ denote the $d$-dimensional unit ball for the Euclidean norm:
  $
    B_d = Big(\{) vect(x) in real^d : norm(vect(x)) <= 1 Big(\}).
  $
  The volume of $B_d$ is defined as the integral of the characteristic function over~$B_d$:
  $
    "vol"(B_d) = underbrace(integral_real dots integral_real, d "times") chi(vect(x)) thin dd x_1 dots dd x_d,
    wide
    chi(vect(x)) :=
    cases(
      1 & "if " vect(x) in B_d,
      0 & "otherwise."
    )
  $
  Complete the following tasks:
  - Write a function
    ```julia
    function hyperball_volume(dim, n)
        # Your code comes here
        return vol, σ
    end
    ```
    that calculates the volume of the unit ball $B_d$ with $d = #raw("dim")$
    using a Monte Carlo approach with~`n` samples drawn from an appropriate distribution.
    Your function should return an estimation of the volume
    together with the standard deviation of the estimator
    (which you should estimate from the samples).

  - Using the function `hyperball_volume`,
    plot the volumes for $d$ going from 1 to 15,
    together with a 99% confidence interval.
    See~@fig:hyperballs for an example solution with $n = 10^7$.
  You are allowed to use your knowledge of the fact that
  $"vol"(B_2) = pi$ and $"vol"(B_3) = 4 pi\/3$,
  but do not use the general formula for the volume of $B_d$.
] <exercise:hyperballs>
#figure(
  placement: top,
  image("/figures/hyperballs.pdf", width: 80%),
  caption: [Example solution for~@exercise:hyperballs.],
) <fig:hyperballs>

#compexercise[
  Complete the following tasks:
  - Write a function `legendre(n)` which returns the Legendre polynomial of degree~$n$.
    To this end, you may use the `Polynomials` library and Rodrigues' formula:
    $
      L_n (x) = 1/(2^n n!) (dd^n)/(dd x^n) lr((x^2 - 1))^n.
    $

  - Write a function `get_nodes_and_weights(n)` which returns the nodes and weights of the Gauss--Legendre quadrature with~$n$ nodes.
    In order to construct Lagrange polynomials,
    you may find it useful to use the `fromroots` functions.

  - Write a function `composite_gauss_legendre(u, a, b, n, N)`,
    which returns an approximation of the integral
    $
      integral_a^b u(x) thin dd x
    $
    obtained by partitioning the integration interval $[a, b]$ into $N$ cells,
    and applying the Gauss--Legendre quadrature within each cell.

  - Take $u(x) = cos(x)$, $a = -1$ and $b = 1$.
    Illustrate on the same plot the error for the values~$n in {1, 2, 3}$ and $N$ varying from 1 to 40.
    For each value of~$n$, estimate the order of convergence with respect to $N$, i.e.~find~$alpha(n)$ such that
    $
      |hat(I)_(n, N) - I| prop C N^(-alpha),
    $
    where $I$ denotes the exact value of the Integral
    and $hat(I)_(n, N)$ denotes its approximation.
]

== Discussion and bibliography
The presentation of part of the material follows that in~@Legat,
and some exercises come from~@MR2265914[Chapter 9].
The main advantage of probabilistic integration approaches is that
they generalize naturally to high-dimensional and infinite-dimensional settings.

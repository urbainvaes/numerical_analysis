#import "../template.typ": *

// \ip*{...} (delimiters scaled to the content) and \ip[big]{...}
#let ipstar(..args) = $lr(chevron.l #args.pos().join($,$) chevron.r)$
#let ipbig(..args) = big($chevron.l #args.pos().join($,$) chevron.r$)

= Interpolation and approximation <cha:interpolation_and_approximation>
#minitoc()

#pagebreak()
#section-star[Introduction]
In this chapter,
we study numerical methods for interpolating and approximating functions.
The Cambridge dictionary defines interpolation as _the addition of something different in the middle of a text, piece of music, etc.~or the thing that is added_.
The concept of interpolation in mathematics is consistent with this definition;
interpolation consists in finding, given a set of points~$(x_i, y_i)$,
a function~$f$ in a finite-dimensional space that goes through these points.
Throughout this course, you use the~`plot` function in Julia,
which performs piecewise linear interpolation for drawing functions,
but there are a number of other standard interpolation methods.
Our first goal in this chapter is to present an overview of these methods and the associated error estimates.

In the second part of this chapter,
we focus on _function approximation_,
which is closely related to the subject of mathematical interpolation.
Indeed, a simple manner for approximating a general function by another one in a finite-dimensional space is to select a set of real numbers on the $x$ axis,
called _nodes_, and find the associated interpolant.
As we shall demonstrate, not all sets of interpolation nodes are equal,
and special care is required in order to avoid undesired oscillations.
The field of function approximation is vast,
so our aim in this chapter is to present only an introduction to the subject.
In order to quantify the quality of an approximation,
a metric on the space of functions,
or a subset thereof, must be specified in order to measure errors.
Without a metric, saying that two functions are close is almost meaningless!

== Interpolation
Assume that we are given $n+1$ nodes $x_0, dots, x_n$ on the $x$ axis,
together with values $u_0, dots, u_n$,
which may be the values taken by an unknown function~$u(x)$ when evaluated at these points.
Suppose that we are looking for an interpolation~$hat(u)(x)$ in a subspace~$Span {phi_0, dots, phi_n}$
of the vector space of continuous functions, i.e.~an interpolating function of the form
$ hat(u)(x) = alpha_0 phi_0(x) + cdots + alpha_n phi_n (x), $
where $alpha_0, dots, alpha_n$ are real coefficients.
In order for~$hat(u)(x)$ to be an interpolating function,
we must require that
$ forall i in {0, dots, n}, wide hat(u)(x_i) = u_i. $
This leads to a linear system of $n+1$ equations and $n+1$ unknowns,
the latter being the coefficients~$alpha_0, dots, alpha_n$.
This system of equations in matrix form reads
#numbered($
  pmat(
    phi_0(x_0), phi_1(x_0), …, phi_n (x_0);
    phi_0(x_1), phi_1(x_1), …, phi_n (x_1);
    ⋮, ⋮, , ⋮;
    phi_0(x_n), phi_1(x_n), …, phi_n (x_n)
  )
  pmat(alpha_0; alpha_1; ⋮; alpha_n)
  =
  pmat(u_0; u_1; ⋮; u_n).
$) <eq:linear_system_interpolation>

=== Vandermonde matrix
Since polynomials are very convenient for evaluation, integration, and differentiation,
they are a natural choice for interpolation purposes.
The simplest basis of the subspace of polynomials of degree less than or equal to $n$ is given by the monomials:
$ phi_0(x) = 1, wide phi_1(x) = x, wide dots, wide phi_n (x) = x^n. $
In this case,
the linear system~@eq:linear_system_interpolation for determining the coefficients of the interpolant reads
#numbered($
  pmat(
    1, x_0, …, x_0^n;
    1, x_1, …, x_1^n;
    ⋮, ⋮, , ⋮;
    1, x_n, …, x_n^n
  )
  pmat(alpha_0; alpha_1; ⋮; alpha_n)
  =
  pmat(u_0; u_1; ⋮; u_n).
$) <eq:linear_system_interpolation_poly>
The matrix on the left-hand side is called a _Vandermonde_ matrix.
If the abcissae $x_0, dots, x_n$ are distinct,
then this is a full rank matrix,
and so~@eq:linear_system_interpolation_poly admits a unique solution,
implying as a corollary that the interpolating polynomial exists and is unique.
It is possible to show that the condition number of the Vandermonde increases dramatically with $n$.
Consequently, solving~@eq:linear_system_interpolation_poly is not a viable method in practice for calculating the interpolating polynomial.

=== Lagrange interpolation formula
One may wonder whether polynomial basis functions $phi_0, dots, phi_n$ can be defined in such a manner that
the matrix in~@eq:linear_system_interpolation is the identity matrix.
The answer to this question is positive;
it suffices to take as a basis the _Lagrange polynomials_,
which are given by
$
  phi_i (x)
  = ((x - x_0) (x - x_1) dots (x - x_(i-1)) (x - x_(i+1)) dots (x - x_n))
  / ((x_i - x_0) (x_i - x_1) dots (x_i - x_(i-1)) (x_i - x_(i+1)) dots (x_i - x_n))
  = product_(j = 0 \ j != i)^n (x - x_j)/(x_i - x_j).
$
It is simple to check that
$
  phi_i (x_j) =
  delta_(i,j) =
  cases(
    1 & quad "if" i = j\,,
    0 & quad "otherwise.",
  )
$
Finding the interpolant in this basis is immediate:
$ hat(u)(x) = u_0 phi_0(x) + cdots + u_n phi_n (x). $
While simple, this approach to polynomial interpolation has a few disadvantages:
- First, evaluating $hat(u)(x)$ is computationally costly when $n$ is large.

- Second, all the basis functions change when adding new interpolation nodes.

- Finally, Lagrange interpolation is numerically unstable because of cancellations between large terms.
  Indeed, it is often the case that Lagrange polynomials take very large values over the interpolation intervals;
  this occurs, for example,
  when many equidistant interpolation nodes are employed,
  as illustrated in~@fig:lagrange.

#figure(
  placement: top,
  stack(
    image("/figures/lagrange_10.pdf", width: 80%),
    image("/figures/lagrange_20.pdf", width: 80%),
  ),
  caption: [Lagrange polynomials associated with equidistant nodes over the $(0, 1)$ interval.],
) <fig:lagrange>

=== Gregory--Newton interpolation
By Taylor's formula,
any polynomial~$p$ of degree $n$ may be expressed as
#numbered($
  p(x) = p(0) + p'(0) x + (p''(0))/2 x^2 + dots + (p^((n))(0))/(n!) x^n.
$) <eq:taylor>
The constant coefficient can be obtained by evaluating the polynomial at 0,
the linear coefficient can be identified by evaluating the first derivative at 0,
and so on.
Assume now that we are given the values taken by~$p$ when evaluated at the integer numbers ${0, dots, n}$.
We ask the following question:
can we find a formula similar in spirit to~@eq:taylor,
but including only evaluations of $p$ and not of its derivatives?
To answer this question, we introduce the difference operator $Delta$ which acts on functions as follows:
$ Delta f(x) = f(x+1) - f(x). $
The operator~$Delta$ is a linear operator on the space of continuous functions.
It maps constant functions to 0,
and the linear function~$x$ to the constant function~$1$,
suggesting a resemblance with the differentiation operator.
In order to further understand this connection,
let us define the _falling power_ of a real number $x$ as
#numbered($
  x^(underline(k)) = x (x-1) (x-2) dots (x-k+1).
$) <eq:falling_powers>
We then have that
#eq-last($
  Delta x^(underline(k))
  &= (x+1) x (x-1) dots (x-k+2) - x (x-1) (x-2) dots (x-k+1) \
  &= big(((x+1) - (x-k+1))) big((x (x-1) dots (x-k+2))) = k x^(underline(k-1))
$) <eq:falling_powers_rule>
In other words,
the action of the difference operator on falling powers mirrors that of the differentiation operator on monomials.
The falling powers form a basis of the space of polynomials,
and so any polynomial in~$poly(n)$, i.e.~of degree less than or equal to $n$, can be expressed as
#numbered($
  p(x) = alpha_0 + alpha_1 x^(underline(1)) + alpha_2 x^(underline(2)) + cdots + alpha_n x^(underline(n)).
$) <eq:gregory_newton>
It is immediate to show that $alpha_i = Delta^i p(0) \/ i!$,
where $Delta^i p$ denotes the function obtained after~$i$ applications of the operator~$Delta$.
Therefore, any polynomial of degree less than or equal to~$n$ may be expressed as
#numbered($
  p(x) = p(0) + Delta p(0) x^(underline(1)) + (Delta^2 p(0))/2 x^(underline(2)) + cdots + (Delta^n p(0))/(n!) x^(underline(n)).
$) <eq:gregory_newton_actually>
An expansion of the form~@eq:gregory_newton_actually is called a _Newton series_,
which is the discrete analog of the continuous Taylor series.
From the definition of~$Delta$,
it is clear that the coefficients in~@eq:gregory_newton_actually depend only on~$p(0), dots, p(n)$.
We conclude that, given points $n+1$ points $(i, u_i)$ for $i in {0, dots, n}$,
the unique interpolating polynomial is given by~@eq:gregory_newton_actually,
after replacing $p(i)$ by $u_i$.

#example[
  Let us use~@eq:gregory_newton in order to calculate the value of
  $ S(n) := sum_(i=0)^n i^2. $
  Since $Delta S(n) = (n+1)^2$,
  which is a second degree polynomial in $n$,
  we deduce that $S(n)$ is a polynomial of degree 3.
  Let us now determine its coefficients.
  #align(center, tabular(
    columns: 5,
    [$n$], [$0$], [$1$], [$2$], [$3$],
    [$Delta^0 S(n)$], [$bold(0)$], [$1$], [$5$], [$14$],
    [$Delta^1 S(n)$], [$bold(1)$], [$4$], [$9$], [],
    [$Delta^2 S(n)$], [$bold(3)$], [$5$], [], [],
    [$Delta^3 S(n)$], [$bold(2)$], [], [], [],
  ))
  We conclude that
  $
    S(n) = bold(1) n + bold(3)/(2!) n(n-1) + bold(2)/(3!) n(n-1)(n-2)
    = (n (2n+1) (n+1))/6
  $
] <example:gregory_newton>
Notice that when falling powers are employed as polynomial basis,
the matrix in~@eq:linear_system_interpolation is lower triangular,
and so the algorithm described in~@example:gregory_newton could be replaced by the forward substitution method.
Whereas the coefficients of the Lagrange interpolant can be obtained immediately from the values of~$u$ at the nodes,
calculating the coefficients of the expansion in~@eq:gregory_newton requires $bigo(n^2)$ operations.
However, Gregory--Newton interpolation has several advantages over Lagrange interpolation:
- If a point~$(n+1, p_(n+1))$ is added to the set of interpolation points,
  only one additional term, corresponding to the falling power~$x^(underline(n+1))$,
  needs to be calculated in~@eq:gregory_newton_actually.
  All the other coefficients are unchanged.
  Therefore, the Gregory--Newton approach is well-suited for incremental interpolation.

- The Gregory--Newton interpolation method is more numerically stable than Lagrange interpolation,
  because the basis functions do not take very large values.

- A polynomial in the form of a Newton series can be evaluated efficiently using Horner's method,
  which is based on rewriting the polynomial as
  $
    p(x) = alpha_0 + x bigg((alpha_1 + (x-1) Big((alpha_2 + (x-2) big((alpha_3 + (x-3) dots)))))).
  $
  Evaluating this expression starting from the innermost bracket leads to an algorithm with a cost scaling linearly with the degree of the polynomial.

==== Non-equidistant nodes
So far,
we have described the Gregory--Newton method in the simple setting where interpolation nodes are just a sequence of successive natural numbers.
The method can be generalized to the setting of nodes $x_0 != dots != x_n$ which are not necessarily equidistant.
In this case, we take as basis the following functions instead of the falling powers:
#numbered($
  phi_i (x) = (x - x_0) (x - x_1) dots (x - x_(i-1)),
$) <eq:basis_newton>
with the convention that the empty product is 1.
By~@eq:linear_system_interpolation,
the coefficients of the interpolating polynomial in this basis solve the following linear system:
#numbered($
  pmat(
    1, , dots, , 0;
    1, x_1-x_0, , , ;
    1, x_2-x_0, (x_2-x_0)(x_2-x_1), , ⋮;
    ⋮, ⋮, , ⋱, ;
    1, x_n-x_0, dots, dots, product_(j=0)^(n-1)(x_n - x_j)
  )
  pmat(alpha_0; alpha_1; alpha_2; ⋮; alpha_n)
  =
  pmat(u_0; u_1; u_2; ⋮; u_n).
$) <eq:matrix_newton>
This system could be solved using, for example, forward substitution.
Clearly $alpha_0 = u_0$ from the first equation,
and then from the second equation we obtain
$ alpha_1 = (u_1 - u_0)/(x_1 - x_0) =: [u_0, u_1], $
which may be viewed as an approximation of the slope of~$u$ at $x_0$.
The right-hand side of this equation is an example of a _divided difference_.
In general, divided differences are defined recursively as follows:
#numbered($
  [u_0, u_1, dots, u_d] := ([u_1, dots, u_d] - [u_0, dots, u_(d-1)])/(x_d - x_0), wide [u_i] = u_i.
$) <eq:definition_divided_difference>

Let us start by observing that our divided differences coincide
with the coefficients in~@eq:gregory_newton.

#lemma[
  Assume that $(0, u_0), dots, (n, u_n) in real^2$.
  Then, provided that $(i, k) in nat^2$ satisfies $i + k <= n$,
  we have
  $ 1/(k!) Delta^k u(i) = [u_i, dots, u_(i + k)]. $
] <lemma:gregory_newton_with_divided_differences>
#proof[
  The result is clearly true for $k = 0$,
  and we shall assume by induction that it holds for~$k - 1$.
  Using our induction hypothesis, we get
  $
    1/(k!) Delta^k u(i)
    &= 1/k lr(((Delta^(k - 1) u(i + 1))/((k - 1)!) - (Delta^(k - 1) u(i))/((k - 1)!))) \
    &= 1/k Big(([u_(i + 1), dots, u_(i + k)] - [u_i, dots, u_(i + k - 1)])).
  $
  By the definition of divided differences in~@eq:definition_divided_difference,
  it follows that we have
  $
    1/(k!) Delta^k u(i)
    = 1/k (underbrace(x_(i + k) - x_i, i + k - i)) [u_i, dots, u_(i + k)]
    = [u_i, dots, u_(i + k)],
  $
  which allows us to conclude the proof.
]

#indent In light of @lemma:gregory_newton_with_divided_differences,
the reader will not be surprised that
the coefficients in the $phi_0, dots, phi_n$ basis in~@eq:basis_newton
are given by the divided differences.

#proposition[
  Assume that $(x_0, u_0), dots, (x_n, u_n)$ are $n+1$ points in the plane with distinct abcissae.
  Then the interpolating polynomial of degree $n$ may be expressed as
  $ p(x) = sum_(i=0)^n [u_0, dots, u_i] phi_i (x), $
  where $phi_i (x)$, for $i = 0, dots, n$, are the basis functions defined in~@eq:basis_newton.
]
#proof[
  The statement is true for $n = 0$.
  Reasoning by induction, we assume that it holds true for polynomials of degree up to $n-1$.
  Let $p_1(x)$ and $p_2(x)$ be the interpolating polynomials at the points
  $x_0, x_1, dots, x_(n-2), x_(n-1)$ and $x_0, x_1, dots, x_(n-2), x_n$, respectively.
  Then
  #numbered($
    p(x) = p_1(x) + (x - x_(n-1))/(x_n - x_(n-1)) big((p_2(x) - p_1(x)))
  $) <eq:interpolating_polynomial>
  is a polynomial of degree~$n$ that interpolates all the data points.
  By the induction hypothesis,
  it holds that
  $
    p_1(x) &= u_0 + [u_0, u_1] (x - x_0) + dots + [u_0, u_1, dots, u_(n-2), bold(upright(u_(n-1)))] product_(i=0)^(n-2) (x - x_i), \
    p_2(x) &= u_0 + [u_0, u_1] (x - x_0) + dots + [u_0, u_1, dots, u_(n-2), bold(upright(u_n))] product_(i=0)^(n-2) (x - x_i).
  $
  Here we used a bold font in order to emphasize the difference between the two expressions.
  Substituting these expressions in~@eq:interpolating_polynomial,
  we obtain
  $
    p(x) =
    & u_0 + [u_0, u_1] (x - x_0) + dots + [u_0, dots, u_(n-2)] product_(i=0)^(n-2) (x - x_i) \
    & + ([u_0, u_1, dots, u_(n-2), u_n] - [u_0, u_1, dots, u_(n-2), u_(n-1)])/(x_n - x_(n-1)) product_(i=0)^(n-1) (x - x_i).
  $
  In @exercise:divided_differences,
  we show that divided differences are invariant under permutations of the data points,
  and so we have that
  $
    ([u_0, u_1, dots, u_(n-2), u_n] - [u_0, u_1, dots, u_(n-2), u_(n-1)])/(x_n - x_(n-1)) = [u_0, dots, u_n],
  $
  which enables to conclude.
]
#example[
  Assume that we are looking for the third degree polynomial going through the following points:
  $ (-1, 10), wide (0, 4), wide (2, -2), wide (4, -40). $
  We have to calculate the divided difference $alpha_i = [u_0, dots, u_i]$ for $i in {0, 1, 2, 3}$.
  To this end,
  it is convenient to use a table:
  #align(center, tabular(
    columns: 5,
    [$i$], [$0$], [$1$], [$2$], [$3$],
    [$[u_i]$], [$bold(10)$], [$4$], [$-2$], [$-40$],
    [$x_(i+1) - x_i$], [$1$], [$2$], [$2$], [],
    [$[u_i, u_(i+1)]$], [$bold(-6)$], [$-3$], [$-19$], [],
    [$x_(i+2) - x_i$], [$3$], [$4$], [], [],
    [$[u_i, u_(i+1), u_(i+2)]$], [$bold(1)$], [$-4$], [], [],
    [$x_(i+3) - x_i$], [$5$], [], [], [],
    [$[u_i, u_(i+1), u_(i+2), u_(i+3)]$], [$bold(-1)$], [], [], [],
  ))
  We deduce that the expression of the interpolating polynomial is
  $
    p(x)
    = bold(10) + (bold(-6))(x+1) + bold(1) (x+1)x + (bold(-1)) (x+1)x(x-2)
    = - x^3 + 2 x^2 + -3x + 4.
  $
]

=== Interpolation error
Assume that~$u(x)$ is a continuous function and denote by $hat(u)(x)$ its interpolation through the points~$(x_i, u_i)$,
where $u_i = u(x_i)$ for $i = 0, dots, n$.
In this section, we study the behavior of the error in the limit as $n -> infinity$.

#theorem(title: [Interpolation error])[
  Assume that~$u colon [a, b] -> real$ is a function in $C^(n+1)([a, b])$ and let~$x_0, dots, x_n$ denote $n+1$ distinct interpolation nodes.
  Then for all $x in [a, b]$,
  there exists~$xi = xi(x)$ in the interval $[a, b]$ such that
  $
    e_n (x) := u(x) - hat(u)(x) = (u^((n+1))(xi))/((n+1)!) (x-x_0) dots (x - x_n).
  $
] <theorem:interpolation_error>
#proof[
  The statement is obvious if $x in {x_0, dots, x_n}$,
  so we assume from now on that $x$ does not coincide with an interpolation node.
  Let us use the compact notation $omega_n = product_(i=0)^n (x - x_i)$ and introduce the function
  #numbered($
    g(t) = e_n (t) omega_n (x) - e_n (x) omega_n (t).
  $) <eq:interpolation_error>
  The function $g$ is smooth and takes the value~0 when evaluated at~$x_0, dots, x_n, x$.
  Since $g$ has $n+2$ roots in the interval $[a, b]$,
  Rolle's theorem implies that $g'$ has at least $n+1$ distinct roots in~$(a, b)$.
  Another application of Rolle's theorem then yields that $g''$ has at least $n$ distinct roots in~$(a, b)$.
  Iterating this reasoning, we deduce that $g^((n+1))$ has one root~$t_*$ in $(a, b)$.
  From~@eq:interpolation_error,
  we calculate that
  #numbered($
    g^((n+1))(t) = e_n^((n+1))(t) omega_n (x) - e_n (x) omega_n^((n+1))(t)
    = u^((n+1))(t) omega_n (x) - e_n (x) (n+1)!.
  $) <eq:interpolation_error_proof>
  Here we employed the fact that $hat(u)^((n+1))(t) = 0$,
  because $hat(u)$ is a polynomial of degree at most $n$.
  Evaluating~@eq:interpolation_error_proof at~$t_*$ and rearranging,
  we obtain that
  $ e_n (x) = (u^((n+1))(t_*))/((n+1)!) omega_n (x), $
  which completes the proof.
]
As a corollary to~@theorem:interpolation_error,
we deduce the following error bound.
#corollary(title: [Upper bound on the interpolation error])[
  Assume that $u$ is smooth in the interval $[a, b]$ and let
  $ C_(n+1) = sup_(x in [a, b]) lr(|u^((n+1))(x)|). $
  Then
  #numbered($
    E_n := sup_(x in [a, b]) big(|e_n (x)|) <= C_(n+1)/(4(n+1)) h^(n+1)
  $) <eq:upper_bound_interp_error>
  where $h$ is the maximum spacing between two successive interpolation nodes.
] <corollary:interpolation_error>
#proof[
  By~@theorem:interpolation_error,
  it holds that
  #numbered($
    forall x in [a, b], wide
    abs(e_n (x)) <= C_(n+1)/((n+1)!) Big(|(x-x_0) dots (x-x_n)|).
  $) <eq:bound_error_from_theorem>
  The product on the right-hand side is bounded from above by
  #numbered($
    h^2/4 times 2h times 3h times 4h times cdots times n h = (n! h^(n+1))/4.
  $) <eq:interpolation_bound_product>
  The first factor comes from the fact that, if $x in [x_i, x_(i+1)]$,
  then
  $ Big(|(x - x_i)(x - x_(i+1))|) <= (x_(i+1) - x_i)^2/4, $
  because the left-hand side is maximized when $x$ is the midpoint of the interval~$[x_i, x_(i+1)]$.
  Substituting~@eq:interpolation_bound_product into~@eq:bound_error_from_theorem,
  we deduce the statement.
]
We now ask the following natural question:
does $E_n$ given in~@eq:upper_bound_interp_error tend to zero as the maximum spacing between successive nodes tends to 0?
By~@corollary:interpolation_error,
the answer to this question is positive when $C_n$ does not grow too quickly as $n -> infinity$.
For example
the interpolation error for the function $u(x) = sin(x)$
decreases very quickly as $n -> infinity$
when equidistant interpolation nodes are employed,
as illustrated in @fig:interpolation_sine.
#figure(
  placement: top,
  grid(
    columns: (49%, 49%), column-gutter: 0.33em,
    image("/figures/interpolation_sine_3.pdf", width: 100%),
    image("/figures/interpolation_sine_4.pdf", width: 100%),
    image("/figures/interpolation_sine_6.pdf", width: 100%),
    image("/figures/interpolation_sine_8.pdf", width: 100%),
  ),
  caption: [Interpolation (in orange) of the function $u(x) = sin(x)$ (in blue) using 3, 4, 6, and 8 equidistant nodes.],
) <fig:interpolation_sine>

In some cases, however,
the constant $C_n$ grows quickly with $n$,
to the extent that $E_n$ may increase with $n$;
in this case, the maximum interpolation error grows when nodes are added!
The classic example illustrating this potential issue is that of the Runge function:
#numbered($
  u(x) = 1/(1 + 25 x^2).
$) <eq:runge_function>
It is possible to show that,
for this function,
the upper bound in~@eq:upper_bound_interp_error tends to $infinity$ in the limit as the number~$n$ of interpolation nodes increases.
We emphasize that this does not prove that~$E_n -> infinity$ in the limit as $n -> infinity$,
because~@eq:upper_bound_interp_error provides only an _upper bound_ on the error.
In fact, the interpolation error for the Runge function can either grow or decrease,
depending on the choice of interpolation nodes.
With equidistant nodes, it turns out that $E_n -> infinity$,
as illustrated in~@fig:interpolation_runge_function.
#figure(
  placement: top,
  grid(
    columns: (49%, 49%), column-gutter: 0.33em,
    image("/figures/interpolation_runge_6.pdf", width: 100%),
    image("/figures/interpolation_runge_10.pdf", width: 100%),
    image("/figures/interpolation_runge_14.pdf", width: 100%),
    image("/figures/interpolation_runge_20.pdf", width: 100%),
  ),
  caption: [Interpolation (in orange) of the Runge function~@eq:runge_function (in blue) using 6, 10, 14, and 20 equidistant nodes.],
) <fig:interpolation_runge_function>

=== Interpolation at Chebyshev nodes
Sometimes,
interpolation is employed as a technique for approximating functions.
The spectral collocation method, for example,
is a technique for solving partial differential equations
where a discrete solution is first calculated,
and then a continuous solution is constructed using polynomial or Fourier interpolation.
When the interpolation nodes are not given a priori as data,
it is natural to wonder whether these can be picked in such a manner that the interpolation error,
measured in a function norm, is minimized.
For example, given a continuous function~$u(x)$ and a number of nodes~$n+1$,
is it possible to choose nodes $x_0, dots, x_n$ such that
$ E := sup_(x in [a, b]) big(|u(x) - hat(u)(x)|) $
is minimized?
Here $hat(u)$ is the polynomial interpolating~$u$ at the nodes.
Achieving this goal in general is a difficult task,
because $xi = xi(x)$ is unknown in the expression of the interpolation error from~@theorem:interpolation_error:
$ e_n (x) = (u^((n+1))(xi))/((n+1)!) (x-x_0) dots (x - x_n). $
In view of this difficulty,
we will focus on the simpler problem of finding interpolation nodes
such that the product $(x-x_0) dots (x - x_n)$ is minimized in the supremum norm.
This problem amounts to finding the optimal interpolation nodes,
in the sense that~$E$ is minimized,
in the particular case where $u$ is a polynomial of degree $n+1$,
because in this case $u^((n+1))(xi)$ is a constant factor.
It turns out that this problem admits an explicit solution,
which we will deduce from the following theorem.
#theorem(title: [Minimum $infinity$ norm])[
  Assume that $p$ is a monic polynomial of degree $n >= 1$:
  $ p(x) = alpha_0 + alpha_1 x + cdots + alpha_(n-1) x^(n-1) + x^n. $
  Then it holds that
  #numbered($
    sup_(x in [-1, 1]) big(|p(x)|) >= 1/(2^(n-1)) =: E.
  $) <eq:chebychev_lower_bound>
  In addition, the lower bound is achieved for~$p_*(x) = 2^(-(n-1)) T_n (x)$,
  where $T_n$ is the Chebyshev polynomial of degree~$n$:
  #numbered($
    T_n (x) = cos(n arccos x) wide (-1 <= x <= 1).
  $) <eq:chebyshev_polynomial>
] <theorem:minimum_infty_norm>
#proof[
  By~@exercise:chebyshev_leading_coefficient,
  the polynomial $x |-> 2^(-(n-1)) T_n (x)$ is indeed monic,
  and it is clear that it achieves the lower bound~@eq:chebychev_lower_bound
  since $abs(T_n (x)) <= 1$ for all $x in [-1, 1]$.

  In order to prove~@eq:chebychev_lower_bound,
  we assume by contradiction that there is a different monic polynomial~$q$ of degree~$n$ such that
  #numbered($
    sup_(x in [-1, 1]) big(|q(x)|) < E.
  $) <eq:chebyshev_minimum>
  Let us introduce $x_i = cos(i pi\/n)$, for $i = 0, dots, n$,
  and observe that
  $ p(x_i) = 2^(-(n-1)) T_n (x_i) = (-1)^i E. $
  The function $h(x) := p(x) - q(x)$ is a polynomial of degree at most $n-1$ which,
  by the assumption~@eq:chebyshev_minimum,
  is strictly positive at $x_0, x_2, x_4, dots$ and strictly negative at $x_1, x_3, x_5, dots$.
  Therefore, the polynomial $h(x)$ changes sign at least $n$ times and so,
  by the intermediate value theorem, it has at least $n$ roots.
  But this is impossible, because $h(x) != 0$ and $h(x)$ is of degree at most $n-1$.
]
#remark(title: [Derivation of Chebyshev polynomials])[
  The polynomial~$p_*$ achieving the lower bound in~@eq:chebychev_lower_bound
  oscillates between the values $-E$ and $E$,
  which are respectively its minimum and maximum values over the interval~$[-1, 1]$.
  It attains the values $E$ or $-E$ at $n+1$ distinct points $x_0 < dots < x_n$,
  with $x_0 = -1$ and $x_n = 1$.
  It turns out that these properties,
  which can be shown to hold a priori using Chebyshev's _equioscillation theorem_,
  are sufficient to derive an explicit expression for the polynomial~$p_*$,
  as we formally demonstrate hereafter.

  The points $x_1, dots, x_(n-1)$
  are local extrema of~$p_*$,
  and so~$p'_*(x) = 0$ at these nodes.
  We therefore deduce that $p_*$ satisfies the differential equation
  #numbered($
    n^2 lr((E^2 - p_*(x)^2)) = p'_*(x)^2 (1 - x^2).
  $) <eq:differential_equation_chebyshev>
  Indeed, both sides are polynomials of degree $2n$ with single roots at -1 and 1,
  with double roots at $x_1, dots, x_(n-1)$,
  and with the same coefficient of the leading power.
  In order to solve~@eq:differential_equation_chebyshev,
  we rearrange the equation and take the square root:
  $
    (p'_*(x)/E)/sqrt(1 - (p_*(x)^2)/E^2) = plus.minus n/sqrt(1 - x^2)
    wide <==> wide
    dd/(dd x) lr((arccos lr(((p_*(x))/E)))) = plus.minus n dd/(dd x) arccos(x).
  $
  Integrating both sides and taking the cosine,
  we obtain
  $ p_*(x) = E cos big((C + n arccos(x))). $
  Requiring that $|p_*(-1)| = E$, we deduce $C = 0$.
] <remark:cheb>

#indent From~@theorem:minimum_infty_norm,
we deduce the following corollary.
#corollary(title: [Chebyshev nodes])[
  Assume that $x_0 < x_1 < dots < x_n$ are in the interval $[a, b]$.
  The supremum norm of the product $omega(x) := (x-x_0) cdots (x-x_n)$ over~$[a, b]$ is minimized when
  #numbered($
    x_i = a + (b-a) (1 + cos lr((((2i + 1) pi)/(2(n+1)))))/2
  $)
] <corollary:chebyshev_nodes>
#proof[
  We consider the affine change of variable
  $
    zeta colon & [-1, 1] -> [a, b]; \
               & y |-> (a + b + y (b - a))/2.
  $
  The function
  $
    p(y) &:= 2^(n+1)/(b-a)^(n+1) omega big((zeta(y))) = 2^(n+1)/(b-a)^(n+1) big((zeta(y)-x_0)) cdots big((zeta(y)-x_n)) \
         &= big((y-y_0)) cdots big((y-y_n)), wide y_i = zeta^(-1)(x_i),
  $
  is a monic polynomial of degree~$n+1$ such that
  #numbered($
    sup_(y in [-1, 1]) abs(p(y)) = 2^(n+1)/(b-a)^(n+1) sup_(x in [a, b]) abs((x-x_0) dots (x - x_n)).
  $) <eq:change_of_variable>
  By~@theorem:minimum_infty_norm,
  the left-hand side is minimized when $p$ is equal to $2^(-n) T_(n+1)$,
  i.e.~when the roots of $p$ coincide with the roots of $T_(n+1)$.
  This occurs when
  $ y_i = zeta^(-1)(x_i) = cos lr((((2i + 1)pi)/(2 (n+1)))). $
  Applying the inverse change of variable $x_i = zeta(y_i)$,
  we deduce the result.
]

#indent @corollary:chebyshev_nodes is useful for interpolation.
The nodes
#numbered($
  x_i = a + (b-a) (1 + cos lr((((2i + 1)pi)/(2 (n+1)))))/2, wide i = 0, dots, n,
$)
are known as Chebyshev nodes and,
more often than not, employing these nodes for interpolation produces much better results than using equidistant nodes,
both in the case where~$u$ is a polynomial of degree $n+1$,
as we just proved,
but also for general~$u$.
As an example we plot in @fig:cheb_runge the interpolation of the Runge function using Chebyshev nodes.
In this case, the interpolating polynomial converges uniformly to the Runge function as we increase the number of interpolation nodes!
#figure(
  placement: top,
  grid(
    columns: (49%, 49%), column-gutter: 0.33em,
    image("/figures/interpolation_cheb_runge_10.pdf", width: 100%),
    image("/figures/interpolation_cheb_runge_15.pdf", width: 100%),
    image("/figures/interpolation_cheb_runge_20.pdf", width: 100%),
    image("/figures/interpolation_cheb_runge_30.pdf", width: 100%),
  ),
  caption: [Interpolation (in orange) of the Runge function~@eq:runge_function (in blue) using 10, 15, 20, and~30 Chebyshev nodes.],
) <fig:cheb_runge>

=== Hermite interpolation
Hermite interpolation,
sometimes also called Hermite--Birkoff interpolation,
generalizes Lagrange interpolation to the case where,
in addition to the function values $u_0, dots, u_n$,
the values of some of the derivatives are given at the interpolation nodes.
For simplicity,
we assume in this section that only the first derivative is specified.
In this case, the aim of Hermite interpolation is to find,
given data~$(x_i, u_i, u'_i)$ for $i in {0, dots, n}$,
a polynomial~$hat(u)$ of degree at most~$2n + 1$ such that
$ forall i in {0, dots, n}, wide hat(u)(x_i) = u_i, wide hat(u)'(x_i) = u'_i. $
In order to construct the interpolating polynomial,
it is useful to define the functions
$
  psi_i (x)
  = product_(j=0, j != i)^n lr(((x-x_j)/(x_i-x_j)))^2,
  wide i = 0, dots, n.
$
The function $psi_i$ is the square of the usual Lagrange polynomials associated with~$x_i$,
and it satisfies
$
  psi_i (x_i) = 1,
  wide psi'_i (x_i) = sum_(j=0,j != i)^n 2/(x_i - x_j),
  wide forall j != i wide psi_i (x_j) = psi'_i (x_j) = 0.
$
We consider the following ansatz for~$hat(u)$:
$ hat(u)(x) = sum_(i=0)^n psi_i (x) q_i (x), $
where $q_i$ are polynomials to be determined of degree at most one,
so that~$hat(u)$ is of degree at most~$2n+1$.
We then require
$ hat(u)(x_i) = q_i (x_i), wide hat(u)'(x_i) = psi'_i (x_i) q(x_i) + q'(x_i). $
From the first equation, we deduce that $q_i (x_i) = u_i$,
and from the second equation we then have~$q'(x_i) = hat(u)'(x_i) - psi'_i (x_i) u_i$.
We conclude that the interpolating polynomial is given by
$
  hat(u)(x) = sum_(i=0)^n psi_i (x) Big((u_i + big((u'_i - psi'_i (x_i) u_i)) (x - x_i))).
$
The following theorem gives an expression of the error.
#theorem(title: [Hermite interpolation error])[
  Assume that~$u colon [a, b] -> real$ is a function in $C^(2n+2)([a, b])$ and let~$hat(u)$
  denote the Hermite interpolation of~$u$ at the nodes~$x_0, dots, x_n$.
  Then for all $x in [a, b]$,
  there exists~$xi = xi(x)$ in the interval $[a, b]$ such that
  $
    u(x) - hat(u)(x) = (u^((2n+2))(xi))/((2n+2)!) (x-x_0)^2 dots (x - x_n)^2.
  $
] <theorem:hermite_interpolation>
#proof[
  See~@exercise:proof_hermite.
]

=== Piecewise interpolation <sub:piecewise_interpolation>
The interpolation methods we discussed so far are in some sense global;
they aim to construct one polynomial that goes through all the data points.
This approach is attractive because the interpolant is infinitely smooth but,
as we observed, it is not always fruitful,
in particular when equidistant interpolation nodes are employed.
An alternative approach is to divide the domain in a number of small intervals
and perform polynomial interpolation within each interval.
Although the resulting interpolating function is usually not smooth over the full domain,
this "local" approach to interpolation is more robust.

Several methods belong in the category of piecewise interpolation.
We mention, for instance, piecewise Lagrange interpolation and cubic splines interpolation.
In this section,
we briefly describe the former method,
which is widely used in the context of the _finite element method_.
Information on the latter method is available in~@MR2265914[Section 8.7.1.].

For simplicity, we illustrate the method in dimension 1,
but piecewise Lagrange interpolation can be extended to several dimensions.
Assume that we wish to approximate a function~$u colon [a, b] -> real$.
We consider a subdivision $a = x_0 < x_1 < dots < x_n = b$ of the interval $[a, b]$
and let~$h$ denote the maximum spacing:
$ h = max_(i in {0, dots, n-1}) abs(x_(i+1) - x_i). $
Within each subinterval $I_i = [x_i, x_(i+1)]$,
we consider a further subdivision
$ x_i = x_i^((0)) < x_i^((1)) < dots < x_i^((m)) = x_(i+1), $
where the nodes $x_i^((0)), dots, x_i^((m))$ are equally spaced with distance $h\/m$.
The idea of piecewise Lagrange interpolation is to calculate,
for each interval $I_i$ in the partition,
the interpolating polynomial~$p_i$ at the nodes $x_i^((0)), dots, x_i^((m))$.
The interpolant is then defined as
#numbered($
  hat(u)(x) = p_iota (x),
$) <eq:piecewise_interpolation>
where $iota = iota(x)$ is the index of the interval to which $x$ belongs.
Since $x_i^((m)) = x_(i+1) = x_(i+1)^((0))$,
the interpolant defined by~@eq:piecewise_interpolation is continuous.
If the function $u$ belongs to $C^(m+1)([a, b])$,
then by~@corollary:interpolation_error the interpolation error within each subinterval may be bounded from above as follows:
#numbered($
  sup_(x in I_i) abs(u(x) - hat(u)(x)) <= (C_(m+1) (h\/m)^(m+1))/(4(m + 1)),
  wide C_(m+1) := sup_(x in [a, b]) abs(u^((m+1))(x)),
$) <eq:error_estimate_piecewise_interpolation>
and so we deduce
$ sup_(x in [a, b]) abs(u(x) - hat(u)(x)) <= C h^(m+1), $
for an appropriate constant~$C$ independent of $h$.
This equation shows that the error is guaranteed to decrease to~0 in the limit as $h -> 0$.
In practice, the number $m$ of interpolation nodes within each interval can be small.
// Several integration methods presented in @cha:quadrature are related to piecewise interpolation.

== Approximation
In this section,
we focus on the subject of _approximation_,
both of discrete data points and of continuous functions.
We begin, in @sub:polynomial_regression with a discussion of least squares approximation for data points,
and in @sub:mean_square_approximation we focus on function approximation in the mean square sense.

=== Least squares approximation of data points <sub:polynomial_regression>
Consider $n + 1$ distinct $x$ values $x_0 < dots < x_n$,
and suppose that we know the values $u_0, dots, u_n$ taken by an unknown function $u$ when evaluated at these points.
We wish to approximate the function~$u$ by a function of the form
#numbered($
  hat(u)(x) = sum_(i=0)^m alpha_i phi_i (x) in Span { phi_0, dots, phi_m },
$) <eq:expansion_approx>
for some $m < n$.
In many cases of practical interest,
the basis functions $phi_0, dots, phi_m$ are polynomials.
In contrast with interpolation,
here we seek a function $hat(u)$ in a finite-dimensional function space of dimension~$m$ strictly lower than the number of data points.
In order for~$hat(u)$ to be a good approximation,
we wish to find coefficients $alpha_0, dots, alpha_m$ such that
the following linear system is approximately satisfied.
$
  mat(A) vect(alpha) :=
  pmat(
    phi_0(x_0), phi_1(x_0), …, phi_m (x_0);
    phi_0(x_1), phi_1(x_1), …, phi_m (x_1);
    phi_0(x_2), phi_1(x_2), …, phi_m (x_2);
    ⋮, ⋮, , ⋮;
    phi_0(x_(n-2)), phi_1(x_(n-2)), …, phi_m (x_(n-2));
    phi_0(x_(n-1)), phi_1(x_(n-1)), …, phi_m (x_(n-1));
    phi_0(x_n), phi_1(x_n), …, phi_m (x_n)
  )
  pmat(alpha_0; alpha_1; ⋮; alpha_m)
  approx
  pmat(u_0; u_1; u_2; ⋮; u_(n-2); u_(n-1); u_n) =: vect(b).
$
In general,
since the matrix on the left-hand side has more lines than columns,
there does not exist an exact solution to this equation.
In order to find an approximate solution,
a natural approach is to find coefficients $alpha_0, dots, alpha_m$ such that
the residual vector $vect(r) = mat(A) vect(alpha) - vect(b)$ is small in some vector norm.
A particularly popular approach,
known as least squares approximation,
is to minimize the Euclidean norm of~$vect(r)$ or,
equivalently, the square of the Euclidean norm:
$
  norm(vect(r))^2 = sum_(i=0)^n big(|u_i - hat(u)(x_i)|)^2 = sum_(i=0)^n lr((u_i - sum_(j=0)^m alpha_j phi_j (x_i)))^2.
$
Let us denote the right-hand side of this equation by~$J(vect(alpha))$,
which we view as a function of the vector of coefficients $vect(alpha)$.
A necessary condition for~$vect(a)_*$ to be a minimizer is that $nabla J(vect(alpha)_*) = 0$.
The gradient of $J$,
written as a column vector, is given by
$
  nabla J(vect(alpha))
  &= nabla Big(((mat(A) vect(alpha) - vect(b))^T (mat(A) vect(alpha) - vect(b)))) \
  &= nabla Big((vect(alpha)^T (mat(A)^T mat(A)) vect(alpha) - vect(b) mat(A) vect(alpha) - vect(alpha)^T mat(A)^T vect(b) + vect(b)^T vect(b)))) \
  &= 2 (mat(A)^T mat(A)) vect(alpha) - 2 mat(A)^T vect(b).
$
We deduce that $vect(alpha)_*$ solves the linear system
#numbered($
  mat(A)^T mat(A) vect(alpha)_* = mat(A)^T vect(b),
$) <eq:normal_equation>
where the matrix on the left-hand side is given by:
$
  mat(A)^T mat(A) :=
  pmat(
    sum_(i=0)^n phi_0(x_i) phi_0(x_i), sum_(i=0)^n phi_0(x_i) phi_1(x_i), …, sum_(i=0)^n phi_0(x_i) phi_m (x_i);
    sum_(i=0)^n phi_1(x_i) phi_0(x_i), sum_(i=0)^n phi_1(x_i) phi_1(x_i), …, sum_(i=0)^n phi_1(x_i) phi_m (x_i);
    ⋮, ⋮, , ⋮;
    sum_(i=0)^n phi_m (x_i) phi_0(x_i), sum_(i=0)^n phi_m (x_i) phi_1(x_i), …, sum_(i=0)^n phi_m (x_i) phi_m (x_i)
  ).
$
Equation~@eq:normal_equation is a system of~$m$ equations with $m$ unknowns,
which admits a unique solution provided that $mat(A)^T mat(A)$ is full rank or,
equivalently, the columns of~$mat(A)$ are linearly independent.
The linear equations~@eq:normal_equation are known as the _normal equations_.
As a side note,
we mention that
the solution $vect(alpha)_* = (mat(A)^T mat(A))^(-1) mat(A)^T vect(b)$
coincides with the maximum likelihood estimator for $alpha$ under the assumption that the data is generated according to
$u_i = u(x_i) + epsilon_i$, for some function $u in Span {phi_0, dots, phi_m}$
and random noise $epsilon_i tilde.op cal(N)(0, 1)$.

#remark[
  From equation~@eq:normal_equation we deduce that
  $ mat(A) vect(alpha)_* = mat(A) (mat(A)^T mat(A))^(-1) mat(A)^T vect(b). $
  The matrix~$Pi_(mat(A)) := mat(A) (mat(A)^T mat(A))^(-1) mat(A)^T$ on the right-hand side is the orthogonal projection operator onto $col(mat(A)) subset real^n$,
  the subspace spanned by the columns of~$mat(A)$.
  Indeed, it holds that $Pi_(mat(A))^2 = Pi_(mat(A))$,
  which is the defining property of a projection operator.
]

#indent To conclude this section,
we note that the matrix $mat(A)^+ = (mat(A)^T mat(A))^(-1) mat(A)^T$ is a left inverse of the matrix $mat(A)$,
because $mat(A)^+ mat(A) = mat(I)$.
It is also called the Moore--Penrose inverse or pseudoinverse of the matrix~$mat(A)$,
which generalizes the usual inverse matrix.
In Julia, the backslash operator silently uses the Moore--Penrose inverse when employed with a rectangular matrix.
Therefore, solving the normal equations~@eq:normal_equation can be achieved by just writing `α = A\b`.

=== Mean square approximation of functions <sub:mean_square_approximation>
The approach described in~@sub:polynomial_regression can be generalized to the setting where the actual function~$u$,
rather than just discrete evaluations of it, is available.
In this section, we seek an approximation of the form~@eq:expansion_approx
such that the error $hat(u)(x) - u(x)$, measured in some function norm,
is minimized.
Of course, the solution to this minimization problem depends in general on the norm employed,
and may in some cases not even be unique.
Instead of specifying a particular norm,
as done in @sub:polynomial_regression,
in this section we retain some generality and assume only that the norm is induced by an inner product on the space of real-valued continuous functions:
#numbered($
  ip(placeholder, placeholder) colon C([a,b]) times C([a, b]) -> real.
$) <eq:inner_product_approx>
In other words,
we seek to minimize
$ J(vect(alpha)) := norm(hat(u) - u)^2 = ip(hat(u) - u, hat(u) - u). $
This is again a function of the $m+1$ variables $alpha_0, dots, alpha_m$.
Before calculating its gradient,
we rewrite the function $J(vect(alpha))$ in a simpler form:
$
  J(vect(alpha))
  &= ipstar(u - sum_(j=0)^m alpha_j phi_j, u - sum_(k=0)^m alpha_k phi_k) \
  &= sum_(j=0)^m sum_(k=0)^m alpha_j alpha_k ip(phi_j, phi_k) - 2 sum_(j=0)^m alpha_j ip(u, phi_j) + ip(u, u)
  = vect(alpha)^T mat(G) vect(alpha) - 2 vect(b)^T vect(alpha) + ip(u, u),
$
where we introduced
#numbered($
  mat(G) :=
  pmat(
    ip(phi_0, phi_0), ip(phi_0, phi_1), …, ip(phi_0, phi_m);
    ip(phi_1, phi_0), ip(phi_1, phi_1), …, ip(phi_1, phi_m);
    ⋮, ⋮, , ⋮;
    ip(phi_m, phi_0), ip(phi_m, phi_1), …, ip(phi_m, phi_m)
  ),
  wide
  vect(b) :=
  pmat(ip(u, phi_0); ip(u, phi_1); ⋮; ip(u, phi_m)).
$) <eq:approx_matrix_m>
Employing the same approach as in the previous section,
we then obtain $nabla J(vect(alpha)) = mat(G) vect(alpha) - vect(b)$,
and so the minimizer of $J(vect(alpha))$ is the solution to the equation
#numbered($
  mat(G) vect(alpha) = vect(b).
$) <eq:normaal_continuous>
The matrix $mat(G)$,
known as the _Gram matrix_,
is positive semi-definite and nonsingular provided that the basis functions are linearly independent,
see~@exercise:matrix_approx.
Therefore, the solution $vect(alpha)_*$ exists and is unique.
In addition, since the Hessian of~$J$ is equal to~$mat(G)$,
the vector~$vect(alpha)_*$ is indeed a minimizer.
Note that if~$ip(placeholder, placeholder)$ is defined as a finite sum of the form
#numbered($
  ip(f, g) = sum_(i=0)^n f(x_i) g(x_i),
$) <eq:inner_finite>
then~@eq:normaal_continuous coincides with the normal equations~@eq:normal_equation from the previous section.
We remark that~@eq:inner_finite is in fact not an inner product on the space of continuous functions,
but it is an inner product on the space of polynomials of degree less than or equal to~$n$.

In practice, the matrix and right-hand side of the linear system~@eq:normaal_continuous can usually not be calculated exactly,
because the inner product~$ip(placeholder, placeholder)$ is defined through an integral; see~@eq:ip_weighted in the next section.

#remark[
  Rewriting the normal equations~@eq:normaal_continuous in terms of~$hat(u)$
  we obtain
  $
    ip(hat(u) - u, phi_0) = 0,
    wide
    dots, wide
    ip(hat(u) - u, phi_m) = 0.
  $
  Therefore, the optimal approximation $hat(u) in Span { phi_0, dots, phi_m }$ satisfies
  $ forall v in Span { phi_0, dots, phi_m }, wide ip(hat(u) - u, v) = 0. $
  This shows that the optimal approximation $hat(u)$,
  in the sense of the norm~$norm(placeholder)$,
  is the orthogonal projection of~$u$ onto $Span { phi_0, dots, phi_m }$.
]
=== Orthogonal polynomials
The Gram matrix~$mat(G)$ in~@eq:normaal_continuous is equal to the identity matrix when the basis functions are orthonormal for the inner product considered.
In this case, the solution to the linear system is
$ alpha_i = ip(u, phi_i), wide i = 0, dots, m, $
and so the best approximation~$hat(u)$ (for the norm induced by the inner product considered!) is simply given by
$ hat(u) = sum_(i=0)^m ip(u, phi_i) phi_i. $
The coefficients $ip(u, phi_i)$ of the basis functions in this expansion are called _Fourier coefficients_.
Given a finite dimensional subspace~$cal(S)$ of the space of continuous functions,
an orthonormal basis can be constructed via the Gram--Schmidt process.
In this section, we focus on the particular case where $cal(S) = poly(n)$
-- the subspace of polynomials of degree less than or equal to $n$.
Another widely used approach,
which we do not explore in this course,
is to use trigonometric basis functions.
We also assume that the inner product~@eq:inner_product_approx is of the form
#numbered($
  ip(f, g)
  = integral_a^b f(x) g(x) w(x) thin dd x,
$) <eq:ip_weighted>
where $w(x)$ is a given nonnegative weight function such that
$ integral_a^b w(x) thin dd x > 0. $
Let $phi_0(x), phi_1(x), phi_2(x) dots$ denote the orthonormal polynomials obtained by applying the Gram--Schmidt procedure to the monomials $1, x, x^2, dots$.
These depend in general on the weight~$w(x)$ and on the interval~$[a, b]$.
A few of the popular classes of orthogonal polynomials are presented in the table below:
#align(center, tabular(
  columns: 3, stretch: 1.5,
  [*Name*], [$w(x)$], [$[a, b]$],
  [Legendre], [$1$], [$[-1, 1]$],
  [Chebyshev], [$1/sqrt(1 - x^2)$], [$(-1, 1)$],
  [Hermite], [$exp lr((- x^2/2))$], [$[-infinity, infinity]$],
  [Laguerre], [$ee^(-x)$], [$[0, infinity]$],
))

Orthogonal polynomials have a rich structure,
and in the rest of this section we prove some of their key properties,
one of which will be very useful in the context of numerical integration in~@cha:quadrature.
We begin by showing that
orthogonal polynomials have distinct real roots.
#proposition[
  Assume for simplicity that $w(x) > 0$ for all $x in [a, b]$,
  and let $phi_0, phi_1, dots$ denote orthonormal polynomials of increasing degree for the inner product~@eq:ip_weighted.
  Then for all $n in nat$,
  the polynomial $phi_n$ has~$n$ distinct roots in the open interval $(a, b)$.
] <proposition:distinct_roots>
#proof[
  Reasoning by contradiction,
  we assume that~$phi_n$ changes sign at only $k < n$ points of the open interval~$(a, b)$,
  which we denote by $x_1, dots, x_k$.
  Then
  $ phi_n (x) times (x - x_0) (x - x_1) dots (x - x_k) $
  is either everywhere nonnegative or everywhere nonpositive over $[a, b]$.
  But then
  $ integral_a^b phi_n (x) times (x - x_1) dots (x - x_k) thin w(x) thin dd x $
  is nonzero,
  which is a contradiction because the product $(x-x_1) dots (x-x_k)$ is a polynomial of degree~$k$,
  which is orthogonal to $phi_n$ by assumption.
  Indeed, being orthogonal to~$phi_0, dots, phi_(n-1)$,
  the polynomial $phi_n$ is also orthogonal to any linear combination of these polynomials.
]

#indent Next, we show that orthogonal polynomials satisfy a three-term recurrence relation.

#proposition[
  Assume that $phi_0, phi_1, dots$ are orthonormal polynomials
  for some inner product of the form~@eq:ip_weighted such that $phi_i$ is of degree~$i$.
  Then
  #numbered($
    forall n in {1, 2, dots}, wide
    alpha_(n+1) phi_(n+1) (x) = (x - beta_n) phi_n (x) - alpha_n phi_(n-1)(x),
  $) <eq:recurrence_ortho>
  where
  $
    alpha_n = ipbig(x phi_n, phi_(n-1)),
    wide beta_n = ipbig(x phi_n, phi_n).
  $
  In addition, $alpha_1 phi_1 (x) = (x - beta_0) phi_0(x)$.
] <proposition:three_term_recurrence>
#proof[
  Since $x phi_n (x)$ is a polynomial of degree $n+1$,
  it may be decomposed as
  #numbered($
    x phi_n (x) = sum_(i=0)^(n+1) gamma_(n,i) phi_i (x).
  $) <eq:expansion_times_x>
  Taking the inner product of both sides of this equation with~$phi_i$
  and employing the orthonormality assumption,
  we obtain an expression for the coefficients:
  $ gamma_(n,i) = ipbig(x phi_n, phi_i). $
  From the expression~@eq:ip_weighted of the inner product,
  it is clear that $ipbig(x phi_n, phi_i) = ipbig(phi_n, x phi_i)$.
  Since~$x phi_i$ is a polynomial of degree $i+1$ and $phi_n$ is orthogonal to all polynomials of degree strictly less than~$n$,
  we deduce that $gamma_(n,i) = 0$ if~$i < n-1$.
  Consequently, we can rewrite the right-hand side of~@eq:expansion_times_x as a sum involving only three terms
  #numbered($
    x phi_n (x) = ipbig(x phi_n, phi_(n-1)) phi_(n-1)(x) + ipbig(x phi_n, phi_n) phi_n (x) + ipbig(x phi_n, phi_(n+1)) phi_(n+1)(x).
  $) <eq:expansion_times_x_good>
  Since $ipbig(x phi_n, phi_(n+1)) = ipbig(x phi_(n+1), phi_n)$,
  we obtain the statement after rearranging.
]
#remark[
  Notice that the polynomials in~@proposition:three_term_recurrence are orthonormal by assumption,
  and so the coefficient~$alpha_(n+1)$ is just a normalization constant.
  We deduce that
  $
    phi_(n+1) (x) = ((x - beta_n) phi_n (x) - alpha_n phi_(n-1)(x))/norm((x - beta_n) phi_n (x) - alpha_n phi_(n-1)(x)),
  $
  which enables to calculate the orthogonal polynomials recursively.
]

=== Orthogonal polynomials and numerical integration: an introduction~#moreinfo <sub:orthogonal_integration>
Equation~@eq:expansion_times_x_good may be rewritten in matrix form as follows:
$
  pmat(
    x phi_0(x);
    x phi_1(x);
    x phi_2(x);
    ⋮;
    x phi_(m-1)(x);
    x phi_m (x)
  )
  =
  pmat(
    beta_0, alpha_1, , , , ;
    alpha_1, beta_1, alpha_2, , , ;
    , alpha_2, beta_2, alpha_3, , ;
    , , ⋱, ⋱, ⋱, ;
    , , , alpha_(m-1), beta_(m-1), alpha_m;
    , , , , alpha_m, beta_m
  )
  pmat(
    phi_0(x);
    phi_1(x);
    phi_2(x);
    ⋮;
    phi_(m-1)(x);
    phi_m (x)
  )
  +
  pmat(
    0;
    0;
    0;
    ⋮;
    0;
    alpha_(m+1) phi_(m+1)(x)
  ).
$
Let $mat(T)$ denote the matrix on the left-hand side of this equation,
and let $r_0, dots, r_m$ denote the roots of~$phi_(m+1)$.
By~@proposition:distinct_roots,
these are distinct and all belong to the interval~$(a, b)$.
The second term on the right-hand side cancels out when~$x$ is a root of $phi_(m+1)$,
and so
$
  forall r in Big(\{ r_0\, dots\, r_m \}), wide
  pmat(
    r phi_0(r);
    r phi_1(r);
    ⋮;
    r phi_(m-1)(r);
    r phi_m (r)
  )
  =
  pmat(
    beta_0, alpha_1, , , ;
    alpha_1, beta_1, alpha_2, , ;
    , ⋱, ⋱, ⋱, ;
    , , alpha_(m-1), beta_(m-1), alpha_m;
    , , , alpha_m, beta_m
  )
  pmat(
    phi_0(r);
    phi_1(r);
    ⋮;
    phi_(m-1)(r);
    phi_m (r)
  ).
$
In other words,
for any root $r$ of~$phi_(m+1)$,
the vector $pmat(phi_0(r), …, phi_m (r))^T$
is an eigenvector of the matrix~$mat(T)$,
with associated eigenvalue equal to~$r$.
Since~$mat(T)$ is a symmetric matrix,
the eigenvectors associated with distinct eigenvalues are orthogonal for the Euclidean inner product of~$real^(m+1)$,
so given that the eigenvalues of~$mat(T)$ are distinct,
we deduce that
#numbered($
  forall i != j, wide
  sum_(i=0)^m phi_i (r_i) phi_i (r_j) = 0.
$) <eq:orthogonality>
Let us construct the matrix
$
  mat(P) =
  pmat(
    phi_0(r_0), phi_1(r_0), …, phi_m (r_0);
    phi_0(r_1), phi_1(r_1), …, phi_m (r_1);
    phi_0(r_2), phi_1(r_2), …, phi_m (r_2);
    ⋮, ⋮, …, ⋮;
    phi_0(r_m), phi_1(r_m), …, phi_m (r_m)
  ).
$
Equation~@eq:orthogonality indicates that the rows of~$mat(P)$ are orthogonal,
and so the matrix~$mat(D) = mat(P) mat(P)^T$ is diagonal with elements given by
$
  d_(i i) = sum_(j=0)^m abs(phi_j (r_i))^2,
  wide i = 0, dots, m.
$
(Here we start counting the rows from 0 for convenience.)
Since $mat(P) mat(P)^T mat(D)^(-1) = mat(I)$,
we deduce that the inverse of~$mat(P)$ is given by $mat(P)^(-1) = mat(P)^T mat(D)^(-1)$.
Consequently,
$ mat(P)^T mat(D)^(-1) mat(P) = mat(P)^(-1) mat(P) = mat(I), $
which means that the _columns_ of $mat(P)$ are orthonormal for the inner product
$(vect(x), vect(y)) |-> vect(x)^T mat(D)^(-1) vect(y)$.
In other words,
the polynomials~$phi_1, dots, phi_m$ are orthonormal for the inner product
$
  ip(placeholder, placeholder)_(m+1) colon & poly(m) times poly(m) -> real; \
                                           & (p, q) |-> vect(sum)_(i=0)^m (p(r_i) q(r_i))/d_(i i).
$
We have thus shown that,
if~$phi_0, phi_1, phi_2, dots$ is a family of orthonormal polynomials for an inner product~$ip(placeholder, placeholder)$,
then these are also orthonormal for the inner product~$ip(placeholder, placeholder)_(m+1)$.
We reformulate our findings in the following result where,
since~$m$ was arbitrary in the previous reasoning,
we add a superscript to indicate when the quantities involved depend on~$m$.
#theorem[
  Orthonormal polynomials~$phi_0, dots, phi_m$ for the inner product
  $ ip(f, g) = integral_a^b f(x) g(x) w(x) thin dd x $
  are also orthonormal for the inner product
  $
    ip(f, g)_(m+1) = sum_(i=0)^m f lr((r^((m+1))_i)) g lr((r^((m+1))_i)) w^((m+1))_i,
  $
  where $r^((m+1))_0, dots, r^((m+1))_m$ are the roots of~$phi_(m+1)$ and
  the weights~$w^((m+1))_i$ are given by
  $
    w^((m+1))_i = 1/(sum_(j=0)^m lr(|phi_j lr((r^((m+1))_i))|)^2), wide i = 0, dots, m.
  $
] <theorem:orthogonal_polynomials>
As an immediate corollary, we deduce that
#numbered($
  forall (p, q) in poly(m) times poly(m), wide
  ip(p, q) = ip(p, q)_(m+1),
$) <eq:orthogonality_equivalence>
Indeed, denoting by $p = alpha_0 phi_0 + cdots + alpha_m phi_m$ and $q = beta_0 phi_0 + cdots + beta_m phi_m$ the expansions of the polynomials~$p$ and~$q$
in the orthonormal basis,
we have
$
  ip(p, q)
  &= ip(alpha_0 phi_0 + cdots + alpha_m phi_m, beta_0 phi_0 + cdots + beta_m phi_m) \
  &= sum_(i=0)^m sum_(j=0)^m alpha_i beta_j ip(phi_i, phi_j)
  = alpha_0 beta_0 + cdots + alpha_m beta_m
  = sum_(i=0)^m sum_(j=0)^m alpha_i beta_j ip(phi_i, phi_j)_(m+1) \
  & = ip(alpha_0 phi_0 + cdots + alpha_m phi_m, beta_0 phi_0 + cdots + beta_m phi_m)_(m+1) = ip(p, q)_(m+1).
$

To conclude this section,
we prove the following statement,
which is another consequence of~@theorem:orthogonal_polynomials and has applications to numerical integration.
#theorem[
  It holds that
  #numbered($
    forall p in poly(2m+1),
    wide integral_a^b p(x) thin w(x) thin dd x
    = sum_(i=0)^m p lr((r_i^((m+1)))) w^((m+1))_i.
  $) <eq:integration_formula>
]
#proof[
  Taking~$q = 1$ in~@eq:orthogonality_equivalence and employing the definitions of~$ip(placeholder, placeholder)$ and $ip(placeholder, placeholder)_(m+1)$,
  we have that~@eq:orthogonality_equivalence is satisfied for any $p in poly(m)$.
  Next, any polynomial $p in poly(2m+1)$ may be decomposed as $p(x) = phi_(m+1)(x) q(x) + rho(x)$,
  for some polynomial~$q$ of degree~$m$ (the quotient of the polynomial division of $p$ by $phi_(m+1)$) and some polynomial~$rho$ of degree lower than or equal to~$m$ (the remainder of the polynomial division).
  Therefore, since~@eq:integration_formula was already shown to hold for polynomials of degree up to~$m$,
  we obtain
  $
    integral_a^b p(x) w(x) thin dd x
    &= integral_a^b phi_(m+1)(x) q(x) w(x) thin dd x + integral_a^b rho(x) w(x) thin dd x \
    &= 0 + integral_a^b rho(x) w(x) thin dd x
    = 0 + sum_(i=0)^m rho lr((r_i)) w_i \
    &= sum_(i=0)^m phi_(m+1) lr((r_i)) q lr((r_i)) w_i + sum_(i=0)^m rho lr((r_i)) w_i = sum_(i=0)^m p lr((r_i)) w_i,
  $
  where we dropped the $""^((m+1))$ superscript for conciseness and
  we used, in the penultimate inequality,
  the fact that~$r_0, dots, r_m$ are the roots of the polynomial~$phi_(m+1)$.
]
Since the left-hand side of~@eq:integration_formula is an integral and the right-hand side is a sum,
we have just constructed an integration formula,
which enjoys a very nice property: it is exact for polynomials of degree up to~$2m + 1$!
A formula of this type is called a _quadrature formula_,
with~$m+1$ _nodes_~$r^((m+1))_0, dots, r^((m+1))_m$ and associated _weights_ $w^((m+1))_0, dots, w^((m+1))_m$.
Note that the nodes and weights of the quadrature depend on the weight~$w(x)$ and on the degree~$m$.
We will revisit this subject in~@cha:quadrature.

== Exercises
#exercise[
  Find the polynomial $p(x) = a x + b$ (a straight line) that goes through the points $(x_0, u_0)$ and $(x_1, u_1)$.
]

#exercise[
  Find the polynomial $p(x) = a x^2 + b x + c$ (a parabola) that goes through the points $(0, 1)$, $(1, 3)$ and $(2, 7)$.
]

#exercise[
  Prove the following recurrence relation for Chebyshev polynomials:
  $ T_(i+1)(x) = 2 x T_i (x) - T_(i-1)(x), wide i = 1, 2, dots. $
]

#exercise[
  Show by recursion that
  #numbered($
    [u_0, u_1, dots, u_n] = sum_(j=0)^n u_j/(product_(k in {0, dots, n} backslash {j}) (x_j - x_k)).
  $) <eq:divided_differences>
  Deduce from this identity that
  $ [u_0, u_1, dots, u_n] = [u_(sigma_1), u_(sigma_2), dots, u_(sigma_n)], $
  for any permutation~$vect(sigma)$ of~$(0, 1, 2, dots, n)$.
] <exercise:divided_differences>
#solution[
  The first statement~@eq:divided_differences is clear when $n = 0$.
  Reasoning by induction, we assume that the statement is true up to $n-1$ and prove that it then also holds for $n$.
  Using the definition~@eq:definition_divided_difference and the induction hypothesis,
  we obtain that
  $
    [u_0, u_1, dots, u_n]
    &= ([u_1, dots, u_n] - [u_0, dots, u_(n-1)])/(x_n - x_0) \
    &= 1/(x_n - x_0) lr((sum_(j=bold(1))^n u_j/(product_(k in {bold(1), dots, n} backslash {j}) (x_j - x_k)) - sum_(j=0)^(bold(upright(n-1))) u_j/(product_(k in {0, dots, bold(upright(n-1))} backslash {j}) (x_j - x_k))))
  $
  Rewriting the fractions with a common denominator leads to
  $
    [u_0, u_1, dots, u_n]
    = 1/(x_n - x_0) sum_(j=0)^n (u_j big(((x_j - x_0) - (x_j - x_n))))/(product_(k in {0, dots, n} backslash {j}) (x_j - x_k))
    = sum_(j=0)^n u_j/(product_(k in {0, dots, n} backslash {j}) (x_j - x_k)),
  $
  which concludes the proof of the first statement.
  The second statement then follows immediately,
  because the right-hand side of~@eq:divided_differences is invariant under permutations.
]

#exercise[
  Using the Gregory--Newton formula,
  find an expression for
  $ sum_(i=1)^n i^4. $
]

#exercise[
  Let $(f_0, f_1, f_2, dots) = (1, 1, 2, dots)$ denote the Fibonacci sequence.
  Prove that there does not exist a polynomial~$p$ such that
  #numbered($
    forall n in nat, wide
    f_n = p(n).
  $) <eq:fibonacci_polynomial>
]
#solution[
  Assume by contradiction that~$p colon real -> real$ is a polynomial such that~@eq:fibonacci_polynomial is satisfied,
  and let~$n$ be the degree of this polynomial.
  Then it holds that $Delta^(n+1) p = 0$~@eq:falling_powers_rule,
  where both sides are viewed as functions from $real$ to $real$.
  On the other hand,
  since $p(n) = f_n$ for all~$n in nat$,
  we can calculate explicitly the values of taken by the function $Delta^m p$ when evaluated at all the natural numbers,
  for all $m in nat$.
  We collate a few values in the following table.
  #align(center, tabular(
    columns: 8,
    [$n$], [0], [1], [2], [3], [4], [5], [6],
    [$Delta^0 p(n)$], [$bold(1)$], [$bold(1)$], [$bold(2)$], [$bold(3)$], [$bold(5)$], [$bold(8)$], [$bold(13)$],
    [$Delta^1 p(n)$], [$0$], [$bold(1)$], [$bold(1)$], [$bold(2)$], [$bold(3)$], [$bold(5)$], [$bold(8)$],
    [$Delta^2 p(n)$], [$1$], [$0$], [$bold(1)$], [$bold(1)$], [$bold(2)$], [$bold(3)$], [$bold(5)$],
    [$Delta^3 p(n)$], [$-1$], [$1$], [$0$], [$bold(1)$], [$bold(1)$], [$bold(2)$], [$bold(3)$],
  ))
  It appears from these calculations that the Fibonacci sequence is shifted one position to the right with each additional application of $Delta$.
  In other words, our calculations suggest that
  #numbered($
    forall (m, n) in nat times nat,
    wide Delta^m p(m+n) = f_n,
  $) <eq:difference_fib>
  which is a contradiction.
  To conclude, let us prove~@eq:difference_fib rigorously.
  This equation is obvious for $m = 0$ by assumption.
  Now, reasoning by contradiction, we assume that~@eq:difference_fib is true up to~$m$.
  Then by definition of the difference operator~$Delta$,
  we have
  $
    Delta^(m+1) p(m+n+1)
    &= Delta^m p(m+n+2) - Delta^m p(m+n+1) \
    &= f_(n+2) - f_(n+1) = f_n.
  $
  Here, we used the induction hypothesis~@eq:difference_fib in the second equality,
  and the definition of the Fibonacci series in the third one.
]

#exercise[
  Using the Gregory--Newton formula,
  show that
  #numbered($
    forall n in nat,
    wide 2^n = 1 + n + n^(underline(2))/(2!) + n^(underline(3))/(3!) + n^(underline(4))/(4!) + cdots
  $) <eq:discrete_exponential>
  // The solution and the remark are inside the exercise in the LaTeX source
  #set text(style: "normal")
  #solution[
    Equation~@eq:discrete_exponential is a particular case of the following more general statement:
    for any function $f in real -> real$,
    it holds that
    #numbered($
      forall n in nat,
      wide f(n) = f(0) + Delta f(0) n + Delta^2 f(0) n^(underline(2))/(2!) + Delta^3 f(0) n^(underline(3))/(3!) + Delta^4 f(0) n^(underline(4))/(4!) + cdots
    $) <eq:more_general_agreement_interpolation>
    In order to show this equation,
    it is sufficient to prove that for any $n_* in nat$,
    the two sides of~@eq:more_general_agreement_interpolation coincide for every $n in {0, dots, n_*}$.
    Since $n^(underline(p)) = 0$ for all $n in {0, dots, p-1}$ by definition~@eq:falling_powers of the falling powers,
    the right-hand side of~@eq:more_general_agreement_interpolation coincides for all $n in {0, dots, n_*}$ with
    $
      g(n) = f(0) + Delta f(0) n + Delta^2 f(0) thin n^(underline(2))/(2!) + cdots + Delta^(n_*) f(0) thin n^(underline(n_*))/(n_*!).
    $
    We recognize on the right-hand side Newton's expression of the interpolating polynomial through the points $big((0, f(0))), dots, big((n_*, f(n_*)))$,
    and so $g(n) = f(n)$ for all $n in {0, dots, n_*}$,
    which concludes the proof.
  ]
  #remark[
    Remarkably, equation~@eq:discrete_exponential holds in fact for _any_ $n in real_(> 0)$.
    However, showing this more general statement is beyond the scope of this course.
  ]
]

#exercise[
  Prove~@theorem:hermite_interpolation.
] <exercise:proof_hermite>

#exercise[
  Show that the matrix $mat(G)$ in @eq:approx_matrix_m is positive definite
  if the basis functions~$phi_0, dots, phi_m$ are linearly independent.
] <exercise:matrix_approx>

#compexercise[
  Write a Julia code for interpolating the following function using a polynomial of degree 20 over the interval $[-1, 1]$.
  $
    f(x) = tanh lr(((x+1\/2)/epsilon)) + tanh lr((x/epsilon)) + tanh lr(((x-1\/2)/epsilon)),
    wide epsilon = .01.
  $
  Use equidistant and then Chebyshev nodes,
  and compare the two approaches in terms of accuracy.
  Plot the function~$f$ together with the approximating polynomials.
]

#compexercise[
  Write from scratch a function to obtain the polynomial interpolating the data points
  $ (x_0, u_0), dots, (x_n, u_n). $
  Your function should return the values taken by the interpolating polynomial
  when evaluated at the points $X_0, dots, X_m$.
  You may use the following code to test your function
  ```julia
  import Plots
  function interp(X, x, u)
      # Your code comes here
  end
  n, m = 10, 100
  f(t) = cos(2π * t)
  x = LinRange(0, 1, n)
  X = LinRange(0, 1, m)
  u = f.(x)
  U = interp(X, x, u)
  Plots.plot(X, f.(X), label="Original function")
  Plots.plot!(X, U, label="Interpolation")
  Plots.scatter!(x, u, label="Data")
  ```
]

#compexercise[
  We wish to use interpolation to approximate the following parametric function,
  called an _epitrochoid_:
  #numbered($
    x (theta) &= (R + r)cos theta + d cos lr(((R + r)/r theta))
  $)
  #v(-16.3pt + 4.5pt)
  #numbered($
    y (theta) &= (R + r)sin theta - d sin lr(((R + r)/r theta)),
  $)
  with $R = 5$, $r = 2$ and $d = 3$, and for $theta in [0, 4pi]$.
  Write a Julia program to interpolate $x(theta)$ and $y(theta)$ using 40 equidistant points.
  Use the `BigFloat` format in order to reduce the impact of round-off errors.
  After constructing the polynomial interpolations $hat(x)(theta)$ and $hat(y)(theta)$,
  plot the parametric curve $theta |-> big((hat(x)(theta), hat(y)(theta)))$.
  Your plot should look similar to @fig:parametric_interpolation.
] <exercise:parametric_interpolation>
// In the LaTeX source, this figure is inside the exercise
#figure(
  placement: top,
  image("/figures/interpolation.pdf", width: 99%),
  caption: [Solution for @exercise:parametric_interpolation.],
) <fig:parametric_interpolation>

#compexercise(title: [Solving the Laplace equation using a spectral method])[
  The classical Laplace equation with homogeneous Dirichlet boundary conditions in dimension 1 reads
  #numbered($
    "Find" u in C^2 big(([0, 1])) "such that"
    wide
    cases(
      - u''(x) = f(x) wide forall x in (0, 1)\,,
      u(0) = u(1) = 0.
    )
  $) <eq:laplace>
  Our goal in this exercise is to approximate the exact solution~$u(x)$ using interpolation.
  Specifically,
  we propose to proceed in two steps:
  - Interpolate the right-hand side using a polynomial with equidistant nodes.
    That is, find a polynomial $hat(f) in poly(n)$ such that
    $
      forall i in {0, dots, n}, wide
      hat(f)(x_i) = f(x_i), wide x_i = i/n.
    $

  - Solve~@eq:laplace with $hat(f)$ instead of $f$.
    Since $hat(f)$ is a polynomial, this can be achieved analytically.
  Implement this program in the case where
  $
    f(x) = exp big((sin(2pi x))) cos(2 pi x)^2 - exp big((sin(2 pi x))) sin(2 pi x),
  $
  and compare for various values of~$n$ the approximate solution you obtain with
  the exact solution to~@eq:laplace,
  which is given by $u(x) = (2pi)^(-2) exp big((big((sin(2pi x))) - 1))$ in this case.
]

#compexercise(title: [Modeling the vapor pressure of mercury])[
  The dataset loaded through the following Julia commands contains data on the vapor pressure of mercury as a function of the temperature.
  ```julia
  import RDatasets
  data = RDatasets.dataset("datasets", "pressure")
  ```
  Find a low-dimensional mathematical model of the form
  #numbered($
    p(T) = exp big((alpha_0 + alpha_1 T + alpha_2 T^2 + alpha_3 T^3))
  $)
  for the pressure as a function of the temperature.
  Plot the approximation together with the data.
  An example solution is given in~@fig:pressure.
] <exercise:approx>
// In the LaTeX source, this figure is inside the exercise
#figure(
  placement: top,
  image("/figures/pressure_model.pdf", width: 75%),
  caption: [Solution for @exercise:approx.],
) <fig:pressure>

#compexercise[
  Let $u colon [0, 2pi] -> real$ and
  #numbered($
    x_k = lr(((2k pi)/(2n+1))), wide k = 0, dots, 2n.
  $) <eq:interp_nodes_fourier>
  We wish to interpolate $u$ at these nodes using complex exponentials:
  $ hat(u) = sum_(k=-n)^n a_k ee^(ii k x). $
  Write a function
  ```julia
  function fourier_interpolate(u, x, X)
      # Your code comes here ...
  end
  ```
  which takes three arguments:
  - `u` is the function to interpolate;

  - `x` are the interpolation nodes, given by~@eq:interp_nodes_fourier in the test code below;
    you can assume that this array contains an odd number of elements.

  - `X` is a one-dimensional array of values on the $x$ axis.
  The function should return a one-dimensional array containing the values that $hat(u)$ takes when evaluated at the points contained in `X`.
  You can use the following code to test your function:
  ```julia
  import Plots
  n, m = 5, 1000
  x = 2π/(2n+1) * (0:2n)
  X = 2π/m * (0:m)

  u(x) = sign(x - π)
  # u(x) = exp(sin(x) + cos(5x))
  # u(x) = x^2 * (x - 2π)^2 / π^4

  @time U = fourier_interpolate(u, x, X)
  Plots.plot(X, u.(X), label="u(x)", legend=:bottomright)
  Plots.plot!(X, U, label="û(x)")
  Plots.scatter!(x, u.(x), label="Interpolation points")
  Plots.xlims!(0, 2π)
  ```
  An example of the output plot is illustrated in~@fig:fourier.
] <exercise:fourier_interp>
#figure(
  placement: top,
  image("/figures/fourier.pdf", width: 80%),
  caption: [Example solution for~@exercise:fourier_interp.],
) <fig:fourier>

== Discussion and bibliography
A comprehensive study of approximation theory would require to cover the $L^infinity$ setting
as well as other functional settings.
A pillar of $L^infinity$ approximation theorem is Chebyshev's equioscillation theorem, which we alluded to in~@remark:cheb.
An excellent introductory reference on approximation theory is~@magnus (in French).
See also~@MR2265914[Chapter 10] and the references therein.

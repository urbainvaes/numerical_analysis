#import "../template.typ": *
#import "@preview/cetz:0.4.2"

// \d and \D: the template versions (upright strings) are followed by a space, as
// for a multi-letter text; upright letters are spaced like ordinary symbols, as in TeX
#let dd = $upright(d)$
#let DD = $upright(D)$

// \subsection*: unnumbered and absent from the table of contents
#let subsection-star(title) = heading(level: 3, numbering: none, outlined: false, title)

// Algorithm without line numbers (algorithmic environment without [1]);
// each line is a pair (indentation level, content)
#let algorithm-unnumbered(caption: none, ..lines) = figure(kind: "algorithm", supplement: [Algorithm],
  numbering: "1", outlined: false, caption: caption,
  grid(row-gutter: float-leading,
    ..lines.pos().map(((level, body)) => pad(left: 1em + level * 1.5em, body))))

// Figures 5.1 and 5.2 (TikZ with yscale = 0.65): the function f(x) = 0.1 x^2 - 1.5
// and its affine approximation y = slope (x - 8) + 4.9 at x_k = 8, plotted over [x0, 9.5]
#let affine-approximation-figure(slope, x0, name) = cetz.canvas({
  import cetz.draw: *
  let s = 0.65
  let p(x, y) = (x, s * y)
  let sep = 0.333em
  let blue = rgb(0, 0, 255)
  let red = rgb(255, 0, 0)
  set-style(stroke: 0.8pt)
  let arrow = (end: "triangle", fill: black, scale: 0.7)
  // Axes
  line(p(-1, 0), p(12, 0), mark: arrow)
  content(p(12, 0), anchor: "south", padding: sep, text(size: large, $x$))
  line(p(0, -2), p(0, 8), mark: arrow)
  content(p(0, 8), anchor: "west", padding: sep, text(size: large, $y$))
  // Function
  let f(x) = 0.1 * x * x - 1.5
  let n = 80
  let xs = range(n + 1).map(i => 1 + 8.5 * i / n)
  line(..xs.map(x => p(x, f(x))), stroke: (paint: blue, thickness: 1.6pt))
  content(p(9.5, f(9.5)), anchor: "east", padding: sep, text(fill: blue, $f(x)$))
  // Affine approximation
  content(p(8.2, 4.9), anchor: "west", padding: sep, text(fill: red, size: large, name))
  line(p(8, 0), p(8, 4.9), stroke: (paint: gray, thickness: 0.4pt, dash: "dotted"))
  circle(p(8, 4.9), radius: 2.5pt, fill: gray, stroke: none)
  let g(x) = slope * (x - 8) + 4.9
  line(p(x0, g(x0)), p(9.5, g(9.5)), stroke: (paint: red, thickness: 0.8pt, dash: "dashed"))
  // Ticks; the zero of the affine approximation is x_{k+1} = 8 - 4.9 / slope
  line(p(8, 0.1), p(8, -0.1))
  content(p(8, -0.1), anchor: "north", padding: sep, $x_k$)
  let root = 8 - 4.9 / slope
  line(p(root, 0.1), p(root, -0.1))
  content(p(root, -0.1), anchor: "north", padding: sep, frame: "rect", fill: white, stroke: none, $x_(k+1)$)
})

= Solution of nonlinear systems <cha:solution_of_nonlinear_systems>
#minitoc()

#section-star[Introduction]

This chapter concerns the numerical solution of nonlinear equations of the general form
#numbered($ vect(f)(vect(x)) = 0, wide vect(f) colon real^n -> real^n. $) <eq:nonlinear_equation>
A solution to this equation is called a _zero_ of the function $f$.
Except in particular cases (for example linear systems),
there does not exist a numerical method for solving~@eq:nonlinear_equation in a finite number of operations,
so iterative methods are required.

In contrast with the previous chapter,
it may not be the case that~@eq:nonlinear_equation admits one and only one solution.
For example, the equation $1 + x^2 = 0$ does not have a (real) solution,
and the equation $cos(x) = 0$ has infinitely many.
Therefore, convergence results usually contain assumptions on the function~$f$ that guarantee the existence and uniqueness of a solution in~$real^n$ or a subset of~$real^n$.

For an iterative method generating approximations $(vect(x)_k)_(k >= 0)$ of a root $vect(x)_*$,
we define the error as $vect(e)_k = vect(x)_k - vect(x)_*$.
If the sequence $(vect(x)_k)_(k >= 0)$ converges to $vect(x)_*$ in the limit as $k -> infinity$
and if
#numbered($ lim_(k -> infinity) norm(vect(e)_(k+1))/norm(vect(e)_k)^q = r, $) <eq:rate_of_convergence>
then we say that $(vect(x)_k)_(k >= 0)$ converges with _order of convergence_ $q$ and
_rate of convergence_ $r$.
In addition, we say that the convergence is linear $q = 1$,
and quadratic if $q = 2$.
The convergence is said to be superlinear if
#numbered($ lim_(k -> infinity) norm(vect(e)_(k+1))/norm(vect(e)_k) = 0. $) <eq:superlinear>
In particular,
the convergence is superlinear if the order of convergence is $q > 1$.

#remark[
  The notion of order of convergence may be defined also when
  the limit in~@eq:rate_of_convergence does not exist.
  A more general definition for the order of convergence of a sequence~$(vect(x)_k)_(k>=0)$ converging to $vect(x)_*$ is the following:
  $
    q(vect(x)_0) = inf lr(\{ p in [1, infinity) : limsup_(k -> infinity) norm(vect(e)_(k+1))/norm(vect(e)_k)^p = infinity \}),
  $
  or $q(vect(x)_0) = infinity$ if the numerator and denominator of the fraction are zero for sufficiently large~$k$.
  It is possible to define similarly the order of convergence of an iterative method
  for an initial guess in a neighborhood $V$ of $vect(x)_*$:
  $
    q = inf lr(\{ p in [1, infinity) : sup_(vect(x)_0 in V) lr((limsup_(k -> infinity) norm(vect(e)_(k+1))/norm(vect(e)_k)^p)) = infinity \}),
  $
  where the fraction should be interpreted as 0 if the numerator and denominator are zero.
  A more detailed discussion of this subject is beyond the scope of this course.
]

#indent The rest of chapter is organized as follows:
- In @sec:bisection_method,
  by way of introduction to the subject of numerical methods for nonlinear equations,
  we present and analyze the bisection method.

- In @sec:fixed_point_methods,
  we present a general method based on a fixed point iteration for solving~@eq:nonlinear_equation.
  The convergence of this method is analyzed in~@sec:convergence_fixed_point_methods.

- In @sec:examples_of_fixed_point_methods,
  two concrete examples of fixed point methods are studied:
  the chord method and the Newton--Raphson method.

== The bisection method <sec:bisection_method>
As an introduction to numerical methods for solving nonlinear equations,
we present the bisection method.
This method applies only in the case of a real-valued function $f colon real -> real$,
and relies on the knowledge of two points $a < b$ such that $f(a)$ and $f(b)$ have different signs.
By the intermediate value theorem,
there necessarily exists $x_* in (a, b)$ such that $f(x_*) = 0$.
The idea of the bisection method it to successively divide the interval in two equal parts,
and to retain, based on the sign of $f$ at the midpoint $x_(1\/2)$,
the one that necessarily contains a root.
If $f(x_(1\/2)) f(a) >= 0$, then $f(x_(1\/2)) f(b) <= 0$ and so there necessarily exists a root of $f$ in the interval $[x_(1\/2), b)$ by the intermediate value theorem.
In contrast, if $f(x_(1\/2)) f(a) < 0$, then there necessarily is a root in the interval $(a, x_(1\/2))$.
The algorithm is presented in @algo:bisection.
#algorithm-unnumbered(
  caption: [Bisection method],
  (0, [Assume that $f(a) f(b) < 0$ with $a < b$.]),
  (0, [Pick $epsilon > 0$.]),
  (0, [$x <- a\/2 + b\/2$]),
  (0, [*while* $|b - a| >= epsilon$ *do*]),
  (1, [*if* $f(x) f(a) >= 0$ *then*]),
  (2, [$a <- x$]),
  (1, [*else*]),
  (2, [$b <- x$]),
  (1, [*end if*]),
  (1, [$x <- a\/2 + b\/2$]),
  (0, [*end while*]),
) <algo:bisection>

#indent The following result establishes the convergence of the method.
#proposition[
  Assume that $f colon real -> real$ is a continuous function and $f(a) f(b) < 0$.
  Let~$[a_j, b_j]$ denote the interval obtained after $j$ iterations of the bisection method,
  and let $x_j = (a_j + b_j)\/2$ denote the midpoint of the interval,
  Then there exists a root~$x_*$ of $f$ such that
  #numbered($ abs(x_j - x_*) <= (b_0 - a_0) 2^(-(j+1)). $) <eq:error_bisection>
] <proposition:convergence_bisection>
#proof[
  By construction, $f(a_j) f(b_j) <= 0$ and $f(b) != 0$.
  Therefore, by the intermediate value theorem,
  there exists a root of $f$ in the interval $[a_j, b_j)$,
  implying that
  $ abs(x_j - x_*) <= (b_j - a_j)/2. $
  Since $b_j - a_j = 2^(-j) (b_0 - a_0)$,
  the statement follows.
]
Although the limit in~@eq:rate_of_convergence may not be well-defined (for example, $x_1$ may be a root of~$f$),
the error $x_j - x_*$ is bounded in absolute value by the sequence $(tilde(e)_j)_(j >= 0)$,
where $tilde(e)_j = (b_0 - a_0) 2^(-(j+1))$ by @proposition:convergence_bisection.
Since the latter sequence exhibits linear convergence to 0,
the convergence of the bisection method is said to be linear,
by a slight abuse of terminology.

== Fixed point methods <sec:fixed_point_methods>

Let $vect(x)_*$ denote a zero of the function $vect(f)$.
The idea of iterative methods for~@eq:nonlinear_equation is to construct,
starting from an initial guess $vect(x)_0$,
a sequence $(vect(x)_k)_(k=0, 1, dots)$ approaching~$vect(x)_*$.
A number of iterative methods for solving ~@eq:nonlinear_equation are based on an iteration of the form
#numbered($ vect(x)_(k+1) = vect(F)(vect(x)_k), $) <eq:fixed_point>
for an appropriate continuous function $F$.
Assume that $vect(x)_k$ converges to some point $vect(x)_* in real^n$ in the limit as $k -> infinity$.
Then, taking the limit $k -> infinity$ in~@eq:fixed_point,
we find that $vect(x)_*$ satisfies
$ vect(F)(vect(x)_*) = vect(x)_*. $
Such a point~$vect(x)_*$ is called a _fixed point_ of the function $vect(F)$.
Several definitions of the function~$vect(F)$ can be employed in order to ensure that
a fixed point of $vect(F)$ coincides with a zero of $vect(f)$.
One may, for example, define $vect(F)(vect(x)) = vect(x) - alpha^(-1) vect(f)(vect(x))$,
for some nonzero scalar coefficient $alpha$.
Then~$vect(F)(vect(x)_*) = vect(x)_*$ if and only if $vect(f)(vect(x)_*) = 0$.
Later in this chapter,
in @sec:examples_of_fixed_point_methods,
we study two instances of numerical methods which can be recast in the form~@eq:fixed_point.
Before this,
we study the convergence of the iteration~@eq:fixed_point for a general function~$vect(F)$.

== Convergence of fixed point methods <sec:convergence_fixed_point_methods>
Equation~@eq:fixed_point may be viewed as a _discrete-time_ dynamical system.
In order to study the behavior of the system as $k -> infinity$,
it is important to understand the concept of stability of a fixed point.
The concept of stability appears also in the field of ordinary differential equations,
which are _continuous-time_ dynamical systems.
Before we define this concept,
we introduce the following notation
for the open ball of radius $delta$ around $vect(x) in real^n$:
$ B_delta (vect(x)) := big(\{ vect(y) in real^n : norm(vect(y) - vect(x)) < delta \}). $
#definition(title: [Stability of fixed points])[
  Let $(vect(x)_k)_(k>=0)$ denote iterates obtained from~@eq:fixed_point when starting from an initial vector~$vect(x)_0$.
  Then we say that a fixed point $vect(x)_*$ is
  - an _attractor_ if there exists a neighborhood $cal(V)$ of $vect(x)_*$ such that
    #numbered($
      forall vect(x)_0 in cal(V), wide
      vect(x)_k stretch(->)_(k -> infinity) vect(x)_*.
    $) <eq:nonlinear_attractor>
    The largest neighborhood for which this is true,
    i.e. the set of values of $vect(x)_0$ such that~@eq:nonlinear_attractor holds true,
    is called the basin of attraction of $vect(x)_*$.

  - stable (in the sense of Lyapunov) if for all $epsilon > 0$,
    there exists $delta > 0$ such that
    $
      forall vect(x)_0 in B_delta (vect(x)_*), wide
      forall k in nat, wide
      norm(vect(x)_k - vect(x)_*) < epsilon.
    $

  - asymptotically stable if it is stable and an attractor.

  - exponentially stable if there exists $C > 0$, $alpha in (0, 1)$, and $delta > 0$ such that
    $
      forall vect(x)_0 in B_delta (vect(x)_*),
      quad forall k in nat, wide
      norm(vect(x)_k - vect(x)_*) <= C alpha^k norm(vect(x)_0 - vect(x)_*).
    $

  - globally exponentially stable if there exists $C > 0$ and $alpha in (0, 1)$ such that
    $
      forall vect(x)_0 in real^n,
      quad forall k in nat, wide
      norm(vect(x)_k - vect(x)_*) <= C alpha^k norm(vect(x)_0 - vect(x)_*).
    $
  - unstable if it is not stable.
]
Clearly, global exponential stability implies exponential stability,
which itself implies asymptotic stability and stability.
If $vect(x)_*$ is globally exponentially stable,
then $vect(x)_*$ is the unique fixed point of~$vect(F)$;
showing this is the aim of~@exercise:global_exponential_stability.
If $vect(x)_*$ is an attractor,
then the dynamical system~@eq:fixed_point is said to be locally convergent to~$vect(x)_*$.
The larger the basin of attraction of $vect(x)_*$,
the less careful we need to be when picking the initial guess~$vect(x)_0$.
Global exponential stability of a fixed point can sometimes be shown
provided that $vect(F)$ satisfies a strong hypothesis.

#definition(title: [Lipschitz continuity])[
  A function $vect(F) colon real^n -> real^n$ is said to be _Lipschitz_ continuous
  with constant $L$ if
  $
    forall (vect(x), vect(y)) in real^n times real^n, wide
    big(||vect(F)(vect(y)) - vect(F)(vect(x))||) <= L norm(vect(y) - vect(x)).
  $
]
A function $vect(F) colon real^n -> real^n$ that is Lipschitz continuous with a constant $L < 1$ is called a _contraction_.
For such a function, it is possible to prove that~@eq:fixed_point has a unique globally exponentially stable fixed point.
#theorem[
  Assume that $vect(F)$ is a contraction.
  Then there exists a unique fixed point of~@eq:fixed_point,
  and it holds that
  #numbered($
    forall vect(x)_0 in real^n,
    quad forall k in nat, wide
    norm(vect(x)_k - vect(x)_*) <= L^k norm(vect(x)_0 - vect(x)_*).
  $) <eq:global_exp>
] <theorem:exponential_convergence_fixed_point>
#proof[
  Existence and uniqueness of the fixed point follows from the _Banach fixed point theorem_,
  see~@theorem:banach_fixed_point,
  so here we show only global exponential convergence.
  Since $vect(F)$ is a contraction,
  it holds that
  #numbered($
    big(||vect(x)_k - vect(x)_*||)
    = big(||vect(F)(vect(x)_(k-1)) - vect(F)(vect(x)_*)||)
    <= L norm(vect(x)_(k-1) - vect(x)_*)
    <= dots <= L^k norm(vect(x)_0 - vect(x)_*),
  $) <eq:contraction>
  which proves~@eq:global_exp.
]

#indent It is possible to prove a weaker, local result under a less restrictive assumptions on the function~$vect(F)$.
#theorem[
  Assume that $vect(x)_*$ is a fixed point of~@eq:fixed_point and that $vect(F) colon real^n -> real^n$ satisfies the local Lipschitz condition
  #numbered($
    forall vect(x) in B_delta (vect(x)_*), wide
    big(||vect(F)(vect(x)) - vect(F)(vect(x)_*)||) <= L norm(vect(x) - vect(x)_*),
  $) <eq:local_lipschitz>
  with $0 <= L < 1$ and $delta > 0$.
  Then~$vect(x)_*$ is the unique fixed point of~$vect(F)$ in $B_delta (vect(x)_*)$ and,
  for all~$vect(x)_0 in B_delta (vect(x)_*)$, it holds that
  - All the iterates $(vect(x)_k)_(k in nat)$ belong to $B_delta (vect(x)_*)$.
  - The sequence $(vect(x)_k)_(k in nat)$ converges exponentially to $vect(x)_*$.
] <theorem:local_convergence>
#proof[
  See~@exercise:prove_local_convergence.
]
It is possible to guarantee that condition~@eq:local_lipschitz holds provided that
we have sufficiently good control of the derivatives of the function~$vect(F)$.
The function $vect(F)$ is said to be differentiable at $vect(x)$ (in the sense of Fréchet) if
there exists a linear operator $DD vect(F)_(vect(x)) colon real^n -> real^n$ such that
#numbered($
  lim_(vect(h) -> 0) norm(vect(F)(vect(x) + vect(h)) - vect(F)(vect(x)) - DD vect(F)_(vect(x)) (vect(h)))/norm(vect(h))
  = 0.
$) <eq:limit_differentiability>
If $vect(F)$ is differentiable,
then all its first partial derivatives~$partial_j F_i$ exist and,
in addition, it holds that $DD vect(F)_(vect(x)) (vect(h)) = mat(J)_F (vect(x)) vect(h)$
where $mat(J)_F (vect(x))$ is the Jacobian matrix of~$vect(F)$ at $vect(x)$:
$
  mat(J)_F (vect(x)) =
  pmat(
    partial_1 F_1(vect(x)), …, partial_n F_1(vect(x));
    ⋮, ⋱, ⋮;
    partial_1 F_n (vect(x)), …, partial_n F_n (vect(x))
  ).
$

#proposition[
  Let $vect(x)_*$ be a fixed point of~@eq:fixed_point,
  and assume that there exists $delta$ and a subordinate matrix norm such that
  $vect(F)$ is differentiable everywhere in $B_delta (vect(x)_*)$ and
  $
    forall vect(x) in B_delta (vect(x)_*), wide
    norm(mat(J)_F (vect(x))) <= L < 1.
  $
  Then condition~@eq:local_lipschitz is satisfied in the associated vector norm,
  and so the fixed point $vect(x)_*$ is locally exponentially stable.
] <proposition:local_convergence_fixed_point>
#proof[
  Let $vect(x) in B_delta (vect(x)_*)$.
  By the fundamental theorem of calculus and the chain rule,
  we have
  $
    vect(F)(vect(x)) - vect(F)(vect(x)_*)
    &= integral_0^1 dd/(dd t) Big((vect(F) big((vect(x)_* + t(vect(x) - vect(x)_*))))) dd t
    = integral_0^1 mat(J)_F big((vect(x)_* + t(vect(x) - vect(x)_*))) lr((vect(x) - vect(x)_*)) dd t.
  $
  Therefore,
  it holds that
  $
    norm(vect(F)(vect(x)) - vect(F)(vect(x)_*))
    <= integral_0^1 lr(||mat(J)_F big((vect(x) + t(vect(x) - vect(x)_*)))||) thin dd t thin norm(vect(x) - vect(x)_*)
    <= integral_0^1 L thin dd t thin norm(vect(x) - vect(x)_*) = L norm(vect(x) - vect(x)_*),
  $
  which is the statement.
]
#remark[
  As a student observed during the lecture,
  in dimension~$n = 1$, @proposition:local_convergence_fixed_point can be proved by using the mean value theorem:
  since $F$ is differentiable in~$(x_* -delta, x_* + delta)$,
  there exists for all $x$ in this interval a $xi = xi(x)$ also in this interval such that
  $ F(x) - F(x_*) = F'(xi) (x - x_*). $
  It then follows immediately that
  $
    big(|F(x) - F(x_*)|)
    = lr(|F'(xi) (x - x_*)|)
    <= L |x - x_*|.
  $
  This proof does not carry over to the multi-dimensional setting, however.
]

#indent In fact, it is possible to prove that a fixed point~$vect(x)_*$ is exponentially locally stable under an even weaker condition,
involving only the derivative of $vect(F)$ at $vect(x)_*$.
#proposition[
  Let $vect(x)_*$ be a fixed point of~@eq:fixed_point and that $F$ is differentiable at $vect(x)_*$ with
  $ norm(mat(J)_F (vect(x)_*)) = L < 1, $
  in a subordinate vector norm.
  Then the fixed point $vect(x)_*$ is locally exponentially stable.
] <proposition:local_convergence>
#proof[
  In this proof, the vector norm used is that associated with the matrix norm in the statement of the proposition.
  By the definition of differentiability~@eq:limit_differentiability,
  there exists for all~$epsilon > 0$ a~$delta > 0$ such that
  $
    forall vect(x) in B_delta (vect(x)_*) backslash {vect(x)_*}, wide
    norm(vect(F)(vect(x)) - vect(F)(vect(x)_*) - mat(J)_F (vect(x)_*)(vect(x) - vect(x)_*))/norm(vect(x) - vect(x)_*) <= epsilon.
  $
  By the triangle inequality,
  this implies that for all $vect(x) in B_delta (vect(x)_*)$,
  $
    norm(vect(F)(vect(x)) - vect(F)(vect(x)_*))
    &<=
    norm(vect(F)(vect(x)) - vect(F)(vect(x)_*) - mat(J)_F (vect(x)_*)(vect(x) - vect(x)_*)) + norm(mat(J)_F (vect(x)_*)(vect(x) - vect(x)_*)) \
    &<=
    epsilon norm(vect(x) - vect(x)_*) + norm(mat(J)_F (vect(x)_*)) norm((vect(x) - vect(x)_*))
    = (L + epsilon) norm(vect(x) - vect(x)_*).
  $
  We have thus shown that for all~$epsilon > 0$,
  there exists $delta > 0$ such that condition~@eq:local_lipschitz is satisfied with constant~$L+epsilon$.
  By taking~$epsilon$ sufficiently small,
  we can ensure that $L + epsilon < 1$,
  and so the fixed point $vect(x)_*$ is locally exponentially stable by~@theorem:local_convergence.
]

#indent The estimate in~@theorem:exponential_convergence_fixed_point suggests that
when the fixed point iteration~@eq:fixed_point converges,
the convergence is linear.
While this is usually the case,
the convergence is superlinear if $mat(J)_F (vect(x)_*) = 0$.

#proposition[
  Assume that $vect(x)_*$ is a fixed point of~@eq:fixed_point and that~$mat(J)_F (vect(x)_*) = 0$.
  Then the convergence to $vect(x)_*$ is superlinear,
  in the sense that if $vect(x)_k -> vect(x)_*$ as $k -> infinity$,
  then
  $ lim_(k -> infinity) norm(vect(x)_(k+1) - vect(x)_*)/norm(vect(x)_k - vect(x)_*) = 0. $
] <proposition:superlinear_convergence>
#proof[
  By @proposition:local_convergence,
  there exists $delta > 0$ such that $(vect(x)_k)_(k>= 0)$ is a sequence converging to $vect(x)_*$ for all $vect(x)_0 in B_delta (vect(x)_*)$.
  It holds that
  $
    norm(vect(x)_(k+1) - vect(x)_*)/norm(vect(x)_k - vect(x)_*)
    = norm(vect(F)(vect(x)_k) - vect(F)(vect(x)_*))/norm(vect(x)_k - vect(x)_*)
    = norm(vect(F)(vect(x)_k) - vect(F)(vect(x)_*) - mat(J)_F (vect(x)_*) (vect(x)_k - vect(x)_*))/norm(vect(x)_k - vect(x)_*).
  $
  Since $vect(x)_k - vect(x)_* -> vect(0)$ as $k -> infinity$,
  the right-hand side converges to 0 by~@eq:limit_differentiability.
]

#indent Similarly, if there exist $delta > 0$, $C > 0$ and $q in (1, infinity)$ such that
#numbered($
  forall vect(x) in B_delta (vect(x)_*), wide
  norm(vect(F)(vect(x)) - vect(F)(vect(x)_*)) <= C norm(vect(x) - vect(x)_*)^q,
$) <eq:condition_quadratic_convergence>
then assuming that $(vect(x)_k)_(k>=0)$ converges to $vect(x)_*$,
it holds for sufficiently large $k$ that
$
  norm(vect(x)_(k+1) - vect(x)_*)/norm(vect(x)_k - vect(x)_*)^q
  = norm(vect(F)(vect(x)_k) - vect(F)(vect(x)_*))/norm(vect(x)_k - vect(x)_*)^q
  <= C.
$
In this case, the order of convergence is at least $q$.

== Examples of fixed point methods <sec:examples_of_fixed_point_methods>
As we mentioned in @sec:fixed_point_methods,
there are several choices for the function $vect(F)$ that guarantee
the equivalence $vect(F)(vect(x)) = vect(x) <=> vect(f)(vect(x)) = vect(0)$.

=== The chord method
In the case where $f$ is a function from $real$ to $real$,
the simplest approach, sometimes called the _chord method_, is to define
$ F(x) = x - alpha^(-1) f(x). $
The fixed point iteration~@eq:error_bisection in this case admits a simple geometric interpretation:
at each step, the function $f$ is approximated by the affine function $x |-> f(x_k) + alpha (x - x_k)$,
and the new iterate is defined as the zero of this affine function,
i.e.
#numbered($ x_(k+1) = x_k - alpha^(-1) f(x_k) = F(x_k). $) <eq:naive_fixed_point>
This is illustrated in @fig:chord_method.
By~@proposition:local_convergence,
a sufficient condition to ensure local convergence is that
#numbered($ abs(F'(x_*)) = abs(1 - alpha^(-1) f'(x_*)) < 1. $) <eq:sufficient_condition_fixed_point>
In order for this condition to hold true,
the slope $alpha$ must be of the same sign as $f'(x_*)$
and the inequality $abs(alpha) >= abs(f'(x_*))\/2$ must be satisfied.
If $f'(x_*) = 0$,
then the sufficient condition~@eq:sufficient_condition_fixed_point is never satisfied;
in this case, the convergence must be studied on a case-by-case basis.
By~@proposition:superlinear_convergence,
the convergence of the chord method is superlinear if~$alpha = f'(x_*)$.
In practice, the solution $x_*$ is unknown,
and so this choice is not realistic.
Nevertheless, the above reasoning suggests that, by letting the slope $alpha$ vary from iteration to iteration in such a manner that~$alpha_k$ approaches $f'(x_*)$ as $k -> infinity$,
fast convergence can be obtained.
This is precisely what the Newton--Raphson method aims to achieve;
see~@sub:newton_raphson
#figure(
  placement: top,
  affine-approximation-figure(0.7, 0, [Affine approximation]),
  caption: [
    Graphical illustration of an iteration of the chord method.
  ],
) <fig:chord_method>

When $vect(f)$ is a function from $real^n$ to $real^n$,
the above approach generalizes to
$
  x_(k+1) = vect(F)(vect(x)_k), wide
  vect(F)(vect(x)) = vect(x) - mat(A)^(-1) vect(f)(vect(x)),
$
where $mat(A)$ is an invertible matrix.
The geometric interpretation of the method in this case is the following:
at each step, the function $vect(f)$ is approximated by the affine function
$vect(x) |-> vect(x)_k + mat(A)(vect(x) - vect(x)_k)$,
and the next iterate is given by the unique zero of the latter function.
Superlinear convergence is achieved when~$mat(A) = mat(J)_f (vect(x)_*)$.
Notice that each iteration requires to calculate $vect(y) := mat(A)^(-1) vect(f)(vect(x)_k)$,
which is generally achieved by solving the linear system $mat(A) vect(y) = vect(f)(vect(x)_k)$.

=== The Newton--Raphson method <sub:newton_raphson>
Let us first consider the case of a function from $real$ to $real$.
A necessary condition for the Newton--Raphson method to apply is that $f$ is differentiable.
At each step, the function $f$ is approximated by the affine function
$x |-> f(x_k) + f'(x_k) (x - x_k)$ and the unique zero of this function is returned.
In other words, one iteration of the Newton--Raphson method reads
#numbered($ x_(k+1) = x_k - f'(x_k)^(-1) f(x_k). $) <eq:newton_raphson>
For this iteration to be well-defined,
it is necessary that $f'(x_k) != 0$.
The Newton--Raphson method may be viewed as a variation on~@eq:naive_fixed_point where the slope~$alpha$ is adapted as the simulation progresses.
If the method converges and $f'$ is continuous,
then $f'(x_k) -> f'(x_*)$ in the limit as~$k -> infinity$,
which is an indication that superlinear convergence could occur in view of our discussion in the previous section.
Equation~@eq:newton_raphson may be recast as a fixed point iteration of the form~@eq:error_bisection with
$ F(x) = x - f(x)/(f'(x)). $
If $x_*$ is a simple root of $f$, that is if $f(x_*) = 0$ and $f'(x_*) != 0$,
then $x_*$ is a fixed point of the function $F$.
If the function~$f$ is twice continuously differentiable,
then the convergence of the Newton--Raphson method is superlinear by~@proposition:superlinear_convergence,
because then
$ F'(x_*) = (f(x_*) f''(x_*))/(f'(x_*)^2) = 0. $
The geometric interpretation of the Newton--Raphson method in dimension 1 is the following:
at each step, the function $vect(f)$ is approximated by the affine function
$x |-> x_k + f'(x_k)(x - x_k)$,
which is _the tangent line to~$f$ at $x_k$_,
and the next iterate is given by the unique zero of the latter function.
This is illustrated in @fig:newton_raphson.
#figure(
  placement: none,
  affine-approximation-figure(1.6, 4.25, [Tangent]),
  caption: [
    Graphical illustration of a Newton--Raphson iteration.
    #text(size: footnotesize)[The code used to create this figure is based on the answer~#link("https://tex.stackexchange.com/a/551205/125558") on LaTeX~stack exchange].
  ],
) <fig:newton_raphson>

The Newton--Raphson method may be generalized to nonlinear equations in $real^n$ of the form~@eq:nonlinear_equation.
In this case $vect(F)(vect(x)) = vect(x) - mat(J)_f (vect(x))^(-1) vect(f)(vect(x))$,
and so an iteration of the method reads
#numbered($ vect(x)_(k+1) = vect(x)_k - mat(J)_f (vect(x)_k)^(-1) vect(f)(vect(x)_k). $) <eq:iteration_newton_raphson_nonlinear>
In the rest of this section,
we show that the iteration~@eq:iteration_newton_raphson_nonlinear is well-defined in a small neighborhood of a root of $vect(f)$ under appropriate assumptions,
and we demonstrate the _second order_ convergence of the method,
first in dimension 1 under simplifying assumption involving the second derivative of~$f$,
and then in the multi-dimensional setting under more general assumptions.

==== Convergence in the one-dimensional setting
We assume in this section that~$(x_k)_(k>= 0)$ is generated from the Newton--Raphson method~@eq:newton_raphson and prove the following result.
#theorem(title: [Quadratic convergence of Newton--Raphson])[
  Assume that $f in C^2(real)$ and that the following assumptions are satisfied:
  - The first derivative of $f$ is uniformly bounded away from zero:
    $ inf_(x in real) abs(f'(x)) = m > 0. $

  - The second derivative of $f$ is uniformly bounded from above in absolute value:
    $ sup_(x in real) abs(f''(x)) = M < infinity. $
  Then $f(x)$ has a unique root~$x_*$ and it holds for all initial $x_0 in real$ that
  #numbered($
    forall k in nat, wide
    abs(x_(k+1) - x_*) <= M/(2m) abs(x_k - x_*)^2.
  $) <eq:convergence_newton_raphson_dim1>
] <theorem:newton_raphson_1d>
#proof[
  By assumption, the function $f$ is continuous and either strictly increasing everywhere or strictly decreasing everywhere.
  Therefore there exists a unique root $x_* in real$ of~$f$.
  In order to prove~@eq:convergence_newton_raphson_dim1,
  we note that
  #numbered($
    x_(k+1) - x_*
    = x_k - f(x_k)/(f'(x_k)) - x_*
    = 1/(f'(x_k)) Big((f'(x_k)(x_k - x_*) - f(x_k))).
  $) <eq:error_nr_dim1>
  By Taylor's theorem, there is $xi in real$ such that
  $ f(x_*) = f(x_k) + f'(x_k) (x_* - x_k) + 1/2 f''(xi) (x_* - x_k)^2. $
  Since~$x_*$ is a root of~$f$,
  the left-hand side of this equation is zero.
  Combining this equation with~@eq:error_nr_dim1,
  we deduce that
  $ x_(k+1) - x_* = (f''(xi) (x_k - x_*)^2)/(2f'(x_k)). $
  Taking absolute values and using the assumptions gives
  $
    |x_(k+1) - x_*|
    <= M/(2m) (x_k - x_*)^2,
  $
  which concludes the proof.
]

#remark[
  As a corollary of~@theorem:newton_raphson_1d,
  we obtain that the Newton--Raphson method is convergent if
  $ abs(x_k - x_*) <= (2m)/M. $
]

#subsection-star[Convergence in the multi-dimensional setting~#moreinfo]
As a first step towards a proof of quadratic convergence for the Newton--Raphson method in the multi-dimensional setting,
we begin by proving the following preparatory lemma,
which we will then employ in the particular case where the matrix-valued function $mat(A)$ is equal to~$mat(J)_f$.

#lemma[
  Let $mat(A) colon real^n -> real^(n times n)$ denote a matrix-valued function on $real^n$ that is both continuous and nonsingular at $vect(x)_*$,
  and let $vect(f)$ be a function that is differentiable at $vect(x)_*$ where $vect(f)(vect(x)_*) = 0$.
  Then the function
  $ vect(G)(vect(x)) = vect(x) - mat(A)(vect(x))^(-1) vect(f)(vect(x)) $
  is well-defined in a neighborhood $B_delta (vect(x)_*)$ of $vect(x)_*$.
  In addition, $vect(G)$ is differentiable at $vect(x)_*$ with
  #numbered($ mat(J)_G (vect(x)_*) = mat(I) - mat(A)(vect(x)_*)^(-1) mat(J)_f (vect(x)_*). $) <eq:jacobian_G>
] <lemma:lemma_newton_raphson>
#proof[
  It holds that
  #numbered($
    mat(A)(vect(x))
    = Big((mat(A)(vect(x)_*) - big((mat(A)(vect(x)_*) - mat(A)(vect(x))))))
    = mat(A)(vect(x)_*) Big((mat(I) - mat(A)(vect(x)_*)^(-1) big((mat(A)(vect(x)_*) - mat(A)(vect(x)))))).
  $) <eq:matrix_function_newton_raphson>
  Let $beta = norm(mat(A)(vect(x)_*)^(-1))$ and $epsilon = (2 beta)^(-1)$.
  By continuity of the matrix-valued function $mat(A)$,
  there exists~$delta > 0$ such that
  $
    forall vect(x) in B_delta (vect(x)_*), wide
    norm(mat(A)(vect(x)) - mat(A)(vect(x)_*)) <= epsilon.
  $
  For $vect(x) in B_delta (vect(x)_*)$ we have $norm(mat(A)(vect(x)_*)^(-1) big((mat(A)(vect(x)_*) - mat(A)(vect(x))))) <= norm(mat(A)(vect(x)_*)^(-1)) norm(mat(A)(vect(x)_*) - mat(A)(vect(x))) <= beta epsilon = 1/2$,
  and so @lemma:linear_inverse_neumann implies that the second factor
  on the right-hand side of~@eq:matrix_function_newton_raphson is invertible with a norm bounded from above by $2$.
  Therefore, we deduce that $mat(A)(vect(x))$ is invertible with
  #numbered($
    forall vect(x) in B_delta (vect(x)_*), wide
    norm(mat(A)(vect(x))^(-1)) <= 2 norm(mat(A)(vect(x)_*)^(-1)) = 2 beta,
  $) <eq:inverse_newton_raphson>
  which shows that~$vect(G)$ is well-defined in~$B_delta (vect(x)_*)$.
  In order to prove~@eq:jacobian_G,
  we need to show that
  $
    lim_(norm(vect(h)) -> 0) big(||vect(G)(vect(x)_* + vect(h)) - vect(G)(vect(x)_*) - big((mat(I) - mat(A)(vect(x)_*)^(-1) mat(J)_f (vect(x)_*))) vect(h)||)/norm(vect(h)) = 0
  $
  By definition of $vect(G)$,
  and using the fact that $vect(f)(vect(x)_*) = vect(0)$,
  we obtain that the argument of the norm in the numerator is equal to
  $
    & mat(A)(vect(x)_*)^(-1) vect(f)(vect(x)_*) - mat(A)(vect(x)_* + vect(h))^(-1) vect(f)(vect(x)_* + vect(h)) + mat(A)(vect(x)_*)^(-1) mat(J)_f (vect(x)_*) vect(h) \
    & wide =
    underbrace(big((mat(A)^(-1)(vect(x)_*) - mat(A)(vect(x)_* + vect(h))^(-1))) mat(J)_f (vect(x)_*) vect(h), =: vect(v)_1)
    - underbrace(mat(A)(vect(x)_* + vect(h))^(-1) big((vect(f)(vect(x)_* + vect(h)) - vect(f)(vect(x)_*) - mat(J)_f (vect(x)_*) vect(h))), =: vect(v)_2).
  $
  Noting that $mat(A)^(-1)(vect(x)_*) - mat(A)(vect(x)_* + vect(h))^(-1) = mat(A)(vect(x)_*)^(-1) big((mat(A)(vect(x)_* + vect(h)) - mat(A)(vect(x)_*))) mat(A)(vect(x)_* + vect(h))^(-1)$,
  we bound the norm of the first term on the right-hand side as follows:
  $
    forall vect(h) in B_delta (vect(0)), wide
    norm(vect(v)_1)
    <= 2 beta^2 norm(mat(A)(vect(x)_* + vect(h)) - mat(A)(vect(x)_*)) norm(mat(J)_f (vect(x)_*)) norm(vect(h)).
  $
  Clearly $norm(vect(v)_1) \/ norm(vect(h)) -> 0$ is the limit as $vect(h) -> vect(0)$ by continuity of the matrix function $mat(A)$.
  It also holds that $norm(vect(v)_2) \/ norm(vect(h)) -> 0$ by differentiability of $vect(f)$ at $vect(x)_*$,
  which concludes the proof.
]

#indent Using this lemma,
we can show the following result on the convergence of the multi-dimensional Newton--Raphson method.
#theorem(title: [Convergence of Newton--Raphson])[
  Let $vect(f) colon real^n -> real^n$ denote a function that is differentiable in a neighborhood~$B_delta (vect(x)_*)$
  of a point $vect(x)_*$ where $vect(f)(vect(x)_*) = 0$.
  Assume that the Jacobian matrix $mat(J)_f (vect(x))$ is nonsingular and continuous at $vect(x)_*$.
  Then $vect(x)_*$ is an attractor of the Newton--Raphson iteration~@eq:iteration_newton_raphson_nonlinear
  and the convergence is superlinear.

  In addition,
  if there is $alpha > 0$ such that the Lipschitz condition
  $
    forall vect(x) in B_delta (vect(x)_*), wide
    norm(mat(J)_f (vect(x)) - mat(J)_f (vect(x)_*)) <= alpha norm(vect(x) - vect(x)_*)
  $
  is satisfied,
  there exists $d in (0, delta)$ and $C > 0$ such that
  $
    forall vect(x)_k in B_d (vect(x)_*), wide
    norm(vect(x)_(k+1) - vect(x)_*) <= C norm(vect(x)_k - vect(x)_*)^2 .
  $
  In other words, the convergence is at least quadratic in~$B_d (vect(x)_*)$.
]
#proof[
  Using~@lemma:lemma_newton_raphson,
  we obtain that the Newton--Raphson update
  $ vect(F)(vect(x)) = vect(x) - mat(J)_f (vect(x))^(-1) vect(f)(vect(x)), $
  is well-defined in a neighborhood $B_delta (vect(x)_*)$ of $vect(x)_*$
  for sufficiently small $delta$.
  In addition, the second statement in~@lemma:lemma_newton_raphson gives that $mat(J)_F (vect(x)_*)^(-1) = mat(I) - mat(J)_F (vect(x)_*)^(-1) mat(J)_F (vect(x)_*) = 0$,
  which establishes the superlinear convergence by~@proposition:superlinear_convergence.

  In order to show that the convergence is quadratic,
  we begin by noticing that,
  since
  $
    vect(f)(vect(x)_k)
    = integral_0^t dd/(dd t) thin vect(f) big((vect(x)_* + t(vect(x)_k - vect(x)_x))) thin dd t
    = integral_0^t mat(J)_f big((vect(x)_* + t(vect(x)_k - vect(x)_x))) (vect(x)_k - vect(x)_*) thin dd t,
  $
  it holds for all $vect(x)_k in B_delta (vect(x)_*)$ that
  #eq-last($
    norm(vect(f)(vect(x)_k) - mat(J)_f (vect(x)_*) (vect(x)_k - vect(x)_*))
    &= lr(||integral_0^1 Big((mat(J)_f big((vect(x)_* + t(vect(x)_k - vect(x)_*))) - mat(J)_f (vect(x)_*))) thin (vect(x)_k - vect(x)_*) thin dd t||) \
    &<= integral_0^1 lr(||mat(J)_f big((vect(x)_* + t(vect(x)_k - vect(x)_*))) - mat(J)_f (vect(x)_*)||) thin norm(vect(x)_k - vect(x)_*) thin dd t \
    &<= integral_0^1 alpha t norm(vect(x)_k - vect(x)_*)^2 thin dd t
    <= alpha/2 norm(vect(x)_k - vect(x)_*)^2.
  $) <eq:inequality_newton_raphson>
  Let $d in (0, delta)$ be sufficiently small to ensure that
  $
    forall vect(x) in B_d (vect(x)_*),
    wide norm(mat(J)_f (vect(x))^(-1)) <= 2 norm(mat(J)_f (vect(x)_*)^(-1)).
  $
  There exists such a~$d$ by~@eq:inverse_newton_raphson.
  Using the inequality~@eq:inequality_newton_raphson,
  we have that for all $vect(x)_k in B_d (vect(x)_*)$,
  $
    & norm(vect(x)_(k+1) - vect(x)_*)
    = norm(vect(F)(vect(x)_k) - vect(x)_*)
    = norm(vect(x)_k - vect(x)_* - mat(J)_f (vect(x)_k)^(-1) vect(f)(vect(x)_k)) \
    & wide = lr(||mat(J)_f (vect(x)_k)^(-1) big((vect(f)(vect(x)_k) - mat(J)_f (vect(x)_k) (vect(x)_k - vect(x)_*)))||)
    <= lr(||mat(J)_f (vect(x)_k)^(-1)||) lr(||vect(f)(vect(x)_k) - mat(J)_f (vect(x)_k) (vect(x)_k - vect(x)_*)||) \
    & wide <= lr(||mat(J)_f (vect(x)_k)^(-1)||) Big((norm(vect(f)(vect(x)_k) - mat(J)_f (vect(x)_*) (vect(x)_k - vect(x)_*)) + norm(mat(J)_f (vect(x)_*) - mat(J)_f (vect(x)_k)) norm(vect(x)_k - vect(x)_*))) \
    & wide <= (3 alpha)/2 lr(||mat(J)_f (vect(x)_k)^(-1)||) norm(vect(x)_k - vect(x)_*)^2
    <= 3 alpha lr(||mat(J)_f (vect(x)_*)^(-1)||) norm(vect(x)_k - vect(x)_*)^2,
  $
  which concludes the proof.
]

=== The secant method~#moreinfo
The Newton--Raphson method exhibits very fast convergence,
but it requires the knowledge of the derivatives of the function $vect(f)$.
To conclude this chapter,
we describe a root-finding algorithm,
known as the secant method,
that enjoys superlinear convergence but does not require the derivatives of~$vect(f)$.
This method applies only when $vect(f)$ is a function from $real$ to $real$,
and so we drop the vector notation in the rest of this section.

Unlike the other methods presented so far in~@sec:fixed_point_methods,
the secant method _can not_ be recast as a fixed point iteration of the form~$x_(k+1) = F(x_k)$.
Instead, it is of the more general form~$x_(k+2) = F(x_k, x_(k+1))$.
The geometric intuition behind the method in the following: given $x_k$ and $x_(k+1)$,
the function $f$ is approximated by the unique linear function that passes through $big((x_k, f(x_k)))$ and
$big((x_(k+1), f(x_(k+1))))$,
and the iterate $x_(k+2)$ is defined as the root of this linear function.
In other words, $f$ is approximated as follows:
$ tilde(f)(x) = f(x_k) + (f(x_(k+1)) - f(x_k))/(x_(k+1) - x_k) (x - x_k). $
Solving $tilde(f)(x) = 0$ gives the following expression for $x_(k+2)$:
#numbered($ x_(k+2) = (f(x_(k+1)) x_k - f(x_k) x_(k+1))/(f(x_(k+1)) - f(x_k)), $) <eq:secant_method>
Showing the convergence of the secant method rigorously under general assumptions is tedious,
so in this course we restrict our attention to the case where $f$ is a quadratic function.
Extending the proof of convergence to a more general smooth function can be achieved
by using a quadratic Taylor approximation of~$f$ around the root~$x_*$,
which is accurate in a close neighborhood of~$x_*$.

#theorem(title: [Convergence of the secant method])[
  Assume that $f$ is a convex quadratic polynomial with a simple root at $x_*$
  and that the secant method converges: $lim_(k -> infinity) x_k = x_*$.
  Then the order of convergence is given by the golden ratio
  $ phi = (1 + sqrt(5))/2. $
  More precisely, there exists a positive real number $y_infinity$ such that
  #numbered($ lim_(k -> infinity) abs(x_(k+1) - x_*)/abs(x_k - x_*)^phi = y_infinity. $) <eq:convergence_secant>
]
#proof[
  Equation~@eq:secant_method implies that
  $ x_(k+2) - x_* = (f(x_(k+1)) (x_k - x_*) - f(x_k) (x_(k+1)-x_*))/(f(x_(k+1)) - f(x_k)). $ <eq:secant_error>
  By assumption, the function $f$ may be expressed as
  $ f(x) = lambda (x - x_*) + mu (x - x_*)^2, wide lambda != 0. $
  // In LaTeX, the label of this unnumbered equation (equation*) refers to the current subsection
  Substituting this expression in~#context link(<eq:secant_error>, text(fill: darkblue)[(#numbering("1.1.1", ..counter(heading).at(<eq:secant_error>)))]) and letting $e_k = x_k - x_*$,
  we obtain
  $
    e_(k+2)
    = (mu e_k e_(k+1) (e_(k+1) - e_k))/(lambda (e_(k+1) - e_k) + mu (e_(k+1)^2 - e_k^2))
    = (mu e_k e_(k+1))/(lambda + mu (e_(k+1) + e_k)).
  $
  Rearranging this equation,
  we have
  #numbered($
    e_(k+2)/e_(k+1)
    = (mu e_k)/(lambda + mu (e_(k+1) + e_k)).
  $) <eq:secant_method_rearranged>
  By assumption, the right-hand side converges to zero,
  and so the left-hand side must also converge to zero;
  the convergence is superlinear.

  To conclude the proof,
  we first reason formally in order to guess the order convergence,
  and then give a rigorous proof that our guess is correct.
  If~$e_k$ is small, then it holds approximately by~@eq:secant_method_rearranged that
  #numbered($
    e_(k+2)/e_(k+1)
    approx mu e_k.
  $) <eq:recursion_error_secant>
  Assume that there exists $q > 0$ such that the equation $e_(k+1) = C e_k^q$ is valid for all~$k$.
  Then it holds that $e_(k+2) = C e_(k+1)^q = C (C e_k^q)^q$ and~@eq:recursion_error_secant enables to determine~$q$:
  $
    (C (C e_k^q)^q) / (C e_k^q) = mu/lambda e_k
    wide => wide
    C^q e_k^(q^2 - q) = mu/lambda e_k
    wide => wide
    q^2 - q - 1 = 0.
    wide => wide
    q = phi.
  $
  Now comes the rigorous justification.
  Take absolute values in~@eq:secant_method_rearranged to obtain,
  after rearranging,
  $
    abs(e_(k+2))/abs(e_(k+1))^phi
    = lr((abs(e_(k+1))/abs(e_k)^(1/(phi - 1))))^(1 - phi) mu/abs(lambda + mu (e_(k+1) + e_k))
    = lr((abs(e_(k+1))/abs(e_k)^phi))^(1 - phi) abs(mu)/abs(lambda + mu (e_(k+1) + e_k)),
  $
  where we used that $phi = 1/(phi - 1)$,
  since $phi$ is a root of the equation $phi^2 - phi - 1 = 0$.
  Thus, introducing the ratio $y_k = abs(e_(k+1))\/abs(e_k)^phi$,
  we have
  $ y_(k+1) = y_k^(1 - phi) abs(mu)/abs(lambda + mu (e_(k+1) + e_k)). $
  Taking logarithms in this equation,
  we deduce
  $
    log(y_(k+1)) = (1 - phi) log(y_k) + c_k,
    wide c_k := log lr((abs(mu)/abs(lambda + mu (e_(k+1) + e_k)))).
  $
  This is a recurrence equation for $log(y_k)$,
  whose explicit solution can be obtained from the variation-of-constants formula:
  $ log(y_k) = (1 - phi)^(k-1) log(y_1) + sum_(i=1)^(k-1) (1 - phi)^(k-1-i) c_i. $
  Since $(c_k)_(k >= 0)$ converges to the constant $c_infinity = log abs(mu\/lambda)$ by the assumption that $e_k -> 0$,
  the sequence $big((log(y_k)))_(k>= 0)$ converges to $c_infinity \/ phi$ (prove this!).
  Therefore, by continuity of the exponential function,
  it holds that
  $
    y_k = exp big((log(y_k))) stretch(->)_(k -> infinity) exp lr((c_infinity/phi))
    = lr(|mu/lambda|)^(1/phi)
  $
  and so we deduce~@eq:convergence_secant.
]

== A numerical experiment
To conclude this chapter,
we present the results of a numerical experiment.
Specifically,
we consider four different fixed point methods for calculating the square root of 2,
i.e.~for solving the nonlinear equation
#numbered($ f(x) := x^2 - 2 = 0. $) <eq:numerical_experiment>
The unique positive solution to this equation is $x_* = sqrt(2)$.
The methods we consider are the following:
- The chord method with large~$alpha = 10$.

- The chord method with the optimal parameter~$alpha$,
  which is such that $F'(x_*) = 0$.
  The optimum value for $alpha$ for solving~@eq:numerical_experiment is given by~$alpha_* = 2 sqrt(2)$.

- The Newton--Raphson method,
  where each iteration is given by
  $
    x_(k+1) = x_k - f(x_k)/(f'(x_k)) = x_k - (x_k^2 - 2)/(2x_k) = 1/2 lr((x_k + 2/x_k))
    =: F(x_k),
  $
  with
  $ F(x) = 1/2 lr((x + 2/x)). $
  Notice that $F'(x_*) = 0$ and that $F in C^2 big(((0, infinity)))$.
  Therefore, by Taylor's theorem it holds for all $x in (x_*-1, x_* + 1)$ that
  $
    lr(|F(x) - F(x_*)|)
    = lr(|F'' big((xi(x)))|)
    <= L (x - x_*)^2,
    wide L := sup_(abs(x - x_*) <= 1) abs(F''(x)).
  $
  We deduce that the convergence is at least quadratic by~@eq:condition_quadratic_convergence.
#remark[
  Note that the ancient Babylonian method coincides with the Newton--Raphson method applied to~@eq:numerical_experiment.
]
The following code implements these methods.
Note that we use the arbitrary precision `BigFloat` format with a precision we manually set to $2000$ bits,
which enables using a very small~$epsilon$ in the stopping criterion.
```julia
function count_digits(x, y)
    xdigits = split(string(x), "")
    ydigits = split(string(y), "")
    len = min(length(xdigits), length(ydigits))
    for i in 1:len
        xdigits[i] != ydigits[i] && return i-2
    end
end

function my_sqrt(a)
    exact = sqrt(a)
    f(x) = x*x - a
    fp(x) = 2x

    # Uncomment desired line
    F(x) = x - f(x)/10      # Chord method
    # F(x) = x - f(x)/(2√a)   # Chord method with optimal α
    # F(x) = 1/2 * (x + a/x)  # Babylonian / Newton Raphson

    r, ε = 1, 1e-200
    while abs(f(r)) > ε
        r = F(r)
        digits = ceil(Int, -log10(abs(r - exact)))
        println("Number of correct digits: $digits")
    end
end

# Sets the precision of BigFloats to 1000 bits
setprecision(2000)
my_sqrt(BigFloat(2))
```
For each of the methods,
the number of correct digits of the approximation as the iterations progress is illustrated in~@table:numerical_experiment_square_root.
Observe that for all the methods except the first one,
the number of correct digits is approximately doubled at each iteration,
which is consistent with quadratic convergence.
#figure(
  placement: none,
  kind: table,
  tabular(
    columns: 4,
    [Method], [Chord $alpha = 10$], [Chord $alpha = 2 sqrt(2)$], [Newton--Raphson],
    [\# Iterations], [*1357*], [*8*], [*9*],
    [\# Correct digits $i = 1$], [1], [1], [1],
    [\# Correct digits $i = 2$], [1], [3], [3],
    [\# Correct digits $i = 3$], [1], [6], [6],
    [\# Correct digits $i = 4$], [1], [12], [12],
    [\# Correct digits $i = 5$], [1], [26], [24],
    [\# Correct digits $i = 6$], [1], [53], [48],
    [\# Correct digits $i = 7$], [1], [107], [97],
    [\# Correct digits $i = 8$], [1], [214], [196],
    [\# Correct digits $i = 9$], [1], [n/a], [392],
  ),
  caption: [
    Comparison of different fixed point methods for calculating $sqrt(2)$.
    Here $i$ denotes the iteration index.
  ],
) <table:numerical_experiment_square_root>

== Exercises

#compexercise[
  Implement the bisection method for finding the solution(s) to the equation
  $ x = cos(x). $
]

#exercise[
  Find a discrete-time dynamical system over $real$ of the form
  $ x_(k+1) = F(x_k) $
  for which $0$ is an attractor but is not stable.

  #par(first-line-indent: 0pt)[*Hint:* Use a function $F$ that is discontinuous.]
]

#exercise[
  Show that if $vect(x)_*$ is a globally exponentially stable fixed point of~$F$,
  then~$F$ does not have any other fixed point: $vect(x)_*$ is the unique fixed point.
] <exercise:global_exponential_stability>

#exercise[
  Prove @theorem:local_convergence.
] <exercise:prove_local_convergence>

#exercise[
  Let $vect(x)_*$ be a fixed point of~@eq:fixed_point.
  Show that if
  $ rho big((mat(J)_F (vect(x)_*))) < 1, $
  then $vect(x)_*$ is locally exponentially stable.
  It is sufficient by~@proposition:local_convergence to find a subordinate matrix norm such that $norm(mat(J)_F (vect(x)_*)) < 1$.
  In other words, this exercise amounts to showing that for any matrix $mat(A) in real^(n times n)$ with $rho(mat(A)) < 1$,
  there exists a matrix norm such that $norm(mat(A)) < 1$.

  #par(first-line-indent: 0pt)[*Hint:* One may employ a matrix norm of the form $norm(mat(A))_(mat(T)) := norm(mat(T)^(-1) mat(A) mat(T))_2$,
  which is a subordinate norm by~@exercise:induced_matrix_norm.
  The Jordan normal form is useful for constructing the matrix $mat(T)$,
  and equation~@eq:preliminary_equation is also useful.]
] <exercise:spectral_and_norm>

#solution[
  Let $mat(J) = mat(P)^(-1) mat(A) mat(P)$ denote the Jordan normal form of~$mat(A)$,
  and let
  $
    mat(E)_epsilon =
    pmat(
      epsilon, , , ;
      , epsilon^2, , ;
      , , ⋱, ;
      , , , epsilon^n
    )
  $
  By~Eq.~@eq:preliminary_equation,
  the matrix $mat(J)_epsilon := mat(E)_epsilon^(-1) mat(J) mat(E)_epsilon$ coincides with~$mat(J)$,
  except that the first superdiagonal is multiplied by~$epsilon$.
  Let $mat(D)$ denote the diagonal part of $mat(J)_epsilon$.
  We have that
  $ norm(mat(J)_epsilon - mat(D))_2 = sqrt(lambda_max (mat(E)_epsilon^T mat(E)_epsilon)). $
  The matrix $mat(E)_epsilon^T mat(E)_epsilon$ is diagonal with entries equal to either 0 or $epsilon^2$,
  and so $norm(mat(J)_epsilon - mat(D))_2 < epsilon$.
  By the triangle inequality, we have
  #numbered($
    norm(mat(J)_epsilon) <= norm(mat(D)) + norm(mat(J)_epsilon - mat(D))_2
    <= rho(mat(A)) + epsilon.
  $) <eq:alternative_norm_spectral_radius>
  Let $norm(mat(A))_epsilon := norm(mat(E)_epsilon^(-1) mat(P)^(-1) mat(A) mat(P) mat(E)_epsilon)$.
  By~@exercise:induced_matrix_norm with~$mat(T) = mat(P) mat(E)_epsilon$,
  this is indeed a subordinate matrix norm.
  By~@eq:alternative_norm_spectral_radius and the assumption that~$rho(mat(A)) < 1$,
  it is clear that $norm(mat(A))_epsilon < 1$ provided that~$epsilon$ is sufficiently small.
]
#remark[
  A corollary of @exercise:induced_matrix_norm is that the spectral radius of a matrix~$mat(A)$ is _the infimum_ of $norm(mat(A))$ over all subordinate matrix norms.
]

#exercise[
  Calculate $x = root(3, 3 + root(3, 3 + root(3, 3 + sqrt(dots))))$ using the bisection method.
]

#exercise[
  Solve the equation $f(x) = ee^x - 2 = 0$ using a fixed point iteration of the form
  $
    x_(k+1) = F(x_k), wide
    F(x) = x - alpha^(-1) f(x).
  $
  Using your knowledge of the exact solution $x_* = log 2$,
  write a sufficient condition on $alpha$ to guarantee that $x_*$ is locally exponentially stable.
  Verify your findings numerically and plot,
  using a logarithmic scale for the $y$ axis,
  the error in absolute value as a function of $k$.
]

#exercise[
  Implement the Newton--Raphson method for solving $f(x) = ee^x - 2 = 0$,
  and plot the error in absolute value as a function of the iteration index $k$.
]

#exercise[
  Find the point $(x, y)$ on the parabola $y = x^2$ that is closest to the point $(3, 1)$.
]

#exercise[
  Consider the linear system
  $
    cases(
      & y = (x-1)^2,
      & x^2 + y^2 = 4
    )
  $
  By drawing these two constraints in the $x y$ plane,
  find an approximation of the solution(s).
  Then calculate the solution(s) using a fixed-point method.
]

#exercise[
  Find solutions $(psi, lambda)$, with $lambda > 0$,
  to the following eigenvalue problem:
  $ psi'' = - lambda^2 psi, wide psi(0) = 0, wide psi'(1) = psi(1). $
]

#exercise[
  Suppose that we have $n$ data points $(x_i, y_i)$ of an unknown function~$y = f(x)$.
  We wish to approximate $f$ by a function of the form
  $ tilde(f)(x) = a/(b + x) $
  by minimizing the sum of squares
  $ sum_(i=1)^n big(|tilde(f)(x_i) - y_i|)^2. $
  Write a system of nonlinear equations that the minimizer $(a, b)$ must satisfy,
  and solve this system using the Newton--Raphson method starting from $(1, 1)$.
  The data is given below:
  ```julia
  x = [0.0; 0.1; 0.2; 0.3; 0.4; 0.5; 0.6; 0.7; 0.8; 0.9; 1.0]
  y = [0.6761488864859304; 0.6345697680852508; 0.6396283580587062; 0.6132010027973919;
       0.5906142598705267; 0.5718728461471725; 0.5524549902830562; 0.538938885654085;
       0.5373495476994958; 0.514904589752926; 0.49243437874655027]
  ```
  Plot the data points together with the function~$tilde(f)$ over the interval $[0, 1]$.
  Your plot should look like @fig:solution_exercise_approximation.
  #figure(
    placement: top,
    image("/figures/approx.pdf", width: 60%),
    caption: [Solution to @exercise:nonlinear_approximation.],
  ) <fig:solution_exercise_approximation>
] <exercise:nonlinear_approximation>

#exercise(title: [Nonlinear least-squares])[
  Suppose that we are given $n$ data points $(x_i, y_i)$ of an unknown function~$y = f(x)$.
  We wish to approximate $f$ by a straight line
  $ tilde(f)(x) = a x + b $
  by minimizing the sum of squared _Euclidean distances_ between the data points and the straight line~$tilde(f)$.
  Since the distance between a point $(x_i, y_i)$ and the straight line is given by
  $ |y_i - a x_i - b|/sqrt(1+a^2), $
  the objective function to minimize is given by
  $ J(a, b) := sum_(i=1)^n lr((y_i - a x_i - b))^2/(1+a^2). $
  This is a smooth function of~$a$ and $b$,
  and so a necessary condition for a pair $(a_*, b_*) in real^2$ to be a minimizer is that
  $ nabla J(a_*, b_*) = vect(0), $
  which is a nonlinear equation for the unknowns~$a_*$ and $b_*$.
  Solve this equation by using the Newton--Raphson method initialized at~$(1, 1)$,
  and then plot the data points together with the function~$tilde(f)$ over the interval~$[0, 1]$.
  Your plot should look like @fig:solution_exercise_approximation_least_squares.
  The data is given hereafter:
  ```julia
  x = [0.0; 0.1; 0.2; 0.3; 0.4; 0.5; 0.6; 0.7; 0.8; 0.9; 1.0]
  y = [-0.9187980789440975; -0.6159791344678258; -0.25568734869121856;
       -0.14269370171581808; 0.3094396057228459; 0.6318327173549161;
       0.8370437988106428; 1.0970402798788812; 1.6057799131867696;
       1.869090784869698; 2.075369730726694]
  ```
  #figure(
    placement: top,
    image("/figures/approx_line.pdf", width: 60%),
    caption: [Solution to @exercise:nonlinear_least_squares.],
  ) <fig:solution_exercise_approximation_least_squares>
] <exercise:nonlinear_least_squares>

== Discussion and bibliography

The content of this chapter is largely based on the lecture notes~@VanDooren.
Several of the exercises are taken or inspired from~@Legat.
The proof of convergence of the secant method is inspired from the general proof presented in the short paper~@MR1186462.
For a detailed treatment of iterative methods for nonlinear equations,
see the book~@MR1744713.

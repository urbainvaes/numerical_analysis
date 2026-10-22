#import "../template.typ": *

// \ip[\big]{...} and \ip[\Big]{...}
#let ipbig(..args) = big($chevron.l #args.pos().join($,$) chevron.r$)
#let ipBig(..args) = Big($chevron.l #args.pos().join($,$) chevron.r$)

= Optimization <cha:Optimization>
#minitoc()

In this chapter, we focus on optimization problems of the following form:
#numbered($
  "Find" quad vect(x)_* in argmin_(vect(x) in cal(K)) J(vect(x)),
$) <eq:optimization>
where $cal(K)$ is a given subset of $real^n$ and $J colon cal(K) -> real$ is a given _objective function_.
We came across several examples of such problems earlier in these notes:
- In @cha:interpolation_and_approximation,
  in the context of least-squares approximation,
  we considered the problem of minimizing
  $ J(vect(alpha)) = 1/2 norm(mat(A) vect(alpha) - vect(b))^2. $

- In @cha:solution_of_linear_systems,
  we observed that,
  if $mat(A)$ is a symmetric and positive definite matrix,
  then solving the linear system~$mat(A) vect(x) = vect(b)$
  amounts to finding the minimizer of the functional
  $ J(vect(x)) = 1/2 vect(x)^T mat(A) vect(x) - vect(b)^T vect(x). $
When $cal(K) = real^n$,
equation~@eq:optimization is an _unconstrained_ optimization problem,
and when~$cal(K) subset.neq real^n$,
equation~@eq:optimization is a _constrained_ optimization problem.
In practice,
the set~$cal(K)$ is often an intersection of sets of the form
$
  big({vect(x) in real^n : phi.alt(vect(x)) <= 0}),
  wide
  " or "
  wide
  big({vect(x) in real^n : phi.alt(vect(x)) = 0}),
$
for appropriate $phi.alt colon real^n -> real$.
Constraints of the former form are called _inequality constraints_,
while constraints of the latter form are called _equality constraints_.
Our aim in this chapter is to give a brief introduction to numerical optimization.
We focus on the simplest method,
namely the _steepest descent method_ with fixed step.
The rest of this chapter is organized as follows:
- We begin in~@sec:convexity by defining the notions of~_convexity_, _strict convexity_ and _strong convexity_,
  which play an important role in optimization.

- Then, in~@sec:unconstrained_optimization,
  we analyze the steepest descent method with fixed step in the setting of unconstrained optimization.
  To this end, we first establish conditions under which~@eq:optimization is well posed.

- Finally, in~@sec:constrained_optimization,
  we extend the steepest descent method to the case of optimization with constraints.

#remark[
  For generality, we could consider the setting where the set~$cal(K)$ in~@eq:optimization is a subset of some finite dimensional or infinite dimensional vector space~$V$.
  An optimization problem over (a subset of) a finite dimensional vector space of dimension~$n$ can always be recast as an optimization problem over (a subset of) $real^n$
  -- the type we study in this chapter -- by fixing a basis.
  The case of an infinite dimensional vector space, however,
  is more delicate, and we do not address it here.
]

== Definition and characterization of convexity <sec:convexity>

#definition(title: [Convexity])[
  Assume that $J colon cal(K) -> real$.
  - The function~$J$ is said to be _convex_ if
    #numbered($
      forall (vect(x), vect(y)) in cal(K) times cal(K),
      quad forall theta in [0, 1],
      quad
      J big((theta vect(x) + (1 - theta) vect(y)))
      <= theta J(vect(x)) + (1 - theta) J(vect(y)).
    $) <eq:convexity>

  - The function~$J$ is called _strictly convex_ if~@eq:convexity holds with strict inequality
    if $vect(x) != vect(y)$ and~$theta in (0, 1)$.

  - The function~$J$ is called _strongly convex_ with parameter~$alpha > 0$ if
    for all $(vect(x), vect(y)) in cal(K) times cal(K)$ and for all~$theta in [0, 1]$,
    #numbered($
      J big((theta vect(x) + (1 - theta) vect(y)))
      <= theta J(vect(x)) + (1 - theta) J(vect(y))
      - alpha/2 theta (1 - theta) norm(vect(x) - vect(y))^2.
    $) <eq:strong_convexity>
]

#indent If the function~$J$ is differentiable,
then convexity, strict convexity and strong convexity can be characterized in terms of the gradient~$nabla J$.
We illustrate this for strong convexity,
noting that a characterization of convexity is obtained by substituting~$alpha = 0$ in the following result.

#proposition[
  A differentiable function~$J colon real^n -> real$ is strongly convex with parameter~$alpha$ if and only if
  #numbered($
    forall (vect(x), vect(y)) in real^n times real^n,
    wide J(vect(x)) >= J(vect(y)) + ipbig(nabla J(vect(y)), vect(x) - vect(y)) + alpha/2 norm(vect(x) - vect(y))^2,
  $) <eq:strong_convex_eq1>
  or, equivalently,
  #numbered($
    forall (vect(x), vect(y)) in real^n times real^n,
    wide ipbig(nabla J(vect(x)) - nabla J(vect(y)), vect(x) - vect(y))
    >= alpha norm(vect(x) - vect(y))^2.
  $) <eq:strong_convex_eq2>
]
#proof[
  For clarity, we divide the proof into items and prove one implication per item.
  - @eq:strong_convexity $=>$~@eq:strong_convex_eq1.
    Rearranging~@eq:strong_convexity, we have
    $
      (J big((vect(y) + theta (vect(x) - vect(y)))) - J(vect(y)))/theta
      <= J(vect(x)) - J(vect(y)) - alpha/2 (1 - theta) norm(vect(x) - vect(y))^2.
    $
    Taking the limit~$theta -> 0$,
    we deduce that
    $
      ipbig(nabla J(vect(y)), vect(x) - vect(y))
      <= J(vect(x)) - J(vect(y)) - alpha/2 norm(vect(x) - vect(y))^2.
    $
    This gives~@eq:strong_convex_eq1 after rearranging.

  - @eq:strong_convex_eq1 $=>$~@eq:strong_convexity.
    To prove this implication, suppose that~@eq:strong_convex_eq1 holds,
    take $(vect(x), vect(y)) in real^n times real^n$
    and let $vect(z) = theta vect(x) + (1 - theta) vect(y)$.
    Using~@eq:strong_convex_eq1 successively with $(vect(x), vect(z))$ and~$(vect(y), vect(z))$,
    we deduce
    $
      J(vect(x)) & >= J(vect(z)) + ipbig(nabla J(vect(z)), vect(x) - vect(z)) + alpha/2 norm(vect(x) - vect(z))^2, \
      J(vect(y)) & >= J(vect(z)) + ipbig(nabla J(vect(z)), vect(y) - vect(z)) + alpha/2 norm(vect(y) - vect(z))^2.
    $
    Combining these inequalities,
    we deduce that
    $
      theta J(vect(x)) + (1 - theta) J(vect(y))
      & >= J(vect(z)) + ipBig(nabla J(vect(z)), theta vect(x) + (1 - theta) vect(y) - vect(z)) \
      & wide + (alpha theta)/2 norm(vect(x) - vect(z))^2 + (alpha(1 - theta))/2 norm(vect(y) - vect(z))^2 \
      & = J(vect(z)) + 0 + alpha/2 theta (1-theta) norm(vect(x) - vect(y))^2.
    $
    Rearranging gives~@eq:strong_convexity.

  - @eq:strong_convex_eq1 $=>$~@eq:strong_convex_eq2.
    Assuming that~@eq:strong_convex_eq1 holds and applying this inequality first to~$(vect(x), vect(y))$ and then to~$(vect(y), vect(x))$,
    we obtain
    $
      J(vect(x)) & >= J(vect(y)) + ipbig(nabla J(vect(y)), vect(x) - vect(y)) + alpha/2 norm(vect(x) - vect(y))^2 \
      J(vect(y)) & >= J(vect(x)) + ipbig(nabla J(vect(x)), vect(y) - vect(x)) + alpha/2 norm(vect(x) - vect(y))^2.
    $
    Adding these equations and rearranging,
    we deduce~@eq:strong_convex_eq2.

  - @eq:strong_convex_eq2 $=>$~@eq:strong_convex_eq1.
    Suppose that~@eq:strong_convex_eq2 holds and
    take $(vect(x), vect(y)) in real^n times real^n$.
    Using the fundamental theorem of analysis and~@eq:strong_convex_eq2,
    we have
    $
      J(vect(x))
      & = J(vect(y)) + integral_0^1 ipBig(nabla J big((vect(y) + theta (vect(x) - vect(y)))), vect(x) - vect(y)) thin dd theta \
      & >= J(vect(y)) + integral_0^1 ipBig(nabla J(vect(y)), vect(x) - vect(y)) + alpha theta norm(vect(x) - vect(y))^2 thin dd theta \
      & = J(vect(y)) + ipBig(nabla J(vect(y)), vect(x) - vect(y)) + alpha/2 norm(vect(x) - vect(y))^2,
    $
    which gives~@eq:strong_convex_eq1.
  We have proved all the implications required to conclude the proof.
]

== Unconstrained optimization <sec:unconstrained_optimization>
Throughout this section $cal(K) = real^n$.
We begin by establishing conditions under which the optimization problem~@eq:optimization admits a unique solution in this setting.
We first prove existence of a global minimizer under appropriate conditions.
#proposition(title: [Existence of a global minimizer])[
  Suppose that $J colon real^n -> real$ is continuous and coercive,
  the latter meaning that $J(vect(x)) -> infinity$ when $norm(vect(x)) -> infinity$.
  Then there exists a global minimizer of~$J$ in $real^n$.
] <proposition:existence_minimizer>
#proof[
  Let $(vect(x)_n)_(n in nat)$ be a minimizing sequence of~$J$,
  i.e.~a sequence in $real^n$ such that
  $ J(vect(x)_n) -> inf_(vect(x) in real^n) J(vect(x)) quad " as " n -> infinity. $
  The sequence~$(vect(x)_n)$ is bounded,
  because otherwise it would hold that $J(vect(x)_n) -> infinity$ by coercivity.
  Therefore, since closed bounded sets in~$real^n$ are compact,
  there is a subsequence $(vect(x)_(n_k))_(k in nat)$ converging to some $vect(x)_* in real^n$.
  Since $J$ is continuous, we have that
  $ J(vect(x)_*) = lim_(k -> infinity) J(vect(x)_(n_k)) = inf_(vect(x) in real^n) J(vect(x)). $
  We conclude that $vect(x)_*$ is a minimizer of $J$.
]
#remark[
  We relied crucially in the proof of~@proposition:existence_minimizer on the fact that closed bounded sets in $real^n$ are compact.
  In the infinite-dimensional setting,
  coercivity and continuity alone are not sufficient to guarantee the existence of a minimizer.
]

#indent Uniqueness of the minimizer can be established under a strict convexity assumption.

#proposition(title: [Uniqueness of the minimizer])[
  If $J$ is strictly convex,
  then there exists at most one global minimizer.
] <proposition:uniqueness_minimizer>
#proof[
  Suppose for contradiction that there were two minimizers~$vect(x)_*$ and~$vect(y)_*$.
  Then by strict convexity we have
  $
    J lr(((vect(x)_* + vect(y)_*)/2))
    < 1/2 big((J(vect(x)_*) + J(vect(y)_*))) = J(vect(x)_*),
  $
  which contradicts the minimality of~$J(vect(x)_*)$.
]

#indent Finally, before introducing the steepest descent algorithm,
we recall the following standard result from analysis,
the proof of which is left as an exercise.

#theorem(title: [Euler condition])[
  Suppose that $J colon real^n -> real$ is differentiable.
  - If $vect(x)_*$ is a local minimizer of~$J$,
    then
    $nabla J(vect(x)_*) = 0.$

  - If $J$ is convex,
    then~$nabla J(vect(x)_*) = 0$ if and only if~$vect(x)_*$ is a global minimizer.
] <theorem:euler>

#paragraph[Steepest descent method.][
  In this section, we study the more general version of the steepest descent with _fixed step_ given in~@algo:steepest_descent_method_optimization.
]
#algorithm(
  caption: [Steepest descent method],
  (0, [Pick $lambda$, and initial $vect(x)_0$.]),
  (0, [*for* $k in {0, 1, dots}$ *do*]),
  (1, [$vect(d)_k <- nabla J(vect(x)_k)$]),
  (1, [$vect(x)_(k+1) <- vect(x)_k - lambda vect(d)_k$]),
  (0, [*end for*]),
) <algo:steepest_descent_method_optimization>

#remark[
  We encountered the steepest descent with fixed step for a quadratic objective function when we
  analyzed Richardson's method for solving linear equations in~@cha:solution_of_linear_systems.
]

#indent In practice, @algo:steepest_descent_method_optimization must be supplemented with an appropriate stopping criterion.
This could be, for example, a criterion of the form $norm(vect(x)_(k+1) - vect(x)_k) <= epsilon$,
or $big(|J(vect(x)_(k+1)) - J(vect(x)_k)|) <= epsilon$.
It is sometimes also useful to use a normalized criterion of the form $norm(vect(x)_(k+1) - vect(x)_k) <= epsilon norm(vect(x)_0).$
The steepest descent method may be viewed as a fixed point iteration for the function
#numbered($
  vect(F)_lambda (vect(x)) = vect(x) - lambda nabla J(vect(x)).
$) <eq:fixed_point_gradient_descent>
A point $vect(x)_* in real^n$ is a fixed point of this function if and only if~$vect(x)_*$ is a solution to the nonlinear equation $nabla J(vect(x)_*) = 0$.
We shall now prove the convergence of the steepest descent under appropriate assumptions on the function~$J$.
#theorem(title: [Convergence of the steepest descent method])[
  Suppose that $J$ is differentiable, strongly convex with parameter~$alpha$,
  and that its gradient $nabla J colon real^n -> real^n$ is Lipschitz continuous with parameter~$L$:
  #numbered($
    forall (vect(x), vect(y)) in real^n times real^n, wide
    norm(nabla J(vect(x)) - nabla J(vect(y))) <= L norm(vect(x) - vect(y)).
  $) <eq:lipschitz_gradient>
  Then provided that
  #numbered($
    0 < lambda < (2 alpha)/L,
  $) <eq:condition_time_step>
  the steepest descent method with fixed step is convergent.
  More precisely,
  there exists $rho in (0, 1)$ such that for all~$k >= 0$
  #numbered($
    norm(vect(x)_k - vect(x)_*) <= rho^k norm(vect(x)_0 - vect(x)_*).
  $) <eq:exponential_convergence>
] <theorem:convergence_steepest_descent_optimization>
#proof[
  Under the assumptions of the theorem,
  there exists a unique global minimizer of~$J$,
  which is the unique fixed point of~$vect(F)_lambda$.
  We begin by proving that $vect(F)_lambda$ defined in~@eq:fixed_point_gradient_descent is globally Lipschitz continuous.
  We have
  $
    norm(vect(F)_lambda (vect(x)) - vect(F)_lambda (vect(y)))^2
    & = lr(||vect(x) - vect(y) - lambda big((nabla J(vect(x)) - nabla J(vect(y))))||)^2 \
    & = norm(vect(x) - vect(y))^2 - 2 lambda ipbig(vect(x) - vect(y), nabla J(vect(x)) - nabla J(vect(y))) + lambda^2 norm(nabla J(vect(x)) - nabla J(vect(y)))^2 \
    & <= (1 - 2 alpha lambda + lambda^2 L) norm(vect(x) - vect(y))^2,
  $
  where we employed~@eq:strong_convex_eq2 for the second term and~@eq:lipschitz_gradient for the third term.
  Thus, $vect(F)_lambda$ is globally Lipschitz continuous with constant $rho = sqrt(1 - 2 alpha lambda + lambda^2 L)$,
  which is less than 1 if and only~@eq:condition_time_step is satisfied.
  The bound~@eq:exponential_convergence then follows by noting that
  $
    norm(vect(x)_k - vect(x)_*) = norm(vect(F)_lambda (vect(x)_(k-1)) - vect(F)_lambda (vect(x)_*))
    <= rho norm(vect(x)_(k-1) - vect(x)_*) <= dots <= rho^k norm(vect(x)_0 - vect(x)_*),
  $
  which concludes the proof.
  (Note that~@eq:exponential_convergence also follows from~@theorem:exponential_convergence_fixed_point.)
]

#remark(title: [Convergence speed])[
  The choice of~$lambda$ minimizing the Lipschitz constant~$rho$ is given by $lambda_* = alpha/L^2$,
  which corresponds to $rho_* = 1 - lr((alpha/L))^2$.
  Often, in practice, it holds that~$alpha << L$,
  in which case the convergence of the steepest descent with fixed step is slow.
]

== Constrained optimization <sec:constrained_optimization>

In this section,
we assume that $cal(K) subset real^n$.
We begin by establishing well-posedness of the optimization problem~@eq:optimization in this setting.

#proposition(title: [Well posedness of~@eq:optimization in the constrained setting])[
  The two items below concern existence and uniqueness, respectively.
  - Suppose that~$cal(K) subset real^n$ is _closed_ and that $J colon cal(K) -> real$ is continuous and coercive.
    Then there exists a global minimizer of~$J$ in $cal(K)$.

  - Suppose that $cal(K) subset real^n$ is _convex_ and that $J colon cal(K) -> real$ is strictly convex.
    Then there exists at most one global minimizer.
] <proposition:well-posedness_constrained>
#proof[
  The proof is very similar to those of~@proposition:existence_minimizer and~@proposition:uniqueness_minimizer,
  and so we leave it to the reader.
  Note that the set~$cal(K)$ must be closed to ensure existence,
  and convex to guarantee uniqueness.
  These assumptions are clearly satisfied when $cal(K) = real^n$,
  so~@proposition:well-posedness_constrained indeed generalizes @proposition:existence_minimizer and @proposition:uniqueness_minimizer.
]

#indent The following theorem,
which generalizes~@theorem:euler,
establishes a characterization of the minimizer when $J$ is differentiable.
#theorem(title: [Euler--Lagrange conditions])[
  Suppose that $J colon cal(K) -> real$ is differentiable and
  that $cal(K) subset real^n$ is closed and convex.
  Then the following statements hold.
  - If $vect(x)_*$ is a local minimizer of~$J$,
    then
    #numbered($
      forall vect(x) in cal(K),
      wide ip(nabla J(vect(x)_*), vect(x) - vect(x)_*) >= 0.
    $) <eq:euler_lagrange>

  - Conversely, if~@eq:euler_lagrange is satisfied and~$J$ is convex,
    then~$vect(x)_*$ is a global minimizer of~$J$.
] <theorem:euler_lagrange>
#proof[
  Suppose that $vect(x)_*$ is a local minimizer of $J$.
  This means that there exists~$delta > 0$ such that
  $
    forall vect(x) in B_delta (vect(x)_*) inter cal(K), wide
    J(vect(x)_*) <= J(vect(x)).
  $
  Therefore $J(vect(x)_*) <= J big(((1-t)vect(x)_* + t vect(x)))$ for all $t in [0, 1]$ sufficiently small.
  But then
  $
    ip(nabla J(vect(x)_*), vect(x) - vect(x)_*)
    = lim_(t -> 0) (J big(((1-t)vect(x)_* + t vect(x))) - J(vect(x)_*))/t >= 0.
  $
  Conversely, suppose that~@eq:euler_lagrange is satisfied and that~$J$ is convex.
  Since~$J$ is convex, equation~@eq:strong_convex_eq1 holds with $alpha = 0$,
  and applying this equation with~$vect(y) = vect(x)_*$,
  we deduce that $vect(x)_*$ is a global minimizer.
]

#indent The steepest descent~@algo:steepest_descent_method_optimization can be extended to optimization problems with constraints by introducing an additional projection step.
In order to precisely formulate the algorithm,
we begin by introducing the projection operator~$Pi_cal(K)$.

#proposition(title: [Projection on a closed convex set])[
  Suppose that $cal(K)$ is a closed convex subset of~$real^n$.
  Then for all $vect(x) in real^n$ there a unique $Pi_cal(K) vect(x) in cal(K)$,
  called the orthogonal projection of~$vect(x)$ onto~$cal(K)$,
  such that
  $ norm(Pi_cal(K) vect(x) - vect(x)) = inf_(vect(y) in cal(K)) norm(vect(y) - vect(x)). $
]
#proof[
  The functional $J_vect(x) (vect(y)) = norm(vect(y) - vect(x))^2$ is strongly convex,
  and so~@proposition:well-posedness_constrained immediately implies the existence and uniqueness of~$Pi_cal(K) vect(x)$.
]

#remark[
  In view of~@theorem:euler_lagrange,
  the projection~$Pi_cal(K) vect(x)$ is the unique element of~$cal(K)$ which satisfies
  #numbered($
    forall vect(y) in cal(K),
    wide ip(Pi_cal(K) vect(x) - vect(x), vect(y) - Pi_cal(K) vect(x)) >= 0.
  $) <eq:euler_lagrange_projection>
]

#indent We are now ready to present the steepest descent method with projection:
see~@algo:steepest_descent_method_optimization_projection.
Like~@algo:steepest_descent_method_optimization,
the steepest descent with projection may be viewed as a fixed point iteration,
this time for the function
#numbered($
  vect(F)_lambda (vect(x)) := Pi_cal(K) big((vect(x) - lambda nabla J(vect(x)))).
$) <eq:constrained_fixed_point>
We now prove the convergence of the method.
#algorithm(
  caption: [Steepest descent with projection],
  (0, [Pick $lambda$, and initial $vect(x)_0$.]),
  (0, [*for* $k in {0, 1, dots}$ *do*]),
  (1, [$vect(d)_k <- nabla J(vect(x)_k)$]),
  (1, [$vect(x)_(k+1) <- Pi_cal(K) (vect(x)_k - lambda vect(d)_k)$]),
  (0, [*end for*]),
) <algo:steepest_descent_method_optimization_projection>

#theorem(title: [Convergence of steepest descent with projection])[
  Suppose that $J$ is differentiable, strongly convex with parameter~$alpha$,
  and that its gradient $nabla J colon real^n -> real^n$ is Lipschitz continuous with parameter~$L$.
  Assume also that $cal(K) subset real^n$ is closed and convex.
  Then provided that
  #numbered($
    0 < lambda < (2 alpha)/L,
  $) <eq:condition_time_step_constroined>
  the steepest descent method with fixed step is convergent.
  More precisely,
  there exists $rho in (0, 1)$ such that for all~$k >= 0$
  $ norm(vect(x)_k - vect(x)_*) <= rho^k norm(vect(x)_0 - vect(x)_*). $
] <theorem:convergence_steepest_descent_optimization_constraint>
#proof[
  Under the assumptions,
  there exists a unique global minimizer $vect(x)_* in cal(K)$.
  We already showed in the proof of~@theorem:convergence_steepest_descent_optimization that the mapping $vect(x) |-> vect(x) - lambda nabla J(vect(x))$ is a contraction if and only if~$lambda$ satisfies~@eq:condition_time_step_constroined.
  In order to prove that~$vect(F)_lambda$ given in~@eq:constrained_fixed_point is a contraction under the same condition,
  it is sufficient to prove that $Pi_cal(K) colon real^n -> cal(K)$ satisfies the following estimate:
  $
    forall (vect(x), vect(y)) in real^n times real^n, wide
    norm(Pi_cal(K) vect(x) - Pi_cal(K) vect(y)) <= norm(vect(x) - vect(y)).
  $
  To this end, take $(vect(x), vect(y)) in real^n times real^n$ and let $vect(delta) = Pi_cal(K) vect(x) - Pi_cal(K) vect(y)$.
  By~@eq:euler_lagrange_projection,
  it holds that
  $
    norm(vect(delta))^2
    & = ip(vect(delta), Pi_cal(K) vect(x) - vect(x)) + ip(vect(delta), vect(x) - vect(y)) + ip(vect(delta), vect(y) - Pi_cal(K) vect(y)) \
    & <= 0 + ip(vect(delta), vect(x) - vect(y)) + 0 <= norm(vect(delta)) norm(vect(x) - vect(y)),
  $
  which yields the required inequality.
  Therefore $vect(F)_lambda$ in~@eq:constrained_fixed_point is a contraction and so,
  by the Banach fixed point theorem, it admits a unique fixed point~$vect(y)_* in cal(K)$.
  To show that $vect(y)_* = vect(x)_*$, note that if~$vect(F)_lambda (vect(y)_*) = vect(y)_*$,
  then by~@eq:euler_lagrange_projection it holds that
  $
    forall vect(y) in cal(K),
    wide ip(lambda nabla J(vect(y)_*), vect(y) - vect(y)_*) >= 0.
  $
  Therefore, using~@theorem:euler_lagrange, we obtain that $vect(y)_*$ is a global minimizer of~$J$,
  so~$vect(y)_* = vect(x)_*$.
]

#remark[
  The applicability of~@algo:steepest_descent_method_optimization_projection is limited in practice,
  as computing~$Pi_cal(K) (vect(x))$ analytically is possible only in simple settings.
]

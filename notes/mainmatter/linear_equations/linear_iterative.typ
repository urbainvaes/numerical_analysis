#import "../../template.typ": *

// \ddots and \vdots (`dots` is redefined in the template)
#let ddots = sym.dots.down
#let vdots = sym.dots.v
// \dotsb before a binary operator: the template `cdots` loses the space after it
#let cdotsb = sym.dots.c
// \eqref applied to a non-equation label: the number in parentheses, e.g. "(4.19)"
#let eqref-number(lbl) = context {
  let el = query(lbl).first()
  let loc = el.location()
  link(loc, text(fill: darkblue)[(#chapter-number(loc).#el.counter.at(loc).first())])
}

== Iterative methods for linear systems <sec:iterative_methods>
Iterative methods enjoy more flexibility than direct methods,
because they can be stopped at any point if the residual is deemed sufficiently small.
This generally enables to obtain a good solution at a computational cost that is significantly lower than that of direct methods.
In this section,
we present and study two classes of iterative methods:
basic iterative methods based on a splitting of the matrix $mat(A)$,
and the so-called Krylov subspace methods.

=== Basic iterative methods <sub:basic_iterative_methods>
The basic iterative methods are particular cases of a general _splitting method_.
Given a splitting of the matrix of the linear system as $mat(A) = mat(M) - mat(N)$,
for a nonsingular matrix $mat(M) in complex^(n times n)$ and a matrix $mat(N) in complex^(n times n)$,
together with an initial guess $vect(x)^((0))$ of the solution,
one step of this general method reads
#numbered($ mat(M) vect(x)^((k+1)) = mat(N) vect(x)^((k)) + vect(b). $) <eq:linear_system_iterative_method>
For any choice of splitting,
the exact solution $vect(x)_*$ to the linear system is a fixed point of this iteration,
in the sense that if $vect(x)^((0)) = vect(x)_*$, then $vect(x)^((k)) = vect(x)_*$ for all $k >= 0$.
Equation~@eq:linear_system_iterative_method is a linear system with matrix~$mat(M)$,
unknown $vect(x)^((k+1))$, and right-hand side $mat(N) vect(x)^((k)) + vect(b)$.
There is a compromise between the cost of a single step and the speed of convergence of the method.
In the extreme case where $mat(M) = mat(A)$ and $mat(N) = 0$,
the method converges to the exact solution in one step,
but performing this step amounts to solving the initial problem.
In practice, in order for the method to be useful,
the linear system~@eq:linear_system_iterative_method should be relatively simple to solve.
Concretely, this means that the matrix~$mat(M)$ should be diagonal, triangular, block diagonal, or block triangular.
The error $vect(e)^((k))$ and residual $vect(r)^((k))$ at iteration $k$ are defined as follows:
$ vect(e)^((k)) = vect(x)^((k)) - vect(x)_*, wide vect(r)^((k)) = mat(A) vect(x)^((k)) - vect(b). $

==== Convergence of the splitting method <ssub:convergence_of_the_basic_splitting_method>
Before presenting concrete examples of splitting methods,
we obtain a necessary and sufficient condition for the convergence of~@eq:linear_system_iterative_method
for any initial guess~$vect(x)^((0))$.

#proposition(title: [Convergence])[
  The splitting method~@eq:linear_system_iterative_method converges for any initial~$vect(x)^((0))$
  if and only if $rho(mat(M)^(-1) mat(N)) < 1$.
  In addition, for any $epsilon > 0$ there exists $K > 0$ such that
  #numbered($
    forall k >= K, wide
    norm(vect(e)^((k))) <= big((rho(mat(M)^(-1) mat(N)) + epsilon))^k norm(vect(e)^((0))).
  $) <eq:additional_claim>
] <proposition:linear_convergence>
#proof[
  Let $vect(x)_*$ denote the solution to the linear system.
  Since $mat(M) vect(x)_* - mat(N) vect(x)_* = vect(b)$,
  we have
  $ mat(M) (vect(x)^((k+1)) - vect(x)_*) = mat(N) (vect(x)^((k)) - vect(x)_*). $
  Using the assumption that $mat(M)$ is nonsingular,
  we obtain that the error satisfies the equation
  $ vect(e)^((k+1)) = (mat(M)^(-1) mat(N)) vect(e)^((k)). $
  Applying this equality repeatedly,
  we deduce that
  #numbered($ vect(e)^((k)) = (mat(M)^(-1) mat(N)) vect(e)^((k-1)) = dots = (mat(M)^(-1) mat(N))^k vect(e)^((0)). $) <eq:error_equation_basic_iterative_method>

  #paragraph[Proof that $rho(mat(M)^(-1) mat(N)) < 1$ is necessary for convergence.][
    We prove the equivalent claim that if $rho(mat(M)^(-1) mat(N)) >= 1$,
    then there exists $vect(x)^((0))$ such that the method is not convergent.
    Indeed, assume that $vect(x)^((0)) = vect(x)_* + vect(v)^((0))$,
    where $vect(v)^((0))$ is the eigenvector of $mat(M)^(-1) mat(N)$ associated with the eigenvalue of largest modulus.
    Then $vect(e)^((0)) = vect(v)^((0))$ and the right-hand side of~@eq:error_equation_basic_iterative_method does not converge to 0 in the limit as $k -> 0$,
    because
    $ norm((mat(M)^(-1) mat(N))^k vect(e)^((0))) = rho(mat(M)^(-1) mat(N))^k norm(vect(v)^((0))) >= norm(vect(v)^((0))). $
    Thus, the condition $rho(mat(M)^(-1) mat(N)) < 1$ is necessary to ensure convergence for all initial guess $vect(x)^((0))$.
  ]

  #paragraph[Proof that $rho(mat(M)^(-1) mat(N)) < 1$ is sufficient for convergence.][
    In order to show that the condition is also sufficient,
    note that by~@eq:error_equation_basic_iterative_method
    $
      forall k >= 0, wide
      norm(vect(e)^((k)))
      <= big(||(mat(M)^(-1) mat(N))^k||) norm(vect(e)^((0))).
    $
    By Gelfand's formula,
    proved in~@proposition:matrices_gelfands of @cha:vectors_and_matrices,
    it holds that
    #numbered($
      lim_(k -> infinity) big(||(mat(M)^(-1) mat(N))^k||)^(1/k) = rho(mat(M)^(-1) mat(N)).
    $) <eq:gelfand>
    Therefore, we deduce that if $rho(mat(M)^(-1) mat(N)) < 1$,
    then $big(||(mat(M)^(-1) mat(N))^k||) -> 0$ and so $vect(e)^((k)) -> vect(0)$.
    In addition, it follows from~@eq:gelfand that for all $epsilon > 0$ there is $K in nat$ such that
    $
      forall k >= K, wide
      big(||(mat(M)^(-1) mat(N))^k||)^(1/k) <= rho(mat(M)^(-1) mat(N)) + epsilon.
    $
    Rearranging this inequality gives~@eq:additional_claim.
  ]
]

#indent At this point,
it is natural to wonder whether there exist sufficient conditions on the matrix~$mat(A)$ such that the inequality $rho(mat(M)^(-1) mat(N)) < 1$ is satisfied,
which is best achieved on a case by case basis.
In the next sections,
we present four instances of splitting methods.
For each of them,
we obtain a sufficient condition for convergence.
We are particularly interested in the case where the matrix~$mat(A)$ is Hermitian and positive definite,
which often arises in applications,
and in the case where $mat(A)$ is strictly row or column diagonally dominant.
We recall that a matrix $mat(A)$ is said to be row or column diagonally dominant
if, respectively,
$
  |a_(i i)| >= sum_(j != i) |a_(i j)| quad forall i
  wide " or " wide
  |a_(j j)| >= sum_(i != j) |a_(i j)| quad forall j.
$

==== Richardson's method <ssub:richardson_s_method>

Arguably the simplest splitting of the matrix $mat(A)$ is given by $mat(A) = 1/omega mat(I) - big((1/omega mat(I) - mat(A)))$,
for some parameter~$omega in real$,
which leads to _Richardson's method_:
#numbered($ vect(x)^((k+1)) = vect(x)^((k)) + omega (vect(b) - mat(A) vect(x)^((k))). $) <eq:linear_richardson>
In this case the spectral radius which enters in the asymptotic rate of convergence is given by
$
  rho(mat(M)^(-1) mat(N))
  = rho lr((omega lr((1/omega mat(I) - mat(A)))))
  = rho big((mat(I) - omega mat(A)))
$
The eigenvalues of the matrix $mat(I) - omega mat(A)$ are given by $1 - omega lambda_i$,
where $(lambda_i)_(1 <= i <= L)$ are the eigenvalues of~$mat(A)$.
Therefore, the spectral radius is given by
$
  rho(mat(M)^(-1) mat(N))
  = max_(1 <= i <= L) abs(1 - omega lambda_i).
$
If the eigenvalues of the matrix~$mat(A)$ do not either (i) all have a positive real part or (ii) all have a negative real part,
then
$
  forall omega in real,
  wide max_(1 <= i <= L) abs(1 - omega lambda_i) >= 1.
$
In other words,
by~@proposition:linear_convergence,
there is for any choice of~$omega in real$ some $vect(x)^((0))$ such that Richardson's method is non-convergent.
Therefore, in order for convergence to hold for all~$vect(x)^((0))$,
it is necessary that the eigenvalues of $mat(A)$ either all have a negative real part,
or all have a positive real part.
We focus in the next paragraph on the latter case and we also assume,
for simplicity, that $mat(A)$ is Hermitian.

#paragraph[Case of symmetric positive definite $mat(A)$.][
  If the matrix $mat(A)$ is Hermitian and positive definite,
  the eigenvalues of $mat(A)$ are all real and positive,
  and it is possible to explicitly calculate the optimal value of $omega$ for convergence.
  In order for convergence to be as fast as possible,
  the spectral radius of $mat(M)^(-1) mat(N)$ should be as small as possible,
  in view of~@proposition:linear_convergence.
  Denoting by $lambda_min$ and $lambda_max$ the minimum and maximum eigenvalues of $mat(A)$,
  it is not difficult to show that
  #numbered($
    rho(mat(M)^(-1) mat(N))
    = max_(1 <= i <= L) abs(1 - omega lambda_i)
    = max big(\{ abs(1 - omega lambda_min)\, abs(1 - omega lambda_max) \}).
  $) <eq:spectral_richardson>
  The right-hand side is minimized $1 - omega lambda_min = omega lambda_max - 1$,
  in which case the two arguments of the maximum in~@eq:spectral_richardson are equal.
  From this we deduce the optimal value of $omega$ and the associated spectral radius:
  $
    omega_"opt" = 2/(lambda_max + lambda_min),
    wide
    rho_"opt"
    = 1 - (2 lambda_min)/(lambda_max + lambda_min)
    = (lambda_max - lambda_min)/(lambda_max + lambda_min)
    = (kappa_2(mat(A)) - 1)/(kappa_2(mat(A)) + 1).
  $
  We observe that the smaller the condition number of the matrix $mat(A)$,
  the better the asymptotic rate of convergence.
] <par:case_of_symmetric_positive_definite_mat_a_>

#remark(title: [Link to optimization])[
  In the case where $mat(A) in real^(n times n)$ is symmetric and positive definite,
  the Richardson update~@eq:linear_richardson may be viewed as a step of the steepest descent algorithm,
  // which we study carefully in~\cref{cha:Optimization},
  for the function $f(vect(x)) = 1/2 vect(x)^T mat(A) vect(x) - vect(b)^T vect(x)$:
  #numbered($ vect(x)^((k+1)) = vect(x)^((k)) - omega nabla f(vect(x)^((k))). $) <eq:linear_richardson_gradient>
  The gradient of this function is $nabla f(vect(x)) = mat(A) vect(x) - vect(b)$,
  and its Hessian matrix is $mat(A)$.
  Since the Hessian matrix is positive definite, the function is convex
  and attains its global minimum when $nabla f$ is zero,
  i.e.~when $mat(A) x = vect(b)$.
] <remark:linear_link_optimization>

==== Jacobi's method <ssub:jacobi_s_method>
In Jacobi's method, the matrix $mat(M)$ in the splitting is the diagonal matrix $mat(D)$ with the same entries as those of $mat(A)$ on the diagonal.
We denote by $mat(L)$ and $mat(U)$ the lower and upper triangular parts of~$mat(A)$,
without the diagonal.
One step of the method reads
#numbered($
  mat(D) vect(x)^((k+1))
  = (mat(D) - mat(A)) vect(x)^((k)) + vect(b)
  = - (mat(L) + mat(U)) vect(x)^((k)) + vect(b)
$) <eq:linear_jacobi>
Since the matrix $mat(D)$ on the left-hand side is diagonal,
this linear system with unknown $vect(x)^((k+1))$ is very simple to solve.
The equation~@eq:linear_jacobi can be rewritten as
$
  cases(gap: #0.5em,
    a_11 x^((#text(fill: darkblue, $k+1$)))_1 + a_12 x^((k))_2 + cdotsb + a_(1n) x^((k))_n = b_1,
    a_21 x^((k))_1 + a_22 x^((#text(fill: darkblue, $k+1$)))_2 + cdotsb + a_(2n) x^((k))_n = b_2,
    vdots,
    a_(n 1) x^((k))_1 + a_(n 2) x^((k))_2 + cdotsb + a_(n n) x^((#text(fill: darkblue, $k+1$)))_n = b_n.
  )
$
The updates for each of the entries of $vect(x)^((k+1))$ are independent,
and so the Jacobi method lends itself well to parallel implementation.
The computational cost of one iteration,
measured in number of floating point operations required,
scales as $bigo(n^2)$ if~$mat(A)$ is a full matrix,
or $bigo(n k)$ if~$mat(A)$ is a sparse matrix with $k$ nonzero elements per row on average.
It is simple to prove the convergence of Jacobi's method is the case where $mat(A)$ is diagonally dominant.
#proposition[
  Assume that $mat(A)$ is strictly (row or column) diagonally dominant.
  Then it holds that $rho(mat(M)^(-1) mat(N)) < 1$ for the Jacobi splitting.
] <proposition:linear_convergence_jacobi>
#proof[
  Assume that $lambda$ is an eigenvalue of $mat(M)^(-1) mat(N)$
  and $vect(v)$ is the associated unit eigenvector.
  Then
  $
    mat(M)^(-1) mat(N) vect(v) = lambda vect(v)
    quad <=> quad
    mat(N) vect(v) = lambda mat(M) vect(v)
    quad <=> quad
    (mat(N) - lambda mat(M)) vect(v) = 0.
  $
  In the case of Jacobi's splitting,
  this is equivalent to
  $ - (mat(L) + lambda mat(D) + mat(U)) vect(v) = 0. $
  If $abs(lambda) > 1$,
  then the matrix on the left-hand side of this equation is diagonally dominant and thus invertible
  (see @exercise:invertibility_diagonal_dominant).
  Therefore $vect(v) = 0$, but this is a contradiction because $vect(v)$ is vector of unit norm.
  Consequently, all the eigenvalues are bounded from above strictly by 1 in modulus.
]

==== Gauss--Seidel's method <ssub:gauss_seidel_s_method>

In Gauss Seidel's method, the matrix $mat(M)$ in the splitting is the lower triangular part of $mat(A)$,
including the diagonal.
One step of the method then reads
#numbered($ (mat(L) + mat(D)) vect(x)^((k+1)) = - mat(U) vect(x)^((k)) + vect(b) $) <eq:linear_gauss_seidel>
The system is solved by forward substitution.
The equation~@eq:linear_gauss_seidel can be rewritten equivalently as
$
  cases(gap: #0.5em,
    a_11 x^((#text(fill: darkblue, $k+1$)))_1 + a_12 x^((k))_2 + a_13 x^((k))_3 + cdotsb + a_(1n) x^((k))_n = b_1,
    a_21 x^((#text(fill: darkblue, $k+1$)))_1 + a_22 x^((#text(fill: darkblue, $k+1$)))_2 + a_23 x^((k))_3 + cdotsb + a_(2n) x^((k))_n = b_2,
    a_32 x^((#text(fill: darkblue, $k+1$)))_1 + a_32 x^((#text(fill: darkblue, $k+1$)))_2 + a_33 x^((#text(fill: darkblue, $k+1$)))_3 + cdotsb + a_(3n) x^((k))_n = b_3,
    vdots,
    a_(n 1) x^((#text(fill: darkblue, $k+1$)))_1 + a_(n 2) x^((#text(fill: darkblue, $k+1$)))_2 + a_(n 3) x^((#text(fill: darkblue, $k+1$)))_3 + cdotsb + a_(n n) x^((#text(fill: darkblue, $k+1$)))_n = b_n.
  )
$
Given $vect(x)^((k))$,
the first entry of $vect(x)^((k+1))$ is obtained from the first equation.
Then the value of the second entry is obtained from the second equation, etc.
Unlike Jacobi's method,
the Gauss--Seidel method is sequential and the entries of $vect(x)^((k+1))$ cannot be updated in parallel.

It is possible to prove the convergence of the Gauss--Seidel method in particular cases.
For example, the method converges if $mat(A)$ is strictly diagonally dominant.
Proving this,
using an approach similar to that in the proof of~@proposition:linear_convergence_jacobi,
is the goal of~@exercise:linear_convergence_gauss_seidel.
It is also possible to prove convergence when $mat(A)$ is Hermitian and positive definite.
We show this in the next section for the relaxation method,
which generalizes the Gauss--Seidel method.

==== Relaxation method

The relaxation method generalizes the Gauss--Seidel method.
It corresponds to the splitting
#numbered($
  mat(A) = lr(((mat(D))/omega + mat(L))) - lr(((1 - omega)/omega mat(D) - mat(U))).
$)
When $omega = 1$,
this is simply the Gauss--Seidel splitting.
The idea of the relaxation method is that,
by letting $omega$ be a parameter that can differ from 1,
faster convergence can be achieved.
This intuition will be verified later.
The equation~@eq:linear_system_iterative_method for this splitting can be rewritten equivalently as
$
  cases(gap: #0.5em,
    a_11 big((x^((#text(fill: darkblue, $k+1$)))_1 - x^((k))_1)) = -omega lr((a_11 x^((k))_1 + a_12 x^((k))_2 + cdotsb + a_(1n) x^((k))_n - b_1)),
    a_22 big((x^((#text(fill: darkblue, $k+1$)))_2 - x^((k))_2)) = -omega lr((a_21 x^((#text(fill: darkblue, $k+1$)))_1 + a_22 x^((k))_2 + cdotsb + a_(2n) x^((k))_n - b_2)),
    vdots,
    a_(n n) big((x^((#text(fill: darkblue, $k+1$)))_n - x^((k))_n)) = -omega lr((a_(n 1) x^((#text(fill: darkblue, $k+1$)))_1 + a_(n 2) x^((#text(fill: darkblue, $k+1$)))_2 + cdotsb + a_(n n) x^((k))_n - b_n)).
  )
$
The coefficient on the right-hand side is larger than in the Gauss--Seidel method if $omega > 1$,
and smaller if $omega < 1$.
These regimes are called _over-relaxation_ and _under-relaxation_, respectively.

To conclude this section,
we establish a sufficient condition for the convergence of the relaxation method,
and also of the Gauss--Seidel method as a particular case when $omega = 1$,
when the matrix $mat(A)$ is Hermitian and positive definite.
To this end, we begin by showing the following preparatory result,
which concerns a general splitting $mat(A) = mat(M) - mat(N)$.
#proposition[
  Let $mat(A)$ be Hermitian and positive definite.
  If the Hermitian matrix $mat(M)^* + mat(N)$ is positive definite,
  then $rho(mat(M)^(-1) mat(N)) < 1$.
] <proposition:criterion_convergence>
#proof[
  First, notice that $mat(M)^* + mat(N)$ is indeed Hermitian because
  $ (mat(M)^* + mat(N))^* = mat(M) + mat(N)^* = mat(A) + mat(N) + mat(N)^*. $
  We will show that $norm(mat(M)^(-1) mat(N))_mat(A) < 1$,
  where $norm(placeholder)_mat(A)$ is the matrix norm induced by the following norm on vectors:
  $ norm(vect(x))_mat(A) := sqrt(vect(x)^* mat(A) vect(x)). $
  Showing that this indeed defines a vector norm is the goal of~@exercise:linear_norm_induced_A.
  Since $mat(N) = mat(M) - mat(A)$, it holds that
  $
    norm(mat(M)^(-1) mat(N))_mat(A)
    = norm(mat(I) - mat(M)^(-1) mat(A))_mat(A),
  $
  and so
  $
    norm(mat(M)^(-1) mat(N))_mat(A)
    = sup big(\{ norm(vect(x) - mat(M)^(-1) mat(A) vect(x))_mat(A) colon norm(vect(x))_mat(A) <= 1 \}).
  $
  Take $vect(x) in complex^n$ with $norm(vect(x))_mat(A) <= 1$
  and let $vect(y) = mat(M)^(-1) mat(A) vect(x)$.
  We calculate
  $
    norm(vect(x) - mat(M)^(-1) mat(A) vect(x))_mat(A)^2
    &= vect(x)^* mat(A) vect(x) - vect(y)^* mat(A) vect(x) - vect(x)^* mat(A) vect(y) + vect(y)^* mat(A) vect(y) \
    &= vect(x)^* mat(A) vect(x) - vect(y)^* mat(M) mat(M)^(-1) mat(A) vect(x) - (mat(M)^(-1) mat(A) vect(x))^* mat(M)^* vect(y) + vect(y)^* mat(A) vect(y) \
    &= vect(x)^* mat(A) vect(x) - vect(y)^* mat(M) vect(y) - vect(y)^* mat(M)^* vect(y) + vect(y)^* (mat(M) - mat(N)) vect(y) \
    &= vect(x)^* mat(A) vect(x) - vect(y)^* (mat(M)^* + mat(N)) vect(y) <= 1 - vect(y)^* (mat(M)^* + mat(N)) vect(y) < 1,
  $
  where we used in the last inequality the assumption that $mat(M)^* + mat(N)$ is positive definite.
  This inequality holds true for all $vect(x) in complex^n$ with $norm(vect(x))_mat(A) = 1$,
  and so we conclude that $norm(mat(M)^(-1) mat(N))_mat(A) < 1$,
  which implies that $rho(mat(M)^(-1) mat(N)) < 1$.
]

#indent As a corollary,
we obtain a sufficient condition for the convergence of the relaxation method.

#corollary[
  Assume that $mat(A)$ is Hermitian and positive definite.
  Then the relaxation method converges if $omega in (0, 2)$.
] <corollary:convergence_relaxation>
#proof[
  For the relaxation method,
  we have
  $
    mat(M) + mat(N)^*
    = lr(((mat(D))/omega + mat(L))) + lr(((1 - omega)/omega mat(D) - mat(U)))^*.
  $
  Since $mat(A)$ is Hermitian,
  it holds that $mat(D)^* = mat(D)$ and $mat(U)^* = mat(L)$.
  Therefore,
  $
    mat(M) + mat(N)^*
    = (2 - omega)/omega mat(D).
  $
  The diagonal elements of $mat(D)$ are all positive,
  because $mat(A)$ is positive definite.
  (Indeed, if there was an index $i$ such that $d_(i i) <= 0$,
  then it would hold that $vect(e)_i^T mat(A) vect(e)_i = d_(i i) <= 0$,
  contradicting the assumption that $mat(A)$ is positive definite.)
  We deduce that $mat(M) + mat(N)^*$ is positive definite if and only if $omega in (0, 2)$.
  We can then conclude the proof by using @proposition:criterion_convergence.
]
Note that @corollary:convergence_relaxation implies as a particular case
the convergence of the Gauss--Seidel method when $mat(A)$ is Hermitian and positive definite.
The condition $omega in (0, 2)$ is in fact necessary for the convergence of the relaxation method,
not only in the case of a Hermitian positive definite matrix $mat(A)$ but in general.

#proposition(title: [Necessary condition for the convergence of the relaxation method])[
  Let $mat(A) in complex^(n times n)$ be an invertible matrix,
  and let $mat(A) = mat(M)_omega - mat(N)_omega$ denote the splitting of the relaxation method with parameter~$omega$.
  It holds that
  $
    forall omega != 0, wide
    rho(mat(M)_omega^(-1) mat(N)_omega)
    >= abs(omega - 1).
  $
] <proposition:necessary_sor>
#proof[
  We recall the following facts:
  - the determinant of a product of matrices is equal to the product of the determinants.
  - the determinant of a triangular matrix is equal to the product of its diagonal entries;
  - the determinant of a matrix is equal to the product of its eigenvalues,
    to the power of their algebraic multiplicity.
    This can be shown from the previous two items,
    by passing to the Jordan normal form.
  Therefore,
  we have that
  $
    det (mat(M)_omega^(-1) mat(N)_omega) = det (mat(M)_omega)^(-1) det (mat(N)_omega)
    = (det lr(((1 - omega)/omega mat(D) - mat(U))))/(det lr(((mat(D))/omega + mat(L))))
    = (1 - omega)^n.
  $
  Since the determinant on the left-hand side is the product of the eigenvalues of $mat(M)_omega^(-1) mat(N)_omega$,
  it is bounded from above in modulus by $rho(mat(M)_omega^(-1) mat(N)_omega)^n$,
  and so we deduce $rho(mat(M)_omega^(-1) mat(N)_omega)^n >= abs(1 - omega)^n$.
  The statement then follows by taking the $n$-th root.
]

==== Comparison between Jacobi and Gauss--Seidel for tridiagonal matrices~#moreinfo <ssub:comparison_between_jacobi_and_gauss_seidel_for_tridiagonal_matrices>
For tridiagonal matrices,
the convergence rate of the Jacobi and Gauss--Seidel methods satisfy an explicit relation,
which we prove in this section.
We denote the Jacobi and Gauss--Seidel splittings by $mat(M)_cal(J) - mat(N)_cal(J)$ and $mat(M)_cal(G) - mat(N)_cal(G)$,
respectively,
and use the following notation for the entries of the matrix $mat(A)$:
$
  pmat(
    a_1, b_1;
    c_1, ddots, ddots;
    , ddots, ddots, b_(n-1);
    , , c_(n-1), a_n
  ).
$
Before presenting and proving the result,
notice that for any $mu != 0$ it holds that
#numbered($
  pmat(
    mu;
    , mu^2;
    , , ddots;
    , , , mu^n
  )
  mat(A)
  pmat(
    mu^(-1);
    , mu^(-2);
    , , ddots;
    , , , mu^(-n)
  )
  =
  pmat(
    a_1, mu^(-1) b_1;
    mu c_1, ddots, ddots;
    , ddots, ddots, mu^(-1) b_(n-1);
    , , mu c_(n-1), a_n
  ).
$) <eq:preliminary_equation>
#proposition[
  Assume that $mat(A)$ is tridiagonal with nonzero diagonal elements,
  so that both $mat(M)_cal(J) = mat(D)$ and $mat(M)_cal(G) = mat(L) + mat(D)$ are invertible.
  Then
  $
    rho(mat(M)_cal(G)^(-1) mat(N)_cal(G))
    = rho(mat(M)_cal(J)^(-1) mat(N)_cal(J))^2
  $
] <proposition:tridiagonal>
#proof[
  If $lambda$ is an eigenvalue of $mat(M)_cal(G)^(-1) mat(N)_cal(G)$ with associated unit eigenvector $vect(v)$,
  then
  $
    mat(M)_cal(G)^(-1) mat(N)_cal(G) vect(v) = lambda vect(v)
    quad <=> quad
    mat(N)_cal(G) vect(v) = lambda mat(M)_cal(G) vect(v)
    quad <=> quad
    (mat(N)_cal(G) - lambda mat(M)_cal(G)) vect(v) = 0.
  $
  For fixed $lambda$,
  there exists a nontrivial solution $vect(v)$ to the last equation if and only if
  $
    p_cal(G) (lambda) := det (mat(N)_cal(G) - lambda mat(M)_cal(G))
    = det (-lambda mat(L) - lambda mat(D) - mat(U)) = 0.
  $
  Likewise, $lambda$ is an eigenvalue of $mat(M)_cal(J)^(-1) mat(N)_cal(J)$ if and only if
  $
    p_cal(J) (lambda) := det (mat(N)_cal(J) - lambda mat(M)_cal(J))
    = det (-mat(L) - lambda mat(D) - mat(U)) = 0.
  $
  Now notice that, for $lambda != 0$,
  $
    p_cal(G) (lambda^2)
    = det big((-lambda^2 mat(L) - lambda^2 mat(D) - mat(U)))
    = lambda^n det big((-lambda mat(L) - lambda mat(D) - lambda^(-1) mat(U))).
  $
  Applying~@eq:preliminary_equation with $mu = lambda != 0$,
  we deduce
  $
    p_cal(G) (lambda^2)
    = lambda^n det big((-mat(L) - lambda mat(D) - mat(U)))
    = lambda^n p_cal(J) (lambda)
  $
  It is clear that this relation is true also if $lambda = 0$.
  Consequently, it holds that if $lambda$ is an eigenvalue of the matrix $mat(M)_cal(J)^(-1) mat(N)_cal(J)$ then
  $lambda^2$ is an eigenvalue of $mat(M)_cal(G)^(-1) mat(N)_cal(G)$.
  Conversely, if $lambda$ is a nonzero eigenvalue of~$mat(M)_cal(G)^(-1) mat(N)_cal(G)$,
  then the two square roots of $lambda$ are eigenvalues of $mat(M)_cal(J)^(-1) mat(N)_cal(J)$.
]

#indent If a matrix $mat(A)$ is tridiagonal and Toeplitz,
i.e.~if it is of the form
$
  pmat(
    a, b;
    c, ddots, ddots;
    , ddots, ddots, b;
    , , c, a
  ),
$
then it is possible to prove that the eigenvalues of~$mat(A)$ are given by
#numbered($
  lambda_k =
  a + 2 sqrt(b c) cos lr(((k pi)/(n+1))),
  wide k = 1, dots, n.
$) <eq:linear_eigenvalues_toeplitz>
In this case, the spectral radius of $mat(M)_cal(J)^(-1) mat(N)_cal(J)$ can be determined explicitly.

==== Monitoring the convergence <ssub:monitoring_the_convergence>
In practice,
we have access to the residual $vect(r)^((k)) = mat(A) vect(x)^((k)) - vect(b)$ at each iteration,
but not to the error $vect(e)^((k)) = vect(x)^((k)) - vect(x)_*$,
as calculating the latter would require to know the exact solution of the problem.
Nevertheless,
the two are related by the equation
$
  vect(r)^((k)) = mat(A) vect(e)^((k))
  quad <=> quad vect(e)^((k)) = mat(A)^(-1) vect(r)^((k)).
$
Therefore, it holds that $norm(vect(e)^((k))) <= norm(mat(A)^(-1)) norm(vect(r)^((k)))$.
Likewise, the relative error satisfies
$
  norm(vect(e)^((k)))/norm(vect(x)_*)
  = norm(mat(A)^(-1) vect(r)^((k)))/norm(mat(A)^(-1) vect(b)),
$
and since $norm(vect(b)) = norm(mat(A) mat(A)^(-1) vect(b)) <= norm(mat(A)) norm(mat(A)^(-1) vect(b))$,
we deduce
$
  norm(vect(e)^((k)))/norm(vect(x)_*)
  <= kappa(mat(A)) norm(vect(r)^((k)))/norm(vect(b)).
$
The fraction on the right-hand side is the _relative residual_.
If the system is well conditioned,
that is if $kappa(A)$ is close to one,
then controlling the relative residual enables a good control of the relative error.

==== Stopping criterion <ssub:stopping_criterion>

In practice,
several criteria can be employed in order to decide when to stop iterating.
Given a small number $epsilon$ (unrelated to the machine epsilon in @cha:rounding_errors),
the following alternatives are available:
- Stop when $norm(vect(r)^((k))) <= epsilon$.
  The downside of this approach is that
  it is not _scaling invariant_:
  when used for solving the following rescaled system
  $ k mat(A) vect(x) = k vect(b), wide k != 1, $
  a splitting method with rescaled initial guess $k vect(x)^((0))$
  will require a number of iterations that depends on $k$:
  fewer if $k << 1$ and more if $k >> 1$.
  In practice, controlling the relative residual and the relative error is often preferable.

- Stop when $norm(vect(r)^((k))) \/ norm(vect(r)^((0))) <= epsilon$.
  This criterion is scaling invariant,
  but the number of iterations is dependent on the quality of the initial guess $vect(x)^((0))$.

- Stop when $norm(vect(r)^((k))) \/ norm(vect(b))$.
  This criterion is generally the best,
  because it is both scaling invariant and the quality of the final iterate is independent of that of the initial guess.

=== The conjugate gradient method <sub:the_conjugate_gradient_method>

As already mentioned in~@remark:linear_link_optimization,
when the matrix $mat(A) in complex^(n times n)$ in the linear system $mat(A) vect(x) = vect(b)$ is symmetric and positive definite,
the system can be interpreted as a minimization problem for the function
#numbered($ f(vect(x)) = 1/2 vect(x)^T mat(A) vect(x) - vect(b)^T vect(x). $) <eq:function_optimization>
The fact that
the exact solution $vect(x)_*$ to the linear system is the unique minimizer of this function appears clearly
when rewriting $f$ as follows:
#numbered($
  f(vect(x))
  = 1/2 (vect(x) - vect(x)_*)^T mat(A) (vect(x) - vect(x)_*) - 1/2 vect(x)_*^T mat(A) vect(x)_*.
$) <eq:rewritten_expression_of_f>
The second term is constant with $vect(x)$,
and the first term is strictly positive if $vect(x) - vect(x)_* != vect(0)$,
because~$mat(A)$ is positive definite.
We saw that Richardson's method can be interpreted as a steepest descent with fixed step size,
$ vect(x)^((k+1)) = vect(x)^((k)) - omega nabla f(vect(x)^((k))). $
In this section,
we will present and study other methods
for solving the linear system~@eq:linear_system
which can be viewed as optimization methods.
Since $mat(A)$ is symmetric,
it is diagonalizable and the function $f$ can be rewritten as
$
  f(vect(x))
  &= 1/2 (vect(x) - vect(x)_*)^T mat(Q) mat(D) mat(Q)^T (vect(x) - vect(x)_*) - 1/2 vect(x)_*^T mat(A) vect(x)_* \
  &= 1/2 (mat(Q)^T vect(e))^T mat(D) (mat(Q)^T vect(e)) - 1/2 vect(x)_*^T mat(A) vect(x)_*, wide vect(e) = vect(x) - vect(x)_*.
$
Therefore, we have that
$
  f(vect(x)) = 1/2 sum_(i=1)^n lambda_i eta_i^2 - 1/2 vect(x)_*^T mat(A) vect(x)_*, wide vect(eta) = vect(Q)^T (vect(x) - vect(x)_*),
$
where $(lambda_i)_(1 <= i <= n)$ are the diagonal entries of $mat(D)$.
This shows that $f$ is a paraboloid
after a change of coordinates.

==== Steepest descent method <ssub:steepest_descent_method>

The steepest descent method is more general than Richardson's method in the sense that
the step size changes from iteration to iteration
and the method is not restricted to quadratic functions of the form~@eq:function_optimization.
Each iteration is of the form
$ vect(x)^((k+1)) = vect(x)^((k)) - omega_k nabla f(vect(x)^((k))). $
It is natural to wonder whether the step size $omega_k$ can be fixed in such a way that
$f(vect(x)^((k+1)))$ is as small as possible.
For the case of the quadratic function~@eq:function_optimization,
this value of $omega_k$ can be calculated explicitly for a general search direction $vect(d)$,
and in particular also when $vect(d) = nabla f(vect(x)^((k)))$.
We calculate that
#eq-last($
  f(vect(x)^((k+1)))
  &= f big((vect(x)^((k)) - omega_k vect(d))) = 1/2 lr((vect(x)^((k)) - omega_k vect(d)))^T mat(A) lr((vect(x)^((k)) - omega_k vect(d))) - lr((vect(x)^((k)) - omega_k vect(d)))^T vect(b) \
  &= f big((vect(x)^((k)))) + omega_k^2/2 vect(d)^T mat(A) vect(d) - omega_k vect(d)^T vect(r)^((k)).
$) <eq:linear_steepest>
When viewed as a function of the real parameter $omega_k$,
the right-hand side is a convex quadratic function.
It is minimized when its derivative is equal to zero,
i.e.~when
#numbered($
  omega_k vect(d)^T mat(A) vect(d) - vect(d)^T (mat(A) vect(x)_k - vect(b)) = 0
  wide => wide omega_k = (vect(d)^T vect(r)^((k)))/(vect(d)^T mat(A) vect(d)).
$) <eq:linear_optimal_step_size_steepest_descent>
The steepest descent algorithm with step size obtained from this equation is summarized in~@algo:steepest_descent_method below.
By construction,
the function value $f(vect(x)^((k)))$ is nonincreasing with $k$,
which is equivalent to saying that
the error $vect(x) - vect(x)_*$ is nonincreasing in the norm $vect(x) |-> sqrt(vect(x)^T mat(A) vect(x))$.
In order to quantify more precisely the decrease of the error in this norm,
we introduce the notation
$
  E_k = norm(vect(x) - vect(x)_*)^2_mat(A) := (vect(x)^((k)) - vect(x)_*)^T mat(A) (vect(x)^((k)) - vect(x)_*)
  = (mat(A) vect(x)^((k)) - vect(b))^T mat(A)^(-1) (mat(A) vect(x)^((k)) - vect(b)).
$
We begin by showing the following auxiliary lemma.
#lemma(title: [Kantorovich inequality])[
  Let $mat(A) in real^(n times n)$ be a symmetric and positive definite matrix,
  and let $lambda_0 <= dots <= lambda_n$ denote its eigenvalues.
  Then for all nonzero $vect(z) in real^n$ it holds that
  $
    (vect(z)^T vect(z))^2/((vect(z)^T mat(A) vect(z)) (vect(z)^T mat(A)^(-1) vect(z)))
    >= (4 lambda_1 lambda_n)/(lambda_1 + lambda_n)^2.
  $
]
#proof[
  By the AM-GM (arithmetic mean-geometric mean) inequality,
  it holds for all $t > 0$ that
  $
    sqrt((vect(z)^T mat(A) vect(z)) (vect(z)^T mat(A)^(-1) vect(z)))
    &= sqrt((t vect(z)^T mat(A) vect(z)) (t^(-1) vect(z)^T mat(A)^(-1) vect(z)))
    <= 1/2 lr((t vect(z)^T mat(A) mat(z) + 1/t vect(z)^T mat(A)^(-1) mat(z))) \
    &= 1/2 vect(z)^T lr((t mat(A) + 1/t mat(A)^(-1))) vect(z).
  $
  The matrix on the right-hand side is also symmetric and positive definite,
  with eigenvalues equal to $t lambda_i + (t lambda_i)^(-1)$.
  Therefore, we deduce
  #numbered($
    forall t >= 0, wide
    sqrt((vect(z)^T mat(A) vect(z)) (vect(z)^T mat(A)^(-1) vect(z)))
    <= 1/2 lr((max_(i in {1, dots, n}) t lambda_i + (t lambda_i)^(-1))) vect(z)^T vect(z).
  $) <eq:linear_kantorovich_bound>
  The function $x |-> x + x^(-1)$ is convex,
  and so over any closed interval $[x_min, x_max]$
  it attains its maximum either at $x_min$ or at $x_max$.
  Consequently, it holds that
  $
    lr((max_(i in {1, dots, n}) t lambda_i + (t lambda_i)^(-1)))
    = max lr(\{ t lambda_1 + 1/(t lambda_1)\, t lambda_n + 1/(t lambda_n) \}).
  $
  In order to obtain the best possible bound from~@eq:linear_kantorovich_bound,
  we should let $t$ be such that the maximum is minimized,
  which occurs when the two arguments of the maximum are equal:
  $
    t lambda_1 + 1/(t lambda_1) = t lambda_n + 1/(t lambda_n)
    wide => wide t = 1/sqrt(lambda_1 lambda_n).
  $
  For this value of $t$,
  the maximum in~@eq:linear_kantorovich_bound is equal to
  $ sqrt(lambda_1/lambda_n) + sqrt(lambda_n/lambda_1). $
  By substituting this expression in~@eq:linear_kantorovich_bound and rearranging,
  we obtain the statement.
]

#indent We are now able to prove the convergence of the steepest descent method.
#theorem(title: [Convergence of the steepest descent method])[
  It holds that
  $
    E_(k+1) <= lr(((kappa_2(mat(A)) - 1)/(kappa_2(mat(A)) + 1)))^2 E_k.
  $
] <theorem:linear_convergenec_steepest_descent>
#proof[
  Substituting $vect(x)^((k+1)) = vect(x)^((k)) - omega_k vect(d)$ in the expression for $E_(k+1)$,
  we obtain
  $
    E_(k+1)
    &= (vect(x)^((k)) - omega_k vect(d) - vect(x)_*)^T mat(A) (vect(x)^((k)) - omega_k vect(d) - vect(x)_*) \
    &= E_k - 2 omega_k vect(d)^T vect(r)^((k)) + omega_k^2 vect(d)^T mat(A) vect(d) \
    &= E_k - (vect(d)^T vect(d))^2/(vect(d)^T mat(A) vect(d))
    = lr((1 - (vect(d)^T vect(d))^2/((vect(d)^T mat(A) vect(d)) (vect(d)^T mat(A)^(-1) vect(d))))) E_k,
  $
  Using the Kantorovich inequality,
  we have
  $
    E_(k+1)
    <= lr((1 - (4 lambda_1 lambda_n)/(lambda_1 + lambda_n)^2)) E_k
    <= lr(((lambda_1 - lambda_n)/(lambda_1 + lambda_n)))^2 E_k
    = lr(((kappa_2(mat(A)) - 1)/(kappa_2(mat(A)) + 1)))^2 E_k.
  $
  We immediately deduce the statement from this inequality.
]

#algorithm(
  caption: [Steepest descent method],
  (0, [Pick $epsilon$ and initial $vect(x)$]),
  (0, [$vect(r) <- mat(A) vect(x) - vect(b)$]),
  (0, [*while* $norm(vect(r)) >= epsilon norm(vect(b))$ *do*]),
  (1, [$vect(d) <- vect(r)$]),
  (1, [$omega <- vect(d)^T vect(r) \/ vect(d)^T mat(A) vect(d)$]),
  (1, [$vect(x) <- vect(x) - omega vect(d)$]),
  (1, [$vect(r) <- mat(A) vect(x) - vect(b)$]),
  (0, [*end while*]),
) <algo:steepest_descent_method>

#remark[
  The attentive reader will have noticed that the rate of convergence rate of the steepest descent method is not apparently better than that of Richardson's method with optimal~$omega$;
  in both cases, some upper bound on the norm of the error is multiplied by
  $ (kappa_2(mat(A)) - 1)/(kappa_2(mat(A)) + 1) $
  at each iteration.
  For the steepest descent method, this rate of convergence is always guaranteed,
  but for Richardson's method, this rate of convergence holds only for the optimal value of $omega$,
  which itself depends on the condition number~$kappa_2(mat(A))$ and is computationally expensive to approximate.
  Indeed, as we shall see in~@cha:numerical_computation_of_eigenvalues,
  the simplest method for calculating the smallest eigenvalue of a matrix,
  which is necessary for estimating the condition number, requires to solve linear systems with matrix~$mat(A)$.
]

==== Preconditioned steepest descent <ssub:preconditioned_steepest_descent>
We observe from @theorem:linear_convergenec_steepest_descent that the convergence of the steepest descent method is faster when
the condition number of the matrix $mat(A)$ is low.
This naturally leads to the following question:
can we reformulate the minimization of $f(vect(x))$ in~@eq:function_optimization
as another optimization problem which is of the same form
but involves a matrix with a lower condition number,
thereby providing scope for faster convergence?
In order to answer this question,
we consider a linear change of coordinates~$vect(y) = mat(T)^(-1) vect(x)$,
where $mat(T)$ is an invertible matrix,
and we define
#numbered($
  tilde(f)(vect(y)) = f(mat(T) vect(y))
  = 1/2 vect(y)^T (mat(T)^T mat(A) mat(T)) vect(y) - (mat(T)^T vect(b))^T vect(y).
$) <eq:linear_quadratic_function>
This function is of the same form as $f$~in~@eq:function_optimization,
with the matrix $tilde(mat(A)) := mat(T)^T mat(A) mat(T)$ instead of~$mat(A)$ and the vector $tilde(vect(b)) := mat(T)^T vect(b)$ instead of $vect(b)$.
Its minimizer is $vect(y)_* = mat(T)^(-1) vect(x)_*$.
The steepest descent algorithm can be applied to~@eq:linear_quadratic_function and,
from an approximation $vect(y)^((k))$ of the minimizer~$vect(y)_*$,
an approximation $vect(x)^((k))$ of $vect(x)_*$ is obtained by the change of variable $vect(x)^((k)) = mat(T) vect(y)^((k))$.
This approach is called _preconditioning_.
By @theorem:linear_convergenec_steepest_descent,
the steepest descent method
satisfies the following error estimate
when applied to the function~@eq:linear_quadratic_function:
$
  E_(k+1)
  <=
  lr(((kappa_2(mat(T)^T mat(A) mat(T)) - 1)/(kappa_2(mat(T)^T mat(A) mat(T)) + 1)))^2 E_k,
  wide E_k &= (vect(y)^((k)) - vect(y)_*)^T tilde(mat(A)) (vect(y)^((k)) - vect(y)_*), \
  &= (vect(x)^((k)) - vect(x)_*)^T mat(A) (vect(x)^((k)) - vect(x)_*).
$
Consequently, the convergence is faster than that of the usual steepest descent method if $kappa_2(mat(T)^T mat(A) mat(T)) < kappa_2(mat(A))$.
The optimal change of coordinates is given by $mat(T) = mat(C)^(-T)$,
where $mat(C)$ is the factor of the Cholesky factorization of $mat(A)$ as $mat(C) mat(C)^T$.
Indeed, in this case
$
  mat(T)^T mat(A) mat(T) = mat(C)^(-1) mat(C) mat(C)^T mat(C)^(-T) = mat(I)
  wide => wide
  kappa_2(mat(T)^T mat(A) mat(T)) = 1,
$
and the method converges in a single iteration!
However,
this iteration amounts to solving the linear system by direct Cholesky factorization of $mat(A)$.
In practice, it is usual to define $mat(T)$ from an approximation of the Cholesky factorization,
such as the _incomplete Cholesky factorization_.

To conclude this section,
we demonstrate that the change of variable from $vect(x)$ to $vect(y)$ need not be performed explicitly in practice.
Indeed, one step of the steepest descent algorithm applied to function $tilde(f)$ reads
$
  vect(y)^((k+1)) = vect(y)^((k)) - tilde(omega)_k (tilde(mat(A)) vect(y)^((k)) - tilde(vect(b))),
  wide tilde(omega)_k = ((tilde(mat(A)) vect(y)^((k)) - tilde(vect(b)))^T (tilde(mat(A)) vect(y)^((k)) - tilde(vect(b))))/((tilde(mat(A)) vect(y)^((k)) - tilde(vect(b)))^T tilde(mat(A)) (tilde(mat(A)) vect(y)^((k)) - tilde(vect(b)))).
$
Letting $vect(x)^((k)) = mat(T) vect(y)^((k))$,
this equation can be rewritten as the following iteration:
$
  vect(x)^((k+1)) = vect(x)^((k)) - tilde(omega)_k vect(d)_k,
  wide tilde(omega)_k = (vect(d)_k^T vect(r)^((k)))/(vect(d)_k^T mat(A) vect(d)_k), wide vect(d)_k = mat(T) mat(T)^T (mat(A) vect(x)^((k)) - vect(b)).
$
A comparison with~@eq:linear_optimal_step_size_steepest_descent shows that the step size $tilde(omega)_k$ is such that $f(vect(x)^((k+1)))$ is minimized.
This reasoning shows that the preconditioned conjugate gradient method amounts to choosing the direction~$vect(d)_k = mat(T) mat(T)^T vect(r)^((k))$ at each iteration,
instead of just $vect(r)^((k))$,
as is apparent in @algo:preconditioned_steepest_descent_method.
It is simple to check that $- vect(d)_k$ is a descent direction for $f$:
$
  -nabla f(vect(x))^T big((mat(T) mat(T)^T (mat(A) vect(x) - vect(b)))) =
  -big((mat(T)^T (mat(A) vect(x) - vect(b))))^T big((mat(T)^T (mat(A) vect(x) - vect(b)))) <= 0.
$
#algorithm(
  caption: [Preconditioned steepest descent method],
  (0, [Pick $epsilon$, invertible $mat(T)$ and initial $vect(x)$]),
  (0, [$vect(r) <- mat(A) vect(x) - vect(b)$]),
  (0, [*while* $norm(vect(r)) >= epsilon norm(vect(b))$ *do*]),
  (1, [$vect(d) <- mat(T) mat(T)^T vect(r)$]),
  (1, [$omega <- vect(d)^T vect(r) \/ vect(d)^T mat(A) vect(d)$]),
  (1, [$vect(x) <- vect(x) - omega vect(d)$]),
  (1, [$vect(r) <- mat(A) vect(x) - vect(b)$]),
  (0, [*end while*]),
) <algo:preconditioned_steepest_descent_method>

==== Conjugate directions method <ssub:conjugate_directions>
#definition(title: [Conjugate directions])[
  Let $mat(A)$ be a symmetric positive definite matrix.
  Two vectors $vect(d)_1$ and $vect(d)_2$ are called $mat(A)$-orthogonal or conjugate with respect to $mat(A)$
  if $vect(d)_1^T mat(A) vect(d)_2 = 0$,
  i.e.~if they are orthogonal for the inner product $ip(vect(x), vect(y))_mat(A) = vect(x)^T mat(A) vect(y)$.
]

#indent Assume that $vect(d)_0, dots, vect(d)_(n-1)$ are $n$ pairwise $mat(A)$-orthogonal nonzero directions.
By @exercise:linear_independence_conjugate_directions,
these vectors are linearly independent,
and so they form a basis of $real^n$.
Consequently, for any initial guess $vect(x)^((0))$,
the vector $vect(x)^((0)) - vect(x)_*$,
where $vect(x)_*$ is the solution to the linear system $mat(A) vect(x) = vect(b)$,
can be decomposed as
$ vect(x)^((0)) - vect(x)_* = alpha_0 vect(d)_0 + cdotsb + alpha_(n-1) vect(d)_(n-1). $
Taking the $ip(placeholder, placeholder)_mat(A)$ inner product of both sides with $vect(d)_k$,
with $k in {0, dots, n-1}$,
we obtain an expression for the scalar coefficient $alpha_k$:
$
  alpha_k = (vect(d)_k^T mat(A) (vect(x)^((0)) - vect(x)_*))/(vect(d)_k^T mat(A) vect(d)_k)
  = (vect(d)_k^T (mat(A) vect(x)^((0)) - vect(b)))/(vect(d)_k^T mat(A) vect(d)_k).
$
Therefore, calculating the expression of the coefficient does not require to know the exact solution $vect(x)_*$,
but only the residual $vect(r)^((0))$!
Given conjugate directions, the exact solution can be obtained as
#numbered($
  vect(x)_* = vect(x)^((0)) - sum_(k=0)^(n-1) alpha_k vect(d)_k, wide alpha_k = (vect(d)_k^T vect(r)^((0)))/(vect(d)_k^T mat(A) vect(d)_k).
$) <eq:linear_expansion_conjugate_direction>

If $vect(x)^((0)) = vect(0)$,
then $vect(r)^((0)) = -vect(b)$ and this equations gives that
$
  vect(x)_* = sum_(k=0)^(n-1) (vect(d)_k^T vect(b))/(vect(d)_k^T mat(A) vect(d)_k) vect(d)_k
  = lr((sum_(k=0)^(n-1) (vect(d)_k vect(d)_k^T)/(vect(d)_k^T mat(A) vect(d)_k))) vect(b),
$
which implies that that the inverse of~$mat(A)$ is given by
$
  mat(A)^(-1) = sum_(k=0)^(n-1) (vect(d)_k vect(d)_k^T)/nu_k, wide nu_k = vect(d)_k^T mat(A) vect(d)_k.
$
The conjugate directions method is illustrated in~@algo:conjugate_directions.
Its implementation is very similar to that of the steepest descent method,
the only difference being that the descent direction at iteration $k$ is given by $vect(d)_k$ instead of~$vect(r)^((k))$.
In particular, the step size at each iteration is such that $f(vect(x)^((k+1)))$ is minimized.
#algorithm(
  caption: [Conjugate directions method],
  (0, [Assuming $vect(d)_0, dots, vect(d)_(n-1)$ are given.]),
  (0, [Pick initial $vect(x)^((0))$]),
  (0, [*for* $k$ in ${0, dots, n-1}$ *do*]),
  (1, [$vect(r)^((k)) = mat(A) vect(x)^((k)) - vect(b)$]),
  (1, [$omega_k = vect(d)_k^T vect(r)^((k)) \/ vect(d)_k^T mat(A) vect(d)_k$]),
  (1, [$vect(x)^((k+1)) = vect(x)^((k)) - omega_k vect(d)_k$]),
  (0, [*end for*]),
) <algo:conjugate_directions>

Let us now establish the connection between the~@algo:conjugate_directions and~@eq:linear_expansion_conjugate_direction,
which may not be immediately apparent because~@eq:linear_expansion_conjugate_direction involves only the initial residual~$mat(A) vect(x)^((0)) - vect(b)$,
while the residual at the current iteration~$vect(r)^((k))$ is used in the algorithm.
#proposition(title: [Convergence of the conjugate directions method])[
  The vector $vect(x)^((k))$ obtained after $k$ iterations of the conjugate directions method is given by
  #numbered($
    vect(x)^((k)) = vect(x)^((0))
    - sum_(i=0)^(k-1) alpha_i vect(d)_i,
    wide alpha_i = (vect(d)_i^T vect(r)^((0)))/(vect(d)_i^T mat(A) vect(d)_i).
  $) <eq:series_expansion_of_y>
  In particular,
  the method converges in at most $n$ iterations.
] <proposition:conjugate_directions>
#proof[
  Let us denote by $vect(y)^((k))$ the solution obtained after $k$ steps of @algo:conjugate_directions.
  Our goal is to show that $vect(y)^((k))$ coincides with $vect(x)^((k))$ defined in~@eq:series_expansion_of_y.
  The result is trivial for $k = 0$.
  Reasoning by induction,
  we assume that it is true up to~$k$.
  Then performing step $k+1$ of the algorithm gives
  $
    vect(y)^((k+1)) = vect(y)^((k)) - omega_k vect(d)_k, wide omega_k = (vect(d)_k^T vect(r)^((k)))/(vect(d)_k^T mat(A) vect(d)_k).
  $
  On the other hand,
  it holds from~@eq:series_expansion_of_y that
  $
    vect(x)^((k+1)) = vect(x)^((k)) - alpha_k vect(d)_k, wide alpha_k = (vect(d)_k^T vect(r)^((0)))/(vect(d)_k^T mat(A) vect(d)_k).
  $
  By the induction hypothesis, it holds that $vect(y)^((k)) = vect(x)^((k))$,
  so in order to prove that $vect(y)^((k+1)) = vect(x)^((k+1))$,
  it is sufficient to show that $omega_k = alpha_k$, i.e. that
  $
    vect(d)_k^T vect(r)^((k)) = vect(d)_k^T vect(r)^((0))
    quad <=> quad
    vect(d)_k^T (vect(r)^((k)) - vect(r)^((0))) = 0
    quad <=> quad
    vect(d)_k^T mat(A) (vect(x)^((k)) - vect(x)^((0))) = 0.
  $
  The latter equality is obvious from the $mat(A)$-orthonormality of the directions.
]

#indent Since $omega_k$ in~@algo:conjugate_directions coincides with the expression in~@eq:linear_optimal_step_size_steepest_descent,
the conjugate directions algorithm satisfies the following "local optimization" property:
the iterate $vect(x)^((k+1))$ minimizes $f$ on the straight line $omega |-> vect(x)^((k)) - omega vect(d)_k$.
In contrast with the steepest descent method,
however,
the conjugate directions method also satisfies the following stronger property.
#proposition(title: [Optimality of the conjugate directions method])[
  The iterate $vect(x)^((k))$ is the minimizer of $f$ over the set $vect(x)^((0)) + cal(B)_k$,
  where $cal(B)_k = Span {vect(d)_0, dots, vect(d)_(k-1)}$.
] <proposition:optimality_conjugate_directions>
#proof[
  By~@eq:linear_expansion_conjugate_direction,
  it holds that
  $
    vect(x)_* = vect(x)^((0)) - sum_(i=0)^(n-1) alpha_i vect(d)_i,
    wide
    alpha_i = (vect(d)_i^T vect(r)^((0)))/(vect(d)_i^T mat(A) vect(d)_i)
  $
  On the other hand,
  any vector $vect(y) in vect(x)^((0)) + cal(B)_k$ can be expanded as
  $ vect(y) = vect(x)^((0)) - beta_0 vect(d)_0 - cdotsb - beta_(k-1) vect(d)_(k-1). $
  Employing these two expressions, the formula for $f$ in~@eq:rewritten_expression_of_f,
  and the $mat(A)$-orthogonality of the directions,
  we obtain
  $
    f(vect(y))
    &= 1/2 (vect(y) - vect(x)_*)^T mat(A) (vect(y) - vect(x)_*) - 1/2 vect(x)_*^T mat(A) vect(x)_* \
    &= 1/2 sum_(i=0)^(k-1) (beta_i - alpha_i)^2 vect(d)_i^T mat(A) vect(d)_i
    + 1/2 sum_(i=k)^(n-1) alpha_i^2 vect(d)_i^T mat(A) vect(d)_i - 1/2 vect(x)_*^T mat(A) vect(x)_*
  $
  This is minimized when $beta_i = alpha_i$ for all $i in {0, dots, k-1}$,
  in which case $vect(y)$ coincides with the $k$-th iterate $vect(x)^((k))$ of the conjugate directions method in view of~@proposition:conjugate_directions.
]
#remark[
  Let $norm(placeholder)_mat(A)$ denote the norm induced by the inner product $ip(placeholder, placeholder)_mat(A)$.
  Since
  $
    norm(vect(x)^((k)) - vect(x)_*)_mat(A) = sqrt(2 f(vect(x)^((k))) + vect(x)_*^T mat(A) vect(x)_*),
  $
  @proposition:optimality_conjugate_directions shows that $vect(x)^((k))$ minimizes the norm~$norm(vect(x)^((k)) - vect(x)_*)_mat(A)$
  over $vect(x)^((0)) + cal(B)_k$.
  This is not surprising since,
  by construction, the vector $vect(x)^((k)) - vect(x)^((0))$ is the orthogonal projection
  of~$vect(x)_* - vect(x)^((0))$ onto $cal(B)_k$,
  for the inner product $ip(placeholder, placeholder)_mat(A)$.
]

#indent A corollary of~#eqref-number(<proposition:optimality_conjugate_directions>) is that the gradient of~$f$ at $vect(x)^((k))$,
i.e. the residual $vect(r)^((k)) = mat(A) vect(x)^((k)) - vect(b)$,
is orthogonal to any vector in ${vect(d)_0, dots, vect(d)_(k-1)}$ for the usual Euclidean inner product.
This can also be checked directly from the formula
$
  vect(x)^((k)) - vect(x)_* = sum_(i=k)^(n-1) alpha_i vect(d)_i,
  wide alpha_i = (vect(d)_i^T vect(r)^((0)))/(vect(d)_i^T mat(A) vect(d)_i),
$
which follows directly from @proposition:conjugate_directions.
Indeed, it holds that
#numbered($
  forall j in {0, dots, k-1}, wide
  vect(d)_j^T vect(r)^((k)) = vect(d)_j mat(A) (vect(x)^((k)) - vect(x)_*)
  = sum_(i=k)^(n-1) alpha_i vect(d)_j^T mat(A) vect(d)_i = 0.
$) <eq:orthogonality_residual>

==== The conjugate gradient method <ssub:the_conjugate_gradient_method>
In the previous section,
we showed that, given $n$ conjugate directions,
the solution to the linear system $mat(A) vect(x) = vect(b)$ can be obtained in a finite number of iterations using~@algo:conjugate_directions.
The conjugate gradient method can be viewed as a particular case of the conjugate directions method.
Instead of assuming that the conjugate directions are given,
they are constructed iteratively as part of the algorithm.
Given an initial guess $vect(x)^((0))$,
the first direction is the residual~$vect(r)^((0))$,
which coincides with the gradient of~$f$ at~$vect(x)^((0))$.
The directions employed for the next iterations are obtained by applying the Gram-Schmidt process
to the residuals.
More precisely, given conjugate directions~$vect(d)_0, dots, vect(d)_(k-1)$,
and letting $vect(x)^((k))$ denote the $k$-th iterate of the conjugate directions method,
the direction $vect(d)_k$ is obtained by
#numbered($
  vect(d)_k = vect(r)^((k)) - sum_(i=0)^(k-1) (vect(d)_i^T mat(A) vect(r)^((k)))/(vect(d)_i^T mat(A) vect(d)_i) vect(d)_i,
  wide vect(r)^((k)) = mat(A) vect(x)^((k)) - vect(b).
$) <eq:conjugate_gradient_gram_schmidt>
It is simple to check that $vect(d)_k$ is indeed $mat(A)$-orthogonal to $vect(d)_i$ for $i in {0, dots, k-1}$,
and that $vect(d)_k$ is nonzero if $vect(r)^((k))$ is nonzero.
To prove the latter claim,
we can take the Euclidean inner product of both sides with $vect(r)^((k))$ and use~@proposition:optimality_conjugate_directions
to deduce that
#numbered($ vect(d)_k^T vect(r)^((k)) = (vect(r)^((k)))^T vect(r)^((k)) > 0. $) <eq:conjugate_gradient_positive_omega>
Note also that since the directions are obtained by applying the Gram--Schmidt process to the residuals,
it holds that
#numbered($
  forall k in {0, dots, n-1}, wide
  cal(B)_(k+1) := Span lr(\{ vect(d)_0\, dots\, vect(d)_k \}) =
  Span lr(\{ vect(r)^((0))\, dots\, vect(r)^((k)) \}).
$) <eq:conjugate_gradient_equality_span>
The following result characterizes precisely the subspace $cal(B)_(k+1)$.
#proposition[
  Assume that $norm(vect(r)^((k))) != 0$ for all $k < m <= n$.
  Then it holds that
  #numbered($
    forall k in {0, dots m}, quad
    Span lr(\{ vect(r)^((0))\, vect(r)^((1))\, dots\, vect(r)^((k)) \})
    = Span lr(\{ vect(r)^((0))\, mat(A) vect(r)^((0))\, dots\, mat(A)^k vect(r)^((0)) \})
  $) <eq:conjugate_gradient_krylov>
  The subspace on the right-hand side is called a Krylov subspace.
] <proposition:krylov_subspace>
#proof[
  The result is clear for $k = 0$.
  Reasoning by induction,
  we prove that if the result is true up to $k < m$,
  then it is also true for $k+1$.
  A simple calculation gives that
  #eq-last($
    vect(r)^((k+1))
    &= mat(A) lr((vect(x)^((k)) - omega_k vect(d)_k)) - vect(b) \
    &= vect(r)^((k)) - omega_k mat(A) vect(d)_k.
  $) <eq:residual_iteration>
  From~@eq:conjugate_gradient_gram_schmidt,
  we deduce that
  $
    vect(r)^((k+1))
    = vect(r)^((k)) - omega_k mat(A) lr((vect(r)^((k)) - sum_(i=0)^(k-1) (vect(d)_i^T mat(A) vect(r)^((k)))/(vect(d)_i^T mat(A) vect(d)_i) vect(d)_i)).
  $
  By~@eq:conjugate_gradient_equality_span and the induction hypothesis,
  the bracketed expression on the right-hand side belongs to $cal(B)_(k+1)$,
  so the inclusion~$subset$ in~@eq:conjugate_gradient_krylov is clear.
  The inclusion $supset$ then follows from the fact the dimension of the subspace
  $
    Span lr(\{ vect(d)_0\, dots\, vect(d)_k \}) =
    Span lr(\{ vect(r)^((0))\, dots\, vect(r)^((k)) \})
  $
  is equal to $k+1$.
]

#indent It appears from~@eq:conjugate_gradient_gram_schmidt that the cost of calculating a new direction grows linearly with the iteration index.
In fact, it turns out that only the last term in the sum is nonzero,
and so the cost of calculating a new direction does not grow with the iteration index~$k$.
Indeed, notice that if $i <= k-2$, then
$
  vect(d)_i^T mat(A) vect(r)^((k))
  = (mat(A) vect(d)_i)^T mat(A) vect(r)^((k))
  = 0,
$
because $mat(A) vect(d)_i in cal(B)_(i+2) subset cal(B)_k$ by~@proposition:krylov_subspace,
and $vect(r)^((k))$ orthogonal to~$cal(B)_k$ for the Euclidean inner product by~@eq:orthogonality_residual.
This observation leads to~@algo:conjugate_gradient.
#algorithm(
  caption: [Conjugate gradient method],
  (0, [Pick initial $vect(x)^((0))$]),
  (0, [$vect(d)_0 = vect(r)^((0)) = mat(A) vect(x)^((0)) - vect(b)$]),
  (0, [*for* $k$ in ${0, dots, n-1}$ *do*]),
  (1, [*if* $norm(vect(r)^((k))) = 0$ *then*]),
  (2, [Stop]),
  (1, [*end if*]),
  (1, [$omega_k = vect(d)_k^T vect(r)^((k)) \/ vect(d)_k^T mat(A) vect(d)_k$]),
  (1, [$vect(x)^((k+1)) = vect(x)^((k)) - omega_k vect(d)_k$]),
  (1, [$vect(r)^((k+1)) = mat(A) vect(x)^((k+1)) - vect(b)$]),
  (1, [$beta_k = vect(d)_k^T mat(A) vect(r)^((k+1)) \/ vect(d)_k^T mat(A) vect(d)_k$.]),
  (1, [$vect(d)_(k+1) = vect(r)^((k+1)) - beta_k vect(d)_k$.]),
  (0, [*end for*]),
) <algo:conjugate_gradient>

Although the conjugate gradient method converges in a finite number of iterations,
performing $n$ iterations for very large systems would require an excessive computational cost,
and so it is sometimes desirable to stop iterating when the residual is sufficiently small.
To conclude this section, we study the convergence of the method.
#theorem(title: [Convergence of the conjugate gradient method])[
  The error for the conjugate gradient method,
  measured as
  $ E_k := (vect(x)^((k)) - vect(x)_*)^T mat(A) (vect(x)^((k)) - vect(x)_*), $
  satisfies the following inequality:
  #numbered($
    forall q_k in poly(k), wide
    E_(k+1) <= max_(1 <= i <= n) big((1 + lambda_i q_k (lambda_i)))^2 E_0.
  $) <equation:convergence_conjugate_gradient>
  Here $poly(k)$ is the vector space of polynomials of degree less than or equal to $k$.
] <theorem:convergence_conjugate_gradient>
#proof[
  In view of @proposition:krylov_subspace,
  the iterate $vect(x)^((k+1))$ can be written as
  $
    vect(x)^((k+1)) = vect(x)^((0)) + sum_(i=0)^k alpha_i mat(A)^i vect(r)^((0))
    = vect(x)^((0)) + p_k (mat(A)) vect(r)^((0)),
  $
  where $p_k$ is a polynomial of degree $k$.
  By~@proposition:optimality_conjugate_directions,
  $p_k$ is in fact the polynomial of degree~$k$ such that $f(vect(x)^((k+1)))$ is minimized,
  and thus also $E_(k+1)$ by~@eq:rewritten_expression_of_f.
  Noting that
  $
    vect(x)^((k+1)) - vect(x)_*
    &= vect(x)^((0)) - vect(x)_* + p_k (mat(A)) vect(r)^((0))
    = vect(x)^((0)) - vect(x)_* + p_k (mat(A)) mat(A) (vect(x)^((0)) - vect(x)_*) \
    &= big((mat(I) + mat(A) p_k (mat(A)))) (vect(x)^((0)) - vect(x)_*),
  $
  we deduce that
  $
    forall q_k in poly(k), wide
    E_(k+1) <= (vect(x)^((0)) - vect(x)_*)^T mat(A) big((mat(I) + mat(A) q_k (mat(A))))^2 (vect(x)^((0)) - vect(x)_*).
  $
  In order to exploit this inequality,
  it is useful to diagonalize $mat(A)$ as $mat(A) = mat(Q) mat(D) mat(Q)^T$,
  for an orthogonal matrix $mat(Q)$ and a diagonal matrix $mat(D)$.
  Since $q_k (mat(A)) = mat(Q) q_k (mat(D)) mat(Q)^T$ for all $q_k in poly(k)$,
  it holds that
  $
    forall q_k in poly(k), wide
    E_(k+1)
    &= big((mat(Q)^T (vect(x)^((0)) - vect(x)_*)))^T mat(D) big((mat(I) + mat(D) q_k (mat(D))))^2 big((mat(Q)^T (vect(x)^((0)) - vect(x)_*))) \
    &<= max_(1 <= i <= n) big((1 + lambda_i q_k (lambda_i)))^2 underbrace(big((mat(Q)^T (vect(x)^((0)) - vect(x)_*)))^T mat(D) big((mat(Q)^T (vect(x)^((0)) - vect(x)_*))), E_0),
  $
  which completes the proof.
]

#indent A corollary of~@theorem:convergence_conjugate_gradient is that,
if $mat(A)$ has $m <= n$ distinct eigenvalues,
then the conjugate gradient method converges in at most $m$ iterations.
Indeed, in this case we can take
$
  q_(m-1) (lambda) = 1/lambda lr((((lambda_1 - lambda) dots (lambda_m - lambda))/(lambda_1 dots lambda_m) - 1)).
$
It is simple to check that the right-hand side is indeed a polynomial,
and that $1 + lambda_i q_(m-1) (lambda_i) = 0$ for all eigenvalues of $mat(A)$.

In general, finding the polynomial that minimizes the right-hand side of~@equation:convergence_conjugate_gradient is not possible,
because the eigenvalues of $mat(A)$ are unknown.
However, it is possible to derive from this equation an error estimate with an explicit dependence on the condition number~$kappa = kappa_2(mat(A))$.
#theorem[
  It holds that
  $
    forall k >= 0, wide
    E_k <= 4 lr(((sqrt(kappa) - 1)/(sqrt(kappa) + 1)))^(2k) E_0,
  $
]
#proof[
  @theorem:convergence_conjugate_gradient implies that
  $
    forall q_k in poly(k), wide
    E_(k+1) <= max_(lambda in [lambda_1, lambda_n]) big((1 + lambda q_k (lambda)))^2 E_0,
  $
  where $lambda_1$ and $lambda_n$ are the minimum and maximum eigenvalues of $mat(A)$.
  Notice that
  $
    Big(\{) 1 + lambda q_k colon q_k in poly(k) Big(\})
    = Big(\{) p_k colon p_k in poly(k+1) " and " p_k (0) = 1 Big(\})
  $
  Therefore, it follows from~@exercise:conjugate_gradient_chebyshev that the right-hand side is minimized when
  #numbered($
    1 + lambda q_k (lambda) = (T_(k+1) lr(((lambda_n + lambda_1 - 2 lambda)/(lambda_n - lambda_1))))/(T_(k+1) lr(((lambda_n + lambda_1)/(lambda_n - lambda_1)))),
  $) <eq:conjugate_gradient_optimal_polynom>
  where $T_(k+1)$ is the _Chebyshev_ polynomial of degree $k+1$,
  see~@eq:chebyshev.
  We recall that $|T_(k+1) (x)| <= 1$ for all $x in [-1, 1]$.
  Consequently, by the expression of Chebyshev polynomials given in~@exercise:yet_another_expression_cheby,
  the following inequality holds true for all $lambda in [lambda_1, lambda_n]$:
  $
    big(|1 + lambda q_k (lambda)|)
    <= 1/(T_(k+1) lr(((lambda_n + lambda_1)/(lambda_n - lambda_1))))
    &= 2 lr((Big((r + sqrt(r^2 - 1)))^(k+1) + Big((r - sqrt(r^2 - 1)))^(k+1)))^(-1), \
    &= 2 lr((lr(((sqrt(kappa) + 1)/(sqrt(kappa) - 1)))^(k+1) + lr(((sqrt(kappa) - 1)/(sqrt(kappa) + 1)))^(k+1)))^(-1).
  $
  where $r = (lambda_n + lambda_1)/(lambda_n - lambda_1)$.
  Since the first term in the bracket converges to zero as $k -> infinity$,
  it is natural to bound this expression by keeping only the second term,
  which after simple algebraic manipulations leads to
  $
    forall lambda in [lambda_1, lambda_n], wide
    big(|1 + lambda q_k (lambda)|)
    <= 2 lr(((sqrt(kappa) - 1)/(sqrt(kappa) + 1)))^(k+1).
  $
  From this inequality,
  the statement of the theorem follows immediately.
]

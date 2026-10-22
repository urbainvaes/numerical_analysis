#import "../template.typ": *

// \vdots and \ddots (`dots` is redefined by the template)
#let vdots = sym.dots.v
#let ddots = sym.dots.down
// Environments are figures, which are not breakable by default: without this,
// a list inside an exercise cannot break across pages (should go in the template)
#show selector.or(..theorem-kinds.map(k => figure.where(kind: k))): set block(breakable: true)

= Numerical computation of eigenvalues <cha:numerical_computation_of_eigenvalues>
#minitoc()

#section-star[Introduction]

Calculating the eigenvalues and eigenvectors of a matrix is a task often encountered in scientific and engineering applications.
Eigenvalue problems naturally arise in quantum physics,
solid mechanics, structural engineering and molecular dynamics,
to name just a few applications.
The aim of this chapter is to present an overview of the standard methods for calculating eigenvalues and eigenvectors numerically.
We focus predominantly on the case of a Hermitian matrix $mat(A) in complex^(n times n)$,
which is technically simpler and arises in many applications.
The reader is invited to go through the background material in~@sec:diagonalization before reading this chapter.
The rest of this chapter is organized as follows
- In @sec:general_remarks,
  we make general remarks concerning the calculation of eigenvalues.

- In @sec:simple_vector_iterations,
  we present standard methods based on a simple vector iteration.

- In @sec:subspace_iteration,
  we present a method for calculating several eigenvectors simultaneously,
  based on iterating a subspace.

- In @sec:projection_methods,
  we present method for constructing an approximation of the eigenvectors in a given subspace of~$complex^n$.

== Numerical methods for eigenvalue problems: general remarks <sec:general_remarks>

As mentioned in~@sec:diagonalization,
a complex number $lambda in complex$ is an eigenvalue of $mat(A) in complex^(n times n)$
if and only if $lambda$ is a root of the characteristic polynomial $p_A colon complex -> complex$ of $mat(A)$,
which is given by
$ p_A (lambda) = det(mat(A) - lambda mat(I)). $
One may,
therefore,
calculate the eigenvalues of $mat(A)$ by calculating the roots of the polynomial~$p_A$ using,
for example, one of the methods presented in~@cha:solution_of_nonlinear_systems.
While feasible for small matrices,
this approach is not viable for large matrices,
because the number of floating point operations required for calculating calculating the coefficients of the characteristic polynomial scales as the factorial of $n$.

In view of the prohibitive computational cost required for calculating the characteristic polynomial,
other methods are required for solving large eigenvalue problems numerically.
All the methods that we study in this chapter are of iterative nature.
While some of them are aimed at calculating all the eigenpairs of the matrix~$mat(A)$,
other methods enable to calculate only a small number of eigenpairs at a lower computational cost,
which is often desirable.
Indeed, calculating all the eigenvalues of a large matrix is computationally expensive;
on a personal computer, the following Julia code takes well over a second to terminate:
```julia
import LinearAlgebra
A = rand(2000, 2000)
LinearAlgebra.eigen(A)
```

In many applications,
the matrix $mat(A)$ is sparse,
and in this case it is important to use algorithms for eigenvalue problems that do not destroy the sparsity structure.
Note that the eigenvectors of a sparse matrix are generally not sparse.

To conclude this section,
we introduce some notation used throughout this chapter.
For a diagonalizable matrix $mat(A)$,
we denote the eigenvalues by $lambda_1, dots, lambda_n$,
with $abs(lambda_1) >= abs(lambda_2) >= dots >= abs(lambda_n)$.
The associated normalized eigenvectors are denoted by $vect(v)_1, dots, vect(v)_n$.
Therefore, it holds that
$
  mat(V)^(-1) mat(A) mat(V) = mat(D) = diag(lambda_1, dots, lambda_n),
  wide "where " mat(V) = pmat(vect(v)_1, dots, vect(v)_n).
$

== Simple vector iterations <sec:simple_vector_iterations>

In this section,
we present simple iterative methods aimed at calculating just one eigenvector of the matrix~$mat(A)$,
which we assume to be diagonalizable for simplicity.

=== The power iteration
The power iteration is the simplest method for calculating the eigenpair associated with the eigenvalue of~$mat(A)$ with largest modulus.
Since the eigenvectors of $mat(A)$ span $complex^n$,
any vector $vect(x)_0$ may be decomposed as
#numbered($ vect(x)_0 = alpha_1 vect(v)_1 + cdots + alpha_n vect(v)_n. $) <eq:eigen_decomposition_x0>
The idea of the power iteration is to repeatedly left-multiply this vector by the matrix $mat(A)$,
in order to amplify the coefficient of $vect(v)_0$ relative to the other ones.
Indeed, notice that
$ mat(A)^k vect(x)_0 = lambda_1^k alpha_1 vect(v)_1 + cdots + lambda_n^k alpha_n vect(v)_n. $
If $lambda_1$ is strictly greater in modulus than the other eigenvalues,
and if $alpha_1 != 0$,
then for large $k$ the vector $mat(A)^k vect(x)_0$ is approximately aligned,
in a sense made precise below,
with the eigenvector~$vect(v)_1$.
In order to avoid overflow errors at the numerical level,
the iterates are normalized at each iteration.
The power iteration is presented in~@algo:power_iteration.
#algorithm(
  caption: [Power iteration],
  (0, [$vect(x) <- vect(x)_0$]),
  (0, [*for* $i in {1, 2, dots}$ *do*]),
  (1, [$vect(x) <- mat(A) vect(x)$]),
  (1, [$vect(x) <- vect(x) \/ norm(vect(x))$]),
  (0, [*end for*]),
) <algo:power_iteration>

#indent To precisely quantify the convergence of the power method,
we introduce the notion of _acute angle_ between vectors of $complex^n$.
$
  angle(vect(x), vect(y))
  &= arccos lr(((|vect(x)^* vect(y)|)/(sqrt(vect(x)^* vect(x)) sqrt(vect(y)^* vect(y))))) \
  &= arcsin lr(((norm((mat(I) - mat(P)_vect(y)) vect(x)))/(norm(vect(x))))),
  wide mat(P)_vect(y) := (vect(y) vect(y)^*)/(vect(y)^* vect(y)).
$
This definition generalizes the familiar notion of angle for vectors in $real^2$ or $real^3$,
and we note that the angle function satisfies $angle(ee^(ii theta_1) vect(x), ee^(ii theta_2) vect(y)) = angle(vect(x), vect(y))$ as well as $angle(vect(x), vect(y)) in [0, pi\/2]$.
We can then prove the following convergence result.
#proposition(title: [Convergence of the power iteration])[
  Suppose that $mat(A)$ is diagonalizable and that $abs(lambda_1) > abs(lambda_2)$.
  Then, for every initial guess with $alpha_1 != 0$,
  the sequence $(vect(x)_k)_(k >= 0)$ generated by the power iteration satisfies
  $ lim_(k -> infinity) angle(vect(x)_k, vect(v)_1) = 0. $
] <proposition:convergence_of_the_power_iteration>
#proof[
  By construction,
  it holds that
  #numbered($
    vect(x)_k
    = (lambda_1^k alpha_1 vect(v)_1 + cdots + lambda_n^k alpha_n vect(v)_n)/norm(lambda_1^k alpha_1 vect(v)_1 + cdots + lambda_n^k alpha_n vect(v)_n)
    = ee^(ii theta_k) (vect(v)_1 + (lambda_2^k alpha_2)/(lambda_1^k alpha_1) vect(v)_2 + cdots + (lambda_n^k alpha_2)/(lambda_1^k alpha_1) vect(v)_n)/lr(||vect(v)_1 + (lambda_2^k alpha_2)/(lambda_1^k alpha_1) vect(v)_2 + cdots + (lambda_n^k alpha_n)/(lambda_1^k alpha_1) vect(v)_n||),
  $) <eq:power_iteration_iterate>
  where
  $ ee^(ii theta_k) := (lambda_1^k alpha_1)/abs(lambda_1^k alpha_1). $
  It follows from~@eq:power_iteration_iterate that $ee^(-ii theta_k) vect(x)_k -> vect(v)_1 \/ norm(vect(v)_1) = vect(v)_1$ in the limit as $k -> infinity$,
  where we employed the fact that $norm(vect(v)_1) = 1$.
  Using the definition of the angle between two vectors in $complex^n$,
  and the continuity with respect to either argument of the $complex^n$ Euclidean inner product and of the $arccos$ function,
  we obtain that
  $
    angle(vect(x)_k, vect(v)_1)
    &= arccos lr(((|vect(v)_1^* vect(x)_k|)/(sqrt(vect(v)_1^* vect(v)_1) sqrt(vect(x)_k^* vect(x)_k)))) = arccos lr((abs(vect(v)_1^* vect(x)_k))) \
    &= arccos lr((lr(|vect(v)_1^* lr((ee^(-ii theta_k) vect(x)_k))|)))
    stretch(->)_(k -> infinity) arccos(1) = 0,
  $
  which concludes the proof.
]

#indent An inspection of the proof also reveals that the dominant term in the error,
asymptotically in the limit as $k -> infinity$,
is the one with coefficient $(lambda_2^k alpha_2)/(lambda_1^k alpha_1)$.
Therefore, we deduce that
$ angle(vect(x)_k, vect(v)_1) = bigo lr((lr(|lambda_2/lambda_1|)^k)). $
The convergence is slow if $abs(lambda_2\/lambda_1)$ is close to one,
and fast if $abs(lambda_2) << abs(lambda_1)$.
Once an approximation of the eigenvector $vect(v)_1$ has been calculated,
the corresponding eigenvalue $lambda_1$ can be estimated from the _Rayleigh quotient:_
#numbered($
  rho_mat(A) colon complex^n_* -> complex colon vect(x) |-> (vect(x)^* mat(A) vect(x))/(vect(x)^* vect(x)).
$) <eq:rayleigh_quotient>
For any eigenvector~$vect(v)$ of $mat(A)$,
the corresponding eigenvalue is equal to~$rho_mat(A) (vect(v))$.
In order to study the error on the eigenvalue~$lambda_1$ for the power iteration,
we assume for simplicity that $mat(A)$ is Hermitian
and that the eigenvectors $vect(v)_1, dots, vect(v)_n$ are orthonormal.
Substituting~@eq:power_iteration_iterate in the Rayleigh quotient~@eq:rayleigh_quotient,
we obtain
$
  rho_mat(A) (vect(x)_k)
  = (lambda_1 + lr(|(lambda_2^k alpha_2)/(lambda_1^k alpha_1)|)^2 lambda_2 + cdots + lr(|(lambda_n^k alpha_n)/(lambda_1^k alpha_1)|)^2 lambda_n)
  /(1 + lr(|(lambda_2^k alpha_2)/(lambda_1^k alpha_1)|)^2 + cdots + lr(|(lambda_n^k alpha_n)/(lambda_1^k alpha_1)|)^2).
$
Therefore,
by reducing to a common denominator we deduce
$
  abs(rho_mat(A) (vect(x)_k) - lambda_1)
  &= lr(|(lambda_1 + lr(|(lambda_2^k alpha_2)/(lambda_1^k alpha_1)|)^2 lambda_2 + cdots + lr(|(lambda_n^k alpha_n)/(lambda_1^k alpha_1)|)^2 lambda_n)
  /(1 + lr(|(lambda_2^k alpha_2)/(lambda_1^k alpha_1)|)^2 + cdots + lr(|(lambda_n^k alpha_n)/(lambda_1^k alpha_1)|)^2) - lambda_1|) \
  &<= lr(|(lambda_2^k alpha_2)/(lambda_1^k alpha_1)|)^2 abs(lambda_2 - lambda_1) + cdots + lr(|(lambda_n^k alpha_n)/(lambda_1^k alpha_1)|)^2 abs(lambda_n - lambda_1)
  = bigo lr((lr(|lambda_2/lambda_1|)^(2k))).
$

The convergence of the eigenvalue in the particular case of a Hermitian matrix is faster than for a general matrix in~$complex^(n times n)$.
For general matrices,
it is possible to show using a similar argument
that the error is of order~$bigo big((abs(lambda_2\/lambda_1)^k))$
in the limit as $k -> infinity$.

#paragraph[Essential convergence.][
  It is useful at this point to introduce the concept of _essential convergence_.
  A sequence~$(vect(x)_k)$ in $complex^n$ is said to _converge essentially_ to a vector~$vect(x)_infinity$ if
  there exists a sequence of complex numbers~$lr((ee^(ii phi_k)))$ of modulus 1 such that
  the sequence $lr((ee^(ii phi_k) vect(x)_k))$ converges to $vect(x)_infinity$.
  For a sequence $(vect(x)_k)$ of normalized vectors,
  the essential convergence of $(vect(x)_k)$ to $vect(x)_infinity$ is equivalent to the convergence of $angle(vect(x)_k, vect(x)_infinity)$ to 0.
  Proving this equivalence is the goal of~@exercise:essential_convergence.
  Reformulated in this new terminology,
  @proposition:convergence_of_the_power_iteration states that the sequence~$(vect(x)_k)$ obtained from the power iteration converges essentially to $vect(v)_1$.
]

=== Inverse iteration
The power iteration is simple but enables to calculate only the dominant eigenvalue of the matrix~$mat(A)$,
i.e.~the eigenvalue of largest modulus.
In addition, the convergence of the method is slow when $abs(lambda_2) approx abs(lambda_1)$.

The inverse iteration enables a more efficient calculation of not only the dominant eigenvalue
but also the other eigenvalues of~$mat(A)$.
It is based on applying the power iteration to $(mat(A) - mu mat(I))^(-1)$,
where $mu in complex$ is a shift.
The eigenvalues of $(mat(A) - mu mat(I))^(-1)$ are given by $(lambda_1 - mu)^(-1), dots, (lambda_n - mu)^(-1)$,
with associated eigenvectors $vect(v)_1, dots, vect(v)_n$.
If $0 < abs(lambda_J - mu) < abs(lambda_j - mu)$ for all $j != J$,
then the dominant eigenvalue of the matrix $(mat(A) - mu mat(I))^(-1)$ is $(lambda_J - mu)^(-1)$,
and so the power iteration applied to this matrix yields an approximation of the eigenvector $vect(v)_J$.
In other words,
the inverse iteration with shift~$mu$ enables to calculate an approximation of the eigenvector of~$mat(A)$ corresponding to the eigenvalue nearest~$mu$.
The inverse iteration is presented in~@algo:inverse_iteration.
In practice, the inverse matrix $(mat(A) - mu mat(I))^(-1)$ need not be calculated,
and it is often preferable to solve a linear system at each iteration.
#algorithm(
  caption: [Inverse iteration],
  (0, [$vect(x) <- vect(x)_0$]),
  (0, [*for* $i in {1, 2, dots}$ *do*]),
  (1, [Solve $(mat(A) - mu mat(I)) vect(y) = vect(x)$]),
  (1, [$vect(x) <- vect(y) \/ norm(vect(y))$]),
  (0, [*end for*]),
  (0, [$lambda <- vect(x)^* mat(A) vect(x) \/ vect(x)^* vect(x)$]),
  (0, [*return* $vect(x), lambda$]),
) <algo:inverse_iteration>

#indent An application of~@proposition:convergence_of_the_power_iteration immediately gives the following convergence result for the inverse iteration.
#proposition(title: [Convergence of the inverse iteration])[
  Assume that $mat(A) in complex^n$ is diagonalizable
  and that there exist $J$ and $K$ such that
  $ 0 < abs(lambda_J - mu) < abs(lambda_K - mu) <= abs(lambda_j - mu) wide forall j != J. $
  Assume also that $alpha_J != 0$,
  where $alpha_J$ is the coefficient of $vect(v)_J$ in the expansion of $vect(x)_0$ given in~@eq:eigen_decomposition_x0.
  Then the iterates of the inverse iteration satisfy
  $ lim_(k -> infinity) angle(vect(x)_k, vect(v)_J) = 0. $
  More precisely,
  $ angle(vect(x)_k, vect(v)_J) = bigo lr((lr(|(lambda_J - mu)/(lambda_K - mu)|)^k)). $
] <proposition:inverse_iteration>
@proposition:inverse_iteration states that $vect(x)_k$ converges essentially to~$vect(v)_J$.
Notice that the closer $mu$ is to~$lambda_J$,
the faster the inverse iteration converges.
Note also that with $mu = 0$,
the inverse iteration enables to calculate the eigenvalue of~$mat(A)$ of smallest modulus.

=== Rayleigh quotient iteration
Since the inverse iteration is fast when $mu$ is close to an eigenvalue $lambda_J$,
it is natural to wonder whether the method can be improved by progressively updating $mu$ as the simulation progresses.
Specifically, an approximation of the eigenvalue associated with the current vector may be employed in place of~$mu$.
This leads to the Rayleigh quotient iteration,
presented in~@algo:rayleigh_quotient.
#algorithm(
  caption: [Inverse iteration],
  (0, [$vect(x) <- vect(x)_0$]),
  (0, [*for* $i in {1, 2, dots}$ *do*]),
  (1, [$mu <- vect(x)^* mat(A) vect(x) \/ vect(x)^* vect(x)$]),
  (1, [Solve $(mat(A) - mu mat(I)) vect(y) = vect(x)$]),
  (1, [$vect(x) <- vect(y) \/ norm(vect(y))$]),
  (0, [*end for*]),
  (0, [$lambda <- vect(x)^* mat(A) vect(x) \/ vect(x)^* vect(x)$]),
  (0, [*return* $vect(x), lambda$]),
) <algo:rayleigh_quotient>

#indent It is possible to show that, when~$mat(A)$ is Hermitian,
the Rayleigh quotient iteration converges to an eigenvector for almost every initial guess~$vect(x)_0$.
Furthermore, if convergence to an eigenvector occurs,
then $mu$ converges cubically to the corresponding eigenvalue.
See~@MR3396212 and the references therein for more details.

== Methods based on a subspace iteration <sec:subspace_iteration>
The subspace iteration resembles the power iteration but it is more general:
not just one but several vectors are updated at each iteration.

=== Simultaneous iteration <sub:simultaneous_iteration>

Let $mat(X)_0 = pmat(vect(x)_1, dots, vect(x)_p)$ denote an initial set of linearly independent vectors.
Before we present the simultaneous iteration,
we recall a statement concerning the $mat(Q) mat(R)$ decomposition of a matrix,
which is related to the Gram--Schmidt orthonormalization process.
We recall that the Gram--Schmidt method enables to construct,
starting form an ordered set of vectors ${vect(x)_1, dots, vect(x)_p}$ in~$complex^n$,
a new set of vectors ${vect(q)_1, dots, vect(q)_p}$
which are _orthonormal_ and span the same subspace of $complex^n$ as the original vectors.
#proposition(title: [Reduced QR decomposition])[
  Assume that $mat(X) in complex^(n times p)$ has linearly independent columns.
  Then there exist a matrix $mat(Q) in complex^(n times p)$ with orthonormal columns and
  an upper triangular matrix~$mat(R) in complex^(p times p)$ such that the following factorization holds:
  #numbered($ mat(X) = mat(Q) mat(R). $) <eq:qr>
  This decomposition is known as a _reduced $mat(Q) mat(R)$ decomposition_ if $p < n$,
  or simply _$mat(Q) mat(R)$ decomposition_ if $p = n$,
  in which case $mat(X)$ is a square matrix and $mat(Q)$ is a unitary matrix.
  The decomposition is unique if we require that the diagonal elements of~$mat(R)$ are real and positive.
]
#proof[
  The statement is clear when $p=1$.
  Reasoning by induction,
  we assume that the result is true up to $p-1$, and prove that it then also holds true for~$p$.
  We wish to show that there is a unique matrix $mat(Q) in complex^(n times p)$ with orthonormal columns and
  a unique upper triangular matrix~$mat(R) in complex^(p times p)$ with real and positive diagonal elements such that~@eq:qr is satisfied.
  To this end, we decompose the matrices $mat(Q)$ and $mat(R)$ as follows:
  #numbered($
    mat(Q) = pmat(mat(Q)_(p-1), vect(q)),
    wide
    mat(R) = pmat(mat(R)_(p-1), vect(r); vect(0)_(p-1)^T, r).
  $) <eq:qr_decomposition_second_expression>
  Here $mat(Q)_(p-1) in complex^(n times (p-1))$ is a matrix with orthonormal columns,
  $mat(R) in complex^((p-1) times (p-1))$ is an upper triangular matrix with positive real diagonal elements,
  $vect(q) in complex^n$ is a normalized vector orthogonal to all the columns of~$mat(Q)_(p-1)$,
  $vect(r) in complex^(n-1)$ is a vector
  and $r in real_(>0)$ is a scalar.
  Let us also denote by~$mat(X)_(p-1) in complex^(n times (p-1))$ the matrix containing the $p-1$ first columns of~$mat(X)$,
  and by $vect(x)_p in complex^n$ the $p$-th column of $mat(X)$.
  Substituting~@eq:qr_decomposition_second_expression into~@eq:qr,
  we then obtain
  #numbered($
    pmat(mat(X)_(p-1), vect(x)_p)
    =
    pmat(mat(Q)_(p-1) mat(R)_(p-1), mat(Q)_(p-1) vect(r) + vect(q) r),
  $) <eq:qr_decomposition>
  By the induction hypothesis,
  there exist a unique choice of matrices $mat(Q)_(p-1)$ and $mat(R)_(p-1)$ with the required structure such that $mat(X)_(p-1) = mat(Q)_(p-1) mat(R)_(p-1)$.
  Comparing the last column of both sides in~@eq:qr_decomposition,
  we obtain
  #numbered($ vect(x)_p = mat(Q)_(p-1) vect(r) + vect(q) r. $) <eq:qr_decomposition_intermediate>
  Left-multiplying both sides by~$mat(Q)_(p-1)^*$ and employing the orthogonality between~$vect(q)$ and the columns of~$mat(Q)_(p-1)$,
  we deduce that necessarily
  $vect(r) = mat(Q)_(p-1)^* vect(x)_p.$
  It then follows from~@eq:qr_decomposition_intermediate that
  $
    vect(q) = 1/r lr((vect(x)_p - mat(Q)_(p-1) mat(Q)_(p-1)^* vect(x)_p)), wide
    r = lr(||vect(x)_p - mat(Q)_(p-1) mat(Q)_(p-1)^* vect(x)_p||).
  $
  It is simple to check that $vect(q)$ is indeed orthogonal to the columns of~$mat(Q)$,
  which concludes the proof.
  Note that~$mat(Q)_(p-1) mat(Q)_(p-1)^* vect(x)_p$ is the orthogonal projection of~$vect(x)_p$ onto the subspace spanned by the columns of~$mat(Q)_(p-1)$.
]
Note that the columns of the matrix~$mat(Q)$ of the decomposition coincide with the vectors that
would be obtained by applying the Gram--Schmidt method to the columns of the matrix~$mat(X)$.
In fact, the Gram--Schmidt process is one of several methods by which the $mat(Q) mat(R)$ decomposition can be calculated in practice.

// [ht!] in LaTeX: floated to the top of the next page
#place(top, float: true)[#algorithm(
  caption: [Simultaneous iteration],
  (0, [$mat(X) <- mat(X)_0$]),
  (0, [*for* $k in {1, 2, dots}$ *do*]),
  (1, [$mat(Q)_k mat(R)_k = mat(A) mat(X)_(k-1)$ (QR decomposition).]),
  (1, [$mat(X)_k <- mat(Q)_k$.]),
  (0, [*end for*]),
) <algo:subspace_iteration>]
The simultaneous iteration method is presented in~@algo:subspace_iteration.
Like the normalization in the power iteration~@algo:power_iteration,
the $mat(Q) mat(R)$ decomposition performed at each step in~@algo:subspace_iteration enables to avoid overflow errors.
Notice that when $p = 1$,
the simultaneous iteration reduces to the power iteration.
We emphasize that
the factorization step at each iteration does not influence the subspace spanned by the columns of~$mat(X)$.
Therefore, this subspace after~$k$ iterations
coincides with that spanned by the columns of the matrix $mat(A)^k mat(X)_0$.
In fact, in exact arithmetic, it would be equivalent to perform the $mat(Q) mat(R)$ decomposition only once as a final step,
after the `for` loop.
Indeed, denoting by $mat(Q)_k mat(R)_k$ the $mat(Q) mat(R)$ decomposition of $mat(A) mat(X)_(k-1)$,
we have
$
  & mat(X)_k = mat(A) mat(X)_(k-1) mat(R)_k^(-1) = mat(A)^2 mat(X)_(k-2) mat(R)_(k-1)^(-1) mat(R)_k^(-1) = cdots
  = mat(A)^k mat(X)_0 mat(R)_1^(-1) dots mat(R)_k^(-1) \
  & <=> quad mat(X)_k (mat(R)_k dots mat(R)_1) = mat(A)^k mat(X)_0.
$
Since $mat(X)_k$ has orthonormal columns and $mat(R)_k dots mat(R)_1$ is an upper triangular matrix (see~@exercise:linear_product_of_lower_triangular)
with real positive elements on the diagonal (check this!),
it follows that~$mat(X)_k$ can be obtained by~$mat(Q) mat(R)$ factorization of $mat(A)^k mat(X)_0$.
In order to show the convergence of the simultaneous iteration,
we begin by proving the following preparatory lemma.
#lemma(title: [Continuity of the reduced QR decomposition])[
  If $mat(Q)_k mat(R)_k -> mat(Q) mat(R)$,
  where $mat(Q) in complex^(n times p)$ has orthonormal columns and~$mat(R) in complex^(p times p)$ is upper triangular with positive real entries on the diagonal,
  then $mat(Q)_k -> mat(Q)$.
] <lemma:continuity_qr>
#proof[
  We reason by contradiction and assume
  there is $epsilon > 0$ and a subsequence $(mat(Q)_(k_n))_(n >= 0)$ such that $norm(mat(Q)_(k_n) - mat(Q)) >= epsilon$ for all~$n$.
  Since the set of matrices with normalized columns is a compact subset of $mat(C)^(n times p)$,
  there exists a further subsequence $(mat(Q)_(k_(n_m)))_(m >= 0)$ that converges to a limit $mat(Q)_infinity$
  which also has orthonormal columns and is at least $epsilon$ away in norm from $mat(Q)$.
  But then
  $
    mat(R)_(k_(n_m))
    = mat(Q)_(k_(n_m))^* (mat(Q)_(k_(n_m)) mat(R)_(k_(n_m)))
    stretch(->)_(m -> infinity) mat(Q)_infinity^* (mat(Q) mat(R)) =: mat(R)_infinity.
  $
  Since $mat(R)_k$ is upper triangular with positive diagonal elements for all $k$,
  clearly $mat(R)_infinity$ is also upper triangular with positive diagonal elements.
  But then $mat(Q)_infinity mat(R)_infinity = mat(Q) mat(R)$,
  and by uniqueness of the decomposition we deduce that $mat(Q) = mat(Q)_infinity$,
  which is a contradiction.
]

#indent Before presenting the convergence theorem,
we introduce the following terminology: we say that $mat(X)_k in complex^(n times p)$ converges essentially to a matrix $mat(X)_infinity$
if each column of $mat(X)_k$ converges essentially to the corresponding column of $mat(X)_infinity$.
We prove the convergence in the Hermitian case for simplicity.
In the general case of~$mat(A) in complex^(n times n)$,
it cannot be expected that~$mat(X)_k$ converges essentially to $mat(V)$,
because the columns of $mat(X)_k$ are orthogonal but eigenvectors may not be orthogonal.
In this case, the columns of $mat(X)_k$ converge not to the eigenvectors
but to the so-called Schur vectors of~$mat(A)$;
see~@MR3396212 for more information.
#theorem(title: [Convergence of the simultaneous iteration~#moreinfo])[
  Assume that $mat(A) in complex^(n times n)$ is Hermitian,
  that~$mat(X)_0 in complex^(n times p)$ has linearly independent columns,
  and finally that the subspace spanned by the column of~$mat(X)_0$ satisfies
  #numbered($
    col(mat(X)_0) inter Span {vect(v)_(p+1), dots, vect(v)_n} = nothing.
  $) <eq:assumption_subspace_iteration>
  If it holds that
  #numbered($
    abs(lambda_1) > abs(lambda_2) > cdots > abs(lambda_p) > abs(lambda_(p+1)) >= abs(lambda_(p+2)) >= dots >= abs(lambda_n),
  $) <eq:order_eigenvalues>
  then $mat(X)_k$ converges _essentially_ to $mat(V)_1 := pmat(vect(v)_1, dots, vect(v)_p)$.
] <theorem:convergence_subspace_iteration>
#proof[
  Let $mat(B) = mat(V)^(-1) mat(X)_0 in complex^(n times p)$, so that $mat(X)_0 = mat(V) mat(B)$,
  and note that $mat(A)^k mat(X)_0 = mat(V) mat(D)^k mat(B)$.
  We denote by~$mat(B)_1 in complex^(p times p)$ and $mat(B)_2 in complex^((n-p) times p)$ the upper $p times p$ and lower $(n-p) times p$ blocks of $mat(B)$,
  respectively.
  The matrix~$mat(B)_1$ is nonsingular,
  otherwise the assumption~@eq:assumption_subspace_iteration would not hold.
  Indeed, if there was a nonzero vector $vect(z) in complex^p$ such that $mat(B)_1 vect(z) = 0$,
  then
  $
    mat(X)_0 vect(z)
    = mat(V) pmat(mat(B)_1; mat(B)_2) vect(z)
    = pmat(mat(V)_1, mat(V)_2) pmat(vect(0); mat(B)_2 vect(z))
    = mat(V)_2 mat(B)_2 vect(z).
  $
  implying that $mat(X)_0 vect(z) in col(mat(X)_0)$ is a linear combination of the vectors in
  $mat(V)_2 = pmat(vect(v)_(p+1), dots vect(v)_n),$
  which contradicts the assumption.
  We also denote by $mat(D)_1$ and $mat(D)_2$ the $p times p$ upper-left and the~$(n-p) times (n-p)$ lower-right blocks of $mat(D)$, respectively.
  From the expression of $mat(A)^k mat(X)_0$,
  we have
  #eq-last($
    mat(A)^k mat(X)_0
    &= pmat(mat(V)_1, mat(V)_2)
    pmat(mat(D)_1^k, ""; "", mat(D)_2^k)
    pmat(mat(B)_1; mat(B)_2)
    = mat(V)_1 mat(D)_1^k mat(B)_1 + mat(V)_2 mat(D)_2^k mat(B)_2, \
    &= lr((mat(V)_1 + mat(V)_2 mat(D)_2^k mat(B)_2 mat(B)_1^(-1) mat(D)_1^(-k))) mat(D)_1^k mat(B)_1.
  $) <eq:factorization>
  The second term in the bracket on the right-hand side converges to zero in the limit as $k -> infinity$ by~@eq:order_eigenvalues.
  Let $tilde(mat(Q)_k) tilde(mat(R)_k)$~denote the reduced $mat(Q) mat(R)$ decomposition of the bracketed term.
  By~@lemma:continuity_qr,
  we deduce from $tilde(mat(Q)_k) tilde(mat(R)_k) -> mat(V)_1$ that $tilde(mat(Q)_k) -> mat(V)_1$.
  Rearranging~@eq:factorization,
  we have
  $ mat(A)^k mat(X)_0 = tilde(mat(Q)_k) (tilde(mat(R)_k) mat(D)_1^k mat(B)_1). $
  Since the matrix between brackets is a $p times p$ square invertible matrix,
  this equation implies that $col(mat(A)^k mat(X)_0) = col(tilde(mat(Q)_k))$.
  Denoting by $mat(Q)_k mat(R)_k$ the $mat(Q) mat(R)$ decomposition of~$mat(A)_k mat(X)_0$,
  we therefore have~$col(mat(Q)_k) = col(tilde(mat(Q)_k))$,
  and so the projectors on these subspaces are equal.
  We recall that, for a set of orthonormal vectors $vect(y)_1, dots, vect(y)_p$
  gathered in a matrix $mat(Y) = pmat(vect(y)_1, dots, vect(y)_p)$,
  the projector on $col(mat(Y)) = Span {vect(y)_1, dots, vect(y)_p} subset complex^n$
  is the square $n times n$ matrix
  $ mat(Y) mat(Y)^* = vect(y)_1 vect(y)_1^* + cdots + vect(y)_p vect(y)_p^*. $
  Consequently,
  the equality of the projectors implies
  $mat(Q)_k mat(Q)_k^* = tilde(mat(Q)_k) tilde(mat(Q)_k)^*,$
  and so $mat(Q)_k mat(Q)_k^* -> mat(V)_1 mat(V)_1^*$.
  Similarly, noting that $mat(Q)_k#tt("[:,1:i]") mat(R)#tt("[1:i,1:i]")$ is the $mat(Q) mat(R)$ decompostion of the matrix $mat(A)^k mat(X)_0#tt("[:,1:i]")$ for all~$i in {1, dots, p}$,
  we obtain the convergence
  #numbered($
    forall i in {1, dots, p}, wide
    mat(Q)_k#tt("[:,1:i]") mat(Q)_k#tt("[:,1:i]")^* stretch(->)_(k -> infinity) vect(v)_1 vect(v)_1^* + cdots + vect(v)_i vect(v)_i^*.
  $) <eq:convergence_projectors>
  This is not surprising given that the first $k$ columns of $mat(X)_0$ undergo a simultaneous iteration independent of the other columns.

  Next, we establish the essential convergence of $mat(Q)_k$ to $mat(V)_1$.
  To this end, we denote the columns of $mat(Q)_k$ by $vect(q)_1^((k)), dots, vect(q)_p^((k))$ and first show by induction that
  $vect(q)_i^((k)) #box($vect(q)_i^((k))$)^* -> vect(v)_i vect(v)_i^*$ in the limit~$k -> infinity$.
  For~$i = 1$ this follows from~@eq:convergence_projectors.
  Assume now that $vect(q)_placeholder^((k)) #box($vect(q)_placeholder^((k))$)^* -> vect(v)_placeholder vect(v)_placeholder^*$ for $placeholder$ up to~$i-1$.
  Then
  $
    vect(q)_i^((k)) #box($vect(q)_i^((k))$)^*
    &= mat(Q)_k#tt("[:,1:i]") mat(Q)_k#tt("[:,1:i]")^* - vect(q)_1^((k)) #box($vect(q)_1^((k))$)^* - cdots - vect(q)_(i-1)^((k)) #box($vect(q)_(i-1)^((k))$)^* \
    &stretch(->)_(k -> infinity) mat(V)_1#tt("[:,1:i]") mat(V)_1#tt("[:,1:i]")^* - vect(v)_1 vect(v)_1^* - cdots - vect(v)_(i-1) vect(v)_(i-1)^* = vect(v)_i vect(v)_i^*.
  $
  It remains to show that the convergence $vect(q)_i^((k)) #box($vect(q)_i^((k))$)^* -> vect(v)_i vect(v)_i^*$ implies the desired essential convergence.
  Noting that~$abs(a) = sqrt(a overline(a))$ for every $a in complex$,
  we have
  $
    abs(vect(v)_i^* vect(q)_i^((k))) = sqrt(vect(v)_i^* vect(q)_i^((k)) #box($vect(q)_i^((k))$)^* vect(v)_i)
    stretch(->)_(k -> infinity) sqrt(vect(v)_i^* vect(v)_i vect(v)_i^* vect(v)_i) = 1,
  $
  Finally, observing that
  $
    lr(||ee^(-ii theta_k) vect(q)_i^((k)) - vect(v)_i||)^2
    = 2 - 2 abs(vect(v)_i^* vect(q)_i^((k))) stretch(->)_(k -> infinity) 0, wide
    ee^(ii theta_k) = (vect(v)_i^* vect(q)_i^((k)))/abs(vect(v)_i^* vect(q)_i^((k))),
  $
  we conclude that $vect(q)_i^((k))$ converges essentially to $vect(v)_i$.
]

#indent In addition to this convergence result,
it is possible to show that the error satisfies
$
  angle Big((col(mat(X)_k), col(mat(V)_1)))
  = bigo lr((lr(|lambda_(p+1)/lambda_p|)^k)).
$
Here,
the angle between two subspaces $cal(A)$ and $cal(B)$ of $complex^n$ is defined as
$
  angle(cal(A), cal(B)) = max_(vect(a) in cal(A) backslash {vect(0)}) bigg((min_(vect(b) in cal(B) backslash {vect(0)}) angle(vect(a), vect(b)))).
$

=== The $mat(Q) mat(R)$ algorithm
The $mat(Q) mat(R)$ algorithm,
which is based on the $mat(Q) mat(R)$ decomposition,
is one of the most famous algorithms for calculating _all_ the eigenpairs of a matrix.
We first present the algorithm and then relate it to the simultaneous iteration in~@sub:simultaneous_iteration.
The method is presented in~@algo:qr.
#algorithm(
  caption: [QR algorithm],
  (0, [$mat(X)_0 = mat(A)$]),
  (0, [*for* $i in {1, 2, dots}$ *do*]),
  (1, [$mat(Q)_k mat(R)_k = mat(X)_(k-1)$ (QR decomposition)]),
  (1, [$mat(X)_k = mat(R)_k mat(Q)_k$]),
  (0, [*end for*]),
) <algo:qr>

#indent Successive iterates of the QR algorithm are related by the equation
#numbered($
  mat(X)_k = mat(Q)_k^(-1) mat(X)_(k-1) mat(Q)_k
  = mat(Q)_k^* mat(X)_(k-1) mat(Q)_k
  = dots = (mat(Q)_1 dots mat(Q)_k)^* mat(X)_0 (mat(Q)_1 dots mat(Q)_k)
$) <eq:qr_iteration_equation>
Therefore, all the iterates are related by a unitary similarity transformation,
and so they all have the same eigenvalues as $mat(X)_0 = mat(A)$.
Rearranging~@eq:qr_iteration_equation,
we have
$ (mat(Q)_1 dots mat(Q)_k) mat(X)_k = mat(A) (mat(Q)_1 dots mat(Q)_k), $
and so, introducing $tilde(mat(Q))_k = mat(Q)_1 dots mat(Q)_k$
and noting that $mat(X)_k = mat(Q)_(k+1) mat(R)_(k+1)$ by the algorithm,
we deduce
$ tilde(mat(Q))_(k+1) mat(R)_(k+1) = mat(A) tilde(mat(Q))_k. $
This reveals that
the matrix sequence $(tilde(mat(Q))_k)_(k >= 1)$ undergoes a simultaneous iteration and so,
assuming that $mat(A)$ is Hermitian with $n$~distinct nonzero eigenvalues,
we deduce that $tilde(mat(Q))_k -> mat(V)$ essentially in the limit as $k -> infinity$,
by~@theorem:convergence_subspace_iteration.
As a consequence, by~@eq:qr_iteration_equation,
it holds that~$mat(X)_k -> mat(V)^* mat(X)_0 mat(V) = mat(D)$;
in other words, the matrix $mat(X)_k$ converges to a diagonal matrix with the eigenvalues of~$mat(A)$ on the diagonal.

== Projection methods <sec:projection_methods>
In this section,
we begin by presenting a general method for constructing an approximation of the eigenvectors of~$mat(A)$ in a given subspace~$cal(U)$ of~$complex^n$.
We then discuss a particular choice for the subspace $cal(U)$
as a Krylov subspace, which is very useful in practice.

Assume that ${vect(u)_1, dots, vect(u)_p}$ is an orthonormal basis of $cal(U)$.
Then for any vector $vect(v) in complex^n$,
the vector of $cal(U)$ that is closest to $vect(v)$ in the Euclidean distance is given by the orthogonal projection
$ mat(P)_cal(U) vect(v) := mat(U) mat(U)^* vect(v) = (vect(u)_1 vect(u)_1^* + cdots + vect(u)_p vect(u)_p^*) vect(v). $
In practice,
the eigenvectors of~$mat(A)$ are unknown,
and so it is impossible to calculate approximations using this formula.
The Rayleigh--Ritz method,
which we present hereafter,
is an alternative and practical method for constructing approximations of the eigenvectors and eigenvalues.
In general, the subspace $cal(U)$ does not contain any eigenvector of $mat(A)$,
and so the problem
#numbered($ mat(A) vect(v) = lambda vect(v), wide vect(v) in cal(U) $) <eq:eigenproblem_subspace>
does not admit a solution.
Let us denote by $mat(U)$ the matrix with columns~$vect(u)_1, dots, vect(u)_p$.
Since any vector $vect(v) in cal(U)$ is equal to $mat(U) vect(z)$ for some vector $vect(z) in complex^p$,
equation~@eq:eigenproblem_subspace is equivalent to the problem
$ mat(A) mat(U) vect(z) = lambda mat(U) vect(z), $
which is a system of $n$ equations with $p < n$ unknowns.
The Rayleigh--Ritz method is based on the idea that,
in order to obtain a problem with as many unknowns as there are equations,
we can multiply this equation by~$mat(U)^*$,
which leads to the problem
#numbered($ mat(B) vect(z) := (mat(U)^* mat(A) mat(U)) vect(z) = lambda vect(z). $) <eq:rayleigh_ritz>
This is standard eigenvalue problem for the matrix~$mat(U)^* mat(A) mat(U) in complex^(p times p)$,
which is much easier to solve than the original problem if $p << n$.
Equivalently, equation~@eq:rayleigh_ritz may be formulated as follows: find $vect(v) in cal(U)$
such that
#numbered($ vect(u)^* (mat(A) vect(v) - lambda vect(v)), wide forall vect(u) in cal(U). $) <eq:equivalent_rayleigh_ritz>
The solutions to~@eq:rayleigh_ritz and~@eq:equivalent_rayleigh_ritz are related by the equation~$vect(v) = mat(U) vect(z)$.
Of course, the eigenvalues of~$mat(B)$ in problem~@eq:rayleigh_ritz,
which are called the Ritz values of~$mat(A)$ relative to $cal(U)$,
are in general different from those of~$mat(A)$.
Once an eigenvector $vect(y)$ of~$mat(B)$ has been calculated,
an approximate eigenvector of~$mat(A)$,
called a _Ritz vector_ of~$mat(A)$ relative to~$cal(U)$,
is obtained from the equation $hat(vect(v)) = mat(U) vect(y)$.
The Rayleigh--Ritz algorithm is presented in full in~@algo:rayleigh_ritz.
#algorithm(
  caption: [Rayleigh--Ritz],
  (0, [Choose $cal(U) subset complex^n$]),
  (0, [Construct a matrix~$mat(U)$ whose columns are orthonormal and span $cal(U)$]),
  (0, [Find the eigenvalues $hat(lambda)_i$ and eigenvectors $vect(y)_i in complex^p$ of $mat(B) := mat(U)^* mat(A) mat(U)$]),
  (0, [Calculate the corresponding Ritz vectors $hat(vect(v))_i = mat(U) vect(y)_i in complex^n$.]),
) <algo:rayleigh_ritz>

#indent It is clear that if $vect(v)_i in cal(U)$,
then~$lambda_i$ is an eigenvalue of $mat(B)$ in~@eq:rayleigh_ritz.
In fact,
we can show the following more general statement.
#proposition[
  If $cal(U)$ is an invariant subspace of~$mat(A)$,
  meaning that $mat(A) cal(U) subset cal(U)$,
  then each Ritz vector of~$mat(A)$ relative to~$cal(U)$ is an eigenvector of~$mat(A)$.
] <proposition:invariant_subpsace_ritz>
#proof[
  Let $mat(U) in complex^(n times p)$ and $mat(W) in complex^(n times (n-p))$ be matrices whose columns form orthonormal bases of~$cal(U)$ and~$cal(U)^perp$,
  respectively.
  Here~$cal(U)^perp$ denotes the orthogonal complement of~$cal(U)$ with respect to the Euclidean inner product.
  Then, since~$mat(W)^* mat(A) mat(U) = mat(0)$ by assumption,
  it holds that
  $
    mat(Q)^* mat(A) mat(Q) =
    pmat(
      mat(U)^* mat(A) mat(U), mat(U)^* mat(A) mat(W);
      mat(W)^* mat(A) mat(U), mat(W)^* mat(A) mat(W)
    )
    =
    pmat(
      mat(U)^* mat(A) mat(U), mat(U)^* mat(A) mat(W);
      mat(0), mat(W)^* mat(A) mat(W)
    ),
    wide
    mat(Q) = pmat(mat(U), mat(W)).
  $
  If $(vect(y), hat(lambda))$ is an eigenvector of $mat(U)^* mat(A) mat(U)$,
  then
  $
    mat(Q)^* mat(A) mat(Q) pmat(vect(y); vect(0))
    =
    pmat((mat(U)^* mat(A) mat(U)) vect(y); vect(0))
    =
    hat(lambda) pmat(vect(y); vect(0))
    =:
    hat(lambda) vect(x),
  $
  and so $(vect(x), hat(lambda))$ is an eigenpair of $mat(Q)^* mat(A) mat(Q)$.
  But then $(mat(Q) vect(x), hat(lambda)) = (mat(U) vect(y), hat(lambda))$ is an eigenpair of~$mat(A)$,
  which proves the statement.
]
#remark[
  @proposition:invariant_subpsace_ritz can be proved more directly from~@eq:equivalent_rayleigh_ritz by noting that,
  if~$cal(U)$ is an invariant subspace of~$mat(A)$,
  then this equation implies that $mat(A) vect(v) - lambda vect(v)$ belongs to both $cal(U)$ and its orthogonal complement $cal(U)^perp$,
  and so this vector is $vect(0)$.
]

#indent If $cal(U)$ is close to being an invariant subspace of~$mat(A)$,
then it is expected that the Ritz vectors and Ritz values of~$mat(A)$ relative to~$cal(U)$ will provide good approximations of some of the eigenpairs of~$mat(A)$.
Quantifying this approximation is difficult,
so we only present without proof the following error bound.
See~@MR1990645 for more information.
#proposition[
  Let~$mat(A)$ be a full rank Hermitian matrix and $cal(U)$ a $p$-dimensional subspace of~$complex^n$.
  Then there exists eigenvalues $lambda_(i_1), dots, lambda_(i_p)$ of~$mat(A)$ which satisfy
  $
    forall j in {1, dots, p}, wide
    abs(lambda_(i_j) - hat(lambda)_j) <= norm((mat(I) - mat(P)_cal(U)) mat(A) mat(P)_cal(U))_2.
  $
]

#indent In the case where~$mat(A)$ is Hermitian,
it is possible to show that the Ritz values are bounded from above by the eigenvalues of~$mat(A)$.
The proof of this result relies on the Courant--Fisher theorem for characterizing the eigenvalues of a Hermitian matrix,
which is recalled in~@theorem:courant-fisher in the appendix.
#proposition[
  If $mat(A) in complex^(n times n)$ is Hermitian,
  then
  $
    forall i in {1, dots, p}, wide
    hat(lambda)_i <= lambda_i
  $
]
#proof[
  By the Courant--Fisher theorem,
  it holds that
  $
    hat(lambda)_i
    = max_(cal(S) subset complex^p, "dim"(cal(S)) = i) lr((min_(vect(x) in cal(S) backslash {0}) (vect(x)^* mat(B) vect(x))/(vect(x)^* vect(x))))
  $
  Letting $vect(y) = mat(U) vect(x)$ and then $cal(R) = mat(U) cal(S)$,
  we deduce that
  $
    hat(lambda)_i
    &= max_(cal(S) subset complex^p, "dim"(cal(S)) = i) lr((min_(vect(y) in mat(U) cal(S) backslash {0}) (vect(y)^* mat(A) vect(y))/(vect(y)^* vect(y)))) \
    &= max_(cal(R) subset cal(U), "dim"(cal(R)) = i) lr((min_(vect(y) in cal(R) backslash {0}) (vect(y)^* mat(A) vect(y))/(vect(y)^* vect(y))))
    <= max_(cal(R) subset complex^n, "dim"(cal(R)) = i) lr((min_(vect(y) in cal(R) backslash {0}) (vect(y)^* mat(A) vect(y))/(vect(y)^* vect(y)))) = lambda_i,
  $
  where we used the Courant--Fisher theorem for the matrix~$mat(A)$ in the last equality.
]

#indent This projection approach is sometimes combined with a simultaneous subspace iteration:
an approximation $mat(X)_k$ of the $p$ first eigenvector is first calculated using~@algo:subspace_iteration,
and then the matrix $mat(X)_k$ is used in place of~$mat(U)$ in~@algo:rayleigh_ritz.

=== Projection method in a Krylov subspace
The power iteration constructs at iteration $k$ an approximation of $vect(v)_1$ in the one-dimensional subspace spanned by the vector $mat(A)^k vect(x)_0$,
and only the previous iteration $vect(x)_k$ is employed to construct $vect(x)^(k+1)$.
One may wonder whether,
by employing all the previous iterates rather than only the previous one,
a better approximation of~$vect(v)_1$ can be constructed.
More precisely, instead of looking for an approximation in the subspace $Span {mat(A)^k vect(x)_0}$,
would it be useful to extend the search area to the Krylov subspace
$
  cal(K)_(k+1)(mat(A), vect(x)_0) := Span big(\{) vect(x)_0, mat(A) vect(x)_0, dots, mat(A)^k vect(x)_0 big(\})?
$
The answer to this question is positive,
and the resulting method is often much faster than the power iteration.
This is achieved by employing the Rayleigh--Ritz projection method~@algo:rayleigh_ritz with the choice $cal(U) = cal(K)_(k+1)(mat(A), vect(x)_0)$.
Applying this method requires to calculate an orthonormal basis of the Krylov subspace and to calculate the reduced matrix $mat(U)^* mat(A) mat(U)$.
The _Arnoldi method_ enables to achieve these two goals simultaneously.

=== The Arnoldi iteration
This Arnoldi iteration is based on the Gram--Schmidt process and presented in~@algo:arnoldi_method.
#algorithm(
  caption: [Arnoldi iteration for constructing an orthonormal basis of $cal(K)_p (mat(A), vect(u)_1)$],
  (0, [Choose $vect(u)_1$ with unit norm.]),
  (0, [*for* $j in {1, dots p}$ *do*]),
  (1, [$vect(u)_(j+1) <- mat(A) vect(u)_j$]),
  (1, [*for* $i in {1, dots, j}$ *do*]),
  (2, [$h_(i,j) <- vect(u)_i^* vect(u)_(j+1)$]),
  (2, [$vect(u)_(j+1) <- vect(u)_(j+1) - h_(i,j) vect(u)_i$]),
  (1, [*end for*]),
  (1, [$h_(j+1,j) <- norm(vect(u)_(j+1))$]),
  (1, [$vect(u)_(j+1) <- vect(u)_(j+1) \/ h_(j+1,j)$]),
  (0, [*end for*]),
) <algo:arnoldi_method>
The iteration breaks down if $h_(j+1,j) = 0$,
which indicates that $cal(A) vect(u)_j$ belongs to the Krylov subspace $Span {vect(u)_1, dots, vect(u)_j} = cal(K)_j (mat(A), vect(u)_1)$,
implying that $cal(K)_(j+1)(mat(A), vect(u)_1) = cal(K)_j (mat(A), vect(u)_1)$.
In this case, the subspace $cal(K)_j (mat(A), vect(u)_1)$ is an invariant subspace of~$mat(A)$ because, by~@exercise:krylov,
we have
$
  mat(A) cal(K)_j (mat(A), vect(u)_1)
  subset cal(K)_(j+1)(mat(A), vect(u)_1)
  = cal(K)_j (mat(A), vect(u)_1).
$
Therefore, applying the Rayleigh--Ritz with $cal(U) = Span {vect(u)_1, dots, vect(u)_j}$ yields exact eigenpairs in view of~@proposition:invariant_subpsace_ritz.
If the iteration does not break down then,
by construction,
the vectors ${vect(u)_1, dots, vect(u)_p}$ at the end of the algorithm are orthonormal.
It is also simple to show by induction that they form a basis of~$cal(K)_p (mat(A), vect(u)_1)$.
The scalar coefficients $h_(i,j)$ can be collected in a matrix square $p times p$ matrix
$
  mat(H)
  =
  pmat(
    h_(1,1), h_(1,2), h_(1,3), cdots, h_(1,p);
    h_(2,1), h_(2,2), h_(2,3), cdots, h_(2,p);
    0, h_(3,2), h_(3,3), cdots, h_(3,p);
    vdots, ddots, ddots, ddots, vdots;
    0, cdots, 0, h_(p,p-1), h_(p,p)
  ).
$
This matrix contains only zeros under the first subdiagonal;
such a matrix is called a~_Hessenberg_ matrix.
Inspecting the algorithm,
we notice that the $j$-th column contains the coefficients of the projection of the vector~$mat(A) vect(u)_j$ onto the basis ${vect(u)_1, dots, vect(u)_p}$.
In other words,
#numbered($ mat(U)^* mat(A) mat(U) = mat(H), $) <eq:reduced_matrix_arnoldi>
We have thus shown that the Arnoldi algorithm enables to construct both an orthonormal basis of a Krylov subspace and
the associated reduced matrix.
In fact, we have the following equation
#numbered($
  mat(A) mat(U) = mat(U) mat(H) + h_(p+1,p) (vect(v)_(p+1) vect(e)_p^*),
  wide vect(e)_p = pmat(0; vdots; 1) in complex^p.
$) <eq:arnoldi_equation>
The Arnoldi algorithm, coupled with the Rayleigh--Ritz method,
has very good convergence properties in the limit as $p -> infinity$,
in particular for eigenvalues with a large modulus.
The following result shows that the residual $vect(r) = mat(A) hat(vect(v)) - hat(lambda) vect(v)$
associated with a Ritz vector can be estimated inexpensively.
Specifically, the norm of the residual is equal to the last component of the associated eigenvector of~$mat(H)$ multiplied by $h_(p+1,p)$.
#proposition(title: [Formula for the residual~#moreinfo])[
  Let $vect(y)_i$ be an eigenvector of $mat(H)$ associated with the eigenvalues~$hat(lambda)_i$,
  and let $hat(vect(v))_i = mat(U) vect(y)_i$ denote the corresponding eigenvector.
  Then
  $ mat(A) hat(vect(v))_i - hat(lambda) vect(v)_i = h_(p+1,p) (vect(y)_i)_p vect(v)_(p+1). $
  Consequently, it holds that
  $ norm(mat(A) hat(vect(v))_i - hat(lambda) vect(v)_i) = abs(h_(p+1,p) (vect(y)_i)_p). $
]
#proof[
  Multiplying both sides of~@eq:arnoldi_equation by $vect(y)_i$,
  we obtain
  $ mat(A) mat(U) vect(y)_i = mat(U) mat(H) vect(y)_i + h_(p+1,p) (vect(v)_(p+1) vect(e)_p^*) vect(y)_i. $
  Using the definition of $hat(vect(v))_i$ and rearranging the equation,
  we have
  $ mat(A) hat(vect(v)_i) - hat(lambda)_i vect(y)_i = h_(p+1,p) (vect(v)_(p+1) vect(e)_p^*) vect(y)_i, $
  which immediately gives the result.
]

#indent In practice, the larger the dimension~$p$ of the subspace~$cal(U)$ employed in the Rayleigh--Ritz method,
the more memory is required for storing an orthonormal basis of~$cal(U)$.
In addition, for large values of~$p$,
computing the reduced matrix~@eq:reduced_matrix_arnoldi and its eigenpairs becomes computationally expensive;
the computational cost of computing the matrix~$mat(H)$ scales as $bigo(p^2)$.
To remedy these potential issues,
the algorithm can be restarted periodically.
For example, @algo:arnoldi_power can be employed as an alternative to the power iteration
in order to find the eigenvector associated with the eigenvalue with largest modulus.

#algorithm(
  caption: [Restarted Arnoldi iteration],
  (0, [Choose $vect(u)_1 in complex^n$ and $p << n$]),
  (0, [*for* $i in {1, 2, dots}$ *do*]),
  (1, [Perform $p$ iterations of the Arnoldi iteration and construct~$cal(U)$;]),
  (1, [Calculate the Ritz vector $hat(vect(v))_1$ associated with the largest Ritz value relative to~$cal(U)$;]),
  (1, [If this vector is sufficiently accurate, then stop. Otherwise, restart with $vect(u)_1 = hat(vect(v))_1$.]),
  (0, [*end for*]),
) <algo:arnoldi_power>

=== The Lanczos iteration
The Lanczos iteration may be viewed as a simplified version of the Arnoldi iteration in the case where the matrix~$mat(A)$ is Hermitian.
Let us denote by ${vect(u)_1, dots, vect(u)_p}$ the orthonormal vectors generated by the Arnoldi iteration.
When~$mat(A)$ is Hermitian, it holds that
$ h_(i,j) = vect(u)_i^* (mat(A) vect(u)_j) = (mat(A) vect(u)_i)^* vect(u)_j = overline(h_(j,i)). $
Therefore, the matrix $mat(H)$ is Hermitian.
This is not surprising, since we showed that $mat(H) = mat(U)^* mat(A) mat(U)$ and the matrix $mat(A)$ is Hermitian.
Since~$mat(H)$ is also of Hessenberg form,
we deduce that~$mat(H)$ is tridiagonal.
An inspection of~@algo:arnoldi_method shows that the subdiagonal entries of~$mat(H)$ are real.
Since $mat(A)$ is Hermitian, the diagonal entries $h_(i,i) = vect(u)_i^* (mat(A) vect(u)_j)$ are also real,
and so we conclude that all the entries of the matrix~$mat(H)$ are in fact real.
This matrix if of the form
$
  mat(H)
  =
  pmat(
    alpha_1, beta_2, "", "", "";
    beta_2, alpha_2, beta_3, "", "";
    "", beta_3, ddots, ddots, "";
    "", "", ddots, ddots, beta_p;
    "", "", "", beta_p, alpha_p
  )
$
Adapting the Arnoldi iteration to this setting leads to~@algo:lanczos_method.
#algorithm(
  caption: [Lanczos iteration for constructing an orthonormal basis of $cal(K)_p (mat(A), vect(u)_1)$],
  (0, [Choose $vect(u)_1$ with unit norm.]),
  (0, [$beta_1 <- 0$, $vect(u)_0 <- vect(0) in complex^n$]),
  (0, [*for* $j in {1, dots p}$ *do*]),
  (1, [$vect(u)_(j+1) <- mat(A) vect(u)_j - beta_j vect(u)_(j-1)$]),
  (1, [$alpha_j <- vect(u)_j^* vect(u)_(j+1)$]),
  (1, [$vect(u)_(j+1) <- vect(u)_(j+1) - alpha_j vect(u)_j$]),
  (1, [$beta_(j+1) <- norm(vect(u)_(j+1))$]),
  (1, [$vect(u)_(j+1) <- vect(u)_(j+1) \/ beta_(j+1)$]),
  (0, [*end for*]),
) <algo:lanczos_method>

== Exercises
#exercise[
  PageRank is an algorithm for assigning a rank to the vertices of a directed graph.
  It is used by many search engines, notably Google,
  for sorting search results.
  In this context,
  the directed graph encodes the links between pages of the World Wide Web:
  the vertices of the directed graph are webpages,
  and there is an edge going from page~$i$ to page~$j$ if page~$i$ contains a hyperlink to page~$j$.

  Let us consider a directed graph $G(V, E)$ with vertices $V = {1, dots, n}$ and edges $E$.
  The graph can be represented by its adjacency matrix $mat(A) in {0, 1}^(n times n)$,
  whose entries are given by
  $
    a_(i j) =
    cases(
      1 & quad #emph[if there is an edge from $i$ to $j$,],
      0 & quad #emph[otherwise.]
    )
  $
  Let $r_i$ denote the "value" assigned to vertex $i$.
  The idea of PageRank, in its simplest form,
  is to assign values to the vertices
  by solving the following system of equations;
  #numbered($
    forall i in V, wide
    r_i
    = sum_(j in cal(N)(i)) r_j/o_j.
  $) <eq:equation_page_rank>
  where $o_j$ is the outdegree of vertex $j$,
  i.e.~the number of edges leaving from~$j$.
  Here the sum is over the set of nodes $cal(N)(i)$, which denotes all the "incoming" neighbors of $i$,
  i.e.~those that have an edge pointing towards node $i$.
  - Read the Wikipedia page on _PageRank_ to familiarize yourself with the algorithm.

  - Let~$vect(r) = pmat(r_1, dots, r_n)^T$.
    Show using~@eq:equation_page_rank that $vect(r)$ satisfies
    $
      vect(r) =
      mat(A)^T
      pmat(
        1/o_1, "", "";
        "", ddots, "";
        "", "", 1/o_n
      )
      vect(r) =: mat(A)^T mat(O)^(-1) vect(r).
    $
    In other words, $vect(r)$ is an eigenvector with eigenvalue $1$ of the matrix $mat(M) = mat(A)^T mat(O)^(-1)$.

  - Show that $mat(M)$ is a left-stochastic matrix,
    i.e.~that each column sums to 1.

  - Prove that the eigenvalues of any matrix $mat(B) in real^(n times n)$ coincide with those of $mat(B)^T$.
    You may use the fact that $det(mat(B)) = det(mat(B)^T)$.

  - Using the previous items,
    show that 1 is an eigenvalue and that $rho(mat(M)) = 1$.
    For the second part, find a subordinate matrix norm such that $norm(mat(M)) = 1$.

  - Implement PageRank in order to rank pages from a 2013 snapshot of English Wikipedia.
    You can use either the simplified version of the algorithm given in~@eq:equation_page_rank or the improved version with a damping factor described on Wikipedia.
    In the former case,
    the following are both sensible stopping criteria:
    $
      norm(mat(M) hat(vect(r)) - hat(vect(r)))_1/norm(hat(vect(r)))_1 < 10^(-15)
      wide " or " wide
      norm(mat(M) hat(vect(r)) - hat(lambda) hat(vect(r)))_1/norm(vect(r))_1 < 10^(-15),
      wide hat(lambda) = (hat(vect(r))^T mat(M) hat(vect(r)))/(hat(vect(r))^T hat(vect(r))),
    $
    where $hat(vect(v))$ is an approximation of the eigenvector corresponding to the dominant eigenvalue.
    A dataset is available on the course website to complete this part.
    This dataset contains a subset of the data publicly available #link("https://snap.stanford.edu/data/enwiki-2013.html")[here],
    and was generated from the full dataset by retaining only the 5% best rated articles.
    After decompressing the archive,
    you can load the dataset into Julia by using the following commands:

    ```julia
    import CSV
    import DataFrames

    # Data (nodes and edges)
    nodes = CSV.read("names.csv", DataFrames.DataFrame)
    edges = CSV.read("edges.csv", DataFrames.DataFrame)

    # Convert data to matrices
    nodes = Matrix(nodes)
    edges = Matrix(edges)
    ```

    After you have assigned a rank to all the pages,
    print the 10 pages with the highest ranks.
    My code returns the following entries:
    // multicols{3}: columns balanced as in LaTeX
    #let item(n, body) = [#box(width: 2em, align(right)[#n.])#h(0.5em)#body]
    #block(above: leading + 9pt, below: leading + 9pt, grid(
      columns: (1fr, 1fr, 1fr), row-gutter: leading,
      item(1)[United States], item(5)[France], item(9)[Canada],
      item(2)[United Kingdom], item(6)[Germany], item(10)[India],
      item(3)[World War~II], item(7)[English language], [],
      item(4)[Latin], item(8)[China], [],
    ))
  - *Extra credit:*
    Write a function `search(keyword)` that can be employed for searching the database.
    Here is an example of what it could return:

    ```julia
    julia> search("New York")
    481-element Vector{String}:
     "New York City"
     "New York"
     "The New York Times"
     "New York Stock Exchange"
     "New York University"
     …
    ```
]

#exercise[
  Show the following properties of the Krylov subspace $cal(K)_p (mat(A), vect(x))$.
  - $cal(K)_p (mat(A), vect(x)) subset cal(K)_(p+1)(mat(A), vect(x)).$

  - $mat(A) cal(K)_p (mat(A), vect(x)) subset cal(K)_(p+1)(mat(A), vect(x)).$

  - The Krylov subspace $cal(K)_p (mat(A), vect(x))$ is invariant under rescaling:
    for all $alpha in complex$,
    $ cal(K)_p (mat(A), vect(x)) = cal(K)_p (alpha mat(A), vect(x)) = cal(K)_p (mat(A), alpha vect(x)). $

  - The Krylov subspace $cal(K)_p (mat(A), vect(x))$ is invariant under shift of the matrix~$mat(A)$:
    for all $alpha in complex$,
    $ cal(K)_p (mat(A), vect(x)) = cal(K)_p (mat(A) - alpha mat(I), vect(x)). $

  - Similarity transformation: If $mat(T) in complex^(n times n)$ is nonsingular,
    then
    $
      cal(K)_p (mat(T)^(-1) mat(A) mat(T), mat(T)^(-1) vect(x)) =
      mat(T)^(-1) cal(K)_p (mat(A), vect(x)).
    $
] <exercise:krylov>

#exercise[
  The minimal polynomial of a matrix $mat(A) in complex^(n times n)$ is the monic polynomial~$p$ of lowest degree such that $p(mat(A)) = 0$.
  Prove that, if~$mat(A)$ is Hermitian with $m <= n$ distinct eigenvalues,
  then the minimal polynomial is given by
  $ p(t) = product_(i=1)^m (t - lambda_i). $
]

#exercise[
  The minimal polynomial for a general matrix $mat(A) in complex^(n times n)$ is given by
  $ p(t) = product_(i=1)^m (t - lambda_i)^(s_i). $
  where $s_i$ is the size of the largest Jordan block associated with the eigenvalue~$lambda_i$ in the normal Jordan form of~$mat(A)$.
  Verify that $p(mat(A)) = 0$.
]

#exercise[
  Let $d$ denote the degree of the minimal polynomial of~$mat(A)$.
  Show that
  $
    forall p >= d, wide
    cal(K)_(p+1)(mat(A), vect(x)) =
    cal(K)_p (mat(A), vect(x)).
  $
  Deduce that, for $p >= n$,
  the subspace $cal(K)_p (mat(A), vect(x))$ is an invariant subspace of~$mat(A)$.
]

#exercise[
  Let $mat(A) in complex^(n times n)$.
  Show that $cal(K)_n (mat(A), vect(x))$ is the smallest invariant subspace of~$mat(A)$ that contains $vect(x)$.
]

#compexercise[
  Consider the matrix
  $
    mat(M) =
    pmat(
      0, 1, 2, 0;
      1, 0, 1, 0;
      2, 1, 0, 2;
      0, 0, 2, 0
    )
  $
  - Find the dominant eigenvalue of~$mat(M)$ by using the power iteration.
  - Find the eigenvalue of~$mat(M)$ closest to 1 by using the inverse iteration.
  - Find the other two eigenvalues of~$mat(M)$ by using a method of your choice.
]

#exercise(title: [A posteriori error bound])[
  Assume that $mat(A) in complex^(n times n)$ is Hermitian,
  and that $hat(vect(v))$ is a normalized approximation of an eigenvector which satisfies
  $
    norm(hat(vect(z))) := norm(mat(A) hat(vect(v)) - hat(lambda) hat(vect(v))) = delta,
    wide hat(lambda) = (hat(vect(v))^* mat(A) hat(vect(v)))/(hat(vect(v))^* hat(vect(v))).
  $
  Prove that there is an eigenvalue~$lambda$ of~$mat(A)$ such that
  $ abs(hat(lambda) - lambda) <= delta. $

  *Hint:* Consider first the case where $mat(A)$ is diagonal.
] <exercise:a_posteriori_1>

#exercise(title: [Bauer--Fike theorem])[
  Assume that $mat(A) in complex^(n times n)$ is diagonalizable: $mat(A) mat(V) = mat(V) mat(D)$.
  Show that, if $hat(vect(v))$ is a normalized approximation of an eigenvector which satisfies
  $ norm(vect(r)) := norm(mat(A) hat(vect(v)) - hat(lambda) hat(vect(v))) = delta $
  for some $hat(lambda) in complex$,
  then there is an eigenvalue~$lambda$ of~$mat(A)$ such that
  $ abs(hat(lambda) - lambda) <= kappa_2(mat(V)) delta. $

  *Hint*:
  Rewrite
  $
    norm(hat(vect(v))) = norm((mat(A) - hat(lambda) mat(I))^(-1) vect(r))
    = norm(mat(V) (mat(D) - hat(lambda) mat(I))^(-1) mat(V)^(-1) vect(r))).
  $
] <exercise:a_posteriori_2>

#exercise[
  In @exercise:a_posteriori_1 and @exercise:a_posteriori_2,
  we saw examples~_a posteriori_ error estimates
  which guarantee the existence of an eigenvalue of~$mat(A)$ within a certain distance of the approximation~$hat(lambda)$.
  In this exercise,
  our aim is to provide an answer to the following different question:
  given an approximate eigenpair $(hat(vect(v)), hat(lambda))$,
  what is the smallest perturbation~$mat(E)$ that we need to apply to $mat(A)$ in order to guarantee that
  $(hat(vect(v)), hat(lambda))$ is an exact eigenpair, i.e.~that
  $ (mat(A) + mat(E)) hat(vect(v)) = hat(lambda) hat(vect(v))? $
  Assume that $hat(vect(v))$ is normalized and
  let $cal(E) = Big(\{) mat(E) in complex^(n times n): (mat(A) + mat(E)) hat(vect(v)) = hat(lambda) hat(vect(v)) Big(\})$.
  Prove that
  #numbered($
    min_(mat(E) in cal(E)) norm(mat(E))_2 = norm(vect(r))_2 := norm(mat(A) hat(vect(v)) - hat(lambda) hat(vect(v))).
  $) <eq:kahan_parlett_jiang>
  To this end,
  you may proceed as follows:
  - Show first that any $mat(E) in cal(E)$ satisfies $mat(E) hat(vect(v)) = - vect(r)$.

  - Deduce from the first item that
    $ inf_(mat(E) in cal(E)) norm(mat(E))_2 >= norm(vect(r))_2. $

  - Find a rank one matrix $mat(E)_*$ such that
    $norm(mat(E)_*)_2 = norm(vect(r))_2,$
    and then conclude.
    Recall that any rank 1 matrix can be written in the form $mat(E)_* = vect(u) vect(w)^*$,
    with norm $norm(vect(u))_2 norm(vect(w))_2$.
  Equation~@eq:kahan_parlett_jiang is a simplified version of the Kahan--Parlett--Jiang theorem
  and is an example of a _backward error estimate_.
  Whereas forward error estimates quantify the distance between an approximation and the exact solution,
  backward error estimates give the size of the perturbation that must be applied to a problem so that an approximation is exact.
]

#exercise[
  Assume that $(vect(x)_k)_(k >= 0)$ is a sequence of normalized vectors in $complex^n$.
  Show that the following statements are equivalent:
  - $(vect(x)_k)_(k >= 0)$ converges essentially to $vect(x)_infinity$ in the limit as $k -> infinity$.
  - The angle $angle(vect(x)_k, vect(x)_infinity)$ converges to zero in the limit as $k -> infinity$.
  - The projector $mat(P)_(vect(x)_k)$ converges to $mat(P)_(vect(x)_infinity)$ in the limit as $k -> infinity$.
] <exercise:essential_convergence>

#exercise[
  Assume that $mat(A) in complex^(n times n)$ is skew-Hermitian.
  Derive a Lanczos-like algorithm for constructing an orthonormal basis of $cal(K)_p (mat(A), vect(x))$
  and calculating the reduced matrix
  $ mat(U)^* mat(A) mat(U), $
  where $mat(U) in complex^(n times p)$ contains the vectors of the basis as columns.
  Implement your algorithm with~$p = 20$ in order to approximate the dominant eigenvalue of the matrix~$mat(A)$ constructed by the following piece of code:
  ```julia
  n = 5000
  A = rand(n, n) + im * randn(n, n)
  A = A - A'  # A is now skew-Hermitian
  ```
]

#exercise[
  Assume that ${vect(u)_1, dots, vect(u)_p}$ and ${vect(w)_1, dots, vect(w)_n}$
  are orthonormal bases of the same subspace $cal(U) subset complex^n$.
  Show that the projectors $mat(U) mat(U)^*$ and $mat(W) mat(W)^*$ are equal.
]

#exercise[
  Assume that $mat(A) in complex^(n times n)$ is a Hermitian matrix with distinct eigenvalues,
  and suppose that we know the dominant eigenpair $(vect(v)_1, lambda_1)$,
  with $vect(v)_1$ normalized.
  Let
  $ mat(B) = mat(A) - lambda_1 vect(v)_1 vect(v)_1^*. $
  If we apply the power iteration to this matrix,
  what convergence can we expect?
]

#exercise[
  Assume that $hat(vect(v)_1)$ and $hat(vect(v)_2)$ are two Ritz vectors of a Hermitian matrix~$mat(A)$ relative to a subspace~$cal(U) subset complex^n$.
  Show that, if the associated Ritz values are distinct,
  then $hat(vect(v)_1)^* hat(vect(v)_2) = 0$.
]

== Discussion and bibliography

The content of this chapter is inspired mainly from~@VanDooren and also from~@MR3396212.
The latter volume contains a detailed coverage of the standard methods for eigenvalue problems.
Some of the exercises are taken from~@Vuik,
and others are inspired from~@MR3396212.

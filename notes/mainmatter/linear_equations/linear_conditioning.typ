#import "../../template.typ": *

// \begin{psmallmatrix} ... \end{psmallmatrix}
#let psmallmat(..rows) = $lr((script(#math.mat(delim: none, row-gap: 0.2em, column-gap: 0.4em, ..rows))))$

== Conditioning <sec:conditioning>
Before studying the properties of numerical methods for solving linear system,
it is crucial to understand the impact of _round-off errors_ on the solution,
which impose a bound on the accuracy we can hope to achieve.
We begin with a motivating example.
#example[
  Suppose that we wish to solve the following equation
  $
    mat(A) vect(x) :=
    pmat(1, 1; 1, 1 - 10^(-12))
    pmat(x_1; x_2)
    =
    pmat(0; 10^(-12)) =: vect(b).
  $
  The exact solution is given by $psmallmat(1, -1)^T$.
  In Julia, this equation can be solved as follows:
  ```julia
  A = [1 1; 1 (1-1e-12)]
  b = [0; 1e-12]
  x = A\b
  ```
  The solution returned by the program is the following:
  #raw("  1.0000221222095027\n -1.0000221222095027", block: true)
  The relative error on the solution is of the order of $10^(-5)$,
  which is about 12 orders of magnitude larger than the machine epsilon for the `Float64` type.
  In order to understand this,
  note that the final error on~`x` arises from three sources:
  - First, the machine representation `b` of the right-hand side is only _an approximation_ of the true right-hand side~$vect(b)$.

  - Similarly, the machine representation `A` of the matrix is only _an approximation_ of the true matrix~$mat(A)$.

  - Finally, the operation `A\b` itself leads to additional round-off errors,
    as the computer implementation of the backslash operator is based on elementary arithmetic operations.
    Understanding the magnitude of the error introduced at this step is delicate,
    and we shall not address this question.
  It follows from the first two items that, when writing `A\b`,
  we are in fact asking the computer for a solution~$vect(x) + Delta vect(x)$ to a perturbed equation
  $ (mat(A) + Delta mat(A)) (vect(x) + Delta vect(x)) = vect(b) + Delta vect(b), $
  where $Delta mat(A)$ and~$Delta vect(b)$ represent the rounding errors on the matrix and right-hand side, respectively.
  Understanding the impact of the perturbations~$Delta vect(b)$ and $Delta mat(A)$ is the goal of this section.
  We shall prove,
  in particular,
  that the magnitude of the error~$Delta vect(x)$ is related to the perturbations~$Delta vect(b)$ and~$Delta mat(A)$ through a quantity known as the _condition number_ of the matrix~$mat(A)$.
]

#indent In general, the condition number for a given problem measures the sensitivity of the solution to the input data.
In order to define this concept precisely,
we consider a general problem of the form~$F(x, d) = 0$,
with unknown $x$ and data $d$.
The linear system~@eq:linear_system can be recast in this form,
with the input data equal to $vect(b)$ or $mat(A)$ or both.
We denote the solution corresponding to perturbed input data $d + Delta d$ by $x + Delta x$.
The absolute and relative condition numbers are defined as follows.

#definition(title: [Condition number for the problem $F(x, d) = 0$])[
  The absolute and relative condition numbers with respect to perturbations of $d$ are defined as
  $
    K_("abs")(d) = lim_(epsilon -> 0) lr((sup_(norm(Delta d) <= epsilon) norm(Delta x)/norm(Delta d))),
    wide
    K(d) = lim_(epsilon -> 0) lr((sup_(norm(Delta d) <= epsilon) (norm(Delta x) \/ norm(x))/(norm(Delta d) \/ norm(d)))).
  $
  The short notation $K$ is reserved for the relative condition number,
  which is often more useful in applications.
]

#indent In the rest of this section,
we obtain an upper bound on the relative condition number for the linear system~@eq:linear_system with respect to perturbations first of $vect(b)$,
and then of $mat(A)$.
We use the notation~$norm(placeholder)$ to denote both a vector norm on $complex^n$ and the induced operator norm on matrices.

#proposition(title: [Perturbation of the right-hand side])[
  Let $vect(x) + Delta vect(x)$ denote the solution to the perturbed equation $mat(A) (vect(x) + Delta vect(x)) = vect(b) + Delta vect(b)$.
  Then it holds that
  #numbered($
    norm(Delta vect(x))/norm(vect(x)) <= norm(mat(A)) norm(mat(A)^(-1)) thin norm(Delta vect(b))/norm(vect(b)),
  $) <eq:linear_perturbation_rhs>
] <proposition:linear_perturbation_rhs>
#proof[
  It holds by definition of $Delta vect(x)$ that $mat(A) Delta vect(x) = Delta vect(b)$.
  Therefore, we have
  #numbered($
    norm(Delta vect(x))
    = norm(mat(A)^(-1) Delta vect(b))
    <= norm(mat(A)^(-1)) norm(Delta vect(b))
    = norm(mat(A) vect(x))/norm(vect(b)) norm(mat(A)^(-1)) norm(Delta vect(b))
    <= (norm(mat(A)) norm(vect(x)))/norm(vect(b)) norm(mat(A)^(-1)) norm(Delta vect(b)).
  $) <eq:linear_perturbation_rhs_to_rearrange>
  Here we employed~@eq:submultiplicative_mat_vec,
  proved in @cha:vectors_and_matrices,
  in the first and last inequalities.
  Rearranging the inequality~@eq:linear_perturbation_rhs_to_rearrange,
  we obtain~@eq:linear_perturbation_rhs.
]
@proposition:linear_perturbation_rhs implies that
the relative condition number of~@eq:linear_system with respect to perturbations of the right-hand side is bounded from above by $norm(mat(A)) norm(mat(A)^(-1))$.
@exercise:linear_sharp_inequality shows that there are values of $vect(x)$ and $Delta vect(b)$ for which the inequality~@eq:linear_perturbation_rhs is an equality,
indicating that the inequality is sharp.

Studying the impact of perturbations of the matrix~$mat(A)$ is slightly more difficult,
because this time the variation~$Delta vect(x)$ of the solution does not depend linearly on the perturbation of the data.
Before stating and proving the main result,
we show an ancillary lemma.
#lemma[
  Let $mat(B) in complex^(n times n)$ be such that $norm(mat(B)) < 1$ in some submultiplicative matrix norm~$norm(placeholder)$.
  Then $mat(I) - mat(B)$ is invertible and
  #numbered($
    norm((mat(I) - mat(B))^(-1))
    <= 1/(1 - norm(mat(B))),
  $) <eq:linear_bound_inverse_perturbation_identity>
  where $mat(I) in complex^(n times n)$ is the identity matrix.
] <lemma:linear_inverse_neumann>
#proof[
  It holds for any matrix $mat(B) in complex^(n times n)$ that
  #numbered($
    mat(I) - mat(B)^(n+1) = (mat(I) - mat(B))(mat(I) + mat(B) + cdots + mat(B)^n).
  $) <eq:first_equation_neumann>
  Since $norm(mat(B)) < 1$ in a submultiplicative matrix norm,
  both sides of the equation are convergent in the limit as $n -> infinity$.
  The left-hand side converges to the identity matrix $mat(I)$,
  and the right-hand side converges as $n -> infinity$ because ${mat(S)_0, mat(S)_1, dots}$ with
  $ mat(S)_n := mat(I) + mat(B) + cdots + mat(B)^n $
  is a Cauchy sequence in the vector space of matrices endowed with the norm for which~$norm(mat(B)) < 1$.
  Indeed, by the triangle inequality and the submultiplicative property of the norm,
  it holds that
  $
    norm(mat(S)_(n+m) - mat(S)_n)
    &= norm(mat(B)^(n+1) + cdots + mat(B)^(n+m)) \
    &<= norm(mat(B)^(n+1)) + cdots + norm(mat(B)^(n+m))
    <= norm(mat(B))^(n+1) + cdots + norm(mat(B))^(n+m) \
    &<= norm(mat(B))^(n+1)/(1 - norm(mat(B))) stretch(->)_(n -> infinity) mat(0),
  $
  where we employed the formula for a geometric series in the last inequality.
  Equating the limits of both sides of~@eq:first_equation_neumann,
  we obtain
  $ mat(I) = (mat(I) - mat(B)) sum_(i=0)^infinity mat(B)^i. $
  This implies that $(mat(I) - mat(B))$ is invertible with inverse
  given by a so-called _Neumann_ series
  $ (mat(I) - mat(B))^(-1) = sum_(i=0)^infinity mat(B)^i. $
  Applying the triangle inequality repeatedly,
  and then using the submultiplicative property of the norm,
  we obtain
  $
    forall n in nat,
    wide
    lr(||sum_(i=0)^n mat(B)^i||)
    <= sum_(i=0)^n norm(mat(B)^i)
    <= sum_(i=0)^n norm(mat(B))^i
    = 1/(1 - norm(mat(B))).
  $
  where we used the summation formula for geometric series in the last equality.
  Letting $n -> infinity$ in this equation and
  using the continuity of the norm enables to conclude the proof.
]

#proposition(title: [Perturbation of the matrix])[
  Let $vect(x) + Delta vect(x)$ denote the solution to the perturbed equation $(mat(A) + Delta mat(A)) (vect(x) + Delta vect(x)) = vect(b)$.
  If $mat(A)$ is invertible and $norm(Delta mat(A)) < norm(mat(A)^(-1))^(-1)$,
  then
  #numbered($
    norm(Delta vect(x))/norm(vect(x))
    <= norm(mat(A)) norm(mat(A)^(-1)) norm(Delta mat(A))/norm(mat(A))
    lr((1/(1 - norm(mat(A)^(-1) Delta mat(A))))).
  $) <eq:linear_perturbation_matrix>
] <proposition:linear_perturbation_matrix>
#proof[
  Left-multiplying both sides of the perturbed equation with $mat(A)^(-1)$,
  we obtain
  #numbered($
    lr((mat(I) + mat(A)^(-1) Delta mat(A))) (vect(x) + Delta vect(x)) = vect(x)
    quad <=> quad
    lr((mat(I) + mat(A)^(-1) Delta mat(A))) Delta vect(x) = - mat(A)^(-1) Delta mat(A) vect(x).
  $) <eq:linear_perturbation_matrix_initial>
  Since $norm(mat(A)^(-1) Delta mat(A)) <= norm(mat(A)^(-1)) norm(Delta mat(A)) < 1$ by assumption,
  we deduce from~@lemma:linear_inverse_neumann that the matrix on the left-hand side is invertible
  with a norm bounded as in~@eq:linear_bound_inverse_perturbation_identity.
  Consequently,
  using in addition the assumed submultiplicative property of the norm,
  we obtain that
  $
    norm(Delta vect(x))
    = norm((mat(I) + mat(A)^(-1) Delta mat(A))^(-1) mat(A)^(-1) Delta mat(A) vect(x))
    <= norm(mat(A)^(-1) Delta mat(A))/(1 - norm(mat(A)^(-1) Delta mat(A))) norm(vect(x)).
  $
  which enables to conclude the proof.
]
Using @proposition:linear_perturbation_matrix,
we deduce that the relative condition number of~@eq:linear_system with respect to perturbations of the matrix $mat(A)$ is also bounded from above by $norm(mat(A)) norm(mat(A)^(-1))$,
because the term between brackets on the right-hand side of~@eq:linear_perturbation_matrix converges to 1 as $norm(Delta mat(A)) -> 0$.

@proposition:linear_perturbation_rhs[Propositions] and @proposition:linear_perturbation_matrix[] show that
the condition number of the linear system~@eq:linear_system,
with respect to perturbations of either~$vect(b)$ or~$mat(A)$,
depends only on $mat(A)$.
This motivates the following definition of the condition number~_of a matrix_.
#definition(title: [Condition number of a matrix])[
  The condition number of a matrix $mat(A)$ associated with a vector norm $norm(placeholder)$ is defined as
  $ kappa(mat(A)) = norm(mat(A)) norm(mat(A)^(-1)). $
  The condition number for the $p$-norm,
  defined in @definition:pnorm_vector,
  is denoted by $kappa_p (mat(A))$.
]
Note that the condition number $kappa(mat(A))$ associated with an induced norm,
i.e.~a matrix norm induced by a vector norm,
is at least one.
Indeed, since the identity matrix has induced norm~1,
it holds that
$ 1 = norm(mat(I)) = norm(mat(A) mat(A)^(-1)) <= norm(mat(A)) norm(mat(A)^(-1)). $

Since the 2-norm of an invertible matrix $mat(A) in complex^(n times n)$ coincides with the spectral radius $rho(mat(A)^T mat(A))$,
the condition number $kappa_2$ corresponding to the $2$-norm is equal to
$ kappa_2 (mat(A)) = sqrt((lambda_(max)(mat(A)^T mat(A)))/(lambda_(min)(mat(A)^T mat(A)))), $
where $lambda_(max)(mat(A)^T mat(A))$ and $lambda_(min)(mat(A)^T mat(A))$ are the maximal and minimal (both real and positive) eigenvalues of the matrix $mat(A)^T mat(A)$.
#example(title: [Perturbation of the matrix])[
  Consider the following linear system
  with perturbed matrix
  $
    (mat(A) + Delta mat(A))
    pmat(x_1; x_2)
    = pmat(0; 0.01),
    wide
    mat(A)
    = pmat(1, 0; 0, 0.01),
    wide
    Delta mat(A) =
    pmat(0, 0; 0, epsilon),
  $
  where $0 < epsilon << 0.01$.
  Here the eigenvalues of $mat(A)$ are given by $lambda_1 = 1$ and $lambda_2 = 0.01$.
  The solution when $epsilon = 0$ is given by $(0, 1)^T$,
  and the solution to the perturbed equation is
  $
    pmat(x_1 + Delta x_1; x_2 + Delta x_2)
    =
    pmat(0; 1/(1 + 100 epsilon)).
  $
  Consequently, we deduce that, in the 2-norm,
  $
    norm(Delta vect(x))/norm(vect(x))
    = lr(|(100 epsilon)/(1 + 100 epsilon)|)
    approx 100 epsilon
    = 100 norm(Delta mat(A))/norm(mat(A)).
  $
  In this case,
  the relative impact of perturbations of the matrix is close to $kappa_2 (mat(A)) = 100$.
]

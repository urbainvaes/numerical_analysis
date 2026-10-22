#import "../template.typ": *

= Background material <cha:vectors_and_matrices>
#minitoc()

In this chapter,
we collect basic results that are useful for this course.

== Inner products and norms <sec:inner_product_and_norm>

We begin by recalling the definitions of the fundamental concepts of _norm_ and _inner product_.
For generality,
we consider the case of a _complex_ vector space,
i.e.~a vector space for which the scalar field is $complex$.
#definition[
  A norm on a complex vector space $cal(X)$ is a function $norm(placeholder): cal(X) -> real$ satisfying the following axioms:
  - *Positivity*:
    $forall vect(x) in cal(X) backslash {vect(0)}, quad
      norm(vect(x)) > 0.$

  - *Homogeneity*:
    $forall (c, vect(x)) in complex times cal(X), quad
      norm(c vect(x)) = abs(c) thin norm(vect(x)) .$

  - *Triangular inequality*:
    $forall (vect(x), vect(y)) in cal(X) times cal(X), quad
      norm(vect(x) + vect(y)) <= norm(vect(x)) + norm(vect(y)).$
]
For example,
the Euclidean norm on~$complex^n$ is given by
$ norm(vect(x)) = sqrt(|x_1|^2 + dots + |x_n|^2). $

#definition[
  An inner product on a _complex_ vector space $cal(X)$ is a function
  $ ip(placeholder, placeholder): cal(X) times cal(X) -> complex $
  satisfying the following axioms:
  - *Conjugate symmetry*:
    Here $overline(placeholder)$ denotes the complex conjugate.
    $
      forall (vect(x), vect(y)) in cal(X) times cal(X), wide
      ip(vect(x), vect(y)) = overline(ip(vect(y), vect(x))).
    $
  - *Linearity*:
    For all $(alpha, beta) in complex^2$ and all $(vect(x), vect(y), vect(z)) in cal(X)^3$,
    it holds that
    $
      ip(alpha vect(x) + beta vect(y), vect(z))
      = alpha ip(vect(x), vect(z)) + beta ip(vect(y), vect(z)).
    $

  - *Positive-definiteness*:
    $
      forall vect(x) in cal(X) backslash {0}, wide
      ip(vect(x), vect(x)) > 0.
    $
]
For example, the familiar Euclidean inner product on $complex^n$ is given by
$ ip(vect(x), vect(y)) := sum_(i=1)^n x_i overline(y)_i. $
A vector space with an inner product is called an _inner product space_.
Any inner product on~$cal(X)$ induces a norm via the formula
#numbered($ norm(vect(x)) = sqrt(ip(vect(x), vect(x))). $) <eq:induced_norm>
The Cauchy--Schwarz inequality enables to bound inner products using norms.
It is also useful for showing that the functional defined in~@eq:induced_norm satisfies the triangle inequality,
which is the goal of @exercise:triangle_inequality.
#proposition(title: [Cauchy--Schwarz inequality])[
  Let $cal(X)$ be an inner product space.
  Then
  #numbered($
    forall (vect(x), vect(y)) in cal(X) times cal(X), wide
    big(|ip(vect(x), vect(y))|) <= norm(vect(x)) norm(vect(y)).
  $) <eq:cauchy_shwartz>
] <proposition:cauchy_shwartz>
#proof[
  The statement is obvious if $vect(y) = vect(0)$,
  so we assume in the rest of the proof that $vect(y) != vect(0)$.
  Let us define $p colon real in.rev t |-> norm(vect(x) + t vect(y))^2$.
  Using the bilinearity of the inner product,
  we have
  $ p(t) = norm(vect(x))^2 + 2 t ip(vect(x), vect(y)) + t^2 norm(vect(y))^2. $
  This shows that $p$ is a convex second-order polynomial with a minimum at $t_* = - ip(vect(x), vect(y)) \/ norm(vect(y))^2$.
  Substituting this value in the expression of $p$,
  we obtain
  $
    p(t_*) = norm(vect(x))^2 - 2 (big(|ip(vect(x), vect(y))|)^2)/(norm(vect(y))^2) + (big(|ip(vect(x), vect(y))|)^2)/(norm(vect(y))^2)
    = norm(vect(x))^2 - (big(|ip(vect(x), vect(y))|)^2)/(norm(vect(y))^2).
  $
  Since $p(t_*) >= 0$ by definition of $p$,
  we obtain~@eq:cauchy_shwartz.
]

#indent Several norms can be defined on the same vector space $cal(X)$.
Two norms $norm(placeholder)_alpha$ and $norm(placeholder)_beta$ on $cal(X)$ are said to be equivalent if
there exist positive real numbers~$c_ell$ and $c_u$ such that
#numbered($
  forall vect(x) in cal(X),
  wide c_ell norm(vect(x))_alpha
  <= norm(vect(x))_beta
  <= c_u norm(vect(x))_alpha.
$) <eq:norm_equilavence>
As the terminology indicates, norm equivalence is an _equivalence relation_.
When working with norms on finite-dimensional vector spaces,
it is important to keep in mind the following result.
The proof is provided for information purposes only.
#proposition[
  Assume that $cal(X)$ is a finite-dimensional vector space.
  Then all the norms defined on~$cal(X)$ are pairwise equivalent.
]
#proof[
  Let $(vect(e)_1, dots, vect(e)_n)$ be a basis of $cal(X)$.
  Any $vect(x) in cal(X)$ admits a unique representation in this basis as $vect(x) = lambda_1 vect(e)_1 + cdots + lambda_n vect(e)_n$.
  We will show that any norm $norm(placeholder)$ on~$cal(X)$ is equivalent to the norm $norm(placeholder)_*$ given by
  #numbered($ norm(vect(x))_* = abs(lambda_1) + dots + abs(lambda_n). $) <eq:norm_star>
  By the triangle inequality,
  it holds that
  #eq-last($
    norm(vect(x)) <= abs(lambda_1) norm(vect(e)_1) + cdots + abs(lambda_n) norm(vect(e)_n)
    &<= Big(\() abs(lambda_1) + cdots + abs(lambda_n) Big(\)) thin
    max Big(\{) norm(vect(e)_1), dots, norm(vect(e)_n) Big(\}) \
    &= norm(vect(x))_* max Big(\{) norm(vect(e)_1), dots, norm(vect(e)_n) Big(\}).
  $) <eq:norm_equivalence_one>
  It remains to show that
  there exists a positive constant $ell$ such that
  #numbered($
    forall vect(x) in cal(X), wide
    norm(vect(x))
    >= ell Big(\() abs(lambda_1) + cdots + abs(lambda_n) Big(\)).
  $) <eq:bound2>
  To this end, we reason by contradiction.
  If this inequality were not true,
  then there would exist a sequence $(vect(x)^((i)))_(i in nat)$ such
  that~$norm(vect(x)^((i))) -> 0$ as $i -> infinity$ and $norm(vect(x)^((i)))_* = 1$ for all~$i in nat$.
  Since~$lambda_1^((i)) in [-1, 1]$ for all~$i in nat$,
  we can extract a subsequence, still denoted by~$(vect(x)^((i)))_(i in nat)$ for simplicity,
  such that the corresponding coefficient $lambda_1^((i))$ satisfies $lambda_1^((i)) -> lambda_1^* in [-1, 1]$,
  by compactness of the interval~$[-1, 1]$.
  Repeating this procedure for $lambda_2, lambda_3, dots$,
  taking a new subsequence every time,
  we obtain a subsequence $(vect(x)^((i)))_(i in nat)$ such that $lambda^((i))_j -> lambda_j^*$ in the limit as~$i -> infinity$, for all $j in {1, dots, n}$.
  Therefore, it holds that $vect(x)^((i)) -> vect(x)^* := lambda_1^* vect(e)_1 + cdots lambda_n^* vect(e)_n$ in the~$norm(placeholder)_*$ norm,
  and thus also in the~$norm(placeholder)$ norm by~@eq:norm_equivalence_one.
  Since $vect(x)^((i)) -> vect(0)$ in the latter norm by assumption,
  we deduce that $vect(x)^* = vect(0)$.
  But the vectors $vect(e)_1, dots, vect(e)_n$ are linearly independent,
  and so this implies that $lambda_1^* = dots = lambda_n^* = 0$,
  which is a contradiction because we also have that
  $
    abs(lambda_1^*) + cdots + abs(lambda_n^*)
    = lim_(i -> infinity) big(|lambda_1^((i))|) + cdots + big(|lambda_n^((i))|) = 1.
  $
  This concludes the proof of~@eq:bound2.
]

#exercise[
  Show that $norm(placeholder)_* colon cal(X) -> real$ defined in~@eq:norm_star is indeed a norm.
]

#exercise[
  Using~@proposition:cauchy_shwartz,
  show that the function $norm(placeholder)$ defined by~@eq:induced_norm satisfies the triangle inequality.
] <exercise:triangle_inequality>

== Completeness
Assume that $cal(X)$ is a vector space with a norm~$norm(placeholder)$.
Together, $(cal(X), norm(placeholder))$ form a _normed vector space_.
A sequence $(vect(x)_n)_(n >= 0)$ in~$cal(X)$ is convergent in this space
if there exists $vect(x)_* in cal(X)$ such that
$ norm(vect(x)_n - vect(x)_*) -> 0 quad " in the limit " n -> infinity". " $
In this case,
we write $lim_(n -> infinity) vect(x)_n = vect(x)_*$ or $vect(x)_n -> vect(x)_*$.
#definition(title: [Cauchy sequence])[
  A sequence $(vect(x)_n)_(n >= 0)$ in~$cal(X)$ is said to be Cauchy if
  $ lim_(m, n -> +infinity) norm(vect(x)_n - vect(x)_m) = 0. $
  The normed vector space $(cal(X), norm(placeholder))$ is called _complete_ if every Cauchy sequence is convergent.
]
Every convergent sequence is Cauchy,
but the converse is not always true.

#example[
  Consider the case where $cal(X) = C([-1, 1])$,
  the space of continuous functions from~$[-1, 1]$ to~$real$,
  endowed with the norm
  $ norm(f) = integral_(-1)^1 abs(f(x)) thin dd x. $
  The sequence of functions $(f_n)_(n >= 0)$ in~$cal(X)$ given by
  #numbered($ f_n (x) = x^(1/(2n+1)) $) <eq:example_non_complete>
  is Cauchy but not convergent.
  Indeed, assume for contradiction that there existed~$f_* in cal(X)$ such~that
  #numbered($ norm(f_n - f_*) stretch(->)_(n -> infinity) 0. $) <eq:example_non_complete_convergence>
  Then $f_*$ necessarily coincides with sign function:
  $
    upright("sgn")(x)
    := cases(
      -1 & "if " x < 0\,,
      0 & "if " x = 0\,,
      1 & "if " x > 0.
    )
  $
  However, this function is discontinuous,
  which contradicts the statement that~$f_* in cal(X)$.
]

#indent In this course,
all the vector spaces encountered are complete.
For example,
- $real^n$ with any vector norm is complete;
- $complex^(m times n)$ with any matrix norm is complete.
In order to show that a sequence is convergent in a complete normed vector space,
it is sufficient to show that the sequence is Cauchy.
This approach is used in~@lemma:linear_inverse_neumann and~@theorem:banach_fixed_point.

#exercise[
  Prove that every convergent sequence is Cauchy.
]

== Contraction mappings and the Banach fixed point theorem
Let $(cal(X), norm(placeholder))$ denote a normed vector space.
A map $phi.alt colon cal(X) -> cal(X)$ is called a contraction mapping if
there is a constant $L in (0, 1)$ such that
$
  forall (x, y) in cal(X) times cal(X), wide
  norm(phi.alt(x) - phi.alt(y)) <= L norm(x - y).
$
The importance of contraction mappings in this course stems from the following theorem.

#theorem(title: [Banach fixed point theorem])[
  Let $(cal(X), norm(placeholder))$ be a complete normed space,
  and let~$phi.alt colon cal(X) -> cal(X)$ be a contraction mapping.
  Then $phi.alt$ has a unique fixed point in~$cal(X)$.
] <theorem:banach_fixed_point>
// The paragraphs of this proof are not indented (\noindent)
#[#set par(first-line-indent: 0pt)
#proof[
  We prove first existence and then uniqueness.

  *Existence.*
  Take $x_0 in cal(X)$,
  and define the sequence $(x_k)_(k in nat)$ inductively by
  #numbered($ x_(k+1) = phi.alt(x_k). $) <eq:inductive_sequence>
  It holds that
  $
    big(||x_(k+1) - x_k||)
    = big(||phi.alt(x_k) - phi.alt(x_(k-1))||)
    <= L norm(x_k - x_(k-1))
    <= dots <= L^k norm(x_1 - x_0).
  $
  Therefore, for any $n >= m$,
  we have by the triangle inequality
  $
    norm(x_n - x_m)
    &<= norm(x_n - x_(n-1)) + cdots + norm(x_(m+1) - x_m) \
    &<= (L^(n-1) + cdots + L^m) norm(x_1 - x_0) \
    &<= L^m (1 + L + L^2 + cdots) norm(x_1 - x_0)
    = L^m/(1-L) norm(x_1 - x_0).
  $
  It follows that the sequence $(x_k)_(k in nat)$ is Cauchy in $cal(X)$,
  implying by completeness that $x_k -> x_*$ in the limit as $k -> infinity$,
  for some limit $x_* in cal(X)$.
  Being a contraction, the mapping $phi.alt$ is continuous,
  and so taking the limit~$k -> infinity$ in~@eq:inductive_sequence, we obtain that
  $
    x_* = lim_(k -> infinity) x_(k+1)
    = lim_(k -> infinity) phi.alt(x_k)
    = phi.alt lr((lim_(k -> infinity) x_k))
    = phi.alt(x_*).
  $
  In other words, $x_*$ is a fixed point of~$phi.alt$.

  #v(.2cm)
  *Uniqueness.*
  Assume that $y_* in cal(X)$ is a fixed point.
  Then,
  $
    big(||y_* - x_*||)
    = big(||phi.alt(y_*) - phi.alt(x_*)||)
    <= L norm(y_* - x_*),
  $
  which implies that $y_* = x_*$ since $L < 1$.
]]
#remark[
  // Notice that the vector space structure is not employed in the proof of @theorem:banach_fixed_point.
  The Banach fixed point theorem holds also in complete metric spaces.
]

== Vector norms <sub:vector_norms>

In the vector space $complex^n$,
the most commonly used norms are particular cases of the $p$-norm, also called Hölder norm.

#definition[
  Given $p in [1, infinity]$,
  the $p$-norm of a vector $vect(x) in complex^n$ is defined as follows:
  $
    norm(vect(x))_p :=
    cases(
      lr((sum_(i=1)^n abs(x_i)^p))^(1/p) & "if " p < infinity",",
      max Big(\{) abs(x_1)\, dots\, abs(x_n) Big(\}) & "if " p = infinity".",
    )
  $
] <definition:pnorm_vector>
The values of $p$ most commonly encountered in applications are $1$, $2$ and $infinity$.
The $1$-norm is sometimes called the _taxicab_ or _Manhattan_ norm,
and the $2$-norm is usually called the _Euclidean norm_.
The explicit expressions of these norms are
$
  norm(vect(x))_1 = sum_(i=1)^n abs(x_i),
  wide
  norm(vect(x))_2 = sqrt(sum_(i=1)^n abs(x_i)^2).
$
Notice that the infinity norm $norm(placeholder)_infinity$ may be defined as the limit of the $p$-norm as $p -> infinity$:
$
  norm(vect(x))_infinity
  := lim_(p -> infinity) norm(vect(x))_p.
$

In the rest of this chapter,
the notations $ip(placeholder, placeholder)$ and $norm(placeholder)$ without subscript always refer to the Euclidean inner product~@eq:induced_norm and induced norm,
unless specified otherwise.

== Matrix norms <sec:matrix_norms>

Given two norms $norm(placeholder)_alpha$ and $norm(placeholder)_beta$ on $complex^m$ and $complex^n$, respectively,
we define the _operator norm_ induced by $norm(placeholder)_alpha$ and $norm(placeholder)_beta$ of the matrix $mat(A)$ as
#numbered($
  norm(mat(A))_(alpha, beta) = sup big(\{) norm(mat(A) vect(x))_alpha : vect(x) in complex^n, norm(vect(x))_beta <= 1 big(\}).
$) <eq:subordinate_norm>
The term _operator norm_ is motivated by the fact that,
to any matrix $mat(A) in complex^(m times n)$, there naturally corresponds the linear operator from $complex^n$ to $complex^m$ with action~$vect(x) |-> mat(A) vect(x)$.
// Equation @eq:subordinate_norm comes from the general definition of the norm of a bounded linear operator between normed spaces.
Matrix norms of the type~@eq:subordinate_norm are also called _subordinate_ matrix norms.
An immediate corollary of the definition~@eq:subordinate_norm is that,
for all $vect(x) in complex^n$,
#numbered($
  norm(mat(A) vect(x))_alpha
  = norm(mat(A) hat(vect(x)))_alpha norm(vect(x))_beta
  <= sup Big(\{) norm(mat(A) vect(y))_alpha: norm(vect(y))_beta <= 1 Big(\}) norm(vect(x))_beta
  = norm(mat(A))_(alpha, beta) norm(vect(x))_beta,
  wide hat(vect(x)) = vect(x)/norm(vect(x))_beta.
$) <eq:submultiplicative_mat_vec>

#exercise[
  Show that equation~@eq:subordinate_norm defines a norm on $complex^(m times n)$.
]

#indent The matrix $p$-norm is defined as the operator norm~@eq:subordinate_norm in the particular case
where $norm(placeholder)_alpha$ and $norm(placeholder)_beta$ are both Hölder norms with the same value of $p$.
#definition[
  Given $p in [1, infinity]$,
  the $p$-norm of a matrix $mat(A) in complex^(m times n)$ is given by
  #numbered($
    norm(mat(A))_p := sup big(\{) norm(mat(A) vect(x))_p : vect(x) in complex^n, norm(vect(x))_p <= 1 big(\}).
  $) <eq:matrix_p_norm>
]
Not all matrix norms are induced by vector norms.
For example, the Frobenius norm,
which is widely used in applications,
is not induced by a vector norm.
It is, however, induced by an inner product on~$complex^(m times n)$.
#definition[
  The Frobenius norm of $mat(A) in complex^(m times n)$ is given by
  #numbered($
    norm(mat(A))_"F" = lr((sum_(i=1)^m sum_(j=1)^n abs(a_(i j))^2))^(1/2).
  $) <eq:frobenius_norm>
]

#indent A matrix norm $norm(placeholder)$ is said to be submultiplicative if,
for any two matrices $mat(A) in complex^(m times n)$ and~$mat(B) in complex^(n times ell)$,
it holds that
$ norm(mat(A) mat(B)) <= norm(mat(A)) norm(mat(B)). $
All subordinate matrix norms,
for example the $p$-norms, are submultiplicative,
and so is the Frobenius norm.

#exercise[
  Write down the inner product on $complex^(m times n)$ corresponding to~@eq:frobenius_norm.
]
#exercise[
  Show that the matrix $p$-norm is submultiplicative.
]

#pagebreak()
== Diagonalization and spectral theorem <sec:diagonalization>

#definition[
  A square matrix $mat(A) in complex^(n times n)$ is said to be diagonalizable if there exists an invertible matrix $mat(P) in complex^(n times n)$
  and a diagonal matrix $mat(D) in complex^(n times n)$ such that
  #numbered($ mat(A) mat(P) = mat(P) mat(D). $) <eq:eigen_decomposition>
  In this case,
  the diagonal elements of $mat(D)$ are called the eigenvalues of~$mat(A)$,
  and the columns of $mat(P)$ are called the eigenvectors of~$mat(A)$.
] <definition:diagonalizable>
Denoting by $vect(e)_i$ the $i$-th column of $mat(P)$ and by $lambda_i$ the $i$-th diagonal element of $mat(D)$,
we have by~@eq:eigen_decomposition that $mat(A) vect(e)_i = lambda_i vect(e)_i$ or,
equivalently, $(mat(A) - lambda_i mat(I)_n) vect(e)_i = vect(0)$.
Here $mat(I)_n$ is the~$complex^(n times n)$ identity matrix.
Therefore, a complex number $lambda$ is an eigenvalue of $mat(A)$ if and only if $det(mat(A) - lambda mat(I)_n) = 0$.
In other words, the eigenvalues of $mat(A)$ are the roots of $det(mat(A) - lambda mat(I)_n)$,
which is called the _characteristic polynomial_.

#subsection-star[Symmetric matrices and spectral theorem]
The transpose of a matrix $mat(A) in complex^(m times n)$ is denoted by $mat(A)^T in complex^(n times m)$
and defined as the matrix with entries~$a^T_(i j) = a_(j i)$.
The conjugate transpose of $mat(A)$, denoted by $mat(A)^*$,
is the matrix obtained by taking the transpose and taking the complex conjugate of all the entries.
A real matrix that is equal to its transpose is necessarily square and called _symmetric_,
and a complex matrix that is equal to its conjugate transpose is called _Hermitian_.
Hermitian matrices, of which real symmetric matrices are a subset,
enjoy many nice properties,
the main one being that they are diagonalizable with a matrix $mat(Q)$ that is unitary,
i.e. such that $mat(Q)^(-1) = mat(Q)^*$.
This is the content of the _spectral theorem_,
a pillar of linear algebra with important generalizations to infinite-dimensional operators.

#theorem(title: [Spectral theorem for Hermitian matrices])[
  If $mat(A) in complex^(n times n)$ is Hermitian,
  then there exists an unitary matrix $mat(Q) in complex^(n times n)$ and a diagonal matrix $mat(D) in real^(n times n)$ such that
  $ mat(A) mat(Q) = mat(Q) mat(D). $
] <theorem:spectral_theorem>
#proof(title: [Sketch of the proof])[
  The result is trivial for $n = 1$.
  Reasoning by induction,
  we assume that the result is true for Hermitian matrices in $complex^(n-1 times n-1)$
  and prove that it then also holds for~$mat(A) in complex^(n times n)$.

  #v(.3cm)
  *Step 1. Existence of a real eigenvalue*.
  By the fundamental theorem of algebra,
  there exists at least one solution $lambda_1 in complex$ to the equation $det(mat(A) - lambda mat(I)_n) = 0$,
  to which there corresponds at least one solution $vect(q)_1 in complex^n$ of norm 1 to the equation $(mat(A) - lambda_1 mat(I)_n) vect(q)_1 = vect(0)$.
  The eigenvalue $lambda_1$ is necessarily real because
  $
    lambda_1 ip(vect(q)_1, vect(q)_1)
    = ip(lambda_1 vect(q)_1, vect(q)_1)
    = ip(mat(A) vect(q)_1, vect(q)_1)
    = ip(vect(q)_1, mat(A) vect(q)_1)
    = ip(vect(q)_1, lambda_1 vect(q)_1)
    = overline(lambda)_1 ip(vect(q)_1, vect(q)_1).
  $

  #v(.3cm)
  #indent *Step 2. Using the induction hypothesis*.
  Next,
  take an orthonormal basis $(vect(e)_2, dots, vect(e)_n)$ of the orthogonal complement $Span{vect(q)_1}^perp$
  and construct the unitary matrix
  $ V = pmat(vect(q)_1, vect(e)_2, dots, vect(e)_n), $
  i.e.~the matrix with columns $vect(q)_1$, $vect(e)_2$, etc.
  A calculation gives,
  $
    mat(V)^* mat(A) mat(V) =
    pmat(
      ip(vect(q)_1, mat(A) vect(q)_1), ip(vect(q)_1, mat(A) vect(e)_2), dots, ip(vect(q)_1, mat(A) vect(e)_n);
      ip(vect(e)_2, mat(A) vect(q)_1), ip(vect(e)_2, mat(A) vect(e)_2), dots, ip(vect(e)_2, mat(A) vect(e)_n);
      ⋮, ⋮, ⋱, ⋮;
      ip(vect(e)_n, mat(A) vect(q)_1), ip(vect(e)_n, mat(A) vect(e)_2), dots, ip(vect(e)_n, mat(A) vect(e)_n);
    )
    =
    pmat(
      lambda_1, 0, dots, 0;
      0, ip(vect(e)_2, mat(A) vect(e)_2), dots, ip(vect(e)_2, mat(A) vect(e)_n);
      ⋮, ⋮, ⋱, ⋮;
      0, ip(vect(e)_n, mat(A) vect(e)_2), dots, ip(vect(e)_n, mat(A) vect(e)_n);
    ).
  $
  Let us denote the $n-1 times n-1$ lower right block of this matrix by $mat(V)_(n-1)$.
  This is a Hermitian matrix of size $n-1$ so,
  using the induction hypothesis,
  we deduce that $mat(V)_(n-1) = mat(Q)_(n-1) mat(D)_(n-1) mat(Q)_(n-1)^*$
  for appropriate matrices $mat(Q)_(n-1) in complex^(n-1 times n-1)$ and $mat(D)_(n-1) in real^(n-1 times n-1)$
  which are unitary and diagonal, respectively.

  #v(.3cm)
  *Step 3. Constructing $mat(Q)$ and $mat(D)$*.
  Define now
  $
    mat(Q) =
    mat(V)
    pmat(
      1, vect(0)^T;
      vect(0), mat(Q)_(n-1);
    ).
  $
  It is not difficult to verify that $mat(Q)$ is a unitary matrix,
  and we have
  $
    mat(Q)^* mat(A) mat(Q) =
    pmat(
      1, vect(0)^T;
      vect(0), mat(Q)_(n-1)^*;
    )
    mat(V)^* mat(A) mat(V)
    pmat(
      1, vect(0)^T;
      vect(0), mat(Q)_(n-1);
    )
    =
    pmat(
      1, vect(0)^T;
      vect(0), mat(Q)_(n-1)^*;
    )
    pmat(
      lambda_1, vect(0)^T;
      vect(0), mat(V)_(n-1);
    )
    pmat(
      1, vect(0)^T;
      vect(0), mat(Q)_(n-1);
    ).
  $
  Developing the last expression, we obtain
  $
    mat(Q)^* mat(A) mat(Q)
    =
    pmat(
      lambda_1, vect(0)^T;
      vect(0), mat(D)_(n-1);
    ),
  $
  which concludes the proof.
]

#indent We deduce, as a corollary of the spectral theorem,
that if $vect(e)_1$ and $vect(e)_2$ are eigenvectors of a Hermitian matrix associated with different eigenvalues,
then they are necessarily orthogonal for the Euclidean inner product.
Indeed, since $mat(A) = mat(A)^*$ and the eigenvalues are real, it holds that
$
  (lambda_1 - lambda_2) ip(vect(e)_1, vect(e)_2)
  &= ip(lambda_1 vect(e)_1, vect(e)_2) - ip(vect(e)_1, overline(lambda)_2 vect(e)_2) \
  &= ip(mat(A) vect(e)_1, vect(e)_2) - ip(vect(e)_1, mat(A) vect(e)_2)
  = ip(mat(A) vect(e)_1, vect(e)_2) - ip(mat(A)^* vect(e)_1, vect(e)_2) = 0.
$

The largest eigenvalue of a matrix, in modulus,
is called the _spectral radius_ and denoted by~$rho$.
The following result relates the 2-norm of a matrix to the spectral radius of $mat(A) mat(A)^*$.
#proposition[
  It holds that $norm(mat(A))_2 = sqrt(rho(mat(A)^* mat(A)))$.
]
#proof[
  Since $mat(A)^* mat(A)$ is Hermitian,
  it holds by the spectral theorem that $mat(A)^* mat(A) = mat(Q) mat(D) mat(Q)^*$ for some unitary matrix $mat(Q)$ and real diagonal matrix $mat(D)$.
  Therefore, denoting by $(mu_i)_(1 <= i <= n)$ the (positive) diagonal elements of $mat(D)$
  and introducing $vect(y) := mat(Q)^* vect(x)$,
  we have
  #eq-last($
    norm(mat(A) vect(x))
    &= sqrt(vect(x)^* mat(A)^* mat(A) vect(x))
    = sqrt(vect(x)^* mat(Q) mat(D) mat(Q)^* vect(x)) \
    &= sqrt(sum_(i=1)^n mu_i y_i^2)
    <= sqrt(rho(mat(A)^* mat(A))) sqrt(sum_(i=1)^n y_i^2)
    = sqrt(rho(mat(A)^* mat(A))) norm(vect(x)),
  $) <eq:matrices_link_norm_spectral_radius>
  where we used in the last equality the fact that $vect(y)$ has the same norm as $vect(x)$,
  because $mat(Q)$ is unitary.
  It follows from~@eq:matrices_link_norm_spectral_radius that $norm(mat(A)) <= sqrt(rho(mat(A)^* mat(A)))$,
  and the converse inequality also holds true since $norm(mat(A) vect(x)) = sqrt(rho(mat(A)^* mat(A))) norm(vect(x))$
  if $vect(x)$ is the eigenvector of $mat(A)^* mat(A)$ corresponding to an eigenvalue of modulus~$rho(mat(A)^* mat(A))$.
]

#indent To conclude this section,
we recall and prove the Courant--Fisher theorem.
#theorem(title: [Courant--Fisher Min-Max theorem])[
  The eigenvalues $lambda_1 >= lambda_2 >= cdots >= lambda_n$ of a Hermitian matrix are characterized by the relation
  #numbered($
    lambda_k = max_(cal(S), upright("dim")(cal(S)) = k) lr((min_(vect(x) in cal(S) backslash {0}) (vect(x)^* mat(A) vect(x))/(vect(x)^* vect(x)))).
  $) <eq:courant_fischer>
] <theorem:courant-fisher>
#proof[
  Let $vect(v)_1, dots, vect(v)_n$ be normalized and pairwise orthogonal eigenvectors associated with the eigenvalues $lambda_1, dots, lambda_n$,
  and let $cal(S)_k = Span {vect(v)_1, dots, vect(v)_k}$.
  Any $vect(x) in cal(S)_k$ may be expressed as a linear combination~$vect(x) = alpha_1 vect(v)_1 + cdots + alpha_k vect(v)_k$,
  and so
  $
    forall vect(x) in cal(S)_k, wide
    (vect(x)^* mat(A) vect(x))/(vect(x)^* vect(x)) =
    (sum_(i=1)^k lambda_i abs(alpha_i)^2)/(sum_(i=1)^k abs(alpha_i)^2)
    >= lambda_k.
  $
  Therefore,
  it holds that
  $ min_(vect(x) in cal(S)_k backslash {0}) (vect(x)^* mat(A) vect(x))/(vect(x)^* vect(x)) >= lambda_k, $
  which proves the $>=$ direction of~@eq:courant_fischer.
  For the $<=$ direction,
  let $cal(U)_k = Span {vect(v)_k, dots vect(v)_n}$.
  Using a well-known result from linear algebra,
  we calculate that, for any subspace $cal(S) subset complex^n$ of dimension~$k$,
  $
    dim(cal(S) inter cal(U)_k)
    &= dim(cal(S)) + dim(cal(U)_k) - dim(cal(S) + cal(U)_k) \
    &>= k + (n - k + 1) - n = 1.
  $
  Therefore, any $cal(S) subset complex^n$ of dimension $k$
  has a nonzero intersection with $cal(U)_k$.
  But since any vector in $cal(U)_k$ can be expanded as $beta_1 vect(v)_k + cdots + beta_n vect(v)_n$,
  we have
  $
    forall vect(x) in cal(U)_k, wide
    (vect(x)^* mat(A) vect(x))/(vect(x)^* vect(x)) =
    (sum_(i=k)^n lambda_i abs(alpha_i)^2)/(sum_(i=k)^n abs(alpha_i)^2)
    <= lambda_k.
  $
  This shows that
  $
    forall cal(S) subset complex^n " with " dim(cal(S)) = k,
    wide
    min_(vect(x) in cal(S) backslash {0}) (vect(x)^* mat(A) vect(x))/(vect(x)^* vect(x)) <= lambda_k,
  $
  which enables to conclude the proof.
]

#exercise[
  Prove that if $mat(A) in real^(n times n)$ is diagonalizable as in~@eq:eigen_decomposition,
  then $mat(A)^n = mat(P) mat(D)^n mat(P)^(-1)$.
]

== Similarity transformation and Jordan normal form <sec:similarity_transformation_and_jordan_normal_form>
In this section, we work with matrices in $complex^(n times n)$.
A _similarity transformation_ is a mapping of the type~$complex^(n times n) in.rev mat(A) |-> mat(P)^(-1) mat(A) mat(P) in complex^(n times n)$,
where $mat(P) in complex^(n times n)$ is a nonsingular matrix.
If two matrices are related by a similarity transformation,
then they are called _similar_,
because they may be viewed as two representations of the same linear mapping in different bases.

#definition(title: [Jordan block])[
  A Jordan block with dimension $n$ is a matrix of the form
  $
    mat(J)_n (lambda) =
    pmat(
      lambda, 1, , , ;
      , lambda, 1, , ;
      , , ⋱, ⋱, ;
      , , , lambda, 1;
      , , , , lambda;
    )
  $
  The parameter $lambda in complex$ is called the eigenvalue of the Jordan block.
]
A Jordan block is diagonalizable if and only if it is of dimension 1.
The only eigenvector of a Jordan block is $pmat(1, 0, dots, 0)^T$.
The power of a Jordan block admits an explicit expression.
#lemma[
  It holds that
  #numbered($
    mat(J)_n (lambda)^k =
    pmat(
      lambda^k, binom(k, 1) lambda^(k-1), binom(k, 2) lambda^(k-2), cdots, cdots, binom(k, n-1) lambda^(k-n+1);
      , lambda^k, binom(k, 1) lambda^(k-1), cdots, cdots, binom(k, n-2) lambda^(k-n+2);
      , , ⋱, ⋱, , ⋮;
      , , , ⋱, ⋱, ⋮;
      , , , , lambda^k, binom(k, 1) lambda^(k-1);
      , , , , , lambda^k;
    ).
  $) <eq:matrices_power_jordan_block>
  // In addition, for any $mu > abs(lambda)$, there exists $C > 0$ such that
  // $forall k in nat, norm(J_n (lambda)^k) <= C mu^k$.
] <lemma:matrices_power_jordan_block>
#proof[
  The explicit expression of the Jordan block can be obtained by decomposing the block as $mat(J)_n (lambda) = lambda mat(I) + mat(N)$
  and using the binomial formula:
  $
    (lambda mat(I) + mat(N))^k =
    sum_(i=1)^k binom(k, i) (lambda mat(I))^(k-i) mat(N)^i.
  $
  To conclude the proof,
  we use the fact that $mat(N)^i$ is a matrix with zeros everywhere except for $i$-th super-diagonal,
  which contains only ones.
  Moreover $mat(N)^i = mat(0)_(n times n)$ if $i >= n$.
]

#indent A matrix is said to be of _Jordan normal form_ if it is block-diagonal with Jordan blocks on the diagonal.
In other words, a matrix~$mat(J) in complex^(n times n)$ is of Jordan normal form if
$
  mat(J) =
  pmat(
    J_(n_1)(lambda_1), , , , ;
    , J_(n_2)(lambda_2), , , ;
    , , ⋱, , ;
    , , , J_(n_(k-1))(lambda_(k-1)), ;
    , , , , J_(n_k)(lambda_k);
  )
$
with $n_1 + cdots + n_k = n$.
Note that $lambda_1, dots, lambda_k$ are the eigenvalues of $mat(A)$.
We state without proof the following important result.

#proposition(title: [Jordan normal form])[
  Any matrix $mat(A) in complex^(n times n)$ is _similar_ to a matrix in Jordan normal form.
  In other words, there exists an invertible matrix $mat(P) in complex^(n times n)$
  and a matrix in normal Jordan form $mat(J) in complex^(n times n)$ such that
  $ mat(A) = mat(P) mat(J) mat(P)^(-1) $
] <proposition:matrices_jordan_normal_form>

== Oldenburger's theorem and Gelfand's formula
The following result establishes a necessary and sufficient condition for the convergence of $norm(mat(A)^k)$ to 0 in terms of the spectral radius of $mat(A)$,
and for any matrix norm $norm(placeholder)$.
#proposition(title: [Oldenburger])[
  Let $rho(mat(A))$ denote the spectral radius of $mat(A) in complex^(n times n)$ and $norm(placeholder)$ be a matrix norm.
  Then
  - $norm(mat(A)^k) -> 0$ in the limit as $k -> infinity$ if and only if $rho(mat(A)) < 1$.
  - $norm(mat(A)^k) -> infinity$ in the limit as $k -> infinity$ if and only if $rho(mat(A)) > 1$.
] <proposition:matrices_convergence_power_of_matrix>
#proof[
  Since all matrix norms are equivalent,
  we can assume without loss of generality that~$norm(placeholder)$ is the 2-norm.
  We prove only the equivalence $norm(mat(A)^k) -> 0 <=> rho(mat(A)) < 1$.
  The other statement can be proved similarly.
  By~@proposition:matrices_jordan_normal_form,
  there exists a nonsingular matrix $mat(P)$ such that~$mat(A) = mat(P) mat(J) mat(P)^(-1)$,
  for a matrix $mat(J) in complex^(n times n)$ which is in normal Jordan form.
  Since $rho(mat(A)) = rho(mat(J))$ and $norm(mat(A)^k) -> 0$ if and only if $norm(mat(J)^k) -> 0$,
  it is sufficient to show that $norm(mat(J)^k) -> 0 <=> rho(mat(J)) < 1$.
  The latter statement follows from the expression of the power of a Jordan block given in~@lemma:matrices_power_jordan_block.
]

#indent With this result,
we can prove Gelfand's formula,
which relates the spectral radius to the asymptotic growth of $norm(mat(A)^k)$,
and is used in @cha:solution_of_linear_systems.
#proposition(title: [Gelfand's formula])[
  Let $mat(A) in complex^(n times n)$.
  It holds for any norm that
  $ lim_(k -> infinity) norm(mat(A)^k)^(1/k) = rho(mat(A)) $
] <proposition:matrices_gelfands>
#proof[
  Let $0 < epsilon < rho(mat(A))$
  and define $mat(A)^+ = mat(A)/(rho(mat(A)) + epsilon)$ and $mat(A)^- = mat(A)/(rho(mat(A)) - epsilon)$.
  It holds by construction that~$rho(mat(A)^+) < 1$ and $rho(mat(A)^-) > 1$.
  Using~@proposition:matrices_convergence_power_of_matrix, we deduce that
  $
    lim_(k -> infinity) norm((mat(A)^+)^k) = 0, wide lim_(k -> infinity) norm((mat(A)^-)^k) = infinity.
  $
  Therefore, it holds that
  $
    limsup_(k -> infinity) lr(||(mat(A)^+)^k||)^(1/k) <= 1,
    wide
    liminf_(k -> infinity) lr(||(mat(A)^-)^k||)^(1/k) >= 1.
  $
  Substituting the expressions of $mat(A)^+$ and $mat(A)^-$,
  we deduce that
  $
    limsup_(k -> infinity) Big(||mat(A)^k||)^(1/k) <= rho(mat(A)) + epsilon,
    wide
    liminf_(k -> infinity) big(||mat(A)^k||)^(1/k) >= rho(mat(A)) - epsilon.
  $
  Since $epsilon$ was arbitrary, we obtain that
  $
    rho(mat(A))
    <= liminf_(k -> infinity) big(||mat(A)^k||)^(1/k)
    <= limsup_(k -> infinity) big(||mat(A)^k||)^(1/k)
    <= rho(mat(A)),
  $
  which implies the statement.
]

#import "../../template.typ": *
#import "@preview/cetz:0.4.2"

// \ddots and \vdots (`dots` is redefined in the template)
#let ddots = sym.dots.down
#let vdots = sym.dots.v

// algorithmic environment outside of an algorithm float: no caption, no rules
// and no line numbers. Each line is a pair (indentation level, content).
#let algorithmic(..lines) = block(width: 100%, above: leading + 9pt, below: leading + 9pt, {
  set par(first-line-indent: 0pt, leading: float-leading, spacing: float-leading)
  grid(row-gutter: float-leading, ..lines.pos().map(((level, body)) => pad(left: level * 1.5em, body)))
})

// algorithm float with an algorithmic environment without line numbers
#let algorithm-unnumbered(caption: none, ..lines) = figure(kind: "algorithm", supplement: [Algorithm],
  numbering: "1", outlined: false, caption: caption,
  grid(row-gutter: float-leading,
    ..lines.pos().map(((level, body)) => pad(left: 1em + level * 1.5em, body))))

// minted with linenos and xleftmargin=\parindent (the code lines start with a space)
#let numbered-code(body) = pad(left: parindent + 0.525em, {
  show raw.line: it => {
    box(width: 0pt, move(dx: -parindent - 0.525em, box(width: 0.4cm,
      align(right, text(font: serif, size: 6pt, fill: black, str(it.number))))))
    it
  }
  body
})

// Cycle graph with 8 nodes, as drawn by TikZ (graph [spring layout]): a regular
// octagon, node A on the right just above the horizontal, then clockwise.
// `labels` gives the label of each node A, ..., H (none for an empty node), and
// `red` the indices of the nodes filled with red!60 instead of blue!40.
#let octagon-graph(labels, red: (), factor: 1) = box(cetz.canvas({
  import cetz.draw: *
  let side = 0.93 * factor
  let R = side / (2 * calc.sin(22.5deg))
  let pts = range(8).map(k => (R * calc.cos(22.5deg - k * 45deg), R * calc.sin(22.5deg - k * 45deg)))
  line(..pts, close: true, stroke: 0.4pt)
  for (k, p) in pts.enumerate() {
    let l = labels.at(k)
    let fill = if k in red { rgb(255, 102, 102) } else { rgb(153, 153, 255) }
    circle(p, radius: if l == none { 0.13 } else { 0.29 }, fill: fill, stroke: none)
    if l != none { content(p, text(size: footnotesize, [#l])) }
  }
}))

// \resizebox{.8\textwidth}{!}{...} inside a minipage
#let resized(body) = layout(size => {
  let natural = measure(body).width
  scale(0.8 * size.width / natural * 100%, reflow: true, body)
})

#let red1 = text(fill: red)[$1$]

== Direct solution method <sec:direct_solution_method>
In this section,
we present the _direct method_ for solving linear systems of the form~@eq:linear_system
with a general invertible matrix~$mat(A) in complex^(n times n)$.
The direct method can be decomposed into three steps:
- First calculate the so-called $mat(L) mat(U)$ decomposition of $mat(A)$,
  i.e.~find an upper triangular matrix~$mat(U)$ and a _unit_ lower triangular matrix~$mat(L)$ such that
  $mat(A) = mat(L) mat(U).$
  A unit lower triangular matrix is a lower triangular matrix with only ones on the diagonal.

- Then solve
  $mat(L) vect(y) = vect(b)$
  using a method called _forward substitution_.

- Finally, solve
  $mat(U) vect(x) = vect(y)$
  using a method called _backward substitution_.
By construction, the solution $vect(x)$ thus obtained is a solution to~@eq:linear_system.
Indeed, we have that
$ mat(A) vect(x) = mat(L) mat(U) vect(x) = mat(L) vect(y) = vect(b). $

=== LU decomposition <sub:lu_decomposition>

In this section,
we first discuss the existence and uniqueness of the $mat(L) mat(U)$ factorization of a matrix.
We then describe a numerical algorithm for calculating the factors $mat(L)$ and $mat(U)$,
based on _Gaussian elimination_.

==== Existence and uniqueness of the decomposition
We present a necessary and sufficient condition for the existence of a unique $mat(L) mat(U)$ decomposition of a matrix.
To this end, we define the principal submatrix of order $i$ of a matrix $mat(A) in complex^(n times n)$
as the matrix $mat(A)_i = mat(A)[1:i, 1:i]$, in Julia notation.

#proposition[
  The $mat(L) mat(U)$ factorization of a matrix $mat(A) in complex^(n times n)$,
  with $mat(L)$ unit lower triangular and $mat(U)$ upper triangular,
  exists and is unique if and only if
  the principal submatrices of $mat(A)$ of all orders are nonsingular.
] <proposition:linear_existence_lu>
#proof[
  We prove only the "if" direction; see~@MR2265914[Theorem 3.4] for the "only if" implication.

  The statement is clear if $n = 1$.
  Reasoning by induction,
  we assume that the result is proved up to $n - 1$.
  Since the matrix $mat(A)_(n-1)$ and all its principal submatrices are nonsingular by assumption,
  it holds that $mat(A)_(n-1) = mat(L)_(n-1) mat(U)_(n-1)$
  for a unit lower triangular matrix $mat(L)_(n-1)$ and an upper triangular matrix $mat(U)_(n-1)$.
  These two matrices are nonsingular,
  for if either of them were singular then the product $mat(A)_(n-1) = mat(L)_(n-1) mat(U)_(n-1)$ would be singular as well.
  Let us decompose~$mat(A)$ as follows:
  $
    mat(A) =
    pmat(mat(A)_(n-1), vect(c); vect(d)^T, a_(n n)).
  $
  Let $vect(ell)$ and $vect(u)$ denote the solutions to $mat(L)_(n-1) vect(u) = vect(c)$ and $mat(U)_(n-1)^T vect(ell) = vect(d)$.
  These solutions exist and are unique,
  because the matrices $mat(L)_(n-1)$ and $mat(U)_(n-1)$ are nonsingular.
  Letting $u_(n n) = a_(n n) - (vect(ell)^T vect(u))^(-1)$,
  we check that $mat(A)$ factorizes as
  $
    pmat(mat(A)_(n-1), vect(c); vect(d)^T, a_(n n))
    =
    pmat(mat(L)_(n-1), vect(0)_(n-1); vect(ell)^T, 1)
    pmat(mat(U)_(n-1), vect(u); vect(0)_(n-1)^T, u_(n n)).
  $
  This completes the proof of the existence of the decomposition.
  The uniqueness of the factors follows from the uniqueness of $vect(ell)$, $vect(u)$ and $u_(n n)$.
]

#indent @proposition:linear_existence_lu raises the following question:
are there classes of matrices whose principal matrices are all nonsingular?
The answer is positive, and we mention,
as an important example, the class of positive definite matrices.
Proving this is the aim of @exercise:linear_positive_definite_matrix_nonsingular_principal_components.

#subsection-star[Gaussian elimination algorithm for computing $mat(L)$ and $mat(U)$] <sub:gaussian_elimination_algorithm_for_computing_mat_l_and_mat_u_>

So far we have presented a condition under which the $mat(L) mat(U)$ decomposition of a matrix exists and is unique,
but not a practical method for calculating the matrices $mat(L)$ and $mat(U)$.
We describe in this section an algorithm,
known as _Gaussian elimination_,
for calculating the $mat(L) mat(U)$ decomposition of a matrix.
We begin by introducing the concept of _Gaussian transformation_.
#definition[
  A Gaussian transformation is a matrix of the form $mat(M)_k = mat(I) - vect(c)^((k)) vect(e)_k^T$,
  where~$vect(e)_k$ is the column vector with entry at index $k$ equal to 1 and all the other entries equal to zero,
  and $vect(c)^((k))$ is a column vector of the following form:
  $
    vect(c)^((k)) =
    pmat(0, 0, dots, 0, c^((k))_(k+1), c^((k))_(k+2), dots, c^((k))_n)^T.
  $
]

#indent The action of a Gaussian transformation $mat(M)_k$ left-multiplying a matrix $mat(A) in complex^(n times n)$ is
to replace the rows from index $k + 1$ to index $n$ by a linear combination involving themselves and the $k$-th row.
To see this, let us denote by $(vect(r)^((i)))_(1 <= i <= n)$ the rows of a matrix $mat(T) in complex^(n times n)$.
Then, we have
$
  mat(M)_k mat(T)
  = big((mat(I) - vect(c)^((k)) vect(e)_k^T)) mat(T)
  =
  pmat(
    1, , , , , , ;
    , 1, , , , , ;
    , , ddots, , , , ;
    , , , 1, , , ;
    , , , - c^((k))_(k+1), 1, , ;
    , , , vdots, , ddots, ;
    , , , - c^((k))_n, , , 1;
  )
  pmat(
    vect(r)^((1));
    vect(r)^((2));
    vdots;
    vect(r)^((k));
    vect(r)^((k+1));
    vdots;
    vect(r)^((n))
  )
  =
  pmat(
    vect(r)^((1));
    vect(r)^((2));
    vdots;
    vect(r)^((k));
    vect(r)^((k+1)) - c^((k))_(k+1) vect(r)^((k));
    vdots;
    vect(r)^((n)) - c^((k))_(n) vect(r)^((k))
  )
$
We show in @exercise:inverse_gaussian_transformation that
the inverse of a Gaussian transformation matrix is given by
#numbered($
  (mat(I) - vect(c)^((k)) vect(e)_k^T)^(-1) = mat(I) + vect(c)^((k)) vect(e)_k^T.
$) <eq:inverse_gaussian_transformation>
The idea of the Gaussian elimination algorithm is to successively left-multiply $mat(A)$
with Gaussian transformation matrices $mat(M)_1$, then $mat(M)_2$, etc.~appropriately chosen in such a way that the matrix~$mat(A)^((k))$,
obtained after $k$ iterations,
is upper triangular up to column $k$.
That is to say, the Gaussian transformations are constructed so that
all the entries in columns~1 to~$k$ under the diagonal of the matrix $mat(A)^((k))$ are equal to zero.
The resulting matrix $mat(A)^((n-1))$ after $n-1$ iterations is then upper triangular
and satisfies
$ mat(A)^((n-1)) = mat(M)_(n-1) dots mat(M)_1 mat(A). $
Rearranging this equation,
we deduce that
$ mat(A) = (mat(M)_1^(-1) dots mat(M)_(n-1)^(-1)) mat(A)^((n-1)). $
The first factor is lower triangular by~@eq:inverse_gaussian_transformation and @exercise:linear_product_of_lower_triangular.
The product in the definition of the matrix $mat(L)$ admits a simple explicit expression.
#lemma[
  It holds that
  $
    mat(M)_1^(-1) cdots mat(M)_(n-1)^(-1)
    = (mat(I) + vect(c)^((1)) vect(e)_1^T) cdots (mat(I) + vect(c)^((n-1)) vect(e)_(n-1)^T)
    = mat(I) + sum_(i=1)^(n-1) vect(c)^((i)) vect(e)_i^T.
  $
] <lemma:linear_inverse_product_gaussian_transformations>
#proof[
  Notice that, for $i < j$,
  $
    vect(c)^((i)) vect(e)_i^T vect(c)^((j)) vect(e)_j^T
    = vect(c)^((i)) (vect(e)_i^T vect(c)^((j))) vect(e)_j^T
    = vect(c)^((i)) 0 vect(e)_j^T = 0.
  $
  The statement then follows easily by expanding the product.
]
A corollary of~@lemma:linear_inverse_product_gaussian_transformations is that all the diagonal entries of the lower triangular matrix $mat(L)$ are equal to 1;
the matrix $mat(L)$ is _unit lower triangular_.
The full expression of the matrix $mat(L)$ given the Gaussian transformations is
#numbered($
  mat(L)
  = mat(I) +
  pmat(vect(c)^((1)), dots, vect(c)^((n-1)), vect(0)_n)
  =
  pmat(
    1, , , , , ;
    c^((1))_2, 1, , , , ;
    c^((1))_3, c^((2))_3, 1, , , ;
    c^((1))_4, c^((2))_4, c^((3))_4, 1, , ;
    vdots, vdots, vdots, , ddots, ;
    c^((1))_n, c^((2))_n, c^((3))_n, dots, c^((n-1))_n, 1;
  )
$) <eq:linear_matrix_L>
Therefore, the Gaussian elimination algorithms, if all the steps are well-defined,
correctly gives the $mat(L) mat(U)$ factorization of the matrix $mat(A)$.
Of course, the success of the strategy outlined above for the calculation of the $mat(L) mat(U)$ factorization hinges on
the existence of an appropriate Gaussian transformation at each iteration.
It is not difficult to show that,
if the $(k+1)$-th diagonal entry of the matrix $mat(A)^((k))$ is nonzero for all~$k in {1, dots, n-2}$,
then the Gaussian transformation matrices exist and are uniquely defined.
#lemma[
  Assume that $mat(A)^((k))$ is upper triangular up to column $k$ included,
  with $k <= n-2$.
  If~$a^((k))_(k+1,k+1) > 0$,
  then there is a unique Gaussian transformation matrix $mat(M)_(k+1)$ such that~$mat(M)_(k+1) mat(A)^((k))$ is upper triangular up to column $k + 1$.
  This transformation matrix is given by~$mat(I) - vect(c)^((k+1)) vect(e)_(k+1)^T$, where
  $
    vect(c)^((k+1)) =
    pmat(0, 0, dots, 0, a^((k))_(k+2,k+1)/a^((k))_(k+1,k+1), a^((k))_(k+3,k+1)/a^((k))_(k+1,k+1), dots, a^((k))_(n,k+1)/a^((k))_(k+1,k+1))^T.
  $
] <lemma:linear_expression_gaussian_transformations>
#proof[
  We perform the multiplication explicitly.
  Denoting denote by $(vect(r)^((i)))_(1 <= i <= n)$ the rows of $mat(A)^((k))$,
  we have
  $
    mat(M)_(k + 1) mat(A)^((k))
    =
    pmat(
      1, , , , , , ;
      , 1, , , , , ;
      , , ddots, , , , ;
      , , , 1, , , ;
      , , , - c^((k+1))_(k+2), 1, , ;
      , , , vdots, , ddots, ;
      , , , - c^((k+1))_n, , , 1;
    )
    pmat(
      vect(r)^((1));
      vect(r)^((2));
      vdots;
      vect(r)^((k+1));
      vect(r)^((k+2));
      vdots;
      vect(r)^((n))
    )
    =
    pmat(
      vect(r)^((1));
      vect(r)^((2));
      vdots;
      vect(r)^((k+1));
      vect(r)^((k+2)) - c^((k+1))_(k+2) vect(r)^((k+1));
      vdots;
      vect(r)^((n)) - c^((k+1))_(n) vect(r)^((k+1))
    ).
  $
  We need to show that the matrix on the right-hand side is upper triangular up to column $k+1$ included.
  This is clear by definition of $vect(c)^((k+1))$ and from the fact that $mat(A)^((k))$ is upper triangular up to column $k$ by assumption.
]

#indent The diagonal elements $a^((k))_(k+1,k+1)$, where $k in {0, dots, n-2 }$, are called the pivots.
We now prove that,
if an invertible matrix $mat(A)$ admits an $mat(L) mat(U)$ factorization,
then the pivots are necessarily nonzero and the Gaussian elimination algorithm is successful.
#proposition(title: [Gaussian elimination works~#moreinfo])[
  If $mat(A)$ is invertible and admits an $mat(L) mat(U)$ factorization,
  then the Gaussian elimination algorithm is well-defined and successfully terminates.
] <proposition:linear_uniqueness_lu>
#proof[
  We denote by $vect(c)^((1)), dots, vect(c)^((n-1)),$ the columns of the matrix $mat(L) - mat(I)$.
  Then the matrices given by $mat(M)_k = mat(I) - vect(c)^((k)) vect(e)_(k)^T$,
  for $k in {1, dots, n-1 }$,
  are Gaussian transformations and it holds that
  $ mat(L) = mat(M)_1^(-1) cdots mat(M)_(n-1)^(-1) $
  in view of @lemma:linear_inverse_product_gaussian_transformations.
  Since $mat(A) = mat(L) mat(U)$ by assumption,
  the result of the product
  $ mat(M)_(n-1) cdots mat(M)_1 mat(A) = mat(U) $
  is upper triangular.
  Let us use the notation $mat(A)^((k)) = mat(M)_(k) cdots mat(M)_1 mat(A)$.
  Our goal is to show that the matrices~$mat(M)_1, dots, mat(M)_(n-1)$ are the same as those obtained by the Gaussian elimination algorithm.

  Of all the Gaussian transformations $mat(M)_1, dots, mat(M)_(n-1)$,
  only $mat(M)_1$ acts on the second row of the matrix it multiplies.
  Therefore, the entry $(2, 1)$ of~$mat(U)$ coincides with the entry $(2, 1)$ of $mat(A)^((1))$,
  which implies that $a^((1))_(2,1) = 0$.
  Then notice that $a^((k))_(3,1) = a^((1))_(3, 1)$ for all $k >= 1$,
  because the entry $(3, 1)$ of the matrix $mat(M)_2 mat(A)^((1))$ is given by $a^((1))_(3,1) - c^((2))_3 a^((1))_(2,1) = a^((1))_(3,1)$,
  and all the other transformation matrices $mat(M)_3, dots, mat(M)_(n-1)$ leave the third row invariant.
  Consequently, it holds that $a^((1))_(3,1) = u_(3,1) = 0$.
  Continuing in this manner,
  we deduce that $mat(A)^((1))$ is upper triangular in the first column and that,
  since $mat(A)$ is invertible by assumption,
  the first pivot~$a_(11)$ is nonzero.
  Since this pivot is nonzero,
  the matrix $mat(M)_1$ is uniquely defined by~@lemma:linear_expression_gaussian_transformations.

  The reasoning can then be repeated with other columns,
  in order to deduce that $mat(A)^((k))$ is upper triangular up to column $k$ and that all the pivots $a^((k-1))_(k k)$ are nonzero.
  Therefore, all the Gaussian transformation matrices are uniquely defined given~@lemma:linear_expression_gaussian_transformations.
]

#subsection-star[Computer implementation] <sub:computer_implementation>

The Gaussian elimination procedure is summarized as follows.

#algorithmic(
  (0, [$mat(A)^((0)) <- mat(A), mat(L) <- mat(I)$]),
  (0, [*for* $i in {1, dots, n-1}$ *do*]),
  (1, [Construct $mat(M)_(i)$ as in~@lemma:linear_expression_gaussian_transformations.]),
  (1, [$mat(A)^((i)) <- mat(M)_(i) mat(A)^((i-1)), mat(L) <- mat(L) mat(M)_i^(-1)$]),
  (0, [*end for*]),
  (0, [$mat(U) <- mat(A)^((n-1))$.]),
)

#indent In practice, it is not necessary to explicitly create the Gaussian transformation matrices,
or to perform full matrix multiplications.
A more realistic version of the algorithm in Julia is given below.
The code exploits the relation~@eq:linear_matrix_L between $mat(L)$ and the parameters of the Gaussian transformations.
#numbered-code[
```julia
# A is an invertible matrix of size n x n
L = [i == j ? 1.0 : 0.0 for i in 1:n, j in 1:n]
U = copy(A)
for i in 1:n-1
    for r in i+1:n
        U[i, i] == 0 && error("Pivotal entry is zero!")
        ratio = U[r, i] / U[i, i]
        L[r, i] = ratio
        U[r, i:end] -= U[i, i:end] * ratio
    end
end
# L is unit lower triangular and U is upper triangular
```
]

==== Computational cost <ssub:computational_cost>
The computational cost of the algorithm,
measured as the number of floating point operations (flops) required,
is dominated by the Gaussian transformations,
in line 9 in the above code.
All the other operations amount to a computational cost scaling as $bigo(n^2)$,
which is negligible compared to the cost of the $mat(L) mat(U)$ factorization when $n$ is large.
This factorization requires
$
  overbrace(2 times, #[`-` and `*`]) underbrace(sum_(i=1)^(n-1), #[`for i in 1:n-1`]) overbrace((n - i), #[`for r in i+1:n`]) underbrace((n - i + 1), "indices" #[`[i:end]`]) quad "flops"
  = 2/3 n^3 + bigo(n^2) quad "flops".
$

=== Backward and forward substitution <sub:backward_and_forward_substitution>
Once the $mat(L) mat(U)$ factorization has been completed,
the solution to the linear system can be obtained by first using forward, and then backward substitution,
which are just bespoke methods for solving linear systems with lower and upper triangular matrices, respectively.
Let us consider the case of a lower triangular system:
$ mat(L) vect(y) = vect(b) $
Notice that the unknown $y_1$ may be obtained from the first equation of the system.
Then, since $y_1$ is known, the value of $y_2$ can be obtained from the second equation, etc.
A simple implementation of this algorithm is as follows:
```julia
# L is unit lower triangular
y = copy(b)
for i in 2:n
    for j in 1:i-1
        y[i] -= L[i, j] * y[j]
    end
end
```

=== Gaussian elimination with pivoting~#moreinfo <sub:pivoting>
The Gaussian elimination algorithm that
we presented in @sub:lu_decomposition relies on the existence of an $mat(L) mat(U)$ factorization.
In practice,
this assumption may not be satisfied,
and in this case a modified algorithm,
called Gaussian elimination _with pivoting_,
is required.

In fact, pivoting is useful even if the usual $mat(L) mat(U)$ decomposition of $mat(A)$ exists,
as it enables to reduce the condition number of the matrices matrices~$mat(L)$ and~$mat(U)$.
There are two types of pivoting:
partial pivoting, where only the rows are rearranged through a permutation at each iteration,
and complete pivoting, where both the rows and the columns are rearranged at each iteration.

Showing rigorously why pivoting is useful is beyond the scope of this course.
In this section, we only present the partial pivoting method.
Its influence on the condition number of the factors $mat(L)$ and $mat(U)$ is studied empirically in @exercise:linear_lu_with_partial_pivoting.
It is useful at this point to introduce the concept of a row permutation matrix.

==== Row permutation matrix
#definition[
  Let $sigma : {1, dots, n} -> {1, dots, n}$ be a permutation,
  i.e.~a bijection on the set~${1, dots, n}$.
  The row permutation matrix associated with $sigma$ is the matrix with entries
  $
    p_(i j) =
    cases(
      1 & "if" i = sigma(j)",",
      0 & "otherwise."
    )
  $
] <definition:row_permutation_matrix>
When a row permutation $mat(P)$ left-multiplies a matrix $mat(B) in complex^(n times n)$,
row~$i$ of matrix $mat(B)$ is moved to row index $sigma(i)$ in the resulting matrix,
for all $i in {1, dots, n}$.
A permutation matrix has a single entry equal to~1 per row and per column,
and its inverse coincides with its transpose:~$mat(P)^(-1) = mat(P)^T$.

==== Partial pivoting
Gaussian elimination with partial pivoting applies for any invertible matrix $mat(A)$,
and it outputs 3 matrices: a row permutation $mat(P)$, a unit triangular matrix $mat(L)$,
and an upper triangular matrix $mat(U)$.
These are related by the relation
$ mat(P) mat(A) = mat(L) mat(U). $
This is sometimes called a $mat(P) mat(L) mat(U)$ decomposition of the matrix $mat(A)$.
It is not unique in general but, unlike the usual $mat(L) mat(U)$ decomposition,
it always exists provided that $mat(A)$ is invertible.
We take this for granted in this course.

The idea of partial pivoting is to rearrange the rows at each iteration of the Gaussian elimination procedure in such a manner that
the pivotal entry is as large as possible in absolute value.
One step of the procedure reads
#numbered($
  mat(A)^((k+1)) = mat(M)_(k+1) mat(P)_(k+1) mat(A)^((k)).
$) <eq:step_partial_pivoting>
Here $mat(P)_(k+1)$ is a simple row permutation matrix which,
when acting on $mat(A)^((k))$,
interchanges row~$k+1$ and row $ell$,
for some index $ell >= k+1$.
The row index $ell$ is selected in such a way that the absolute value of the pivotal entry,
in position $(k+1, k+1)$ of the product~$mat(P)_(k+1) mat(A)^((k))$, is maximum.
The matrix $mat(M)_(k+1)$ is then the unique Gaussian transformation matrix ensuring that~$mat(A)^((k+1))$ is upper triangular up to column~$k+1$,
obtained as in~@lemma:linear_expression_gaussian_transformations.
The resulting matrix $mat(A)^((n-1))$ after~$n-1$ steps of the form~@eq:step_partial_pivoting is upper triangular and satisfies
$
  mat(A)^((n-1)) = mat(M)_(n-1) mat(P)_(n-1) cdots mat(M)_(1) mat(P)_(1) mat(A)
  quad <=> quad
  mat(A) = (mat(M)_(n-1) mat(P)_(n-1) cdots mat(M)_(1) mat(P)_(1))^(-1) mat(A)^((n-1)).
$
The first factor in the decomposition of $mat(A)$ is not necessarily lower triangular.
However, using the notation $mat(M) = mat(M)_(n-1) mat(P)_(n-1) cdots mat(M)_(1) mat(P)_(1)$ and $mat(P) = mat(P)_(n-1) cdots mat(P)_(1)$,
we have
#numbered($
  mat(P) mat(A) = mat(P) mat(M)^(-1) mat(U) = (mat(P) mat(M)^(-1)) mat(U) =: mat(L) mat(U).
$) <eq:plu_decomposition>
@lemma:linear_matrix_l_pivoting below shows that,
as the notation $mat(L)$ suggests,
the matrix $mat(L) = (mat(P) mat(M)^(-1))$ on the right-hand side is indeed lower triangular.
Before stating and proving the lemma,
we note that $mat(P)$ is a row permutation matrix,
and so the solution to the linear system $mat(A) vect(x) = mat(b)$ can be obtained by solving $mat(L) mat(U) vect(x) = mat(P)^T vect(b)$
by forward and backward substitution.
Since $mat(P)$ is a very sparse matrix,
the right-hand side~$mat(P)^T vect(b)$ can be calculated very efficiently.
#lemma[
  The matrix $mat(L) = mat(P) mat(M)^(-1)$ is unit lower triangular
  with all entries bounded in absolute value from above by 1.
  It admits the expression
  $
    mat(L)
    = mat(I)
    + (mat(P)_(n-1) cdots mat(P)_(2) vect(c)^((1))) vect(e)_(1)^T
    + (mat(P)_(n-1) cdots mat(P)_(3) vect(c)^((2))) vect(e)_(2)^T
    + cdots
    + (mat(P)_(n-1) vect(c)^((n-2))) vect(e)_(n-2)^T
    + vect(c)^((n-1)) vect(e)_(n-1)^T.
  $
] <lemma:linear_matrix_l_pivoting>
#proof[
  Let $mat(M)^((k)) = mat(M)_(k) mat(P)_(k) cdots mat(M)_(1) mat(P)_(1)$ and $mat(P)^((k)) = mat(P)_(k) cdots mat(P)_(1)$.
  It is sufficient to show that
  #numbered($
    mat(P)^((k)) big((mat(M)^((k))))^(-1)
    = mat(I)
    + (mat(P)_(k) cdots mat(P)_(2) vect(c)^((1))) vect(e)_(1)^T
    + (mat(P)_(k) cdots mat(P)_(3) vect(c)^((2))) vect(e)_(2)^T
    + cdots
    + (mat(P)_(k) vect(c)^((k-1))) vect(e)_(k-1)^T
    + vect(c)^((k)) vect(e)_(k)^T
    // TeX moves an equation to the left when its number does not fit beside it
    #h(1cm)
  $) <eq:linear_lower_triangular_l>
  for all $k in {1, dots, n-1 }$.
  The statement is clear for $k = 1$,
  and we assume by induction that it is true up to~$k-1$.
  Then notice that
  $
    mat(P)^((k)) big((mat(M)^((k))))^(-1)
    &= mat(P)_(k) lr(( mat(P)^((k-1)) big((mat(M)^((k-1))))^(-1) )) mat(P)_k^(-1) mat(M)_k^(-1) \
    &= mat(P)_(k) Big(\() mat(I)
    + (mat(P)_(k-1) cdots mat(P)_(2) vect(c)^((1))) vect(e)_(1)^T
    + cdots
    + (mat(P)_(k-1) vect(c)^((k-2))) vect(e)_(k-2)^T
    + vect(c)^((k-1)) vect(e)_(k-1)^T Big(\)) mat(P)_k^(-1) mat(M)_k^(-1) \
    &= Big(\() mat(I)
    + (mat(P)_k mat(P)_(k-1) cdots mat(P)_(2) vect(c)^((1))) vect(e)_(1)^T
    + cdots
    + (mat(P)_k mat(P)_(k-1) vect(c)^((k-2))) vect(e)_(k-2)^T
    + (mat(P)_k vect(c)^((k-1))) vect(e)_(k-1)^T Big(\)) mat(M)_k^(-1).
  $
  In the last equality,
  we used that $vect(e)_(i)^T mat(P)_k^(-1) = (mat(P)_k vect(e)_i)^T = vect(e)_i^T$ for all $i in {1, dots, k-1}$,
  because the row permutation $mat(P)_k$ does not affect rows~$1$ to~$k-1$.
  Using the expression $mat(M)_k^(-1) = mat(I) + vect(c)^((k)) vect(e)_k^T$,
  expanding the product and noting that $vect(e)_j^T vect(c)^((k)) = 0$ if $j <= k$,
  we obtain~@eq:linear_lower_triangular_l.
  The statement that the entries are bounded in absolute value from above by 1 follows from the choice of the pivot at each iteration.
]
The expression of $mat(L)$ in @lemma:linear_matrix_l_pivoting suggests the iterative procedure given in~@algo:LU_decomposition_with_partial_pivoting for performing the $mat(L) mat(U)$ factorization with partial pivoting.
A Julia implementation of this algorithm is presented in~@julia:lu_factorization_with_pivoting.
#algorithm-unnumbered(
  caption: [$mat(L) mat(U)$ decomposition with partial pivoting],
  (0, [Assign $mat(A)^((0)) <- mat(A)$ and $mat(P) <- mat(I)$]),
  (0, [*for* $i in {1, dots, n-1}$ *do*]),
  (1, [Find the row index $k >= i$ such that $mat(A)^((i-1))_(k,i)$ is maximum in absolute value.]),
  (1, [Interchange the rows $i$ and $k$ of matrices $mat(A)^((i-1))$ and $mat(P)$, and of vectors $vect(c)^((1)), dots, vect(c)^((i-1))$.]),
  (1, [Construct $mat(M)_(i)$ with corresponding column vector $vect(c)^((i))$ as in~@lemma:linear_expression_gaussian_transformations.]),
  (1, [Assign $mat(A)^((i)) <- mat(M)_(i) mat(A)^((i-1))$]),
  (0, [*end for*]),
  (0, [Assign $mat(U) <- mat(A)^((n-1))$.]),
  (0, [Assign $mat(L) <- mat(I) + display(pmat(vect(c)^((1)), dots, vect(c)^((n-1)), vect(0)_n))$.]),
) <algo:LU_decomposition_with_partial_pivoting>

#listing(caption: [$mat(L) mat(U)$ factorization with partial pivoting.])[#block(width: 100%, inset: (left: 13pt))[
```julia
# Auxiliary function
function swap_rows!(i, j, matrices...)
    for M in matrices
        M_row_i = M[i, :]
        M[i, :] = M[j, :]
        M[j, :] = M_row_i
    end
end

n = size(A)[1]
L, U = zeros(n, 0), copy(A)
P = [i == j ? 1.0 : 0.0 for i in 1:n, j in 1:n]
for i in 1:n-1
    # Pivoting
    index_row_pivot = i - 1 + argmax(abs.(U[i:end, i]))
    swap_rows!(i, index_row_pivot, U, L, P)

    # Usual Gaussian transformation
    c = [zeros(i-1); 1.0; zeros(n-i)]
    for r in i+1:n
        ratio = U[r, i] / U[i, i]
        c[r] = ratio
        U[r, i:end] -= U[i, i:end] * ratio
    end
    L = [L c]
end
L = [L [zeros(n-1); 1.0]]
# It holds that P*A = L*U
```
]] <julia:lu_factorization_with_pivoting>

#remark[
  It is possible to show that,
  if the matrix $mat(A)$ is column diagonally dominant in the sense that
  $
    forall j in {1, dots, n}, wide
    abs(a_(j j)) >= sum_(i=1, i != j)^(n) abs(a_(i j)),
  $
  then pivoting does not have an effect:
  at each iteration,
  the best pivot is already on the diagonal.
]

=== Direct method for Hermitian positive definite matrices <sub:direct_method_for_symmetric_positive_definite_matrices>
The $mat(L) mat(U)$ factorization with partial pivoting applies to any matrix $mat(A) in complex^(n times n)$ that is invertible.
If~$mat(A)$ is Hermitian positive definite,
however, it is possible to compute a factorization into lower and upper triangular matrices at half the computational cost,
using the so-called _Cholesky decomposition_.
#lemma(title: [Cholesky decomposition])[
  If~$mat(A)$ is Hermitian positive definite,
  then there exists a lower-triangular matrix $mat(C) in complex^(n times n)$ such that
  #numbered($
    mat(A) = mat(C) mat(C)^*.
  $) <eq:linear_cholesky_decomposition>
  Equation~@eq:linear_cholesky_decomposition is called the Cholesky factorization of $mat(A)$.
  The matrix $mat(C)$ is unique if we require that all its diagonal entries are positive.
]
#proof[
  Since $mat(A)$ is positive definite,
  its $mat(L) mat(U)$ decomposition exists and is unique by~@proposition:linear_existence_lu and~@proposition:linear_uniqueness_lu.
  Let $mat(D)$ denote the diagonal matrix with the same diagonal as that of $mat(U)$.
  Then
  $ mat(A) = mat(L) mat(D) (mat(D)^(-1) mat(U)). $
  Note that the matrix $mat(D)^(-1) mat(U)$ is unit upper triangular.
  Since $mat(A)$ is Hermitian,
  we have
  $ mat(A) = mat(A)^* = (mat(D)^(-1) mat(U))^* (mat(L) mat(D))^*. $
  The first and second factors on the right-hand side are respectively unit lower triangular and upper triangular,
  and so we deduce, by uniqueness of the $mat(L) mat(U)$ decomposition,
  that $mat(L) = (mat(D)^(-1) mat(U))^*$ and $mat(U) = (mat(L) mat(D))^*$.
  But then
  $ mat(A) = mat(L) mat(U) = mat(L) mat(D) mat(L)^* = (mat(L) sqrt(mat(D))) (sqrt(mat(D)) mat(L))^*. $
  Here $sqrt(mat(D))$ denotes the diagonal matrix whose diagonal entries are obtained by taking the square root of those of $mat(D)$,
  which are necessarily real and positive because $mat(A)$ is positive definite.
  This implies the existence of a Cholesky factorization with $mat(C) = mat(L) sqrt(mat(D))$.
]

==== Calculation of the Choleski factor

The matrix $mat(C)$ can be calculated from~@eq:linear_cholesky_decomposition.
For example, developing the matrix product gives that $a_(1,1) = c_(1, 1)^2$ and so $c_(1, 1) = sqrt(a_(1,1))$.
It is then possible to calculate $c_(2,1)$ from the equation~$a_(2, 1) = c_(2, 1) c_(1, 1)$, and so on.
Implementing the Cholesky factorization is the goal of @exercise:linear_cholesky.

=== Direct methods for banded matrices

In applications related to partial differential equations,
the matrix $mat(A) in complex^(n times n)$ very often has a bandwidth which is small in comparison with $n$.
#definition[
  The bandwidth of a matrix $mat(A) in complex^(n times n)$ is the smallest number $k in nat$ such that
  $a_(i j) = 0$ for all $(i, j) in {1, dots, n}^2$ with $abs(i - j) > k$.
]
It is not difficult to show that,
if $mat(A)$ is a matrix with bandwidth $k$,
then so are $mat(L)$ and $mat(U)$ in the absence of pivoting.
This can be proved by equaling the entries of the product $mat(L) mat(U)$ with those of the matrix $mat(A)$.
We emphasize, however, that the sparsity structure _within_ the band of $mat(A)$ may be destroyed in $mat(L)$ and $mat(U)$;
this phenomenon is called _fill-in_.

==== Reducing the bandwidth: the Cuthill--McKee algorithm~#moreinfo <ssub:reducing_the_bandwidth_the_cuthill_mckee_algorithm>

The computational cost of calculating the $mat(L) mat(U)$ or Cholesky decomposition of a matrix with bandwidth $k$ scales as $bigo(n k^2)$,
which is much better than the general scaling $bigo(n^3)$ if $k << n$.
In applications, the bandwidth $k$ is often related to the matrix size $n$.
For example, if~$mat(A)$ arises from the discretization of the Laplacian operator, then $k = bigo(sqrt(n))$
provided that a good ordering of the vertices is employed.
In this case, the computational cost scales as $bigo(n^2)$.

Since a narrow band is associated with a lower computational cost of the $mat(L) mat(U)$ decomposition,
it is natural to wonder whether the bandwidth of a matrix $mat(A)$ can be reduced.
A possible strategy to this end is to use permutations.
More precisely, is it possible to identify a row permutation matrix $mat(P)$ such that
$mat(P) mat(A) mat(P)^T$
has minimal bandwidth?
Given such a matrix,
the solution to the linear system~@eq:linear_system can be obtained by first solving $(mat(P) mat(A) mat(P)^T) vect(y) = mat(P) b$,
and then letting $vect(x) = mat(P)^T vect(y)$.

The Cuthill--McKee algorithm is a heuristic method for finding a good,
but sometimes not optimal,
permutation matrix~$mat(P)$ in the particular case where $mat(A)$ is a _Hermitian_ matrix.
It is based on the fact that,
to a Hermitian matrix~$mat(A)$, we can associate a unique undirected graph whose adjacency matrix~$mat(A)_*$ has the same sparsity structure as that of $mat(A)$,
i.e.~zeros in the same places.
For any row permutation matrix $mat(P)_(sigma)$ with corresponding permutation~$sigma colon {1, dots, n} -> {1, dots, n}$ (see @definition:row_permutation_matrix),
the matrices~$mat(P)_(sigma) mat(A) mat(P)_(sigma)^T$ and~$mat(P)_(sigma) mat(A)_* mat(P)_(sigma)^T$ also have the same sparsity structure.
Therefore, minimizing the bandwidth of~$mat(P)_(sigma) mat(A) mat(P)_(sigma)^T$ is equivalent to minimizing the bandwidth of~$mat(P)_(sigma) mat(A)_* mat(P)_(sigma)^T$.
The key insight for understanding the Cuthill--McKee method is that $mat(P)_(sigma) mat(A)_* mat(P)_(sigma)^T$ is the adjacency matrix of the graph obtained by renumbering the nodes
according to the permutation $sigma$,
i.e.~by changing the number of the nodes from~$i$ to~$sigma(i)$.
Consider, for example,
the following graph and renumbering:

#block(width: 100%, move(dx: parindent - 0.5cm, grid(
  columns: (40%, 10%, 40%),
  column-gutter: (0.1cm, 0.7cm),
  align: center + horizon,
  pad(y: 0.5cm, resized(octagon-graph((1, 2, 3, 4, 5, 6, 7, 8)))),
  text(size: huge)[$->$],
  pad(y: 0.5cm, resized(octagon-graph((1, 2, 4, 6, 8, 7, 5, 3)))),
)))
The associated adjacency matrices are given by:
$
  pmat(
    1, 1, , , , , , red1;
    1, 1, 1, , , , , ;
    , 1, 1, 1, , , , ;
    , , 1, 1, 1, , , ;
    , , , 1, 1, 1, , ;
    , , , , 1, 1, 1, ;
    , , , , , 1, 1, 1;
    red1, , , , , , 1, 1;
  )
  wide
  #text(size: huge)[$->$]
  wide
  pmat(
    1, 1, 1, , , , , ;
    1, 1, , 1, , , , ;
    1, , 1, , 1, , , ;
    , 1, , 1, , 1, , ;
    , , 1, , 1, , 1, ;
    , , , 1, , 1, , 1;
    , , , , 1, , 1, 1;
    , , , , , 1, 1, 1;
  )
$
We assume that the nodes are all self-connected,
although this is not depicted,
and so the diagonal entries of the adjacency matrices are equal to 1.
This renumbering corresponds to the permutation
$
  pmat(
    i class("relation", :), 1, 2, 3, 4, 5, 6, 7, 8;
    sigma(i) class("relation", :), 1, 2, 4, 6, 8, 7, 5, 3
  ),
$
and we may verify that the adjacency matrix of the renumbered graph can be obtained from the associated row permutation matrix:
$
  mat(P) mat(A)_* mat(P)^T
  =
  pmat(
    1, , , , , , , ;
    , 1, , , , , , ;
    , , , , , , , 1;
    , , 1, , , , , ;
    , , , , , , 1, ;
    , , , 1, , , , ;
    , , , , , 1, , ;
    , , , , 1, , , ;
  )
  pmat(
    1, 1, , , , , , red1;
    1, 1, 1, , , , , ;
    , 1, 1, 1, , , , ;
    , , 1, 1, 1, , , ;
    , , , 1, 1, 1, , ;
    , , , , 1, 1, 1, ;
    , , , , , 1, 1, 1;
    red1, , , , , , 1, 1;
  )
  pmat(
    1, , , , , , , ;
    , 1, , , , , , ;
    , , , 1, , , , ;
    , , , , , 1, , ;
    , , , , , , , 1;
    , , , , , , 1, ;
    , , , , 1, , , ;
    , , 1, , , , , ;
  )
$

#indent In this example,
renumbering the nodes of the graph enables a significant reduction of the bandwidth,
from 7 to 2.
The Cuthill--McKee algorithm,
which was employed to calculate the permutation,
is an iterative method that produces an ordered $n$-tuple $R$ containing the nodes in the new order;
in other words, it returns $big((sigma^(-1)(1), dots, sigma^(-1)(n)))$.
The first step of the algorithm is to find the node $i$ with the lowest _degree_,
i.e.~with the smallest number of connections to other nodes,
and to initialize $R = (i)$.
Then the following steps are repeated until $R$ contains all the nodes of the graph:
- Define $A_i$ as the set containing all the nodes which are adjacent to a node in $R$ but not themselves in $R$;

- Sort the nodes in $A_i$ according to the following rules:
  a node~$i in A_i$ comes before~$j in A_i$ if~$i$ is connected to a node in $R$ that comes before all the nodes in $R$ to which $j$ is connected.
  As a tiebreak,
  precedence is given to the node with highest degree.

- Append the nodes in $A_i$ to $R$,
  in the order determined in the previous item.

#figure(
  placement: none,
  [#octagon-graph((1, none, none, none, none, none, none, none), factor: 0.6)
    #h(.5cm)
    #octagon-graph((1, 2, none, none, none, none, none, 3), red: (1, 7), factor: 0.7)
    #h(.5cm)
    #octagon-graph((1, 2, 4, none, none, none, 5, 3), red: (2, 6), factor: 0.7)
    #h(.5cm)
    #octagon-graph((1, 2, 4, 6, none, 7, 5, 3), red: (3, 5), factor: 0.7)
    #h(.5cm)
    #octagon-graph((1, 2, 4, 6, 8, 7, 5, 3), red: (4,), factor: 0.7)],
  caption: [
    Illustration of the Cuthill--McKee algorithm.
    The new numbering of the nodes is illustrated.
    The first node was chosen randomly since all the nodes have the same degree.
    In this example, the ordered tuple $R$ evolves as follows:
    $(1) -> (1, 2, 8) -> (1, 2, 8, 3, 7) -> (1, 2, 8, 3, 7, 4, 6) -> (1, 2, 8, 3, 7, 4, 6, 5)$.
  ],
) <figure:linear_cuthill_mckee>
The steps of the algorithm for the example above are depicted in~@figure:linear_cuthill_mckee.
Another example,
taken from the original paper by Cuthill and McKeen~@cuthill1969reducing,
is presented in~@figure:linear_example_from_cuthill_mckee.

#figure(
  placement: none,
  [#box(image("/figures/linear_cuthill-mckee_original.png", width: 49%)) #box(image("/figures/linear_cuthill-mckee_improved.png", width: 49%))],
  caption: [Example from the original Cuthill--McKee paper~@cuthill1969reducing.],
) <figure:linear_example_from_cuthill_mckee>

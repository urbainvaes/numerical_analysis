#import "../../template.typ": *

// \ddots and \vdots (`dots` is redefined in the template)
#let ddots = sym.dots.down
#let vdots = sym.dots.v
// algorithmic environment outside of an algorithm float: no caption, no rules
// and no line numbers. Each line is a pair (indentation level, content).
#let algorithmic(..lines) = block(width: 100%, above: leading + 9pt, below: leading + 9pt, {
  set par(first-line-indent: 0pt, leading: float-leading, spacing: float-leading)
  grid(row-gutter: float-leading, ..lines.pos().map(((level, body)) => pad(left: level * 1.5em, body)))
})

== Exercises

#exercise[
  In the simple case where $mat(A)$ is symmetric,
  find values of $vect(x)$, $vect(b)$ and $Delta vect(b)$ for which the inequality~@eq:linear_perturbation_rhs is in fact an equality?
] <exercise:linear_sharp_inequality>

#exercise(title: [Inverse of Gaussian transformation])[
  Prove the formula~@eq:inverse_gaussian_transformation.
] <exercise:inverse_gaussian_transformation>

#exercise[
  Prove that the product of two lower triangular matrices is lower triangular.
] <exercise:linear_product_of_lower_triangular>

#exercise[
  Assume that $mat(A) in real^(n times n)$ is positive definite,
  i.e.~that
  $ forall vect(x) in real^n backslash {vect(0)_n}, wide vect(x)^T mat(A) vect(x) > 0. $
  Show that all the principal submatrices of $mat(A)$ are nonsingular.
] <exercise:linear_positive_definite_matrix_nonsingular_principal_components>

#compexercise[
  Implement the backward substitution algorithm for solving $mat(U) x = y$.
  What is the computational cost of the algorithm?
]

#compexercise[
  Compare the condition number of the matrices $mat(L)$ and $mat(U)$ with and without partial pivoting.
  For testing, use a matrix with pseudo-random entries generated as follows
  ```julia
  import Random
  # Set the seed so that the code is deterministic
  Random.seed!(0)
  n = 1000 # You can change this parameter
  A = randn(n, n)
  ```
] <exercise:linear_lu_with_partial_pivoting>
#solution[
  See the Jupyter notebook for this chapter.
]

#compexercise[
  Write a code for calculating the Cholesky factorization of a symmetric positive definite matrix $mat(A)$ by comparing the entries of the product $mat(C) mat(C)^T$ with those of the matrix $mat(A)$.
  What is the associated computational cost,
  and how does it compare with that of the $mat(L) mat(U)$ factorization?

  #par(first-line-indent: 0pt)[*Extra credit:* ... if your code is able to exploit the potential banded structure of the matrix passed as argument for better efficiency.
  Specifically, your code will be tested with a matrix is of the type `BandedMatrix` defined in the #tt[BandedMatrices.jl] package,
  which you will need to install.
  The following code can be useful for testing purposes.]
  ```julia
  import BandedMatrices
  import LinearAlgebra

  function cholesky(A)
      m, n = size(A)
      m != n && error("Matrix must be square")
      # Convert to banded matrix
      B = BandedMatrices.BandedMatrix(A)
      B.u != B.l && error("Matrix must be symmetric")
      # --> Your code comes here <--
  end
  n, u, l = 20000, 2, 2
  A = BandedMatrices.brand(n, u, l)
  A = A*A'
  # so that A is symmetric and positive definite (with probability 1).
  C = @time cholesky(A)
  LinearAlgebra.norm(C*C' - A, Inf)
  ```
  For information, my code takes about 1 second to run with the parameters given here.
] <exercise:linear_cholesky>

#exercise(title: [Matrix square root])[
  Let $mat(A) in real^(n times n)$ be a symmetric positive definite matrix.
  Show that $mat(A)$ has a positive definite square root,
  i.e.~that there exists a symmetric matrix $mat(B)$ such that $mat(B) mat(B) = mat(A)$.
] <exercise:matrix_square_root>
#solution[
  Since $mat(A)$ is symmetric,
  there exist a diagonal matrix $mat(D)$ and an orthogonal matrix $mat(Q)$ such that $mat(A) = mat(Q) mat(D) mat(Q)^T$.
  Let $mat(D)^(1\/2)$ denote the diagonal matrix obtained by applying the square root function to the entries of~$mat(D)$,
  and notice that $mat(D)^(1\/2) mat(D)^(1\/2) = mat(D)$.
  Then it holds that
  $ mat(A) = big((mat(Q) mat(D)^(1\/2) mat(Q)^T)) big((mat(Q) mat(D)^(1\/2) mat(Q)^T)). $
  The matrix $mat(A)^(1\/2) := mat(Q) mat(D)^(1\/2) mat(Q)^T$ is a square root of the matrix~$mat(A)$,
  in the sense that $mat(A)^(1\/2) mat(A)^(1\/2)$,
  and it is positive definite because the diagonal elements of $mat(D)^(1\/2)$ are strictly positive.
]

#exercise[
  Show that if $mat(A)$ is row or column diagonally dominant,
  then $mat(A)$ is invertible.
] <exercise:invertibility_diagonal_dominant>

#exercise[
  Let $mat(T)$ be a nonsingular matrix.
  Show that
  $ norm(mat(A))_mat(T) := norm(mat(T)^(-1) mat(A) mat(T))_2 $
  defines a matrix norm induced by a vector norm.
] <exercise:induced_matrix_norm>

#exercise[
  Let $mat(A) in real^(n times n)$ be a symmetric positive definite matrix.
  Show that the functional
  $ norm(placeholder)_mat(A) colon vect(x) |-> sqrt(vect(x)^T mat(A) vect(x)) $
  defines a norm on $real^n$.
] <exercise:linear_norm_induced_A>
#solution[
  We need to prove that the three axioms of a norm are satisfied:
  - (*Positivity*)
    Since $mat(A)$ is positive definite,
    it holds that $norm(vect(x))_mat(A) > 0$ for any $vect(x) in real^n without {vect(0)}$.

  - (*Homogeneity*)
    It is clear that $norm(c vect(x))_mat(A) = |c| norm(vect(x))_mat(A)$ for any $c in real$.

  - (*Triangle inequality*)
    Let $mat(A)^(1\/2)$ denote the positive definite square root of~$mat(A)$,
    which exists by~@exercise:matrix_square_root.
    Then
    $ norm(vect(x))_mat(A) = norm(mat(A)^(1\/2) vect(x))_2. $
    The triangle inequality for $norm(placeholder)_mat(A)$ then follows from that for~$norm(placeholder)_2$:
    $
      norm(vect(x) + vect(y))_mat(A)
      = norm(mat(A)^(1\/2) vect(x) + mat(A)^(1\/2) vect(y))_2
      <= norm(mat(A)^(1\/2) vect(x)) + norm(mat(A)^(1\/2) vect(y))_2
      = norm(vect(x))_mat(A) + norm(vect(y))_mat(A).
    $
  Another option for solving this exercise is to show that
  $ ip(vect(x), vect(y))_mat(A) := vect(x)^T mat(A) vect(y) $
  defines an inner product,
  with induced norm given by~$norm(placeholder)_mat(A)$.
]

#exercise[
  Show that the residual satisfies the equation
  $ vect(r)^((k+1)) = mat(N) mat(M)^(-1) vect(r)^((k)) = (mat(I) - mat(A) mat(M)^(-1)) vect(r)^((k)). $
]

#exercise[
  Show that, if $mat(A)$ and $mat(B)$ are two square matrices,
  then $rho(mat(A) mat(B)) = rho(mat(B) mat(A))$.
]

#exercise[
  Is $rho(placeholder)$ a norm? Prove or disprove.
]

#exercise[
  Prove that, if $mat(A)$ is a diagonal matrix, then
  $ norm(mat(A))_1 = norm(mat(A))_2 = norm(mat(A))_infinity = rho(mat(A)). $
]

#exercise[
  Show that, for any matrix norm $norm(placeholder)$ induced by a vector norm,
  $ rho(mat(A)) <= norm(mat(A)). $
]

#exercise[
  Let $norm(placeholder)$ denote the Euclidean vector norm on $real^n$.
  We define in~@cha:vectors_and_matrices the induced matrix norm as
  $ norm(mat(A)) = sup big(\{ norm(mat(A) vect(x)) colon norm(vect(x)) <= 1 \}). $
  Show from this definition that, if $mat(A)$ is symmetric and positive definite,
  then
  $ norm(mat(A)) = norm(mat(A))_* := sup big(\{ abs(vect(x)^T mat(A) vect(x)) colon norm(vect(x)) <= 1 \}). $
  // What is the induced matrix norm on $\real^{n \times n}$?
]
#solution[
  By the Cauchy--Schwarz inequality and the definition of $norm(mat(A))$,
  it holds that
  $
    forall vect(x) in real^n "with" norm(vect(x)) <= 1, wide
    abs(vect(x)^T mat(A) vect(x))
    <= norm(vect(x)) norm(mat(A) vect(x))
    <= norm(vect(x)) norm(mat(A)) norm(vect(x)) <= norm(mat(A)).
  $
  This shows that $norm(mat(A))_* <= norm(mat(A))$.
  Conversely,
  letting $mat(B)$ denote a matrix square root of $mat(A)$ (see @exercise:matrix_square_root),
  we have
  $
    forall vect(x) in real^n "with" norm(vect(x)) <= 1, wide
    norm(mat(A) vect(x))
    &= sqrt(vect(x)^T mat(A)^T mat(A) vect(x))
    = sqrt((mat(B) vect(x))^T mat(B) mat(B) (mat(B) vect(x)))
    = sqrt((mat(B) vect(x))^T mat(A) (mat(B) vect(x))) \
    &= norm(mat(B) vect(x)) sqrt(vect(y)^T mat(A) vect(y)),
    wide vect(y) = (mat(B) vect(x))/norm(mat(B) vect(x)).
  $
  It holds that $norm(mat(B) vect(x)) = sqrt(vect(x)^T mat(A) vect(x)) <= sqrt(norm(mat(A))_*)$.
  In addition $norm(vect(y)) = 1$,
  so the expression inside the square root is bounded from above by $norm(mat(A))_*$,
  which enables to conclude the proof.
]


#exercise[
  Prove that, if the matrix $mat(A)$ is strictly diagonally dominant (by rows or columns),
  then the Gauss--Seidel method converges, i.e.~$rho(mat(M)^(-1) mat(N)) < 1$.
  You can use the same approach as in the proof of @proposition:linear_convergence_jacobi.
] <exercise:linear_convergence_gauss_seidel>

#exercise[
  Let $mat(A) in real^(n times n)$ denote a symmetric positive definite matrix,
  and assume that the vectors $vect(d)_1, dots, vect(d)_n$ are pairwise $mat(A)$-orthogonal directions.
  Show that $vect(d)_1, dots, vect(d)_n$ are linearly independent.
] <exercise:linear_independence_conjugate_directions>

#compexercise(title: [Steepest descent algorithm])[
  Consider the linear system
  #numbered($
    mat(A) vect(x) :=
    pmat(3, 1; 1, 3)
    pmat(x_1; x_2)
    =
    pmat(1; 1) =: vect(b).
  $) <eq:exercise_linear_system>
  - Show that $mat(A)$ is positive definite.

  - Draw the contour lines of the function
    $ f(vect(x)) = 1/2 vect(x)^T mat(A) vect(x) - vect(b)^T vect(x). $

  - Plot the contour lines of~$f$ in Julia using the function `contourf` from the package `Plots`.

  - Using @theorem:linear_convergenec_steepest_descent,
    estimate the number~$K$ of iterations of the steepest descent algorithm required in order to guarantee that $E_K <= 10^(-8)$,
    when starting from the vector $vect(x)^((0)) = (2 #h(0.33em) 3)^T$.

  - Implement the steepest descent method for finding the solution to~@eq:exercise_linear_system,
    and plot the iterates as linked dots over the filled contour of~$f$.

  - Plot the error $E_k$ as a function of the iteration index,
    using a linear scale for the $x$ axis and a logarithmic scale for the $y$ axis.
]

#exercise[
  Compute the number of floating point operations required for performing one iteration of the conjugate gradient method,
  assuming that the matrix $mat(A)$ contains $alpha << n$ nonzero elements per row.
]

#compexercise(title: [Solving the Poisson equation over a rectangle])[
  We consider in this exercise Poisson's equation in the domain $Omega = (0, 2) times (0, 1)$,
  equipped with homogeneous Dirichlet boundary conditions:
  $
    - laplacian f(x, y) &= b(x, y), wide x in Omega, \
    f(x) &= 0, quad wide x in partial Omega.
  $
  The right-hand side is
  $ b(x, y) = sin(4 pi x) + sin(2 pi y). $
  A number of methods can be employed in order to discretize this partial differential equation.
  After discretization, a finite-dimensional linear system of the form~$mat(A) vect(x) = vect(b)$ is obtained.
  // and solving this system is often the most computationally expensive part of approximating the solution~$f$.
  A Julia function for calculating the matrix $mat(A)$ and the vector $vect(b)$ using the finite difference method is given to you on the course website,
  as well as a function to plot the solution.
  The goal of this exercise is to solve the linear system using the conjugate gradient method.
  Use the same stopping criterion as in~@exercise:one_dim_finite_difference.
]

#exercise[
  Show that if $mat(A) in real^(n times n)$ is nonsingular,
  then the solution to the equation $mat(A) vect(x) = vect(b)$ belongs to the Krylov subspace
  $
    cal(K)_n (mat(A), vect(b))
    = Span Big(\{) vect(b), mat(A) vect(b), mat(A)^2 vect(b), dots, mat(A)^(n-1) vect(b) Big(\}).
  $
]

#exercise[
  Write a function `lu(A)` for calculating the $mat(L) mat(U)$ decomposition of a square matrix $mat(A) in real^(n times n)$,
  with $mat(L)$ unit lower triangular and $mat(U)$ upper triangular,
  not by Gaussian elimination but by comparing the entries of the product $mat(L) mat(U)$ with those of $mat(A)$.
  To this end, one option is to compare the entries one by one
  in the order $(1, 1)$, $(1, 2)$, #sym.dots, $(1, n)$, $(2, 1)$, $(2, 2)$, #sym.dots,
  i.e.~row by row starting from the top.
  For example,
  - Comparing the entry $(1, k)$ with $k in {1, dots, n}$ gives
    $ ell_11 u_(1k) = a_(1k). $
    Since $ell_11 = 1$ as $mat(L)$ is unit lower triangular,
    this implies that $u_(1k) = a_(1k)$.

  - Comparing the entry $(2, 1)$ gives
    $ ell_21 u_11 = a_21 $
    and so $ell_21 = a_21 \/ u_11$.

  - Comparing the entry $(2, k)$ with $k in {2, dots, n}$ gives
    $ ell_21 u_(1k) + ell_22 u_(2k) = a_(2k). $
    Given the previous items,
    the only unknown in this equation is~$u_(2k)$.

  - Comparing the entry $(3, 1)$ gives
    $ ell_31 u_11 = a_31, $
    and so $ell_31 = a_31 \/ u_11$.

  - Comparing the entry $(3, 2)$ gives
    $ ell_31 u_12 + ell_32 u_22 = a_32, $
    Given the previous items,
    the only unknown in this equation is~$ell_32$.

  - Comparing the entry $(3, k)$ with $k in {3, dots, n}$ gives
    $ ell_31 u_(1k) + ell_32 u_(2k) + ell_33 u_(3k) = a_(3k), $
    Given the previous items,
    the only unknown in this equation is~$u_(3k)$.
  Notice that a pattern seems to be emerging:
  when going through the entries row by row starting from the top left corner of the matrix,
  comparing the entry $(i, j)$ provides an equation for $ell_(i j)$ if~$j < i$,
  and an equation for $u_(i j)$ if $j >= i$.
  Do not use any external package for this exercise.

  #par(first-line-indent: 0pt)[*Extra credit:* ... if your code is able to exploit the potential banded structure of the matrix passed as argument for better efficiency.
  Specifically, your code will be tested with a matrix created as follows]
  ```julia
  b, n = 5, 10000
  A = [abs(i-j) <= b ? rand() : 0.0 for i in 1:n, j in 1:n]
  ```
]

#compexercise[
  Implement an iterative method based on a splitting for finding a solution to the following linear system on $real^n$.
  $
    1/h^2
    pmat(
      2, -1;
      -1, 2, -1;
      , -1, 2, -1;
      , , ddots, ddots, ddots;
      , , , -1, 2, -1;
      , , , , -1, 2;
    )
    pmat(x_1; x_2; x_3; vdots; x_(n-1); x_n)
    =
    pmat(1; 1; 1; vdots; 1; 1),
    wide
    h = 1/(n+1).
  $
  Plot the norm of the residual as a function of the iteration index.
  Use as stopping criterion the condition
  $
    norm(vect(r)^((k))) <= epsilon norm(vect(b)),
    wide epsilon = 10^(-8).
  $
  As initial guess, use a vector of zeros.
  The code will be tested with $n = 500$.
  Do not use any library (except for plotting),
  and do not use the backslash operator.
] <exercise:one_dim_finite_difference>

#exercise[
  Find a formula for the optimal value of $omega$ in the relaxation method given~$n$,
  for the linear system in~@exercise:one_dim_finite_difference.
  The proof of @proposition:tridiagonal,
  as well as the formula~@eq:linear_eigenvalues_toeplitz for the eigenvalues of a tridiagonal matrix,
  are useful to this end.
]

#solution[
  @corollary:convergence_relaxation and @proposition:necessary_sor imply that a sufficient and necessary condition for convergence,
  when $mat(A)$ is Hermitian and positive definite, is that $omega in (0, 2)$.
  Let $M_omega = 1/omega mat(D) + mat(L)$ and $N_omega = (1-omega)/omega mat(D) - mat(U)$.
  A nonzero scalar $lambda in complex$ is an eigenvalue of $mat(M)_omega^(-1) mat(N)_omega$ if and only if
  $
    det(mat(M)_omega^(-1) mat(N)_omega - lambda mat(I)) = 0
    quad <=> quad
    det(mat(M)_omega^(-1)) det(mat(N)_omega - lambda mat(M)_omega) = 0
    quad <=> quad
    det(lambda mat(M)_omega - mat(N)_omega) = 0.
  $
  Substituting the expressions of $mat(M)_omega$ and $mat(N)_omega$,
  we obtain that this condition can be equivalently rewritten as
  $
    det lr((lambda mat(L) + lr(((lambda + omega - 1)/omega)) mat(D) + mat(U))) = 0
    quad <=> quad
    det lr((sqrt(lambda) mat(L) + lr(((lambda + omega - 1)/omega)) mat(D) + sqrt(lambda) mat(U))) = 0
  $
  where we used~@eq:preliminary_equation for the last equivalence.
  The equality of the determinants in these two equations is valid for
  $sqrt(lambda)$ denoting either of the two complex square roots of~$lambda$.
  This condition is equivalent to
  $
    det lr((mat(L) + lr(((lambda + omega - 1)/(sqrt(lambda) omega))) mat(D) + mat(U))) = 0.
  $
  We recognize from the proof of @proposition:tridiagonal
  that this condition is equivalent to
  $
    (lambda + omega - 1)/(sqrt(lambda) omega) in spectrum(mat(M)_cal(J)^(-1) mat(N)_cal(J)).
  $
  In other words, for any $(lambda, mu) in complex^2$ such that
  #numbered($
    (lambda + omega - 1)^2/(lambda omega^2) = mu^2,
  $) <eq:relationship_lambda_mu>
  it holds that $mu in spectrum(mat(M)_cal(J)^(-1) mat(N)_cal(J))$ if and only if $lambda in spectrum(mat(M)_omega^(-1) mat(N)_omega)$.
  By~@eq:linear_eigenvalues_toeplitz,
  the eigenvalues of $mat(M)_cal(J)^(-1) mat(N)_cal(J)$ are real and given by
  #numbered($
    mu_j = cos lr(((j pi)/(n+1))),
    wide 1 <= j <= n.
  $) <eq:explicit_eigenvalues_jacobi>
  Rearranging~@eq:relationship_lambda_mu,
  we find
  $
    lambda^2 + lambda lr((2(omega - 1) - omega^2 mu^2)) + (omega - 1)^2 = 0.
  $
  For given $omega in (0, 2)$ and $mu in real$,
  this is a quadratic equation for $lambda$ with solutions
  $
    lambda_plus.minus = lr(((omega^2 mu^2)/2 + 1 - omega)) plus.minus omega mu sqrt((omega^2 mu^2)/4 + 1 - omega),
  $
  Since the first bracket is positive when the argument of the square root is positive,
  it is clear that
  $
    max big(\{ abs(lambda_-)\, abs(lambda_+) \})
    = lr(|(omega^2 mu^2)/2 + 1 - omega + omega abs(mu) sqrt((omega^2 mu^2)/4 + 1 - omega)|).
  $
  Combining this with~@eq:explicit_eigenvalues_jacobi,
  we deduce that the spectral radius of $mat(M)_omega^(-1) mat(N)_omega$ is given by
  #numbered($
    rho(mat(M)_omega^(-1) mat(N)_omega) =
    max_(j in {1, dots, n}) lr(|(omega^2 mu_j^2)/2 + 1 - omega + omega abs(mu_j) sqrt((omega^2 mu_j^2)/4 + 1 - omega)|).
  $) <eq:sor_spectral_radius>
  We wish to minimize this expression over the interval $omega in (0, 2)$.
  While this can be achieved by algebraic manipulations,
  we content ourselves here with graphical exploration.
  @fig:modulus_lambdas depicts the amplitude of the modulus in~@eq:sor_spectral_radius
  for different values of $mu$.
  It is apparent that, for given $omega$,
  the modulus increases as $mu$ increases,
  which suggests that
  #numbered($
    rho(mat(M)_omega^(-1) mat(N)_omega) =
    lr(|(omega^2 mu_*^2)/2 + 1 - omega + omega abs(mu_*) sqrt((omega^2 mu_*^2)/4 + 1 - omega)|),
    wide mu_* = rho(mat(M)_cal(J)^(-1) mat(N)_cal(J)).
  $)
  The figure also suggests that for a given value of $mu$,
  the modulus is minimized at the discontinuity of the first derivative,
  which occurs when the argument of the square root is zero.
  We conclude that the optimal $omega$ satisfies
  $
    (omega_"opt"^2 mu_*^2)/4 + 1 - omega_"opt" = 0
    quad stretch(=>)_(omega < 2) quad
    omega_"opt" = 2 (1 - sqrt(1 - mu_*^2))/mu_*^2
    = 2/(1 + sqrt(1 - mu_*^2))
    = 2/(1 + sin lr((pi/(n+1)))).
  $
]
#figure(
  placement: top,
  image("/figures/optimal_omega.pdf", width: 80%),
  caption: [
    Modulus of $abs(lambda_+)$ as a function of $omega$,
    for different eigenvalues of $mu$.
  ],
) <fig:modulus_lambdas>


#exercise(title: [Midterm 2022])[
  Let $mat(A) in real^(n times n)$ be a symmetric positive definite matrix and let~$vect(b) in real^n$.
  The steepest descent algorithm for solving $mat(A) vect(x) = vect(b)$ is given hereafter:

  #algorithmic(
    (0, [Pick $epsilon > 0$ and initial $vect(x)$]),
    (0, [$vect(r) <- mat(A) vect(x) - vect(b)$]),
    (0, [*while* $norm(vect(r)) >= epsilon norm(vect(b))$ *do*]),
    (1, [$omega <- vect(r)^T vect(r) \/ vect(r)^T mat(A) vect(r)$]),
    (1, [$vect(x) <- vect(x) - omega vect(r)$]),
    (1, [$vect(r) <- mat(A) vect(x) - vect(b)$]),
    (0, [*end while*]),
  )

  - Why is this method called the _steepest descent_ algorithm?

  - How many floating point operations does an iteration of this algorithm require?

  - Are the following statements true of false? (2 marks)
    + There exists a unique solution~$vect(x)_*$ to the linear system~$mat(A) vect(x) = vect(b)$.

    + The iterates converge to $vect(x)_*$ in at most $n$ iterations.

    + We consider the following modification of the algorithm:
      #algorithmic(
        (0, [Pick $epsilon > 0$, $omega > 0$ and initial $vect(x)$]),
        (0, [$vect(r) <- mat(A) vect(x) - vect(b)$]),
        (0, [*while* $norm(vect(r)) >= epsilon norm(vect(b))$ *do*]),
        (1, [$vect(x) <- vect(x) - omega vect(r)$]),
        (1, [$vect(r) <- mat(A) vect(x) - vect(b)$]),
        (0, [*end while*]),
      )
      If $omega$ is sufficiently small, then this algorithm converges.

    + Here we no longer assume that $mat(A)$ is positive definite.
      Instead, we consider that
      $ mat(A) = pmat(-1, 0; 0, -2). $
      In this case, the steepest descent algorithm is convergent for any initial $vect(x)$.
]

#exercise(title: [Final exam Spring 2022])[
  Assume that $mat(A) in real^(n times n)$ is a nonsingular matrix and that $vect(b) in real^n$.
  We wish to solve the linear system~@eq:linear_system
  using an iterative method where each iteration is of the form
  #numbered($ mat(M) vect(x)_(k+1) = mat(N) vect(x)_k + vect(b). $) <eq:iterative_scheme>
  Here $mat(A) = mat(M) - mat(N)$ is a splitting of~$mat(A)$ such that $mat(M)$ is nonsingular,
  and $vect(x)_k in real^n$ denotes the $k$-th iterate of the numerical scheme.

  + Let $vect(e)_k := vect(x)_k - vect(x)_*$,
    where $vect(x)_*$ is the exact solution to~@eq:linear_system.
    Prove that
    $ vect(e)_(k+1) = mat(M)^(-1) mat(N) vect(e)_k. $

  + Let $L = norm(mat(M)^(-1) mat(N))_infinity$.
    Prove that
    #numbered($
      forall k in nat, wide
      norm(vect(e)_k)_infinity <= L^k norm(vect(e)_0)_infinity.
    $) <eq:inequality_interpolation>

  + Is the condition $norm(mat(M)^(-1) mat(N))_infinity < 1$ necessary
    for convergence when $vect(x)_0 != vect(x)_*$?

  + Assume that $mat(A)$ is strictly row diagonally dominant, in the sense that
    $
      forall i in {1, dots, n}, wide
      |a_(i i)| > sum_(j=1, j != i)^n |a_(i j)|.
    $
    Show that, in this case, the inequality $norm(mat(M)^(-1) mat(N))_infinity < 1$ holds for the Jacobi method,
    i.e.~when $mat(M)$ contains just the diagonal of~$mat(A)$.
    You may take for granted the following expression for the~$infinity$-norm of a matrix $mat(X) in real^(n times n)$:
    $ norm(mat(X))_infinity = max_(1 <= i <= n) sum_(j=1)^n abs(x_(i j)). $

  + Write down a few iterations of the Jacobi method when
    $
      mat(A) =
      pmat(1, 2; 0, 1),
      wide
      vect(b)
      pmat(1; 1),
      wide
      vect(x)_0 =
      pmat(0; 0).
    $
    Is the method convergent?
]

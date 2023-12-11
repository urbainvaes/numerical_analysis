# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.5
# ---

# # NYU Paris - Numerical Analysis

# ## Final Exam
#
# - Submit this notebook on Brightspace before 6:00 PM
#
# - The exam comprises six independent exercises. In each exercise, the
#   cells may depend on the previous cells.
#
# - To facilitate correction of your code,
#   do not change the signatures of the functions to implement.
#
# - Run the cell below to import necessary packages and
#   load a macro used for unit tests thoughout the notebook.

# +
using ForwardDiff
using LinearAlgebra
using Polynomials
using SpecialFunctions
using TestImages
using LaTeXStrings
using Plots

macro mark(bool_expr)
    return :( print($bool_expr ? "✔️" : "❌"))
end
# -

# ### 1. Floating point arithmetic
#
# Throughout this exercise,
# $(\bullet)_\beta$ denotes base $\beta$ representation.
# True or false? (+1/0/-1)
# 1. It holds that
#    $
#      (1.0111)_2 + (0.1001)_2 = (10)_2.
#    $
# 1. It holds that
#    $
#      (1000)_5 \times (0.003)_5 = 3.
#    $
# 1. It holds that
#    $
#      (0.\overline{4})_5 = 4
#    $
# 1. In base 16, all the natural numbers from 1 to 1000 can be represented using 3 digits.
#
# 1. Machine multiplication according to the IEEE754 standard is a commutative operation.
#
# 1. Machine addition according to the IEEE754 standard is an associative operation.
#
# 1. Only finitely many real numbers can be represented exactly in the `Float64` format.
#
# 1. The spacing (in absolute value) between successive `Float64` numbers is strictly increasing.
#
# 1. In Julia, `eps(Float32)` is the smallest positive number representable in the `Float32` format.
#
# 1. In Julia, `nextfloat(1f0) - 1f0` equals the machine epsilon for the `Float32` format.
#
# 1. Assume that $x \in \real$ belongs to $\mathbf F_{64}$. Then $2x \in \mathbf F_{64}$.
#
# 1. The real number $\sqrt{2}$ can be represented exactly in the single-precision format.

# +
answers_q1 = zeros(Int, 12)

# Use -1 for false, 1 for true
answers_q1[1] = 0
answers_q1[2] = 0
answers_q1[3] = 0
answers_q1[4] = 0
answers_q1[5] = 0
answers_q1[6] = 0
answers_q1[7] = 0
answers_q1[8] = 0
answers_q1[9] = 0
answers_q1[10] = 0
answers_q1[11] = 0
answers_q1[12] = 0
# -


# ### 2. Exercise on interpolation and approximation

# 1. Find the quadratic polynomial $p(x) = ax^2 + b x + c$ that goes through the points $(0, 1)$, $(1, 3)$ and $(2, 7),$
#    and then plot the polynomial $p$ together with the data points on the same graph.

# +
# Calculate a, b and c here

@mark begin p(x) = a*x^2 + b*x + c; p(0) == 1 end
@mark begin p(x) = a*x^2 + b*x + c; p(1) == 3 end
@mark begin p(x) = a*x^2 + b*x + c; p(2) == 7 end

# Create the plot here
# -

# 2. Find the function $f$ of the form $f(x) = a \sin(x) + b \cos(x) + c \sin(2x)$ that goes through the points $(x_1, y_1)$, $(x_2, y_2)$ and $(x_3, y_3)$,
#    where the data is given below.
#    Plot the function $f$ together with the data points on the same graph.

# +
x = [1, 2, 3]
y = [2.8313730233698577, -0.6797987415765313, -2.11828048333995]

# Calculate a, b and c here

@mark begin f(z) = a*sin(z) + b*cos(z) + c*sin(2z); f(x[1]) == y[1] end
@mark begin f(z) = a*sin(z) + b*cos(z) + c*sin(2z); f(x[2]) == y[2] end
@mark begin f(z) = a*sin(z) + b*cos(z) + c*sin(2z); f(x[3]) == y[3] end

# Create the plot here
# -

# 3. Find the function $g$ of the form $g(x) = a \sin(x) + b \cos(x)$
#    such that the sum of squared residuals
#    $$
#    \sum_{i=1}^{n} |g(x_i) - y_i|^2, \qquad n = 10,
#    $$
#    is minimized for the data given below.
#    Plot on the same graph the function $g$ and the data points.

# +
x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
y = [0.46618950237090034, 0.9365762460801551, 0.5907432672662498,
     -0.21370874419278404, -0.8869010982313607, -0.6906040605442342,
     0.1784931250350807, 1.024036713535387, 0.7837248688922458,
     -0.1544083454499539]

# Calculate a and b and create the plot here
# -

# ### 3. Exercise on numerical integration
#
# Our goal in this exercise is to write a program to calculate integrals of the form
# $$
# I[f] := \int_{-\infty}^{\infty} f(x) \, e^{-x^2} \, d x.
# $$
# To this end, we will use Hermite polynomials,
# which are orthogonal polynomials for the following inner product:
#    $$
#     \langle f, g \rangle := \int_{-\infty}^{\infty} f(x) \, g(x) \, e^{-x^2} \, d x
#    $$
# 1. Using that Hermite polynomials satisfy the recurrence relation
#    $$
#        H_{k + 1}(x) = 2x H_k(x) - H_k'(x), \qquad H_0(x) = 1,
#    $$
#    write a function `hermite(n)` which returns the Hermite polynomial of degree $n$.

# +
function hermite(n)
    # Write your code here
end

@mark hermite(2) == Polynomial([-2, 0, 4])
@mark hermite(3) == Polynomial([0, -12, 0, 8])
@mark hermite(4) == Polynomial([12, 0, -48, 0, 16])
@mark degree(hermite(10)) == 10
# -

# 2. Write a function `integrate_monomial(n)`, which return the value of the integral
#    $$
#    \int_{-\infty}^{\infty} x^n e^{-x^2} \, dx =
#    \frac{1}{2} \Bigl((-1)^n + 1\Bigr) \Gamma\left( \frac{n+1}{2} \right),
#    $$
#    where $\Gamma$ is the Gamma function,
#    available as `gamma` in the `SpecialFunctions` library.

# +
function integrate_monomial(n)
    # Write your code here
end

@mark integrate_monomial(0) ≈ √π
@mark integrate_monomial(2) ≈ √π/2
@mark integrate_monomial(3) ≈ 0
@mark integrate_monomial(4) ≈ 3√π/4
@mark integrate_monomial(5) == 0
# -

# 3. Write a function `get_nodes_and_weights(n)` which returns the nodes and weights of the Gauss-Hermite quadrature with $n$ nodes.
#
#    **Hints:**
#    - Remember that the nodes of the quadrature are given by the roots of $H_{n}$.
#    - In order to construct Lagrange polynomials, you may find it useful to use the `fromroots` function.

# +
function get_nodes_and_weights(n)
    # Write your code here
    return nodes, weights
end

@mark length(get_nodes_and_weights(1)) == 2
@mark length(get_nodes_and_weights(3)[1]) == 3
@mark length(get_nodes_and_weights(5)[2]) == 5
# -

# 4. Write a function `integrate_hermite(f, n)`, which returns an approximation of $I[f]$ obtained by Gauss-Hermite integration with $n$ nodes.

# +
function integrate_hermite(f, n)
    # Write your code here
end

@mark integrate_hermite(x -> 1, 10) ≈ integrate_monomial(0)
@mark abs(integrate_hermite(x -> x, 10)) <= 1e-10
@mark integrate_hermite(x -> x^2, 10) ≈ integrate_monomial(2)
@mark integrate_hermite(x -> x^4, 10) ≈ integrate_monomial(4)
@mark integrate_hermite(x -> x^6, 10) ≈ integrate_monomial(6)
@mark abs(integrate_hermite(sin, 10)) <= 1e-10
# -

# 5. Set $n = 4$, and calculate numerically the degree of precision of the corresponding quadrature rule.

# +
# Write your code here
# -

# 6. Set $f(x) = \cos(x)$, and plot the integration error as a function of $n$.
#    You may use that
#    $$
#    \int_{-\infty}^{\infty} \cos(x) \, e^{-x^2} \, d x = \frac{\sqrt{\pi}}{\sqrt[4]{e}}
#    $$

# +
# Write your code here
# -

# ### 4. Exercise on linear systems
#
# This exercise focuses on solving the Euler-Bernoulli beam equation in one dimension,
# with homogeneous Dirichlet boundary conditions:
# $$
# u''''(x) = 1, \qquad u(0) = u'(0) = u'(1) = u(1) = 0.
# $$
# This equation models the deflection of a beam clamped at both ends,
# under a uniform load.
# A discretization of this equation on a uniform grid using the finite difference method leads to the following linear system:
# $$
# \begin{pmatrix}
#     6  & -4 & 1 \\
#     -4 & 6  & -4 & 1 \\
#      1 & -4 & 6      & -4 & 1 \\
#      & 1 & -4 & 6      & -4 & 1 \\
#      &  &  \ddots  & \ddots & \ddots & \ddots & \ddots  \\
#      & &    &   1    & -4    & 6      & -4 & 1 \\
#      & & &    &   1    & -4    & 6      & -4 \\
#      & & &    &        &  1  & -4      & 6 \\
# \end{pmatrix}
# \begin{pmatrix}
#     u_2 \\
#     u_3 \\
#     u_4 \\
#     u_5 \\
#     \vdots \\
#     u_{n-4} \\
#     u_{n-3} \\
#     u_{n-2}
# \end{pmatrix}
# =
# h^4
# \begin{pmatrix}
#     1 \\
#     1 \\
#     1 \\
#     1 \\
#     \vdots \\
#     1 \\
#     1 \\
#     1
# \end{pmatrix},
# \quad
# h := \frac{1}{n},
# \quad
# x_i := ih.
# $$
# where $h$ is the spacing between the discretization points, and $(u_2, u_3, \cdots, u_{n-3}, u_{n-2})$ are the values of the unknown function $u$ at the points $(x_2, x_3, \cdots, x_{n-3}, x_{n-2})$.
#
# 1. Write a function `build_matrix(n)`, which returns the matrix of the linear system.
#
#    **Hint:** the function `diagm` is useful here.

# +
function build_matrix(n)
    # Write your code here
end

@mark size(build_matrix(20)) == (17, 17)
@mark build_matrix(20)[5, 3] ≈ 1
@mark build_matrix(20)[5, 4] ≈ -4
@mark build_matrix(20)[5, 5] ≈ 6
@mark build_matrix(20)[5, 6] ≈ -4
@mark build_matrix(20)[5, 7] ≈ 1
@mark build_matrix(20)[5, 8] ≈ 0
# -

# 2. Write a function `build_rhs(n)`, which returns the right-hand side of the linear system.

# +
function build_rhs(n)
    # Write your code here
end

@mark length(build_rhs(20)) == 17
@mark build_rhs(20)[end] == 1 / 20^4
@mark build_rhs(20)[1] == 1 / 20^4
# -

# 3. Write a function `forward_substitution(L, y)`
#    that solves the linear system $\mathsf L \mathbf x = \mathbf y$,
#    in the case where `L` is a lower-triangular matrix,
#    by using forward substitution.

# +
function forward_substitution(L, y)
    # Write your code here
end

@mark begin n = 10; L = [1 0; 2 1]; forward_substitution(L, [1; 0]) ≈ [1; -2] end
@mark begin n = 10; L = LowerTriangular(ones(n, n)); y = fill(2, n); forward_substitution(L, L*y) ≈ y end
@mark begin n = 10; L = LowerTriangular(randn(n, n)); y = randn(n); forward_substitution(L, L*y) ≈ y end
@mark begin n = 20; L = LowerTriangular(2ones(n, n)); y = rand(n); forward_substitution(L, L*y) ≈ y end
# -

# 4. The successive over-relaxation method is a splitting method for solving linear systems of the form $\mathsf A \mathbf x = \mathbf b$.
#    It is based on the iteration
#    $$
#    \mathsf M \mathbf x_{k+1} = \mathsf N \mathbf x_{k} + \mathbf{b}, \qquad \text{where} \quad \mathsf M = \frac{\mathsf D}{\omega} + \mathsf L \quad \text{and} \quad \mathsf N = \mathsf M - \mathsf A.
#    \tag{SOR}
#    $$
#    <a id="SOR"></a>
#    Here $\omega \in (0, 2)$ is a parameter,
#    $\mathsf D$ is diagonal matrix containing the diagonal of $\mathsf A$,
#    and $\mathsf L$ is a strictly lower triangular matrix containing the strict lower triangular part of $\mathsf A$,
#    not including the diagonal.
#    Write a function `sor(A, b, ω)` that
#    implements this iteration, without using Julia's `\` and `inv` functions.
#    Initialize the iteration with $\mathbf x_0 = \mathbf 0$,
#    and use as stopping criterion that
#    $$
#    \lVert \mathsf A \mathbf x - \mathbf b \rVert
#    \leq \varepsilon \lVert \mathbf b \rVert,
#    \qquad \varepsilon := 10^{-10}.
#    $$
#    If after `maxiter` iterations,
#    a solution that satisfies this stopping criterion has not been found,
#    return `nothing`.
#
#    **Hints**:
#    - The functions `Diagonal` and `LowerTriangular` are useful to construct the matrices $\mathsf D$ and $\mathsf L$.
#    - In order to solve <a href="#NR">(SOR)</a> at each iteration,
#      use the function `forward_substitution` you wrote above.

# +
function sor(A, b; ω = 1, ε = 1e-10, maxiter = 10^5)
    # Write your code here
end

# Test code
n = 20
X = qr(randn(n, n)).Q
A = X*Diagonal(rand(n) .+ 5)*X'
b = randn(n)
@mark begin ε = 1e-10; sol = sor(A, b; ω = 1.5, ε = ε); norm(A*sol - b) ≤ ε*norm(b) end
@mark begin ε = 1e-10; sol = sor(A, b; ω = 1.5, ε = ε); norm(A*sol - b) ≥ 1e-15*norm(b) end
@mark begin ε = 1e-12; sol = sor(A, b; ω = 1.5, ε = ε); norm(A*sol - b) ≤ ε*norm(b) end
@mark begin ε = 1e-12; sor(A, b; ω = 2, ε = ε) == nothing end
@mark begin ε = 1e-12; sor(A, b; ω = 1, ε = ε) ≈ A\b end
# -

# 5. Use the relaxation method implemented in the previous item
#    with parameters $\omega = 1.5$ and $\varepsilon = 10^{-8}$ to solve the linear system of this exercise with $n = 40$.
#    Compare on a graph the solution you obatined with the exact solution,
#    which is given by
#    $$ u(x) = \frac{1}{24} x^2(x - 1)^2. $$

# +
u(x) = 1/24 * x^2 * (x - 1)^2
# Write your code here
# -


# ### 5. Exercise on nonlinear equations
#
# We wish to find all the solutions to the following transcendental equation for $x \in [-5, -5]$.
# $$
# x = 5 \sin(\pi x)
# \tag{1}
# $$
#
# 1. Plot the functions $x \mapsto x$ and $x \mapsto 5 \sin(\pi x)$ on the same graph,
#    for $x$ in the range $[-5, 5]$,
#    and count the number of solutions

# +
# Write your code here
# -

# 2. Using the bisection method, calculate precisely the only solution in the interval $[\frac{1}{2}, 1]$.

# +
function bisection(f, a, b, ε = 1e-10)
    # Write your code here
end

# Calculate solution here
# -

# 3. Write a function `newton_raphson(f, x₀; maxiter = 100, ε = 1e-12)` that returns the result of the Newton-Raphson iteration applied to the equation $f(x) = 0$ and initialized at `x₀`,
#    or `nothing` if a solution is not found after `maxiter` iterations.
#    Use the `ForwardDiff` library to compute derivatives,
#    and use the following stopping criterion
#    $$
#    |f(x_k)| ≤ \varepsilon.
#    $$

# +
function newton_raphson(f, x₀; maxiter = 100, ε = 1e-12)
    # Write your code here
end

@mark newton_raphson(x -> x^2 - 2, 1) ≈ √2
@mark newton_raphson(x -> x^2 - 2, -1) ≈ -√2
@mark newton_raphson(x -> x^3 - 2, 1) ≈ cbrt(2)
@mark newton_raphson(x -> x^2 + 1, 1.5) == nothing
@mark newton_raphson(x -> x^2 + 1, 1) == nothing
# -

# 4. Using the Newton-Raphson method you implemented,
#    calculate all the solutions to the transcendental equation <a href="#NR">(1)</a>.

# +
# Write your code here.
# -

# 5. We now consider another approach,
#    called the *secant method*.
#    An iteration of this method is given by
#    $$
#    x_{k+2} = \frac{x_{k} f(x_{k+1}) - x_{k+1}f(x_k)}{f(x_{k+1}) - f(x_{k})}.
#    $$
#    Write a function `secant(f, x₀, x₁; maxiter = 100, ε = 1e-12)`
#    that returns the result of the secant iteration applied to the equation $f(x) = 0$, and initialized with `x₀` and `x₁`,
#    or `nothing` if a solution is not found after `maxiter` iterations.
#    Use the same stopping criterion as above.

# +
function secant(f, x₀, x₁; maxiter = 100, ε = 1e-12)
    # Write your code here
end

@mark secant(x -> x^2 - 2, 1., 2.) ≈ √2
@mark secant(x -> x^2 - 2, -1., -2.) ≈ -√2
@mark secant(x -> x^2 + 1, 1.,  2.) == nothing
# -

# 6. Use the secant method you implemented to solve the transcendental equation
# $$x = e^{-x}$$

# +
# Write your code here.
# -

# ### 6. Exercise on eigenvalue problems

# Our goal in this exercise is to implement,
# without using the functions `eigen` and `svd` from the `LinearAlgebra` library,
# a simple algorithm for calculating the dominant *singular values* and associated *singular vectors* of a matrix,
# and to use this algorithm for image compression.
# Recall that any matrix $\mathsf A \in \mathbf R^{n \times n}$ admits a *singular value decomposition* of the form
# $$
# \mathsf A = \mathsf U \mathsf \Sigma \mathsf V^\top,
# $$
# where $\mathsf U \in \mathbf R^{n \times n}$ and $\mathsf V \in \mathbf R^{n \times n}$ are orthogonol matrices,
# and $\mathsf \Sigma \in \mathbf R^{n \times n}$ is a diagonal matrix.
# The diagonal entries of $\mathsf \Sigma$ are called *singular value*,
# while the columns of $\mathsf U$ and $\mathsf V$ are the associated left and right *singular vectors*, respectively.
# We restrict our attention to square matrices, for simplicity.
#
# 1. Write a function `myEigen(B, p, niter)` which returns the `p` dominant eigenvalues,
#    together with the associated eigenvectors,
#    of a symmetric matrix `B`.
#    To this end, implement the simultaneous iteration seen in class.
#    You can use the `myQR` function given below in order to calculate reduced QR decompositions.

# +
function myQR(B)
    Q, R =  qr(B)
    return Matrix(Q), R
end

function myEigen(B, p, niter)
    # Write your code here
end

@mark begin n = 2; A = randn(n, n); myEigen(A'A, n, 100)[1] ≈ reverse(eigen(A'A).values) end
@mark begin n = 4; A = randn(n, n); myEigen(A'A, n, 100)[1] ≈ reverse(eigen(A'A).values) end
@mark begin A = randn(5, 5); q, r = qr(A); B = q*Diagonal(1:5)*q'; myEigen(B, 3, 100)[1] ≈ [5; 4; 3] end
# -

# 2. Write a function `mySVD(B, p, niter)`
#    that returns the `p` dominant singular values of square matrix `B` (in a vector `σs`),
#    together with the associated left and right singular vectors (in matrices `Up` and `Vp`).
#    To this end,
#    notice that
#    $$
#    \mathsf A \mathsf A^\top = \mathsf U \mathsf \Sigma^2 \mathsf U^\top, \qquad
#    \mathsf A^\top \mathsf A = \mathsf V \mathsf \Sigma^2 \mathsf V^\top.
#    $$
#    Therefore, the left singular vectors of $\mathsf A$
#    are the eigenvectors of $\mathsf A \mathsf A^\top$,
#    while the right singular vectors of $\mathsf A$
#    are the eigenvectors of $\mathsf A^\top \mathsf A$,
#
#    **Hints**:
#
#    - In order to calculate singular vectors,
#      apply the function `myEigen` you wrote above to the matrices $\mathsf A \mathsf A^\top$
#      (in order to obtain left singular vectors),
#      and $\mathsf A^\top \mathsf A$ (in order to obtain right singular vectors).
#
#    - Once you have calculated the left and right singular vectors
#      associated with the `p` dominant singular values,
#      the singular values themselves can be obtained by extracting the diagonal from the matrix
#      $$
#      \Sigma_p = \mathsf U_p^\top \mathsf B \mathsf V_p.
#      $$
#      Here $\mathsf U_p$ and $\mathsf V_p$ are the matrices containing as columns the left and right singular vectors associated with the `p` dominant singular values,
#      respectively.

# +
function mySVD(B, p, niter)
    # Write your code here
end
# -

# 3. Write a function to test your implementation of `mySVD`.
#    Proceed as follows:
#    - Create a 10 x 10 matrix `B` with random entries using the `randn` function.
#    - Calculate its full singular value decomposition using `mySVD`.
#    - Check that `U` and `V` are orthogonal,
#      and that `U*Diagonal(σs)*V'` is close to `B`.

# +
# Write you code here
# -

# 4. The singular value decomposition is very useful for compressing matrices.
#    The idea of matrix compression based on SVD is the following:
#    given $p \leqslant n$, a matrix $B \in \mathbf R^{n\times n}$
#    can be approximated by
#    $$
#    \widetilde {\mathsf B} := \mathsf U_p \mathsf \Sigma_p \mathsf V_p^\top,
#    $$
#    where $\Sigma_p \in \mathbf R^{p \times p}$ is a diagonal matrix containing the $p$ dominant singular values of $\mathsf B$ on its diagonal,
#    and where $\mathsf U_p \in \mathbf R^{n \times p}$ and $\mathsf V_p \in \mathbf R^{n \times p}$ are rectangular matrices containing the associated left and right singular vectors,
#    respectively.
#
#    Since a grayscale image can be represented by a matrix containing the intensity values of all the pixels,
#    this approach for compressing matrices can be used for compressing grayscale images.
#    Use this method, i.e. calculate $\widetilde {\mathsf B}$, for $p \in \{5, 10, 20, 30\}$,
#    in order to compress the image `fabio_gray_256` given below,
#    and plot the compressed image for these values of `p`.

#    **Remarks**:
#    - (For information only) In practice, instead of storing the full matrix $\widetilde {\mathsf B}$, which contains $n^2$ entries,
#      we can store only the matrices $\mathsf U_p$, $\mathsf \Sigma_p$ and $\mathsf V_p$,
#      which together contain only $2np + p^2$ entries.
#      If $p \ll n$,
#      then the memory required to store these matrices is much smaller than the memory required to store the initial matrix $\mathsf B$.
#
#    - A function for drawing images based on the matrix of pixel intensity values is provided below.

# +
# Do not change the code in this cell
A = testimage("fabio_gray_256")

# Convert image to matrix of Float64
M = Float64.(A)

# Function to plot a grayscale image from the matrix
# containing the intensity values of all the pixels
function plot_matrix(B)
    plot(Gray.(B), ticks=false, showaxis=false)
end

plot_matrix(M)
# -

# +
p = 5
# Write you code here
# -

# +
p = 10
# Write you code here
# -

# +
p = 20
# Write you code here
# -

# +
p = 30
# Write you code here
# -

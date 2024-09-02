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

# ### 1. Exercise on Interpolation
#
# We aim to implement a numerical method to approximately solve the Euler-Bernoulli beam equation with homogeneous Dirichlet boundary conditions:
#
# $$
#   u\in C^4([0,1]),\quad\left\{\begin{aligned}  u''''(x) &= \varphi(x) \qquad \forall\, x\in(0,1),\\
#   u(0) &= u'(0) = u'(1) = u(1) = 0, \end{aligned}\right.
# $$
# where $\varphi(x) = (2\pi)^4\cos(2\pi x)$.
# In order to solve the equation approximately, we approximate the right-hand side $\varphi$ by its interpolating polynomial $\widehat \varphi$, and then we solve the equation exactly with the right-hand side $\widehat \varphi$ instead of $\varphi$.
#
# 1. Write a function `fit_values_and_slopes(u₀, up₀, u₁, up₁)` which returns the unique polynomial $p$ of degree 3 such that
#    $$
#    p(0) = u_0, \qquad p'(0) = up_0, \qquad p(1) = u_1, \qquad p'(1) = up_1.
#    $$
# 1. Write a function `approx(n)` implementing the approach described above for solving the PDE. The function should return a polynomial approximation of the solution based on an interpolation of **degree** $n$ of the right-hand side at equidistant points between 0 and 1, inclusive.
#
# 1. Write a function `estimate_error(n)` that approximates the error,
#    in $L^\infty$ norm,
#    between the exact and approximate solutions.
#    Note that the exact solution is given by
#    $$
#       \varphi(x) = \cos(2\pi x) - 1.
#    $$
#
# 1. Plot the error for $n$ in the range $\llbracket 5, 50 \rrbracket$. Use a logarithmic scale for the $y$ axis.
#
# *Hints:*
# - The Wikipedia page on *cubic Hermite splines* may be useful for the first exercise.
# - You can use the `fit` function from the `Polynomials.jl` library to obtain the interpolating polynomial:
#
#     ```julia
#     p = fit(x, y)
#     ```
#
#     where `x` are the interpolation nodes, and `y` are the values of the function to interpolate.
#
# - To calculate the analytical solution with a polynomial right-hand side, notice that all solutions are polynomials, and without boundary conditions, the solution is unique modulo a polynomial of degree at most 3.
#
# - You can use the `integrate` function from the `Polynomials.jl` library, which calculates an antiderivative of a polynomial:
#
#     ```julia
#     P = integrate(p)
#     ```
#
# - Use the `BigFloat` format to limit rounding errors.

# +
using LaTeXStrings
using Polynomials
using Plots

function fit_values_and_slopes(u₀, up₀, u₁, up₁)
    # We look for polynomials p(x) = a₀ + a₁ x + a₂ x² + a₃ x³
    A = [1 0 0 0; 0 1 0 0; 1 1 1 1; 0 1 2 3]
    α = A\[u₀; up₀; u₁; up₁]
    return Polynomial(α)
end

# Test our code
p = fit_values_and_slopes(-1, -1, 1, 1)
plot(p, xlims=(0, 1))
# -

# +
φ(x) = (2π)^4 * cospi(2*x)
u(x) = cospi(2*x) - 1

function approx(n)
    X = LinRange{BigFloat}(0, 1, n + 1)
    Y = φ.(X)
    p = fit(X, Y)
    uh = integrate(integrate(integrate(integrate(p))))
    ∂uh = derivative(uh)
    uh -= fit_values_and_slopes(uh(0), ∂uh(0), uh(1), ∂uh(1))
    return uh
end

plot(approx(3), xlims=(0, 1), label="Exact solution")
plot!(u, xlims=(0, 1), label="Approximate solution")
# -


# +
function estimate_error(n)
    un = approx(n)
    x_fine = LinRange(0, 1, 1000)
    un_fine, u_fine = un.(x_fine), u.(x_fine)
    return maximum(abs.(u_fine - un_fine))
end
# -

# +
ns = 5:50
errors = estimate_error.(ns)
plot(ns, errors, marker = :circle, label=L"$L^{\infty}$ Error")
plot!(yaxis=:log, lw=2)
# -

# ### 2. Exercise on numerical integration
#
# Our goal in this exercise is to write a program in order to calculate integrals of the form
# $$
# I[f] := \int_0^{\infty} f(x) e^{-x} \, d x
# $$
# To this end, we will use Laguerre polynomials,
# which are orthogonal polynomials for the following inner product:
#    $$
#     \langle f, g \rangle := \int_0^{\infty} f(x) g(x) e^{-x} \, d x
#    $$
# 1. Using that Laguerre polynomials satisfy the recurrence relation
#    $$
#        L_{k + 1}(x) = \frac{(2k + 1 - x)L_k(x) - k L_{k - 1}(x)}{k + 1}, \qquad L_0(x) = 1, \qquad L_1(x) = 1-x,
#    $$
#    write a function `laguerre(n)` which returns the Laguerre polynomial of degree $n$.
#
# 1. Write a function `get_nodes_and_weights(n)` which returns the nodes and weights of the Gauss-Laguerre quadrature with $n$ nodes. In order to construct Lagrange polynomials, you may find it useful to use the `fromroots` function from `Polynomials.jl`.
#    You will need to use that
#    $$
#    \int_0^{\infty} x^n e^{-x} \, dx = n!
#    $$
#
# 1. Write a function `integrate_laguerre(f, n)`, which returns an approximation of $I[f]$ obtained by Gauss-Laguerre integration with $n$ nodes.
#
# 1. Set $n = 3$, and calculate numerically the degree of precision of the corresponding quadrature rule.
#
# 1. Set $f(x) = \sin(x)$, and plot the integration error as a function of $n$.

# +
function laguerre(n)
    if n == 0
        return Polynomial([1])
    elseif n == 1
        return Polynomial([1, -1])
    else
        k = n-1
        x = Polynomial([0, 1])
        return ((2k + 1 - x) * laguerre(k) - k*laguerre(k-1))/(k+1)
    end
end
# -

# +
function get_nodes_and_weights(n)
    nodes = roots(laguerre(n))
    weights = zero(nodes)
    for i in 1:n
        ℓ = fromroots(nodes[1:end .!= i])
        ℓ /= ℓ(nodes[i])
        weights[i] = factorial.(0:n-1)'ℓ.coeffs
    end
    return nodes, weights
end
# -

# +
function integrate_laguerre(f, n)
    nodes, weights = get_nodes_and_weights(n)
    return f.(nodes)'weights
end
# -

# +
n = 5
for i in 1:15
    correct = integrate_laguerre(x -> x^i, n) ≈ factorial(i)
    println("f = x^$i, Rule exact? ", correct)
    if !correct
        println("Degree of precision = ", i - 1)
        break
    end
end
# -

# +
ns = 1:20
f(x) = sin(x)
In = 1/2
Ih = integrate_laguerre.(f, ns)
plot(ns, abs.(Ih .- In), yscale=:log10, xlabel=L"n", ylabel="Error")
scatter!(ns, abs.(Ih .- In))
# -

# ### 3. Exercise on linear systems
#
# We consider in this exercise Poisson's equation in the domain $\Omega = (0, 2) \times (0, 1)$,
# equipped with homogeneous Dirichlet boundary conditions:
# $$
#     \begin{aligned}
#         - \Delta f(x, y) &= b(x, y), \qquad x \in \Omega, \\
#         f(x) &= 0, \quad \qquad x \in \partial \Omega.
#     \end{aligned}
# $$
# The right-hand side is
# $$
#     b(x, y) = \sin(4\pi x) + \sin(2\pi y).
# $$
# A number of methods can be employed in order to discretize this partial differential equation.
# After discretization, a finite-dimensional linear system of the form $\mathsf A \mathbf x = \mathbf b$ is obtained.
# A Julia function for calculating the matrix $\mathsf A$ and the vector $\mathbf b$ using the finite difference method is given to you below,
# as well as a function to plot the solution.
# The goal of this exercise is to solve the linear system using the conjugate gradient method.
#
# 1. Based on the conjugate gradient algorithm given in the lecture notes,
#    write a function `conjugate_gradient(A, b)` to solve the linear system $\mathsf A \mathbf x = \mathbf b$.
#    Use a vector of zeros for $\mathbf x_0$, and the stopping criterion
#    $$
#    \lVert \mathsf A \mathbf x_k - \mathbf b \rVert \leq \varepsilon \lVert \mathbf b \rVert.
#    $$
#
# 1. Make a graph illustrating the evolution of the residual $\lVert \mathbf r_k \rVert$, the error $\lVert \mathbf x_k - \mathbf x_* \rVert$,
#      and the error $\lVert \mathbf x_k - \mathbf x_* \rVert_{\mathsf A}$ as a function of $k$.
#      Choose appropriate scales for the graph (linear-linear, log-log, linear-log, or log-linear).

# +
using LinearAlgebra
using Plots
import SparseArrays

# Domain size
Lx, Ly = 2, 1

# Number of discretization points along x and y, including the boundary points
nx, ny = 101, 101

function discretize(nx, ny)
    hx, hy = Lx/(nx - 1), Ly/(ny - 1)
    Dxx = (1/hx^2) * SymTridiagonal(2ones(nx-2), -ones(nx-3))
    Dyy = (1/hy^2) * SymTridiagonal(2ones(ny-2), -ones(ny-3))
    A = kron(Dxx, I(ny-2)) + kron(I(nx-2), Dyy)
    xgrid = Lx/(nx-1) * (1:nx-2)
    ygrid = Ly/(ny-1) * (1:ny-2)
    x_2d = reshape([x for y in ygrid, x in xgrid], (nx-2)*(ny-2))
    y_2d = reshape([y for y in ygrid, x in xgrid], (nx-2)*(ny-2))
    b = sin.(4π*x_2d) + sin.(2π*y_2d)
    return SparseArrays.SparseMatrixCSC(A), b
end

function plot_solution(f)
    f = reshape(f, ny-2, nx-2)
    z = [zeros(nx)'; zeros(ny-2) f zeros(ny-2); zeros(nx)']  # Add boundary
    xgrid = Lx/(nx-1) * (0:nx-1)
    ygrid = Ly/(ny-1) * (0:ny-1)
    heatmap(xgrid, ygrid, z, c=:viridis, levels=50)
end

# +
function conjugate_gradients(A, b)
    n = length(b)
    x = ones(n)
    d = r = A*x - b
    ε = 1e-10
    xs = [x]
    while √(r'r) ≥ ε * √(b'b)
        ω = d'r / (d'A*d)
        x -= ω*d
        r = A*x - b
        β = r'A*d/(d'A*d)
        d = r - β*d
        push!(xs, x)
    end
    xs
end

# +
A, b = discretize(nx, ny)
x_ = A\b
xs = conjugate_gradients(A, b)
plot_solution(xs[end])
# -

# +
rs = [A*x - b for x in xs]
es = [x - x_ for x in xs]
plot(xlabel=L"k", yaxis=:log)
plot!([√(r'r) for r in rs], label=L"\Vert Ax - b \Vert")
plot!([√(e'r) for (e, r) in zip(es, rs)], label=L"\Vert x - x_{*} \Vert_A")
plot!([√(e'e) for e in es], label=L"\Vert x - x_{*} \Vert")
# -

# ### 4. Exercise on nonlinear equations
#
# Data from a biological experiment are collated in the cell below.
# We wish to fit to this data a model of the form
# $$
#     R = \frac{\beta_1 C}{\beta_2 + C},
# $$
# where $R$ is the rate, $C$ is the concentration,
# and $\beta_1$ and $\beta_2$ are parameters to be determined.
# To this end, we shall minimize the sum of square residuals:
# $$
#     J(\mathbf \beta) := \sum_{i=1}^{N} \lvert r_i(\mathbf \beta) \rvert^2, \qquad r_i(\mathbf \beta) := R_i - \frac{\beta_1 C_i}{\beta_2 + C_i}
# $$
#
# 1. Find the optimal $\mathbf \beta$ using the chord method, applied to the nonlinear equation $\nabla J(\beta) = 0$.
#
# 1. Find the optimal $\mathbf \beta$ using the Newton-Raphson method.
#
# 1. In applications, it is often desirable solve nonlinear least squares problems of the type considered here without calculating second derivatives, i.e. without resorting to the Newton-Raphson method. The most famous algorithm to this end is the so-called Gauss-Newton method.
#
#    The idea of the Gauss-Newton method is the following:
#    given an approximation $\beta_k$ of the optimal vector of parameters $\beta$,
#    we approximate each residual $r_i(\beta)$ by its linearization around $\beta_k$,
#    which is given by
#    $$
#       \widetilde r_i(\beta) = r_i(\beta_k) + \nabla r_i(\beta_k)^\top (\beta - \beta_k).
#    $$
#    We then calculate $\beta_{k+1}$ as the minimizer of
#    $$
#    \widetilde J(\mathbf \beta) = \sum_{i=1}^{N} \lvert \widetilde r_i(\mathbf \beta) \rvert^2.
#    $$
#    This is now a **linear** least squares problem,
#    the solution of which is given by the normal equations
#    $$
#    \beta_{k+1} - \beta_{k} = (\mathsf A^\top \mathsf A)^{-1} \mathsf A^\top \mathbf b,
#    $$
#    where
#    $$
#    \mathsf A =
#    \begin{pmatrix}
#    \nabla r_1(\beta_k)^\top \\
#    \vdots \\
#    \nabla r_N(\beta_k)^\top
#    \end{pmatrix},
#    \qquad
#    \mathbf b =
#    -\begin{pmatrix}
#    r_1(\beta_k) \\
#    \vdots \\
#    r_N(\beta_k)
#    \end{pmatrix}
#    $$
#
# *Hints*:
#
# - You may use `ForwardDiff` in order to compute derivatives.

# +
concentrations = [0.038, 0.194, 0.425, 0.626, 1.253, 2.500, 3.740]
rates = [0.050, 0.127, 0.094, 0.2122, 0.2729, 0.2665, 0.3317]

# +
using ForwardDiff
r(i, β) = rates[i] - β[1] * concentrations[i] / (β[2] + concentrations[i])
r(β) = [r(i, β) for i in 1:length(rates)]
J(β) = r(β)'r(β)
∇J(β) = ForwardDiff.gradient(J, β)
HJ(β) = ForwardDiff.hessian(J, β)   

# Chord method
β = [0.3; 0.5]
α = 1000
maxiter = 100_000
for i in 1:maxiter
    β -= ∇J(β) / α
    if norm(∇J(β)) < 1e-12
        println("Converged in $i iterations!")
        break
    end
end
println(β, " ", norm(∇J(β)))

# Newton-Raphson
β = [0.3; 0.5]
α = 1000
maxiter = 100_000
for i in 1:maxiter
    β -= HJ(β)\∇J(β)
    if norm(∇J(β)) < 1e-12
        println("Converged in $i iterations!")
        break
    end
end
println(β, " ", norm(∇J(β)))

# Gauss-Newton
β = [0.3; 0.5]
α = 1000
maxiter = 100_000
for i in 1:maxiter
    A = ForwardDiff.jacobian(r, β)
    β -= (A'A)\(A'r(β))
    if norm(∇J(β)) < 1e-6
        println("Converged in $i iterations!")
        break
    end
end
println(β)
# -

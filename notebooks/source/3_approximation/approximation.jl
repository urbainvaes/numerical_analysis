# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.1
#   kernelspec:
#     display_name: Julia 1.10.4
#     language: julia
#     name: julia-1.10
# ---

# ## Notebook 3: Approximation

# + nbgrader={"grade": false, "grade_id": "cell-4d638a879ba6f86e", "locked": true, "schema_version": 3, "solution": false, "task": false}
using LaTeXStrings
using LinearAlgebra
using Plots
using Polynomials

# + [markdown] nbgrader={"grade": false, "grade_id": "cell-b493a5b381036ad7", "locked": true, "schema_version": 3, "solution": false, "task": false}
# ### <font color='orange'>[Exercise 1]</font> Least squares approximation

# + [markdown] nbgrader={"grade": false, "grade_id": "cell-9024d5e738948ab8", "locked": true, "schema_version": 3, "solution": false, "task": false}
# The objective of this exercise is to approximate some given data `x`, `y` using a
# polynomial of a given degree,
# potentially much lower than the number of data points.
#
# 1. Without using any library,
#    write a function `polyfit(x, y, d)` which,
#    given data points
#    $
#        (x_1, y_1), \dotsc, (x_N, y_N)
#    $
#    and an integer $0 \leq d \leq N-1$,
#    calculates the unique polynomial $p \in \mathbb R[x]$ of degree at most $d$ minimizing the total error
#
#    $$
#        E := \frac{1}{2} \sum_{n=0}^{N} \bigl\lvert p(x_n) - y_n \bigr\rvert^2.
#    $$
#    Your function should return a vector containing the $d+1$ coefficients of $p$,
#    starting from the constant term.
#    <details>
#        <summary>
#            <em><font color='gray'>Hint (click to display)</font></em>
#        </summary>
#
#    Within the function, you can proceed as follows.
#    First, create the following matrix and vector:
#    $$
#        \mathbf{A} =
#        \begin{pmatrix}
#            1 & x_0 & \dots & x_0^d \\
#            \vdots & \vdots & & \vdots \\
#            1 & x_{N} & \dots & x_N^d
#        \end{pmatrix},
#        \qquad
#        \mathbf{b} =
#        \begin{pmatrix}
#            y_0 \\
#            \vdots \\
#            y_N
#        \end{pmatrix}.
#    $$
#    Note that the error $E$ rewrites as follows:
#    $$
#    E(\boldsymbol \alpha) = \frac{1}{2} \bigl\lVert \mathsf A \boldsymbol \alpha - \boldsymbol b \bigr\rVert^2,
#    $$
#    where $\boldsymbol \alpha$ is a vector containing the coefficients
#    of the polynomial, in order of increasing degree, and $\lVert \cdot \rVert$
#    is the Euclidean norm.
#    In other words, the function $E$ is a quadratic function of the
#    vector of coefficients of the polynomial.
#    Writing $\nabla E(\boldsymbol \alpha) = 0$ leads to the so-called **normal equations**:
#    $$
#        \mathsf{A}^\top \mathsf{A} \boldsymbol{\alpha} = \mathsf{A}^\top \mathbf{b}.
#    $$
#    This is a linear system with a square invertible matrix on the left-hand side,
#    which can be solved using the backslash operator `\`;
#    in fact you can write just `A\b` instead of `(A'A)\(A'b)`,
#    because the operator `\` performs a least squares approximation by default.
#    </details>

# + nbgrader={"grade": false, "grade_id": "cell-16c5d2765a6be02d", "locked": false, "schema_version": 3, "solution": true, "task": false}
function polyfit(x, y, d)
    # Your code here
end

# + nbgrader={"grade": true, "grade_id": "cell-ba35e3fde06b7198", "locked": true, "points": 1, "schema_version": 3, "solution": false, "task": false}
n = 10 ; x = 1:n ; y = randn(n)
@assert polyfit([0.], [0.], 0) ≈ [0.]
@assert polyfit(1:5, 1:5, 1) ≈ [0., 1.]
@assert polyfit(x, y, 0) ≈ [sum(y)/n]
@assert polyfit(x, y, 0) ≈ [sum(y)/n]
@assert polyfit(x, y, 2) ≈ fit(x, y, 2).coeffs

# + [markdown] nbgrader={"grade": false, "grade_id": "cell-abfabcd8fa09fad4", "locked": true, "schema_version": 3, "solution": false, "task": false}
# 2. Write also a function `polyval(α, X)`
#    to evaluate the polynomial
#    $$
#        p(x) = \alpha_0 + \alpha_1 x + \dotsc + \alpha_d x^d
#    $$
#    at all the points in `X`.
#    The function should return the result in a vector.

# + nbgrader={"grade": false, "grade_id": "cell-ad0df3d23dd2927a", "locked": false, "schema_version": 3, "solution": true, "task": false}
function polyval(α, X)
    # Your code here
end

# + nbgrader={"grade": true, "grade_id": "cell-e4b2a8fabeacbd21", "locked": true, "points": 1, "schema_version": 3, "solution": false, "task": false}
n = 10 ; α = randn(n)
@assert polyval([0.], [1., 2., 3.]) == [0., 0., 0.]
@assert polyval(α, [0., 1.]) ≈ [α[1], sum(α)]

# + [markdown] nbgrader={"grade": false, "grade_id": "cell-107ce28d10b5c285", "locked": true, "schema_version": 3, "solution": false, "task": false}
# Use the data given below,
# of the altitude of a marble in free fall as a function of time on a remote planet, to test your code.
# The experiment was performed on a different planet. Can you find which one? See [Gravitational acceleration](https://en.wikipedia.org/wiki/Gravitational_acceleration).

# +
# Time since dropping the marble
x = [0., 1., 2., 3., 4., 5.]

# Altitude of the marble
y = [100., 98.26, 93.56, 81.79, 71.25, 53.22]

# Fit using polyfit
α = polyfit(x, y, 2)

# Evalutate at X
X = LinRange(0, 5, 200)
Y = polyval(α, X)

plot(X, Y, label="My approximation")
scatter!(x, y, label="Data")

# + [markdown] nbgrader={"grade": false, "grade_id": "cell-32309e8ba988d9c6", "locked": true, "schema_version": 3, "solution": false, "task": false}
#    <details>
#        <summary>
#            <em><font color='gray'>Least squares approximation using `Polynomials.fit` (click to display)</font></em>
#        </summary>
#
#    We learnt in the previous notebook that `Polynomials.fit` could be employed for
#    polynomial interpolation. This function can also be used for fitting
#    data by minimizing the sum of squared residuals,
#    which can be achieved as follows:
#    </details>

# + nbgrader={"grade": false, "grade_id": "cell-aa104b090b54806d", "locked": true, "schema_version": 3, "solution": false, "task": false}
# This returns structure of type `Polynomial`, with associated degree 2
p = fit(x, y, 2)
@show p

# The coefficients can be obtained as follows
@show p.coeffs

# The structure `p` behaves like a function
@show p(0)

X = LinRange(0, 5, 200)
plot(X, p.(X), label="Polynomials.jl approximation")

# + [markdown] nbgrader={"grade": false, "grade_id": "cell-b493a5b381036ad7", "locked": true, "schema_version": 3, "solution": false, "task": false}
# ### <font color='orange'>[Exercise 2]</font> Least squares approximation for data analysis

# The dataset loaded through the following Julia commands contains data on the vapor pressure of mercury as a function of the temperature.

# +
import RDatasets
data = RDatasets.dataset("datasets", "pressure")
# -

# Using `Polynomials.jl`,
# find a low-dimensional mathematical model of the form
# \begin{equation}
#     p(T) = \exp \bigl( \alpha_0 + \alpha_1 T + \alpha_2 T^2 + \alpha_3 T^3 \bigr)
# \end{equation}
# for the pressure as a function of the temperature.
# Plot the approximation together with the data.

# +
# -

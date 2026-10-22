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
#     display_name: Julia 1.10.5
#     language: julia
#     name: julia-1.10
# ---

# # NYU Paris - Numerical Analysis

# ## Final exam 
#
# All three questions are graded out of 10. The total grade is out of 40, with only your best 4 questions counting.

# +
using ForwardDiff
using LinearAlgebra
using Polynomials
using Plots

macro mark(bool_expr)
    return :( print($bool_expr ? "✔️" : "❌"))
end
# -

# ### 1. Floating point arithmetic
#
# 1. <span style="color: gray;">T/F:</span> Let $(\cdot)_2$ denote base 2 representation. It holds that  
#    $$
#    3 \times (0.0101)_2 = (0.1111)_2.
#    $$
# 
# 1. <span style="color: gray;">T/F:</span>  Does the following equality hold? Explain your reasoning.  
#    $$
#    (0.\overline{011})_2 = \frac{3}{4}.
#    $$
# 
# 1. <span style="color: gray;">T/F:</span>  In Julia, `Float64(x) == Float32(x)` is `true` if `x` is a rational number.
# 
# 1. <span style="color: gray;">T/F:</span>  The value of the machine epsilon for the double precision format is the same in Julia and Python.
# 
# 1. <span style="color: gray;">T/F:</span>  The spacing (in absolute value) between successive double-precision (`Float64`) floating point numbers is equal to the machine epsilon.
# 
# 1. <span style="color: gray;">T/F:</span>  All the natural numbers can be represented exactly in the double precision floating point format `Float64`.
# 
# 1. <span style="color: gray;">T/F:</span>  Machine addition in the `Float64` format is associative but not commutative.
# 
# 1. <span style="color: gray;">T/F:</span>  In Julia, let `f(x) = (x == x/100.0) ? x : f(x/100.0)`. Then `f(a)` returns `0.0` for all finite numbers `a` representable in the `Float64` format.
# 
# 1. In Julia, `exp(eps()) == 1 + eps()` evaluates to `true`. Explain briefly why.
# 
# 1. In Julia, `sqrt(1 + eps()) == 1 + eps()` evaluates to `false`. Explain briefly why.

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

# ### 3. Exercise on nonlinear equations
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

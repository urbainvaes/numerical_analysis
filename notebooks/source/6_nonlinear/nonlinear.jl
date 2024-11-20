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

# ## Notebook 6 : Nonlinear equations

using LinearAlgebra
using Plots
using Polynomials

# ### <font color='orange'>[Exercise 1]</font> Newton-Raphson in dimension 2
#
# We consider the following linear system:
# $$
# \left \{
#     \begin{aligned}
#         &y = (x-1)^2 \\
#         &x^2 + y^2 = 4
#     \end{aligned}
# \right.
# $$
#
# 1. Plot appropriate graphs to roughly visualize the zone(s) containing solution(s).
#

# + nbgrader={"grade": true, "grade_id": "cell-87d5b9021583a6dc", "locked": false, "points": 1, "schema_version": 3, "solution": true, "task": false}
### BEGIN SOLUTION
plot(x->(x-1)^2, xlim=(-1,2), xlabel="x", ylabel="y", label="(x-1)²")
plot!(x->√(4-x^2), label="√(4-x²)", aspect_ratio=:equal)
### END SOLUTION
# -

# 2. Implement the chord method to compute precise approximations the solutions.

# + nbgrader={"grade": true, "grade_id": "cell-210b54914b0fd6c7", "locked": false, "points": 0, "schema_version": 3, "solution": true, "task": false}
### BEGIN SOLUTION
nothing
### END SOLUTION
# -

# 3. Implement the Newton-Raphson method to compute precise approximations the solutions,
#    and print the number of iterations that are required.

# + nbgrader={"grade": true, "grade_id": "cell-45afdb05a9da8e1d", "locked": false, "points": 0, "schema_version": 3, "solution": true, "task": false}
# ### BEGIN SOLUTION
nothing
# ### END SOLUTION
# -

# 4. Estimate the order of convergence, i.e., $q$ such that
#    $$
#    \lim_{k \to \infty} \frac{\lVert \mathbf{x}_{k+1} - \mathbf{x}_* \rVert}{\lVert \mathbf{x}_k - \mathbf{x}_* \rVert^q} \in \mathbb{R}^*_+
#    $$
#    for the solutions.
#
#    <details>
#        <summary>
#            <em><font color='gray'>Hint (click to display)</font></em>
#        </summary>
#
#    Let $y_k := -\log(\lVert \mathbf{x}_k - \mathbf{x}_* \rVert^q)$.
#    The limit given implies that
#    $$
#    \lim_{k \to \infty} y_{k+1} - q y_k = C \in \mathbb{R}.
#    $$
#    From this, we deduce that
#    $$
#    q = \lim_{k \to \infty} \frac{y_{k+1}}{y_k}.
#    $$
#    This equation allows estimating $q$ from $y_{k+1}$ and $y_k$ for sufficiently large $k$.

# + nbgrader={"grade": true, "grade_id": "cell-aaa07e2866c41119", "locked": false, "points": 1, "schema_version": 3, "solution": true, "task": false}
### BEGIN SOLUTION
nothing
### END SOLUTION
# -

# ### <font color='green'> Introduction to automatic differentiation</font>

# ### <font color='orange'>[Exercise 2]</font> Calculation of the square root using the Babylonian method
#
# Let a real parameter $a > 0$ and the sequence defined by
# <a id="baby"></a>
# $$
# \tag{1}
# x_0 > 0 \qquad ; \qquad ∀k ∈ \mathbb{N}, \quad x_{k+1} = \frac{1}{2}\left(x_k + \frac{a}{x_k}\right)
# $$
#
#   > *Preliminary questions (to be done on scratch paper but not required for the submission)*
#   >
#   > i) By writing $x_{k+1} - \sqrt{a}$ as a function of $x_k - \sqrt{a}$ and then $x_{k+1} - x_k$, show that $(x_k)$ converges quadratically to $x_* = \sqrt{a}$ for any $x_0 > 0$.
#   >
#   >    <details>
#   >        <summary>
#   >            <em><font color='gray'>Help (click to reveal)</font></em>
#   >        </summary>
#   >
#   >    - First note that if $x_0 > 0$, then $x_k > 0$ for all $k$.
#   >    - Show that $x_{k+1} - \sqrt{a} = \frac{(x_k - \sqrt{a})^2}{2 x_k}$ and that $x_{k+1} - x_k = \frac{a - x_k^2}{2 x_k}$.
#   >    - Deduce that $(x_k)_{k \geq 1}$ is bounded below by $\sqrt{a}$ and is decreasing (be careful to only consider the reasoning for $k \geq 1$), so it converges.
#   >    - Conclude that the limit is necessarily $\sqrt{a}$ and that the convergence is quadratic.
#   >    </details>
#   >
#   > ii) Show that the recurrence formulation <a href="#baby">(1)</a> is nothing but the Newton-Raphson algorithm applied to a function to be identified that vanishes at $x_* = \sqrt{a}$.
#
# 1. Construct a function `Babylonian` that takes `a` and an integer `n` (defaulting to `10`) as arguments and returns the vector $[x_0, x_1, \ldots, x_n]$, initializing the sequence with $x_0 = \frac{1 + a}{2}$.

# + nbgrader={"grade": false, "grade_id": "cell-51aeec58821a4395", "locked": false, "schema_version": 3, "solution": true, "task": false}
function Babylonian(a; n = 10)
    ### BEGIN SOLUTION
    x = [(1+a)/2]
    for i = 1:n push!(x, (x[end]+a/x[end])/2) end
    return x
    ### END SOLUTION
end

# + nbgrader={"grade": true, "grade_id": "cell-77f4b908283dc59e", "locked": true, "points": 1, "schema_version": 3, "solution": false, "task": false}
for a in (0.1, 2, 25, 100)
    @assert Babylonian(a)[end] ≈ √a
end
# -

# 2. Plot the error $|x_k - x_*|$ as a function of the index $k$ for $a = 2$.

# + nbgrader={"grade": true, "grade_id": "cell-928ca89fac43223b", "locked": false, "points": 1, "schema_version": 3, "solution": true, "task": false}
### BEGIN SOLUTION
plot(abs.(Babylonian(2) .- √2), yaxis=:log10, xlabel="k", ylabel="|xₖ-x*|", label="")
### END SOLUTION
# -

# The idea behind the remainder of the exercise is to apply the `Babylonian` function defined earlier to an argument `a` not of type `Float64` but of a new type that allows us to estimate both the value of $\sqrt{a}$ and the derivative of $a \mapsto \sqrt{a}$, which is $\frac{1}{2\sqrt{a}}$. For this, we introduce new numbers called **dual numbers**. These are defined similarly to complex numbers, based on the definition of a special number denoted $\varepsilon$, such that a dual number is written as $x = a + b\varepsilon$, where $a$ and $b$ are real numbers. In a sense, $\varepsilon$ plays a role analogous to the complex $i$, with the difference that we set $\varepsilon^2 = 0$. The purpose of such numbers is to be able to store both the value of a function and its derivative by writing:
#
# <a id="fdual"></a>
# $$
# \tag{2}
# f(a + b\varepsilon) = f(a) + f'(a) b \varepsilon
# $$
#
# This means that the derivative of $f$ at $a$ can be obtained by extracting the component on $\varepsilon$ of $f(a + \varepsilon)$ (i.e., by setting $b = 1$).
#
# In practice, it is necessary to redefine the behavior of common functions in accordance with <a href="#fdual">(2)</a>. However, in the current application, only the operations `+`, `-`, `*`, and `/` will be needed and must be overloaded to allow dual numbers as arguments. Additionally, it will be necessary to implement the `convert` function to convert a real number to a dual number and the `promote_rule` to express that in the presence of an operation involving two numbers, one of which is dual, both must first be expressed as dual numbers before the operation is performed. Note that operator and function overloading is only possible if they are explicitly imported using, for example, `import Base: +, -, ...`. It is also possible to define the `Base.show` function so that a dual number is displayed in the explicit form `a + bɛ`.
#
# The overloading of operators is mathematically expressed as follows:
#
# $$
# \begin{align*}
# (a + b\varepsilon) + (c + d\varepsilon) &= (a + c) + (b + d)\varepsilon \\
# (a + b\varepsilon) - (c + d\varepsilon) &= (a - c) + (b - d)\varepsilon \\
# (a + b\varepsilon) \cdot (c + d\varepsilon) &= ac + (bc + ad)\varepsilon \\
# \frac{(a + b\varepsilon)}{(c + d\varepsilon)} &= \frac{a}{c} + \frac{bc - ad}{c^2} \varepsilon
# \end{align*}
# $$
#
# Alternatively, the last operation can be defined as $\mathrm{inv}(a + b\varepsilon) = \frac{1}{a} - \frac{b}{a^2} \varepsilon$, and then `u/v = u * inv(v)`.
#
# 3. Study the `struct D` defined below to represent a dual number, along with the associated lines of code. Complete the missing parts of the code, specifically the implementations of `/` and `inv`.

# + nbgrader={"grade": false, "grade_id": "cell-6a89a36f52efd269", "locked": false, "schema_version": 3, "solution": true, "task": false}
import Base: +, -, *, /, inv, isapprox, convert, promote_rule
using LinearAlgebra

struct D <: Number
    f::Tuple{Float64, Float64}
end
D(a::Real, b::Real) = D((a, b))
+(x::D, y::D) = D(x.f .+ y.f)
-(x::D, y::D) = D(x.f .- y.f)
*(x::D, y::D) = D(x.f[1]*y.f[1], x.f[2]*y.f[1] + x.f[1]*y.f[2])
### BEGIN SOLUTION
/(x::D, y::D) = D(x.f[1]/y.f[1], (y.f[1]*x.f[2] - x.f[1]*y.f[2])/y.f[1]^2)
inv(x::D) = D(1/x.f[1], -x.f[2]/x.f[1]^2)
### END SOLUTION
-(x::D) = D(.-(x.f))
isapprox(x::D, y::D; kwargs...) = all(isapprox.(x.f, y.f ; kwargs...))
convert(::Type{D}, x::Real) = D((x,zero(x)))
promote_rule(::Type{D}, ::Type{<:Real}) = D
Base.show(io::IO,x::D) = print(io,x.f[1],x.f[2]<0 ? " - " : " + ",abs(x.f[2])," ε")

# Construction of a dual number
x = D(0.1, -1.6)
# -

# 4. Define an instance of the number `ɛ` (use `\varepsilon` and press TAB to display ε), in other words the number `0 + 1ɛ`, and perform some operations to verify the implementations (use the `@show` macro to display a dual number), for example:
#
#    ```julia
#    @show (1+2ɛ)*(3+4ɛ)
#    @show 1/(1+ɛ)
#    @show (1+2ɛ)/(2-ɛ)
#    ```

# + nbgrader={"grade": true, "grade_id": "cell-ec62ba59865b5495", "locked": false, "points": 1, "schema_version": 3, "solution": true, "task": false}
### BEGIN SOLUTION
ε = D((0,1))
@show (1+2ɛ)*(3+4ɛ)
@show 1/(1+ɛ)
@show (1+2ɛ)/(2-ɛ)
### END SOLUTION

# + nbgrader={"grade": true, "grade_id": "cell-c9da9f0b525d17d9", "locked": true, "points": 1, "schema_version": 3, "solution": false, "task": false}
@assert (1+2ɛ)*(3+4ɛ) == 3+10ɛ "error"
@assert 1/(1+ɛ) == 1-ɛ "error"
@assert (1+2ɛ)/(2-ɛ) == 1/2+5ɛ/4 "error"

### BEGIN HIDDEN TESTS
a, b, c, d = rand(4)
@assert 1/(a+b*ɛ) == inv(a+b*ɛ) == 1/a-b/a^2*ε
@assert (a+b*ɛ)/(c+d*ɛ) == a/c + (b*c-a*d)/c^2*ε
### END HIDDEN TESTS
# -

# 5. Use the dual number structure to estimate the derivative of the square root function from the Babylonian method (by directly using the `Babylonian` function without rewriting it).

# + nbgrader={"grade": false, "grade_id": "cell-b0b7e3c217eea6a2", "locked": false, "schema_version": 3, "solution": true, "task": false}
function derivative_sqrt(a; n = 10)
    ### BEGIN SOLUTION
    return Babylonian(a+ε)[end].f[2]
    ### END SOLUTION
end

# + nbgrader={"grade": true, "grade_id": "cell-f3a9a91b390a3270", "locked": true, "points": 1, "schema_version": 3, "solution": false, "task": false}
for a in (0.1, 2, 25, 100)
    @assert derivative_sqrt(a) ≈ 1/2√a
end
# -

# 6. Overlay on a graph the derivative of the square root obtained by the Babylonian method using dual numbers and the analytical expression.

# + nbgrader={"grade": true, "grade_id": "cell-b10d808d11cb751d", "locked": false, "points": 1, "schema_version": 3, "solution": true, "task": false}
### BEGIN SOLUTION
xplot = LinRange(0.1,10,200)
plot(xplot, x -> Babylonian(x+ε)[end].f[2], label="Méth. babylonienne")
plot!(xplot, x -> 1/2√x, linestyle=:dashdot, linewidth=3, label="1/2√x")
### END SOLUTION
# -

# 7. Propose an analogous method to calculate the $p^\textrm{th}$ root of a number $a$, i.e., $\sqrt[p]{a}$. Verify that the derivative of the $p^\textrm{th}$ root can also be obtained using dual numbers without any additional lines of code.

# + nbgrader={"grade": false, "grade_id": "cell-fe8183299ca40337", "locked": false, "schema_version": 3, "solution": true, "task": false}
function nthrt(a, p=2; x=1, n=100)
    ### BEGIN SOLUTION
    for i = 1:n x = ((p-1)*x+a/x^(p-1))/p end
    return x
    ### END SOLUTION
end

# + nbgrader={"grade": true, "grade_id": "cell-611b2e70df58e157", "locked": true, "points": 1, "schema_version": 3, "solution": false, "task": false}
for a in (0.1, 2, 25, 100), p in (2, 3, 5)
    @assert nthrt(a+ε, p) ≈ a^(1/p) + a^(1/p-1)/p*ε "error for (a,p)=($a,$p)"
end
# -

# ### <font color='orange'>[Exercise 3]</font> Application to optimization
#
# In Julia, automatic differentiation based on dual numbers is implemented in the `ForwardDiff` library.
# For example, the following code snippet
#    ```julia
#        using ForwardDiff
#        f(t) = 9.81 * t^2/2
#        ForwardDiff.derivative(f, 1)
#    ```
# returns the value of $f'(1)$.
# Unlike our simple implementation of dual numbers above,
# `ForwardDiff` is able to calculate derivatives of higher orders.
#    For example,
#    ```julia
#        using ForwardDiff
#        f(x) = x -> x[1]^2 + 3x[2]^2
#        ForwardDiff.hessian(f, [1, 2])
#    ```
# returns the Hessian of the function $f(x, y) = x^2 + 3y^2$ at $(1, 2)$.
# Using the `ForwardDiff` library, solve the following exercises.
#
#
# 1. We have $n$ points $(x_i, y_i)$ of an unknown function $y = f(x)$. We want to approximate the data with a function of the form
#    $$
#    \widetilde f(x) = \frac{a}{b + x}
#    $$
#    by minimizing
#    $$
#    E(a,b) := \sum_{i=1}^{n} |\widetilde f(x_i) - y_i|^2.
#    $$
#    To find the minimizer, we propose to solve the nonlinear system of equations given by
#    $$
#    \begin{aligned}
#    \partial_a E(a, b) &= 0, \\
#    \partial_b E(a, b) &= 0.
#    \end{aligned}
#    $$
#    Use the Newton-Raphson method to solve these equations, and
#    then plot on the same graph the data points given and the approximating function.

# + nbgrader={"grade": true, "grade_id": "cell-e7bf77d0fa7b9458", "locked": false, "points": 0, "schema_version": 3, "solution": true, "task": false}
x = [0.0; 0.1; 0.2; 0.3; 0.4; 0.5; 0.6; 0.7; 0.8; 0.9; 1.0]
y = [0.6761488864859304; 0.6345697680852508; 0.6396283580587062; 0.6132010027973919;
     0.5906142598705267; 0.5718728461471725; 0.5524549902830562; 0.538938885654085;
     0.5373495476994958; 0.514904589752926; 0.49243437874655027]
f(a,b) = x -> a / (b+x)

### BEGIN SOLUTION
nothing
### END SOLUTION
# -

# 2. We have $n$ new points $(x_i, y_i)$ of an unknown function $y = f(x)$, and we want to approximate $f$ with an affine function
#    $$
#    \widetilde f(x) = ax+b
#    $$
#    by minimizing the sum of Euclidean distances between the points and the line defined by $\widetilde f$. Given that the distance between a point $(x_i, y_i)$ and the straight line is given by
#    $$
#    \frac{\lvert y_i - a x_i - b \rvert}{\sqrt{1+a^2}},
#    $$
#    the objective function to minimize is
#    $$
#    J(a, b) := \sum_{i=1}^{n} \frac{ \left( y_i - a x_i - b \right)^2 }{1+a^2}
#    $$
#    Find the optimal parameters $a$ and $b$ using the Newton-Raphson method and plot on the same graph the straight line $\tilde f$ along with the data points.

# + nbgrader={"grade": true, "grade_id": "cell-391efd40e50910f9", "locked": false, "points": 0, "schema_version": 3, "solution": true, "task": false}
x = [0.0; 0.1; 0.2; 0.3; 0.4; 0.5; 0.6; 0.7; 0.8; 0.9; 1.0]
y = [-0.9187980789440975; -0.6159791344678258; -0.25568734869121856;
     -0.14269370171581808; 0.3094396057228459; 0.6318327173549161;
     0.8370437988106428; 1.0970402798788812; 1.6057799131867696;
     1.869090784869698; 2.075369730726694]
f(a,b) = x -> a*x+b

### BEGIN SOLUTION
nothing
### END SOLUTION

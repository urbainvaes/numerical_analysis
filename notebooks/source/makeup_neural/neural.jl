# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.19.6
#   kernelspec:
#     display_name: Julia 1.12
#     language: julia
#     name: julia-1.12
# ---

# ## Make-up session: a neural network from scratch
#
# In this session, we build and train neural networks using only `LinearAlgebra` and the numerical tools seen in the course:
# no machine learning library.
# Along the way, we will meet familiar friends:
#
# - **Dual numbers** (automatic differentiation) to compute gradients;
# - **Gradient descent**, which is nothing but Richardson's iteration in disguise;
# - **The Runge phenomenon**, to compare neural networks with polynomial interpolation;
# - **Fourier series**, to help a network paint a picture.
#
# | Part | Topic | Suggested time |
# |------|-------|----------------|
# | 1 | Dual numbers and gradients | 25 min |
# | 2 | A first network: fitting the Runge function | 35 min |
# | 3 | Backpropagation | 40 min |
# | 4 | Learning to separate two spirals | 25 min |
# | 5 | Painting with a neural network | 30 min |
#
# **On using AI assistants.**
# You are welcome to use an AI assistant for the "plumbing" (plotting, animations, Julia syntax).
# However, the numerical core of each exercise (dual number arithmetic, gradient descent, backpropagation, Adam)
# should be written by you, and you should be able to explain every line of it.
# A good use of an assistant: once your code works, ask it to *critique* your implementation, and check whether its claims are correct.

# +
using LaTeXStrings
using LinearAlgebra
using Plots
using Polynomials
using Random
# -

# ### <font color='orange'>[Exercise 1]</font> Dual numbers and gradients
#
# Recall that a dual number is an expression of the form $a + b \varepsilon$, with $\varepsilon^2 = 0$.
# For any smooth function $f$, it holds that
# $$
# f(a + b \varepsilon) = f(a) + b f'(a) \varepsilon,
# $$
# so evaluating a function on a dual number gives its derivative for free.
#
# 1. Complete the definition of the `Dual` type below by implementing
#    the operations `-` (unary), `*`, `/`, `exp` and `tanh`.

# + nbgrader={"grade": false, "grade_id": "dual-ops", "locked": false, "schema_version": 3, "solution": true, "task": false}
import Base: +, -, *, /, exp, tanh, convert, promote_rule

struct Dual <: Number
    v::Float64  # Value
    d::Float64  # Derivative
end

+(x::Dual, y::Dual) = Dual(x.v + y.v, x.d + y.d)
-(x::Dual, y::Dual) = Dual(x.v - y.v, x.d - y.d)
### BEGIN SOLUTION
-(x::Dual) = Dual(-x.v, -x.d)
*(x::Dual, y::Dual) = Dual(x.v * y.v, x.d * y.v + x.v * y.d)
/(x::Dual, y::Dual) = Dual(x.v / y.v, (x.d * y.v - x.v * y.d) / y.v^2)
exp(x::Dual) = Dual(exp(x.v), exp(x.v) * x.d)
tanh(x::Dual) = Dual(tanh(x.v), (1 - tanh(x.v)^2) * x.d)
### END SOLUTION

# These two lines make it possible to combine dual numbers with real numbers, e.g. `2 * Dual(1, 1)`
convert(::Type{Dual}, x::Real) = Dual(x, 0)
promote_rule(::Type{Dual}, ::Type{<:Real}) = Dual
# -

# +
@assert tanh(Dual(0.5, 1)).d ≈ 1 / cosh(0.5)^2
@assert (Dual(2, 1) / Dual(3, 0) * exp(Dual(1, 0))).d ≈ exp(1) / 3
@assert (-Dual(1, 2)).d == -2
# -

# 2. Write a function `gradient(f, θ)` that returns the gradient of a function $f\colon \mathbb R^n \to \mathbb R$ at $θ$.
#    <details>
#        <summary>
#            <em><font color='gray'>Hint (click to display)</font></em>
#        </summary>
#
#    The $i$-th partial derivative can be obtained by evaluating $f$ at the point
#    $θ + \varepsilon e_i$, where $e_i$ is the $i$-th vector of the canonical basis.
#    How many evaluations of $f$ are required to compute the full gradient?
#    </details>

# + nbgrader={"grade": false, "grade_id": "dual-gradient", "locked": false, "schema_version": 3, "solution": true, "task": false}
function gradient(f, θ)
    ### BEGIN SOLUTION
    n = length(θ)
    return [f([Dual(θ[j], i == j) for j in 1:n]).d for i in 1:n]
    ### END SOLUTION
end

rosenbrock(x) = (1 - x[1])^2 + 100(x[2] - x[1]^2)^2
∇rosenbrock(x) = [-2(1 - x[1]) - 400x[1]*(x[2] - x[1]^2), 200(x[2] - x[1]^2)]
@assert gradient(rosenbrock, [0.5, 0.3]) ≈ ∇rosenbrock([0.5, 0.3])
# -

# ### <font color='orange'>[Exercise 2]</font> A first network: fitting the Runge function
#
# A neural network with one hidden layer of width $m$ is a function of the form
# $$
# N_θ(x) = c + \sum_{j=1}^{m} a_j \tanh(w_j x + b_j),
# $$
# where $θ = (w_1, \dotsc, w_m, b_1, \dotsc, b_m, a_1, \dotsc, a_m, c) \in \mathbb R^{3m+1}$ are the parameters of the network.
# Compare with polynomial approximation, where the basis functions $1, x, x^2, \dotsc$ are *fixed*:
# here, the basis functions $x \mapsto \tanh(w_j x + b_j)$ themselves are learnt from the data.
#
# Given data points $(x_i, y_i)_{1 \leq i \leq n}$, *training* the network means minimizing the mean squared error
# $$
# J(θ) = \frac{1}{n} \sum_{i=1}^{n} \bigl(N_θ(x_i) - y_i\bigr)^2.
# $$
# The simplest approach is gradient descent:
# $$
# θ_{k+1} = θ_k - η \nabla J(θ_k),
# $$
# where $η > 0$ is called the *learning rate*.
#
# We consider the Runge function $f(x) = 1/(1 + 25x^2)$, sampled at $n = 30$ equispaced points in $[-1, 1]$.
#
# 1. Complete the function `network(θ, x)`.

# + nbgrader={"grade": false, "grade_id": "network-1d", "locked": false, "schema_version": 3, "solution": true, "task": false}
m = 10
function network(θ, x)
    w, b, a, c = θ[1:m], θ[m+1:2m], θ[2m+1:3m], θ[3m+1]
    ### BEGIN SOLUTION
    return c + sum(a[j] * tanh(w[j]*x + b[j]) for j in 1:m)
    ### END SOLUTION
end

f(x) = 1 / (1 + 25x^2)
xs = LinRange(-1, 1, 30)
ys = f.(xs)
J(θ) = sum((network(θ, x) - y)^2 for (x, y) in zip(xs, ys)) / length(xs)
# -

# 2. Train the network by gradient descent with $η = 0.2$ and $5000$ iterations,
#    starting from random parameters `θ = randn(3m + 1)`.
#    Store the value of $J(θ_k)$ at each iteration in a vector `losses`,
#    and plot it using a logarithmic scale for the y-axis.

# + nbgrader={"grade": false, "grade_id": "train-1d", "locked": false, "schema_version": 3, "solution": true, "task": false}
Random.seed!(0)
θ = randn(3m + 1)
η, niter = 0.2, 5000

### BEGIN SOLUTION
losses = Float64[]
for k in 1:niter
    global θ -= η * gradient(J, θ)
    push!(losses, J(θ))
end
plot(losses, yscale=:log10, xlabel=L"k", ylabel=L"J(\theta_k)", label=nothing)
### END SOLUTION
# -

# 3. On the same graph, plot the Runge function, the trained network,
#    and the polynomial interpolating the data points.
#    Which approximation is better? Is the comparison fair?
#    <details>
#        <summary>
#            <em><font color='gray'>Hint (click to display)</font></em>
#        </summary>
#
#    The interpolating polynomial can be computed with `Polynomials.fit(xs, ys)`.
#    </details>

# + nbgrader={"grade": false, "grade_id": "runge-plot", "locked": false, "schema_version": 3, "solution": true, "task": false}
### BEGIN SOLUTION
xfine = LinRange(-1, 1, 500)
p = fit(xs, ys)
plot(xfine, f.(xfine), label="Runge function", ylims=(-0.5, 1.5))
plot!(xfine, network.(Ref(θ), xfine), label="Neural network")
plot!(xfine, p.(xfine), label="Interpolating polynomial")
scatter!(xs, ys, label="Data")
### END SOLUTION
# -

# 4. **Learning rate and Richardson's iteration.**
#    For a quadratic function $J(θ) = \frac{1}{2} θ^T \mathsf A θ - \mathsf b^T θ$ with $\mathsf A$ symmetric positive definite,
#    show that gradient descent coincides with Richardson's iteration for the system $\mathsf A θ = \mathsf b$.
#    For which values of $η$ does it converge?
#
#    Our function $J$ is not quadratic, but close to a minimizer it is approximately quadratic with $\mathsf A = \nabla^2 J(θ_*)$.
#    Rerun the training for $η \in \{0.2, 0.5, 1, 2\}$ and plot the loss curves.
#    What do you observe?

# + nbgrader={"grade": false, "grade_id": "learning-rate", "locked": false, "schema_version": 3, "solution": true, "task": false}
### BEGIN SOLUTION
plot(yscale=:log10, xlabel=L"k", ylabel=L"J(\theta_k)")
for η in (0.2, 0.5, 1, 2)
    local θ = randn(MersenneTwister(0), 3m + 1)
    local losses = Float64[]
    for k in 1:1000
        θ -= η * gradient(J, θ)
        push!(losses, J(θ))
    end
    plot!(filter(isfinite, losses), label="η = $η")
end
plot!()
### END SOLUTION
# -

# ### <font color='orange'>[Exercise 3]</font> Backpropagation
#
# Computing the gradient with dual numbers requires one evaluation of $J$ *per parameter*.
# This is fine for 31 parameters, but modern networks have billions of them!
# In this exercise, we implement *reverse-mode* automatic differentiation, known as **backpropagation** in this context,
# which computes the full gradient at a cost comparable to a few evaluations of $J$.
#
# We now consider deep networks, also called *multilayer perceptrons*,
# acting on a batch of inputs stored as columns of a matrix $\mathsf X \in \mathbb R^{d_0 \times n}$:
# $$
# \mathsf A_0 = \mathsf X, \qquad
# \mathsf Z_\ell = \mathsf W_\ell \mathsf A_{\ell-1} + \mathbf b_\ell, \qquad
# \mathsf A_\ell = \tanh(\mathsf Z_\ell) \quad \text{for } \ell = 1, \dotsc, L-1,
# \qquad \mathsf A_L = \mathsf Z_L.
# $$
# Here $\mathsf W_\ell \in \mathbb R^{d_\ell \times d_{\ell-1}}$ and $\mathbf b_\ell \in \mathbb R^{d_\ell}$,
# and the vector $\mathbf b_\ell$ is added to every column.
# The parameters are stored as a vector of tuples `params = [(W₁, b₁), …, (W_L, b_L)]`.
#
# Given targets $\mathsf Y \in \mathbb R^{d_L \times n}$, the loss is
# $J = \frac{1}{n} \lVert \mathsf A_L - \mathsf Y \rVert_{\rm F}^2$.
# Using the chain rule, show that the matrices $\mathsf G_\ell := \partial J / \partial \mathsf Z_\ell$ satisfy
# $$
# \mathsf G_L = \frac{2}{n} (\mathsf A_L - \mathsf Y), \qquad
# \mathsf G_{\ell-1} = (\mathsf W_\ell^T \mathsf G_\ell) \odot (1 - \mathsf A_{\ell-1}^{\odot 2}),
# $$
# where $\odot$ denotes the elementwise product, and that
# $$
# \frac{\partial J}{\partial \mathsf W_\ell} = \mathsf G_\ell \mathsf A_{\ell-1}^T, \qquad
# \frac{\partial J}{\partial \mathbf b_\ell} = \mathsf G_\ell \mathbf 1.
# $$
# Notice that the gradients are propagated *backwards*, from the last layer to the first.

# +
init(sizes) = [(randn(sizes[l+1], sizes[l]) / sqrt(sizes[l]), zeros(sizes[l+1]))
               for l in 1:length(sizes)-1]

function forward(params, X)
    A = X
    for (ℓ, (W, b)) in enumerate(params)
        Z = W*A .+ b
        A = ℓ < length(params) ? tanh.(Z) : Z
    end
    return A
end

# Conversion between `params` and a flat vector of parameters
flatten(params) = vcat([vcat(vec(W), b) for (W, b) in params]...)
function unflatten(θ, params)
    result, i = [], 0
    for (W, b) in params
        Wnew = reshape(θ[i+1:i+length(W)], size(W))
        i += length(W)
        bnew = θ[i+1:i+length(b)]
        i += length(b)
        push!(result, (Wnew, bnew))
    end
    return result
end
# -

# 1. Complete the function `loss_and_gradient`, which returns the loss and its gradient
#    in the same format as `params`.

# + nbgrader={"grade": false, "grade_id": "backprop", "locked": false, "schema_version": 3, "solution": true, "task": false}
function loss_and_gradient(params, X, Y)
    L, n = length(params), size(X, 2)

    # Forward pass: store all the activations
    As = [X]
    for (ℓ, (W, b)) in enumerate(params)
        Z = W*As[end] .+ b
        push!(As, ℓ < L ? tanh.(Z) : Z)
    end
    loss = sum(abs2, As[end] - Y) / n

    # Backward pass
    grads = similar(params)
    ### BEGIN SOLUTION
    G = 2(As[end] - Y) / n
    for ℓ in L:-1:1
        W, b = params[ℓ]
        grads[ℓ] = (G * As[ℓ]', vec(sum(G, dims=2)))
        if ℓ > 1
            G = (W' * G) .* (1 .- As[ℓ].^2)
        end
    end
    ### END SOLUTION
    return loss, grads
end
# -

# 2. Check your implementation against the gradient computed with dual numbers,
#    for a small network with layer sizes $[1, 4, 4, 1]$ and random data.

# + nbgrader={"grade": false, "grade_id": "backprop-check", "locked": false, "schema_version": 3, "solution": true, "task": false}
Random.seed!(1)
params = init([1, 4, 4, 1])
X, Y = randn(1, 20), randn(1, 20)

### BEGIN SOLUTION
Jflat(θ) = sum((forward(unflatten(θ, params), X) - Y).^2) / size(X, 2)
_, grads = loss_and_gradient(params, X, Y)
g_backprop = flatten(grads)
g_dual = gradient(Jflat, flatten(params))
@show norm(g_backprop - g_dual)
### END SOLUTION
# -

# 3. For networks with layer sizes $[1, w, w, 1]$ with $w \in \{4, 8, 16, 32, 64\}$,
#    measure the time required to compute the gradient with both methods,
#    and plot it as a function of the number of parameters on a log-log scale.
#    What are the observed scalings?

# + nbgrader={"grade": false, "grade_id": "backprop-timing", "locked": false, "schema_version": 3, "solution": true, "task": false}
### BEGIN SOLUTION
widths = [4, 8, 16, 32, 64]
nparams, t_dual, t_backprop = Int[], Float64[], Float64[]
for w in widths
    local params = init([1, w, w, 1])
    Jw(θ) = sum((forward(unflatten(θ, params), X) - Y).^2) / size(X, 2)
    local θ = flatten(params)
    gradient(Jw, θ); loss_and_gradient(params, X, Y)  # Compilation
    push!(nparams, length(θ))
    push!(t_dual, @elapsed gradient(Jw, θ))
    push!(t_backprop, @elapsed loss_and_gradient(params, X, Y))
end
plot(nparams, t_dual, xscale=:log10, yscale=:log10, marker=:o, label="Dual numbers")
plot!(nparams, t_backprop, marker=:o, label="Backpropagation")
plot!(xlabel="Number of parameters", ylabel="Time [s]", legend=:topleft)
### END SOLUTION
# -

# ### <font color='orange'>[Exercise 4]</font> Learning to separate two spirals
#
# We now train a network to *classify* points in the plane: each point belongs to one of two intertwined spirals,
# labelled $+1$ and $-1$.
# The network outputs a real number, and we predict the label by its sign.

# +
function spirals(n; noise=0.1, rng=MersenneTwister(1))
    X, Y = zeros(2, 2n), zeros(1, 2n)
    for i in 1:n, (j, s) in enumerate((1, -1))
        t = sqrt(rand(rng)) * 3π
        X[:, 2(i-1)+j] = s * t * [cos(t), sin(t)] / 3π + noise/3π * randn(rng, 2)
        Y[1, 2(i-1)+j] = s
    end
    return X, Y
end

X, Y = spirals(300)
scatter(X[1, :], X[2, :], zcolor=vec(Y), c=:coolwarm, aspect_ratio=:equal, label=nothing, colorbar=false)
# -

# Plain gradient descent is slow on this problem.
# In practice, neural networks are trained with the **Adam** algorithm,
# which uses *momentum* (an exponential moving average of the gradients)
# and a *per-parameter* step size (based on a moving average of the squared gradients):
# $$
# \begin{aligned}
# m_k &= β_1 m_{k-1} + (1-β_1) g_k, \\
# v_k &= β_2 v_{k-1} + (1-β_2) g_k^{2}, \\
# θ_k &= θ_{k-1} - η \frac{m_k / (1 - β_1^k)}{\sqrt{v_k / (1 - β_2^k)} + 10^{-8}},
# \end{aligned}
# $$
# where $g_k$ is the gradient at $θ_{k-1}$, all operations are performed elementwise,
# and $m_0 = v_0 = 0$.
# The divisions by $1 - β^k$ correct the bias introduced by initializing the averages at zero.
#
# 1. Complete the function `train!`, which updates `params` *in place* and returns the vector of losses.
#    The optional `callback` is a function called at the end of each iteration with argument $k$,
#    which will be useful for animations.

# + nbgrader={"grade": false, "grade_id": "adam", "locked": false, "schema_version": 3, "solution": true, "task": false}
function train!(params, X, Y; η=1e-2, nsteps=1000, β₁=0.9, β₂=0.999, callback=nothing)
    ms = [(zero(W), zero(b)) for (W, b) in params]
    vs = [(zero(W), zero(b)) for (W, b) in params]
    losses = Float64[]
    for k in 1:nsteps
        loss, grads = loss_and_gradient(params, X, Y)
        push!(losses, loss)
        ### BEGIN SOLUTION
        for ℓ in eachindex(params), i in 1:2
            θ, g, m, v = params[ℓ][i], grads[ℓ][i], ms[ℓ][i], vs[ℓ][i]
            @. m = β₁*m + (1 - β₁)*g
            @. v = β₂*v + (1 - β₂)*g^2
            @. θ -= η * (m / (1 - β₁^k)) / (sqrt(v / (1 - β₂^k)) + 1e-8)
        end
        ### END SOLUTION
        callback === nothing || callback(k)
    end
    return losses
end
# -

# 2. Train a network with layer sizes $[2, 32, 32, 1]$ for 3000 iterations with $η = 10^{-2}$.
#    Compute the fraction of correctly classified points.

# + nbgrader={"grade": false, "grade_id": "spirals-train", "locked": false, "schema_version": 3, "solution": true, "task": false}
Random.seed!(2)
params = init([2, 32, 32, 1])

### BEGIN SOLUTION
losses = train!(params, X, Y; η=1e-2, nsteps=3000)
accuracy = sum(sign.(forward(params, X)) .== Y) / length(Y)
println("Accuracy: $accuracy")
plot(losses, yscale=:log10, xlabel=L"k", label="Loss")
### END SOLUTION
# -

# 3. Create an animation showing how the *decision boundary*
#    (the zero level set of the network) evolves during training.
#    <details>
#        <summary>
#            <em><font color='gray'>Hint (click to display)</font></em>
#        </summary>
#
#    Use `Plots.Animation()`, and call `frame(anim)` in the callback every 50 iterations.
#    The network can be evaluated on a grid at once by building a $2 \times N$ matrix of grid points.
#    Use `gif(anim, "spirals.gif")` to save the animation.
#    </details>

# + nbgrader={"grade": false, "grade_id": "spirals-anim", "locked": false, "schema_version": 3, "solution": true, "task": false}
### BEGIN SOLUTION
grid = LinRange(-1.2, 1.2, 100)
G = hcat([[x, y] for y in grid for x in grid]...)

function plot_boundary(params, k)
    Z = reshape(forward(params, G), length(grid), length(grid))'
    contourf(grid, grid, clamp.(Z, -1.5, 1.5), c=:coolwarm, levels=20, colorbar=false)
    contour!(grid, grid, Z, levels=[0], c=:black, lw=2)
    scatter!(X[1, :], X[2, :], zcolor=vec(Y), c=:coolwarm, ms=3, label=nothing)
    plot!(aspect_ratio=:equal, title="Iteration $k", xlims=extrema(grid), ylims=extrema(grid))
end

Random.seed!(2)
params = init([2, 32, 32, 1])
anim = Animation()
train!(params, X, Y; η=1e-2, nsteps=3000,
       callback = k -> k % 50 == 1 && (plot_boundary(params, k); frame(anim)))
gif(anim, "spirals.gif", fps=10)
### END SOLUTION
# -

# 4. Try a network with a single hidden layer of width 32, and a network with layer sizes $[2, 8, 8, 8, 1]$.
#    How do the number of parameters and the final accuracy compare?

# + nbgrader={"grade": false, "grade_id": "spirals-arch", "locked": false, "schema_version": 3, "solution": true, "task": false}
### BEGIN SOLUTION
for sizes in ([2, 32, 1], [2, 8, 8, 8, 1], [2, 32, 32, 1])
    Random.seed!(2)
    local params = init(sizes)
    train!(params, X, Y; η=1e-2, nsteps=3000)
    local accuracy = sum(sign.(forward(params, X)) .== Y) / length(Y)
    println("$sizes: $(length(flatten(params))) parameters, accuracy $accuracy")
end
### END SOLUTION
# -

# ### <font color='orange'>[Exercise 5]</font> Painting with a neural network
#
# A neural network can also represent an image:
# we look for a function $N_θ\colon [-1, 1]^2 \to \mathbb R$ such that $N_θ(x, y)$ is the gray level of the pixel at $(x, y)$.
# This idea, called *implicit neural representation*, is at the heart of recent techniques for 3D scene reconstruction.

# +
function smiley(x, y)
    face = hypot(x, y) < 0.9
    eyes = hypot(x - 0.35, y - 0.3) < 0.13 || hypot(x + 0.35, y - 0.3) < 0.13
    mouth = 0.45 < hypot(x, y + 0.05) < 0.6 && y < -0.15
    return float(face && !eyes && !mouth)
end

npix = 64
pixels = LinRange(-1, 1, npix)
P = hcat([[x, y] for y in pixels for x in pixels]...)  # 2 × npix² matrix of pixel coordinates
Yimg = [smiley(x, y) for x in pixels, y in pixels] |> vec |> permutedims
show_image(v, n=npix) = heatmap(reshape(v, n, n)', c=:grays, aspect_ratio=:equal, axis=false, colorbar=false, yflip=false)
show_image(Yimg)
# -

# 1. Train a network with layer sizes $[2, 64, 64, 1]$ on the image, using Adam with $η = 3 \times 10^{-3}$ and 1000 iterations.
#    Display the resulting image. What do you observe?

# + nbgrader={"grade": false, "grade_id": "paint-plain", "locked": false, "schema_version": 3, "solution": true, "task": false}
### BEGIN SOLUTION
Random.seed!(3)
params = init([2, 64, 64, 1])
losses_plain = train!(params, P, Yimg; η=3e-3, nsteps=1000)
show_image(forward(params, P))
### END SOLUTION
# -

# The result is blurry: neural networks have a *spectral bias*, meaning that they learn low frequencies first,
# and struggle to represent sharp features.
# A simple remedy is to feed the network with **Fourier features** of the input instead of the input itself:
# $$
# γ(\mathbf x) = \begin{pmatrix} \sin(2π \mathsf B \mathbf x) \\ \cos(2π \mathsf B \mathbf x) \end{pmatrix} \in \mathbb R^{2q},
# $$
# where $\mathsf B \in \mathbb R^{q \times 2}$ is a fixed random matrix with i.i.d. $\mathcal N(0, σ^2)$ entries.
#
# 2. Train networks with layer sizes $[128, 64, 64, 1]$ on the Fourier features, with $q = 64$ and $σ \in \{1, 3, 10\}$.
#    Compare the loss curves.

# + nbgrader={"grade": false, "grade_id": "paint-fourier", "locked": false, "schema_version": 3, "solution": true, "task": false}
fourier_features(P, B) = vcat(sin.(2π*B*P), cos.(2π*B*P))

### BEGIN SOLUTION
Random.seed!(3)
results = Dict()
plot(losses_plain, label="No Fourier features", yscale=:log10, xlabel=L"k")
for σ in (1, 3, 10)
    local B = σ * randn(64, 2)
    local params = init([128, 64, 64, 1])
    local losses = train!(params, fourier_features(P, B), Yimg; η=3e-3, nsteps=1000)
    results[σ] = (B, params)
    plot!(losses, label="σ = $σ")
end
plot!()
### END SOLUTION
# -

# 3. A network is a function defined on the whole square, not only at the pixels.
#    Evaluate each of the trained networks on a finer grid of $256 \times 256$ points, and display the results.
#    Which value of $σ$ gives the best image? Relate your observations to what you know about interpolation.

# + nbgrader={"grade": false, "grade_id": "paint-superres", "locked": false, "schema_version": 3, "solution": true, "task": false}
### BEGIN SOLUTION
nfine = 256
fine = LinRange(-1, 1, nfine)
Pfine = hcat([[x, y] for y in fine for x in fine]...)
plots = [show_image(forward(params, fourier_features(Pfine, B)), nfine) for (B, params) in (results[σ] for σ in (1, 3, 10))]
plot(plots..., layout=(1, 3), size=(900, 300), title=["σ = 1" "σ = 3" "σ = 10"])
### END SOLUTION
# -

# With large $σ$, the network fits the pixels perfectly but oscillates in between:
# this is overfitting, a close cousin of the Runge phenomenon!
#
# 4. Create an animation of the network learning the image (for $σ = 3$), and enjoy.

# + nbgrader={"grade": false, "grade_id": "paint-anim", "locked": false, "schema_version": 3, "solution": true, "task": false}
### BEGIN SOLUTION
Random.seed!(4)
B = 3 * randn(64, 2)
F = fourier_features(P, B)
params = init([128, 64, 64, 1])
anim = Animation()
train!(params, F, Yimg; η=3e-3, nsteps=500,
       callback = k -> (k ≤ 50 || k % 10 == 0) && (show_image(forward(params, F)); frame(anim)))
gif(anim, "smiley.gif", fps=15)
### END SOLUTION
# -

# ### Going further
#
# If you finish early, here are a few ideas. These are good opportunities to work with an AI assistant,
# but remember to check that you understand (and agree with!) what it produces.
#
# - **Your own picture.** Load a photograph with the `Images` and `FileIO` packages and learn it with a network with 3 outputs (RGB).
# - **Classification with cross-entropy.** Replace the mean squared error by the logistic loss $\log(1 + e^{-y N_θ(x)})$. Which lines of `loss_and_gradient` change?
# - **Handwritten digits.** Download the MNIST dataset with the `MLDatasets` package, and train a network with layer sizes $[784, 64, 10]$. Draw a digit yourself and see whether your network recognizes it.
# - **Stochastic gradient descent.** Instead of using all the data at each iteration, use a random subset (a *mini-batch*). How does this affect the cost per iteration and the convergence?
# - **Reverse-mode automatic differentiation in general.** Our backpropagation is specific to multilayer perceptrons. How would you design a `Dual`-like type that records operations, so that gradients of *any* function can be computed backwards?

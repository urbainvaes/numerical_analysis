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

# ## Make-up session: real-time smoke
#
# In this session, we write a smoke simulator based on the *stable fluids* method of Jos Stam (1999),
# which is used in many video games and visual effects.
# The method combines three ingredients seen in the course:
#
# - **Interpolation**, to transport quantities along the flow (semi-Lagrangian advection);
# - **ODE solvers**, to trace the trajectories of fluid particles backwards in time;
# - **Iterative methods for linear systems** (Jacobi, Gauss–Seidel, SOR), to make the flow incompressible.
#
# In the last part, we take the simulation to the browser, in real time.
#
# | Part | Topic | Suggested time |
# |------|-------|----------------|
# | 1 | Iterative solvers for the Poisson equation | 45 min |
# | 2 | Semi-Lagrangian advection | 30 min |
# | 3 | Projection and the full simulator | 40 min |
# | 4 | Going real-time in the browser | 45 min |
#
# **On using AI assistants.**
# You are welcome to use an AI assistant for the "plumbing" (plotting, animations, and in Part 4, JavaScript and HTML).
# However, the numerical core of Parts 1–3 (solvers, interpolation, projection) should be written by you,
# and you should be able to explain every line of it.

# +
using LaTeXStrings
using LinearAlgebra
using Plots
using Random
using SparseArrays
# -

# ### The model
#
# We model the air as an incompressible, inviscid fluid with velocity field $\mathbf u = (u, v)$,
# governed by the Euler equations
# $$
# \partial_t \mathbf u + (\mathbf u \cdot \nabla) \mathbf u = - \nabla p + \mathbf f,
# \qquad
# \nabla \cdot \mathbf u = 0,
# $$
# where $p$ is the pressure and $\mathbf f$ an external force (here, buoyancy: hot smoke rises).
# The smoke density $s$ is simply transported by the flow:
# $$
# \partial_t s + \mathbf u \cdot \nabla s = 0.
# $$
# Stam's method splits each time step into three substeps:
#
# 1. **Forces:** $\mathbf u \leftarrow \mathbf u + Δt\, \mathbf f$;
# 2. **Advection:** transport $\mathbf u$ and $s$ along the flow, ignoring pressure;
# 3. **Projection:** find $φ$ such that $\mathbf u - \nabla φ$ is divergence-free, i.e. solve the *Poisson equation*
#    $$
#    \Delta φ = \nabla \cdot \mathbf u,
#    $$
#    and set $\mathbf u \leftarrow \mathbf u - \nabla φ$.
#
# Step 3 is by far the most expensive one, and it is where numerical linear algebra comes in.

# ### <font color='orange'>[Exercise 1]</font> Iterative solvers for the Poisson equation
#
# We discretize the square $[0, 1]^2$ into $n \times n$ cells of size $h = 1/n$.
# The unknown $φ_{i,j}$ lives at the center of cell $(i, j)$, and we use the boundary condition $φ = 0$ outside the domain
# (the box is surrounded by air at rest).
# The five-point discretization of $\Delta φ = r$ reads
# $$
# \frac{φ_{i-1,j} + φ_{i+1,j} + φ_{i,j-1} + φ_{i,j+1} - 4 φ_{i,j}}{h^2} = r_{i,j},
# \qquad 1 \leq i, j \leq n,
# $$
# with the convention that $φ_{i,j} = 0$ if $(i, j)$ is outside $\{1, \dotsc, n\}^2$.
# This is a linear system with $n^2$ unknowns, which we never assemble:
# all the solvers below are *matrix-free* and work directly on $n \times n$ arrays.

# +
# Value of φ[i, j], or 0 outside the domain
@inline at(φ, i, j) = (1 ≤ i ≤ size(φ, 1) && 1 ≤ j ≤ size(φ, 2)) ? φ[i, j] : 0.0

function residual(φ, r, h)
    n₁, n₂ = size(φ)
    return [r[i, j] - (at(φ, i-1, j) + at(φ, i+1, j) + at(φ, i, j-1) + at(φ, i, j+1) - 4φ[i, j]) / h^2
            for i in 1:n₁, j in 1:n₂]
end
# -

# 1. Write functions `jacobi!(φ, r, h, niter)` and `sor!(φ, r, h, niter; ω)`
#    that perform `niter` iterations of the Jacobi and SOR methods, respectively, *in place*.
#    Recall that the Gauss–Seidel method corresponds to SOR with $ω = 1$.
#    <details>
#        <summary>
#            <em><font color='gray'>Hint (click to display)</font></em>
#        </summary>
#
#    Solving the $(i, j)$ equation for $φ_{i,j}$ gives
#    $$
#    φ_{i,j} = \frac{1}{4} \left(φ_{i-1,j} + φ_{i+1,j} + φ_{i,j-1} + φ_{i,j+1} - h^2 r_{i,j}\right).
#    $$
#    The Jacobi method uses the old values on the right-hand side, so it needs a copy of $φ$;
#    the Gauss–Seidel method overwrites $φ_{i,j}$ immediately.
#    </details>

# + nbgrader={"grade": false, "grade_id": "solvers", "locked": false, "schema_version": 3, "solution": true, "task": false}
function jacobi!(φ, r, h, niter)
    ### BEGIN SOLUTION
    n₁, n₂ = size(φ)
    ψ = similar(φ)
    for _ in 1:niter
        for j in 1:n₂, i in 1:n₁
            ψ[i, j] = (at(φ, i-1, j) + at(φ, i+1, j) + at(φ, i, j-1) + at(φ, i, j+1) - h^2 * r[i, j]) / 4
        end
        φ .= ψ
    end
    ### END SOLUTION
    return φ
end

function sor!(φ, r, h, niter; ω=1)
    ### BEGIN SOLUTION
    n₁, n₂ = size(φ)
    for _ in 1:niter, j in 1:n₂, i in 1:n₁
        φ_gs = (at(φ, i-1, j) + at(φ, i+1, j) + at(φ, i, j-1) + at(φ, i, j+1) - h^2 * r[i, j]) / 4
        φ[i, j] = (1 - ω) * φ[i, j] + ω * φ_gs
    end
    ### END SOLUTION
    return φ
end
# -

# +
Random.seed!(0)
n = 32
h = 1/n
r = randn(n, n)
@assert norm(residual(jacobi!(zeros(n, n), r, h, 5000), r, h)) < 1e-6 * norm(r)
@assert norm(residual(sor!(zeros(n, n), r, h, 2000), r, h)) < 1e-6 * norm(r)
# -

# 2. The eigenvalues of the Jacobi iteration matrix for this problem are
#    $$
#    μ_{k, \ell} = \frac{1}{2} \left( \cos \frac{k π}{n+1} + \cos \frac{\ell π}{n+1} \right), \qquad 1 \leq k, \ell \leq n.
#    $$
#    Deduce the spectral radius $ρ_J$ of the Jacobi method, and recall from the course
#    the spectral radius $ρ_{GS}$ of the Gauss–Seidel method,
#    the optimal relaxation parameter $ω_*$ of the SOR method, and the corresponding spectral radius $ρ_{SOR}$.
#
#    Then, starting from $φ = 0$, plot the norm of the residual as a function of the iteration number $k \leq 1000$
#    for the Jacobi, Gauss–Seidel and optimal SOR methods, using a logarithmic scale for the y-axis.
#    Add the theoretical rates $ρ^k$ to the plot.

# + nbgrader={"grade": false, "grade_id": "rates", "locked": false, "schema_version": 3, "solution": true, "task": false}
### BEGIN SOLUTION
ρ_J = cos(π/(n+1))
ρ_GS = ρ_J^2
ω_opt = 2 / (1 + sqrt(1 - ρ_J^2))
ρ_SOR = ω_opt - 1

function residuals(method!, niter)
    φ = zeros(n, n)
    res = [norm(r)]
    for k in 1:niter
        method!(φ)
        push!(res, norm(residual(φ, r, h)))
    end
    return res
end

K = 1000
ks = 0:K
plot(yscale=:log10, xlabel=L"k", ylabel="Residual", legend=:bottomleft, ylims=(1e-12, 1e2))
plot!(ks, residuals(φ -> jacobi!(φ, r, h, 1), K), label="Jacobi", c=1)
plot!(ks, residuals(φ -> sor!(φ, r, h, 1), K), label="Gauss-Seidel", c=2)
plot!(ks, residuals(φ -> sor!(φ, r, h, 1; ω=ω_opt), K), label="SOR, ω = $(round(ω_opt, digits=3))", c=3)
plot!(ks, norm(r) * ρ_J.^ks, ls=:dash, c=1, label=L"\rho_J^k")
plot!(ks, norm(r) * ρ_GS.^ks, ls=:dash, c=2, label=L"\rho_{GS}^k")
plot!(ks, norm(r) * ρ_SOR.^ks, ls=:dash, c=3, label=L"\rho_{SOR}^k")
### END SOLUTION
# -

# 3. For $ω \in [1, 1.99]$, compute the number of SOR iterations needed to reduce the residual by a factor $10^6$,
#    and plot it as a function of $ω$. Compare the empirical optimum with $ω_*$.

# + nbgrader={"grade": false, "grade_id": "omega-scan", "locked": false, "schema_version": 3, "solution": true, "task": false}
### BEGIN SOLUTION
function iterations_needed(ω; tol=1e-6, maxiter=5000)
    φ = zeros(n, n)
    for k in 1:maxiter
        sor!(φ, r, h, 1; ω)
        norm(residual(φ, r, h)) < tol * norm(r) && return k
    end
    return maxiter
end

ωs = 1:0.01:1.99
plot(ωs, iterations_needed.(ωs), xlabel=L"\omega", ylabel="Iterations", label=nothing, yscale=:log10)
vline!([ω_opt], ls=:dash, label=L"\omega_*")
### END SOLUTION
# -

# 4. In two dimensions, sparse direct solvers are very efficient.
#    Assemble the matrix of the system using `kron` and `spdiagm`, and compare the time needed by `A \ b`
#    with the time needed by optimal SOR to reach a relative residual of $10^{-6}$.
#    Can you think of reasons to use an iterative method in a real-time simulation anyway?
#    <details>
#        <summary>
#            <em><font color='gray'>Hint (click to display)</font></em>
#        </summary>
#
#    The matrix of the one-dimensional problem is `D = spdiagm(-1 => ones(n-1), 0 => -2ones(n), 1 => ones(n-1)) / h^2`,
#    and the two-dimensional matrix is $\mathsf I \otimes \mathsf D + \mathsf D \otimes \mathsf I$.
#    Be careful with the ordering of the unknowns: `vec(φ)` stacks the columns of `φ`.
#    </details>

# + nbgrader={"grade": false, "grade_id": "direct", "locked": false, "schema_version": 3, "solution": true, "task": false}
### BEGIN SOLUTION
D = spdiagm(-1 => ones(n-1), 0 => -2ones(n), 1 => ones(n-1)) / h^2
A = kron(I(n), D) + kron(D, I(n))
φ_direct = reshape(A \ vec(r), n, n)
@assert norm(residual(φ_direct, r, h)) < 1e-8 * norm(r)

A \ vec(r); sor!(zeros(n, n), r, h, 1)  # Compilation
t_direct = @elapsed A \ vec(r)
t_sor = @elapsed sor!(zeros(n, n), r, h, iterations_needed(ω_opt); ω=ω_opt)
println("Direct: $t_direct s, SOR: $t_sor s")
### END SOLUTION
# -

# For later use, here is a function returning the matrix of the system.

# +
function poisson_matrix(n)
    D = spdiagm(-1 => ones(n-1), 0 => -2ones(n), 1 => ones(n-1)) * n^2
    return kron(I(n), D) + kron(D, I(n))
end
# -

# In a simulation, the solution $φ$ changes little from one time step to the next.
# An iterative method can be *warm-started* from the previous solution,
# and a fixed (small) number of iterations per frame gives a predictable cost, which matters for real-time applications.
# Iterative methods are also easy to parallelize on GPUs.

# ### <font color='orange'>[Exercise 2]</font> Semi-Lagrangian advection
#
# We store the velocity on a *staggered grid* (also called MAC grid):
# the horizontal velocity $u$ is stored at the centers of the vertical cell faces,
# and the vertical velocity $v$ at the centers of the horizontal cell faces.
# The smoke density $s$ and the potential $φ$ are stored at cell centers.
# For each field, we store the coordinates of the first grid point, called the *offset*:
#
# | Field | Array size | Position of entry $(i, j)$ | Offset |
# |-------|------------|-----------------------------|--------|
# | $u$ | $(n+1) \times n$ | $\bigl((i-1)h, (j-\frac12)h\bigr)$ | $(0, \frac h2)$ |
# | $v$ | $n \times (n+1)$ | $\bigl((i-\frac12)h, (j-1)h\bigr)$ | $(\frac h2, 0)$ |
# | $s$ | $n \times n$ | $\bigl((i-\frac12)h, (j-\frac12)h\bigr)$ | $(\frac h2, \frac h2)$ |
#
# Staggering makes the discrete divergence and gradient fit together exactly, as we will see in Exercise 3.

# +
mutable struct Fluid
    n::Int
    h::Float64
    u::Matrix{Float64}  # Horizontal velocity, (n+1) × n
    v::Matrix{Float64}  # Vertical velocity, n × (n+1)
    s::Matrix{Float64}  # Smoke density, n × n
    φ::Matrix{Float64}  # Potential, n × n (kept for warm starts)
    t::Float64          # Current time
end

Fluid(n) = Fluid(n, 1/n, zeros(n+1, n), zeros(n, n+1), zeros(n, n), zeros(n, n), 0)

offset_u(f) = (0, f.h/2)
offset_v(f) = (f.h/2, 0)
offset_s(f) = (f.h/2, f.h/2)
# -

# **Semi-Lagrangian advection.**
# The transport equation $\partial_t q + \mathbf u \cdot \nabla q = 0$ states that $q$ is constant along the trajectories of fluid particles.
# Therefore, the value of $q$ at a grid point $\mathbf x$ at time $t + Δt$ equals the value of $q$ at time $t$
# at the point $\mathbf x_{\rm back}$ where the particle came from.
# To compute $\mathbf x_{\rm back}$, we solve the ODE $\dot{\mathbf x} = \mathbf u(\mathbf x)$ *backwards* in time over $[0, Δt]$,
# and since $\mathbf x_{\rm back}$ is generally not a grid point, we evaluate $q$ there by *interpolation*.
#
# 1. Complete the function `sample(q, x, y, offset, h)`, which returns the bilinear interpolation of the field `q` at point $(x, y)$.
#    Points outside the grid should be projected onto the nearest boundary of the grid (use `clamp`).
#    <details>
#        <summary>
#            <em><font color='gray'>Hint (click to display)</font></em>
#        </summary>
#
#    Compute the fractional indices $σ = (x - o_x)/h + 1$ and $τ = (y - o_y)/h + 1$,
#    clamp them to the valid range,
#    and let $i = \lfloor σ \rfloor$, $a = σ - i$ (and similarly for $j$ and $b$). Then
#    $$
#    q(x, y) \approx (1-a)(1-b)\, q_{i,j} + a(1-b)\, q_{i+1,j} + (1-a)b\, q_{i,j+1} + ab\, q_{i+1,j+1}.
#    $$
#    Be careful with the case $σ = n_1$, where $i + 1$ is out of bounds.
#    </details>

# + nbgrader={"grade": false, "grade_id": "sample", "locked": false, "schema_version": 3, "solution": true, "task": false}
function sample(q, x, y, offset, h)
    n₁, n₂ = size(q)
    ### BEGIN SOLUTION
    σ = clamp((x - offset[1]) / h + 1, 1, n₁)
    τ = clamp((y - offset[2]) / h + 1, 1, n₂)
    i, j = min(floor(Int, σ), n₁ - 1), min(floor(Int, τ), n₂ - 1)
    a, b = σ - i, τ - j
    return (1-a)*(1-b)*q[i, j] + a*(1-b)*q[i+1, j] + (1-a)*b*q[i, j+1] + a*b*q[i+1, j+1]
    ### END SOLUTION
end

velocity(f, x, y) = (sample(f.u, x, y, offset_u(f), f.h), sample(f.v, x, y, offset_v(f), f.h))

# Test: bilinear interpolation is exact for bilinear functions
q = [2 + 3x - y + x*y for x in 0:0.1:1, y in 0:0.1:1]
@assert sample(q, 0.234, 0.789, (0, 0), 0.1) ≈ 2 + 3*0.234 - 0.789 + 0.234*0.789
@assert sample(q, 1.5, -1, (0, 0), 0.1) ≈ q[end, 1]
# -

# 2. Complete the function `advect(q, offset, f, Δt; rk2)`, which returns the advected field.
#    The backward trajectory should be computed with one step of the explicit Euler method if `rk2 = false`,
#    and with one step of the midpoint (RK2) method otherwise.

# + nbgrader={"grade": false, "grade_id": "advect", "locked": false, "schema_version": 3, "solution": true, "task": false}
function advect(q, offset, f, Δt; rk2=true)
    h = f.h
    result = similar(q)
    for j in axes(q, 2), i in axes(q, 1)
        x, y = offset[1] + (i-1)*h, offset[2] + (j-1)*h
        ### BEGIN SOLUTION
        u, v = velocity(f, x, y)
        if rk2
            u, v = velocity(f, x - Δt/2 * u, y - Δt/2 * v)
        end
        result[i, j] = sample(q, x - Δt*u, y - Δt*v, offset, h)
        ### END SOLUTION
    end
    return result
end
# -

# 3. We test the method on a rigid rotation $\mathbf u(x, y) = 2π \bigl(-(y - \tfrac12), x - \tfrac12\bigr)$, with period 1.
#    Advect the image below for one full turn with $Δt = 1/50$, with the Euler and RK2 backward trajectories.
#    Compute the Courant number $\max |\mathbf u| Δt / h$: an explicit finite-difference scheme would be unstable here!
#    What do you observe? Repeat with $Δt = 1/200$.

# +
n = 100
f = Fluid(n)
f.u .= [-2π * (y - 0.5) for x in (0:n) * f.h, y in ((1:n) .- 0.5) * f.h]
f.v .= [ 2π * (x - 0.5) for x in ((1:n) .- 0.5) * f.h, y in (0:n) * f.h]
cells = ((1:n) .- 0.5) / n
s₀ = [float(abs(x - 0.5) < 0.3 && abs(y - 0.7) < 0.1 || abs(x - 0.5) < 0.1 && 0.25 < y < 0.7) for x in cells, y in cells]
show_field(s; kwargs...) = heatmap(cells, cells, s', aspect_ratio=:equal, c=:inferno,
                                   xlims=(0, 1), ylims=(0, 1), colorbar=false; kwargs...)
show_field(s₀, title="Initial condition")
# -

# + nbgrader={"grade": false, "grade_id": "rotation", "locked": false, "schema_version": 3, "solution": true, "task": false}
### BEGIN SOLUTION
plots = []
for nsteps in (50, 200), rk2 in (false, true)
    local s = copy(s₀)
    for _ in 1:nsteps
        s = advect(s, offset_s(f), f, 1/nsteps; rk2)
    end
    courant = round(2π * sqrt(0.5) / nsteps / f.h, digits=1)
    push!(plots, show_field(s, title="$(rk2 ? "RK2" : "Euler"), Δt = 1/$nsteps, C ≈ $courant", titlefontsize=9))
end
plot(plots..., layout=(2, 2), size=(700, 700))
### END SOLUTION
# -

# The method is stable for any time step: this is why it is popular in video games.
# The price to pay is *numerical diffusion*: each interpolation step smooths the field a little.

# ### <font color='orange'>[Exercise 3]</font> Projection and the full simulator
#
# On the staggered grid, the discrete divergence at the center of cell $(i, j)$ is
# $$
# (\nabla_h \cdot \mathbf u)_{i,j} = \frac{u_{i+1,j} - u_{i,j} + v_{i,j+1} - v_{i,j}}{h},
# $$
# and the discrete gradient of $φ$ at the face between cells $(i-1, j)$ and $(i, j)$ is $(φ_{i,j} - φ_{i-1,j})/h$,
# with the convention $φ = 0$ outside the domain (and similarly in the vertical direction).
#
# 1. Show that the discrete divergence of the discrete gradient is exactly the five-point Laplacian of Exercise 1.
#    Deduce that if $φ$ solves $\Delta_h φ = \nabla_h \cdot \mathbf u$, then $\mathbf u - \nabla_h φ$ is exactly divergence-free.
#
# 2. Complete the functions `divergence(f)` and `project!(f, solve!)`,
#    where `solve!(φ, r, h)` is a function that approximately solves the Poisson equation in place.
#    Use `f.φ` as the initial guess (warm start).

# + nbgrader={"grade": false, "grade_id": "project", "locked": false, "schema_version": 3, "solution": true, "task": false}
function divergence(f)
    n, h, u, v = f.n, f.h, f.u, f.v
    ### BEGIN SOLUTION
    return [(u[i+1, j] - u[i, j] + v[i, j+1] - v[i, j]) / h for i in 1:n, j in 1:n]
    ### END SOLUTION
end

function project!(f, solve!)
    n, h, φ = f.n, f.h, f.φ
    solve!(φ, divergence(f), h)
    ### BEGIN SOLUTION
    for j in 1:n, i in 1:n+1
        f.u[i, j] -= (at(φ, i, j) - at(φ, i-1, j)) / h
    end
    for j in 1:n+1, i in 1:n
        f.v[i, j] -= (at(φ, i, j) - at(φ, i, j-1)) / h
    end
    ### END SOLUTION
    return f
end

# Test: after an accurate projection, the velocity field is divergence-free
f = Fluid(32)
f.u .= randn(33, 32); f.v .= randn(32, 33)
A₃₂ = lu(poisson_matrix(32))
project!(f, (φ, r, h) -> φ .= reshape(A₃₂ \ vec(r), 32, 32))
@assert norm(divergence(f)) < 1e-8
# -

# We can now put everything together.
# Smoke and an upward velocity are injected in a small disk at the bottom of the box,
# with a small horizontal velocity oscillating in time to break the symmetry.
# The smoke is lighter than the surrounding air, which creates a buoyancy force $\mathbf f = (0, β s)$.

# +
function add_sources!(f, Δt; center=(0.5, 0.1), radius=0.05, speed=1.0, β=4.0, wobble=0.3)
    n, h = f.n, f.h
    for j in 1:n, i in 1:n
        x, y = (i - 0.5)*h, (j - 0.5)*h
        if hypot(x - center[1], y - center[2]) < radius
            f.s[i, j] = 1
            f.v[i, j] = f.v[i, j+1] = speed
            f.u[i, j] = f.u[i+1, j] = wobble * sin(10f.t)
        end
    end
    # Buoyancy: average the smoke density on horizontal faces
    for j in 2:n, i in 1:n
        f.v[i, j] += Δt * β * (f.s[i, j-1] + f.s[i, j]) / 2
    end
    return f
end

function step!(f, Δt, solve!)
    add_sources!(f, Δt)
    # All three fields are advected with the same (old) velocity field
    u, v, s = advect(f.u, offset_u(f), f, Δt), advect(f.v, offset_v(f), f, Δt), advect(f.s, offset_s(f), f, Δt)
    f.u, f.v, f.s = u, v, s
    project!(f, solve!)
    f.t += Δt
    return f
end
# -

# 3. Simulate the smoke with $n = 96$, $Δt = 0.01$ and 300 time steps,
#    using 20 iterations of optimal SOR for each projection,
#    and create an animation of the smoke density.

# + nbgrader={"grade": false, "grade_id": "smoke", "locked": false, "schema_version": 3, "solution": true, "task": false}
### BEGIN SOLUTION
n = 96
cells = ((1:n) .- 0.5) / n
ω_opt = 2 / (1 + sin(π/(n+1)))
solve_sor!(φ, r, h) = sor!(φ, r, h, 20; ω=ω_opt)

f = Fluid(n)
anim = @animate for k in 1:300
    step!(f, 0.01, solve_sor!)
    show_field(f.s, axis=false, ticks=false, clims=(0, 1))
end every 2
gif(anim, "smoke.gif", fps=25)
### END SOLUTION
# -

# 4. Now compare four versions of the simulation, side by side:
#    no projection at all, 20 iterations of Jacobi, 20 iterations of optimal SOR, and an exact solve.
#    For each version, also plot the norm of the divergence of the velocity field after each time step.
#    Which solver would you choose for a video game?

# + nbgrader={"grade": false, "grade_id": "smoke-compare", "locked": false, "schema_version": 3, "solution": true, "task": false}
### BEGIN SOLUTION
A = lu(poisson_matrix(n))

solvers = [
    "No projection" => (φ, r, h) -> φ .= 0,
    "Jacobi (20 it.)" => (φ, r, h) -> jacobi!(φ, r, h, 20),
    "SOR (20 it.)" => solve_sor!,
    "Exact" => (φ, r, h) -> φ .= reshape(A \ vec(r), n, n),
]

fluids = [Fluid(n) for _ in solvers]
divergences = [Float64[] for _ in solvers]
anim = @animate for k in 1:300
    for (f, (name, solve!), d) in zip(fluids, solvers, divergences)
        step!(f, 0.01, solve!)
        push!(d, norm(divergence(f)) * f.h)
    end
    plot([show_field(f.s, axis=false, ticks=false, clims=(0, 1), title=name, titlefontsize=10)
          for (f, (name, _)) in zip(fluids, solvers)]..., layout=(1, 4), size=(1000, 280))
end every 3
display(gif(anim, "smoke_comparison.gif", fps=20))

plot(xlabel="Time step", ylabel=L"\|\nabla_h \cdot \mathbf{u}\|", yscale=:log10, legend=:outerright)
for ((name, _), d) in zip(solvers, divergences)
    plot!(max.(d, 1e-16), label=name)
end
plot!()
### END SOLUTION
# -

# ### <font color='orange'>[Exercise 4]</font> Going real-time in the browser
#
# Our Julia simulation is not interactive.
# In this last part, the goal is to port the simulator to JavaScript, so that it runs in real time in a web browser,
# and to let the user stir the smoke with the mouse.
# This is a good task to carry out **with the help of an AI assistant**, for instance as follows:
#
# 1. Give your Julia code to the assistant and ask it to translate `sample`, `advect`, `divergence`, `project!` and `sor!`
#    to JavaScript, using `Float64Array` buffers of size $(n+1) n$, $n (n+1)$ and $n^2$ instead of matrices.
#    **Check the translation carefully:** Julia is 1-based and column-major, and off-by-one errors are easy to make
#    (and also easy for an assistant to make!).
# 2. Ask for an HTML page with a `<canvas>` that runs one time step per animation frame (`requestAnimationFrame`)
#    and draws the smoke density, pixel by pixel, with `ImageData`.
# 3. Add mouse interaction: when the mouse is dragged, add smoke and set the velocity near the cursor
#    to the velocity of the mouse.
# 4. Add controls to choose the solver (none, Jacobi, Gauss–Seidel, SOR), the number of iterations per frame, and $ω$.
#    Add a view showing the divergence of the velocity field.
#
# Then play with it, and answer the following questions:
#
# - With 5 iterations per frame, how do Jacobi and SOR compare visually? What about 50 iterations?
# - What happens with SOR when $ω$ is close to 2? When $ω > 2$?
# - What is the largest grid size for which your simulation still runs at 60 frames per second?
#
# A reference implementation is provided in `smoke.html`: open it in your browser to see what you are aiming for,
# but try to write your own first!
#
# ### Going further
#
# - **Obstacles.** Add a solid obstacle (a disk) in the box. The velocity must vanish on the faces inside the obstacle,
#   and the Poisson equation must be modified accordingly (Neumann boundary conditions on the obstacle).
#   You should observe vortex shedding behind the obstacle.
# - **Vorticity confinement.** Numerical diffusion damps the small vortices. Look up *vorticity confinement*
#   (Fedkiw, Stam & Jensen, 2001) and implement it.
# - **Multigrid.** Jacobi and Gauss–Seidel quickly damp high-frequency errors but are slow on low frequencies.
#   The *multigrid* method addresses this by solving the problem on a hierarchy of coarser grids.
#   Implement a simple V-cycle and compare it with SOR.
# - **GPU.** Port the solver to WebGL or WebGPU shaders. Which of Jacobi and Gauss–Seidel is easier to parallelize? Look up *red-black Gauss–Seidel*.

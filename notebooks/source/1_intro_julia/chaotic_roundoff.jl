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
#     display_name: Julia 1.10.2
#     language: julia
#     name: julia-1.10
# ---

# ## Roundoff errors in a chaotic differential equation

using Plots

# The Lorenz system
# $$
# \begin{aligned}
# \dot x &= \sigma(y-x), \\
# \dot y &= x(\rho-z)-y, \\
# \dot z &= xy-\beta z,
# \end{aligned}
# $$
# is chaotic for the classical parameters
# $\sigma=10$, $\rho=28$ and $\beta=8/3$.
# Nearby trajectories separate exponentially, so tiny roundoff errors
# eventually produce macroscopically different numerical solutions.
#
# To isolate the effect of arithmetic precision, we use the same
# fixed-step fourth-order Runge--Kutta method, nominal time step and
# initial condition with `Float32`, `Float64` and 256-bit `BigFloat` arithmetic.

# +
lorenz(u) = begin
    T = eltype(u)
    σ = T(10)
    ρ = T(28)
    β = T(8)/T(3)
    (σ*(u[2] - u[1]),
     u[1]*(ρ - u[3]) - u[2],
     u[1]*u[2] - β*u[3])
end

function rk4_step(f, u, h)
    k₁ = f(u)
    k₂ = f(u .+ (h/2) .* k₁)
    k₃ = f(u .+ (h/2) .* k₂)
    k₄ = f(u .+ h .* k₃)
    return u .+ (h/6) .* (k₁ .+ 2 .* k₂ .+ 2 .* k₃ .+ k₄)
end

function lorenz_trajectory(::Type{T}; h=one(T)/1000, tmax=T(40)) where {T<:AbstractFloat}
    nsteps = round(Int, tmax/h)
    t = Vector{T}(undef, nsteps + 1)
    trajectory = Matrix{T}(undef, 3, nsteps + 1)
    u = (T(1), T(1), T(1))

    for n in 0:nsteps
        t[n+1] = n*h
        trajectory[1, n+1], trajectory[2, n+1], trajectory[3, n+1] = u
        n < nsteps && (u = rk4_step(lorenz, u, h))
    end
    return t, trajectory
end

# +
t32, u32 = lorenz_trajectory(Float32)
t64, u64 = lorenz_trajectory(Float64)
t_big, u_big = setprecision(256) do
    lorenz_trajectory(BigFloat)
end

# Compare all states in a common arithmetic before taking Euclidean distances.
distance32 = vec(sqrt.(sum(abs2, Float64.(u32) .- u64; dims=1)))
distance64 = setprecision(256) do
    Float64.(vec(sqrt.(sum(abs2, BigFloat.(u64) .- u_big; dims=1))))
end

function divergence_time(t, distance; threshold=1)
    index = findfirst(>(threshold), distance)
    return isnothing(index) ? Inf : t[index]
end

@show divergence_time(t64, distance32)
@show divergence_time(t64, distance64);

# +
p_trajectory = plot(t64, u64[1, :], label="Float64", xlabel="t", ylabel="x(t)")
plot!(p_trajectory, Float64.(t32), Float64.(u32[1, :]), label="Float32", alpha=.8)
plot!(p_trajectory, Float64.(t_big), Float64.(u_big[1, :]),
      label="BigFloat (256 bits)", alpha=.8)

p_separation = plot(t64[2:end], distance32[2:end], yscale=:log10,
                    label="Float32 versus Float64", xlabel="t",
                    ylabel="Euclidean separation")
plot!(p_separation, t64[2:end], distance64[2:end],
      label="Float64 versus BigFloat")
hline!(p_separation, [1], color=:black, linestyle=:dash, label="separation = 1")

plot(p_trajectory, p_separation, layout=(2, 1), size=(900, 700))
# -

# The trajectories initially overlap, but their separation grows by many
# orders of magnitude. `Float64` follows the high-precision discrete trajectory
# much longer than `Float32`, yet it eventually diverges too.
#
# This does **not** mean that the `BigFloat` curve is the exact solution of the
# differential equation: all three curves also have time-discretization error.
# Rather, the comparison isolates the roundoff accumulated by the *same discrete
# algorithm*. After divergence, pointwise prediction has lost meaning, although
# every trajectory may still reproduce statistical properties of the Lorenz attractor.

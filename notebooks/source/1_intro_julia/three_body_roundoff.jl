# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.19.5
#   kernelspec:
#     display_name: Julia 1.12
#     language: julia
#     name: julia-1.12
# ---

# ## Roundoff errors in the three-body problem

using Plots

# The gravitational three-body problem is a good illustration of sensitivity
# to roundoff, provided that the initial condition is chosen carefully: not
# every three-body orbit is chaotic. Here we use **Burrau's problem**, also
# called the Pythagorean three-body problem. Three bodies with masses 3, 4 and 5
# start from rest at the vertices of a 3--4--5 triangle:
# $$
# \begin{aligned}
# m_1 &= 3, & r_1(0) &= (1,3), \\
# m_2 &= 4, & r_2(0) &= (-2,-1), \\
# m_3 &= 5, & r_3(0) &= (1,-1).
# \end{aligned}
# $$
# Their center of mass is at the origin. With the gravitational constant set
# to one, their positions satisfy
# $$
# \ddot r_i = \sum_{j\ne i}m_j
# \frac{r_j-r_i}{\lVert r_j-r_i\rVert^3}.
# $$
# Close encounters strongly amplify small numerical perturbations.

# We integrate the same discrete problem with `Float16`, `Float32`, `Float64`
# and 256-bit `BigFloat` arithmetic. All four calculations use the same fourth-order
# Runge--Kutta method, nominal step size and initial data.

# +
function threebody_rhs(u, masses)
    T = eltype(u)
    du = similar(u)
    du[1:6] .= u[7:12]
    du[7:12] .= zero(T)

    for i in 1:3, j in 1:3
        i == j && continue
        position_i = 2i - 1
        position_j = 2j - 1
        velocity_i = 6 + position_i
        dx = u[position_j] - u[position_i]
        dy = u[position_j+1] - u[position_i+1]
        distance_squared = dx^2 + dy^2
        factor = masses[j] / (distance_squared * sqrt(distance_squared))
        du[velocity_i] += factor*dx
        du[velocity_i+1] += factor*dy
    end
    return du
end

function rk4_step(u, h, masses)
    k₁ = threebody_rhs(u, masses)
    k₂ = threebody_rhs(u .+ (h/2) .* k₁, masses)
    k₃ = threebody_rhs(u .+ (h/2) .* k₂, masses)
    k₄ = threebody_rhs(u .+ h .* k₃, masses)
    return u .+ (h/6) .* (k₁ .+ 2 .* k₂ .+ 2 .* k₃ .+ k₄)
end

function threebody_trajectory(::Type{T}; step_denominator=512,
                              tmax=1.88, save_every=1) where {T<:AbstractFloat}
    h = one(T)/T(step_denominator)
    masses = T[3, 4, 5]
    u = T[1, 3, -2, -1, 1, -1, 0, 0, 0, 0, 0, 0]
    nsteps = round(Int, tmax*step_denominator)
    nsaved = div(nsteps, save_every) + 1
    t = Vector{Float64}(undef, nsaved)
    trajectory = Matrix{T}(undef, 12, nsaved)
    t[1] = 0
    trajectory[:, 1] = u

    saved = 1
    for step in 1:nsteps
        u = rk4_step(u, h, masses)
        if step % save_every == 0
            saved += 1
            t[saved] = step/step_denominator
            trajectory[:, saved] = u
        end
    end
    return t, trajectory
end

# +
t16, u16 = threebody_trajectory(Float16)
t32, u32 = threebody_trajectory(Float32)
t64, u64 = threebody_trajectory(Float64)
t_big, u_big = setprecision(256) do
    threebody_trajectory(BigFloat)
end

# Compare the six position coordinates directly with the `BigFloat` reference.
function position_error(trajectory, reference)
    setprecision(256) do
        Float64.(vec(sqrt.(sum(abs2,
            BigFloat.(trajectory[1:6, :]) .- reference[1:6, :]; dims=1))))
    end
end

error16 = position_error(u16, u_big)
error32 = position_error(u32, u_big)
error64 = position_error(u64, u_big)

function divergence_time(t, distance; threshold=0.1)
    index = findfirst(>(threshold), distance)
    return isnothing(index) ? Inf : t[index]
end

@show divergence_time(t64, error16)
@show divergence_time(t64, error32)
@show divergence_time(t64, error64)
@show error16[end]
@show error32[end]
@show error64[end];

# +
p_distance = plot(t64[2:end], error16[2:end], yscale=:log10,
                  label="Float16", xlabel="t",
                  ylabel="Position error against BigFloat")
plot!(p_distance, t64[2:end], error32[2:end], label="Float32")
plot!(p_distance, t64[2:end], error64[2:end], label="Float64")
hline!(p_distance, [0.1], color=:black, linestyle=:dash, label="error = 0.1")
p_distance

# The following cell renders a synchronized animation. All three trajectory
# panels use the same fixed, tightly zoomed camera.

# +
plot_trajectories = (Float64.(u16), Float64.(u32), u64)
format_names = ("Float16", "Float32", "Float64")
reference_trajectory = Float64.(u_big)
body_colors = (:dodgerblue, :darkorange, :seagreen)
marker_sizes = (6, 8, 10)

function trajectory_panel(trajectory, frame, format_name; trail=120)
    first_frame = max(1, frame - trail)
    p = plot(xlims=(-2.5, 2.5), ylims=(-1.5, 3.5), aspect_ratio=:equal,
             xlabel="x", ylabel="y", title=format_name,
             legend=false, framestyle=:box)
    for body in 1:3
        xindex = 2body - 1
        plot!(p, reference_trajectory[xindex, first_frame:frame],
             reference_trajectory[xindex+1, first_frame:frame],
             color=:gray, linewidth=1, linestyle=:dash)
        plot!(p, trajectory[xindex, first_frame:frame],
             trajectory[xindex+1, first_frame:frame],
             color=body_colors[body], linewidth=2)
        scatter!(p, [trajectory[xindex, frame]], [trajectory[xindex+1, frame]],
                 color=body_colors[body], markersize=marker_sizes[body])
    end
    return p
end

function error_panel(frame)
    indices = frame < 2 ? (1:0) : (2:frame)
    p = plot(t64[indices], error16[indices], yscale=:log10,
             ylims=(1e-16, 1), label="Float16", xlabel="t",
             ylabel="Position error against BigFloat", legend=:topleft)
    plot!(p, t64[indices], error32[indices], label="Float32")
    plot!(p, t64[indices], error64[indices], label="Float64")
    hline!(p, [0.1], color=:black, linestyle=:dash, label="error = 0.1")
    return p
end

frame_indices = unique(round.(Int, LinRange(1, length(t64), 240)))
animation = @animate for frame in frame_indices
    panels = map(zip(plot_trajectories, format_names)) do (trajectory, format_name)
        trajectory_panel(trajectory, frame, format_name)
    end
    plot(panels..., error_panel(frame), layout=(2, 2), size=(1000, 850),
         plot_title="Burrau three-body problem, t = $(round(t64[frame], digits=2))")
end
gif(animation, "three_body_precision.gif", fps=24)

# ![Comparison of Float16, Float32 and Float64 trajectories](three_body_precision.gif)
#
# The step size $h=1/512$ is exactly representable in all three binary formats
# and is large enough for every `Float16` position coordinate to change. The
# `Float16` trajectory visibly drifts from the reference, while `Float32` and
# `Float64` remain much closer over this interval.
# At $t=1.88$, their respective position errors are approximately
# $2.73\times 10^{-1}$, $5.31\times 10^{-5}$ and $2.08\times 10^{-13}$.
# We stop just before the first near-singular close encounter, which a fixed-step
# method cannot resolve smoothly at these precisions.
#
# The high-precision curve is still not the exact solution of the continuous
# equations: all three trajectories have time-discretization error. The
# comparison instead approximates the roundoff error of a fixed discrete
# algorithm. A responsible numerical study should also repeat the calculation
# with a smaller time step and monitor conserved quantities such as total energy.

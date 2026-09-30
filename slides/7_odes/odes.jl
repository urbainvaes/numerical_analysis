using LinearAlgebra
using Plots
using DifferentialEquations
using LaTeXStrings
Plots.default(label=nothing,
              fontfamily="Computer Modern",
              background_color=colorant"#f5ede9",
              # background_color=colorant"#f5ede9",
              foreground_color=:black,
              linewidth=2,
              xlabelfontsize=18,
              ylabelfontsize=18,
              legendfontsize=18,
              xtickfontsize=18,
              ytickfontsize=18);


N = 10
T = 10.
ts = LinRange(0, T, N)
x₀ = 0
xlims=(ts[j], T + 1)
Δ = ts[2]
xs = zero(ts)
xs[1] = x₀

α = .15
f(t, x) = α * (x - sin(t)) + cos(t)
sol(t₀, x₀) = t -> sin(t) + exp(α*(t - t₀)) * (x₀ - sin(t₀))
ylims = (-1.2, 3)

anim = Animation()
for i in 1:length(ts) - 1
    plot(size=(900, 600), xlabel=L"t", bottom_margin=5Plots.mm)
    plot!(sol(ts[1], xs[1]), xlims=(0, T + 1), color=:green)
    for j in 1:i-1
        plot!(sol(ts[j], xs[j]), xlims=(ts[j], T + 1), color=:lightgray)
    end
    scatter!(ts[1:i], xs[1:i], color=:red)
    plot!(ts[1:i], xs[1:i], color=:blue)
    plot!(sol(ts[i], xs[i]), xlims=(ts[i], T + 1), color=:gray)
    plot!(xlims=(0, T), ylims=ylims)
    frame(anim)
    xs[i + 1] = xs[i] + Δ*f(ts[i], xs[i])
    scatter!(ts[1:i+1], xs[1:i+1], color=:red)
    plot!(ts[1:i+1], xs[1:i+1], color=:blue)
    frame(anim)
end
gif(anim, "forward_euler.gif", fps=N/10)


α = -1
x0s = -2:.2:2
plot(size=(900, 600))
for x0 in x0s
    plot!(sol(0, x0), xlims=(0, 10), color=:gray)
end
annotate!(5, 1.5, text(L"\alpha = " * string(α), :left))
plot!(sol(0, 0), xlims=(0, 10), color=:green, label=L"\sin(t)")
plot!(xlabel=L"t")

α = 1
x0s = LinRange(-.0001, .0001, 20)
plot(size=(900, 600))
for x0 in x0s
    plot!(sol(0, x0), xlims=(0, 10), color=:gray)
end
annotate!(5, 1.5, text(L"\alpha = " * string(α), :left))
plot!(sol(0, 0), xlims=(0, 10), color=:green, label=L"\sin(t)")
plot!(xlabel=L"t")


# f(t, x) = α * x
# sol(t₀, x₀) = t -> x₀ * exp(α * (t - t₀))
# x₀ = 1
# ylims = (.5, 5)
# N = 4

α = .15
function forward_euler(f, x₀, Δ, T)
    ts = 0:Δ:T
    xs = zero(ts)
    xs[1] = x₀
    for i in 1:length(ts) - 1
        xs[i + 1] = xs[i] + Δ*f(ts[i], xs[i])
    end
    return ts, xs
end

T = 10
Δs = 2.0.^(0:-1:-5)
anim = @animate for Δ in Δs
    ts, xs = forward_euler(f, x₀, Δ, T)
    plot(size=(900, 600), xlabel=L"t", bottom_margin=5Plots.mm)
    plot!(sol(0, x₀), xlims=(0, T), ylims=ylims, color=:green)
    scatter!(ts, xs, color=:red)
    plot!(ts, xs, color=:blue)
    annotate!(.5, 2.5, text(L"\Delta = " * string(Δ), :left))
end
gif(anim, "convergence.gif", fps=1)

α = .15
Δs = 2.0.^(0:-1:-10)
errors = zero(Δs)
for (i, Δ) in enumerate(Δs)
    ts, xs = forward_euler(f, x₀, Δ, T)
    xexact = sin
    errors[i] = maximum(abs.(xs - xexact.(ts)))
end
plot(Δs, errors, seriestype=[:line,:scatter],
     xscale=:log, yscale=:log, bg=:transparent,
     size=(900, 600), xlabel=L"\Delta", title=L"\max_n|x(t_n) - x_n|",
     bottom_margin=5Plots.mm)

using Polynomials
fit(log.(Δs), log.(errors),1)

xlims = (-6, 2)
ylims = (-4, 4)
Plots.plot(size=(900, 600))
x = LinRange(xlims..., 400)
y = LinRange(ylims..., 400) 
is_stable(z, p) = abs(sum(z^i/factorial(i) for i in 0:p)) < 1
stable = sum(is_stable.(x' .+ im*y, p) for p in 1:4)
Plots.contourf!(x, y, stable, c=:Pastel1_3,
                aspect_ratio=:equal, levels=[0, 1, 2, 3, 4, 5],
                xlabel=L"\Re(\Delta \lambda)",
                ylabel=L"\Im(\Delta \lambda)",
                colorbar=:none, xlims=xlims, ylims=ylims)
Plots.contour!(x, y, stable, levels=5, c=:green)
Plots.vline!([0], color=:gray) 
Plots.hline!([0], color=:gray) 
Plots.plot!(bg=:transparent)



α = -10
T = 2
Δ = .18
x₀ = 1.
f(t, x) = α * (x - sin(t)) + cos(t)
sol(t₀, x₀) = t -> sin(t) + exp(α*(t - t₀)) * (x₀ - sin(t₀))
ts, xs = forward_euler(f, x₀, Δ, T)
plot(ts, xs, seriestype=[:line,:scatter])
for (t, x) in zip(ts, xs)
    plot!(sol(t, x), xlims=(t, T), c=:lightgray)
end
plot!(sol(0, x₀), xlims=(0, T), c=:green)

xlims, ylims = (-5, 1), (-2, 2)
x = LinRange(xlims..., 100)
y = LinRange(ylims..., 100)
stable = @. clamp(abs(1 + x' + im*y), 0, 2)
Plots.contourf(x, y, stable, c=:Pastel1_3,
               aspect_ratio=:equal, levels=[0., 1., 2],
               xlabel=L"\Re(\Delta \lambda)",
               ylabel=L"\Im(\Delta \lambda)",
               xlims=xlims, ylims=ylims, bg=:transparent)
Plots.contour!(x, y, stable, color=:green, 
               levels=[0., 1., 1e5])

Plots.vline!([0], color=:gray)
Plots.hline!([0], color=:gray)

xlims, ylims = (-1, 5), (-2, 2)
x = LinRange(xlims..., 200)
y = LinRange(ylims..., 200)
stable = @. clamp(abs(1 / (1 - x' - im*y)), 0, 2)
Plots.contourf(x, y, stable, c=:Pastel1_3,
               aspect_ratio=:equal, levels=[0., 1., 2.],
               xlabel=L"\Re(\Delta \lambda)",
               ylabel=L"\Im(\Delta \lambda)",
               xlims=xlims, ylims=ylims, bg=:transparent)
Plots.contour!(x, y, stable, levels=[0., 1., 2.], color=:green)
Plots.vline!([0], color=:gray)
Plots.hline!([0], color=:gray)

λ = .3
expf(t, x) = - λ*x + .35
expsol(t) = exp(- λ*t)
ts = [0, 1, 2]
xs = expsol.(ts)
fs = expf.(ts, xs)
using Polynomials
p = fit(ts, fs)
plot(grid=false, showaxis=false, size=(900, 600), bg=:transparent) 
plot!(p, xlims=(0, 3.5), label=L"Quadratic interpolant $\widehat f$")
scatter!(ts, expf.(ts, xs), xlims=(-2.5, 4))
annotate!(ts[1] - 1.7, fs[1], text(L"f(t_{n-2}, x_{n-2})", :left, 20))
annotate!(ts[2] - 1.7, fs[2], text(L"f(t_{n-1}, x_{n-1})", :left, 20))
annotate!(ts[3] - 1.2, fs[3], text(L"f(t_n, x_{n})", :left, 20))
newts = LinRange(2, 3, 20)
plot!(newts, zero(newts), fillrange=p.(newts), fillalpha=.5)
annotate!(2, .05, text(L"\int_{t_n}^{t_{n+1}} \widehat f(t) \, \rm d t", :left, 18))
hline!([0], color=:gray, legend=:outertop)
plot!([ts[1], ts[1]], [0, fs[1]], color=:lightgray, ls=:dash)
plot!([ts[2], ts[2]], [0, fs[2]], color=:lightgray, ls=:dash)
plot!([ts[3], ts[3]], [0, fs[3]], color=:lightgray, ls=:dash)
annotate!(ts[1], -.01, text(L"t_{n-2}", :left, 20))
annotate!(ts[2], -.01, text(L"t_{n-1}", :left, 20))
annotate!(ts[3], -.01, text(L"t_{n}", :left, 20))
annotate!(3, -.01, text(L"t_{n+1}", :left, 20))

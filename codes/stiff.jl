import DifferentialEquations as ode
import Plots
using LaTeXStrings

Plots.default(fontfamily="Computer Modern",
              titlefontsize=14,
              xlabelfontsize=14,
              ylabelfontsize=14,
              legendfontsize=12,
              xtickfontsize=12,
              ytickfontsize=12,
              framestyle=:box,
              label=nothing,
              grid=true)

α = 10
function stiff_equation!(dx, x, p, t)
    global α
    dx[1] = - α * (x[1] - sin(t)) + cos(t)
end

T = 10
x0, tspan = [1], (0, T)

Plots.plot(title=L"Solutions to a stiff differential equation for $\alpha = %$α$")
for x0 in [exp.(-3:10); -exp.(-3:10)]
    stiff_problem = ode.ODEProblem(stiff_equation!, [x0], tspan)
    sol = ode.solve(stiff_problem)
    Plots.plot!(sol, color=:gray)
end
Plots.plot!(legend=:none, ylims=(-2, 2))
Plots.plot!(t -> sin(t), color=:green, linewidth=2, xlabel=L"t")
Plots.savefig("ode_stiff_$α.pdf")

α = 10
stiff_problem = ode.ODEProblem(stiff_equation!, [1], tspan)
sol = ode.solve(stiff_problem)

T = 10
f(x, t) = - α * (x - sin(t)) + cos(t)
Δ_critical = 2/α
factor = 2
Δ = factor * Δ_critical
ts = 0:Δ:T
n = length(ts) - 1
xs = zeros(n+1)
xs[1] = 1
# function forward_euler(Δ, ts)
for i in 1:n
    # xs[i+1] = xs[i] + Δ * f(xs[i], ts[i])
    xs[i+1] = (xs[i] + α*Δ*sin(ts[i+1]) + Δ*cos(ts[i+1])) / (1 + α*Δ)
end

Plots.plot(ts, xs, linewidth=2)
Plots.scatter!(ts, xs, color=:black, markersize=2)
Plots.plot!(sol, color=:green, linewidth=2, xlabel=L"t", legend=:none)
Plots.plot!(title=L"Backward Euler method with $\Delta = %$factor \Delta_*$")
# Plots.savefig("forward_euler_$factor.pdf")
Plots.savefig("backward_euler_$factor.pdf")

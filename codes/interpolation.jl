import Polynomials
import Plots
import LaTeXStrings

Plots.default(fontfamily="Computer Modern",
              titlefontsize=20,
              xlabelfontsize=20,
              ylabelfontsize=20,
              legendfontsize=12,
              xtickfontsize=12,
              ytickfontsize=12,
              linewidth=2,
              markersize=7,
              framestyle=:box,
              label=nothing,
              grid=false)

lims = (-1, 1)

n = 30
k = 0:(n-1)
x = cos.(π*(k.+1/2)/n) |> sort
# x = LinRange(0, 2π, n)
# x = LinRange(lims[1], lims[2], n)

# ε = .01
f(x) = sin(x)
f(x) = 1 / (1 + 25x^2)
y = f.(x)

# Interpolation
p = Polynomials.fit(x, y)

# Grid for the plots
# x_plot = LinRange(0, 2π, 100)
x_plot = LinRange(lims[1], lims[2], 600)
Plots.plot(x_plot, f.(x_plot), xlims=lims, size=(900,450))
Plots.plot!(x_plot, p.(x_plot))
Plots.scatter!(x, f.(x), legend=:bottomright)
Plots.xlims!(lims[1], lims[2])
# Plots.xlims!(lims[1]*1.001, lims[2]*1.001)
Plots.savefig("interpolation_sine_$n.pdf")

# Plots.plot!(title="Interpolation of the sine function with $n nodes")

n = 15
k = 0:(n-1)
x = cos.(π*(k.+1/2)/n) |> sort
# x = LinRange(-1, 1, n)

# ε = .01
f(x) = sin(x)
f(x) = 1 / (1 + 25x^2)
y = f.(x)

# Interpolation
p = Polynomials.fit(x, y)

# Grid for the plots
x_plot = LinRange(-1, 1, 500)
Plots.plot(x_plot, f.(x_plot), xlims=(-1, 1), size=size=(900,450))
Plots.plot!(x_plot, p.(x_plot))
Plots.scatter!(x, f.(x), legend=:bottomright)
# Plots.plot!(title="Interpolation of the Runge function with $n nodes")
Plots.savefig("interpolation_cheb_runge_$n.pdf")


ε = .01
f(x) = tanh((x+1/2)/ε) + tanh(x/ε) + tanh((x-1/2)/ε)
# f(x) = x^(n)
y = f.(x)

# Interpolation
p = Polynomials.fit(x, y)

# Grid for the plots
x_plot = range(-1, 1, length=200)

using Symbolics
@variables x

runge = 1 / (1 + 25x^2)

n = 40
a, b = 0, 4π
θ = LinRange(a, b, n)
# θ = a .+ (b-a)*(cos.(π*(0:n .+ 1/2)/(n+1)) .+ 1)/2

# ε = .01
R, r, d = 5, 2, 3
f1(θ) = (R+r)*cos(θ) + d*cos((R+r)/r*θ)
f2(θ) = (R+r)*sin(θ) - d*sin((R+r)/r*θ)
p1 = Polynomials.fit(θ, BigFloat.(f1.(θ)))
p2 = Polynomials.fit(θ, BigFloat.(f2.(θ)))

# Grid for the plots
θ_plot = LinRange(a, b, 500)
Plots.plot(f1.(θ_plot), f2.(θ_plot), label="Exact function", axis=nothing, xlims=(-10, 25))
Plots.plot!(p1.(θ_plot), p2.(θ_plot), size=size=(500,450), aspect_ratio=:equal, label="Polynomial interpolation")
Plots.scatter!(p1.(θ), p2.(θ), size=(900,450), label="Interpolation nodes")
Plots.savefig("interpolation.pdf")

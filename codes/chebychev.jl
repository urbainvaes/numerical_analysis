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
              framestyle=:box,
              label=nothing,
              grid=true)

# For information on any of these, type
# Plots.plotattr("framestyle")

n = 20
k = 0:(n-1)
x = cos.(π*(k.+1/2)/n) |> sort

ε = .01
f(x) = tanh((x+1/2)/ε) + tanh(x/ε) + tanh((x-1/2)/ε)
# f(x) = x^(n)
y = f.(x)

# Interpolation
p = Polynomials.fit(x, y)

# Grid for the plots
x_plot = range(-1, 1, length=200)

# Figure for cover page
Plots.plot(x_plot, f.(x_plot), linewidth=4, ticks=false)
Plots.plot!(x_plot, p.(x_plot), linewidth=4)
Plots.scatter!(x, y, markersize=4)
Plots.xlims!(-1, 1)
Plots.savefig("figures/chebychev_cover.pdf")

# Full figure
Plots.plot(x_plot, f.(x_plot), linewidth=4, label="Function", ticks=true)
Plots.plot!(x_plot, p.(x_plot), linewidth=4, label="Interpolating polynomial")
Plots.scatter!(x, y, markersize=4, legend=:bottomright)
Plots.title!("Chebychev interpolation")
Plots.xlims!(-1, 1)
Plots.savefig("figures/chebychev.pdf")

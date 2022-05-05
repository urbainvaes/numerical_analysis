# +
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

n = 10
k = 0:(n-1)
x = cos.(π*(k.+1/2)/n) |> sort
# x = LinRange(-1, 1, n)

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
# Plots.savefig("figures/chebychev_cover.pdf")

anim = Plots.Animation()
for n in 2:50
    k = 0:(n-1)
    x = cos.(π*(k.+1/2)/n) |> sort
    # x = LinRange(-1, 1, n)
    ε = .01
    f(x) = BigFloat.(tanh((x+1/2)/ε) + tanh(x/ε) + tanh((x-1/2)/ε))
    # f(x) = x^(n)
    y = f.(x)
    # Interpolation
    p = Polynomials.fit(x, y)
    # Grid for the plots
    x_plot = range(-1, 1, length=200)
    # Figure for cover page
    Plots.plot(x_plot, f.(x_plot), linewidth=4, ticks=false, size=(900, 400))
    Plots.plot!(x_plot, p.(x_plot), linewidth=4)
    Plots.scatter!(x, y, markersize=4)
    Plots.xlims!(-1, 1)
    Plots.ylims!(-4, 4)
    Plots.frame(anim)
end
Plots.gif(anim, "animation.webm", fps=1)

# Full figure
Plots.plot(x_plot, f.(x_plot), linewidth=4, label="Function", ticks=true)
Plots.plot!(x_plot, p.(x_plot), linewidth=4, label="Interpolating polynomial")
Plots.scatter!(x, y, markersize=4, legend=:bottomright)
Plots.title!("Chebychev interpolation")
Plots.xlims!(-1, 1)
# Plots.savefig("figures/chebychev.pdf")

anim = Plots.Animation()
series = decompose(x -> exp(log(exact0(x)) + x^2/2), dmax, quad)
dmax = 100
for d in 0:dmax
    approx_h = eval(series[1:d+1], xplot)
    Plots.plot(xplot, exact0.(xplot), size=(1200,800))
    Plots.plot!(xplot, approx_h .* factor.(xplot))
    Plots.ylims!((-.1, 1))
    Plots.title!("d=$d")
    Plots.frame(anim)
end
Plots.gif(anim, "animation.gif", fps=4)

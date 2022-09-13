import Polynomials
import Plots
using LaTeXStrings

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


f(x) = exp(x)
d(δ) = (f(1+δ) - f(1))/δ
δs = 10 .^(0:-.1:-17)
err = abs.(d.(δs) .- exp(1))
Plots.plot(δs, err, xscale=:log10, yscale=:log10, label=L"|d(\delta) - f'(1)|", size=(900,450))
Plots.plot!([eps()], seriestype=:vline, label = L"\varepsilon")
Plots.plot!([sqrt(eps())], seriestype=:vline, label = L"\sqrt{\varepsilon}")
Plots.savefig("figures/numerical_differentiation.pdf")

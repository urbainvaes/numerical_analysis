import Plots
using LaTeXStrings
Plots.resetfontsizes()
Plots.scalefontsizes(1.5)
Plots.default(titlefont = ("computer modern"), 
              legendfont = ("computer modern"))

eps_64 = eps(Float64)
x_64 = Float64.(2.0.^(-4:.01:4))
successors_64 = nextfloat.(x_64)
distance_64 = successors_64 - x_64

eps_32 = eps(Float32)
x_32 = Float32.(2.0.^(-4:.01:4))
successors_32 = nextfloat.(x_32)
distance_32 = successors_32 - x_32

Plots.plot(x_64, 1 ./ distance_64, label=L"\frac{1}{\Delta(x)}~\mathrm{(Float64)}")
Plots.plot!(xaxis=:log, yaxis=:log, top_margin=2Plots.mm, bottom_margin=2Plots.mm)
Plots.xticks!(2.0.^(-4:4), [latexstring("2^{"* string(i) * "}") for i in -4:4])
Plots.title!("Density of Float64 numbers", fontsize=12)
Plots.xlabel!(L"x")
Plots.savefig("figures/float64_density.pdf")

Plots.plot(x_64, (distance_64./x_64) / eps_64, xaxis=:log, label=L"\frac{\Delta(x)}{x \, \varepsilon_M}")
Plots.xticks!(2.0.^(-4:4), [latexstring("2^{"* string(i) * "}") for i in -4:4])
Plots.xlabel!(L"x")
Plots.title!("Relative spacing of Float64 numbers", fontsize=12)
Plots.savefig("figures/float64_spacing.pdf")

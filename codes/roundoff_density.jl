import Plots
using LaTeXStrings
Plots.resetfontsizes()
Plots.scalefontsizes(1.5)
Plots.default(titlefont = "computer modern",
              legendfont = "computer modern",
              linewidth = 2,
              top_margin=20Plots.mm, 
              bottom_margin=20Plots.mm)

eps_64 = eps(Float64)
exp_min, exp_max = -4, 4
exponents = exp_min:.01:exp_max
x_64 = Float64.(2.0.^exponents)
successors_64 = nextfloat.(x_64)
distance_64 = successors_64 - x_64

eps_32 = eps(Float32)
x_32 = Float32.(2.0.^exponents)
successors_32 = nextfloat.(x_32)
distance_32 = successors_32 - x_32

Plots.plot(x_64, distance_64, label=L"\Delta(x)~\mathrm{(Float64)}")
Plots.plot!(xaxis=:log, yaxis=:log, size=(900,450))
Plots.xticks!(2.0.^(exp_min:exp_max), [latexstring("2^{"* string(i) * "}") for i in exp_min:exp_max])
Plots.title!("Absolute spacing between Float64 numbers")
Plots.xlabel!(L"x")
Plots.savefig("figures/float64_density.pdf")

Plots.plot(x_64, (distance_64./x_64) / eps_64, xaxis=:log, label=L"\frac{\Delta(x)}{x}",
           bottom_margin=2Plots.mm, size=(900,450))
Plots.xticks!(2.0.^(-4:4), [latexstring("2^{"* string(i) * "}") for i in exp_min:exp_max])
yticks = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
Plots.yticks!(yticks, [latexstring(string(i) * "\\, \\varepsilon_M") for i in yticks])
Plots.xlabel!(L"x")
Plots.title!("Relative spacing between Float64 numbers")
Plots.savefig("figures/float64_spacing.pdf")

eps_64 = eps(Float64)
exp_min, exp_max = -1024, -1016
exponents = exp_min:.01:exp_max
x_64 = Float64.(2.0.^exponents)
successors_64 = nextfloat.(x_64)
distance_64 = successors_64 - x_64

Plots.plot(x_64, (distance_64./x_64) / eps_64, xaxis=:log, label=L"\frac{\Delta(x)}{x}",
           bottom_margin=2Plots.mm, size=(900,450))
Plots.xticks!(2.0.^(exp_min:exp_max), [latexstring("2^{"* string(i) * "}") for i in exp_min:exp_max])
Plots.plot!([2^(-1022)], seriestype=:vline, label = false)
yticks = [1, 2, 3, 4]
Plots.yticks!(yticks, [latexstring(string(i) * "\\, \\varepsilon_M") for i in yticks])
Plots.xlabel!(L"x")
Plots.title!("Relative spacing between Float64 numbers")
Plots.savefig("figures/float64_spacing_denormalized.pdf")

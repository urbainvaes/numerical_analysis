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

# name = "backward Euler"
# fun(z) = abs(1 / (1 - z))
# xlims = (-1, 5)
# ylims = (-2, 2)

# name = "Crank-Nicolson"
# fun(z) = abs((1 + z/2)/(1 - z/2))
# xlims = (-3, 3)
# ylims = (-2, 2)

name = "forward Euler"
xlims = (-5, 1)
ylims = (-2, 2)

name = "Taylor2"
xlims = (-6, 2)
ylims = (-4, 4)

name = "Taylor3"
xlims = (-6, 2)
ylims = (-4, 4)

name = "Taylor4"
xlims = (-6, 2)
ylims = (-3.5, 3.5)

fun = Dict(
    "forward Euler" => z -> abs(1 + z),
    "Taylor2" => z -> abs(1 + z + z^2/2),
    "Taylor3" => z -> abs(1 + z + z^2/2 + z^3/6),
    "Taylor4" => z -> abs(1 + z + z^2/2 + z^3/6 + z^4/24))

x = LinRange(xlims..., 100)
y = LinRange(ylims..., 100) 
z = x' .+ im*y
fun_plot = fun["Taylor4"].(z)

Plots.contourf(x, y, fun_plot, 
               title="Stability region for the $name method",
               color=Plots.cgrad(:Pastel1_3, rev=true),
               aspect_ratio=:equal,
               levels=[-10_000, 1, 10_000],
               xlabel=L"\Re(\Delta \lambda)",
               ylabel=L"\Im(\Delta \lambda)",
               colorbar=:none)

for key in keys(fun)
    fun_plot = fun[key].(z)
    Plots.contour!(x, y, fun_plot, color=:green, levels=[1])
end

Plots.vline!([0], color=:gray) 
Plots.hline!([0], color=:gray) 
Plots.title!("Stability regions for the Taylor methods") 
Plots.plot!(xlims=xlims, ylims=ylims)
Plots.savefig("stability_$name.pdf")

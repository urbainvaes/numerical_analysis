import Plots

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

np = 20
xs = LinRange(0, 1, np)

x_interp = LinRange(0, 1, 200)
result = zeros(length(x_interp))

pl = Plots.plot(size=(900,450), bottom_margin=4Plots.mm)
for d in 1:length(xs)
    temp = xs[[1:d-1; d+1:length(xs)]]
    p(x) = prod(x .- temp) / prod(xs[d] .- temp)
    result += p.(x_interp)
    Plots.plot!(p, 0, 1, label=nothing)
end

Plots.scatter!(xs, 0*xs .+ 1, label=nothing)
Plots.plot!(xlims=(0, 1), xlabel="x", title="Lagrange polynomials for $(np+1) equidistant nodes")
Plots.savefig("lagrange_$(np).pdf")
display(pl)

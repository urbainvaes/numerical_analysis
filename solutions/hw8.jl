import RDatasets
import Polynomials
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

data = RDatasets.dataset("datasets", "pressure")
temperature, pressure = data[:, 1], data[:, 2] 

poly = Polynomials.fit(temperature, log.(pressure), 4)
Plots.scatter(temperature, pressure, label="Data", size=(900,450), legend=:topleft,
              bottom_margin=5Plots.mm, left_margin=5Plots.mm)
xplot = LinRange(temperature[1], temperature[end], 100)
Plots.plot!(xplot, exp.(poly.(xplot)), label="Model")
Plots.xlabel!("Temperature")
Plots.ylabel!("Pressure")
Plots.savefig("pressure_model.pdf")

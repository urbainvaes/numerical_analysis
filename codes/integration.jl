using LaTeXStrings
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

import Plots
Plots.default(fontfamily="Computer Modern",
              titlefontsize=20,
              xlabelfontsize=20,
              ylabelfontsize=20,
              legendfontsize=12,
              xtickfontsize=12,
              ytickfontsize=12,
              framestyle=:box,
              linewidth=3,
              label=nothing,
              grid=false)

# Composite Simpson's rule
function composite_simpson(u, a, b, n)
    # Integration nodes
    x = LinRange(a, b, n + 1)
    # Evaluation of u at the nodes
    ux = u.(x)
    # Step size
    h = x[2] - x[1]
    # Approximation of the integral
    return (h/3) * sum([ux[1]; ux[end]; 4ux[2:2:end-1]; 2ux[3:2:end-2]])
end
# Function to integrate
u(x) = cos(x)
# Integration bounds
a, b = 0, π/2
# Exact integral
I = 1.0
# Number of subintervals
ns = [8; 16; 32; 64; 128]
# Approximate integrals
Î = composite_simpson.(u, a, b, ns)
# Calculate exact and approximate errors
for i in 2:length(ns)
    println("Exact error: $(I - Î[i]), ",
            "Approx error: $((Î[i] - Î[i-1])/15)")
end

J(h) = 1 + sin(h)
J₁(h) = (2J(h) - J(2h))
J₂(h) = (4J₁(h) - J₁(2h))/3
J₃(h) = (8J₂(h) - J₂(2h))/7
Plots.plot(bottom_margin=4Plots.mm, right_margin=4Plots.mm,
           xlims=(0, .2), xlabel="h", size=(900,450), xaxis=:log, yaxis=:log)
Plots.plot!(h -> abs(J(h) - 1), 0, .2, label=L"J(h)")
Plots.plot!(h -> abs(J₁(h) - 1), 0, .2, label=L"J_1(h)")
Plots.plot!(h -> abs(J₂(h) - 1), 0, .2, label=L"J_2(h)")
Plots.plot!(h -> abs(J₃(h) - 1), 0, .2, label=L"J_3(h)")
Plots.xlims!(10^(-10), 10)
Plots.ylims!(10^(-10), 10)
Plots.savefig("figures/richardson.pdf")


function composite_trapezoidal(u, a, b, h)
    # Number of intervals
    n = Int((b-a)/h)
    # Integration nodes
    x = LinRange(a, b, n + 1)
    # Evaluation of u at the nodes
    ux = u.(x)
    # Approximation of the integral
    return (h/2) * sum([ux[1]; ux[end]; 2ux[2:end-1]])
end
a, b = 0, π/2
J(h) = composite_trapezoidal(u, a, b, h)
J₁(h) = (4J(h) - J(2h))/3
J₂(h) = (16J₁(h) - J₁(2h))/15
J₃(h) = (64J₂(h) - J₂(2h))/63

hs = (b-a) * 2.0.^(-8:-3)
colors = Plots.palette(:default)
Plots.plot(hs, .1hs.^2, color=colors[1])
Plots.scatter!(hs, abs.(1 .- J.(hs)), color=colors[1], label=L"|I - J(h)|~~~O(h^2)")
Plots.plot!(hs, .005hs.^4, color=colors[2])
Plots.scatter!(hs, abs.(1 .- J₁.(hs)), color=colors[2], label=L"|I - J_1(h)|~~~O(h^4)")
Plots.plot!(hs, .0025hs.^6, color=colors[3])
Plots.scatter!(hs, abs.(1 .- J₂.(hs)), color=colors[3], label=L"|I - J_2(h)|~~~O(h^6)")
Plots.plot!(hs, .003hs.^8, color=colors[4])
Plots.scatter!(hs, abs.(1 .-J₃.(hs)), color=colors[4], label=L"|I - J_3(h)|~~~O(h^8)")
Plots.plot!(xlabel="h", xaxis=:log, yaxis=:log, size=(900,450), color=colors[1])
Plots.plot!(bottom_margin=10Plots.mm, right_margin=4Plots.mm, size=(900,450), legend=:bottomright)
Plots.xticks!((b-a)*2.0.^(-8:-3), [latexstring("\\frac{b-a}{2^{"* string(i) * "}}") for i in 8:-1:3])
Plots.savefig("figures/romberg.pdf")

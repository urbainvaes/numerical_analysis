import Symbolics
import LinearAlgebra
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


a, b = π, -1

x = 0:.1:1 |> collect
y = a*x .+ b + .2randn(length(x))

# Data
# x = [0.0; 0.1; 0.2; 0.3; 0.4; 0.5; 0.6; 0.7; 0.8; 0.9; 1.0]
# y = [-0.9187980789440975; -0.6159791344678258; -0.25568734869121856;
#      -0.14269370171581808; 0.3094396057228459; 0.6318327173549161;
#      0.8370437988106428; 1.0970402798788812; 1.6057799131867696;
#      1.869090784869698; 2.075369730726694]

# @variables a b x y
# f = (y -a*x - b)^2/(1+a^2)
# G = Symbolics.gradient(f, [a; b])
# H = Symbolics.hessian(f, [a; b])
# simplify(G[1])
# simplify(H[1, 1])

function grad(p)
    a, b = p
    ∂a = sum(@. (2b*x + 2a*x^2 + 4a*b*y + 2x*y*a^2 - 2a*b^2 - 2x*y - 2a*y^2 - 2b*x*a^2) / ((1 + a^2)^2))
    ∂b = sum(@. (2b + 2a*x - 2y) / (1 + a^2))
    return [∂a; ∂b]
end

function my_hessian(p)
    a, b = p
    ∂aa = sum(@. (2(x^2)) / (1 + a^2) + (-2((y - b - a*x)^2)) / ((1 + a^2)^2) - 2a*((-x*(2y - 2b - 2a* x)) / ((1 + a^2)^2) - 2a*(((y - b - a*x)^2) / ((1 + a^2)^4))*(2 + 2(a^2))) - 2a*((-x* (2y - 2b - 2a*x)) / ((1 + a^2)^2)))
    ∂ab = sum(@. (2x) / (1 + a^2) + (-2a*(2b + 2a*x - 2y)) / ((1 + a^2)^2))
    ∂bb = sum(@. 0*x + 2 / (1 + a^2))
    return [∂aa ∂ab; ∂ab ∂bb]
end

p = [1, 1]
ε = 1e-15
iterates = zeros(2, 0)
while LinearAlgebra.norm(grad(p)) > ε
    p = p - my_hessian(p) \ grad(p)
    println("$(p[1]) $(p[2])")
    # println(LinearAlgebra.norm(grad(p)))
    iterates = [iterates p]
end
print(p)

# error = iterates[:,1:end-1] .- iterates[:,end]
# norm = sqrt.(error[1,:].^2 + error[2,:].^2)
# Plots.plot(norm, yscale=:log10)

xplot = collect(0:.01:1)
Plots.plot(xplot, p[1]*xplot .+ p[2], label="Approximation", linewidth=2)
Plots.scatter!(x, y, markersize=2, label="Data", legend=:topleft)
Plots.savefig("approx_line.pdf")

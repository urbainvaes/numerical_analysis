import Plots
import LinearAlgebra

function hyperball_volume(dim, n)
    number = 0
    for i in 1:n
        x = rand(dim)
        number += LinearAlgebra.norm(x, p) <= 1
    end
    average = number/n
    var = average*(1-average)
    return average * 2^dim, sqrt(var/n) * 2^dim
end

n = 10^7

n_dims = 15
dims, vols, vars = 1:n_dims, zeros(n_dims), zeros(n_dims)
for dim in dims
    vols[dim], vars[dim] = hyperball_volume(dim, 2, n)
end

conf = vars / sqrt(.01)
Plots.scatter(dims, vols, label="Volume estimation")
Plots.plot!(dims, vols, ribbon=conf, fillalpha=0.35, label="99% confidence interval", xlabel="d")
Plots.savefig("hyperballs.pdf")

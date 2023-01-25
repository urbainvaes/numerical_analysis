import Plots
import LinearAlgebra as la
import SparseArrays

# The easy way
function solve(n)
    ρ = cos(π/(n+1))
    ω = 2 / (1 + sqrt(1 - ρ^2))
    b = ones(n)
    L = SparseArrays.spdiagm(-1 => -1*ones(n-1))
    D = SparseArrays.spdiagm(0 => 2*ones(n))
    U = SparseArrays.spdiagm(1 => -1*ones(n-1))
    A = L + D + U
    M = D/ω + L
    N = (1/ω - 1)*D - U 
    x = zeros(n)
    iter = 0
    while la.norm(A*x - b) > 1e-8 * la.norm(b)
        x = M\(N*x + b)
        iter += 1
    end
    return iter
end
@time solve(5000)

# The "manual" way
function solve(n)
    ω = 2 / (1 + sin(π/(n+1)))
    d, u, l = 2n^2, -n^2, -n^2
    x, residual = zeros(n+2), - ones(n)
    iter = 0
    while la.norm(residual) > 1e-8 * sqrt(n)
        println(iter)
        for i in 2:n+1
            x[i] -= ω/d*(l*x[i-1] + d*x[i] + u*x[i+1] - 1)
        end
        residual = d*x[2:n+1] + l*x[1:n] + u*x[3:n+2] .- 1
        iter += 1
    end
    return x[2:end-1]
end

@time solve(5000)

# n = 5000
# ω = 2 / (1 + sin(π/(n+1)))
# d, u, l = 2n^2, -n^2, -n^2
# x, residual = zeros(n+2), - ones(n)
# iter = 0
# while la.norm(residual) > 1e-8 * sqrt(n)
#     println(iter)
#     for i in 2:n+1
#         x[i] -= ω/d*(l*x[i-1] + d*x[i] + u*x[i+1] - 1)
#     end
#     residual = d*x[2:n+1] + l*x[1:n] + u*x[3:n+2] .- 1
#     iter += 1
# end

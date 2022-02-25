# We return "Inf" if the problem is impossible
function T(n, r)
    r == 1 && return 0
    n == 1 && return 1
    r == 2 && return Inf
    optimum = minimum(2T(k, r) + T(n-k, r-1) for k ∈ 1:(n-1))
    return (optimum < Inf) ? Int(optimum) : Inf
end

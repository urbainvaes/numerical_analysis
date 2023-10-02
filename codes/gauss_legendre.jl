using Polynomials
using Plots

function legendre(n)
    p = Polynomial([-1, 0, 1])
    factor = 1 // (2^n * factorial(n))
    return factor * derivative(p^n, n)
end

function get_nodes_and_weights(n)
    nodes = sort(roots(legendre(n)))
    weights = zero(nodes)
    for i in 1:n
        ℓ = fromroots(nodes[1:end .!= i])
        ℓ = ℓ / ℓ(nodes[i])
        weights[i] = integrate(ℓ, -1, 1)
    end
    return nodes, weights
end

function composite_gauss_legendre(u, a, b, n, N)
    h = 2/N
    result = 0.
    X = LinRange(a, b, N + 1)
    x, w = get_nodes_and_weights(n)
    x, w = (x .+ 1)/2, w/2
    for i in 1:N
        result += h * w'u.(X[i] .+ h*x)
    end
    return result
end

# Function to integrate
u(x) = cos(x)

# Integration interval
a, b = -1, 1

# Number of nodes
n = 3

# Exact value of the integral
I = 2sin(1)

N = 1:40
errors = composite_gauss_legendre.(u, a, b, 2, N) .- I
scatter(N, abs.(errors), scale=:log10)
fit(log.(N), log.(abs.(errors)), 1)

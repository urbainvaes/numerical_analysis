using Polynomials
using LaTeXStrings
using Plots

function legendre(n)
    p = Polynomial([-1, 0, 1])
    return 1 / (2^n * factorial(n)) * derivative(p^n, n)
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
    h = (b-a)/N
    X = LinRange(a, b, N + 1)
    z, w = get_nodes_and_weights(n)
    result = 0.
    for i in 1:N
        nodes = X[i] + h/2 .+ z*h/2
        result += h/2 * w'u.(nodes)
    end
    return result
end

# Function to integrate
u(x) = cos(x)

# Integration interval
a, b = -1, 1

# Exact value of the integral
I = 2sin(1)

# Number of nodes
ns = [1, 2, 3]

# Number of cells
N = 1:40

p = plot(title="Convergence of Gauss Legendre quadrature")
for n in ns
    errors = composite_gauss_legendre.(u, a, b, n, N) .- I
    polyfit = fit(log.(N), log.(abs.(errors)), 1)
    α = round(- polyfit[1], digits=2)
    scatter!(N, abs.(errors), scale=:log10, label="n=$n, α=$α")
    xlabel!(L"N")
    ylabel!(L"|I - \widehat I|")
end
display(p)

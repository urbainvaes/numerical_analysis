import FFTW

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

function fourier_interpolate(u, x, X)
    n = (length(x) - 1) ÷ 2
    ξ = FFTW.fft(u.(x))
    result = 0*X
    for k in 0:n
        result += ξ[k+1] * exp.(im*k*X)
    end
    for k in n+1:2n
        result += ξ[k+1] * exp.(im*(k-2n-1)*X)
    end
    return real.(result/(2n+1))
end

import Plots
n, m = 5, 1000
x = 2π/(2n+1) * (0:2n)
X = 2π/m * (0:m)

u(x) = sign(x-π)
# u(x) = exp(sin(x) + cos(5x))
# u(x) = x^2 * (x - 2π)^2 / π^4

@time U = fourier_interpolate(u, x, X)
Plots.plot(X, u.(X), label="u(x)", legend=:bottomright)
Plots.plot!(X, U, label="û(x)")
Plots.scatter!(x, u.(x), label="Interpolation points")
Plots.xlims!(0, 2π)
Plots.savefig("fourier.pdf")

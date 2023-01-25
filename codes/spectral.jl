import Polynomials
import Plots

n = 2
exact(x) = 1/(n*π)^2 * (exp(sin(n*π * x)) - 1)
f(x) = exp(sin(n*π*x)) * cos(n*π*x)^2 - exp(sin(n*π*x)) * sin(n*π*x)

x_interp = LinRange(0, 1, 20)
p = Polynomials.fit(x_interp, f.(x_interp))

inner = Polynomials.integrate(Polynomials.integrate(p))
boundary = Polynomials.Polynomial([0; -inner(1)])
sol = inner + boundary

x_plot = LinRange(0, 1, 200)
Plots.plot(x_plot, sol.(x_plot), legend=:none)
Plots.plot!(x_plot, exact.(x_plot), legend=:none)

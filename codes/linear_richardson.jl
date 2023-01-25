import Plots
import LinearAlgebra

n = 10

A = LinearAlgebra.diagm(0 => 2*ones(n), 1 => -ones(n-1), -1 => -ones(n-1))
b = ones(n)

eigvals, eigvecs = LinearAlgebra.eigen(A)
λmax = maximum(eigvals)
λmin = minimum(eigvals)

ω = .1
M = (1/ω) * LinearAlgebra.I(n)
N = M - A

tol = 1e-10

x = zeros(n)
norm_residue = LinearAlgebra.norm(A*x - b)
while norm_residue > tol
    x = M\(N*x + b)
    norm_residue = LinearAlgebra.norm(A*x - b)
    println(norm_residue)
    sleep(.1)
end

# # Example of an ill-conditioned system

# +
import LinearAlgebra
import Random
Random.seed!(0)

# Shorthand notation
norm = LinearAlgebra.norm
cond = LinearAlgebra.cond
# -

# # Example of an ill-conditioned system
# +
# Parameters of linear system
A = [1 1; 1 (1-1e-12)]
b = [0; 1e-12]

# Use much more accurate BigFloat format
exact_A = [1 1; 1 (1- BigFloat("1e-12"))]
exact_b = [0, BigFloat("1e-12")]

# Relative error on data
relative_error_b = norm(b - exact_b) / norm(b)
relative_error_A = norm(A - exact_A) / norm(A)
println("Relative error on b: $relative_error_b")
println("Relative error on A: $relative_error_A")
# -

# +
# Exact and approximate solutions
x_exact = [1; -1]
x_approx = A\b

relative_error_x = norm(x_exact - x_approx)/norm(x_approx)
println("Relative error on x: $relative_error_x")
# -

# +
# Condition number of A
cond_A = cond(A)
# -

# The relative error on the solution is much larger than the relative error on the data. The ratio between the two is of the same order of magnitude as the condition number of the matrix $A$.

# # LU decomposition

# +
function lu(A)
    n = size(A)[1]
    L = [i == j ? 1.0 : 0.0 for i in 1:n, j in 1:n]
    U = copy(A)
    for i in 1:n-1
        for r in i+1:n
            U[i, i] == 0 && error("Pivotal entry is zero!")
            ratio = U[r, i] / U[i, i]
            L[r, i] = ratio
            U[r, i:end] -= U[i, i:end] * ratio
        end
    end
    return L, U
end

ns = 2 .^ collect(5:9)
times = zeros(length(ns))
for (i, n) in enumerate(ns)
    A = randn(n, n)
    result = @timed lu(A)
    L, U = result.value
    times[i] = result.time
    error = norm(L*U - A)
    println("Error: $(LinearAlgebra.norm(L*U - A)), time: $(result.time)")
end
# -

# # LU decomposition with pivoting

# +
function lu_pivot(A)
    n = size(A)[1]
    L, U = zeros(n, 0), copy(A)
    P = [i == j ? 1.0 : 0.0 for i in 1:n, j in 1:n]
    function swap_rows!(i, j, matrices...)
        for M in matrices
            M_row_i = M[i, :]
            M[i, :] = M[j, :]
            M[j, :] = M_row_i
        end
    end
    for i in 1:n-1

        # Pivoting
        index_row_pivot = i - 1 + argmax(abs.(U[i:end, i]))
        swap_rows!(i, index_row_pivot, U, L, P)

        # Usual Gaussian transformation
        c = [zeros(i-1); 1.0; zeros(n-i)]
        for r in i+1:n
            ratio = U[r, i] / U[i, i]
            c[r] = ratio
            U[r, i:end] -= U[i, i:end] * ratio
        end
        L = [L c]
    end
    L = [L [zeros(n-1); 1.0]]
    return P, L, U
end

ns = 2 .^ collect(5:9)
times = zeros(length(ns))
for (i, n) in enumerate(ns)
    A = randn(n, n)
    L, U = lu(A)
    P, Lpivot, Upivot = lu_pivot(A)
    println("\n-- n = $n --")
    println("Condition numbers without pivoting: κ(L) = $(cond(L)), κ(U) = $(cond(U))")
    println("Condition numbers with    pivoting: κ(L) = $(cond(Lpivot)), κ(U) = $(cond(Upivot))")
end

# -

# +
import LinearAlgebra

# Shorthand notation
norm = LinearAlgebra.norm
cond = LinearAlgebra.cond
# -

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

# The relative error on the solution is much larger than the relative error on the data.
# The ratio between the two is of the same order of magnitude as the condition number of the matrix $A$.

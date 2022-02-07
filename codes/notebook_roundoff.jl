# # Exercise 1.5

# +
# Convert number of the form 0.a₁a₂... to binary and truncate the result after
# a number 'nbits' of bits. The function the bits b₁, b₂, b₃... of the binary
# representation, from the most to the least significative

function to_binary(x, nbits)
    result = zeros(Bool, nbits)
    for i in 1:nbits
        result[i] = (2*x >= 1)
        x = 2*x - result[i]
    end
    return result
end

x, nbits = .1, 60
bits = to_binary(x, nbits)
println(bits)
# -

# +
# Check
approx = BigFloat(0)
for (i, b) in enumerate(bits)
    approx += BigFloat(2.0)^(-i) * b
end
println("Approximation: $approx")
println("\nFloat64 approximation of 0.1:")
println("$(BigFloat(0.1))")
# -

# # Exercise 1.14
# A better formula for computing the sample variance is the following:
# $$
# s^2 = \frac{1}{N-1} \sum_{n=1}^N \bigl( x_n - \bar x \bigr)^2, \qquad \bar x = \frac{1}{N} \sum_{n=1}^N x_n
# $$

# +
import Random

# This ensures that the same random numbers are generated every time the code is run
Random.seed!(0)

N, L = 10^6, 10^9
x = L .+ rand(N)

average = (sum(x)/N)

s1 = 1/(N-1) * (sum(x.^2) - N*average^2)
s2 = 1/(N-1) * (sum((x .- average).^2))
println("Method 1: $s1\nMethod 2: $s2")

# We use 'BigFloat' to calculate a very precise result
x = BigFloat.(x)
exact_average = (sum(x)/N)
exact_value_1 = 1/(N-1) * (sum(x.^2) - N*(sum(x)/N)^2)
exact_value_2 = 1/(N-1) * (sum(x.^2) - N*(sum(x)/N)^2)
println("Exact value: $exact_value_1")
# -

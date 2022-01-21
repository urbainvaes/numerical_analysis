# Define function x
f(x) = exp(x)

# Machine epsilon for Float64
ε = eps()

δs = ε * 2.0 .^ (-2:0)

for δ in δs
    derivative = (f(δ) - f(0))/δ
    println(derivative)
end

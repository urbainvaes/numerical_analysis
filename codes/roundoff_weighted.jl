N, L = 10^6, 10^9
x = L .+ rand(N)

1/(N-1) * (sum(x.^2) - N*(sum(x)/N)^2)

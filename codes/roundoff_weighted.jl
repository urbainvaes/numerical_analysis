N, L = 10^6, 10^9
x = L .+ rand(N)

average = (sum(x)/N)

s1 = 1/(N-1) * (sum(x.^2) - N*(sum(x)/N)^2)
s2 = 1/(N-1) * (sum((x .- average).^2))

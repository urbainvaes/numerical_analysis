function get_primes(N)
    is_prime = [false; ones(Bool, N-1)]
    for i ∈ 2:ceil(Int, sqrt(N))
        is_prime[2i:i:end] .= 0
    end
    return [i for i in 1:N if is_prime[i]]
end

function to_binary(n)
    result = []
    while n != 0
        append!(result, n % 2)
        n ÷= 2
    end
    return result
end

# One-liner
to_binary(n) = n == 0 ? [] : [n%2; to_binary(n÷2)]

# To check that it works
number = 123456789
bits = to_binary(number)
pows2 = 2 .^ range(0, length(bits) - 1)
@assert sum(bits'pows2) == number

struct Float8
    sign :: Bool
    exponent :: Vector{Bool}
    fraction :: Vector{Bool}
end

# function Float8(Float64::x)

# end
function round8(x::Float64)
    exponent = log2(x)
end

function +(x::Float8, y::Float8)
    return 5
end

function bitstring(x::Float8)
    sign = x.sign
    exponent = x.exponent
    fraction = x.fraction
    return "$sign"
end

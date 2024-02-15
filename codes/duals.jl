import Base: +, -, *, /, cos, sin, convert, promote_rule, show
struct D <: Number
    v::Tuple{Number,Number}
end
+(x::D, y::D) = D(x.v .+ y.v)
-(x::D, y::D) = D(x.v .- y.v)
*(x::D, y::D) = D((x.v[1]*y.v[1], (x.v[2]*y.v[1] + x.v[1]*y.v[2])))
/(x::D, y::D) = D((x.v[1]/y.v[1], (y.v[1]*x.v[2] - x.v[1]*y.v[2])/y.v[1]^2))
-(x::D) = D(.-(x.v))
cos(x::D) = D((cos(x.v[1]), -sin(x.v[1]) * x.v[2]))
sin(x::D) = D((sin(x.v[1]), cos(x.v[1]) * x.v[2]))
convert(::Type{D}, x::Real) = D((x,zero(x)))
promote_rule(::Type{D}, ::Type{<:Real}) = D
show(io::IO,x::D) = print(io, "(", x.v[1], ")" ," + (", x.v[2],") ε")
ε = D((0, 1))

# p-th derivative of f at x
function derivative(f, p, x)
    p == 0 && return f(x)
    return derivative(f, p - 1, D((x, 1))).v[2]
end

using Plots
x = LinRange(0, 2, 200)
f(x) = sin(x)^2 / (x^2 - x + 1)
plot(x, f.(x))
plot!(x, derivative.(f, 1, x))
plot!(x, derivative.(f, 2, x))
plot!(x, derivative.(f, 3, x))

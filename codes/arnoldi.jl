import SparseArrays
using LinearAlgebra

n = 500
A = Tridiagonal(-1*ones(n-1), 2*ones(n), -1*ones(n-1))
x = BigFloat.(rand(n))

function power(A, v)
    v /= norm(v)
    Av = A*v
    λ = v'A*v
    error = [norm(A*v-λ*v)]
    while error[end] > 1e-12
        println(error[end])
        v = A*v
        v /= norm(v)
        Av = A*v
        λ = v'A*v
        append!(error, norm(A*v-λ*v))
    end
    return λ, v, error
end

function lanczos(v, p)
    v /= norm(v)
    λ = v'A*v
    error = [norm(A*v-λ*v)]
    U = zeros(n, p + 1)
    U[:, 1] = v/sqrt(v'v)
    α, β = zeros(p), zeros(p)
    for j in 1:p
        println(error[end])
        w = A*U[:, j] - (j == 1 ? zeros(n) : β[j-1]*U[:, j-1])
        α[j] = U[:, j]'w
        w -= α[j]*U[:, j]
        β[j] = sqrt(w'w)
        U[:, j+1] = w/β[j]
        H = Tridiagonal(β[1:j-1], α[1:j], β[1:j-1])
        λ, y, _ = power(H, ones(1))
        v = U[:, 1:j]*y
        λ = v'*A*v / v'v
        append!(error, norm(A*v-λ*v))
    end
    return λ, v, error
end

function arnoldi(v, m)
    v /= norm(v)
    λ = v'A*v
    error = [norm(A*v-λ*v)]
    U = zeros(Complex, n, m + 1)
    U[:, 1] = v/sqrt(v'v)
    H = zeros(Complex, m + 1, m)
    for j in 1:m
        println(error[end])
        w = A*U[:, j]
        H[1:j, j] = U[:, 1:j]'w
        w -= U[:, 1:j] * H[1:j, j]
        H[j+1, j] = sqrt(w'w)
        U[:, j+1] = w/H[j+1, j]
        λ, y, _ = power(H[1:j, 1:j], ones(j))
        v = U[:, 1:j]*y
        λ = v'*A*v / v'v
        append!(error, norm(A*v-λ*v))
    end
    return λ, v, error
end

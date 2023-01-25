b, n = 5, 1000
A = [abs(i-j) <= b ? rand() : 0.0 for i ∈ 1:n, j ∈ 1:n]

function lu(A)
    n, n = size(A)
    L, U = zeros(n, n), zeros(n, n)
    for i in 1:n
        L[i, i] = 1
        for j in 1:n
            if j < i
                L[i, j] = A[i, j]
                for k in 1:j-1
                    L[i, j] -= L[i, k]*U[k, j]
                end
                L[i, j] /= U[j, j]
            else
                U[i, j] = A[i, j]
                for k in 1:i-1
                    U[i, j] -= L[i, k]*U[k, j]
                end
            end
        end
    end
    return L, U
end

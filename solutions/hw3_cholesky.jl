import BandedMatrices
import LinearAlgebra
using BenchmarkTools

function cholesky1(A)
    _, n = size(A)
    C = copy(LinearAlgebra.LowerTriangular(A))
    b = A.u
    for i ∈ 1:n, j ∈ max(1, i-b):i
        C[i, j] -= @view C[i, max(1, i-b):j-1]'C[j, max(1, i-b):j-1]
        C[i, j] = (j < i) ? C[i, j] / C[j, j] : sqrt(C[i, j])
    end
    return C
end

function cholesky2(A)
    m, n = size(A)
    m != n && error("Matrix must be square")
    # Convert to banded matrix
    A = BandedMatrices.BandedMatrix(A)
    A.u != A.l && error("Matrix must be symmetric")
    C = copy(LinearAlgebra.LowerTriangular(A))
    for i ∈ 1:n, j ∈ max(1,i-A.u):i
        for k ∈ max(1,i-A.u):j-1
            C[i, j] -= C[i, k]*C[j, k]
        end
        C[i, j] += A[i, j]
        if j < i
            C[i, j] /= C[j, j]
        else
            C[i, j] = sqrt(C[i, j])
        end
    end
    return C
end

n, l, u = 2000, 2, 2;
A = BandedMatrices.brand(n, l, u);
A = A'*A;

@belapsed cholesky1(A)
C = @time cholesky2(A);





























# LinearAlgebra.norm(C*C' - A, Inf)

# import BandedMatrices
# import LinearAlgebra
# bm = BandedMatrices


# n, l, u = 20000, 2, 2
# A = BandedMatrices.brand(n, l, u)
# A = A'*A

# C = @time cholesky(A)
# LinearAlgebra.norm(C*C' - A, Inf)

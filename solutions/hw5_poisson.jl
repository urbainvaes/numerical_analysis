using LinearAlgebra
import SparseArrays
import Plots
import PlotlyJS

# Domain size
Lx, Ly = 2, 1

# Number of discretization points along x and y, including the boundary points
nx, ny = 101, 101

function discretize(nx, ny)
    hx, hy = Lx/(nx - 1), Ly/(ny - 1)
    Dxx = (1/hx^2) * Tridiagonal(-ones(nx-3), 2ones(nx-2), -ones(nx-3))
    Dyy = (1/hy^2) * Tridiagonal(-ones(ny-3), 2ones(ny-2), -ones(ny-3))
    A = kron(Dxx, I(ny-2)) + kron(I(nx-2), Dyy)
    xgrid = Lx/(nx-1) * (1:nx-2)
    ygrid = Ly/(ny-1) * (1:ny-2)
    x_2d = reshape([x for y in ygrid, x in xgrid], (nx-2)*(ny-2))
    y_2d = reshape([y for y in ygrid, x in xgrid], (nx-2)*(ny-2))
    b = sin.(4π*x_2d) + sin.(2π*y_2d)
    return SparseArrays.SparseMatrixCSC(A), b
end

function plot_solution(f)
    f = reshape(f, ny-2, nx-2)

    # Boundary condition
    z = [zeros(nx)'; zeros(ny-2) f zeros(ny-2); zeros(nx)']
    xgrid = Lx/(nx-1) * (0:nx-1)
    ygrid = Ly/(ny-1) * (0:ny-1)

    # Plots.contourf(xgrid, ygrid, z, c=:viridis, levels=50)
    PlotlyJS.plot(PlotlyJS.contour(x=xgrid, y=ygrid, z=z))
end

# Calculate matrix and right-hand side of linear system
A, b = discretize(nx, ny)

########################
# Your code comes here #

f = A\b
########################

# Plot the solution
plot_solution(f)
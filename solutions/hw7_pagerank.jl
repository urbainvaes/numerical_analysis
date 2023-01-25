import CSV
import SparseArrays
import LinearAlgebra
import DataFrames
norm = LinearAlgebra.norm

# data
nodes = CSV.read("names.csv", DataFrames.DataFrame)
edges = CSV.read("edges.csv", DataFrames.DataFrame)

# nodes = CSV.read("new_nodes.csv", DataFrames.DataFrame)
# edges = CSV.read("new_edges.csv", DataFrames.DataFrame)







































# Convert data to matrices
nodes = Matrix(nodes)
edges = Matrix(edges)
nn, ne = length(nodes), size(edges)[1]

# Count the number of outbound and inbound edges for each node
nn, ne = length(nodes), size(edges)[1]
n_outbound = zeros(Int, nn)
for e in eachrow(edges)
    n_outbound[e[1]] += 1
end

# The data no longer contains sinks
I, J = edges[:, 1], edges[:, 2]
V = 1 ./ n_outbound[I]
adjacency = SparseArrays.sparse(J, I, V, nn, nn)


x = ones(nn) / nn
my_prod = adjacency*x
while norm(my_prod - x, 1) / norm(x, 1) > 1e-15
    x = my_prod
    my_prod = adjacency*my_prod
    iter += 1
end

x = x ./ LinearAlgebra.norm(x)
x = sortperm(x, rev=true)
nodes[x]


function search(keyword)
    return [n for n in nodes if occursin(keyword, n)]
    return filter(s -> occursin(keyword, s), nodes)
end

function clean(keep)
    nodes = nodes[keep]
    edges = edges[keep[edges[:, 1]], :]
    edges = edges[keep[edges[:, 2]], :]
    conversion = cumsum(keep .== 1)
    edges = conversion[edges]
    df_nodes = DataFrames.DataFrame([nodes], :auto)
    df_edges = DataFrames.DataFrame(edges, :auto)
    CSV.write("new_nodes.csv", df_nodes)
    CSV.write("new_edges.csv", df_edges)
end

import CSV
import SparseArrays
import DataFrames

# data
names = CSV.read("names.csv", DataFrames.DataFrame)
edges = CSV.read("edges.csv", DataFrames.DataFrame)

# Convert data to matrices
names = Matrix(names)
edges = Matrix(edges)
nn, ne = length(names), size(edges)[1]

# Count the number of outbound and inbound edges for each node
nn, ne = length(names), size(edges)[1]
n_outbound = zeros(Int, nn)
n_inbound = zeros(Int, nn)
for e in eachrow(edges)
    n_outbound[e[1]] += 1
    n_inbound[e[2]] += 1
end

# The data no longer contains sinks
I, J = edges[:, 1], edges[:, 2]
V = 1 ./ n_outbound[I]
adjacency = SparseArrays.sparse(I, J, V, nn, nn)

x0 = ones(1, nn) / nn

for i in 1:500
    print(i)
    x0 = x0*adjacency
end
x0 = x0[1, :]

keep = sort(collect(zip(1:nn, x0)), by = x -> x[2], rev=true)
keep_names = sort([k[1] for k in keep])
df = DataFrames.DataFrame([new_names], :auto)
CSV.write("new_names.csv", df)
new_names = names[keep_names]

threshold = sort(x0, rev=true)[2*10^5]
new_names = findall(x0 .>= threshold)

conversion = cumsum(x0 .>= threshold)

keep_edges = edges[x0[edges[:, 1]] .>= threshold, :]
keep_edges = keep_edges[x0[keep_edges[:, 2]] .>= threshold, :]
new_edges = conversion[keep_edges]


df = DataFrames.DataFrame([names], :auto)
CSV.write("names.csv", df)

df = DataFrames.DataFrame(edges, :auto)
CSV.write("edges.csv", df)


# adjacency = SparseArrays.spzeros(n, n)

function clean(names, edges)
    nn, ne = length(names), size(edges)[1]

    # Count the number of outbound and inbound edges for each node
    nn, ne = length(names), size(edges)[1]
    n_outbound = zeros(Int, nn)
    n_inbound = zeros(Int, nn)
    for e in eachrow(edges)
        n_outbound[e[1]] += 1
        n_inbound[e[2]] += 1
    end

    # For simplicity, we will remove the sinks
    not_sinks = findall(n_outbound .> 0)

    # Remove edges to sinks
    edges = edges[n_outbound[edges[:, 2]] .> 0, :]

    # Remove sinks from the nodes
    names = names[not_sinks]

    # Conversion array from old node indices to new ones
    conversion = cumsum(n_outbound .> min_outbound)

    # Update edges to new node indices
    edges = conversion[edges]

    return names, edges
end

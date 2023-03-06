import CSV
import SparseArrays
import DataFrames

# data
names = CSV.read("enwiki-2013-names.csv", DataFrames.DataFrame, escapechar='\\')
edges = CSV.read("enwiki-2013.txt", DataFrames.DataFrame, header=["From"; "To"], skipto=5, delim=' ')

# Convert data to matrices
names = names[:, 2]
edges = Matrix(edges) .+ 1

nn, ne = length(names), size(edges)[1]

# Count the number of outbound and inbound edges for each node
nn, ne = length(names), size(edges)[1]
n_outbound = zeros(Int, nn)
n_inbound = zeros(Int, nn)
for e in eachrow(edges)
    n_outbound[e[1]] += 1
    n_inbound[e[2]] += 1
end

# Find nodes without any outbound edges
sinks = findall(n_outbound .== 0)

# Add 'n' artificial outbound edges from sinks
n = 10
for i in 1:n
    to = floor.(Int, nn*rand(length(sinks)))
    new_edges = [sinks to]
    edges = [edges; new_edges]
    n_outbound[sinks] .+= 1
    n_inbound[to] .+= 1
    ne += length(sinks)
end


# The data no longer contains sinks

I, J, V = edges[:, 1], edges[:, 2], 1 ./ n_outbound[I]
adjacency = SparseArrays.sparse(I, J, V, nn, nn)

x0 = ones(1, nn) / nn
for i in 1:100
    print(i)
    x0 = x0*adjacency
end

sort(collect(zip(names, x0)), by = x -> x[2], rev=true)

# adjacency = SparseArrays.spzeros(n, n)


x = rand(nn)
A = adjacency

B = SparseArrays.SparseMatrixCSC(A)
for r in 1:nn
    if r % 100 == 0
        println(r)
    end
    if length(B[r, :]) == 0
        print("zero")
        # adjacency[e[1] +  1, e[2] + 1] = 1
    end
end


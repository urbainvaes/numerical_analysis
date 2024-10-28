# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.1
#   kernelspec:
#     display_name: Julia 1.10.2
#     language: julia
#     name: julia-1.10
# ---

# ## Notebook 1: Introduction to Julia

using Plots
using LaTeXStrings

# + [markdown] jp-MarkdownHeadingCollapsed=true
# ### <font color='green'>Understanding the output of Jupyter cells</font>
# -

# - The output of a series of Julia commands is the last expression:

# +
# Jupyter will display the result of `a + b` after running the cell

a = 5
b = 4
c = a + b
# -

# - Using `;` at the end of a command suppresses output

# +
# Jupyter will not display anything

a = 5
b = 4
a + b;
# -

# - This holds also for plots. Compare the following two situations:

plot(cos)

plot(sin);

# - You can also use `@show` or `print`:

a = 5
b = 4
@show a
println(b)
a + b

# + [markdown] jp-MarkdownHeadingCollapsed=true
# ### <font color='green'>Manipulation of strings</font>
# -

# - Strings are defined in Julia between double quotes `"`. In contrast with Python, single quote `'` are used for single characters only.

s = "Hello world!"

'o' ∈ s # ∈ is an alias for in

# - Unlike many other languages using `+`, Julia uses `*` for string concatenation. The idea is to emphasize the non-commutativity of this operation

s1 = "Hello " ; s2 = "world!"
s = s1 * s2
println(s)

print(s1^5)

# - Interpolation within a string can be obtained using `$`

col = "red" ; x = 20_000 ; discount = 5
s = "My $col car cost $(x*(1-discount/100)) €"
print(s)

s = "Hello world!"
ordinal_suffix(i) = i == 1 ? "st" : i == 2 ? "nd" : "th"
for (i, c) ∈ enumerate(s) 
    println("The character at $i$(ordinal_suffix(i)) place is $c") 
end

# > Note that LaTeX strings can be built thanks to the [LaTeXStrings](https://github.com/JuliaStrings/LaTeXStrings.jl) package, imported above and used further in this notebook.

# + [markdown] jp-MarkdownHeadingCollapsed=true
# ### <font color='green'>Using and defining functions</font>
# -

# #### <font color='orange'>[Exercise 1]</font>
# Write a function that returns the sum of the `n` first integers using a `for` loop.
#
# > *Hints:*
# >  - read the documentation of `function` (open a REPL and type `?` followed by the searched word, here `function`),
# >  - a range of integers from `m` to `n` simply writes `m:n` <mark>(unlike Python bounds are included here in the range)</mark>.

# + nbgrader={"grade": false, "grade_id": "squares", "locked": false, "schema_version": 3, "solution": true, "task": false}
function sum_of_squares(n)
    ### BEGIN SOLUTION
    result = 0
    for i in 1:n
        result += i*i
    end
    return result
    ### END SOLUTION
end;

# + nbgrader={"grade": true, "grade_id": "squares-tests", "locked": true, "points": 2, "schema_version": 3, "solution": false, "task": false}
@assert sum_of_squares(1) == 1 "❌"
@assert sum_of_squares(2) == 5 "❌"
@assert sum_of_squares(5) == 55 "❌"
@assert sum_of_squares(10) == 385 "❌"
# -

# **Remarks:**
# - The same function can be written in short form

# +
short_sum_of_squares(n) = (1:n)'*(1:n)

@show short_sum_of_squares(2)
@show short_sum_of_squares(4);
# -

# - Read the documentation of the function `sum` to understand why the following code also works.

# +
other_sum_of_squares(n) = sum(x -> x^2, 1:n)

@show other_sum_of_squares(2)
@show other_sum_of_squares(4);
# -

# - The function can be applied elementwise to an array by using the broadcasting operator `.`

# +
@show short_sum_of_squares.([2, 4, 5])

# The above command is syntactic sugar for
@show broadcast(short_sum_of_squares, [2, 4, 5]);
# -

# - Broadcasting works even for functions with multiple arguments

geometric_mean(a, b) = √(a*b)
@show geometric_mean.(2, [8, 18]);

# - Sometimes, it is convenient to define *anonymous* functions with the `->` notation:

# Plot the sum of cubes for different n
scatter(n -> sum(x -> x^3, 1:n), 1:5, label="sum of cubes")

# #### <font color='orange'>[Exercise 2]</font>
#
# The ancient Babylonian method for calculating the square root of $a > 0$ is given by the iteration
# $$
# ∀k∈\mathbb{N},\qquad x_{k+1}=\frac{1}{2}\left(x_k+\frac{a}{x_k}\right)
# $$
#
# Write a function `my_sqrt(a, x₀)` that returns an approximation of the square root of `a` based on 10 iterations of this method, starting from `x₀`.

# + nbgrader={"grade": false, "grade_id": "baby", "locked": false, "schema_version": 3, "solution": true, "task": false}
function my_sqrt(a, x₀)
    ### BEGIN SOLUTION
    x = x₀
    for i in 1:10
        x = (x + a/x)/2
    end
    return x
    ### END SOLUTION
end;

# + nbgrader={"grade": true, "grade_id": "baby-tests", "locked": true, "points": 2, "schema_version": 3, "solution": false, "task": false}
@assert my_sqrt(4, 1) ≈ 2 "❌"
@assert my_sqrt(144, 1) ≈ 12 "❌"

# + [markdown] jp-MarkdownHeadingCollapsed=true
# ### <font color='green'>Plotting and creating animations</font>
# -

# - The function `plot` can be called in a variety of ways.
#   One approach is to pass it a function:

f(x) = cos(x) * exp(-abs(x))
plot(f, xlims=(-5, 5))

# - Another, often more flexible, option is pass it $x$ and $y$ arrays:

f(x) = cos(x) * exp(-abs(x))
x = -5:.01:5
plot(x, f.(x))

# - Use commands followed by `!` to modify an existing plot.

plot(cos, label = "f1")
plot!(sin, label = "f2")
xlims!(-5, 5)
title!("title")
xlabel!("x")
ylabel!("y")

# - Options can be grouped together within a `plot!` call:

plot(cos, label = "f1")
plot!(sin, label = "f2")
plot!(xlims = (-5, 5), title = "title", xlabel = "x", ylabel = "y")

# - Animations can be created as follows:

anim = Animation()
for θ in LinRange(0, 2π, 50)
    plot(x -> cos(x - θ), label=nothing)
    title!("θ = $(round(θ, digits=3))")
    frame(anim)
end
gif(anim, fps=10)

# - The `Plots` library defines the `@animate` macro to simplify the creation of animations:

anim = @animate for θ in LinRange(0, 2π, 50)
    plot(x -> cos(x - θ), label=nothing)
    title!("θ = $(round(θ, digits=3))")
end
gif(anim, fps=10)

# #### <font color='orange'>[Exercise 3]</font>
# Create an animation illustrating the motion a golf ball whose position in the $x,y$ plane with position
# $$
# \begin{aligned}
# x(t) &= u_0 t, \\
# y(t) &= v_0 t - g \frac{t^2}{2},
# \end{aligned}
# $$
# Use the parameters provided as well as the function `scatter`,
# and make the animation as good as possible.

# + nbgrader={"grade": true, "grade_id": "golf", "locked": false, "points": 4, "schema_version": 3, "solution": true, "task": false}
u₀ = 50
v₀ = 50
g = 9.81
### BEGIN SOLUTION
X(t) = u₀*t
Y(t) = v₀*t - g*t^2/2
tf = 2v₀/g
nframes = 100

anim = @animate for t in LinRange(0, tf, nframes)
    scatter([X(t)], [Y(t)], label=nothing, markersize=10)
    title!("t = $(round(t, digits=3))")
    xlims!(0, X(tf))
    ylims!(0, Y(tf/2) * 1.1)
end
gif(anim, fps=nframes/tf)
### END SOLUTION

# + [markdown] jp-MarkdownHeadingCollapsed=true
# ### <font color='green'>Floating point formats</font>
# -

# Only finitely many numbers can be represented exactly on a computer.
# Consequently, some real numbers can be represented exactly,
# and some cannot.
# The IEEE 754 standard,
# with which virtually all programming languages comply,
# specifies a small number of **floating point formats**.
# The set of real numbers representable in each format is of the form
# $$
# \begin{align*}
#     F(p, E_{\min}, E_{\max})
#     = \Bigl\{ & (-1)^s 2^E (b_0. b_1 b_2 \dots b_{p-1})_2 \colon \\
#               & \qquad s \in \{0, 1\}, b_i \in \{0, 1\} \, \text{and} \, E_{\min} \leq E \leq E_{\max} \Bigr\}.
# \end{align*}
# $$
# Here $(b_0.b_1 b_2 \dots b_{p-1})_2$ is a number in binary representation, i.e.
# $$
# (b_0.b_1 b_2 \dots b_{p-1})_2 = b_0 + \frac{1}{2} b_1 + \dotsc + \frac{1}{2^{p-1}} b_{p-1}.
# $$
# For an element $x \in F(p, E_{\min}, E_{\max})$, $(b_0. b_1 b_2 \dots b_{p-1})_2$ is called the *significand*, $s$ is the *sign* and $E$ is the *exponent*.
# Three parameters appear in the set definition given above:
# - $p$ is the number of significant bits (also called the **precision**),
# - $E_{\min}$ is the minimum exponent.
# - $E_{\max}$ is the maximum exponent.
#
# From the precision, the **machine epsilon** associated with the format is defined as
# $$
# \varepsilon = 2^{-(p-1)}
# $$
#
# The parameters of the most commonly-used formats specified by the standard are the following:
#
# | Parameter | Half-precision (`Float16`) | Single precision (`Float32`) | Double precision (`Float64`) |
# | --- | --- | --- | --- |
# | $p$ | 11 | 24 | 53 |
# | $E_{\min}$ | -14 | -126 | -1022 |
# | $E_{\max}$ | 15 | 127 | 1023 |
#
# In Julia, half, single and double precision floating point formats are given by
# `Float16`, `Float32`, and `Float64`, the default one being `Float64`.

# +
a = √2
b = Float32(√2)
c = Float16(√2)

@show a
@show b
@show c

@show typeof(a)
@show typeof(b)
@show typeof(c);
# -

# The machine epsilon should be understood as a measure of the relative spacing between a floating point number and the next one, within a format.
# Indeed, the number 1 is representable exactly in all formats,
# and the next representable number is precisely $1 + \varepsilon$

a = Float16(1.0)
@show nextfloat(a) - a
@show eps(Float16);

# Likewise, the next representable number from $2$ is $2 + 2\varepsilon$
# and the next representable number from $4$ is $4 + 4ε$:

ε = eps(Float64)
@show nextfloat(2.0) - 2.0
@show 2ε;
@show nextfloat(4.0) - 4.0
@show 4ε;

# As expected from the names,
# a `Float16` number is stored using 16 bits, 
# a `Float32` is stored using 32 bits, 
# and a `Float64` number is stored using 64 bits.
# Of these bits, one is used to store the sign, 
# and the rest are used to store the significand and the exponent.
# Additionally, a few combinations of the bits are reserved to represent the special values
# `Inf`, `-Inf` and `NaN`, the latter being an acronym for "not a number".

# A **byte**, abbreviated with a capital B, is a unit that consists of 8 bits.
# This is the unit that is usually used to describe the amount of memory available in commercial products;
# for example, a laptop typically has around 8GB or 16GB of RAM memory.
# It is often useful to have a rough idea of how much space mathematical objects take in memory.
# For example, a 1000 x 1000 matrix of `Float64` numbers requires $8 \times 10^6$ bytes to store,
# and a 10000 x 10000 matrix of `Float64` numbers requires almost 800MB of RAM.

A = zeros(10_000, 10_000)
@show Base.summarysize(A);

# To develop our understanding of floating point format,
# let us plot on the same graph the spacing between
# successive `Float16`, `Float32` and `Float64` numbers
# in the range $[1, 10^4]$. Notice that the density of
# representable floating point numbers around $x$ decreases
# with $x$.

n = 1000
x_16 = (^).(10, LinRange{Float16}(0, 4, n))
x_32 = (^).(10, LinRange{Float32}(0, 4, n))
x_64 = (^).(10, LinRange{Float64}(0, 4, n))
spacing16 = nextfloat.(x_16) - x_16
spacing32 = nextfloat.(x_32) - x_32
spacing64 = nextfloat.(x_64) - x_64
plot(x_16, spacing16, label="Float16", )
plot!(x_32, spacing32, label="Float32")
plot!(x_64, spacing64, label="Float64")
plot!(xscale=:log10, yscale=:log10,
      xlabel="x", ylabel="Δx", legend=:bottomright)

# #### <font color='orange'>[Exercise 4]</font>
# To further develop your understanding,
# - Read the documentation of the `nextfloat` function.
# - Using this function,
#   plot the **relative** spacing, normalized by the machine epsilon of the `Float32` format,
#   between successive `Float32` numbers in the range $[10^{-2}, 10^2]$.
#   More precisely, plot the function
#   $$
#   x \mapsto \frac{\text{nextfloat}(x) - x}{x \varepsilon}
#   $$
#
# Use a logarithmic scale for the $x$ axis only.
# You may find it useful to use `LinRange{Type}(a, b, n)`
# to create a vector of $n$ equidistant numbers of type `Type` between $a$ and $b$.

# + nbgrader={"grade": true, "grade_id": "spacing", "locked": false, "points": 4, "schema_version": 3, "solution": true, "task": false}
### BEGIN SOLUTION
x = (^).(10, LinRange{Float32}(-2, 2, 10000))
spacing = (nextfloat.(x) - x) ./ (x*eps(Float32))
plot(x, spacing, label="Relative spacing")
plot!(xscale=:log2, xlabel="x", ylabel="Δx/x (in ε units)")
### END SOLUTION
# -

# Let us now address a natural question stemming from this discussion:
# what happens if the result of a mathematical operation is not exactly representable in the floating point format?
# The answer is simple: the result is rounded to the nearest representable number.
# This explains the following surprising results:

@show 1 + eps()/3      # the result is 1
@show 1 + 2*eps()/3    # the result is 1 + ε
@show .1 + .2 == .3    # this is False

# The errors caused by this rounding process are called **roundoff errors**.
# They are particularly impactful for the calculation of derivatives.
# Indeed, suppose that we wish to calculate,
# for a differentiable function $f\colon \mathbb R \to \mathbb R$, an approximation of $f'(1)$.
# The most natural approach to this end is to calculate for small $\delta > 0$
# $$
# d(\delta) = \frac{f(1+ \delta) - f(1)}{\delta}.
# $$
# In exact arithmetic, this would provide an approximation of $f'(1)$ with an error scaling as $\mathcal O(\delta)$,
# assuming $f''(1) \neq 0$.
# In computer arithmetic, however, an additional error occurs due to rounding,
# which can be shown to be of the order $\mathcal O(\varepsilon/\delta)$.
# The best choice is to set $\delta \approx \sqrt{\varepsilon}$,
# so that the errors are approximately of the same order of magnitude.

f(x) = exp(x)
d(δ) = (f(1+δ) - f(1))/δ
δs = 10 .^(0:-.1:-17)
err = abs.(d.(δs) .- exp(1))
plot(δs, err, xscale=:log10, yscale=:log10, label=L"|d(\delta) - f'(1)|", size=(900,450))
plot!([eps()], seriestype=:vline, label = L"\varepsilon")
plot!([sqrt(eps())], seriestype=:vline, label = L"\sqrt{\varepsilon}")

# For very small $\delta$ of the order of the machine $\varepsilon$, 
# the approximation obtained is completely wrong.
# Do you understand the following observations?

δ = eps()/4; @show (exp(δ) - exp(0)) / δ
δ = eps()/2; @show (exp(δ) - exp(0)) / δ
δ = eps(); @show (exp(δ) - exp(0)) / δ;

# **Remark**. In many settings, it is useful to use floating point numbers with more precision than the standard types,
# so that the roundoff errors can be reduced and true error, 
# coming from the inexact nature of many algorithms even in exact arithmetic,
# can be observed.
# The `BigFloat` in Julia is an arbitrary precision format.
# The corresponding number of significant bits can be set by using the `setprecision` function.

# +
@show Float64(π)    

setprecision(100, base=10)  # We want 100 digits in base 10
BigFloat(π)

# + [markdown] jp-MarkdownHeadingCollapsed=true
# ### <font color='green'>Manipulation of arrays</font>
# -

# - Simple explicit definition

u = [1, 2, 3]  # defines a column Vector
display(u)
@show typeof(u)

# - Note the difference with

u = [1 2 3] # defines a Matrix of 1 line
display(u)
@show typeof(u)

# - Explicit definition of a matrix

A = [1 2 3
     4 5 6
     7 8 9]

A = [1 2 3; 4 5 6; 7 8 9]

# - Extraction of a column or a line
# > **Note:** in `Julia`, indices start at `1` as in `Fortran` (**not `0` as in `Python` or `C`**)

A[1,:]

A[:,1]

# - `:` is a shortcut for `begin:end` or `1:end`

A[2,1:end]

A[begin:end,2]

# - Construction of an array by comprehension syntax

A = [10i+j for i in 0:2, j in 0:3]
display(A)
print("\nA = $A\n\n")
@show A
@show A[1,1]
@show A[1,2]
@show length(A)
@show size(A)
@show size(A, 1)
@show size(A, 2)
@show eltype(A)
;

# - Adjoint matrix

A = [1 2 3
     4 5 6]
A'

display(A'A)
display(A*A')

A = [2+3im 4+2im 1+3im
     1+2im 2+3im 5+6im]
A'

display(A'A)
display(A*A')

# - Concatenation of arrays

A = [1 2 3
     4 5 6]
u = [10, 20, 30]
v = [100, 200]
;

[A v]

[A; u']

[A A]

[A; A]

hcat(A, A, A)

vcat(A, A, A)

# - Components are stored column-wise and can be accessed by single index even for multidimensional arrays
# - `eachindex` and `CartesianIndices` allow to iterate over all the indices of an array (single index for the former, multi-index for the latter)
# - A `CartesianIndex` can be converted into a `Tuple`
# - `...` is the splat operator (see the documentation)

# +
A = [1 2 3; 4 5 6]

for k ∈ eachindex(A)
    println("A[$k] = ", A[k])
end
println()

for k ∈ CartesianIndices(A)
    println("A[$k] = A$(Tuple(k)...) = ", A[k])
end
println()

for i in 1:size(A, 1)
    for j in 1:size(A, 2)
        println("A[$i, $j] = ", A[i, j])
    end
end
println()

for i in 1:size(A, 1), j in 1:size(A, 2)
    println("A[$i, $j] = ", A[i, j])
end
# -

# - A matrix (array of order 2) is not a vector of vectors

B = [[10i+j for i in 0:2] for j in 0:3]
display(B)
@show B isa Matrix 
;

# - Arrays are mutable (if the type is consistent), tuple are not!

v = [1, 2]
@show typeof(v)
@show eltype(v)
@show v
v[1] = 3
@show v
# v[1] = 1.5 # error since v is of type Vector{Int64}
# @show v
t = (1, 2)
@show t[1] ;
# t[1] = 3 # error since t is not mutable

# - Definition of a void array and *a posteriori* construction
#
# > Note the difference between `push!` and `append!` (it is recalled that for both `!` expresses that by convention, not obligation, the argument is modified)

v = [] # void vector containing any type (not optimal!)
@show typeof(v)
push!(v, 1)
append!(v, [2.5, 3.6]) # pushes each item of the vector [2.5, 3.6]
@show v
push!(v, [4.1, 5.2]) # pushes a vector as new item
@show v ;
@show v[end]
@show v[end-1] ;

# - For better performance, if possible, it is preferable to inform about the type of the elements of the vector at construction

v = Float64[]
@show typeof(v)
push!(v, 1)
append!(v, [2//3, π]) # pushes each item of the vector [2//3, π] and convert on-the-fly to Float64
@show v
# push!(v, [4.3, 5.7]) # throws an error since [4.3, 5.7] is not of type Float64
;

# - Reference of an array

A = [10i+j for i in 1:2, j in 1:3]
@show A
B = A # B is a reference to A (both point to the same memory location)
@show B == A # same content
@show B === A # same memory location
B .= 99 # .= assigns all the coefficients of the array to 99
@show B
@show A # A is of course also affected
;

# - Copy of an array

A = [10i+j for i in 1:2, j in 1:3]
@show A
B = copy(A) # B and A are two distinct arrays
@show B == A # same content
@show B === A # not the same memory location
B .= 99 # .= assigns all the coefficients of the array to 99
@show B
@show A # A unaffected
;

# - Slicing as a l-value ≡ **reference**

A = [10i+j for i in 1:6, j in 1:8]
display(A)
A[1:2:end, 2:2:end] .= 0 # components of odd line and even column are set to 0
display(A)
;

# - Slicing as a r-value ≡ **copy**

A = [10i+j for i in 1:6, j in 1:8]
display(A)
B = A[1:2:end, 2:2:end] # B is a copy of the submatrix
display(B)
B .= 0
display(B)
display(A) # A remains unaffected
;

# - Reference as r-value by means of `view`

A = [10i+j for i in 1:6, j in 1:8]
display(A)
B = view(A, 1:2:size(A,1), 2:2:size(A,2)) # B is a view of the submatrix (reference not a copy)
display(B)
B .= 0
display(B)
display(A) # A is modified
;

# - Other ways to define slices by selection of indices

A = [10i+j for i in 1:6, j in 1:8]
A[isodd.(1:end), iseven.(1:end)] .= 0 ; display(A)
A = [10i+j for i in 1:6, j in 1:8]
A[findall(isodd,A)] .= 0 ; display(A) ;

# - Explain the following lines

A = [10i+j for i in 1:6, j in 1:8]
B = view(A, collect(1:size(A,1)) .% 3 .== 0, collect(1:size(A,2)) .% 5 .!= 0)
B .= 0
A

A = rand(5, 5) ; display(A)
A[findall(x -> x < 1/2, A)] .= 0
A

n, m = 7, 9
B = rand(n, m)
u = rand(m)
display(B*u)
subi = collect(1:n) .% 3 .!= 0 ; @show subi
subj = collect(1:m) .% 3 .== 0 ; @show subj
B[subi,subj]*u[subj]

# #### <font color='orange'>[Exercise 5]</font>
#
# Write a function which computes the comatrix of a matrix explicitly from the determinant of submatrices.
#
# It is recalled that the component $C_{ij}$ of the comatrix of $A$ is obtained as the determinant of the submatrix obtained from $A$ by removing its $i^\textrm{th}$ row and $j^\textrm{th}$ column multiplied by $(-1)^{i+j}.$ The determinant of a small dense matrix is computed by the function `det` available in the library `LinearAlgebra`.

# + nbgrader={"grade": false, "grade_id": "comatrix", "locked": false, "schema_version": 3, "solution": true, "task": false}
using LinearAlgebra
function comatrix(A)
    C = zero(A)
    ### BEGIN SOLUTION
    for ij ∈ CartesianIndices(A)
        i, j = Tuple(ij)
        C[ij] = (-1)^(i+j)*det(A[1:end .!= i, 1:end .!= j])
    end
    ### END SOLUTION
    return C
end

# + nbgrader={"grade": true, "grade_id": "comatrix-test", "locked": true, "points": 4, "schema_version": 3, "solution": false, "task": false}
n = 10
A = rand(n, n)
C = comatrix(A)
display(C'A)
@assert C'A ≈ det(A)*LinearAlgebra.I

# + [markdown] jp-MarkdownHeadingCollapsed=true
# ### <font color='green'>Creating custom structures</font>
# -

# Julia is not an object-oriented language.
# Instead, it supports *multiple dispatch*,
# which is eloquently motivated in the [Julia documentation](https://docs.julialang.org/en/v1/manual/methods/)
#
# > The choice of which method to execute when a function is applied is called dispatch. Julia allows the dispatch process to choose which of a function's methods to call based on the number of arguments given, and on the types of all of the function's arguments. This is different from traditional object-oriented languages, where dispatch occurs based only on the first argument, which often has a special argument syntax, and is sometimes implied rather than explicitly written as an argument. Using all of a function's arguments to choose which method should be invoked, rather than just the first, is known as multiple dispatch. Multiple dispatch is particularly useful for mathematical code, where it makes little sense to artificially deem the operations to "belong" to one argument more than any of the others: does the addition operation in x + y belong to x any more than it does to y? The implementation of a mathematical operator generally depends on the types of all of its arguments. Even beyond mathematical operations, however, multiple dispatch ends up being a powerful and convenient paradigm for structuring and organizing programs.

# +
struct Dog
  name::String
end

struct Cat
  name::String
end

greet(a::Dog, b::Dog) = "$(a.name): Wouf wouf $(b.name), content de te revoir mon vieux"
greet(a::Dog, b::Cat) = "$(a.name): Wouf wouf wouf wouf !"
greet(a::Cat, b::Dog) = "$(a.name) s'éloigne discrètement de $(b.name)"
greet(a::Cat, b::Cat) = "$(a.name) commence à ronronner"

milou = Dog("Milou")
rantanplan = Dog("Rantanplan")
miaouss = Cat("Miaouss")
mrMoustache = Cat("MrMoustache")

@show greet(milou, rantanplan)
@show greet(milou, mrMoustache)
@show greet(mrMoustache, miaouss)
@show greet(mrMoustache, rantanplan);
# -

# To conclude this session,
# let us illustrate the power of multiple disptach using two other examples.
# First, let us create a specific type for the identity matrix,
# and redefine the operations `+` and `*` on matrices so that they work with this new type.

# +
struct IdentityMatrix end

# This is required to redefine + and *
import Base: +,*

function +(A::Matrix, Id::IdentityMatrix)
    result = copy(A)
    for i in 1:size(A, 1)
        result[i, i] += 1
    end
    return result
end

+(Id::IdentityMatrix, A::Matrix) = +(A, Id)
*(A::Matrix, Id::IdentityMatrix) = copy(A)
*(Id::IdentityMatrix, A::Matrix) = copy(A)

A = ones(3, 3)
Id = IdentityMatrix()

# We can now write operations such as
@show A
@show A + Id
@show Id + A
@show A * Id;
# -

# Note that we did not need to modify the internal implementation of matrices.
# Remark also that, with our implementation of the identity matrix,
# the product `A*I` is virtually free.
# Had the identity matrix been stored as a full matrix,
# a lot of computing power would have been wasted in the matrix multiplication.

# +
n = 5000
A = rand(n, n)
I1 = IdentityMatrix()
I2 = [(i == j)*1.0 for i in 1:n, j in 1:n]

@time A*I1;
@time A*I2;
# -

# #### <font color='orange'>[Exercise 6]</font>
# The aim of this exercise is to create a structure to represent dual numbers.
# Dual numbers are a mathematical construct that extends
# the real numbers by introducing a new element, often
# denoted by $\varepsilon$, such that $\varepsilon^2 = 0$.
# A dual number is composed of a real part and a dual part,
# where the dual part behaves like an infinitesimal quantity.
# They will be very useful later in the course,
# when we study automatic differentiation, a technique
# for efficiently computing derivatives numerically.
# Complete the structure definition below by adding the operations `-`, `*` and `/` to the dual numbers.

# + jupyter={"source_hidden": true} nbgrader={"grade": false, "grade_id": "dual", "locked": false, "schema_version": 3, "solution": true, "task": false}
import Base: +, -, *, /, show

struct Dual
    real::Float64
    dual::Float64
end

+(x::Dual, y::Dual) = Dual(x.real + y.real, x.dual + y.dual)
### BEGIN SOLUTION
*(x::Dual, y::Dual) = Dual(x.real*y.real, (x.dual*y.real + x.real*y.dual))
-(x::Dual, y::Dual) = Dual(x.real - y.real, x.dual - y.dual)
/(x::Dual, y::Dual) = Dual(x.real/y.real, (y.real*x.dual - x.real*y.dual)/y.real^2)
### END SOLUTION

show(io::IO, x::Dual) = print(io, x.real, x.dual ≥ 0 ? " + $(x.dual)" : " - $(-x.dual)", " ε")

# Example usage
a = Dual(3.0, 2.0)  # 3 + 2ε
b = Dual(1.0, 4.0)  # 1 + 4ε

@show a + b
@show a - b
@show a * b
@show a / b;

# + nbgrader={"grade": true, "grade_id": "dual-tests", "locked": true, "points": 4, "schema_version": 3, "solution": false, "task": false}
@assert a + b == Dual(4, 6)
@assert a - b == Dual(2, -2)
@assert a * b == Dual(3, 14)
@assert a / b == Dual(3, -10)

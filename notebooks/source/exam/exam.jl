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
#     display_name: Julia 1.10.4
#     language: julia
#     name: julia-1.10
# ---

# + [markdown] nbgrader={"grade": false, "grade_id": "cell-9546f1bdb893cdcd", "locked": true, "schema_version": 3, "solution": false, "task": false}
# ## Final exam
#
# - Ce notebook est à soumettre sur <a href="https://educnet.enpc.fr/mod/assign/view.php?id=65895">Educnet</a> avant 11h45.
#
# - L’examen comporte trois exercices indépendants. Dans chaque exercice les
#   cellules peuvent éventuellement dependre des cellules précèdentes.
#
# - Afin de faciliter l'évaluation de votre code,
#   ne pas changer les signatures des fonctions à implémenter.
#
# - La cellulle ci-dessous importe les bibliothèques utilisées dans ce notebook. Si une ou plusieurs d'entre elles manquent sur votre machine, vous êtes invités à les installer au préalable dans le gestionnaire de bibliothèques d'une console.

# + nbgrader={"grade": false, "grade_id": "cell-4d638a879ba6f86e", "locked": true, "schema_version": 3, "solution": false, "task": false}
using ForwardDiff
using LaTeXStrings
using LinearAlgebra
using Plots
using Polynomials
using Interpolations
using Random

Plots.default(titlefontsize=14,
              xlabelfontsize=12,
              ylabelfontsize=12,
              legendfontsize=12,
              xtickfontsize=12,
              ytickfontsize=12)
# -

# ## <font color='orange'>[Exercice I]</font> Une méthode de déflation pour calculer plusieurs valeurs propres dominantes

# Le but de cet exercice est de mettre en œuvre une méthode pour approximer plusieurs valeurs propres dominantes d'une matrice hermitienne ${\sf A} \in \mathbb{C}^{n\times n}$.
# Pour décrire la méthode,
# notons $ |\lambda_1|\geq |\lambda_2|\geq \dots \geq |\lambda_n|$ les valeurs propres de $\sf A$,
# et définissons la suite de matrices :
# $$ {\sf A}_1 := {\sf A},\qquad {\sf A}_{k+1} = {\sf A}_{k} - \lambda_{k} v_{k} v_{k}^*,$$
# où $v_{k}$ est un vecteur propre normalisé de $\mathsf A$ associé à la valeur propre $\lambda_k$.
# Il est facile de vérifier par récurrence que $(\lambda_k,v_k)$ est le couple propre dominant de ${\sf A}_k$.
# La méthode de déflation consiste simplement à définir la suite
# $$ \widetilde{\sf A}_1 := {\sf A},\quad \widetilde{\sf A}_{k+1} = \widetilde{\sf A}_{k} - \widetilde\lambda_{k} \widetilde v_{k} \widetilde v_{k}^*,$$
# où $(\widetilde\lambda_{k},\widetilde v_{k})$ est une approximation numérique du couple propre dominant de $ \widetilde{\sf A}_{k} $ calculée par l'itération de la puissance.
#
# 1. Implémenter une fonction `power_iter(A, x₀; ɛ=1e-12, maxiter=5000)` prenant comme arguments une matrice `A` de taille $n\times n$,
#    un vecteur initial de $n$ éléments `x₀`, un seuil de tolérance `ɛ` et un nombre maximal d'itérations `maxiter`,
#    et retournant un tuple `(λ, x)` contenant le résultat de l'itération de la puissance appliquée à `A`.

function power_iter(A, x₀; ε=1e-8, maxiter=100000)
    ### BEGIN SOLUTION
    v = x₀
    for k=1:maxiter
        v = A*v
        normalize!(v)
        λ = v'A*v

        if norm(A*v - λ*v) < ε
            return λ,v
        end

    end
    error("Power iteration failed to converge.")
    ### END SOLUTION
end;

@assert power_iter([1. 2.; 2. 1.], [1., 0.])[1] ≈ 3.
@assert power_iter([1. 0.; 0. .5], [1., 1.])[1] ≈ 1.
@assert [1, -1]'power_iter([1. 2.; 2. 1.], [1., 0.])[2] |> abs < 1e-6
@assert [1, 0]'power_iter([0. 0.; 0. 1.], [1., 1.])[2] |> abs < 1e-6
@assert [0, 1]'power_iter([1. 0.; 0. .5], [1., 1.])[2] |> abs < 1e-6

# 2. Implémenter une fonction `deflation_method(A, x₀, nev; ε=1e-12, maxiter=5000)` prenant comme arguments la matrice `A`, un vecteur initial `x₀`, le nombre de valeurs propres désirées `nev`, la tolérance `ɛ` et le nombre d'itérations `maxiter` pour chaque appel à l'itération de la puissance.
#    La fonction devra retourner un tuple `(λs, vs)` où `λs` est un vecteur de `nev` valeurs propres approximatives,
#    classées de la plus grande à la plus petite en valeur absolue,
#    et `vs` est une matrice dont les colonnes sont constituées des vecteurs propres associés.

function deflation_method(A, x₀, nev; ε=1e-8 ,maxiter=100000)
    n = size(A, 1)

    (nev > n) && error("$nev eigenvalues required for a $n×$n matrix")

    λs = zeros(eltype(A), nev)
    vs = zeros(eltype(A), n, nev)

    ### BEGIN SOLUTION
    for k=1:nev
        λ,v = power_iter(A,x₀ - vs*vs'x₀;ε=ε,maxiter=maxiter)
        λs[k] = λ
        vs[:,k] .= v
        A -= λ*v*v'
    end
    ### END SOLUTION

    return λs,vs
end;

# +
# Tests automatiques
N = 200
seed = 2024
A = randn(Xoshiro(seed), N, N)
A = (A+A')/sqrt(N)
x₀ = randn(N)

nev = 10

@time λs, us = deflation_method(A,x₀,nev)
@assert all(λs .≈ sort(eigvals(A), by=abs, rev=true)[1:nev])

# + [markdown] jp-MarkdownHeadingCollapsed=true nbgrader={"grade": false, "grade_id": "duplicate_id", "locked": true, "schema_version": 3, "solution": false, "task": false}
# ### <font color='orange'>[Exercice 2]</font> Intégration composite de Gauss-Legendre
#
# 1. Écrire une fonction `legendre(n)` qui retourne le polynôme de Legendre de degré $n$,
#    sous forme d'une structure `Polynomial` de la bibliothèque `Polynomials`.
#    Pour ce faire, vous pouvez utiliser la bibliothèque `Polynomials` et la formule de Rodrigues :
#    $$
#    L_n(x) = \frac{1}{2^n n!} \frac{\mathrm{d}^n}{\mathrm{d} x^n} \left(x^2 - 1\right)^n.
#    $$
#
#     <details>
#         <summary>
#             <em><font color='gray'>Indication (cliquer pour afficher)</font></em>
#         </summary>
#
#     - La fonction `factorial(n)` permet de calculer la factiorielle de `n`.
#
#     - La fonction `Polynomials.Polynomial` permet de créer un polynôme à partir de ses coefficients :
#       ```julia
#       p = Polynomial([1, 2, 3])  # p(x) = 1 + 2x + 3x²
#       ```
#     - La fonction `Polynomials.derivative` permet de calculer les dérivées d'un polynôme :
#       ```julia
#       dp = derivative(p)  # dp(x) = 2 + 6x
#       ddp = derivative(p, 2)  # ddp(x) = 6
#       ```
#     </details>
# -
function legendre(n)
    ### BEGIN SOLUTION
    p = Polynomial([-1, 0, 1])
    return 1 / (2^n * factorial(n)) * derivative(p^n, n)
    ### END SOLUTION
end;

# 2. Écrire une fonction `get_nodes_and_weights(n)` qui calcule,
#    sans utiliser d'autres bibliothèques logicielles que celles importées au début du notebook,
#    les nœuds $(x_i)_{i \in \{1, \dots, n\}}$ et poids $(w_i)_{i \in \{1, \dots, n\}}$ de la quadrature de Gauss-Legendre avec $n$ nœuds.
#    Pour rappel, les nœuds et poids doivent être tels que l'approximation
#    $$
#    \int_{-1}^{1} f(x) \, \mathrm d x
#    \approx \sum_{i=1}^{n} w_i f(x_i)
#    $$
#    soit exacte pour tout polynôme $f$ de degré au plus $2n-1$.
#    <details>
#        <summary>
#            <em><font color='gray'>Indication (cliquer pour afficher)</font></em>
#        </summary>
#
#    - On rappelle que les nœuds d'intégration sont donnés par les racines du polynôme de Legendre de degré `n`.
#      Ces racines peuvent être calculées par la fonction `roots` de la biblothèque `Polynomials.jl`.
#
#    - Pour construire les polynômes de Lagrange en vue de calculer les poids,
#      il peut être utile d'utiliser les fonctions `fromroots` et `integrate` de la biblothèque `Polynomials.jl`.
#
#      ```julia
#          p = fromroots([1., 2.])  # Constructs (x - 1)(x - 2) = x² - 3x + 2
#          q = integrate(p)  # q = x^3/3 - 3x^2/2 + 2x
#      ```
#    </details>

function get_nodes_and_weights(n)
    ### BEGIN SOLUTION
    nodes = sort(roots(legendre(n)))
    weights = zero(nodes)
    for i in 1:n
        ℓ = fromroots(nodes[1:end .!= i])
        ℓ = ℓ / ℓ(nodes[i])
        weights[i] = integrate(ℓ, -1, 1)
    end
    ### END SOLUTION
    return nodes, weights
end;

@assert get_nodes_and_weights(5) |> length == 2
@assert get_nodes_and_weights(5)[1] |> length == 5
@assert get_nodes_and_weights(5)[2] |> length == 5
@assert get_nodes_and_weights(1)[1] ≈ [0.]
@assert get_nodes_and_weights(1)[2] ≈ [2.0]
@assert get_nodes_and_weights(3)[1] .|> legendre(3) |> abs∘sum < 1e-10
@assert get_nodes_and_weights(5)[1] .|> legendre(5) |> abs∘sum < 1e-10
@assert get_nodes_and_weights(5)[2] |> sum ≈ 2

# 3. Écrire une fonction `composite_gauss_legendre(u, a, b, n, N)` qui renvoie une approximation de l'intégrale
#     $$
#     \int_{a}^{b} u(x) \, \mathrm{d} x
#     $$
#     obtenue en partitionnant l'intervalle d'intégration $[a, b]$ en $N$ cellules de même taille,
#     et en appliquant la quadrature de Gauss-Legendre avec $n$ nœuds dans chaque cellule.

function composite_gauss_legendre(u, a, b, n, N)
    ### BEGIN SOLUTION
    h = (b-a)/N
    X = LinRange(a, b, N + 1)
    z, w = get_nodes_and_weights(n)
    result = 0.
    for i in 1:N
        nodes = X[i] + h/2 .+ z*h/2
        result += h/2 * w'u.(nodes)
    end
    return result
    ### END SOLUTION
end;

_short(f, n, N) = composite_gauss_legendre(f, 0, 1, n, N)
for d in 1:9
    @assert _short(x -> x^d, 5, 1) ≈ 1/(d+1)
    @assert _short(x -> x^d, 5, 2) ≈ 1/(d+1)
    @assert _short(x -> x^d, 5, 3) ≈ 1/(d+1)
end
@assert !(_short(x -> x^10, 2, 1) ≈ 1/11)
@assert !(_short(x -> x^10, 2, 2) ≈ 1/11)
@assert _short(x -> x^10, 5, 200) ≈ 1/11
@assert _short(x -> exp(x), 5, 200) ≈ ℯ - 1

# 4. Considérons le cas particulier où $u(x) = \cos(x)$, $a = -1$ et $b = 1$,
#     et définissons l'erreur d'intégration,
#     vue comme une fonction de $N$ où $n$ est un paramètre fixé,
#     par la formule
#     $$
#     E_{n}(N) = \lvert \widehat I_{n, N} - I_{\rm exact} \rvert.
#     $$
#     Dans cette équation,
#     $I_{\rm exact}$ est la valeur exacte de l'intégrale
#     tandis que $\widehat I_{n, N}$ est son approximation par la règle de Gauss-Legendre composite.
#     Il est demandé
#
#     - d'estimer, pour chaque valeur de $n \in \{1, 2, 3\}$,
#       l'ordre de convergence de la quadrature de Gauss-Legendre composite par rapport à $N$,
#       c'est-à-dire de trouver $\beta = \beta(n)$ tel que
#       $$
#       E_n(N) \propto C N^{-\beta}.
#       $$
#
#     - d'illustrer sur un même graphique,
#       à l'aide de la fonction `Plots.scatter`,
#       les fonctions $E_1, E_2, E_3$,
#       pour des valeurs de $N$ variant de 1 à 40.
#       Utiliser l'échelle logarithmique pour les deux axes,
#       et inclure l'ordre de convergence `β` trouvé au point précédent dans la légende,
#       en passant par exemple à la fonction `scatter` l'argument `label="n=$n, β=$β"`.

# +
# Function to integrate
uᵢ(x) = cos(x)

# Integration interval
a, b = -1, 1

# Exact value of the integral
I_exact = 2sin(1)

### BEGIN SOLUTION
# Number of nodes
ns = [1, 2, 3]

# Number of cells
N = 1:40

p = plot(title="Convergence of Gauss Legendre quadrature", legend=:bottomleft,
         xticks=([1, 5, 10, 20, 30], ["1", "5", "10", "20", "30"]))
for n in ns
    errors = composite_gauss_legendre.(uᵢ, a, b, n, N) .- I_exact
    polyfit = fit(log.(N), log.(abs.(errors)), 1)
    β = round(- polyfit[1], digits=2)
    scatter!(N, abs.(errors), label="n=$n, β=$β", scale=:log10)
    xlabel!(L"N")
    ylabel!(L"|I - \widehat I_{n,N}|")
end
p
### END SOLUTION
# -

# ### <font color='orange'>[Exercice 3]</font> Trajectoire d'une masse ponctuelle sur un rail

# On considère dans cet exercice un solide de masse $m$ modélisé comme un point susceptible de se déplacer sans frottement le long d'un arc paramétré $\mathbf{M}(x)=\left(x, y(x)\right)^T$ de classe $\mathcal{C}^2$ dans le plan $xy$. On introduit les vecteurs tangent $\mathbf{t}(x)$ et normal $\mathbf{n}(x)$ à l'arc par
#
# $$
# \mathbf{t}(x)=\frac{\mathbf{M}'(x)}{\lVert \mathbf{M}'(x) \rVert}=\frac{(1, y'(x))^T}{\lVert \mathbf{M}'(x) \rVert}
# \quad ; \quad
# \mathbf{n}(x)=\frac{(-y'(x), 1)^T}{\lVert \mathbf{M}'(x) \rVert}
# \quad \textrm{avec} \quad
# \lVert \mathbf{M}'(x) \rVert=\sqrt{1+y'(x)^2}
# $$
#
# Par convention on notera dans la suite $y'$ la dérivation de $y$ par rapport à $x$ et par $\dot{y}$ la dérivation par rapport au temps, ce qui entraîne que $\dot{y}=y' \dot{x}$.
#
# Le solide est soumis à la gravité $\mathbf{P}=-mg\mathbf{e}_y$ et à la réaction du support $\mathbf{R}=R_n\mathbf{n}$ (pas de frottement). On montre alors par projection sur $\mathbf{t}(x)$ du principe fondamental de la dynamique que le mouvement est régi par la résolution d'une équation différentielle d'ordre 2 sur $x(t)$
#
# $$
# \ddot{x}=\frac{-g y'(x) - y'(x) y''(x) \dot{x}^2}{1+y'(x)^2}
# \quad ; \quad
# x(0)=x_0
# \quad ; \quad
# \dot{x}(0)=\dot{x}_0
# \tag{EDOx}
# $$
# <a id="EDOx"></a>
# Cette équation différentielle se ramène classiquement à une équation vectorielle d'ordre 1 d'inconnue $\mathbf{u}=(x, \dot{x})^T$ sous la forme suivante
#
# $$
# \dot{\mathbf{u}}=\mathbf{f}(\mathbf{u})=
# \begin{pmatrix}
# u_2\\
# \frac{-g y'(u_1) - y'(u_1) y''(u_1) u_2^2}{1+y'(u_1)^2}
# \end{pmatrix}
# \quad ; \quad
# \mathbf{u}(0)=\begin{pmatrix} x_0 \\ \dot{x}_0 \end{pmatrix}
# \tag{EDOu}
# $$
# On considère que la fonction $y(x)$ est donnée par
# $$
# y_{\alpha}(x) = \frac{e^{-\alpha x} - e^{-\alpha}}{e^{\alpha} - e^{-\alpha}},
# $$
# où $α$ est un paramètre strictement positif.
# Noter que $y_{\alpha}(-1) = 1$ et $y_{\alpha}(1) = 0$.
# La fonction $y_\alpha$ et ses dérivées sont implémentées ci-dessous sous la forme de fonctions de fonctions, i.e. `y(α)` renvoie la fonction $y_\alpha$ de $x$.

# +
const g = 9.81 ;
const m = 0.03 ;

y(α) = x -> (exp(-α * x) - exp(-α)) / (exp(α) - exp(-α))
dy(α) = x -> -α*exp(-α * x) / (exp(α) - exp(-α))
d2y(α) = x -> α^2*exp(-α * x) / (exp(α) - exp(-α))
# -

# 1. Tracer les courbes correspondant à plusieurs valeurs de $\alpha \in \{1, \dotsc, 10\}$.

### BEGIN SOLUTION
pl = plot(y(0),-1,1, label="")
for α in LinRange(1., 10., 10)
    plot!(pl, y(α), label="")
end
pl
### END SOLUTION

# 2. Écrire la fonction $f$ sous forme d'une fonction Julia `f(u, α)`.

function f(u, α)
    ### BEGIN SOLUTION
    x, dx = u
    yp, ypp = x .|> [dy(α), d2y(α)]
    return [dx, -(g*yp + yp * ypp * dx^2)/(1 + yp^2)]
    ### END SOLUTION
end

# 3. Écrire une fonction `rk4(tₙ, uₙ, f, Δ)` implémentant un pas de temps de taille $\Delta$ de la méthode de Runge-Kutta d'ordre 4 pour une équation différentielle générique de la forme $u' = f(u)$.
#    Cette méthode est basée sur l'itération suivante:
#    $$
#       \mathbf u_{n+1} = \mathbf u_n + \frac{\Delta}{6}\left(\mathbf k_1 + 2\mathbf k_2 + 2\mathbf k_3 + \mathbf k_4 \right),
#    $$
#    où
#    \begin{align*}
#    \mathbf k_1 &= \ \mathbf f(\mathbf u_n), \\
#    \mathbf k_2 &= \ \mathbf f\!\left(\mathbf u_n + \frac{\Delta}{2} \mathbf k_1\right), \\
#    \mathbf k_3 &= \ \mathbf f\!\left(\mathbf u_n + \frac{\Delta}{2} \mathbf k_2\right), \\
#    \mathbf k_4 &= \ \mathbf f\!\left(\mathbf u_n + \Delta\mathbf k_3\right).
#    \end{align*}
#    La fonction devra renvoyer $\mathbf u_{n+1}$. À noter que l'argument générique `f` est considéré ici comme une fonction d'une unique variable vectorielle `u`.

function rk4(uₙ, Δ, f)
    ### BEGIN SOLUTION
    k₁ = f(uₙ)
    k₂ = f(uₙ + Δ/2 * k₁)
    k₃ = f(uₙ + Δ/2 * k₂)
    k₄ = f(uₙ + Δ   * k₃)
    return uₙ + Δ/6 * (k₁ + 2k₂ + 2k₃ + k₄)
    ### END SOLUTION
end

@assert rk4([0.], 1., u -> [1.]) ≈ [1.0]
@assert rk4([3.,12.], 1., u -> u .* [2.,3.]) ≈ [21.,196.5]

# 4. Écrire une fonction `solve_ode(u₀, Δ, α)` pour une condition initiale `u₀` et un paramètre `α`,
#    en utilisant la méthode de Runge-Kutta d'ordre 4 avec pas de temps fixe `Δ`.
#    Votre fonction devra renvoyer un vecteur de temps `T` et un vecteur de vecteurs `U` contenant la solution à ces temps.
#    On demande d'interrompre l'intégration numérique dès que la valeur de la coordonnée $x$ sera devenue supérieure ou égale à 1;
#    il faudra donc que `U[end-1][1]` soit inférieur à 1 et `U[end][1]` soit supérieur ou égal à 1.

function solve_ode(u₀, Δ, α)
    fα(u) = f(u, α)
    U = [u₀ .+ 0.0 * α] # pour éviter les problèmes de type
    T = [0.]
    ### BEGIN SOLUTION
    while U[end][1] < 1.
        # println(length(T))
        push!(U, rk4(U[end], Δ, fα))
        push!(T, T[end]+Δ)
    end
    ### END SOLUTION
    return T, U
end

# 5. Écrire une fonction `final_time(α)`,
#    qui retourne une approximation du temps mis par le solide pour atteindre la position $x = 1$.
#    Pour ce faire, résoudre l'équation différentielle avec un pas de temps $Δ = 0.01$,
#    et estimer le temps requis par interpolation linéaire sur le dernier pas de temps,
#    durant lequel la coordonnée $x$ du solide passe au delà de 1.

function final_time(α)
    Δ = .01
    ### BEGIN SOLUTION
    T, U = solve_ode([-1., 0.], Δ, α)
    t₁, t₂, δx₁, δx₂ = T[end-1], T[end], U[end-1][1]-1, U[end][1]-1
    return (t₁*δx₂ - t₂*δx₁)/(δx₂-δx₁)
    ### END SOLUTION
end


# 6. Faire un plot du temps final en fonction de `α`.
#    Estimer graphiquement, par exemple à l'aide de la fonction `vline`,
#    la valeur de `α` permettant de minimiser le temps final.

### BEGIN SOLUTION
plot(LinRange(0.001, 20, 1000), final_time)
vline!([2.55])
### END SOLUTION

# 7. On se propose ici de déterminer la ou les valeurs de $α$ permettant d'atteindre le point final en un temps donné par une méthode de Newton-Raphson s'appuyant sur les nombres duaux (afin d'être en mesure de dériver une fonction de $\alpha$). Écrire une fonction `newton_raphson_dual(x, f, maxiter=100; ε = 1e-12)` renvoyant une racine de la fonction `f` en partant d'un point initial `x`.
#
#     <details>
#         <summary>
#             <em><font color='gray'> Indication (cliquer pour afficher)</font></em>
#         </summary>
#
#     On rappelle que l'on peut obtenir simultanément la valeur et la dérivée d'une fonction `f` en `x` par la méthode suivante
#
#     ```julia
#       y = f(ForwardDiff.Dual(x, 1.))
#       fx = y.value # renvoie f(x)
#       dfx = y.partials[1] # renvoie f'(x)
#     ```
#
#     **Pour éviter des problèmes de confusion en manipulant plusieurs ordres de dérivation, il est préférable d'utiliser `f(ForwardDiff.Dual(x, 1.))` plutôt que `f(x+ForwardDiff.Dual(0., 1.))`.**
#
#     </details>

function newton_raphson_dual(x, f, maxiter=100; ε = 1e-12)
    ### BEGIN SOLUTION
    for i in 1:maxiter
        y = f(ForwardDiff.Dual(x, 1.))
        x -= y.value/y.partials[1]
        norm(f(x)) < ε && return x
    end
    error("Failed to converge!")
    ### END SOLUTION
end;

# 8. Déterminer la ou les valeurs de $\alpha$ permettant d'atteindre le point final en $t=0.85s$.

### BEGIN SOLUTION
α₁ = newton_raphson_dual(1., x -> final_time(x) - 0.85)
α₂ = newton_raphson_dual(10., x -> final_time(x) - 0.85)
@show α₁, α₂ ;
### END SOLUTION

# 9. Ayant constaté graphiquement l'existence d'un paramètre `α` optimal,
#    calculer précisément ce paramètre en utilisant la méthode de votre choix.

# +
### BEGIN SOLUTION
function bisection(f, a, b; δ = 1e-10)
    @assert f(a) * f(b) ≤ 0
    while abs(b - a) ≥ δ
        x = (a + b) / 2
        a, b = f(a) * f(x) ≤ 0 ? [a, x] : [x, b]
    end
    return (a + b) / 2
end
α_opt = bisection(x -> final_time(ForwardDiff.Dual(x, 1.)).partials[1], 1., 3.)

# ou

α_opt = newton_raphson_dual(1., x -> final_time(ForwardDiff.Dual(x, 1.)).partials[1])
### END SOLUTION
# -

@assert abs(α_opt-2.677)<1e-3

# 10. Écrire une fonction d'animation `animate_list_sol(listα, Δ)` permettant de calculer des solutions correspondant à plusieurs valeurs de $\alpha$ dans `listα` pour un même pas de temps `Δ`,
#     et de superposer les trajectoires.
#
#     Bien prêter attention au fait que les vecteurs temps récupérés lors des différentes simulations ne sont pas de même longueur mais ont le même pas entre éléments consécutifs.
#
#     En notant `α_opt` la valeur définie à la question précédente, on pourra appliquer l'animation à la liste `[α_opt/100, α_opt/2, α_opt, 2α_opt, 4α_opt]` et constater que c'est bien la trajectoire liée à `α_opt` qui arrive en premier.

# +
function animate_list_sol(listα, Δ)
    ### BEGIN SOLUTION
    TUs = [solve_ode([-1, 0], Δ, α) for α in listα]
    Ts, Us = first.(TUs), last.(TUs)
    T = argmax(length, Ts)
    anim = @animate for i in eachindex(T)
        t = T[i]
        plot(title="t = $(round(t, digits=3))", xlims=(-1,1), ylims=(-0.5,1.2),
             aspect_ratio=:equal, xaxis=false, yaxis=false, grid=false, ticks=false)
        for j in eachindex(listα)
            α = listα[j]
            u = Us[j][min(i,length(Us[j]))]
            plot!(y(α), aspect_ratio=:equal, label=nothing, color=:blue)
            scatter!([u[1]], [y(α)(u[1])], label=nothing)
        end
    end
    return anim
    ### END SOLUTION
end

listα = [α_opt/100, α_opt/2, α_opt, 2α_opt, 4α_opt]
gif(animate_list_sol(listα, 0.01), fps=10, show_msg=false)

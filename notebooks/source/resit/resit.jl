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
# ## Examen de rattrapage
#
# - Ce notebook est à soumettre sur <a href="https://educnet.enpc.fr/mod/assign/view.php?id=66064">Educnet</a> avant 16h15.
#
# - L’examen comporte trois exercices indépendants. Dans chaque exercice les
#   cellules peuvent éventuellement dependre des cellules précèdentes.
#
# - Afin de faciliter l'évaluation de votre code,
#   ne pas changer les signatures des fonctions à implémenter.
#
# - La cellulle ci-dessous importe les bibliothèques utilisées dans ce notebook. Si une ou plusieurs d'entre elles manquent sur votre machine, vous êtes invités à les installer au préalable.

# + nbgrader={"grade": false, "grade_id": "cell-4d638a879ba6f86e", "locked": true, "schema_version": 3, "solution": false, "task": false}
using ForwardDiff
using LaTeXStrings
using LinearAlgebra
using Plots
using Polynomials
# -

# ## <font color='orange'>[Exercice 1]</font> Interpolation pour l'équation des poutres

# Le but de cet exercice est d'explorer une application de l'interpolation polynomiale à la résolution d'une équation différentielle.
# Plus précisément, nous allons mettre en œuvre une méthode numérique pour résoudre approximativement l'équation des poutres d'Euler-Bernoulli avec des conditions aux limites de Dirichlet homogènes :
#
# $$
# u \in C^4([0,1]), \qquad \left\{ \begin{aligned} u''''(x) &= \varphi(x) \qquad \forall\, x \in (0,1),\\
# u(0) &= u'(0) = u'(1) = u(1) = 0, \end{aligned} \right.
# $$
#
# où $\varphi(x) = (2\pi)^4\cos(2\pi x)$ est une charge transverse appliquée à la poutre.
# Afin de résoudre l'équation numériquement,
# nous allons approximer le terme de droite $\varphi$ par un polynôme interpolant $\widehat \varphi$,
# puis résoudre l'équation exactement avec $\widehat \varphi$ au lieu de $\varphi$.
#
# 0. Commençons par écrire une fonction `fit_values_and_slopes(u₀, up₀, u₁, up₁)` qui retourne l'unique polynôme de degré 3 tel que
#    $$
#    p(0) = u_0, \qquad p'(0) = up_0, \qquad p(1) = u_1, \qquad p'(1) = up_1.
#    $$


# +
function fit_values_and_slopes(u₀, up₀, u₁, up₁)
    # We look for polynomials p(x) = a₀ + a₁ x + a₂ x² + a₃ x³
    A = [1 0 0 0; 0 1 0 0; 1 1 1 1; 0 1 2 3]
    α = A\[u₀; up₀; u₁; up₁]
    return Polynomial(α)
end

# Test code
p = fit_values_and_slopes(-1, -1, 1, 1)
plot(p, xlims=(0, 1))
# -

# 1. Écrire une fonction `approx(n)` implémentant l'approche décrite ci-dessus pour résoudre l'EDP.
#    La fonction devra retourner une approximation polynomiale de la solution basée sur une interpolation de **degré** $n$ du membre de droite à des points équidistants entre 0 et 1 compris.
#
#    <details>
#       <summary>
#          <em><font color='gray'>Indication (cliquer pour afficher)</font></em>
#       </summary>
#
#    - Utiliser la fonction `fit` de la bibliothèque `Polynomials.jl` pour obtenir le polynôme interpolateur :
#
#        ```julia
#            p = fit(x, y)
#        ```
#
#        où `x` sont les nœuds d'interpolation et `y` sont les valeurs correspondantes de la fonction à interpoler.
#
#    - Pour calculer la solution analytique lorsque le membre de droite est un polynôme,
#      on pourra remarquer que toutes les solutions sont des polynômes,
#      et que sans les conditions aux bords, la solution est unique modulo un polynôme cubique.
#
#    - La fonction `integrate` de la bibliothèque `Polynomials.jl` permet de calculer une primitive d'un polynôme :
#        ```julia
#            P = integrate(p)
#        ```
#
#    - Utiliser le format `BigFloat` pour limiter les erreurs d'arrondi.
#      En particulier, la fonction `LinRange{BigFloat}(a, b, N)` permet de créer un vecteur de `N` nombres au format `BigFloat` également distribués entre `a` et `b` inclus.
#      ```julia
#          X = LinRange{BigFloat}(0, 1, n + 1)
#      ```
#    </details>

# + nbgrader={"grade": false, "grade_id": "interp_pde", "locked": false, "schema_version": 3, "solution": true, "task": false}
# Right-hand side
φ(x) = (2*big(π))^4 * cospi(2*x)

# Exact solution (for comparison purposes)
uₑ(x) = cospi(2*x) - 1

function approx(n)
    X = LinRange{BigFloat}(0, 1, n + 1)
    ### BEGIN SOLUTION
    Y = φ.(X)
    p = fit(X, Y)
    uh = integrate(integrate(integrate(integrate(p))))
    ∂uh = derivative(uh)
    uh -= fit_values_and_slopes(uh(0), ∂uh(0), uh(1), ∂uh(1))
    return uh
    ### END SOLUTION
end;

# + nbgrader={"grade": true, "grade_id": "cell-e79ffe6285cd8198", "locked": true, "points": 2, "schema_version": 3, "solution": false, "task": false}
@assert typeof(approx(3)) == Polynomial{BigFloat, :x}
@assert degree(approx(3)) == 7
@assert approx(3)(0) ≈ 0
@assert approx(3)(1) ≈ 0
@assert derivative(approx(3))(0) ≈ 0
@assert derivative(approx(3))(1) ≈ 0
# -

# 2. Écrire une fonction `estimate_error(n)` qui approxime l'erreur,
#    en norme $L^\infty$,
#    entre la solution approchée par l'approche ci-dessus et la solution exacte.
#    La solution exacte est donnée par
#    $$
#       u_e(x) = \cos(2\pi x) - 1.
#    $$

# + nbgrader={"grade": false, "grade_id": "interp_error", "locked": false, "schema_version": 3, "solution": true, "task": false}
# Exact solution (for comparison purposes)
uₑ(x) = cospi(2*x) - 1

function estimate_error(n)
    ### BEGIN SOLUTION
    un = approx(n)
    x_fine = LinRange{BigFloat}(0, 1, 1000)
    un_fine, u_fine = un.(x_fine), uₑ.(x_fine)
    return maximum(abs.(u_fine - un_fine))
    ### END SOLUTION
end;

# + nbgrader={"grade": true, "grade_id": "cell-9764a7b6ad0837c3", "locked": true, "points": 2, "schema_version": 3, "solution": false, "task": false}
@assert estimate_error(2) > 0
@assert estimate_error(20) < 1e-3
@assert estimate_error(20) > 1e-20
@assert estimate_error(40) < 1e-12
# -

# 3. Tracer un graphique de l'erreur en fonction de $n$ pour $n$ allant de 5 à 30.
#    Utiliser l'échelle par défaut pour l'axe des abcisses et une échelle logarithmique pour l'axe des ordonnées.

# + nbgrader={"grade": true, "grade_id": "interp_plot", "locked": false, "points": 1, "schema_version": 3, "solution": true, "task": false}
# ### BEGIN SOLUTION
ns = 5:50
errors = estimate_error.(ns)
plot(ns, errors, marker = :circle, label=L"$L^{\infty}$ Error")
plot!(yaxis=:log, lw=2)
# ### END SOLUTION
# -

# ### <font color='orange'>[Exercice 2]</font> Différence finies pour Poisson non-linéaire
# On s'intéresse dans cet exercice à la résolution de l'équation de Poisson non-linéaire suivante :
# $$
# \tag{Poisson}
# - \frac{\mathrm d} {\mathrm d x} \left( \kappa(u, \alpha) \frac{\mathrm d u} {\mathrm d x} \right) = 1,
# \qquad u(0) = u(1) = 0, \qquad \kappa(u, α) := 1 + α u^2,
# $$
# où $\alpha \in \mathbb R_{\geq 0}$ un paramètre positif.
# Cette équation modélise la température à l'équilibre dans une barre unidimensionnelle dont les extrémités sont maintenues à température nulle,
# en présence d'une source de chaleur uniforme.
# La non-linéarité de cette équation provient du fait que la conductivité thermique $\kappa$ dépend explicitement de la solution $u$.
# Pour résoudre numériquement cette équation,
# on utilisera la méthode des différences finies.
# Soit une grille d'abcisses équidistantes
# $$
# 0 = x_0 < x_1 < \dots < x_{N-1} < x_N = 1, \qquad x_i = \frac{i}{N},
# $$
# et soit $(u_i)_{i \in \{0, \dots N\}}$ les valeurs de la solution à ces points.
# La condition aux limites implique que $u_0 = u_N = 0$,
# et il ne sera donc pas nécessaire d'inclure ces valeurs dans les inconnues.
# Aux points intérieurs, on utilise l'approximation
# $$
# - \frac{\mathrm d} {\mathrm d x} \left( \kappa(u, α) \frac{\mathrm d u} {\mathrm d x} \right) (x_i)
# \approx \frac{F_{i+\frac{1}{2}} - F_{i-\frac{1}{2}}}{\Delta x},
# \qquad
# i \in \{1, \dotsc, N-1\},
# \tag{Approx}
# $$
# où $\Delta x = 1/N$ et, pour $k \in \{0, \dotsc, N-1\}$,
# $$
# F_{k + \frac{1}{2}} := - \kappa_{k+\frac{1}{2}} \frac{u_{k+1} - u_{k}}{\Delta x}, \qquad
# \kappa_{k + \frac{1}{2}} := \frac{\kappa(u_{i}, α) + \kappa(u_{i+1}, α)}{2}.
# $$
# Le terme $F_{k+ 1/2}$ est une approximation du flux de chaleur vers la droite,
# c'est à dire de la fonction qui à $x \in [0, 1]$ associe $- \kappa\bigl(u(x), α\bigr) u'(x)$,
# au point $x_k + \frac{\Delta x}{2}$.
# En substituant l'approximation <a>(Approx)</a> dans <a>(PDE)</a> et en réarrangant les termes,
# on obtient le système suivant :
# $$
# \forall i \in \{1, \dotsc, N-1\}, \qquad
# \frac{1}{\Delta x^2} \left( - \kappa_{i-\frac{1}{2}} u_{i-1} + \Bigl(\kappa_{i-\frac{1}{2}} + \kappa_{i+\frac{1}{2}}\Bigr) u_{i}  - \kappa_{i+\frac{1}{2}} u_{i+1} \right) = 1,
# $$
# qui constitue un système non-linéaire de $N-1$ équations en les inconnues $\mathbf u = (u_1, \dotsc, u_{N-1})^\top$.
# Notons que les coefficients $(\kappa_{k+1/2})_{k \in \{0, \dots, N-1\}}$ dépendent de manière non-linéaire de $\mathbf u$ et $α$,
# mais cette dépendance est omise dans la notation par souci de concision.
# Ce système peut être réécrit sous forme matricielle comme suit :
# $$
# \tag{NonLin}
# \frac{\mathsf A_N(\mathbf u, α)}{Δx^2}
# \begin{pmatrix}
#     u_1 \\
#     u_2 \\
#     u_3 \\
#     \vdots \\
#     u_{N-2} \\
#     u_{N-1}
# \end{pmatrix}
# =
# \begin{pmatrix}
#     1 \\
#     1 \\
#     1 \\
#     \vdots \\
#     1 \\
#     1
# \end{pmatrix}
# =: \mathbf b
# $$
# où $\mathsf A_N(\mathbf u, α) \in \mathbb R^{(N-1)\times (N-1)}$ est la matrice suivante (dimension <font color=red>N-1 x N-1</font>):
# $$
# A_N(\mathbf u, α) =
# \begin{pmatrix}
#     \kappa_{1/2} + \kappa_{3/2} & - \kappa_{3/2} \\
#     -\kappa_{3/2} & \kappa_{3/2} + \kappa_{5/2}  & - \kappa_{5/2} \\
#        & -\kappa_{5/2} & \kappa_{5/2} + \kappa_{7/2} & - \kappa_{7/2} \\
#        &    & \ddots & \ddots & \ddots & \\
#        &    &        &  -\kappa_{N-5/2}   & \kappa_{N-5/2} + \kappa_{N-3/2}  & -\kappa_{N-3/2} \\
#        &    &        &     & -\kappa_{N-3/2} & \kappa_{N-3/2} + \kappa_{N-1/2} \\
# \end{pmatrix}.
# $$
#
# 0. On fournit ci-dessous une fonction `build_A(N, u, α)` permettant de construire la matrice $\mathsf A_N(\mathbf u, \alpha)$.
#    <details>
#         <summary>
#             <em><font color='gray'>Indication (cliquer pour afficher)</font></em>
#         </summary>
#
#     En vue de pouvoir dans la suite dériver la matrice $\mathsf A_N$ par rapport à $α$ en utilisant la différentiation automatique,
#     il est crucial que la matrice $\mathbf A$ puisse contenir des nombres duaux.
#     C'est pourquoi on utilise la commande
#     ```julia
#     A = zeros(typeof(α), N+1, N+1)
#     ```
#     pour initialiser la matrice `A`.
#     De cette manière, les élément de la matrice sont du même type que le paramètre $α$.
#     </details>

# + nbgrader={"grade": false, "grade_id": "cell-7cc040866d32d6b0", "locked": true, "schema_version": 3, "solution": false, "task": false}
κ(u, α) = 1 + α^2 * u^2

function build_A(N, u, α)
    @assert length(u) == N-1
    A = zeros(typeof(α), N+1, N+1)
    v = [0.; u; 0.]
    for i in 2:N
        A[i, i-1] =  - (κ(v[i-1], α) + κ(v[i], α))  / 2
        A[i, i] = (κ(v[i-1], α) + 2κ(v[i], α) + κ(v[i+1], α)) / 2
        A[i, i+1] = - (κ(v[i], α) + κ(v[i+1], α)) / 2
    end
    return A[2:end-1, 2:end-1]
end

@assert begin N = 20; size(build_A(N, zeros(N-1), 0.)) == (N - 1, N - 1) end
@assert begin N = 20; build_A(N, zeros(N-1), 0.) == SymTridiagonal(fill(2, N-1), fill(-1, N-2)) end
# -

# 1. Pour résoudre le système non-linéaire <a>(NonLin)</a>,
#    on se propose d'utiliser une méthode de point fixe basée sur l'itération
#    $$
#    \frac{\mathsf A_{N}(\mathbf u^{n}, α)}{Δ x^2}  \mathbf u^{n+1} = \mathbf b.
#    $$
#    Cette itération permet de générer,
#    à partir d'une approximation courante de la solution $\mathbf u^n$,
#    une nouvelle approximation $\mathbf u^{n+1}$.
#    Écrire une fonction `solve_nonlin(N, α; maxiter=500, ε=1e-12)` permettant de calculer une solution approximative de <a>(NonLin)</a> par cette approche,
#    en utilisant comme critère d'arrêt que.
#    $$
#    \left\| \frac{\mathsf A_N(\mathbf u, α)}{Δx^2} \mathbf u - \mathbf b \right\| \leq ε
#    $$
#    La fonction devra renvoyer uniquement le résultat final $\mathbf u = (u_1, \dotsc, u_{N-1})^\top$ de l'itération,
#    ou `nothing` si une solution n'a pas été trouvée après `maxiter` itérations.
#    <details>
#         <summary>
#             <em><font color='gray'>Indication (cliquer pour afficher)</font></em>
#         </summary>
#
#     Afin de pouvoir dans la suite dériver la solution par rapport au paramètre $α$ en utilisant la différentiation automaiique,
#     il est important que les éléments `u` puissent contenir des nombres duaux.
#     Il est donc recommandé d'initialiser `u` comme suit :
#     ```julia
#       u = zeros(typeof(α), N-1)
#     ```
#     </details>

# + nbgrader={"grade": false, "grade_id": "solve_nonlin", "locked": false, "schema_version": 3, "solution": true, "task": false}
function solve_nonlin(N, α; maxiter=1000, ε=1e-10)
    b = ones(N-1)
    u = zeros(typeof(α), N-1)
    A = N^2 * build_A(N, u, α)
    ### BEGIN SOLUTION
    for i in 1:maxiter
        u = A\b
        A = N^2*build_A(N, u, α)
        norm(A*u - b) ≤ ε && return u
    end
    return nothing
    ### END SOLUTION
end;

# + nbgrader={"grade": true, "grade_id": "cell-68fd71f3ebf99d3e", "locked": true, "points": 2, "schema_version": 3, "solution": false, "task": false}
@assert size(solve_nonlin(10, 1.)) == (9,)
@assert solve_nonlin(10, 0; ε = 1e-100) == nothing
@assert begin N, α = 20, 1.; u = solve_nonlin(N, α); A_ = build_A(N, u, α); norm(N^2*A_*u - ones(N-1)) ≤ 1e-5 end
# -

# 2. Écrire une fonction `solve_poisson(N, α)` permettant de résoudre approximativement <a>(Poisson)</a>,
#    en faisant appel à la fonction `solve_nonlin` définie au point précédent.
#    Plus précisément, la fonction devra renvoyer les vecteurs $(x_0, x_1, \dotsc, x_N)^\top$ et $(u_0, u_1, \dotsc, u_N)^\top$,
#    où les **points aux limites sont inclus**.
#    <details>
#         <summary>
#             <em><font color='gray'>Hint (click to display)</font></em>
#         </summary>
#
#     Pour ajouter des éléments au début et à la fin d'un vecteur,
#     on pourra utiliser la concaténation verticale grâce à `;`
#
#     ```julia
#         v1 = [1., 2., 3.]
#         v2 = [0.; v1; 0.]  # v2 = [0., 1., 2., 3., 0.]
#     ```
#     </details>

# + nbgrader={"grade": false, "grade_id": "solve_poisson", "locked": false, "schema_version": 3, "solution": true, "task": false}
function solve_poisson(N, α)
    ### BEGIN SOLUTION
    Δx = 1/N
    x = (0:N)*Δx
    u = [0.; solve_nonlin(N, α); 0.]
    return x, u
    ### END SOLUTION
end

# + nbgrader={"grade": true, "grade_id": "cell-0daaf7cc0874d73d", "locked": true, "points": 1, "schema_version": 3, "solution": false, "task": false}
@assert solve_poisson(20, 1.) isa Tuple
@assert length(solve_poisson(20, 1.)) == 2
@assert solve_poisson(20, 1.)[1] == (0:20)/20
# -

# 3. Fixer $α = 30.0$ et illustrer dans ce cas la convergence de la méthode numérique quand $N \to \infty$,
#    en traçant sur un même graphe les solutions numériques obtenues pour les valeurs de $N$ données ci-dessous.

# + nbgrader={"grade": true, "grade_id": "plot_poisson_1", "locked": false, "points": 1, "schema_version": 3, "solution": true, "task": false}
α = 30.
p = plot(title=L"Convergence as $N \to \infty$", xlims=(0, 1),
         legend=:outertopright)
for N in (5, 10, 20, 30, 40)
    # Plot numerical solution for current value of `N`
    ### BEGIN SOLUTION
    x, u = solve_poisson(N, α)
    plot!(x, u, label="N = $N")
    ### END SOLUTION
end
p
# -

# 4. Fixer $N = 50$, et illustrer sur un même graphe
#    - La solution exacte de l'équation de Poisson quand $\alpha = 0$,
#      qui est donnée par
#      $$
#      u(x) = \frac{1}{8} - \frac{1}{2} \left(x - \frac{1}{2} \right)^2.
#      $$
#    - Les solutions numériques obtenues pour $α \in \{0, 10, 20, 50, 100\}$.

# + nbgrader={"grade": true, "grade_id": "plot_poisson_2", "locked": false, "points": 1, "schema_version": 3, "solution": true, "task": false}
N = 50
u_exact(x) = 1/8 - (x - .5)^2/2
p = plot(title="Solution to the nonlinear Poisson equation", xlims=(0, 1),
         legend=:outertopright)

### BEGIN SOLUTION
plot!(u_exact, label="Exact")
for α in (0., 10., 20., 50., 100.)
    x, u = solve_poisson(N, α)
    plot!(x, u, label="α = $α")
end
### END SOLUTION
p
# -

# 5. On suppose maintenant que le paramètre $\alpha$ est inconnu,
#    mais que la température moyenne de la barre a été mesurée:
#    $$
#    \tag{Constraint}
#    \int_{0}^{1} u(x) \, \mathrm dx = \overline u.
#    $$
#    On s'intéresse au problème inverse visant à calculer le paramètre $α$ correspondant à la température moyenne mesurée.
#    Pour ce faire, commencer par écrire une fonction `mean_temperature(N, α)`,
#    donnant la température moyenne dans la barre pour la valeur de $α$ donnée,
#    sur base d'une approximation numérique de la solution de l'équation de Poisson par `solve_poisson`.
#    On pourra supposer que le nombre $N$ est pair,
#    de manière à ce que la méthode de Simpson puisse étre utilisée.

# + nbgrader={"grade": false, "grade_id": "mean_temperature", "locked": false, "schema_version": 3, "solution": true, "task": false}
function mean_temperature(N, α)
    @assert N % 2 == 0
    ### BEGIN SOLUTION
    x, u = solve_poisson(N, α)
    result = u[1] + u[end]
    result += 4sum(u[2:2:end-1])
    result += 2sum(u[3:2:end-2])
    return result/N/3
    ### END SOLUTION
end

# + nbgrader={"grade": true, "grade_id": "2", "locked": true, "points": 0, "schema_version": 3, "solution": false, "task": false}
@assert mean_temperature(10, 10.) isa Float64
@assert mean_temperature(100, 0.) ≈ 1/12
@assert mean_temperature(100, 10.) < 1/12
# -

# 6. Écrire ensuite une fonction `newton_raphson(x, f; maxiter=100, ε = 1e-12)`
#    implémentant la méthode de Newton-Raphson pour résoudre l'équation scalaire $f(x) = 0$,
#    avec comme critère d'arrêt que $| f(x) | \leq \varepsilon$.
#    La fonction devra renvoyer uniquement le résultat final de l'itération,
#    ou `nothing` si une solution n'a pas été trouvée après `maxiter` itérations.
#    On considérera uniquement le cas où `x` est un scalaire et `f` est une fonction de $\mathbb R$ dans $\mathbb R$.
#
#    <details>
#         <summary>
#             <em><font color='gray'>Indication (cliquer pour afficher)</font></em>
#         </summary>
#
#     Pour calculer numériquement la dérivée de $f$,
#     la bibliothèque `ForwardDiff` peut être utilisée.
#     </details>

# + nbgrader={"grade": false, "grade_id": "newton_raphson", "locked": false, "schema_version": 3, "solution": true, "task": false}
const dx = ForwardDiff.Dual(0., 1.)
function newton_raphson(f, x; maxiter=100, ε = 1e-12)
    ### BEGIN SOLUTION
    for i in 1:maxiter
        y = f(x + dx)
        x -= y.partials[1]\y.value
        norm(f(x)) < ε && return x
    end
    error("Failed to converge!")
    ### END SOLUTION
end;

# + nbgrader={"grade": true, "grade_id": "cell-d2e7364ddc277254", "locked": true, "points": 2, "schema_version": 3, "solution": false, "task": false}
@assert newton_raphson(x -> x^2 - 2, 1) ≈ √2
@assert newton_raphson(x -> x^2 - 2, -1) ≈ -√2
@assert newton_raphson(x -> x^3 - 2, 1) ≈ cbrt(2)
@assert newton_raphson(x -> cos(x) - .5, 2) ≈ acos(.5)
# -

# 7. Écrire une fonction `get_alpha(N)` permettant de calculer $α$ de manière à ce que l'équation <a>(Constraint)</a> soit satisfaite avec $\overline u = .07$.
#    Pour ce faire, appliquer la méthode de Newton-Raphson à une fonction appropriée,
#    faisant intervenir un appel à `mean_temperature` avec la valeur de $N$ passée en argument.
#    <details>
#         <summary>
#             <em><font color='gray'>Indications (cliquer pour afficher)</font></em>
#         </summary>
#
#     Il vaut mieux ne pas initialiser l'itération de Newton-Raphson à $α = 0$.
#     En effet, la dérivée partielle par rapport à $α$ de la conductivité thermique $\kappa$ est égale à 0 en $α = 0$.
#     On pourra par exemple initialiser l'itération à $α = 10$
#     </details>

# + nbgrader={"grade": true, "grade_id": "get_alpha", "locked": false, "points": 1, "schema_version": 3, "solution": true, "task": false}
function get_alpha(N)
    ### BEGIN SOLUTION
    u_bar = .07
    α = newton_raphson(a -> mean_temperature(N, a) - u_bar, 10.)
    ### END SOLUTION
    return α
end;
# -

# 8. Illustrer sur un graphe la variation de la valeur trouvée pour $α$ en fonction de $N$,
#    par exemple pour $N \in \{10, \dotsc, 100\}$.

# + nbgrader={"grade": true, "grade_id": "plot_alpha", "locked": false, "points": 1, "schema_version": 3, "solution": true, "task": false}
### BEGIN SOLUTION
plot(title="Solution of inverse problem", xlabel=L"N", ylabel=L"α")
Ns = 10:2:100
plot!(Ns, get_alpha.(Ns))
### END SOLUTION
# -

# ### <font color='orange'>[Exercice 3]</font> Une équation différentielle à la rescousse

# Les plongeurs amateurs disposent parfois d'une bouée de secours autogonflable,
# qu'ils peuvent utiliser pour remonter rapidement à la surface en cas de problème.
# Cet exercice vise à étudier la remontée d'un plongeur à l'aide d'une telle bouée à laquelle il est attaché,
# et à délimiter les conditions d'utilisation de la bouée afin de garantir une ascension suffisamment rapide.
# On fera plusieurs hypothèses simplificatrices :
#
# - La bouée contient du CO₂, modélisé par un gaz idéal.
#
# - La température du CO₂ dans la bouée est constante, égale à 25°C
#
# - Bien que son volume soit variable en fonction de la pression,
#   la bouée reste toujours de forme sphérique.
#
# - Le volume du plongeur est indépendant de sa profondeur.
#
# - Le coefficient du traînée $C_d$ du système "plongeur + bouée" est constant, égal à 0.47.
#
# - On néglige la masse à vide de la bouée; toute sa masse provient du CO₂ qu'elle contient.

# Les paramètres physiques du problème sont donnés ci-dessous:

# +
# Temperature [K]
const T = 298.15  #  = 25°C

# Density of water [kg/m³]
const ρ_w = 1000

# Gravitational acceleration [m/s²]
const g = 9.81

# Mass of the diver [kg]
const m_diver = 80

# Volume of the diver [m³]
const V_diver = .07

# Cross section of diver [m²]
const A_diver = .25

# Drag coefficient
const C_d = 0.47

# Universal gas constant [J/mol·K]
const R = 8.314

# Molar mass of CO₂ [kg/mol]
const M_CO₂ = 0.04401

# Sea level pressure [Pa]
const P₀ = 101325
# -

# **Remarques**:
# - La densité du plongeur est supérieure à la densité de l'eau. Il s'ensuit que, sans bouée, un plongeur inconscient coulerait.
#
# - La masse de la bouée n'a pas encore été définie,
#   car elle a vocation à varier dans cet exercice.
#   On la notera `m_buoy`.

# Notons $z(t)$ la coordonnée verticale du plongeur.
# On supposera, par convention,
# que $z = 0$ correspond au niveau de la surface de l'eau,
# et que $z ≤ 0$ correspond au dessous de la surface.
# Notons aussi $v(t) = \dot z(t)$,
# la vitesse verticale du plongeur vers le haut.
#
# 1. Écrire une fonction `V_buoy(z, m_buoy)`,
#    retournant le volume de la bouée à une profondeur $z ≤ 0$,
#    lorsque celle-ci contient une masse `m_buoy` de CO₂.
#    Pour ce faire, utiliser la loi des gaz parfaits :
#    $$
#    V = \frac{n R T}{P},
#    \qquad n := \frac{m_{\rm buoy}}{M_{\rm CO_2}}.
#    $$
#    Ici $n$ est le nombre de moles de CO₂ contenues dans la bouée.
#    On rappelle que la pression hydrostatique à une profondeur $z ≤ 0$ est donnée par
#    $$
#    P(z) = P_0 - g \rho_w z.
#    $$

# + nbgrader={"grade": false, "grade_id": "V_buoy", "locked": false, "schema_version": 3, "solution": true, "task": false}
# Volume of the buoy at depth `z` (z ≤ 0)
function V_buoy(z, m_buoy)
    ### BEGIN SOLUTION
    P = P₀ - ρ_w * g * z
    return (m_buoy * R * T) / (P * M_CO₂)
    ### END SOLUTION
end;

# + nbgrader={"grade": true, "grade_id": "cell-0eef22c644cfc99f", "locked": true, "points": 1, "schema_version": 3, "solution": false, "task": false}
@assert V_buoy(0, 1.) ≈ (R * T) / (P₀ * M_CO₂)
@assert V_buoy(0, 1.)*P₀ ≈ V_buoy(-1, 1.)*(P₀ + ρ_w * g)
# -

# On modélise les forces sur le système,
# toutes dans la direction $z$, comme suit :
#
# - La **force de gravité** est donnée par
#
#   $$
#   F_g = - m \cdot g, \qquad m := m_{\rm diver} + m_{\rm buoy}.
#   $$
#
# - La **poussée d'Archimède** est donnée par
#
#   $$
#   F_a(z) = \rho_w \cdot V(z) \cdot g, \qquad V(z) := V_{\rm diver} + V_{\rm buoy}(z).
#   $$
#
# - La **force de traînée** est donnée par
#   $$
#   F_d(v) = - \frac{1}{2} \rho_w \cdot A_{\rm total} \cdot C_d \cdot v^2 \cdot {\rm sign}(v).
#   $$
#   Pour la surface de référence $A_{\rm total}$,
#   on prendra pour simplifier le maximum entre celle de la bouée et celle du plongeur.
#   Comme la bouée est supposée sphérique,
#   sa surface de référence peut être calculée explicitement,
#   et on a donc
#   $$
#   A_{\rm total} = \max\left\{ A_{\rm diver}, π\left(\frac{3V_{\rm buoy}}{4π}\right)^{\frac{2}{3}} \right\}.
#   $$
#
# En appliquant la seconde loi de Newton, on obtient
# une équation différentielle décrivant le mouvement du plongeur :
# $$
# \tag{EDO}
# \left\{
# \begin{aligned}
# \dot z &= v, \\
# \dot v &= \frac{F_g + F_a(z) + F_d(v)}{m_{\rm diver} + m_{\rm buoy}}.
# \end{aligned}
# \right.
# $$
#
# 2. Soit $X = (z, v)^\top$.
#    L'équation <a href="#EDO">(EDO)</a> peut-être réécrite sous la forme
#    $$
#    \tag{EDOv}
#    \dot X(t) = f\bigl(X; m_{\rm buoy})
#    $$
#    Écrire la fonction $f$ sous forme d'une fonction Julia `f(X, m_buoy)`,
#    où `m_buoy` est vu comme un paramètre.

# + nbgrader={"grade": false, "grade_id": "f_ode", "locked": false, "schema_version": 3, "solution": true, "task": false}
# Parameter `A` in the drag force (depends on volume `Vb` of the buoy)
A_total(Vb) = max(A_diver, π*(3Vb/4π)^(2/3))

function f(X, m_buoy)
    ### BEGIN SOLUTION
    z, v = X
    Vb = V_buoy(z, m_buoy)
    m = m_diver + m_buoy
    F_a = ρ_w * (V_diver + Vb) * g
    F_g = - m * g
    F_d = - 0.5 * C_d * ρ_w * A_total(Vb) * v^2 * sign(v)
    return [v, (F_g + F_a + F_d) / m]
    ### END SOLUTION
end

# + nbgrader={"grade": true, "grade_id": "cell-e119b00105e73c51", "locked": true, "points": 2, "schema_version": 3, "solution": false, "task": false}
@assert f([0., 0.], 0.) |> length == 2
@assert f([0., 0.], 0.)[1] == 0.
@assert f([0., 5.], 0.)[1] == 5.
@assert f([-1., 0.], 0.)[2] ≈ f([0., 0.], 0.)[2]
@assert f([-1., 0.], .1)[2] ≥ f([0., 0.], 0.)[2]
@assert f([0., 0.], 0.)[2] ≈ -1.22625
@assert f([0., 0.], .1)[2] ≈ 5.5709364718165455
# -

# 3. Écrire une fonction `rkx(Xₙ, h, Δ)` implémentant un pas de temps de taille $\Delta$ de la méthode de Runge-Kutta suivante pour une équation différentielle générique de la forme $X' = h(X)$:
#    $$
#       X_{n+1} = X_n + \frac{\Delta}{9}\left(2k_1 + 3k_2 + 4k_3 \right),
#    $$
#    où
#    \begin{align*}
#    k_1 &= \ h(X_n), \\
#    k_2 &= \ h\!\left(X_n + \frac{\Delta}{2} k_1\right), \\
#    k_3 &= \ h\!\left(X_n + \frac{3\Delta}{4} k_2\right).
#    \end{align*}
#    La fonction devra renvoyer $X_{n+1}$.

# + nbgrader={"grade": false, "grade_id": "rkx", "locked": false, "schema_version": 3, "solution": true, "task": false}
function rkx(Xₙ, h, Δ)
    k₁ = h(Xₙ)
    k₂ = h(Xₙ + Δ/2 * k₁)
    k₃ = h(Xₙ + 3Δ/4 * k₂)
    return Xₙ + Δ/9 * (2k₁ + 3k₂ + 4k₃)
end

# + nbgrader={"grade": true, "grade_id": "cell-35dffc49bd8d021b", "locked": true, "points": 1, "schema_version": 3, "solution": false, "task": false}
@assert rkx([0.], X -> [1.], 1.) ≈ [1]
@assert rkx([1.], X -> X, 1.)  ≈ [2 + 1/2 + 1/6]
# -

# 4. Écrire une fonction `solve_ode(Δ, z₀, m_buoy; tmax=20)` pour résoudre <a href="#EDOv">(EDOv)</a>
#    avec une condition initiale $X(0) = (z₀, 0)^\top$,
#    en utilisant la méthode de Runge-Kutta de point précédent avec pas de temps fixe `Δ`.
#    On supposera que `z₀ < 0`.
#    Votre fonction devra renvoyer un vecteur de temps `ts` et un vecteur de vecteurs `Xs` contenant la solution à ces temps.
#
#    On calculera le mouvement du plongeur jusqu'à ce qu'il ait atteint la surface,
#    ou jusqu'à `tmax` s'il n'a pas atteint la surface après ce temps.
#    Il faudra donc que soit `Xs[end][1] ≥ 0`, soit `ts[end] ≥ tmax`.

# + nbgrader={"grade": false, "grade_id": "solve_ode", "locked": false, "schema_version": 3, "solution": true, "task": false}
function solve_ode(Δ, z₀, m_buoy; tmax=20)
    X₀ = [z₀; 0.]
    ts = [0.]
    Xs = [X₀]
    ### BEGIN SOLUTION
    h(X) = f(X, m_buoy)
    while Xs[end][1] ≤ 0 && ts[end] < tmax
        push!(Xs, rkx(Xs[end], h, Δ))
        push!(ts, ts[end] + Δ)
    end
    ### END SOLUTION
    return ts, Xs
end

# + nbgrader={"grade": true, "grade_id": "cell-362b4d008e1484af", "locked": true, "points": 2, "schema_version": 3, "solution": false, "task": false}
@assert solve_ode(.01, -1., 0.) |> length == 2
@assert solve_ode(.01, -1., 0.)[1][end] ≈ 20
@assert solve_ode(.01, -1., .1)[2][1] |> length == 2
@assert solve_ode(.01, -1., .1)[2][end][1] ≥ 0
@assert solve_ode(.01, -1., .1)[2][end-1][1] ≤ 0
@assert solve_ode(.01, -10., .1)[1][end] > 5
@assert solve_ode(.01, -10., .1)[1][end] < 6
# -

# 5. Écrire une fonction `plot_z(Δ, z₀, ms)` permettant d'illustrer sur un même graphe la coordonnée $z$ du plongeur en fonction du temps,
# pour **une** valeur de $z_0$ donnée et **plusieurs** valeurs de `m_buoy` dans le vecteur `ms`.

# + nbgrader={"grade": true, "grade_id": "plot_z", "locked": false, "points": 1, "schema_version": 3, "solution": true, "task": false}
function plot_z(Δ, z₀, ms)
    p = plot(title="Depth of the diver")
    ### BEGIN SOLUTION
    for m ∈ ms
        ts, Xs = solve_ode(Δ, z₀, m)
        plot!(ts, [x[1] for x in Xs], label="m = $m")
        xlabel!("t [s]")
        ylabel!("z [m]")
    end
    ### END SOLUTION
    return p
end

Δ, z₀, ms = .01, -20., [.1, .2, .3]
plot_z(Δ, z₀, ms)
# -

# 6. Écrire une fonction `plot_v(Δ, z₀, ms)` permettant d'illustrer sur un même graphe la vitesse du plongeur en fonction du temps,
# pour **une** valeur de $z₀$ donnée et **plusieurs** valeurs de `m_buoy` dans le vecteur `ms`.

# + nbgrader={"grade": true, "grade_id": "plot_v", "locked": false, "points": 1, "schema_version": 3, "solution": true, "task": false}
function plot_v(Δ, z₀, ms)
    p = plot(title="Velocity of the diver")
    ### BEGIN SOLUTION
    for m ∈ ms
        ts, Xs = solve_ode(Δ, z₀, m)
        plot!(ts, [x[2] for x in Xs], label="m = $m")
        xlabel!("t [s]")
        ylabel!("v [m/s]")
    end
    ### END SOLUTION
    return p
end

Δ, z₀, ms = .01, -20., [.1, .2, .3]
plot_v(Δ, z₀, ms)
# -

# 7. On fixe à partir de maintenant la masse de la bouée à `m_buoy = .1`.
#    Écrire une fonction `rescue_time(z₀)`,
#    qui retourne une approximation du temps mis par le plongeur pour rejoindre la surface
#    à partir d'une profondeur initiale $z_0$.
#    Pour ce faire, résoudre l'équation différentielle avec un pas de temps fixe $Δ = 0.01$,
#    et estimer le temps de sauvetage par interpolation linéaire sur le dernier pas de temps,
#    durant lequel la coordonnée $z$ passe au dessus de 0.

# + nbgrader={"grade": true, "grade_id": "rescue_time", "locked": false, "points": 1, "schema_version": 3, "solution": true, "task": false}
function rescue_time(z₀)
    Δ, m_buoy = .01, .1
    ### BEGIN SOLUTION
    Ts, Xs = solve_ode(Δ, z₀, m_buoy)
    t₁, t₂, δx₁, δx₂ = Ts[end-1], Ts[end], Xs[end-1][1]-1, Xs[end][1]-1
    return (t₁*δx₂ - t₂*δx₁)/(δx₂-δx₁)
    ### END SOLUTION
end
# -

# 8. Faire un plot du temps de remontée en fonction de `z₀`,
#    pour des valeurs de ce paramètre dans l'intervalle $[-30, -5]$.
#    Estimer graphiquement, par exemple à l'aide de la fonction `vline`,
#    la profondeur maximale telle que le plongeur peut rejoindre la surface en moins de 10 secondes grâce à la bouée.

# + nbgrader={"grade": true, "grade_id": "plot_rescue_time", "locked": false, "points": 1, "schema_version": 3, "solution": true, "task": false}
### BEGIN SOLUTION
plot(-30:.1:-5, rescue_time)
vline!([-16.2])
hline!([10])
xlabel!("z₀ [m]")
ylabel!("Rescue time [s]")
### END SOLUTION
# -

# 9. **Question ouverte**. Afin de déterminer les normes d'utilisation de la bouée,
#    calculer avec précision,
#    par la méthode de votre choix,
#    la valeur de `z₀` permettant au plongeur de rejoindre la surface en exactement 10 secondes.
#    Vous pouvez pour ce faire utiliser la bibliothèque `ForwardDiff`.

# + nbgrader={"grade": true, "grade_id": "bonus", "locked": false, "points": 1, "schema_version": 3, "solution": true, "task": false}
# BEGIN SOLUTION
z = -10
for i in 1:10
    r = rescue_time(z + dx)
    z -= (r.value - 10) / r.partials[1]
end
z
### END SOLUTION

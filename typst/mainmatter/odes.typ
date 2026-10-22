#import "../template.typ": *

// The template's `dd` (math.upright("d")) is followed by a space, e.g. "d s"
#let dd = $upright(d)$

// \left\{ \begin{aligned} a &= b \\ c &= d \end{aligned} \right.: rows separated
// by `;`, alignment points with `&`, first column right-aligned (unlike cases).
// Commas that are not inside parentheses must be escaped: `\,`.
#let system = math.mat.with(delim: ("{", none), align: right)

// \begin{proof}[Title]
#let proof-titled(title, body) = block(width: 100%, above: leading + 6pt, below: leading + 6pt, breakable: true,
  [#emph[#title.]#h(0.5em)#body#h(1fr)$square.stroked$])

= Numerical ordinary differential equations <cha:ode>
#minitoc()

#section-star[Introduction]
This chapter concerns the numerical solution of ordinary differential equations (ODEs) of the following form:
#numbered($
  cases(
    vect(x)'(t) = vect(f) big((t, vect(x)(t)))\,,
    vect(x)(t_0) = vect(x)_0.
  )
$) <eq:ode>
Here $vect(f) colon real times real^n -> real$ and $vect(x)_0$ is the initial condition.
Equations of this type
are the building blocks of a plethora mathematical models in science and engineering.
They have applications in celestial dynamics, molecular simulation and fluid mechanics, to mention just a few.
Ordinary differential equations also arise after discretization of time-dependent partial differential equations,
which are also ubiquitous in science.
More often than not,
it is not possible to find an explicit solution of~@eq:ode,
and so one has to resort to numerical simulation.
The rest of chapter is organized as follows:
- In @sec:existence,
  we define the concepts of local and global solutions for the continuous-time problem~@eq:ode,
  and we recall fundamental results concerning the existence and uniqueness of a solution.

- In~@sec:one-step_methods,
  we analyze the so-called _one-step_ numerical methods to solve~@eq:ode.
  We emphasize in particular the concepts of _consistency_, _stability_ and _convergence_.

- In~@sec:multistep_methods,
  we present _multistep_ methods to solve~@eq:ode,
  and discuss their drawbacks and advantages compared to one-step methods.

- Finally, in~@sec:absolute_stability,
  we introduce the concept of _absolute stability_ and discuss its relevance in the context of _stiff_ differential equations.

== Analysis of the continuous problem <sec:existence>
A differentiable function $vect(x) colon I -> real^n$,
where $I$ denotes an interval of $real$ containing $t_0$,
is a solution of~@eq:ode if $vect(x)(t_0) = vect(x)_0$ and the equation~@eq:ode is satisfied for all~$t in I$.
The solution is called global if~$I = real$,
and local otherwise.

#paragraph[Integral formulation.][
  If $vect(x)$ is a solution to~@eq:ode,
  then it holds that
  #numbered($
    forall t in I, wide
    vect(x)(t) = vect(x)_0 + integral_(t_0)^t vect(f) big((s, vect(x)(s))) thin dd s.
  $) <eq:ode_integral>
  The converse statement is not true in general,
  because a solution to~@eq:ode_integral need not necessarily be differentiable everywhere.
  However, if the integral formulation~@eq:ode_integral holds,
  then necessarily $vect(x)$ is absolutely continuous
  and~@eq:ode is satisfied for almost every~$t$.
  Additionally, if~@eq:ode_integral is satisfied and the function~$vect(f)$ is continuous,
  then the function~$s |-> vect(f) big((s, vect(x)(s)))$ is continuous,
  and so~@eq:ode is satisfied for all $t in I$ by the fundamental theorem of analysis.
  We now focus on the integral formulation~@eq:ode_integral,
  and begin by establishing existence of a local solution.
]

#theorem(title: [Existence of a solution])[
  Let $vect(x)_0 in real^n$ and let $Omega_(cal(T), cal(R))$ denote the set
  $
    big(\{)(t, vect(x)) in real times real^n: abs(t-t_0) <= cal(T) " and " vecnorm(vect(x) - vect(x)_0) <= cal(R) big(\}),
  $
  Assume that the following conditions are satisfied for some $cal(T) > 0$ and $cal(R) > 0$:
  - The function $vect(f)$ is uniformly bounded on~$Omega_(cal(T), cal(R))$:
    #numbered($
      forall (t, vect(x)) in Omega_(cal(T), cal(R)), wide
      vecnorm(vect(f)(t, vect(x))) <= M.
    $) <eq:boundedness_ode>

  - The function $vect(f)$ satisfies the following Lipschitz condition:
    there is $L > 0$ such that
    #numbered($
      forall big(((t, vect(x)_1), (t, vect(x)_2))) in Omega_(cal(T), cal(R)) times Omega_(cal(T), cal(R)), wide
      vecnorm(vect(f)(t, vect(x)_1) - vect(f)(t, vect(x)_2)) <= L vecnorm(vect(x)_1 - vect(x)_2).
    $) <eq:local_lipschitz_ode>
  Then there exists $T in (0, cal(T)]$ depending on $cal(R)$, $M$ and $L$ such that
  the differential equation~@eq:ode_integral has a local solution $vect(x) colon [t_0 - T, t_0 + T] -> real^n$.
] <theorem:existence>
#proof[
  Fix $T in (0, cal(T)]$ and let $I = [t_0 - T, t_0 + T]$.
  Let also $cal(X)$ denote the following subset of continuous functions defined from~$I$ to $real^n$:
  $
    cal(X) := lr({ vect(x) in C lr((I, real^n)) : sup_(t in I) big(||vect(x)(t) - vect(x)_0||) <= cal(R) })
  $
  The set~$cal(X)$ endowed with supremum metric is a closed subset of~$C(I, real^n)$.
  Since $cal(X)$ is a closed subset of a complete metric space,
  it is itself complete.
  Let $Phi colon cal(X) -> C(I, real^n)$ denote the mapping
  $
    Phi(vect(x)) colon t |-> vect(x)_0 + integral_0^t vect(f) big((s, vect(x)(s))) thin dd s.
  $
  The right-hand side, being the integral of a bounded function,
  is indeed a continuous function.
  We will show that, for~$T$ sufficiently small,
  - the mapping $Phi$ maps $cal(X)$ into $cal(X)$;
  - the mapping $Phi$ is a contraction.
  From~@eq:boundedness_ode,
  it follows that
  $
    forall vect(x) in cal(X), wide
    forall t in I, wide
    lr(||Phi(vect(x))(t) - vect(x)_0||) = lr(||integral_(t_0)^t vect(f) big((s, vect(x)(s))) thin dd s||)
    <= M T.
  $
  On the other hand, from the Lipschitz condition~@eq:local_lipschitz_ode,
  it holds that
  $
    forall (vect(x), vect(y)) in cal(X) times cal(X), wide
    norm(Phi(vect(x)) - Phi(vect(y))) <= L T norm(Phi(vect(x)) - Phi(vect(y))).
  $
  Therefore, it suffices to take $T < min lr({ cal(T), cal(R)/M, 1/L })$ to ensure that the above conditions are satisfied.
  For a value of~$T$ in this range,
  the Banach fixed point theorem, @theorem:banach_fixed_point,
  gives the existence of a unique fixed point~$vect(x)_* in cal(X)$ of~$Phi$.
  Since a fixed point of~$Phi$ is a solution to~@eq:ode_integral in view of the definition of~$Phi$,
  the statement is proved.
]
It may seem at first glance that uniqueness of the solution to~@eq:ode_integral follows from the uniqueness of the fixed point guaranteed by~@theorem:banach_fixed_point.
However, this theorem implies uniqueness _only in the set~$cal(X)$_,
a property known as~_conditional uniqueness_.
In order to prove that the solution is unique over the full space $C([t_0-T, t_0+T], real^n)$,
additional assumptions and arguments are required.
A simple approach is to rely on Grönwall's lemma.
#lemma(title: [Grönwall's lemma, simplified integral form])[
  Suppose that $u colon [t_0-T, t_0+T] -> real_(>= 0)$ is continuous, nonnegative, and satisfies
  #numbered($
    forall t in [t_0, t_0+T], wide
    u(t) <= alpha + integral_(t_0)^t beta(s) thin u(s) thin dd s,
  $) <eq:gronwall_condition>
  where $alpha >= 0$ and $beta colon [t_0, t_0 + T] -> real_(>= 0)$ is continuous and nonnegative.
  Then
  #numbered($
    forall t in [t_0, t_0+T], wide
    u(t) <= alpha exp lr((integral_(t_0)^t beta(s) thin dd s)).
  $) <eq:gronwall>
]
#proof[
  Assume first that~$alpha > 0$,
  so that the logarithm in~@eq:intermediate_gronwall is well-defined.
  By the fundamental theorem of calculus and~@eq:gronwall_condition,
  it holds that
  $
    dd/(dd t) lr((alpha + integral_(t_0)^t beta(s) thin u(s) thin dd s))
    <= beta(t) lr((alpha + integral_(t_0)^t beta(s) thin u(s) thin dd s))
  $
  Therefore we have
  #numbered($
    dd/(dd t) log lr((alpha + integral_(t_0)^t beta(s) thin u(s) thin dd s)) <= beta(t),
  $) <eq:intermediate_gronwall>
  and after integrating and exponentiating,
  we obtain
  $
    alpha + integral_(t_0)^t beta(s) thin u(s) thin dd s
    <= alpha exp lr((integral_(t_0)^t beta(s) thin dd s))
  $
  The statement then follows by using~@eq:gronwall_condition again.
  Assume next that $alpha = 0$.
  If~@eq:gronwall_condition is satisfied for~$alpha = 0$,
  then this condition is also satisfied for all $alpha > 0$.
  Therefore the conclusion~@eq:gronwall holds for all~$alpha > 0$,
  and taking the limit $alpha -> 0$ in this equation,
  we obtain the statement.
]
Note that the estimate~@eq:gronwall is sharp,
since the function $v colon [t_0, t_0 + T]$ given by
$
  v(t) = alpha exp lr((integral_(t_0)^t beta(s) thin dd s))
$
satisfies~@eq:gronwall_condition with equality.
We are now ready to prove uniqueness under an appropriate condition.
#theorem(title: [Uniqueness of the solution])[
  Let $vect(x)_0 in real^n$ and let
  $
    Omega_(cal(T), cal(R))
    big(\{)(t, vect(x)) in real times real^n: abs(t - t_0) <= cal(T) " and " vecnorm(vect(x) - vect(x)_0) <= cal(R) big(\}),
  $
  Assume that for all $cal(T) in real_(>0)$ and $cal(R) in real_(>0)$,
  there is $L_(cal(T), cal(R))$ such that
  #numbered($
    forall big(((t, vect(x)_1), (t, vect(x)_2))) in Omega_(cal(T), cal(R)) times Omega_(cal(T), cal(R)), wide
    vecnorm(vect(f)(t, vect(x)_1) - vect(f)(t, vect(x)_2)) <= L_(cal(T), cal(R)) vecnorm(vect(x)_1 - vect(x)_2).
  $) <eq:local_lipschitz_uniqueness>
  Then if $vect(x)_1$ and $vect(x)_2$ in $C big(([t_0 - T, t_0 + T], real^n))$ are local solutions to~@eq:ode_integral,
  it holds that $vect(x)_1 = vect(x)_2$.
] <theorem:uniqueness>
#proof[
  Suppose that $vect(x)_1$ and $vect(x)_2$ are solutions to~@eq:ode_integral.
  Let~$I = [t_0 - T, t_0 + T]$ and
  $
    R := max lr({ sup_(t in I) norm(vect(x)_1(t) - vect(x)_0), sup_(t in I) norm(vect(x)_2(t) - vect(x)_0) }) < infinity,
  $
  Since $vect(x)_1$ and $vect(x)_2$ are solutions,
  it holds that
  $
    forall t in [t_0 - T, t_0 + T], wide
    vect(x)_1(t) - vect(x)_2(t) = integral_(t_0)^t Big((f big((s, vect(x)_1(s))) - f big((s, vect(x)_2(s))))) thin dd s.
  $
  Taking the norm and using~@eq:local_lipschitz_uniqueness,
  we obtain
  $
    forall t in [t_0, t_0 + T], wide
    vecnorm(vect(x)_1(t) - vect(x)_2(t)) <= L_(T,R) integral_(t_0)^t vecnorm(vect(x)_1(s) - vect(x)_2(s)) thin dd s
  $
  Using Grönwall's lemma,
  we deduce that~$vect(x)_1(t) = vect(x)_2(t)$ for all $t in [t_0, t_0 + T]$.
  A similar argument can be employed to show that $vect(x)_1 = vect(x)_2$ on $[t_0 - T, t_0]$.
]

#corollary(title: [Maximal solutions])[
  Assume that $vect(f)$ is continuous in~$t$ and satisfies the local Lipschitz condition~@eq:local_lipschitz_uniqueness.
  Then there exists $- infinity <= T_- < T_+ <= infinity$ such that $t_0 in (T_-, T_+)$ and
  the following properties are satisfied.
  - there exists a solution $vect(x)_* colon (T_-, T_+) -> real^n$ to~@eq:ode_integral;

  - if $vect(x) colon I -> real^n$ is a local solution of~@eq:ode_integral,
    then $I subset (T_-, T_+)$ and $vect(x)(t) = vect(x)_* (t)$ for all~$t in I$.

  - If $T_+$ is finite, then $lim_(t -> T_+) big(||vect(x)(t)||) = infinity$,
    and if $T_-$ is finite,
    then $lim_(t -> T_-) big(||vect(x)(t)||) = infinity$.
  The solution $vect(x)_*$ is called the maximal solution of~@eq:ode_integral.
] <corollary:maximal_solutions>
#proof[
  Let $cal(I)$ denote the union of all the open intervals~$I$ such that there exists a solution in~$C(I, real^n)$ to~@eq:ode_integral.
  The open set $cal(I)$ is connected and, by~@theorem:existence,
  it contains a neighborhood of~$t_0$.
  Therefore $cal(I)$ is of the form $(T_-, T_+)$,
  where $-infinity <= T_- < t_0 < T_+ <= infinity$.
  In view of~@theorem:uniqueness,
  all the local solutions coincide where they are defined,
  and so they can be patched together in order to construct a solution~$vect(x)_* colon (T_-, T_+) -> real$.
  It remains to prove the third item.
  To this end, suppose for contradiction that $T_+$ was finite and that
  there was $(t_n)_(n in nat)$ such that $t_n -> T_+$ in the limit~$n -> infinity$ and~
  $
    K := sup_(n in nat) vecnorm(vect(x)_* (t_n)) < infinity.
  $
  Since $vect(f)$ is continuous,
  there is $M$ such that $|f(t, vect(x))|$ is uniformly bounded from above by~$M$ for all $(t, vect(x)) in [T_- - 1, T_+ + 1] times B_(K+1)(vect(0))$.
  Furthermore, by the assumption~@eq:local_lipschitz_uniqueness,
  there is $L$ such that for all $t in [T_- - 1, T_+ + 1]$,
  the following Lipschitz condition holds:
  $
    forall (vect(x)_1, vect(x)_2) in B_(K+1)(vect(0)) times B_(K+1)(vect(0)), wide
    vecnorm(vect(f)(t, vect(x)_1) - vect(f)(t, vect(x)_2)) <= L vecnorm(vect(x)_1 - vect(x)_2).
  $
  Consequently, @theorem:existence with $cal(T) = cal(R) = 1$ implies for all~$n$ the existence of a solution to
  $
    cases(
      vect(x)'(t) = vect(f) big((t, vect(x)(t)))\,,
      vect(x)(t_n) = vect(x)_* (t_n).
    )
  $
  over the time interval $[t_n - T, t_n + T]$, where $T > 0$ depends only on $M$ and $L$, and not on~$n$.
  But then, for~$n$ sufficiently large,
  this solution extends beyond $T_+$,
  which contradicts the maximality of~$cal(I)$.
  An analogous reasoning can be employed for~$T_-$.
]
#example[
  Consider the ODE
  $
    cases(
      x'(t) = x(t)^2\,,
      x(0) = 1.
    )
  $
  The maximal solution is $x_* colon (-infinity, 1) -> real$ given by
  $
    vect(x)_* (t) = 1/(1-t).
  $
]

#paragraph[Existence of a unique global solution.][
  In certain settings,
  it is possible to prove the maximal solution to~@eq:ode_integral is globally defined for any initial condition.
  We discuss a few important examples.
  - The first case is when $vect(f) colon real times real^n$ is globally Lipschitz in its second argument,
    with a Lipschitz constant that depends continuously on the first argument.

  - The second case,
    generalizing the first,
    is when the growth of~$vect(f)(t, placeholder)$ is at most affine:
    $
      forall (t, vect(x)) in real times real^n, wide
      vecnorm(vect(f)(t, vect(x))) <= C(t) + L(t) vecnorm(vect(x)),
    $
    with continuous constants~$C(t)$ and $L(t)$.

  - The third case is when $vect(f)$ is independent of~$t$ and there is a function~$W in C^1(real^n)$ such that~$W(vect(x)) -> infinity$ in the limit as $vecnorm(vect(x)) -> infinity$
    and
    $
      forall vect(x) in real^n, wide
      nabla W(vect(x)) dot vect(f)(vect(x)) <= c < infinity
    $
    Such a function is called a _Lyapunov function_.
  The strategy of proof for global existence usually relies on an argument by contradiction.
  Consider for example the third setting.
  Since the assumptions of~@corollary:maximal_solutions are satisfied,
  there exists a maximal solution $vect(x)_* colon (T_-, T_+) -> infinity$.
  Assume for contradiction that~$T_+$ is finite.
  Then the third item in~@corollary:maximal_solutions implies that $lim_(t -> T_+) vecnorm(vect(x)_* (t)) -> infinity$,
  and so~$W big((vect(x)_* (t)))$ blows up as~$t$ approaches $T^+$.
  On the other hand, we have
  $
    dd/(dd t) W big((vect(x)_* (t))) = nabla W big((vect(x)_* (t))) dot vect(f) big((vect(x)_* (t))) <= c.
  $
  Therefore $lim_(t -> T_+) W big((vect(x)_* (t))) <= W big((vect(x)_* (t_0))) + |c| (T_+ - t_0)$,
  which is a contradiction.
]

== One-step methods <sec:one-step_methods>
From now on, we assume for simplicity that $t_0 = 0$ and that
the initial value problem~@eq:ode admits a unique solution over the interval $[0, T]$.
Most numerical methods for ODEs construct an approximation of the solution at discrete points:
$
  vect(x)_n approx vect(x)(t_n), wide n = 0, 1, 2, dots.
$
The discretization points $(t_n)_(n in nat)$ are commonly equidistant,
i.e.~$t_n = n Delta$ where $Delta$ is the _discretization step_.
Sometimes, it is useful to employ a variable time step,
but we assume throughout this section that the time step is fixed, for simplicity.
We begin in~@sub:explicit_euler and~@sub:implicit_euler by studying the simplest one-step methods,
namely the forward and backward Euler methods.
Then, in~@sub:one_step_general,
we present a general approach to the analysis of one-step methods.
Finally, in~@sub:other_one_step,
we present other widely used one-step methods in applications.

=== Forward Euler method <sub:explicit_euler>
Assume that~@eq:ode has a unique solution~$vect(x)(t)$ over the interval $[0, T]$.
If $vect(x)(t)$ is twice continuously differentiable,
then by Taylor's formula,
we have
#numbered($
  vect(x)(t + Delta) = vect(x)(t) + Delta vect(f)(t, vect(x)) + Delta^2/2 vect(x)''(tau),
  wide tau in (t, t+Delta).
$) <eq:taylor_forward>
This motivates a method known as the _forward_ or _explicit_ Euler method:
$
  vect(x)_(n+1) = vect(x)_n + Delta vect(f)(t_n, vect(x)_n),
$
with the same initial condition as for the continuous equation~@eq:ode.
The convergence of this method can be proved under a global Lipschitz assumption on the function~$vect(f)$.

#theorem(title: [Convergence of the forward Euler method])[
  Assume that there is $L in real_(>0)$ such that
  #numbered($
    forall (t, vect(x), vect(y)) in real times real^n times real^n, wide
    vecnorm(vect(f)(t, vect(x)) - vect(f)(t, vect(y)))
    <= L vecnorm(vect(x) - vect(y)).
  $) <eq:global_lipschitz>
  Suppose in addition that there exists a unique, twice continuously differentiable of~@eq:ode over the interval~$[0, T]$,
  and let
  $
    M = sup_(t in [0,T]) vecnorm(vect(x)''(t))
  $
  Then the following error estimate holds:
  #numbered($
    forall n in lr({0, 1, dots, lr(floor.l T/Delta floor.r) }),
    wide
    vecnorm(vect(x)(t_n) - vect(x)_n)
    <=
    (Delta M)/2 lr(((ee^(L t_n) - 1)/L)).
  $) <eq:error_bound_forward_euler>
] <theorem:forward_euler>
#proof[
  By Taylor's theorem,
  it holds that
  $
    vect(x)(t_n) = vect(x)(t_(n-1)) + Delta vect(f) big((t_(n-1), vect(x)(t_(n-1)))) + Delta^2/2 vect(alpha)_n,
    wide
    vect(alpha)_n := 2 integral_0^1 (1 - s) vect(x)''(t_n + Delta s) thin dd s.
  $
  Notice that that $norm(vect(alpha)_n) <= M.$
  Therefore, it holds that
  $
    vect(x)(t_n) - vect(x)_n
    &= lr((vect(x)(t_(n-1)) + Delta thin vect(f) big((t_(n-1), vect(x)(t_(n-1)))) + Delta^2/2 vect(alpha)_n))
    - Big((vect(x)_(n-1) + Delta thin vect(f) big((t_(n-1), vect(x)_(n-1))))) \
    &= big((vect(x)(t_(n-1)) - vect(x)_(n-1))) + Delta Big((vect(f) big((t_(n-1), vect(x)(t_(n-1)))) - vect(f) big((t_(n-1), vect(x)_(n-1))))) + Delta^2/2 vect(alpha)_n,
  $
  Let $vect(e)_n = vect(x)(t_n) - vect(x)_n$ and $vect(epsilon)_n = Delta^2/2 vect(alpha)_n$.
  The first term is the error at iteration~$n-1$,
  and the second may be bounded from~@eq:global_lipschitz,
  which gives
  $
    vecnorm(vect(e)_n)
    &<= (1 + Delta L) vecnorm(vect(e)_(n-1))
    + vecnorm(vect(epsilon)_n).
  $
  The structure of this equation is important,
  as it appears in the analysis of all one-step methods for ODEs.
  The first term is an amplification of the error at the previous iteration,
  and the second term is an upper bound on the additional error introduced at step~$n$.
  Applying this inequality to the previous time steps, we obtain
  #eq-last($
    vecnorm(vect(e)_n)
    &<= (1 + Delta L) Big(((1 + Delta L) vecnorm(vect(e)_(n-2)) + vecnorm(vect(epsilon)_(n-1)))) + vecnorm(vect(epsilon)_n) \
    &<= dots
    <= (1 + Delta L)^n vecnorm(vect(e)_0) + sum_(i=1)^n (1 + Delta L)^(n-i) vecnorm(vect(epsilon)_i).
  $) <eq:sum_local_and_accumulation>
  Since $vecnorm(vect(epsilon)_i) <= Delta^2 M\/2$,
  we have by using the formula for geometric series that
  $
    vecnorm(vect(e)_n)
    <= (1 + Delta L)^n vecnorm(vect(e)_0) + ((1+Delta L)^n - 1)/(Delta L) lr(((Delta^2 M)/2)).
  $
  The first term is zero because $vecnorm(vect(e)_0) = 0$.
  Using the bound $(1+Delta L)^n <= big((exp(Delta L)))^n = ee^(L t_n)$ in the second term and rearranging,
  we finally obtain the statement~@eq:error_bound_forward_euler.
]

=== Backward Euler method <sub:implicit_euler>
If we apply the Taylor expansion~@eq:taylor_forward backward around~$t + Delta$,
instead of forward around~$t$,
then we obtain
$
  vect(x)(t) = vect(x)(t + Delta) - Delta vect(f)(t + Delta, vect(x)) + cal(O)(Delta^2).
$ <eq:taylor_backward>
This motivates the so-called _backward_ or _implicit_ Euler method:
#numbered($
  vect(x)_(n+1) = vect(x)_n + Delta vect(f)(t_(n+1), vect(x)_(n+1)).
$) <eq:backward_euler>
Observe that the right-hand side depends on $vect(x)_(n+1)$.
Therefore, given $t_n$ and $vect(x)_n$,
this is a nonlinear equation for the unknown $vect(x)_(n+1)$,
which can be solved by using any of the methods studied in~@cha:solution_of_nonlinear_systems.
Finding a solution to~@eq:backward_euler amounts to finding a fixed point of the function
$
  vect(y) |-> vect(F)(vect(y)) := vect(x)_n + Delta vect(f)(t_n, vect(y)).
$
A priori, the existence and uniqueness of such a fixed point is not guaranteed.
We proved in~@theorem:exponential_convergence_fixed_point that a sufficient condition for these two properties to hold is
that~$vect(F)$ is globally Lipschitz with a constant strictly less than 1,
which holds if and only if the function~$vect(y) |-> vect(f)(t_n, vect(y))$ is globally Lipschitz with a constant strictly less than $1\/Delta$.
If the condition~@eq:global_lipschitz holds, for example,
then the backward Euler method~@eq:backward_euler is guaranteed to be well defined for $Delta < 1/L$.
@theorem:exponential_convergence_fixed_point also ensures that, if~$vect(F)$ is globally Lipschitz with a constant less than~$1$,
then the fixed point can be approximated by using the iteration
#numbered($
  vect(y)_(k+1) = vect(F)(vect(y)_k).
$) <eq:fixed_point_ode>
and there is exponential convergence~$vect(y)_k -> vect(x)_(n+1)$ in the limit as $k -> infinity$.
A natural starting point for~@eq:fixed_point_ode is~$vect(y)_0 = vect(x)_n$.
An alternative approach to the fixed point iteration~@eq:fixed_point_ode is to use the Newton--Raphson method for~@eq:backward_euler,
which is faster in principle but must be initialized sufficiently close to the fixed point.

Using a reasoning similar to that employed for proving~@theorem:forward_euler,
we can prove the following result.
#theorem(title: [Convergence of the backward Euler method])[
  If the assumptions of~@theorem:forward_euler hold and $Delta < 1/L$,
  then the following error estimate holds:
  #numbered($
    forall n in lr({0, 1, dots, lr(floor.l T/Delta floor.r) }),
    wide
    vecnorm(vect(x)(t_n) - vect(x)_n)
    <=
    (Delta M)/2 lr(((lr((1/(1 - Delta L)))^n - 1)/L)).
  $) <eq:estimate_backward_euler>
]
#proof[
  The proof is left as an exercise.
]
#remark[
  Note that, if $Delta < 1/(2L)$,
  then
  $
    1/(1 - Delta L)
    &= 1 + Delta L + (Delta L)^2 + (Delta L)^3 + (Delta L)^4 dots \
    &<= 1 + Delta L + (Delta L)^2 + 1/2 (Delta L)^2 + 1/4 (Delta L)^2 + dots \
    &<= 1 + Delta L + 2 (Delta L)^2 <= exp big((Delta L + (Delta L)^2)),
  $
  and so the error estimate~@eq:estimate_backward_euler gives
  $
    vecnorm(vect(x)(t_n) - vect(x)_n)
    <=
    (Delta M)/2 lr(((exp(L t_n + Delta L^2 t_n) - 1)/L)),
  $
  which makes it clear that the right-hand side of~@eq:estimate_backward_euler is close,
  in absolute and relative terms, to that of~@eq:error_bound_forward_euler when $Delta << 1$.
]
At this point, the reader may be wondering why one would use the backward Euler method instead of the forward Euler method,
given that both methods have same order of convergence but iterations of the former are more computationally costly.
The reason is that the backward Euler method,
like many implicit methods,
is more _stable_ than its forward counterpart.
Implicit methods are especially attractive in the context of~_stiff_ differential equations.
We shall elaborate on this subject in~@sec:absolute_stability.

=== Analysis of general one-step methods <sub:one_step_general>
In general, one-step methods to solve differential equations are of the abstract form
#numbered($
  vect(x)_(n+1) = vect(x)_n + Delta vect(Phi)_Delta (t_n, vect(x)_n).
$) <eq:general_one_step>
where $vect(Phi)_Delta colon real times real^n -> real^n$ is a function such that
#numbered($
  vect(Phi)_Delta (t, vect(x))
  approx 1/Delta integral_t^(t+Delta) vect(f) big((s, vect(x)^(t, vect(x))(s))) thin dd s
  = (vect(x)^(t, vect(x))(t + Delta) - vect(x))/Delta.
$) <eq:one_step_integral>
Here $vect(x)^(t, vect(x))$ denotes the solution to the differential equation~@eq:ode with initial condition $vect(x)(s) = vect(x)$.
The main goal of this section is to establish general conditions,
known as _consistency_ and _stability_,
under which the numerical scheme~@eq:general_one_step is convergent.
As we observed in the proof of~@theorem:forward_euler
-- specifically in equation~@eq:sum_local_and_accumulation --
the error at the final iteration for the forward Euler method is a sum of local errors,
each amplified by a factor depending to the number of iterations left to reach the final time.
_Consistency_ of a numerical method enables to control the size of local errors when they arise,
while~_stability_ enables to control their growth.

We emphasize that both the forward and the backward Euler methods can be recast in the form~@eq:general_one_step.
For the forward Euler method $vect(Phi)_Delta (t, vect(x)) = vect(f)(t, vect(x))$,
while for the backward Euler method, the function $vect(Phi)_Delta$ is defined implicitly as the function which to~$(t, vect(x))$ associates the solution~$vect(phi.alt) in real^n$ to the equation
$
  vect(phi.alt) = vect(f)(t + Delta, vect(x) + Delta vect(phi.alt)).
$
==== Local truncation error and consistency
The local truncation error is the residual error obtained when substituting the exact solution of the differential equation in~@eq:general_one_step:
$
  vect(eta)_(n+1) := (vect(x)(t_(n+1)) - vect(x)(t_n))/Delta - vect(Phi)_Delta big((t_n, vect(x)(t_n))).
$
Since there is a division by~$Delta$, the local truncation error has the same physical dimension as that of $vect(x)'$,
and so it should be viewed as an error _per time unit_.

#definition(title: [Consistency])[
  A numerical method is consistent if
  $
    lim_(Delta -> 0) lr((max_(1 <= n <= N) vecnorm(vect(eta)_n))) = 0, wide N = lr(floor.l T/Delta floor.r).
  $
  It is consistent with order~$p$ if there exists~$C$ such that
  $
    forall Delta > 0, wide
    max_(1 <= n <= N) vecnorm(vect(eta)_n) <= C Delta^p.
  $
]
Proving the consistency of a numerical method is usually achieved on a case-by-case basis by application of Taylor's formula.

==== Stability <ssub:Stability>
The stability of a numerical method qualifies its sensitivity to perturbations.
Roughly speaking, it expresses that small perturbation of the right-hand side of~@eq:general_one_step lead to small perturbations of the numerical solution.
#definition(title: [Stability])[
  A numerical method of the form~@eq:general_one_step is stable if there exists a constant~$S(T) > 0$ independent of $Delta$ such that for all sequence~$(vect(y)_n)_(1 <= n <= N)$ satisfying
  #numbered($
    vect(y)_(n+1) = vect(y)_n + Delta vect(Phi)_Delta (t_n, vect(y)_n) + Delta vect(delta)_(n+1), wide vect(y)_0 = vect(x)_0,
  $) <eq:general_one_step_perturbed>
  it holds that
  #numbered($
    max_(1 <= n <= N) vecnorm(vect(x)_n - vect(y)_n) <= S(T) Delta sum_(n=1)^N vecnorm(vect(delta)_n).
  $) <eq:stability_numerical_method>
]
It is convenient to introduce the following norms for sequences of vectors $(vect(u)_n)_(1 <= n <= N)$:
$
  norm(vect(u)_placeholder)_(ell^1_T) = sum_(n=1)^N vecnorm(vect(u)_n),
  wide
  norm(vect(u)_placeholder)_(ell^infinity_T) = max_(1 <= n <= N) vecnorm(vect(u)_n),
$
With these notations, equation~@eq:stability_numerical_method may be rewritten compactly as follows:
$
  norm(vect(x)_placeholder - vect(y)_placeholder)_(ell^infinity_T) <= S(T) Delta norm(vect(delta)_placeholder)_(ell^1_T)
$
One could argue that this equation is neater than~@eq:stability_numerical_method;
it bounds a norm of just one mathematical object, namely the sequence $(vect(x)_n - vect(y)_n)_(1 <= n <= N)$,
by a norm of another object, namely the sequence $(vect(delta)_n)_(1 <= n <= N)$.
Arguments for proving that a numerical scheme is stable often rely on some form of Lipschitz continuity.
If the function $vect(Phi)_Delta (t, vect(y))$ is globally Lipschitz continuous with respect to~$vect(y)$,
then stability is particularly simple to prove,
as we now demonstrate.

#proposition[
  Assume that there is $L_Phi > 0$ such that for all $t in [0, T]$ and $Delta > 0$,
  the function $vect(Phi)_Delta (t, placeholder)$ is globally Lipschitz continuous with constant~$L_Phi$.
  Then the one-step method~@eq:general_one_step is stable.
]
#proof[
  By~@eq:general_one_step and~@eq:general_one_step_perturbed,
  it holds that
  $
    vect(x)_n - vect(y)_n
    = vect(x)_(n-1) - vect(y)_(n-1)
    + Delta Big((vect(Phi)_Delta (t_(n-1), vect(x)_(n-1)) - vect(Phi)_Delta (t_(n-1), vect(y)_(n-1))))
    - Delta vect(delta)_n.
  $
  Taking the Euclidean norm and employing the Lipschitz continuity assumption,
  we obtain
  $
    vecnorm(vect(x)_n - vect(y)_n)
    <= (1 + Delta L_Phi) vecnorm(vect(x)_(n-1) - vect(y)_(n-1))
    + Delta vecnorm(vect(delta)_n).
  $
  By a reasoning similar to that in the proof of~@theorem:forward_euler,
  we then obtain
  $
    vecnorm(vect(x)_n - vect(y)_n)
    <= (1 + Delta L_Phi)^n vecnorm(vect(x)_0 - vect(y)_0) + sum_(i=1)^n (1 + Delta L_Phi)^(n-i) Delta vecnorm(vect(delta)_i)
    <= 0 + ee^(L_Phi t_n) Delta sum_(i=1)^n vecnorm(vect(delta)_i).
  $
  We conclude that~@eq:stability_numerical_method is satisfied with~$S(T) = ee^(L_Phi T)$.
]

==== Convergence
We are now ready to prove that consistency and stability of the numerical~@eq:general_one_step together imply convergence,
in the sense that
$
  lim_(Delta -> 0) lr((max_(1 <= n <= N) vecnorm(vect(x)(t_n) - vect(x)_n))) = 0, wide N = lr(floor.l T/Delta floor.r).
$
This result is an instance of the _Lax equivalence theorem_,
a pillar of numerical analysis with far-reaching applications.

#theorem(title: [Consistence and stability imply convergence])[
  Assume that the one-step numerical method~@eq:general_one_step is consistent and stable.
  Then the method is also convergent.
] <theorem:lax_convergence>
#proof[
  By definition of the local truncation error, it holds that
  $
    vect(x)(t_(n+1)) = vect(x)(t_n) + Delta vect(Phi)_Delta big((t_n, vect(x)(t_n))) + Delta vect(eta)_(n+1).
  $
  Therefore, the sequence $big((vect(x)(t_n)))_(1 <= n <= N)$ satisfies~@eq:general_one_step_perturbed with $vect(delta)_n = vect(eta)_n$,
  and so the stability estimate~@eq:stability_numerical_method implies that
  $
    max_(1 <= n <= N) vecnorm(vect(x)(t_n) - vect(x)_n) <= S(T) Delta sum_(n=1)^N vecnorm(vect(eta)_n).
  $
  By consistency, the right-hand side converges to zero in the limit as~$Delta -> 0$,
  which concludes the proof.
]

#remark[
  If we assume in @theorem:lax_convergence that the method is consistent with order~$p$,
  then by adapting the proof,
  we find that the error satisfies
  $
    max_(1 <= n <= N) vecnorm(vect(x)(t_n) - vect(x)_n) <= C S(T) Delta^p.
  $
  In this setting,
  the numerical scheme is said to be _convergent with order~$p$_.
]

=== Widely used one-step methods <sub:other_one_step>
In this section, we motivate and describe some of the other widely-used one-step methods,
namely methods of Taylor and Runge--Kutta type.
We assume in this section that the equation~@eq:ode admits a unique smooth solution over the interval~$[0, T]$.

==== Taylor methods
In order to construct a method with a smaller local truncation error than that of the forward Euler method,
a Taylor expansion of higher order than~@eq:taylor_forward can be employed:
#numbered($
  vect(x)(t + Delta) = vect(x)(t)
  + Delta vect(x)'(t)
  + cdots
  + Delta^p/p! vect(x)^((p))(t)
  + cal(O)(Delta^(p+1)).
$) <eq:taylor_high_order>
Since $vect(x) colon [0, T] -> real$ is a smooth solution to~@eq:ode by assumption,
the time derivatives of~$vect(x)$ can be obtained by differentiation of~@eq:ode:
$
  vect(x)'(t) = vect(f) big((t, vect(x)(t))),
  wide
  vect(x)''(t) = partial_t vect(f) big((t, vect(x)(t)))
  + Big((vect(f) lr((t, vect(x)(t))) dot nabla_(vect(x)))) vect(f) big((t, vect(x)(t))),
  wide
  dots
$
In general, it is immediate to show inductively that $vect(x)^((p))(t) = vect(f)^((p-1)) big((t, vect(x)(t)))$,
where the functions~$vect(f)^((p)) colon real times real^n -> real$ are defined recursively from the following equation:
$
  vect(f)^((p+1)) = partial_t vect(f)^((p)) big((t, vect(x)(t))) + Big((vect(f) lr((t, vect(x)(t))) dot nabla_(vect(x)))) vect(f)^((p)) big((t, vect(x)(t))).
$
The Taylor expansion~@eq:taylor_high_order motivates the so-called Taylor methods for integrating~@eq:ode numerically,
which are based on the following iteration:
#numbered($
  vect(x)_(n+1) = vect(x)_n + Delta vect(T)^p_Delta big((t_n, vect(x)_n)),
$) <eq:taylor_scheme>
where
$
  vect(T)^p_Delta (t, vect(x))
  := vect(f)(t, vect(x)) + Delta/(2!) vect(f)^((1))(t, vect(x))
  + cdots + Delta^(p-1)/p! vect(f)^((p-1))(t, vect(x)).
$
Note that, for any~$p$,
the Taylor scheme~@eq:taylor_scheme may be rewritten as
$
  vect(x)_(n+1) = vect(x)_n
  + Delta (dd vect(x)^(t_n, vect(x)_n))/(dd t) (t_n)
  + cdots
  + Delta^p/p! (dd^p vect(x)^(t_n, vect(x)_n))/(dd t^p) (t_n).
$
For~$p=1$, this scheme coincides with the forward Euler scheme.


==== Runge--Kutta methods
// Like Taylor methods, Runge--Kutta methods enjoy a local truncation error of high order,
Runge--Kutta methods resemble Taylor methods,
but they do not require to calculate the derivatives of the function~$vect(f)$.
This is achieved by approximating the derivatives in Taylor methods by difference quotients.
Consider for example the Taylor method of order $p = 2$:
#numbered($
  vect(x)_(n+1) = vect(x)_n
  + Delta (dd vect(x)^(t_n, vect(x)_n))/(dd t) (t_n)
  + Delta^2/(2!) (dd^2 vect(x)^(t_n, vect(x)_n))/(dd t^2) (t_n).
$) <eq:taylor_order_2>
Substituting the approximation
#eq-last($
  (dd^2 vect(x)^(t_n, vect(x)_n))/(dd t^2) (t_n)
  &approx 1/Delta lr(((dd vect(x)^(t_n, vect(x)_n))/(dd t) (t_n + Delta) - (dd vect(x)^(t_n, vect(x)_n))/(dd t)(t_n))) \
  &= 1/Delta Big((vect(f) big((t_n + Delta, vect(x)^(t_n, vect(x)_n)(t_n + Delta))) - vect(f) big((t_n, vect(x)_n)))) \
  &approx 1/Delta Big((vect(f) big((t_n + Delta, vect(x)_n + Delta vect(f)(t_n, vect(x)_n))) - vect(f) big((t_n, vect(x)_n))))
$) <eq:approx_heun>
in~@eq:taylor_order_2,
we obtain an explicit method known as _Heun's method_:
$
  vect(x)_(n+1) = vect(x)_n
  + Delta/2 vect(f)(t_n, vect(x)_n)
  + Delta/2 vect(f) big((t_n + Delta, vect(x)_n + Delta vect(f)(t_n, vect(x)_n))).
$
It is possible to show that
the local truncation error for this method also scales as~$Delta^2$.
Heun's method is a particular instance of a Runge--Kuta method.
In general, an explicit Runge--Kutta method with~$s$ stages is of the form
$
  vect(x)_(n+1) &= vect(x)_n + Delta sum_(i=1)^s b_i vect(k)_i \
  vect(k)_1 &= vect(f)(t_n, vect(x)_n), \
  vect(k)_2 &= vect(f) big((t_n + c_2 Delta, vect(x)_n + Delta(a_(21) vect(k)_1))), \
  vect(k)_3 &= vect(f) big((t_n + c_3 Delta, vect(x)_n + Delta(a_(31) vect(k)_1 + a_(32) vect(k)_2))), \
  &#h(0.56em) ⋮ \
  vect(k)_s &= vect(f) lr((t_n + c_s Delta, vect(x)_n + Delta sum_(j = 1)^(s-1) a_(s j) vect(k)_j)),
$
with appropriate coefficients $c_i$ and $a_(i j)$.
Heun's iteration can be recast in this form as follows:
$
  vect(x)_(n+1) &= vect(x)_n + Delta/2 (vect(k)_1 + vect(k)_2) \
  vect(k)_1 &= vect(f)(t_n, vect(x)_n) \
  vect(k)_2 &= vect(f) big((t_n + Delta, vect(x)_n + Delta vect(k)_1)).
$
The approach we employed to construct Heun's method may be generalized to higher orders.
For example, the most widely known Runge--Kutta method
approximates the Taylor method of order $p=4$ with the following iteration:
$
  vect(x)_(n+1) &= vect(x)_n + Delta/6 (vect(k)_1 + 2 vect(k)_2 + 2 vect(k)_3 + vect(k)_4), \
  vect(k)_1 &= vect(f)(t_n, vect(x)_n), wide
  && vect(k)_2 = vect(f) lr((t_n + Delta/2, vect(x)_n + Delta vect(k)_1/2)), \
  vect(k)_3 &= vect(f) lr((t_n + Delta/2, vect(x)_n + Delta vect(k)_2/2)), wide
  && vect(k)_4 = vect(f) lr((t_n + Delta, vect(x)_n + Delta vect(k)_3)).
$
The local truncation error for this method scales as~$Delta^4$ and,
when~$vect(f)(t, vect(x)) = vect(f)(t)$,
this method coincides with Simpson's formula~@eq:simpsons for the approximation of the integral in~@eq:one_step_integral.
The systematic derivation of Runge--Kutta methods is cumbersome,
and so we do not address this issue in this course.

#remark[
  Explicit Runge--Kutta methods of a given order are not uniquely defined.
  For example, if we employ instead of~@eq:approx_heun the approximation
  $
    (dd^2 vect(x)^(t_n, vect(x)_n))/(dd t^2) (t_n)
    approx 2/Delta lr((vect(f) lr((t_n + Delta/2, vect(x)_n + Delta/2 vect(f)(t_n, vect(x)_n))) - vect(f) big((t_n, vect(x)_n)))),
  $
  then we obtain by substitution in~@eq:taylor_order_2 the so-called _explicit midpoint method_,
  which is also a Runge--Kutta method of order 2:
  $
    vect(x)_(n+1) &= vect(x)_n + Delta vect(f) lr((t_n + 1/2 Delta, vect(x)_n + Delta/2 vect(f)(t_n, vect(x)_n))).
  $
]

==== Implicit methods
To conclude this section,
we mention two common implicit methods with a better order of convergence than that of the backward Euler method.
- The Crank--Nicolson method:
  #numbered($
    vect(x)_(n+1) = vect(x)_n + Delta/2 big((f(t_n, vect(x)_n) + f(t_n + Delta, vect(x)_(n+1)))).
  $) <eq:crank_nicolson>
  When $vect(f)$ is independent of~$vect(x)$ and depends only on~$t$,
  this method coincides with the trapezoidal rule for numerical integration.

- The implicit midpoint method:
  $
    vect(x)_(n+1) = vect(x)_n + Delta f lr((t_n + Delta/2, (vect(x)_n + vect(x)_(n+1))/2)).
  $
Similarly to the backward Euler method,
each iteration of these methods requires the resolution of a nonlinear equation.
Implicit methods often enjoy better stability than their explicit counterparts.
This subject is further discussed in~@sec:absolute_stability.

== Multistep methods <sec:multistep_methods>
The idea of multistep methods is to use,
in the construction of a new iterate,
information from not only the current but also previous iterations.
This degree of freedom enables to construct more economical numerical methods than one-step methods for the same order of convergence,
at the cost of a more difficult initialization.
In this section we focus on~_linear_ multistep methods of the form
#eq-last($
  vect(x)_(n+1) &= a_0 vect(x)_n + a_1 vect(x)_(n-1) + cdots + a_k vect(x)_(n-k) \
  &wide + Delta Big((b_(-1) vect(f)(t_(n+1), vect(x)_(n+1)) + b_0 vect(f)(t_n, vect(x)_n) + cdots + b_k vect(f)(t_(n-k), vect(x)_(n-k)))).
$) <eq:linear_mulitstep>
This equation defines an explicit method if $b_(-1) = 0$,
and an implicit method if $b_(-1) != 0$.
Note that explicit methods of the form~@eq:linear_mulitstep require only one additional evaluation~$f(t_n, vect(x)_n)$ per iteration,
in contrast with Runge--Kutta methods.
When~$b_(-1) != 0$, the iteration~@eq:linear_mulitstep is a nonlinear equation for the unknown~$vect(x)_(n+1)$,
which must itself be solved by resorting to a numerical method.

#paragraph[Initialization.][
  In order to initiate the numerical method~@eq:linear_mulitstep,
  the values $vect(x)_0, dots, vect(x)_k$ are required.
  These can be calculated by using a one-step method with an order of convergence matching that of the multistep method.
]

#paragraph[Local truncation error.][
  Consistently with the setting of one-step methods,
  the local truncation error for~@eq:linear_mulitstep is defined as the residual error left when the exact solution is substituted in the numerical scheme:
  #eq-last($
    Delta vect(eta)_(n+1)
    &:= vect(x)(t_n + Delta) - a_0 vect(x)(t_n) - a_1 vect(x)(t_n - Delta) - cdots - a_k vect(x)(t_n - k Delta) \
    & wide - Delta Big((b_(-1) vect(x)'(t_n + Delta) + b_0 vect(x)'(t_n) + cdots + b_k vect(x)'(t_n - k Delta))).
  $) <eq:truncation_multistep>
  The multistep method~@eq:linear_mulitstep is said to be of order~$p$ if the maximum local truncation error
  over all the discretization points, in norm, scales as~$bigo(Delta^p)$.
  The following result is useful for estimating the order of consistency of a linear multistep method.
]
#proposition[
  The linear multistep method~@eq:linear_mulitstep is consistent with order~$p$ for any smooth $vect(x) colon [0, T] -> real^n$ if and only if
  the local truncation error~@eq:truncation_multistep is everywhere zero when~$vect(x)(t)$ is of the scalar form
  #numbered($
    x(t) = t^q, wide q in {0, dots, p}.
  $) <eq:polynomial_form>
]
#proof[
  Assume that the method is consistent with order~$p$, fix $q in {1, dots, p}$,
  and let $x(t) = t^q$.
  Fix also $t in [0, T]$ and consider the function $xi colon {Delta : t\/Delta in nat_(> 0)} -> real$ given by
  $
    Delta xi(Delta) = Delta eta_((t\/Delta) + 1)
    &= x(t + Delta) - a_1 x(t) - a_2 x(t - Delta) - cdots - a_k x big((t - (k-1) Delta)) \
    & wide - Delta Big((b_0 x'(t + Delta) + b_1 x'(t) + cdots + b_k x' big((t - (k-1) Delta)))).
  $
  The quantity $xi(Delta)$ should be understood as the local truncation error at~$t$ for time step $Delta$.
  It is a polynomial in $Delta$ of degree at most~$p$ and scaling as~$bigo(Delta^(p+1))$.
  Therefore, it holds necessarily that~$xi(Delta) = 0$.

  Conversely, assume that the right-hand side of~@eq:truncation_multistep is equal to zero for any function of the form~@eq:polynomial_form.
  If~$vect(x)(t)$ denotes a smooth solution of~@eq:ode,
  then by Taylor's theorem there is $C > 0$ independent of $t_n$ such that
  $
    forall t in [0, T], wide
    system(
      vecnorm(vect(x)(t) - vect(y)(t)) &<= C|t-t_n|^(p+1);
      vecnorm(vect(x)'(t) - vect(y)'(t)) &<= C|t-t_n|^p
    ) ,
    wide
    vect(y)(t) := vect(x)(t_n) + sum_(i=1)^p vect(e)_i (t-t_n)^i,
  $
  for appropriate vectors $vect(e)_i in real^n$ depending on~$t_n$.
  Substituting $vect(x)(t) = vect(y)(t) + big((vect(x)(t) - vect(y)(t)))$ in the right-hand side of~@eq:truncation_multistep,
  we obtain
  $
    Delta vecnorm(vect(eta)_(n+1)) = bigo(Delta^(p+1)) + Delta bigo(Delta^p) = bigo(Delta^(p+1)),
  $
  with the constant implicit in the big~$bigo$ notation independent of~$n$.
  This concludes the proof.
]
#example[
  In the one-dimensional setting,
  we wish to find parameters $a_0$, $a_1$ and~$b_1$ such that the order of consistency
  of the following multistep scheme is as high as possible:
  $
    x_(n+1) = a_0 x_n + a_2 x_(n-1) + b_0 Delta f(t_n, x_n).
  $
  Substituting $x(t) = 1$ in the formula~@eq:truncation_multistep for the local truncation error,
  we obtain
  $
    eta_(n+1) = x(t_n + Delta) - a_0 x(t_n) - a_1 x(t_n - Delta) - b_0 Delta x'(t_n)
    = 1 - a_0 - a_1.
  $
  Therefore $a_1 = (1 - a_0)$.
  Next, substituting $x(t) = t - t_n$ in~@eq:truncation_multistep,
  we obtain
  $
    eta_(n+1) = Delta (2 - a_0 - b_0),
  $
  which gives $b_0 = 2 - a_0$.
  Finally, substituting $x(t) = (t-t_n)^2$,
  we obtain
  $
    eta_(n+1) = Delta^2 a_0,
  $
  and so $a_0 = 1$.
  We conclude that the best parameters,
  leading to a local truncation error scaling as~$bigo(Delta^2)$,
  are given by~$a_0 = 0$, $a_1 = 1$ and $b_0 = 2$.
  The resulting method reads
  $
    x_(n+1) = x_(n-1) + 2 Delta f(t_n, x_n),
  $
  and is known as the _multistep midpoint method_.
]
We now present two widely used systematic approaches for constructing multistep methods,
known as the Adams--Bashforth and Adams--Moulton approaches.

=== Adams--Bashforth methods
Let $vect(x) colon [0, T] -> real^n$ denote the exact solution to the differential equation~@eq:ode.
Integrating this equation between $t_n$ and $t_(n+1)$,
we obtain
#numbered($
  vect(x)(t_(n+1)) = vect(x)(t_n) + integral_(t_n)^(t_(n+1)) vect(f) big((t, vect(x)(t))) thin dd t.
$) <eq:adams_bashforth_intro>
The key idea of the Adams--Bashforth method is to approximate the function $t |-> vect(f) big((t, vect(x)(t)))$ by the interpolating polynomial~$hat(vect(f))$ of degree $k$ at the nodes
$t_(n-k), dots, t_n$:
#numbered($
  hat(vect(f))(t) = sum_(i=0)^k vect(f) big((t_(n-i), vect(x)(t_(n-i)))) L_i (t),
  wide L_i (t) := product_(j=0 \ j != i)^k (t - t_(n-j))/(t_(n-i) - t_(n-j)).
$) <eq:interpolation>
Substituting this approximation in~@eq:adams_bashforth_intro,
we obtain
$
  vect(x)(t_(n+1)) approx vect(x)(t_n) + sum_(i=0)^k vect(f) big((t_(n-i), vect(x)(t_(n-i)))) integral_(t_n)^(t_(n+1)) L_i (t) thin dd t.
$
This motivates the following _Adams--Bashforth_ numerical scheme:
#numbered($
  vect(x)_(n+1) = vect(x)_n + Delta sum_(i=0)^k b_i vect(f) big((t_(n-i), vect(x)_(n-i))),
  wide b_i := integral_0^1 product_(j=0 \ j != i)^k (s + j)/(-i + j) thin dd s.
$) <eq:adams_bashforth>
Since the Lagrange polynomials $(L_i)_(0 <= i <= k)$ depend on~$k$,
so do the coefficients $b_i$.
However, these are independent of~$Delta$,
and so they can be tabulated.
The value of these coefficients for the first few Adams--Bashforth methods are collated in~@table:adams_bashforth.

#figure(
  placement: top,
  kind: table,
  tabular(
    columns: 5,
    stretch: 1.2,
    [$i$], [$0$], [$1$], [$2$], [$3$],
    [$k=0$], [$1$], [], [], [],
    [$k=1$], [$3/2$], [$-1/2$], [], [],
    [$k=2$], [$23/12$], [$-16/12$], [$5/12$], [],
    [$k=3$], [$55/24$], [$-59/24$], [$37/24$], [$-9/24$],
  ),
  caption: [
    Coefficients $(b_i)_(0 <= i <= k)$ of the Adams--Bashforth methods.
  ],
) <table:adams_bashforth>

#paragraph[Local truncation error.][
  Assuming $vect(x) in C^(k+2) big(([0, T], real^n))$ and applying @theorem:interpolation_error for the interpolation error component-wise,
  we obtain
  $
    forall t in [0, T], wide lr(||vect(x)'(t) - hat(vect(f))(t)||)_infinity
    <= (|t - t_(n-k)| cdots |t - t_n|)/((k+1)!) sup_(t in [0, T]) lr(||vect(x)^((k+2))(t)||)_infinity,
  $
  where $hat(vect(f))$ is the function defined in~@eq:interpolation.
  Since
  $
    Delta vect(eta)_(n+1)
    = vect(x)(t_(n+1)) - vect(x)(t_n) - Delta sum_(i=0)^k b_i vect(f) big((t_(n-i), vect(x)(t_(n-i))))
    = integral_(t_n)^(t_(n+1)) Big((vect(x)'(t) - hat(vect(f))(t))) thin dd t,
  $
  we deduce that
  #numbered($
    vecnorm(vect(eta)_(n+1)) <= C_k M_(k+2) Delta^(k+1),
    wide M_(k+2) := sup_(t in [0, T]) lr(||vect(x)^((k+2))(t)||)_infinity,
  $) <eq:truncation_adams_bashforth>
  for an appropriate numerical constant~$C_k$ independent of~$n$ and of the problem data.
  Therefore the Adams--Bashforth method~@eq:adams_bashforth is consistent with order~$k+1$.
]

#paragraph[Convergence.][
  By using a reasoning similar to that in the proof of~@theorem:forward_euler,
  we can prove a convergence result of the Adams--Bashforth method.
]
#theorem[
  Assume that the solution $vect(x) colon [0, T] -> real^n$ to~@eq:ode is $k+2$ times continuously differentiable
  and that the global Lipschitz condition~@eq:global_lipschitz is satisfied.
  Suppose also that
  $
    forall i in {0, dots, k}, wide
    vecnorm(vect(x)(t_i) - vect(x)_i)
    <= delta.
  $
  Then the following error estimate holds for the Adams--Bashforth method~@eq:adams_bashforth:
  $
    forall n in lr({0, 1, dots, lr(floor.l T/Delta floor.r) }),
    wide
    vecnorm(vect(x)(t_n) - vect(x)_n)
    <=
    delta ee^(L B) + C_k M_(k+2) Delta^(k+1) lr(((ee^(L B t_n) - 1)/(L B))),
  $
  where $C_k$ and $M_(k+2)$ are the constants from~@eq:truncation_adams_bashforth,
  and with $B := abs(b_0) + cdots + abs(b_k)$.
] <theorem:adams_bashforth>
#proof-titled[Sketch of proof][
  Let $vect(e)_n := vect(x)(t_(n+1)) - vect(x)_(n+1)$.
  From the equation
  $
    vect(x)(t_(n+1)) - vect(x)_(n+1) =
    vect(x)(t_n) - vect(x)_n + Delta sum_(i=0)^k b_i Big((vect(f) big((t_(n-i), vect(x)(t_(n-i)))) - vect(f) big((t_(n-i), vect(x)_(n-i)))))
    + Delta vect(eta)_(n+1),
  $
  which is valid for $n >= k$,
  we deduce that
  $
    max Big(\{) vecnorm(vect(e)_0), dots, vecnorm(vect(e)_(n+1)) Big(\}) <= lr((1 + Delta L B)) max Big(\{) vecnorm(vect(e)_0), dots, vecnorm(vect(e)_n) Big(\}) + C_k M_(k+2) Delta^(k+2).
  $
  Since $max Big(\{) vecnorm(vect(e)_0), dots, vecnorm(vect(e)_k) Big(\}) <= delta$ by assumption,
  the statement easily follows.
]

=== Adams--Moulton methods
The Adams--Moulton methods are very similar to their Adams--Bashforth cousins.
The only difference is that the former are obtained by interpolating the function $t |-> vect(f) big((t, vect(x)(t)))$ in~@eq:adams_bashforth_intro at nodes shifted forward by~$Delta$,
i.e.~at the nodes $t_(n-k+1), dots, t_(n+1)$.
This leads to the method
#numbered($
  vect(x)_(n+1) = vect(x)_n + Delta sum_(i=-1)^(k-1) b_i vect(f) big((t_(n-i), vect(x)_(n-i))),
  wide b_i := integral_0^1 product_(j=-1 \ j != i)^(k-1) (s + j)/(-i + j) thin dd s.
$) <eq:adams_moulton>
Unlike the Adams--Bashforth methods, which are _explicit_,
the Adams--Moulton methods are _implicit_.
The value of the coefficients for the first few Adams--Moulton methods are collated in~@table:adams_moulton.
Notice that, for~$k = 0$, the Adams--Moulton method coincides with the backward Euler method,
and for $k = 1$ it coincides with the Crank--Nicolson method.
#figure(
  placement: top,
  kind: table,
  tabular(
    columns: 5,
    stretch: 1.2,
    [$i$], [$-1$], [$0$], [$1$], [$2$],
    [$k=0$], [$1$], [], [], [],
    [$k=1$], [$1/2$], [$1/2$], [], [],
    [$k=2$], [$5/12$], [$8/12$], [$-1/12$], [],
    [$k=3$], [$9/24$], [$19/24$], [$-5/24$], [$1/24$],
  ),
  caption: [
    Coefficients $(b_i)_(0 <= i <= k)$ of the Adams--Moulton methods.
  ],
) <table:adams_moulton>

== Absolute stability <sec:absolute_stability>
To conclude this chapter,
we introduce the notion of _absolute stability_ and explain its relevance.
Absolute stability is a property of a numerical method in relation to the model equation
#numbered($
  system(
    x'(t) &= lambda x(t)\,;
    x(0) &= 1.
  )
$) <eq:model_equation>
A numerical scheme for approximating~@eq:model_equation is called _absolutely stable_ if
#numbered($
  |x_n| -> 0 wide " in the limit as " n -> infinity.
$) <eq:convergence_absolute_stability>
where $(x_n)_(n=0, 1, dots)$ denotes the numerical solution to~@eq:model_equation.
Whether a numerical method is absolutely stable or not depends on the parameters~$lambda$ and $Delta$.
#example[
  The forward Euler method for~@eq:model_equation reads
  $
    x_(n+1) = x_n + Delta lambda x_n = (1 + Delta lambda) x_n.
  $
  Therefore $x_n -> 0$ if and only if~$|1 + Delta lambda| <= 1$.
] <example:absolute_stability_euler>
As @example:absolute_stability_euler illustrates,
whether absolute stability holds for the forward Euler methods depends only the value of the product $Delta lambda in complex$.
This dependence on $lambda$ and $Delta$ only through the product $Delta lambda$ holds in fact generally.
Indeed, all the numerical schemes we considered in this chapter are invariant under linear time rescaling of the ordinary differential equation:
the numerical solution of the rescaled equation, when the time step is rescaled by the same factor,
coincides with the discrete function obtained by linear rescaling of the numerical solution to the original equation.
This motivates the definition of _absolute stability region_ as
$
  cal(A) := { z in complex : #[@eq:convergence_absolute_stability holds when $Delta lambda = z$] } subset complex.
$
The exact solution to the model equation~@eq:model_equation diverges to $infinity$ as $t -> infinity$ if $Re(lambda) > 0$,
and it converges to $0$ if~$Re(lambda) < 0$.
Numerical schemes which exhibit a similar property at the discrete level are called _A-stable_.
More precisely, a numeric method is A-stable if the absolute stability region $cal(A)$ contains the left half-plane,
i.e.~if
$
  { z in complex : Re(z) < 0 } subset cal(A).
$
Before investigating whether the numerical schemes introduced previously in this chapter are absolutely stable,
we address the following natural question:
why focus on the simple model equation~@eq:model_equation?
We provide a couple of motivations:
// and is absolute stability of a numerical method useful outside of more general settings?
- First, note that equations of the form~@eq:model_equation are more relevant in science than might appear at first glance.
  Indeed, when discretizing in space a linear parabolic partial differential equation,
  one often obtains a linear differential equation of the form
  $
    vect(x)'(t) = mat(A) vect(x),
  $
  where $mat(A) in complex^(n times n)$.
  If the matrix $mat(A) = mat(Q) mat(D) mat(Q)^*$ is diagonalizable,
  then the vector $vect(z)(t) := mat(Q)^* vect(x)(t)$ satisfies the differential equation
  $
    vect(z)'(t) = mat(D) vect(z).
  $
  In other words, each component of~$vect(z)$ satisfies an ordinary differential equation of the same form as the model equation~@eq:model_equation.
  In applications, the components of~$vect(z)$ often encode the amplitudes of Fourier modes of the solution to the partial differential equation,
  and for dissipative equations all the eigenvalues of~$mat(A)$ have a negative real part.
  However, the spectral radius of~$mat(A)$ usually diverges as the number of discretization points increases.
  In this context, A-stability is particularly attractive,
  as it ensures that the numerical approximation remains well-behaved as the number of discretization points increases.

- Second, the model equation~@eq:model_equation may be viewed as a linearized approximation of a more interesting equation.
  Consider, for example, the following one-dimensional autonomous differential equation:
  #numbered($
    system(
      x'(t) &= f big((x(t)))\,;
      x(0) &= x_0.
    )
  $) <eq:nonlinear_autonomous>
  Assume that $f(x_*) = 0$ for some $x_* in real$.
  Such a point is called a _critical point_ of the differential equation.
  If $f'(x_*) < 0$,
  then $x_*$ is an attractor of the equation,
  in the sense that $x(t) -> x_*$ provided that $x(0)$ is sufficiently close to $x_*$.
  This result, which is the counterpart of~@proposition:local_convergence for differential equations,
  is a particular case of a theorem due to Poincaré and Lyapunov;
  see~@MR1422255[Theorem 7.1].
  If $|x_0 - x_*|$ is sufficiently small,
  then the solution to~@eq:nonlinear_autonomous is expected to be close
  to that of the linearized equation
  #numbered($
    system(
      y'(t) &= f' big((x_*)) big((y(t) - x_*))\,;
      y(0) &= x_0\,
    )
  $) <eq:linearized>
  which is of the form~@eq:model_equation.
  Often, studying the linearized equation~@eq:linearized enables to gain insight into the behavior of the original equation~@eq:nonlinear_autonomous,
  and analyzing the performance of a numerical method for the linearized equation~@eq:linearized is useful
  to inform the choice of a numerical scheme for~@eq:nonlinear_autonomous.

- More generally, if~$x(t)$ and $x_epsilon (t)$ are respectively the solutions to
  #numbered($
    system(
      x'(t) &= f big((t, x(t)))\,;
      x(0) &= x_0 + epsilon\,
    )
    wide
    " and "
    wide
    system(
      x_epsilon'(t) &= f big((t, x_epsilon (t)))\,;
      x_epsilon (0) &= x_0 + epsilon\,
    )
  $) <eq:multiscale_autonomous>
  then the difference $e(t) := x_epsilon (t) - x(t)$ satisfies the equation
  #numbered($
    e'(t) &= f big((t, x_epsilon (t))) - f big((t, x(t))) approx partial_x f big((t,x(t))) e(t), \
    e(0) &= epsilon,
  $) <eq:ode_error_equation>
  which looks similar to~@eq:model_equation with $partial_x f big((t,x(t)))$ in place of $lambda$.
  At a given time~$t$,
  the solutions tend to converge to each other as time increases
  if~$partial_x f big((t, x(t))) < 0$,
  and diverge from each other if $partial_x f big((t, x(t))) > 0$.
  Testing absolute stability with $lambda = partial_x f big((t, x(t)))$ enables to determine whether this property holds true also at the discrete level.
  Although the latter statement is difficult to state precisely and prove generally,
  we illustrate its validity for the forward Euler method in~@example:perturbed_forward_euler.

#example[
  Let $(x_n)$ and $(x_n^epsilon)$ denote the numerical solutions obtained by applying the forward Euler method to the differential equations in~@eq:multiscale_autonomous.
  If $epsilon << 1$,
  then
  $
    x_(n+1)^epsilon - x_(n+1)
    &= x_n^epsilon - x_n
    + Delta f(t_n, x^epsilon_n) - Delta f(t_n, x_n) \
    &approx x_n^epsilon - x_n + Delta partial_x f(t_n, x_n) (x_n^epsilon - x_n)
    = big((1 + Delta partial_x f(t_n, x_n))) (x_n^epsilon - x_n).
  $
  Therefore, the numerical solutions $(x^epsilon_n)$ and $(x_n)$ tend to become closer as~$n$ increases if
  #numbered($
    Delta partial_x f(t_n, x_n) in cal(A).
  $) <eq:stability_nonlinear>
] <example:perturbed_forward_euler>

#indent The absolute stability regions of the forward and backward Euler methods are illustrated in green in~@fig:stab_region_euler.
For the forward Euler method,
absolute stability holds if and only if~$|1 + Delta lambda| < 1$,
as we proved in~@example:absolute_stability_euler.
A similar reasoning gives that the absolute stability region for backward Euler method is given by $\{z in complex : |1 - z|^(-1) < 1\}$.
The backward Euler method is A-stable but the forward Euler method is not.
Notice that, if the time step is sufficiently large,
then the backward Euler method is absolutely stable even for values of $lambda$ with a positive real part,
for which exact solutions to the model equation are divergent.

#figure(
  placement: top,
  [
    #box(image("/figures/stability_forward Euler.pdf", width: 49%))
    #box(image("/figures/stability_backward Euler.pdf", width: 49%))
  ],
  caption: [Absolute stability regions for the forward (left) and backward (right) Euler methods.],
) <fig:stab_region_euler>

#example(title: [Absolute stability region of the Taylor methods])[
  When applied to~@eq:model_equation,
  the Taylor method of order~$p$ given in~@eq:taylor_scheme reads
  $
    x_(n+1) = lr((1 + Delta lambda + (Delta^2 lambda^2)/2 + cdots + (Delta^p lambda^p)/p!)) x_n.
  $
  Thus, the absolute stability region is given by
  $
    lr({ z in complex : lr(|1 + z + z^2/2 + cdots + z^p/p!|) < 1 }).
  $
  This region is illustrated for various values of~$p$ in~@fig:stab_taylor.
  We observe that the absolute stability region grows as $p$ increases.
]

#figure(
  placement: top,
  image("/figures/stability_Taylor_methods.pdf", width: 70%),
  caption: [Stability regions for the first few Taylor methods.],
) <fig:stab_taylor>

#heading(level: 3, numbering: none, outlined: false)[Stiff differential equations]
In the context of ordinary differential equations,
stiffness is not a precisely defined concept,
but rather rather a generic term employed to describe equations with widely separated time scales.
Roughly speaking, a differential equation of the form~@eq:ode is called stiff if the Jacobian matrix of~$f$,
with respect to the variable~$vect(x)$,
has at least one eigenvalue with a large negative real part.
In the one-dimensional setting,
the solutions to stiff differential equations which are close at the initial time tend to converge quickly to each other,
in view of~@eq:ode_error_equation.
This is illustrated in~@example:stiff_ode.

#example(title: [Stiff differential equation])[
  Consider the following equation~@Legat[Chapter 4]:
  #numbered($
    system(
      x'(t) &= - alpha big((x(t) - sin(t))) + cos(t);
      x(0) &= x_0
    )
  $) <eq:stiff_ode>
  The exact solution to this equation is given by
  $
    x(t) = sin(t) + x_0 thin ee^(-alpha t).
  $
  When $alpha in real$ is large,
  the distance between the solution and the function $t |-> sin(t)$ converges to zero very quickly,
  regardless of the initial condition.
  This behavior is illustrated in~@fig:stiff_ode.
] <example:stiff_ode>
#figure(
  placement: top,
  [
    #box(image("/figures/ode_stiff_2.pdf", width: 49%))
    #box(image("/figures/ode_stiff_5.pdf", width: 49%))
  ],
  caption: [Solutions to~@eq:stiff_ode for various initial conditions when $alpha = 2$ (left) and $alpha = 5$ (right).],
) <fig:stiff_ode>

#indent In the rest of this section,
we use the differential equation~@eq:stiff_ode as a guiding example.
For this problem, we have $partial_x f(t, x) = alpha$.
Therefore, in view of~@eq:stability_nonlinear,
we expect that the forward Euler scheme is non-divergent if $|1 - alpha Delta| < 1$,
i.e.~if
$
  Delta < Delta_* ~ = 2/alpha.
$
It turns out that this prediction is precise,
as depicted in~@fig:forward_euler_stiff.
Note that if the equation is very stiff,
that is to say if $alpha >> 1$,
then a very small time step is required to ensure stability.

In contrast with the forward Euler scheme,
the backward Euler scheme is stable regardless of the time step.
Since the right-hand side of~@eq:stiff_ode is linear in~$x$,
the value of the iterate~$x_(n+1)$ can be calculated explicitly from~$x_n$ for the backward scheme:
$
  x_(n+1) = (x_n + Delta alpha sin(t_(n+1)) + Delta cos(t_(n+1)))/(1 + Delta alpha).
$
Numerical approximations obtained using this scheme are illustrated in~@fig:backward_euler_stiff.
We observe that the method is stable even for the large time step $Delta = 2 Delta_*$.

#figure(
  placement: top,
  [
    #box(image("/figures/stiff_forward_euler_0_7.pdf", width: 48%))
    #box(image("/figures/stiff_forward_euler_0_97.pdf", width: 48%))
    #box(image("/figures/stiff_forward_euler_1_00.pdf", width: 48%))
    #box(image("/figures/stiff_forward_euler_1_03.pdf", width: 48%))
  ],
  caption: [
    Numerical approximations of the solution to~@eq:stiff_ode with $alpha = 10$ obtained with the forward Euler method,
    for four different values of $Delta$.
  ],
) <fig:forward_euler_stiff>

#figure(
  placement: top,
  [
    #box(image("/figures/stiff_backward_euler_0_5.pdf", width: 48%))
    #box(image("/figures/stiff_backward_euler_2_0.pdf", width: 48%))
  ],
  caption: [
    Numerical approximations of the solution to~@eq:stiff_ode with $alpha = 10$ obtained with the backward Euler method,
    for two different values of $Delta$.
  ],
) <fig:backward_euler_stiff>
== Exercises

#exercise[
  Show that the absolute stability region of the Crank--Nicolson method~@eq:crank_nicolson is given by the left half-plane;
  see~@fig:stab_region_crank_nicolson.
]

#figure(
  placement: top,
  image("/figures/stability_Crank-Nicolson.pdf", width: 70%),
  caption: [Absolute stability regions for the Crank Nicolson method.],
) <fig:stab_region_crank_nicolson>

#exercise[
  Calculate the absolute stability region for Gear's method.
]

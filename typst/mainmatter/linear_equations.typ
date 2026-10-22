#import "../template.typ": *

= Solution of linear systems of equation <cha:solution_of_linear_systems>
#minitoc()

#include "linear_equations/linear_introduction.typ"
#include "linear_equations/linear_conditioning.typ"
#include "linear_equations/linear_direct.typ"
#include "linear_equations/linear_iterative.typ"
#include "linear_equations/linear_exercises.typ"

== Discussion and bibliography

In this chapter,
we presented direct methods and some of the standard iterative methods for solving linear systems.
We focused particularly on linear systems with a symmetric positive definite matrix.
@sec:direct_solution_method is based on~@MR2265914 @Vuik
and @sec:iterative_methods roughly follows~@VanDooren[Chapter 2].
The book~@MR1990645 is a very detailed reference on iterative methods for solving sparse linear systems.
The reference~@shewchuk1994introduction is an excellent introduction to the conjugate gradient method.

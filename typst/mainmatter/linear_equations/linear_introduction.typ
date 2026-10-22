#import "../../template.typ": *

#section-star[Introduction]
This chapter is devoted to the numerical solution of linear problems of the following form:
#numbered($
  #[Find $vect(x) in complex^n$ such that] wide
  mat(A) vect(x) = vect(b),
  wide mat(A) in complex^(n times n),
  wide vect(b) in complex^n.
$) <eq:linear_system>
Systems of this type appear in a variety of applications.
They naturally arise in the context of linear partial differential equations,
which we use as main motivating example.
Partial differential equations govern a wide range of physical phenomena including heat propagation, gravity, and electromagnetism,
to mention just a few.
Linear systems in this context often have a particular structure:
the matrix $mat(A)$ is generally very sparse,
which means that most of the entries are equal to 0,
and it is often Hermitian and positive definite,
provided that these properties are satisfied by the underlying operator.

There are two main approaches for solving linear systems:
- Direct methods enable to calculate the exact solution to systems of linear equations,
  up to round-off errors, in a finite number of steps.
  Although this is an attractive property,
  direct methods are usually too computationally costly for large systems:
  The cost of inverting a general $n times n$ matrix,
  measured in number of floating operations,
  scales as $n^3$!

- Iterative methods, on the other hand,
  enable to progressively calculate increasingly accurate approximations of the solution.
  Iterations may be stopped once the _the residual_ is sufficiently small.
  These methods are often preferable when the dimension $n$ of the linear system is very large.

#indent This chapter is organized as follows.
- In @sec:conditioning,
  we introduce the concept of _conditioning_.
  The condition number of a matrix provides information on the sensitivity of the solution to perturbations of the right-hand side $vect(b)$ or matrix $mat(A)$.
  It is useful, for example, in order to determine the potential impact of round-off errors.

- In @sec:direct_solution_method,
  we present the direct method for solving systems of linear equations.
  We study in particular the $mat(L) mat(U)$ decomposition for an invertible matrix,
  as well as its variant for symmetric positive definite matrices,
  which is called the Cholesky decomposition.

- In @sec:iterative_methods,
  we present iterative methods for solving linear systems.
  We focus, in particular, on basic iterative methods based on a splitting of the matrix~$mat(A)$
  and on the conjugate gradient method.

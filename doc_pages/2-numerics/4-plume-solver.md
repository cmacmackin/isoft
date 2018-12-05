Title: Plume Solver
Author: Chris MacMackin
Date: December 2018 

For a one-dimensional domain,
[equations 8-12](./1-equations.html#mjx-eqn-eqplume-cont) become
second order ODEs. If boundary conditions are applied only at the grounding
line then this means that they can be solved using initial value
problem methods such as Runge-Kutta integration. However, the
diffusive terms in
[equations 8-12](./1-equations.html#mjx-eqn-eqnplume-mom-x) require two
boundary conditions each for velocity, salinity, and temperature. The
logical choice would use Dirichlet conditions at the grounding line
and outflow conditions (setting the first derivative to zero) at the
calving front. This turns the plume equations into a boundary
value problem and makes it more difficult to solve.
\( \def\bm#1{{\boldsymbol{#1}}}
\setcounter{27}
\)

One strategy for solving boundary value problems such as this is the
shooting method
[(see, e.g., Press et al.)](../6-bibliog.html#Press2007). With this
technique, unknown boundary values are guessed at the grounding line
to create an initial value problem. This initial value problem is then
integrated and the difference between the values at the calving front
and the required boundary conditions there is noted. A nonlinear
solver is then used to change the guesses at the grounding line in
order to get the correct boundary conditions at the calving
front. However, the diffusion term in these equations leads to
solutions involving exponential growth. If the guesses of grounding
line conditions are not sufficiently good then these exponentials can
lead to overflow and the failure of the solver. It was found that this
made the shooting method unsuitable in practice.

A relaxation method was also tried, wherein a time-dependent version
of the plume model was evolved
forward in time (using an explicit method) until it reached steady
state. However, the weakness of the diffusion coefficient and nearly
hyperbolic character of the time-dependent plume model means that
significant waves often arise, and considerable time is needed to
reach a steady state.

Finally, an approach called the quasilinearisation method
[(Mandelzweig and Tabakin, 2001)](../6-bibliog.html#Mandelzweig2001),
or QLM, was tried. This is a technique, based
on Newton's method, for solving boundary value problems. Though
[Mandelzweig and Tabakin (2001)](../6-bibliog.html#Mandelzweig2001)
present the technique for single-variable
problems, it is trivial to generalise it to the multivariate
case. Consider the differential equation
\begin{equation}
  \label{eq:qlm-equation}
  L^{(n)}\bm{s}(x) = \bm{f}[\bm{s}(x),
  \bm{s}^{(1)}(x), \ldots, \bm{s}^{(n-1)}(x), x], \qquad \bm{s}
  \in \mathbb{R}^{m},
\end{equation}
being solved on the domain \([0,b]\). Here, \(L^{(n)}\) is an \(n^{\rm th}\)
order linear differential operator, \(\bm{f}\) is a nonlinear function,
and \(\bm{s}^{(i)}\) is the \(i^{\rm th}\) derivative of
\(\bm{s}\). Boundary conditions are specified by
\begin{align}
  \label{eq:qlm-boundaries}
  g_k[\bm{s}(0), \bm{s}^{(1)}(0), \ldots, \bm{s}^{(n-1)}(0)] &= 0,& k = 1, \ldots, l; \nonumber \\
  g_k[\bm{s}(b), \bm{s}^{(1)}(b), \ldots, \bm{s}^{(n-1)}(b)] &= 0,& k = l+1, \ldots, mn;
\end{align}
where \(g_1, g_2, \ldots, g_{mn}\) are (potentially) nonlinear
functions. This can be solved iteratively for the \(r+1\) iterate
\(\bm{s}_{r+1}\) using the equation
\begin{multline}
  \label{eq:qlm-iter}
  L^{(n)}\bm{s}_{r+1}(x) = \bm{f}[\bm{s}_{r}(x),
  \bm{s}^{(1)}_{r}(x), \ldots, \bm{s}^{n-1}_{r}(x), x] + \\
  \sum_{i=0}^{n-1}\bm{f}_{\bm{s}^{(i)}}[\bm{s}_{r}(x),
  \bm{s}^{(1)}_{r}(x), \ldots, \bm{s}^{n-1}_{r}(x),
  x]\left[\bm{s}^{(i)}_{r+1}(x) - \bm{s}^{(i)}_{r}(x)\right],
\end{multline}
with boundary conditions for each iteration set by
\begin{align}
  \label{eq:qlm-boundaries-iter}
  \sum_{i=0}^{n-1}
  g_{k,\bm{s}^{(s)}}[\bm{s}_r(0), \bm{s}^{(1)}_r(0), \ldots,
  \bm{s}^{n-1}_r(0), x]\cdot[\bm{s}^{(i)}_{r+1}(0) -
  \bm{s}^{(i)}_r(0)] &= 0,&  k=1,\ldots,l; \nonumber \\
  \sum_{i=0}^{n-1}
  g_{k,\bm{s}^{(i)}}[\bm{s}_r(b), \bm{s}^{(1)}_r(b), \ldots,
  \bm{s}^{n-1}_r(b), x]\cdot[\bm{s}^{(i)}_{r+1}(b) -
  \bm{s}^{(i)}_r(b)] &= 0,&  k=l+1,\ldots,mn.
\end{align}
In these equations,
\(\bm{f}_{\bm{s}^{(i)}} = \partial \bm{f}/\partial \bm{s}^{(i)}\) (i.e.\
the Jacobian of \(\bm{f}\)) and
\(g_{k,\bm{s}^{(i)}} = \partial g_{k}/\partial \bm{s}^{(i)}\). Note
that, for linear boundary conditions,
[equation 31](#mjx-eqn-eqqlm-boundaries-iter) reduces to the boundary
conditions being constant across iterations.  This technique can be
proven to give quadratic convergence to the solution given certain
easily-satisfied assumptions
[(see Mandelzweig and Tabakin, 2001, for details)](../6-bibliog.html#Mandelzweig2001).
Furthermore, convergence is often monotonic.

To apply this method to the plume model, new variables \(U' = U_x\),
\(V' = V_x\), \(S' = S_x\), and \(T' = T_x\) were introduced, where a
subscript \(x\) denotes a derivative.
[Equation 8](./1-equations.html#mjx-eqn-eqplume-cont)
was rewritten so that the left-hand-side is just \(D_x\), while
[equations 9-12](./1-equations.html#mjx-eqn-eqplume-mom-x}) were rewritten
such that the left hand sides are \(U'_{x}\), \(V'_{x}\), etc. This
allowed a linear operator to be constructed with the form
\begin{equation}
  \label{eq:L-form}
  L{[D, U, U', S, S', T, T']}^T = {\left[\frac{d D}{d x}, 
      \frac{d U}{d x} - U',
      \frac{d U'}{d x},
      \frac{d S}{d x} - S',
      \frac{d S'}{d x},
      \frac{d T}{d x} - T',
      \frac{d T'}{d x}\right]}^T.
\end{equation}
The nonlinear operator is zero for \(U\), \(S\), and \(T\) and elsewhere
consists of the right-hand-side of the rearranged version of
[equations 8-12](./1-equations.html#mjx-eqn-eqplume-cont).

In order to find successive iterates, a linear equation must be
solved, consisting of the linear operator and the Jacobian
\(\bm{f}_{\bm{s}^{(i)}} = \partial \bm{f}/\partial \bm{s}^{(i)}\). It
is neither feasible nor efficient to explicitly evaluate the Jacobian,
especially if the solver is to be agnostic to parameterisation
choices. The iterative GMRES solver implemented in NITSOL (very
slightly modified to accept an initial guess of the solution) was used
because it can work knowing only the product of the Jacobian and the
current iterate. Initially these products were calculated using the
finite-difference approximation to the Jacobian
[described for the shelf solver](./3-shelf-solver.html). While this
was sufficiently accurate to run many simple simulations, it proved
unreliable when the plume undergoes a sudden change or when nonlinear
parameterisations are used.  To address this, automatic
differentiation [(Neidinger 2010)](../6-bibliog.html#Neidinger2010)
was applied instead and this proved far more robust. This calculated
the product of the Jacobian with a vector (i.e., the directional
derivative) via operator overloading. See the
[code design section](../3-codeDesign/2-discretisation.html) for
further details of the implementation. All results displayed in this
chapter were obtained using automatic differentiation.

The GMRES algorithm required preconditioning to work properly, as was
the case with the ice shelf solver. The preconditioner was chosen to
be \(\bm{P}^{-1}=L^{-1}\), equivalent to finding the inverse of
[equation 32](#mjx-eqn-eqL-form), which involves integration of the
derivatives. Spectral integration was performed by reversing the steps
for spectral differentiation described on the
[Spatial Discretisation](./2-spatial.html). A similarly modified
version of the NITSOL implementation of the biconjugate gradient
stabilised method (BiCGSTAB) was also found to work when solving the
preconditioned linear system, but it proved less robust.

When solving the linear equation at each iteration, the initial guess
was the previous iterate. The GMRES solver was expected to reduce the
error in the linear system by a factor \(\epsilon\) compared to the
initial guess. It was found that gradually decreasing the magnitude of
\(\epsilon\) over each nonlinear iteration, the final answer could be
determined with a residual norm smaller than \(7N \times 10^{-9}\),
where \(N\) is the number of grid points used and 7 indicates the number
of plume variables being solved for. The QLM proved to be highly
efficient, typically converging within a few iterations of
[equations 30 and 31](#mjx-eqn-eqqlm-iter),
although up to a few hundred iterations would often be needed by the
GMRES solver to perform the necessary intermediate linear solves.

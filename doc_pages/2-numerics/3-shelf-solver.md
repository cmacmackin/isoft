Title: Shelf Solver
Author: Chris MacMackin
Date: December 2018 

Simulating the evolution the ice shelf mass balance using
[equation 1](./1-equations.html#mjx-eqn-eqice-height-nd) requires a time-stepping scheme. In
order to allow numerical stability with large time steps, a
semi-implicit method is used. This works by defining the residual
operator \(
\def\bm#1{{\boldsymbol{#1}}}
\newcommand{\dx}{\mathcal{D}_x}
\setcounter{21}
\)
\begin{equation}
  \label{eq:shelf-resid}
  \bm{f}(h_n) =
  \frac{h_n - h_{n-1}}{\Delta t} +
  \frac{\partial}{\partial x}(h_{n}u_n(h_n)) + \lambda m_{n-1}
\end{equation}
where the subscript \(n\) indicates the value at the time step being
solved for, while subscript \(n - 1\) indicates the value at the
previous time step. This is a semi-implicit scheme (rather than fully
implicit) because melt rate \(m_{n-1}\) is used from the previous time
step, rather using \(m_n\) from the current time step. In this equation,
\(h_n\) and \(u_n\) represent vectors of thickness and velocity values at
each grid point at time step \(n\), while \(h_{n-1}\) is a vector of
thickness values at the previous time step and \(\partial/\partial x\)
is evaluated using the [Chebyshev differentiation procedure](./2-spatial.html)
previously described. This is a nonlinear system
and can be solved using Newton's method, where the root is determined
iteratively by solving the linear equation
\begin{equation}
  \label{eq:newton}
  \bm{J}\delta\bm{s}^k_n = -\bm{f}(\bm{s}^k_n).
\end{equation}
Here, \(\bm{s}_n = h_n\) is the current value of the ice thickness, \(\bm{J}\) is the Jacobian of \(f\),
and the new iterate \(\bm{s}^{k+1}_n = \bm{s}^{k}_n + \delta\bm{s}^k_n\).

In order to avoid having to evaluate the Jacobian of this system, a
Jacobian-free Newton-Krylov method
[(Knoll and Keyes, 2004)](../bibliog.html#Knoll2004) is used. This
solves the linear [equation 23](#mjx-eqn-eqnewton) iteratively via a Krylov
method which only requires the product of the Jacobian with the
iterate, and not the actual Jacobian itself. This product is
approximated as a finite difference:
\begin{equation}
  \label{eq:jac-fin-diff}
  \bm{J}\bm{v} \approx \frac{\bm{f}(\bm{s} + \epsilon\bm{v}) -
    \bm{f}(\bm{s})}{\epsilon}.
\end{equation}
The NITSOL implementation
[(Pernice and Walker, 1998)](../6-bibliog.html#Pernice1998) of a Newton-Krylov
solver was chosen, as it is very flexible and written in Fortran,
which was the language other portions of the code were to be
implemented with.

The spectral discretisation used here corresponds to dense matrices,
making [equation 23](#mjx-eqn-eqnewton) very poorly conditioned. As a
result, the Krylov solvers in NITSOL were unable to converge on a
solution without a preconditioner. Even with relatively-sparse
matrices, preconditioners are often needed for iterative methods
(e.g.: [Pernice and Walker, 1998](../6-bibliog.html#Pernice1998);
[Knoll and Keyes, 2004](../6-bibliog.html#Knoll2004)). A right preconditioner,
\(\bm{P}^{-1}\), is chosen so that the modified problem
\((\bm{J}\bm{P}^{-1})\bm{q} = -\bm{f}(\bm{s})\) is well-conditioned and
can be solved for \(\bm{q}\). It is then easy to find \(\bm{s}\) using
\(\bm{P}\delta\bm{s} = \bm{q} \Rightarrow \delta\bm{s} =
\bm{P}^{-1}\bm{q}\). A good preconditioner will have
\(\bm{P}^{-1} \approx \bm{J}^{-1}\), so that
\(\bm{J}\bm{P}^{-1} \approx \bm{I}\) to a decent approximation. A
tradeoff must be made between a preconditioner which is a sufficiently
good approximation of \(\bm{J}^{-1}\) to be useful and one which is not
too expensive or unstable to apply (e.g. \(P^{-1} = J^{-1}\) would be a
perfect preconditioner but constructing it is of equal difficulty to
solving the original, unpreconditioned problem).

The Jacobian of [equation 1](#mjx-eqn-eqshelf-resid) can be expressed as
\begin{equation}
  \label{eq:shelf-jacob}
  \bm{J} = \frac{1}{\Delta t} + \dx u,
\end{equation}
where we define \(\dx A \equiv \partial A/\partial x + A\Delta_x\), and
\(\Delta_x\) is the differential operator in the \(x\)-direction. Although
\(\Delta_x\) will be a dense matrix when using a spectral method, it can
be approximated as a sparse finite difference operator, as proposed
by [Orszag (1980)](../6-bibliog.html#Orszag1980). In this case
\(\Delta_x\), and thus also \(\dx\),
are tridiagonal matrices. This means that the finite difference form
of the Jacobian can be "inverted" simply by solving the tridiagonal
system, which can be done efficiently using a routine in LAPACK
[(Anderson et al., 1999)](../6-bibliog.html#Anderson1999).
This proved effective at preconditioning the
Krylov solver in NITSOL, whilst maintaining the accuracy of the
underlying pseudospectral method.

The \(u_n(h_n)\) term in [equation 1](#mjx-eqn-eqshelf-resid) can itself be
found by solving a nonlinear system, this time with the form
\begin{equation}
  \label{eq:u-resid}
  \bm{f}(u_n) =
  \frac{\partial}{\partial x}\left(4\eta_{n}h_{n}\frac{\partial
      u_n}{\partial x}\right) - \chi \frac{\partial h_n^2}{\partial x}.
\end{equation}
This has a Jacobian \begin{equation}
  \label{eq:u-jacob}
  \bm{J} = \dx(4\eta h)\Delta_x.
\end{equation}
Every time a new residual is calculated using
[equation 22](#mjx-eqn-eqshelf-resid),
[equation 26](#mjx-eqn-equ-resid) is solved iteratively using NITSOL.

It is also possible to construct a nonlinear system which takes both
\(h_n\) and \(u_n\) as arguments and solve for both simultaneously. While
this avoids the need to repeatedly solve for \(u_n\), it proved to be much
less stable and tended to require smaller time-steps in order to
prevent failure. As such, the approach outlined above proved to be
computationally cheaper overall.

Of the three linear solvers provided with NITSOL, only GMRES proved
reliable, with the BiCGSTAB and TFQMR solvers typically
failing. Consult [Pernice and Walker (1998)](../6-bibliog.html#Pernice1998) for more details on these
routines. It was also found that the "backtracking" globalisation
used in NITSOL, which is meant to prevent the solution from starting
to diverge, had a tendency to make the solver get trapped in local
minima for this problem. As such, it was turned off and this was found
to greatly improve the robustness of the nonlinear solver.

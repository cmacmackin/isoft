Title: Testing and Benchmarking
Author: Chris MacMackin
Date: December 2018 

## Ice Shelf

After programming the nonlinear solver using the algorithms described
in the previous section, various tests were run to ensure that it
would give the correct results. First, the ice shelf component was
tested with a prescribed melt rate matching that of the analytic
steady state solution in \S~\ref{sec:simplified-1-d}. It was
confirmed that when the ice shelf was initialised to the matching
steady state it remained there. Initialising the ice shelf to a
wedge-shape, it was found to evolve to the correct steady state.

As a test of the time-stepping for transient evolution, the 1-D shelf
equations were analysed for the special case where there is no
stretching (\)\chi = 0\)) and the melt rate \(m\) is constant in \(t\) and
\(x\). The velocity of the ice at the grounding line (and thus across
the entire shelf, since there is no stretching) was prescribed to be
\(u(t) = \bar{u} + \tilde{u}_0\sin(\omega t)\). With these assumptions,
equation~\eqref{eq:ice-height-nd} becomes
\begin{equation}
  \label{eq:ice-cont-nostretch}
  \frac{\partial h}{\partial t} - u(t)\frac{\partial h}{\partial x} =
  - \lambda m,
\end{equation}
which can be solved using the method of characteristics. In this
method, a Lagrangian coordinate \(s\) is introduced such that the
thickness of a parcel of ice following the trajectory \(x(s)\), \(t(s)\)
evolves according to \(h(s)\). It can then be shown that
\(\) \frac{d  x}{d  s} = u(t), \quad \frac{d 
  t}{d  s} = 1, \quad \frac{d  h}{d  s} = -\lambda
m. \(\) With the initial conditions \(x = \sigma\), \(t = 0\), and
\(h = h_0(\sigma)\) these equations can be integrated to yield the
transient solution
\begin{equation}
  \label{eq:h-transient-sol}
  h = h_0(\sigma) - \lambda mt,
\end{equation}
where \(\sigma\) can be computed from \(x\) and \(t\) according to
\(\sigma = x - \bar{u}t + \tilde{u}/\omega\left[\cos(\omega t) -
  1\right]\). This solution applies to ice starting in the domain, but
a different form is needed for ice parcels crossing the grounding
line at time \(t = t_g > 0\). Then \(\sigma < 0\) and the initial conditions are
set to \(x = 0\) and \(h = 1\). In this case the method of characteristics
provides the implicit solution
\begin{equation}
  \label{eq:h-implicit-sol}
  h - \frac{\lambda m\tilde{u}_0}{\omega\bar{u}}\cos\left[\omega t +
    \frac{\omega}{\lambda m}(h-1)\right] = 1 - \frac{\lambda
    m}{\bar{u}}\left[x + \frac{\tilde{u}_0}{\omega}\cos(\omega t) \right].
\end{equation}
This algebraic equation can easily be solved numerically for \(h\) using
a bisection-secant method, such as that of \citet[Chapter 4]{Brent1973}.
Possible solutions can be bracketed using the physical insight that
\(h \in [0,1]\).

These solutions provide a way to test accuracy of the ice shelf solver in
time and space. However, the fact the melt rate is constant means
that the semi-implicit approach to time-discretisation is not fully
tested. Using the same technique, a solution can be found for melt
rate \(m = m_t t\), where \(m_t\) is a constant rate of change in the
melt. For \(\sigma > 0 \( (calculated as before) the transient
solution applies:
\begin{equation}
  \label{eq:h-transient-sol2}
  h = h_0(\sigma) - \frac{\lambda m_t}{2}t^2.
\end{equation}
Elsewhere, the solution is again given implicitly:
\begin{equation}
  \label{eq:h-implicit-sol2}
  h - 1 + \lambda m_t (t_g + s/2)s = 0
\end{equation}
where the time since the ice parcel crossed the grounding line and the
time at which it crossed the grounding line are given by
\(\) s = \frac{x}{\bar{u}} + \frac{\tilde{u}_0}{\omega
  \bar{u}}\left(\cos(\omega t) - \cos(\omega t_g)\right), \quad t_g =
\sqrt{\frac{2(h - 1)}{\lambda m_t} + t^2}, \(\)
respectively.
Bracketing this solution is slightly more difficult than in the
constant-melt case, as if the value of \(h\) is too small it will result
in an imaginary value of \(t_g\). As such, the lower bound was set to
the value of \(h = 1 - \lambda m_t t^2/2\), which corresponds to \(t_g=0\) at the time being solved
for, plus some small value \(\epsilon\) to ensure that floating point
error does not become an issue. The upper bound remains set to 1.

A series of simulations were run under these conditions using
different time steps. Parameter values \(\bar{u}=1\), \(\tilde{u}_0=0.5\),
\(\lambda = 100\), \(m_t = 2\e{-4}\), and \(\omega=34.56\) were chosen,
corresponding to the scale choices described in the next section. A
domain of \(x \in [0,6]\) was used, with a wedge-shaped initial ice
profile \(h_0(x) = 1 - 0.1x\). All simulations used 300 grid-points.

\begin{figure}
  \centering
  \includegraphics[width=\textwidth]{benchmarking/err_comp_50.pdf}
  \caption[Comparison of the analytic solution given in
  equations~\eqref{eq:h-transient-sol2} and~\eqref{eq:h-implicit-sol2}
  with numerical solutions at \(t=5\).]{Comparison of the analytic
    solution given in equations~\eqref{eq:h-transient-sol2}
    and~\eqref{eq:h-implicit-sol2} with numerical solutions at \(t=5\),
    using 300 grid-points. These are plotted alongside each other in
    the top panel, while the bottom panel displays the differences
    between the numerical and analytical solutions. Insets offer a
    zoomed-in view of the error, with values at individual grid-points
    indicated by an \(\times\). Parameter values \(\bar{u}=1\),
    \(\tilde{u}_0=0.5\), \(\lambda = 100\), \(m_t = 2\e{-4}\), and
    \(\omega=34.56\) were chosen, corresponding to the scale choices
    described in the next section.}
  \label{fig:err-comp-5}
\end{figure}

\begin{figure}
  \centering
  \includegraphics[width=\fullwidth]{benchmarking/err_comp_100.pdf}
  \caption[Comparison of the analytic solution given in
  equations~\eqref{eq:h-transient-sol2} and~\eqref{eq:h-implicit-sol2}
  with numerical solutions at \(t=10\).]{Comparison of the analytic
    solution given in equation~\eqref{eq:h-implicit-sol2} with
    numerical solutions at \(t=10\), using 300 grid-points. In all
    cases, the plot shows the difference between the solution \(h\) and
    the unforced background state \(\bar{h}\) corresponding to the
    solution to equation~\eqref{eq:h-implicit-sol2} for
    \(\tilde{u}_0 = 0\). Only the first quarter of the domain is
    displayed, to make the plot easier to read; decay continues to
    contribute downstream. Parameter values \(\bar{u}=1\),
    \(\tilde{u}_0=0.5\), \(\lambda = 100\), \(m_t = 2\e{-4}\), and
    \(\omega=34.56\) were chosen, corresponding to the scale choices
    described in the next section.}
  \label{fig:err-comp-10}
\end{figure}

Figure~\ref{fig:err-comp-5} shows the results of two simulations at
time \(t=5\), with time-steps fixed at \(10^{-2}\) and \(10^{-4}\), compared
to the analytical solution given in equations~\eqref{eq:h-transient-sol2}
and~\eqref{eq:h-implicit-sol2}. Both simulations give reasonably good
agreement with the large-scale features of the solution, although
there is more significant error at the transition to the transient
solution. The numerical solutions tend to smooth out those sorts of
discontinuities, although reducing the time-step helps with this
considerably. The main issue, however, is how the numerical solution
handles the ripples which form due to the seasonal forcing of shelf
velocity. These are very small in magnitude, meaning that very high
levels of accuracy are demanded to resolve them. Even the simulation
with the smaller time-step shows signs of diffusion, causing the
ripples to loose amplitude as they move across the domain.

\begin{figure}
  \centering
  \includegraphics[width=\fullwidth]{benchmarking/err_dt.pdf}
  \caption[The error of simulations at times \(t=5\) and \(t=10\) for
  different time-steps.]{The root-mean-square (dashed) and the maximum
    (solid) error of simulations at times \(t=5\) and \(t=10\) for
    different time-steps, using 300 grid-points. Parameter values
    \(\bar{u}=1\), \(\tilde{u}_0=0.5\), \(\lambda = 100\), \(m_t = 2\e{-4}\),
    and \(\omega=34.56\) were chosen, corresponding to the scale choices
    described in the next section.}
  \label{fig:err-dt}
\end{figure}

This can be seen more clearly in figure~\ref{fig:err-comp-10}, which
is produced at \(t=10\) when the transient feature has been advected out
of the domain. All results in this plot are differences between the
time-dependent solution \(h\) with \(\tilde{u}_0 = 0.5\) and the
steady-state result \(\bar{h}\) of equation~\eqref{eq:h-implicit-sol2}
for \(\tilde{u}_0 = 0\), which corresponds to the unforced background
state. In order to make the plot easier to read, the domain only goes
to \(x=1.5\). As can be seen, the amplitude of the ripples decays,
indicating the presence of some numerical diffusion. Smaller
time-steps result in less of this diffusion. The root-mean-square and
the maximum error at times \(t=5\) and \(t=10\) were found for a range of
time-steps (figure~\ref{fig:err-dt}). The error declines fairly slowly
with the time-step. Given that high accuracy is needed for these
simulations, in any future developments of this algorithm it may be
useful to update the time-integration to a second-order or third-order
method, allowing larger time-steps to be used. The memory requirements
of this would be fairly modest, as only the ice thickness would need
to be saved between time-steps.

Similarly, the error was found to fall with an increasing number of
Chebyshev nodes used in the calculation. However, after a certain
point, the error stagnated and adding more nodes did not cause further
improvement. The point at which this stagnation occurs seems to depend
on the time-step, with smaller time-steps permitting higher numbers of
nodes before stagnation. Similarly, the beginnings of stagnation with
any further reductions in the time step can be seen in
figure~\ref{fig:err-dt}. This is consistent with the total error being
the sum of error arising due to temporal discretisation and spatial
discretisation. Increased resolution was found to be a more
computationally expensive means to improve accuracy than reducing the
time-step.

\begin{figure}
  \centering
  \includegraphics[width=\fullwidth]{benchmarking/err_evolution.pdf}
  \caption[Error over the course of two simulations using different
  time-steps.]{The root-mean-square (dashed) and the maximum (solid)
    error over the course of two simulations with different
    time-steps, both using 300 grid-points. Parameter values
    \(\bar{u}=1\), \(\tilde{u}_0=0.5\), \(\lambda = 100\), \(m_t = 2\e{-4}\),
    and \(\omega=34.56\) were chosen, corresponding to the scale choices
    described in the next section.}
  \label{fig:err-evolve}
\end{figure}

\begin{figure}
  \centering
  \includegraphics[width=\fullwidth]{benchmarking/err_evolution_const.pdf}
  \caption[Error over the course of two simulations with constant melt
  rate, using different time-steps.]{The root-mean-square (dashed) and
    the maximum (solid) error over the course of two simulations with
    different time-steps, both using 300 grid-points. In this case a
    constant melt rate was applied, with error found relative to the
    solution in equations~\eqref{eq:h-transient-sol}
    and~\eqref{eq:h-implicit-sol}. Parameter values \(\bar{u}=1\),
    \(\tilde{u}_0=0.5\), \(\lambda = 100\), \(m = 10^{-3}\), and
    \(\omega=34.56\) were chosen, corresponding to the scale choices
    described in the next section.}
  \label{fig:err-evolve-const}
\end{figure}

Plotting the root-mean-square and maximum error over the course of a
simulation shows that both grow approximately linearly
(figure~\ref{fig:err-evolve}), although the latter is very
noisy. Presumably this noise is due to aliasing of small-scale
features of the oscillations onto a discrete grid. The rate of growth
increases with the size of the time-step. There is a spike in error
which occurs as the kink at the transient feature reaches the end of
the domain around \(t=6\). The reason for the continued error growth
after the transient has been advected away is the growth in melt
rate. This means that ripples will tend to be larger and thus display
larger absolute error. Running a simulation with a constant melt rate
of \(m=10^{-3}\) and comparing to the solution in
equations~\eqref{eq:h-transient-sol} and~\eqref{eq:h-implicit-sol}
indicates that the error becomes roughly constant after the transient
feature leaves the domain (figure~\ref{fig:err-evolve-const}).

Experiments with this benchmarking problem showed that using 320
grid-points with a time step of \(10^{-5}\) resulted in absolute error
no larger than \(10^{-4}\). Error in the amplitude of the ripples at the
end of the domain was no more than \(\sim 10\%\), which was felt to be
acceptable when running simulations.


##Plume

Testing the plume solver was more difficult, as the structure of the
solver required a non-zero diffusivity, while the analytic solution in
equation~\eqref{eq:plume-bc-eps0} assumed \(\nu=0\). To avoid this
problem, the equation of state was altered for the benchmark test so
that it would always return the same density profile, regardless of
plume salinity or temperature. The density was chosen so that the
plume would have the same velocity as in the analytic solution. Now
uncoupled from the continuity and momentum
equations~\eqref{eq:plume-cont} and~\eqref{eq:plume-mom-x}, the
salinity and temperature equations~\eqref{eq:plume-salt}
and~\eqref{eq:plume-temp} could be analytically solved individually. A
plume was initialised by giving this analytic solution a sinusoidal
perturbation of amplitude 0.1 and a wavelength twice the size of the
domain. Starting from this initial guess and a prescribed wedge-shaped
ice thickness, the solver was able to converge to the expected result.

The coupled behaviour of the ice shelf/ocean received less rigorous
testing as there are no analytical benchmark solutions available for
the full nonlinear problem. The two components were initialised much
as they were in the plume test (except that the plume density was now
made dependent on salinity) and then allowed to evolve together. As
the resulting steady-state was not known, it was simply ensured that
numerical convergence was achieved as the number of Chebyshev nodes
increased and the time step reduced, and that the results looked
plausible.

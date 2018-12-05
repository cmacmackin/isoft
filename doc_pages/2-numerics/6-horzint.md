Title: Horizontally-Integrated Plume Equations
Author: Chris MacMackin
Date: December 2018 

The [plume solver](./4-plume-solver.html) described earlier assumed
that the plume is uniform in the transverse direction, with no
transverse velocity component. While this is a reasonable
approximation for narrow ice shelves, past modelling has shown that
the Coriolis force steers plume flow within subglacial cavities of
wide ice shelves Although 2-D plume models have been developed and
applied before
[(e.g. Sergienko, 2013)](../6-bibliog.html#Sergienko2013), they are
computationally expensive. Instead, ISOFT implements a
"horizontally_integrated" 1-D model, containing information on the
transverse flow. In addition to its computational simplicity, the
horizontally-integrated model provides a conceptual tool which can be
useful in understanding of the results of observations and more
complex simulations. The full derivation of this model can be found in
Chapter 4 of [MacMackin (2019)](../6-bibliog.html#MacMackin2019). What
follows is an overview of the results.

![A 3-D cartoon diagram of the horizontally-integrated
plume.](|media|/horizontal_integration.svg){: height=250px }
![A planar view of the horizontally-integrated plume.](|media|/horizontal_integration_planar.svg){: height=250px }

In this model, illustrated in figures above,
the plume variables are averaged over both the thickness of the plume
and also some lateral width \(\setcounter{37} \Delta y\) across the shelf. At the lower limit of
this domain in y is a sidewall of the subglacial cavity, through which
there can be no plume flow. The location of the upper limit is a
parameter which can be adjusted, but it is assumed to be an open
boundary through which transverse outflow is allowed. In order for
transverse flow to begin there must be something to break the
horizontal symmetry in the plume equations. This naturally arises due
to the Coriolis force. Simulations indicate that, in a rotational
plume such as this, there would be a narrow longitudinal boundary
current on the opposite side of the cavity. The presence of such a
boundary current is assumed here, rather than being explicitly
modelled; this current would act to drain the transverse flux of water
out from under the ice shelf.

The plume variables are assumed to be separable in \(x\) and \(y\), with the forms
\begin{equation}
  \label{eq:plume-var-sep}
  D(x,y) = \hat{D}(x)f_D(y), \quad U(x,y) = \hat{U}(x)f_U(y), \quad\ldots
\end{equation}
and similar for \(V\), \(T\), and \(S\).  A width-averaging operator,
represented by an over-bar, is defined according to
\begin{equation}
  \label{eq:averaging}
  \overline{G} = \frac{1}{\Delta y}\int^{y_2}_{y_1}G(y)dy,
\end{equation}
where \(G\) is an arbitrary _y_-dependent variable and  \(\Delta y = y_2 - y_1\).
The shape functions \(f(y)\) are defined such that
$$ \overline{f_D} = 1, \quad \overline{f_U} = 1, \quad\ldots $$
There is no general way to relate
\(\widehat{|\vec{U}|}(x)\) to \(\hat{U}(x)\) and \(\hat{V}(x)\), so instead
it is treated as an independent variable with its own shape function:
\begin{equation}
  \label{eq:Ubar-sep}
  |\vec{U}(x,y)| = \widehat{|\vec{U}|}(x)f_{|\vec{U}|}(y).
\end{equation}
However, \(\widehat{|\vec{U}|} = \sqrt{\hat{U}^2 + \hat{V}^2}\) is
exactly true if \(f_U(y) = f_V(y)\) or approximately true if \(U \gg V\)
or \(V \gg U\).

With these definitions in mind, the horizontally-integrated plume equations can be written as
  \begin{align}
    \alpha_{DU}\frac{d}{d x}\left(DU\right) +
    \left.\frac{f_D f_V}{\Delta y}\right\rvert^{y_2}_{y_1}DV &=
                                                               \overline{e} + \overline{m}, \label{eq:plume-cont-int} \\
    \alpha_{DU^2}\frac{d}{d x}\left(DU^2\right) +
    \left.\frac{f_D f_U f_V}{\Delta y}\right\rvert^{y_2}_{y_1}DUV &=
                                                                    D\rho_a\frac{d}{d x}\left(b - \delta\alpha_{D^2} D\right)\label{eq:plume-mom-x-int}\\ \nonumber
                                                             &- D\left(\overline{\rho}\frac{d b}{d x}
                                                               - \delta\alpha_{D^2}\tilde{\rho}\frac{d D}{d x}\right) \\ \nonumber
                                                             &+ \nu\alpha_{DU}\frac{d}{dx}\left(D\frac{dU}{dx}\right)
                                                               + \left.\nu\frac{DU}{\Delta
                                                               y}f_{D}f'_{U}\right\rvert^{y_2}_{y_1} \\ \nonumber
                                                             &- \mu\alpha_{|\vec{U}|U}|\vec{U}|U
                                                               + \frac{\delta\alpha_{D^2} D^2}{2}\frac{d\tilde{\rho}}{d x}
                                                               + \Phi\alpha_{DV}DV, \\
    \alpha_{DUV}\frac{d}{d x}\left(DUV\right) +
    \left.\frac{f_D f_V^2}{\Delta y}\right\rvert^{y_2}_{y_1}DV^2 &=
                                                                   \nu\alpha_{DV}\frac{d}{dx}\left(D\frac{dV}{dx}\right) 
                                                                   + \left.\nu\frac{DV}{\Delta y}f_{D}f'_{V}\right\rvert^{y_2}_{y_1} \label{eq:plume-mom-y-int}\\ \nonumber
                                                             &-\mu\alpha_{|\vec{U}|V}|\vec{U}|V - \left.\frac{\delta D^2}{2\Delta y}f_D^2[\rho_a
                                                               - \rho(x,y)]\right\rvert^{y_2}_{y_1}\\ \nonumber
                                                             &- \Phi\alpha_{DU}DU,
  \end{align}
  \begin{align}
    \alpha_{DUS}\frac{d}{d x}\left(DUS\right) +
    \left.\frac{f_D f_S f_V}{\Delta y}\right\rvert^{y_2}_{y_1}DSV &= \overline{e}S_a
                                                                    + \nu\alpha_{DS}\frac{d}{dx}\left(D\frac{dS}{dx}\right) \label{eq:plume-salt-int}\\ \nonumber
                                                                  &+ \left.\nu\frac{DS}{\Delta y}f_{D}f'_{S}\right\rvert^{y_2}_{y_1}
                                                                    + \overline{m}S_m - \overline{\gamma_S(S-S_m)}, \\
    \alpha_{DUT}\frac{d}{d x}\left(DUT\right) +
    \left.\frac{f_D f_T f_V}{\Delta y}\right\rvert^{y_2}_{y_1}DTV &= \overline{e}T_a
                                                                    + \nu\alpha_{DT}\frac{d}{d x}\left(D\frac{d T}{d x}\right) \label{eq:plume-temp-int}\\ \nonumber
                                                                  &+ \left.\nu\frac{DT}{\Delta y}f_{D}f'_{T}\right\rvert^{y_2}_{y_1}
                                                                    + \overline{m}T_m - \overline{\gamma_T(T-T_m)}.
  \end{align}
  
The constants involving \(\alpha\) are defined below. This result
assumes a linear equation of state, for which
\begin{equation}
  \overline{\rho} = \rho_{\rm ref}[1 + \beta_S(\alpha_{DS}S-S_{\rm
    ref}) - \beta_T(\alpha_{DT}T-T_{\rm ref})],\label{eq:rho-bar}
\end{equation}
and
\begin{equation}
  \tilde{\rho} = \rho_{\rm ref}[1 + \beta_S(\tilde{\alpha}_{DS}S-S_{\rm
    ref}) - \beta_T(\tilde{\alpha}_{DT}T-T_{\rm ref})].\label{eq:rho-tilde}
\end{equation}
The entrainment parameterisation in
[equation 14](./1-equations.html#mjx-eqn-eqentrainment-jenkins-nd) is unchanged. The one-equation
melt formulation of [equation 17](./1-equations.html#mjx-eqn-eqmelt-nondim) becomes
\begin{equation}
  \label{eq:horz-melt}
  \overline{m} = \zeta_1\zeta_2|\vec{U}|(\alpha_{|\vec{U}|T}T - T_m),
\end{equation}
when horizontally integrated. The ice is assumed to be impermeable to
salt, meaning \(\overline{\gamma_S(S - S_m)} = \overline{m}S_m =
0\). After horizontal integration, the thermal transfer term becomes
\begin{equation}
  \label{eq:horz-therm-trans}
  \overline{\gamma_T(T-T_m)} = \zeta_1|\vec{U}|(\alpha_{|\vec{U}|T}T - T_m).
\end{equation}

The \(\alpha\)
coefficients in these equations contain information on the transverse
shape of the plume variables and are defined as
\begin{align}
\newcommand{\horzint}[1]{\frac{1}{\Delta{}y}\int^{y_2}_{y_1}#1 dy}
  \nonumber
  \alpha_{DU} &= \overline{f_{D}f_{U}}, &
                                          \alpha_{DU^2} &= \overline{f_{D}f_{U}^2}, &
                                                                                      \alpha_{D^2} &= \overline{f_{D}^2},\\
  \nonumber
  \alpha_{DV} &= \overline{f_{D}f_{V}}, &
                                          \alpha_{DUV} &= \overline{f_{D}f_{U}f_{V}}, &
                                                                                        \alpha_{|\vec{U}|V} &= \overline{f_{|\vec{U}|}f_{V}},\\
  \label{eq:shape-coefs}
  \alpha_{|\vec{U}|U} &= \overline{f_{|\vec{U}|}f_{U}}, &
                                                          \alpha_{DUS} &= \overline{f_{D}f_{U}f_{S}}, &
                                                                                                        \alpha_{DUT} &= \overline{f_{D}f_{U}f_{T}},\\
  \nonumber
  \alpha_{|\vec{U}|T} &= \overline{f_{|\vec{U}|}f_{T}}, &
                                                          \quad\alpha_{DS} &= \overline{f_{D}f_{S}}, &
                                                                                                       \alpha_{DT} &= \overline{f_{D}f_{T}},\\
  \nonumber
  \tilde{\alpha}_{DS} &= \frac{\overline{f^2_D f_S}}{\alpha_{D^2}}, &
                                                                      \tilde{\alpha}_{DT} &= \frac{\overline{f_{D}^2f_{T}}}{\alpha_{D^2}}. &
\end{align}

The same [numerical methods](4-plume-solver.html) used to solve the
original set of plume equations could be used to solve
[equations 41-45](#mjx-eqn-eqplume-cont-int).  The linear
operator for the plume solver, defined in
[equation 32](./4-plume-solver.html#mjx-eqn-eqL-form), was
modified to contain the Coriolis forcing terms, becoming
\begin{multline}
  \label{eq:L-form-horzint}
  L{[D, U, U', S, S', T, T']}^T = \left[\frac{d D}{d x}, 
    \frac{d U}{d x} - U',
    \frac{d U'}{d x} - \frac{\Phi\alpha_{DV}}{\nu\alpha_{DU}} V,
    \frac{d V}{d x} - V',\right.\\\left.
    \frac{d V'}{d x} + \frac{\Phi\alpha_{DU}}{\nu\alpha_{DV}} U,
    \frac{d S}{d x} - S',
    \frac{d S'}{d x},
    \frac{d T}{d x} - T',
    \frac{d T'}{d x}\right]^T.
\end{multline}
Shape coefficients, drainage terms, and the equation for _y_-momentum
were added to the nonlinear operator. It was found that the existing
preconditioner was adequate to solve the modified equations. The
solver was tested first by running it in the trivial case with
\(\Phi = 0\) and \(V=0\) throughout the domain, with uniform thickness in
\(y\), to ensure that it converged to the solutions used for
[benchmarking](./5-benchmarking.html). It was then further tested by checking
that the values of each variable at the end of the domain agreed with
the asymptotic predictions described in
[MacMackin (2019)](../6-bibliog.html#MacMackin2019) when
\(\Phi \ne 0\) and \(b_x\) was constant.

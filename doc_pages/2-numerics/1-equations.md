Title: Shelf/Plume Equations
Author: Chris MacMackin
Date: December 2018 

The basic equations used in ISOFT to describe ice and plume dynamics
are similar to those of
[Sergienko (2013)](../6-bibliog.html#Sergienko2013), except that the
plume is assumed to be in quasi-steady state.  All equations were
nondimensionalised. Scales for this were chosen to be able to work
with multiple choices of parameterisations. This requirement of
flexibility sometimes resulted in otherwise suboptimal choices.

## Ice Shelf Equations

After rescaling, the dimensionless ice shelf equations have the form
\begin{equation}
    h_{t} + \nabla\cdot(h\vec{u}) = -\lambda m, \label{eq:ice-height-nd}
\end{equation}
\begin{equation}
    {\left[2\eta h\left(2u_x + v_y\right)\right]}_{x} + {\left[\eta
    h\left(u_y + v_x\right)\right]}_{y} - \chi{\left(h^2\right)}_x = 0, \label{eq:ice-mom-x}
\end{equation}
\begin{equation}
    {\left[\eta h\left(u_y + v_x\right)\right]}_{x} + {\left[2\eta
    h\left(u_x + 2v_y\right)\right]}_{y} - \chi{\left(h^2\right)}_y = 0. \label{eq:ice-mom-y}
\end{equation}

In these equations \(h\) is the ice thickness (scaled by reference
thickness \(h_0\)), \(\vec{u} = (u,v)\) is the velocity at which the ice
flows (scaled by reference \(u_0\)), \(m\) is the rate at which the ice
shelf is melting (rescaled by reference \(m_0\), defined below), and
\(\eta\) is the rescaled ice viscosity. The dimensionless parameters
\begin{equation}\label{eq:ice-parameters}
  \lambda \equiv \frac{\rho_0 m_0 x_0}{\rho_i h_0 u_0},
  \quad \chi \equiv \frac{\rho_i gh_0x_0}{2\eta_0 u_0}\left(1 - \frac{\rho_i}{\rho_0}\right)
\end{equation}
represent the ratio of melt versus influx of ice and the stretching rate
(ratio of gravitational stresses that drive stretching versus viscous
stresses resisting), respectively. The ice density is \(\rho_i\), while the reference density for
ocean water is \(\rho_0\). Otherwise, the subscript nought indicates a
typical scale for a variable. The timescale for ice flow is given by
\(t_0 = x_0/u_0\). Gravitational acceleration is \(g\). If viscosity is
modelled as being Newtonian, then the dimensionless \(\eta\) is set to 1.
Alternatively, Glen’s Law can be nondimensionalised to take the form
\begin{equation}\label{eq:glens-nd}
  \eta = \frac{1}{2}\xi D^{1/n - 1}_{2},
  \end{equation} where
\(D_2 = \sqrt{D_{i,j}D_{i,j}/2}\) is the second invariant of the strain
rate and the dimensionless coefficient
\begin{equation}\label{eq:glens-coef-nd}
  \xi \equiv \frac{B}{\eta_0} {\left(\frac{u_0}{x_0}\right)}^{1/n - 1}.
  \end{equation}
Typically, \(n=3\).


## Plume Equations

The plume equations are scaled according to
\begin{equation}\label{eq:scaling}
  \begin{gathered}
    U_0^2 = \frac{h_0g\Delta\rho}{\rho_0},
    \quad e_0=m_0=\gamma_{T0}=\gamma_{S0} = \frac{D_0U_0}{x_0}, \\
    \quad \Delta T_0 = \frac{\Gamma^{*}_T x_0}{D_0}(T_a - T_m),
    \quad \Delta S_0 = \frac{c_0\rho_0\Delta T_0(S_a - S_m)}{\rho_i L}.
  \end{gathered}
\end{equation}
As before, the subscript nought indicates the typical
scale for a variable. The exception to this is \(\rho_0\), which is a
representative value for the water density. Because the density
difference is used in the plume equation and this difference is quite
small, it was found useful to adopt a different scale, \(\Delta\rho\), to
use when nondimensionalising \(\rho_a - \rho\), \(\rho_x\), and \(\rho_y\).
Similarly, temperature and salinity were scaled according to typical
differences rather than by their absolute values. An arbitrary point
could be set to zero for these two variables and it proved convenient to
choose the ambient values as the zero-point. The basal depth, \(b\), is
scaled by \(h_0\). Note that the scale \(m_0\) does not correspond to
typical physical values of melting but is chosen because it is
convenient to have it equal to those of the other variables; as a
result, \(m \ll 1\). This yields the dimensionless system
\begin{equation}
    \nabla\cdot\left(D\vec{U}\right) = e + m, \label{eq:plume-cont}
\end{equation}
\begin{equation}
    \nabla\cdot\left(D\vec{U}U\right) = D(\rho_a - \rho)\left(b_x -
                                        \delta D_x\right) +
                                        \frac{\delta D^2}{2}\rho_x +
                                        \nu\nabla\cdot\left(D\nabla
                                        U\right) -
                                        \mu|\vec{U}|U, \label{eq:plume-mom-x}
\end{equation}
\begin{equation}
    \nabla\cdot\left(D\vec{U}V\right) = D(\rho_a - \rho)\left(b_y -
                                        \delta D_y\right) +
                                        \frac{\delta D^2}{2}\rho_y +
                                        \nu\nabla\cdot\left(D\nabla
                                        V\right) -
                                        \mu|\vec{U}|V, \label{eq:plume-mom-y}
\end{equation}
\begin{equation}
    \nabla\cdot\left(D\vec{U}S\right) = eS_a
                                        + \nu\nabla\cdot\left(D\nabla S\right) + mS_m -
                                        \gamma_S(S-S_m), \label{eq:plume-salt}
\end{equation}
\begin{equation}
    \nabla\cdot\left(D\vec{U}T\right) = eT_a
                                        + \nu\nabla\cdot\left(D\nabla T\right) + mT_m -
                                        \gamma_T(T-T_m). \label{eq:plume-temp}
\end{equation}

These equations were constructed
without making any assumptions about the form of \(m\), \(e\), or \(\rho\).
For this reason different, more generic, scales are adopted for these
values. The velocity scale depends on the density scale, rather than on
buoyancy input from subglacial discharge. Instead of scaling the
salinity in terms of buoyancy input, its scale is based the level of
melt-water input, which is the dominant source of salinity forcing
across most of the domain.

The dimensionless parameters
\begin{equation}\label{eq:plume-parameters}
  \delta \equiv \frac{D_0}{h_0},
  \quad r = \frac{\rho_0}{\rho_i},
  \quad \nu \equiv \frac{\kappa}{x_0U_0},
  \quad \mu \equiv \frac{C_dx_0}{D_0}
\end{equation}
represent the dimensionless
buoyancy correction, density ratio, turbulent eddy diffusivity, and
turbulent drag coefficient, respectively. The latter two depend on the
dimensional eddy diffusivity \(\kappa\), which is assumed to be equal to
the eddy viscosity,
and the unscaled turbulent drag coefficient \(C_d\).

The simple entrainment parameterisation of [Jenkins (1991)](../6-bibliog.html#Jenkins1991) can
be nondimensionalised to have the form
\begin{equation}
\label{eq:entrainment-jenkins-nd}
  e = \frac{E_0}{\delta}|\nabla b||\vec{U}|,
  \end{equation}
suggesting it is
convenient to take \(\delta = E_0\) (or equivalently, \(D_0 = E_0 h\)). The more complex parameterisation of
[Kochergin (1987)](../6-bibliog.html#Kochergin1987) nondimensionalises to
\begin{equation}
\label{eq:entrainment-koch-nd}
  e = \frac{K}{S_m}\sqrt{|\vec{U}|^2 + \frac{\delta(\rho_a - \rho)D}{S_m}},
  \end{equation}
with the dimensionless coefficient
\begin{equation}\label{eq:ent-koch-coef-nd}
  K = \frac{c_L^2 x_0}{D_0}.
\end{equation}
The turbulent Schmidt number, \(S_m\),
depends on the Richardson number:
$$ S_m = \frac{Ri}{0.0725(Ri + 0.186 - \sqrt{Ri^2 - 0.316Ri 0.0346})}. $$
With these scales, the Richardson number is given by
\(Ri = \delta(\rho_a - \rho)D/|\vec{U}|^2\). The simplified melt rate
parameterisation taken from
[Dallaston, Hewitt, and Wells (2015)](../6-bibliog.html#Dallaston2015),
has the dimensionless form
\begin{equation}\label{eq:melt-nondim}
  m = \zeta_1\zeta_2|\vec{U}|(T - T_m),
\end{equation}
where
\begin{equation}\label{eq:thermal-zetas}
  \zeta_1 = \frac{\Gamma^{*}_T x_0}{D_0}, \quad \zeta_2 =
  \frac{c\Delta T_0}{L}.
\end{equation}

Due to the low efficiency of thermal transfer
to the ice shelf, compared to the high rate of entrainment,
\(\zeta_1 \ll 1\). The large latent heat of ice results in
\(\zeta_2 \ll 1\), as well. Together, these results mean \(m \ll e \sim 1\),
so that the mass gain by meltwater input is much smaller than by
entrainment. It can be seen that the thermal transfer coefficient
nondimensionalises to give
\begin{equation}\label{eq:therm-trans-nondim}
  \gamma_T = \zeta_1 |\vec{U}|.
  \end{equation}

The water density is set using a linear equation of state of the form
\begin{equation}
  \label{eq:lin-eos}
  \rho_w = \rho_{\rm ref}[1 + \beta_S(S-S_{\rm ref})-\beta_T(T-T_{\rm ref})].
\end{equation}
Here, \(\beta_S\) is the haline contraction coefficient,
\(\beta_T\) is the thermal expansion coefficient, and \(S_{\rm ref}\),
\(T_{\rm ref}\), and \(\rho_{\rm ref}\) are reference values for salinity,
temperature, and density about which the relation has been
linearised.


## Typical Scales and Parameter Values

Typical scales and values for ice shelf and plume properties are
listed in the table below, along with the values of non-dimensional
parameters which result. “Repr.  val.” stands for “representative
value”. [J11] refers to [Jenkins (2011)](../6-bibliog.html#Jenkins2011),
[B11] to [Bindschadler, Vaughan, and Vornberger (2011)](../6-bibliog.html#Bindschadler2011),
[D15] to [Dallaston, Hewitt, and Wells (2015)](../6-bibliog.html#Dallaston2015),
[J91] to  [Jenkins (1991)](../6-bibliog.html#Jenkins1991),
[J96] to [Jacobs, Hellmer, and Jenkins (1996)](../6-bibliog.html#Jacobs1996),
[K87] to [Kochergin (1987)](../6-bibliog.html#Kochergin1987), and
[S13] to [Sergienko (2013)](../6-bibliog.html#Sergienko2013).
Scales in the third column are chosen to be comparable
to the conditions of the PIG ice shelf and come from the indicated
source. Where a scaling is unconstrained it was chosen to provide
convenient parameter values (e.g. \(x_0\) fixed by \(\chi\)). The value of
\(c_L\) was chosen so that its entrainment parameterisation resulted
in a similar melt rate as that of [@Jenkins1991]. Due to an error,
\(\zeta_2\) is a factor of \(r\) too large.
\(\newcommand{\unit}[1]{{\rm\thinspace #1}}
\newcommand{\m}{\unit{m}}
\newcommand{\km}{\unit{km}}
\newcommand{\kg}{\unit{kg}}
\newcommand{\s}{\unit{s}}
\newcommand{\yr}{\unit{yr}}
\newcommand{\J}{\unit{J}}
\newcommand{\K}{\unit{K}}
\newcommand{\Pa}{\unit{Pa}}
\newcommand{\hPa}{\unit{hPa}}
\newcommand{\kPa}{\unit{kPa}}
\newcommand{\rad}{\unit{rad}}
\newcommand{\psu}{\unit{psu}}
\newcommand{\mps}{\m\s^{-1}}
\newcommand{\mpss}{\m\s^{-2}}
\newcommand{\mpyr}{\m\yr^{-1}}
\newcommand{\kmpyr}{\km\yr^{-1}}
\newcommand{\C}{\rm ^{\circ}\thinspace C}
\newcommand{\kgpmc}{\kg\m^{-3}}
\newcommand{\ppsu}{\psu^{-1}}
\newcommand{\msps}{\m^{2}\s^{-1}}
\newcommand{\Pas}{\Pa\s}
\newcommand{\Jpkg}{\J\kg^{-1}}
\newcommand{\JpkgpK}{\Jpkg\K^{-1}}\)

<table class='table table-sm table-striped'>
<thead>
<tr>
<th>Variable</th>
<th>Description</th>
<th>Value</th>
<th>Source</th>
</tr>
</thead>
<tbody>
<tr>
<td>
<script type="math/tex">\rho_0</script>
</td>
<td>Reference water density</td>
<td>
<script type="math/tex">1030\kgpmc</script>
</td>
<td>Common</td>
</tr>
<tr>
<td>
<script type="math/tex">\rho_i</script>
</td>
<td>Ice density</td>
<td>
<script type="math/tex">916\kgpmc</script>
</td>
<td>Common</td>
</tr>
<tr>
<td>
<script type="math/tex">g</script>
</td>
<td>Acceleration due to gravity</td>
<td>
<script type="math/tex">9.8\mpss</script>
</td>
<td>Common</td>
</tr>
<tr>
<td>
<script type="math/tex">L</script>
</td>
<td>Latent heat of fusion</td>
<td>
<script type="math/tex">3.35\times 10^{5}\Jpkg</script>
</td>
<td>Common</td>
</tr>
<tr>
<td>
<script type="math/tex">c</script>
</td>
<td>Specific heat of water</td>
<td>
<script type="math/tex">3.98\times 10^{3}\JpkgpK</script>
</td>
<td>Common</td>
</tr>
<tr>
<td>
<script type="math/tex">E_0</script>
</td>
<td>Entrainment coefficient (J91)</td>
<td>0.036</td>
<td>[J11]</td>
</tr>
<tr>
<td>
<script type="math/tex">c_L</script>
</td>
<td>Entrainment coefficient (K87)</td>
<td>0.1059</td>
<td></td>
</tr>
<tr>
<td>
<script type="math/tex">\beta_S</script>
</td>
<td>Haline contraction coefficient</td>
<td>
<script type="math/tex">7.86\times 10^{-4}\psu^{-1}</script>
</td>
<td>[J11]</td>
</tr>
<tr>
<td>
<script type="math/tex">\beta_T</script>
</td>
<td>Thermal contraction coefficient</td>
<td>
<script type="math/tex">3.87\times 10^{-5}\K^{-1}</script>
</td>
<td>[J11]</td>
</tr>
<tr>
<td>
<script type="math/tex">C_d</script>
</td>
<td>Turbulent drag coefficient</td>
<td>
<script type="math/tex">2.5\times 10^{-3}</script>
</td>
<td>[J11]</td>
</tr>
<tr>
<td>
<script type="math/tex">\Gamma^{*}_T</script>
</td>
<td>Thermal transfer coefficient</td>
<td>
<script type="math/tex">5.7\times 10^{-5}</script>
</td>
<td>[D15]</td>
</tr>
<tr>
<td>
<script type="math/tex">\kappa</script>
</td>
<td>Turbulent diffusivity/viscosity</td>
<td>
<script type="math/tex">10</script>–<script type="math/tex">100\m^2\s^{-1}</script>
</td>
<td>Repr. val.</td>
</tr>
<tr>
<td>
<script type="math/tex">\eta_0</script>
</td>
<td>Ice viscosity</td>
<td>
<script type="math/tex">2.6\times 10^{13}\Pa\s</script>
</td>
<td>Repr. val.</td>
</tr>
<tr>
<td>
<script type="math/tex">B</script>
</td>
<td>Glen’s Law coefficient</td>
<td>
<script type="math/tex">1.6\times 10^{8}\Pa\s^{1/3}</script>
</td>
<td>[S13]</td>
</tr>
<tr>
<td>
<script type="math/tex">S_g</script>
</td>
<td>Subglacial discharge salinity</td>
<td>
<script type="math/tex">0\psu</script>
</td>
<td>Repr. val.</td>
</tr>
<tr>
<td>
<script type="math/tex">S_a</script>
</td>
<td>Ambient salinity</td>
<td>
<script type="math/tex">34.6\psu</script>
</td>
<td>[J96]</td>
</tr>
<tr>
<td>
<script type="math/tex">T_a - T_m</script>
</td>
<td>Thermal Forcing</td>
<td>
<script type="math/tex">2\C</script>
</td>
<td>[J96]</td>
</tr>
<tr>
<td>
<script type="math/tex">u_0</script>
</td>
<td>Ice velocity scale</td>
<td>
<script type="math/tex">2.5\kmpyr</script>
</td>
<td>[B11]</td>
</tr>
<tr>
<td>
<script type="math/tex">h_0</script>
</td>
<td>Ice thickness scale</td>
<td>
<script type="math/tex">1200\m</script>
</td>
<td>[B11]</td>
</tr>
<tr>
<td>
<script type="math/tex">x_0</script>
</td>
<td>Length scale</td>
<td>
<script type="math/tex">13.8\km</script>
</td>
<td></td>
</tr>
<tr>
<td>
<script type="math/tex">t_0</script>
</td>
<td>Time scale</td>
<td>
<script type="math/tex">5.5\yr</script>
</td>
<td></td>
</tr>
<tr>
<td>
<script type="math/tex">m_0</script>
</td>
<td>Melt scale</td>
<td>
<script type="math/tex">1.94\times 10^{4} \mpyr</script>
</td>
<td></td>
</tr>
<tr>
<td>
<script type="math/tex">X</script>
</td>
<td>Dimensionless domain length</td>
<td>
<script type="math/tex">6</script>
</td>
<td>[B11]</td>
</tr>
<tr>
<td>
<script type="math/tex">Q_g</script>
</td>
<td>Subglacial discharge</td>
<td>
<script type="math/tex">8.5\times 10^{-3}\m^2\s^{-1}</script>
</td>
<td>Repr. val.</td>
</tr>
<tr>
<td>
<script type="math/tex">D_0</script>
</td>
<td>Plume thickness scale</td>
<td>
<script type="math/tex">43.2\m</script>
</td>
<td></td>
</tr>
<tr>
<td>
<script type="math/tex">U_0</script>
</td>
<td>Plume velocity scale</td>
<td>
<script type="math/tex">0.196\mps</script>
</td>
<td></td>
</tr>
<tr>
<td>
<script type="math/tex">\Delta T_0</script>
</td>
<td>Temperature scale</td>
<td>
<script type="math/tex">0.0364\K</script>
</td>
<td></td>
</tr>
<tr>
<td>
<script type="math/tex">\Delta S_0</script>
</td>
<td>Salinity scale</td>
<td>
<script type="math/tex">0.0170\psu</script>
</td>
<td></td>
</tr>
<tr>
<td>
<script type="math/tex">\Delta\rho</script>
</td>
<td>Density variation scale</td>
<td>
<script type="math/tex">3.38\times 10^{-3}\kgpmc</script>
</td>
<td></td>
</tr>
<tr>
<td>
<script type="math/tex">\chi</script>
</td>
<td>Dimensionless stretching rate</td>
<td>
<script type="math/tex">4</script>
</td>
<td></td>
</tr>
<tr>
<td>
<script type="math/tex">\xi</script>
</td>
<td>Dimensionless Glen’s coefficient</td>
<td>
<script type="math/tex">1.919</script>
</td>
<td></td>
</tr>
<tr>
<td>
<script type="math/tex">\lambda</script>
</td>
<td>Dimensionless melt rate</td>
<td>
<script type="math/tex">100</script>
</td>
<td></td>
</tr>
<tr>
<td>
<script type="math/tex">r</script>
</td>
<td>Density ratio</td>
<td>
<script type="math/tex">1.12</script>
</td>
<td></td>
</tr>
<tr>
<td>
<script type="math/tex">\nu</script>
</td>
<td>Dimensionless eddy diffusivity</td>
<td>
<script type="math/tex">3.69\times 10^{-2}</script>
</td>
<td></td>
</tr>
<tr>
<td>
<script type="math/tex">\mu</script>
</td>
<td>Dimensionless drag coefficient</td>
<td>
<script type="math/tex">0.799</script>
</td>
<td></td>
</tr>
<tr>
<td>
<script type="math/tex">\delta</script>
</td>
<td>Buoyancy correction</td>
<td>
<script type="math/tex">0.036</script>
</td>
<td></td>
</tr>
<tr>
<td>
<script type="math/tex">K</script>
</td>
<td>Dimensionless K87 coefficient</td>
<td>
<script type="math/tex">3.58</script>
</td>
<td></td>
</tr>
<tr>
<td>
<script type="math/tex">\zeta_1</script>
</td>
<td>Dimensionless transfer coefficient</td>
<td>
<script type="math/tex">0.0182</script>
</td>
<td></td>
</tr>
<tr>
<td>
<script type="math/tex">\zeta_2</script>
</td>
<td>Dimensionless melt coefficient</td>
<td>
<script type="math/tex">4.86\times 10^{-4}</script>
</td>
<td></td>
</tr>
</tbody>
</table>

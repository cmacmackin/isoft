Title: Spatial Discretisation
Author: Chris MacMackin
Date: December 2018 

The shelf/plume simulation requires computing various derivatives, for
which a pseudospectral method is used. An introduction to this
technique is provided below; for a more thorough explanation, see
\citet{Trefethen2000}. Spectral methods provide a fast and accurate way to numerically
differentiate discrete data. While more
computationally expensive than finite difference methods for the same
number of grid points, spectral methods give exponential convergence
and thus often require significantly fewer grid points to achieve the
same level of accuracy. Numerical accuracy is of particular importance
here, as the purpose of running simulations is to test the stability
of an ice-shelf to potentially small perturbations. Spectral methdos are often used for
problems with periodic boundary conditions, where a Fourier
series, \(f(\theta) = \sum_k a_k e^{ik\theta}\), can be used to
interpolate between grid points. If the grid points are evenly spaced
then the coefficients \(a_k\) can easily be calculated with a discrete
Fourier transform. Typically this would be done using the highly
efficient fast Fourier transform (FFT) algorithm \citep{Cooley1965},
which requires \(O(N\log N)\) operations for N grid-points. The
derivative is then \(f'(\theta) = \sum_k ika_k e^{ik\theta}\) and an
inverse FFT can be used to convert the new coefficients \(ika_k\) to the
values of \(f'\) at each grid point.

However, the boundary conditions for the ice shelf and plume
are not periodic. Instead, say there is an
interpolant \(F(x)\) for data mapped onto \(-1 \le x \le 1\) using a
linear coordinate rescaling. To apply a spectral
method, it is necessary to map the interpolant to a function
\(f(\theta)\), \(0 \le \theta < 2\pi\), where \(x = \cos\theta\). Regardless
of the boundary conditions on \(F\), \(f\) will be periodic and even in
\(\theta\) and can thus be differentiated as before. The results can
then be mapped back onto the grid points in the \(x\)-domain. By
choosing \(x\) grid points to be \textit{Chebyshev collocation points},
defined below, the corresponding grid points in \(\theta\) will be
equally spaced and an FFT can be used to find the Fourier
coefficients. This is known as the \textit{Chebyshev pseudospectral
  method}~\citep{Trefethen2000}. If \(N + 1\) Chebyshev collocation
points are needed, their positions are given by
\begin{equation}
  \label{eq:cheb-colloc}
  x_j = \cos(j\pi/N), \quad j = 0, \ldots,N.
\end{equation}
This approach provides uneven spacing of points in \(x\), with a
clustering of resolution near the domain boundaries, and hence is also
well suited to capture rapid variation near the grounding line.

Following \citet{Trefethen2000}, the practical algorithm used to
differentiate discrete data \(v_j = v(x_j)\), for
\(0 \le j \le N\), corresponding to values at Chebyshev collocation
nodes \(x_0=1,\ldots,x_N=-1\), is as follows:

1. Take a type-I discrete cosine transform of the data, to
  determine the Fourier coefficients
  $$\hat{v}_j = \frac{\pi}{N}\left[v_0 + v_{N}\cos(\pi j) +
    2\sum^{N-1}_{k=1} v_{k}\cos\left(\frac{\pi
        jk}{N-1}\right)\right].$$
2. Let \(\hat{w}_{j} = -j\hat{v}_{j}\) for \(1\le j\le N-1\) and
  \(\hat{w}_N = 0\).
3. Take a type-I discrete sine transform of \(\hat{w}_j\) from \(j=1\)
  to \(j=N-1\), yielding
  $$w_j = \frac{1}{\pi}\sum^{N-1}_{k=1}\hat{w}_k\sin\left(\frac{\pi
      kj}{N}\right).$$
4. Compute
  $$ v_j' =  \begin{cases}
    \frac{1}{2\pi}\left[\frac{N^2}{2}\hat{v}_N + \sum_{k=1}^{N-1}
      k^2\hat{v}_k \right], & j = 0 \\
    \frac{-w_j}{\sqrt{1-x_j^2}}, & 1 \le j \le N - 1 \\
    \frac{1}{2\pi}\left[\frac{{(-1)}^{N+1}N^2}{2}\hat{v}_N +
      \sum_{k=1}^{N-1} {(-1)}^{k+1}k^2\hat{v}_k \right], & j = N
  \end{cases}.
  $$

Discrete sine and cosine transforms are variations of the discrete
Fourier transform which take advantage of data being real and either
even or odd. The FFTW3 package \citep{Frigo2005} was used to compute
these. A more rigorous and detailed explanation of the above methods for
periodic and non-periodic functions is provided
by~\citet{Trefethen2000} in chapters 3 and 8, respectively.

If a domain other than \(-1 \le x \le 1\) is desired then the
Collocation points can be scaled and offset as necessary, giving a
coordinate system \(x^*_j\). The above differentiation algorithm is
applied unchanged, except that the result is scaled by twice the
inverse of the domain-width.

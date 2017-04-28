module specfun_mod
!! Author: William Cody and Laura Stoltz
!!
!! SPECFUN is a FORTRAN90 library which evaluates certain special
!! functions, by William Cody and Laura Stoltz.
!!
!! In particular, SPECFUN can evaluate the I, J, K and Y Bessel
!! functions, of orders 0, 1, or arbitrary positive integer order N,
!! or for any positive non-integer order (an unusual feature).
!!
!! Routines are also available for the Gamma function, the logarithm
!! of the Gamma function, the exponential integrals, the error
!! function, the Psi function, and Dawson's integral.
!!
!! The original, true, correct (FORTRAN77) version of SPECFUN is
!! available through NETLIB:
!! http://www.netlib.org/specfun/index.html".
!!
!! @NOTE The routines used to calculate the real exponential integral
!! have been slightly modified so that they have the `elemental`
!! attribute. To do this, I needed to enclose them in a module. I was
!! having some issues compiling the module (the compiler didn't seem
!! to see the `r8_gamma` function and thus didn't correctly apply
!! name-mangling). To avoid these, I just deleted everything that I
!! didn't need.

contains

elemental subroutine calcei ( arg, result, jint )

!*****************************************************************************80
!
!! CALCEI computes various exponential integrals.
!
!  Discussion:
!
!    This routine computes the exponential integrals Ei(x),
!    E1(x), and  exp(-x)*Ei(x) for real arguments x where
!
!           integral (from t=-oo to t=x) (exp(t)/t),  x > 0,
!    Ei(x) =
!          -integral (from t=-x to t=+oo) (exp(t)/t),  x < 0,
!
!    and where the first integral is a principal value integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody, Henry Thacher,
!    Rational Chebyshev Approximations for the Exponential
!    Integral E1(x),
!    Mathematics of Computation,
!    Volume 22, Number 103, July 1968, pages 641-649.
!
!    William Cody, Henry Thacher,
!    Chebyshev Approximations for the Exponential
!    Integral Ei(x),
!    Mathematics of Computation,
!    Volume 23, Number 106, April 1969, pages 289-303.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG, the argument.  The argument must not
!    be zero.  If JINT = 2, then the argument must be strictly positive.
!
!    Output, real ( kind = 8 ) RESULT, the value of the function,
!    which depends on the input value of JINT:
!    1, RESULT = EI ( ARG );
!    2, RESULT = EONE ( ARG );
!    3, RESULT = EXPEI ( ARG ).
!
!    Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
!    1, Ei(x);
!    2, -Ei(-x);
!    3, exp(-x)*Ei(x).
!
  implicit none

  real ( kind = 8 ), intent(in) :: arg
  real ( kind = 8 ) ei
  real ( kind = 8 ) frac
  integer ( kind = 4 ) i
  integer ( kind = 4 ), intent(in) :: jint
  real ( kind = 8 ) px(10)
  real ( kind = 8 ) qx(10)
  real ( kind = 8 ), intent(out) :: result
  real ( kind = 8 ) sump
  real ( kind = 8 ) sumq
  real ( kind = 8 ) t
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ) xmx0
  real ( kind = 8 ) y
  real ( kind = 8 ) ysq
!
!  Mathematical constants
!  EXP40 = exp(40)
!  X0 = zero of Ei
!  X01/X11 + X02 = zero of Ei to extra precision
!
  real ( kind = 8 ), parameter :: p037 = 0.037D+00
  real ( kind = 8 ), parameter :: three = 3.0D+00
  real ( kind = 8 ), parameter :: six = 6.0D+00
  real ( kind = 8 ), parameter :: two4 = 24.0D+00
  real ( kind = 8 ), parameter :: fourty = 40.0D+00
  real ( kind = 8 ), parameter :: exp40 = 2.3538526683701998541D+17
  real ( kind = 8 ), parameter :: x01 = 381.5D+00
  real ( kind = 8 ), parameter :: x11 = 1024.0D+00
  real ( kind = 8 ), parameter :: x02 = -5.1182968633365538008D-05
  real ( kind = 8 ), parameter :: x0 = 3.7250741078136663466D-01
!
!  Machine-dependent constants
!
  real ( kind = 8 ), parameter :: xinf = 1.79d+308
  real ( kind = 8 ), parameter :: xmax = 716.351d0
  real ( kind = 8 ), parameter :: xbig = 701.84d0
!
!  Coefficients  for -1.0 <= X < 0.0
!
  real ( kind = 8 ), dimension(7), parameter :: a = &
         [1.1669552669734461083368d2, 2.1500672908092918123209d3, &
          1.5924175980637303639884d4, 8.9904972007457256553251d4, &
          1.5026059476436982420737d5,-1.4815102102575750838086d5, &
          5.0196785185439843791020d0]

  real ( kind = 8 ), dimension(6), parameter :: b = &
         [4.0205465640027706061433d1, 7.5043163907103936624165d2, &
          8.1258035174768735759855d3, 5.2440529172056355429883d4, &
          1.8434070063353677359298d5, 2.5666493484897117319268d5]
!
!  Coefficients for -4.0 <= X < -1.0
!
  real ( kind = 8 ), dimension(9), parameter :: c = &
         [3.828573121022477169108d-1, 1.107326627786831743809d+1, &
          7.246689782858597021199d+1, 1.700632978311516129328d+2, &
          1.698106763764238382705d+2, 7.633628843705946890896d+1, &
          1.487967702840464066613d+1, 9.999989642347613068437d-1, &
          1.737331760720576030932d-8]
  real ( kind = 8 ), dimension(9), parameter :: d = &
         [8.258160008564488034698d-2, 4.344836335509282083360d+0, &
          4.662179610356861756812d+1, 1.775728186717289799677d+2, &
          2.953136335677908517423d+2, 2.342573504717625153053d+2, &
          9.021658450529372642314d+1, 1.587964570758947927903d+1, &
          1.000000000000000000000d+0]
!
!  Coefficients for X < -4.0
!
  real ( kind = 8 ), dimension(10), parameter :: e = &
         [1.3276881505637444622987d+2,3.5846198743996904308695d+4, &
          1.7283375773777593926828d+5,2.6181454937205639647381d+5, &
          1.7503273087497081314708d+5,5.9346841538837119172356d+4, &
          1.0816852399095915622498d+4,1.0611777263550331766871d03, &
          5.2199632588522572481039d+1,9.9999999999999999087819d-1]
  real ( kind = 8 ), dimension(10), parameter :: f = &
         [3.9147856245556345627078d+4,2.5989762083608489777411d+5, &
          5.5903756210022864003380d+5,5.4616842050691155735758d+5, &
          2.7858134710520842139357d+5,7.9231787945279043698718d+4, &
          1.2842808586627297365998d+4,1.1635769915320848035459d+3, &
          5.4199632588522559414924d+1,1.0d0]
!
!  Coefficients for rational approximation to ln(x/a), |1-x/a| < .1
!
  real ( kind = 8 ), dimension(4), parameter :: plg = &
           [-2.4562334077563243311d+01,2.3642701335621505212d+02, &
            -5.4989956895857911039d+02,3.5687548468071500413d+02]
  real ( kind = 8 ), dimension(4), parameter :: qlg = &
           [-3.5553900764052419184d+01,1.9400230218539473193d+02, &
            -3.3442903192607538956d+02,1.7843774234035750207d+02]
!
!  Coefficients for  0.0 < X < 6.0,
!  ratio of Chebyshev polynomials
!
  real ( kind = 8 ), dimension(10), parameter :: p = &
         [-1.2963702602474830028590d01,-1.2831220659262000678155d03, &
          -1.4287072500197005777376d04,-1.4299841572091610380064d06, &
          -3.1398660864247265862050d05,-3.5377809694431133484800d08, &
           3.1984354235237738511048d08,-2.5301823984599019348858d10, &
           1.2177698136199594677580d10,-2.0829040666802497120940d11]
  real ( kind = 8 ), dimension(10), parameter :: q = &
          [7.6886718750000000000000d01,-5.5648470543369082846819d03, &
           1.9418469440759880361415d05,-4.2648434812177161405483d06, &
           6.4698830956576428587653d07,-7.0108568774215954065376d08, &
           5.4229617984472955011862d09,-2.8986272696554495342658d10, &
           9.8900934262481749439886d10,-8.9673749185755048616855d10]
!
!  J-fraction coefficients for 6.0 <= X < 12.0
!
  real ( kind = 8 ), dimension(10), parameter :: r = &
        [-2.645677793077147237806d00,-2.378372882815725244124d00, &
         -2.421106956980653511550d01, 1.052976392459015155422d01, &
          1.945603779539281810439d01,-3.015761863840593359165d01, &
          1.120011024227297451523d01,-3.988850730390541057912d00, &
          9.565134591978630774217d00, 9.981193787537396413219d-1]
  real ( kind = 8 ), dimension(9), parameter :: s = &
         [1.598517957704779356479d-4, 4.644185932583286942650d00, &
          3.697412299772985940785d02,-8.791401054875438925029d00, &
          7.608194509086645763123d02, 2.852397548119248700147d01, &
          4.731097187816050252967d02,-2.369210235636181001661d02, &
          1.249884822712447891440d00]
!
!  J-fraction coefficients for 12.0 <= X < 24.0
!
  real ( kind = 8 ), dimension(10), parameter :: p1 = &
         [-1.647721172463463140042d00,-1.860092121726437582253d01, &
          -1.000641913989284829961d01,-2.105740799548040450394d01, &
          -9.134835699998742552432d-1,-3.323612579343962284333d01, &
           2.495487730402059440626d01, 2.652575818452799819855d01, &
          -1.845086232391278674524d00, 9.999933106160568739091d-1]
  real ( kind = 8 ), dimension(9), parameter :: q1 = &
          [9.792403599217290296840d01, 6.403800405352415551324d01, &
           5.994932325667407355255d01, 2.538819315630708031713d02, &
           4.429413178337928401161d01, 1.192832423968601006985d03, &
           1.991004470817742470726d02,-1.093556195391091143924d01, &
           1.001533852045342697818d00]
!
!  J-fraction coefficients for  24 <= X.
!
  real ( kind = 8 ), dimension(10), parameter :: p2 = &
          [1.75338801265465972390d02,-2.23127670777632409550d02, &
          -1.81949664929868906455d01,-2.79798528624305389340d01, &
          -7.63147701620253630855d00,-1.52856623636929636839d01, &
          -7.06810977895029358836d00,-5.00006640413131002475d00, &
          -3.00000000320981265753d00, 1.00000000000000485503d00]
  real ( kind = 8 ), dimension(9), parameter :: q2 = &
          [3.97845977167414720840d04, 3.97277109100414518365d00, &
           1.37790390235747998793d02, 1.17179220502086455287d02, &
           7.04831847180424675988d01,-1.20187763547154743238d01, &
          -7.99243595776339741065d00,-2.99999894040324959612d00, &
           1.99999999999048104167d00]

  x = arg

  if ( x == 0.0D+00 ) then

    ei = - xinf

    if ( jint == 2 ) then
      ei = -ei
    end if
!
!  Calculate EI for negative argument or for E1.
!
  else if ( x < 0.0D+00 .or. jint == 2 ) then

    y = abs ( x )

    if ( y <= 1.0D+00 ) then

      sump = a(7) * y + a(1)
      sumq = y + b(1)
      do i = 2, 6
        sump = sump * y + a(i)
        sumq = sumq * y + b(i)
      end do
      ei = log ( y ) - sump / sumq

      if ( jint == 3 ) then
        ei = ei * exp ( y )
      end if

    else if ( y <= 4.0D+00 ) then

      w = 1.0D+00 / y
      sump = c(1)
      sumq = d(1)
      do i = 2, 9
        sump = sump * w + c(i)
        sumq = sumq * w + d(i)
      end do

      ei = - sump / sumq

      if ( jint /= 3 ) then
        ei = ei * exp ( -y )
      end if

    else

      if ( xbig < y .and. jint < 3 ) then

        ei = 0.0D+00

      else

        w = 1.0D+00 / y
        sump = e(1)
        sumq = f(1)
        do i = 2, 10
          sump = sump * w + e(i)
          sumq = sumq * w + f(i)
        end do

        ei = -w * ( 1.0D+00 - w * sump / sumq )

        if ( jint /= 3 ) then
          ei = ei * exp ( - y )
        end if

      end if

    end if

    if ( jint == 2 ) then
      ei = -ei
    end if
!
!  To improve conditioning, rational approximations are expressed
!  in terms of Chebyshev polynomials for 0 <= X < 6, and in
!  continued fraction form for larger X.
!
  else if ( x < six ) then

    t = x + x
    t = t / three - 2.0D+00
    px(1) = 0.0D+00
    qx(1) = 0.0D+00
    px(2) = p(1)
    qx(2) = q(1)
    do i = 2, 9
      px(i+1) = t * px(i) - px(i-1) + p(i)
      qx(i+1) = t * qx(i) - qx(i-1) + q(i)
    end do
    sump = 0.5D+00 * t * px(10) - px(9) + p(10)
    sumq = 0.5D+00 * t * qx(10) - qx(9) + q(10)
    frac = sump / sumq
    xmx0 = ( x - x01 / x11 ) - x02

    if ( p037 <= abs ( xmx0 ) ) then

      ei = log ( x / x0 ) + xmx0 * frac

      if ( jint == 3 ) then
        ei = exp ( - x ) * ei
      end if
!
!  Special approximation to ln(X/X0) for X close to X0.
!
    else

      y = xmx0 / ( x + x0 )
      ysq = y * y
      sump = plg(1)
      sumq = ysq + qlg(1)
      do i = 2, 4
        sump = sump * ysq + plg(i)
        sumq = sumq * ysq + qlg(i)
      end do
      ei = ( sump / ( sumq * ( x + x0 ) ) + frac ) * xmx0

      if ( jint == 3 ) then
        ei = exp ( - x ) * ei
      end if

    end if

  else if ( x < 12.0D+00 ) then

    frac = 0.0D+00
    do i = 1, 9
      frac = s(i) / ( r(i) + x + frac )
    end do

    ei = ( r(10) + frac ) / x
    if ( jint /= 3 ) then
      ei = ei * exp ( x )
    end if

  else if ( x <= two4 ) then

    frac = 0.0D+00
    do i = 1, 9
      frac = q1(i) / ( p1(i) + x + frac )
    end do

    ei = ( p1(10) + frac ) / x

    if ( jint /= 3 ) then
      ei = ei * exp ( x )
    end if

  else

    if ( xmax <= x .and. jint < 3 ) then

      ei = xinf

    else

      y = 1.0D+00 / x
      frac = 0.0D+00
      do i = 1, 9
        frac = q2(i) / ( p2(i) + x + frac )
      end do
      frac = p2(10) + frac
      ei = y + y * y * frac

      if ( jint /= 3 ) then

        if ( x <= xmax - two4 ) then
          ei = ei * exp ( x )
!
!  Calculation reformulated to avoid premature overflow.
!
        else
          ei = ( ei * exp ( x - fourty ) ) * exp40
        end if

      end if
    end if
  end if

  result = ei

  return
end

elemental function ei ( x )

!*****************************************************************************80
!
!! EI evaluates the exponential integral Ei(X).
!
!  Discussion:
!
!    This routine computes approximate values for the
!    exponential integral Ei(x), where x is real.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) EI, the value of the function.
!
  implicit none

  real ( kind = 8 ) ei
  integer ( kind = 4 ) jint
  real ( kind = 8 ) result
  real ( kind = 8 ), intent(in) :: x

  jint = 1
  call calcei ( x, result, jint )
  ei = result

  return
end

end module specfun_mod

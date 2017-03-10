!
!  nitsol.f90
!  This file is part of ISOFT.
!  
!  Copyright 2016 Chris MacMackin <cmacmackin@gmail.com>
!  
!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!  
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
!  MA 02110-1301, USA.
!  

module nitsol_mod
  !* Author: Christopher MacMackin
  !  Date: July 2016
  !  License: GPLv3
  !
  ! Provides an explicit interface to the
  ! [NITSOL](http://users.wpi.edu/~walker/Papers/nitsol,SISC_19,1998,302-318.pdf)
  ! package. Variables held in common blocks which can be used to
  ! control NITSOL are also provided here. At some point I may produce
  ! a proper object oriented interface for it.
  !
  ! Also present is an interface to the GMRES solver provided by
  ! NITSOL. Both a direct interface and a wrapper are provided. The
  ! wrapper offers a somewhat more general and F90-ish presentation of
  ! the routine.
  !
  use iso_fortran_env, only: r8 => real64
  implicit none

  ! Common blocks used by the legacy NITSOL code. This allows other
  ! code to use them via this module.
  integer  :: iplvl, ipunit
  common /nitprint/ iplvl, ipunit

  integer  :: instep, newstep, krystat
  real(r8) :: avrate, fcurnrm
  common /nitinfo/ avrate, fcurnrm, instep, newstep, krystat

  real(r8) :: choice1_exp, choice2_exp, choice2_coef
  real(r8) :: eta_cutoff, etamax
  real(r8) :: thmin, thmax, etafixed
  common /nitparam/ choice1_exp, choice2_exp, choice2_coef, &
                    eta_cutoff, etamax, thmin, thmax, etafixed

  abstract interface
    subroutine f_intr(n, xcur, fcur, rpar, ipar, itrmf)
      !! Interface for a subroutine which evaluates the function
      !! the zero of which is sought.
      import :: r8
      integer, intent(in)                   :: n
        !! Dimension of the problem
      real(r8), dimension(n), intent(in)    :: xcur
        !! Array of length `n` containing the current \(x\) value
      real(r8), dimension(n), intent(out)   :: fcur
        !! Array of length `n` containing f(xcur) on output
      real(r8), dimension(*), intent(inout) :: rpar
        !! Parameter/work array
      integer, dimension(*), intent(inout)  :: ipar
        !! Parameter/work array
      integer, intent(out)                  :: itrmf
        !! Termination flag. 0 means normal termination, 1 means
        !! failure to produce f(xcur)
    end subroutine f_intr

    subroutine jacv_intr(n, xcur, fcur, ijob, v, z, rpar, ipar, itrmjv)
      !! Interface for a subroutine which optionally evaluates
      !! \(J\vec{v}\) or \(P^{-1}\vec{v}\), where \(J\) is the Jacobian
      !! of \(f\) and \(P\) is a right preconditioning operator.
      import :: r8
      integer, intent(in)                   :: n
        !! Dimension of the problem
      real(r8), dimension(n), intent(in)    :: xcur
        !! Array of length `n` containing the current \(x\) value
      real(r8), dimension(n), intent(in)    :: fcur
        !! Array of length `n` containing the current \(f(x)\) value
      integer, intent(in)                   :: ijob
        !! Integer flag indicating which product is desired. 0
        !! indicates \(z = J\vec{v}\). 1 indicates \(z = P^{-1}\vec{v}\).
      real(r8), dimension(n), intent(in)    :: v
        !! An array of length `n` to be multiplied
      real(r8), dimension(n), intent(out)   :: z
        !! An array of length n containing the desired product on
        !! output.
      real(r8), dimension(*), intent(inout) :: rpar
        !! Parameter/work array 
      integer, dimension(*), intent(inout)  :: ipar
        !! Parameter/work array
      integer, intent(out)                  :: itrmjv
        !! Termination flag. 0 indcates normal termination, 1
        !! indicatesfailure to prodce \(J\vec{v}\), and 2 indicates
        !! failure to produce \(P^{-1}\vec{v}\)
    end subroutine jacv_intr
      
    function dinpr_intr(n, x, sx, y, sy)
      !! Interface for function which calculates vector inner products.
      !! This has the same interace as the BLAS routine
      !! [ddot](http://www.netlib.org/lapack/explore-html/de/da4/group__double__blas__level1_ga75066c4825cb6ff1c8ec4403ef8c843a.html).
      import :: r8
      integer, intent(in)                :: n
        !! The length of the vectors
      real(r8), dimension(*), intent(in) :: x
        !! The first input vector
      integer, intent(in)                :: sx
        !! The stride in memory between successive elements of `x`
      real(r8), dimension(*), intent(in) :: y
        !! The second input vector
      integer, intent(in)                :: sy
        !! The stride in memory between successive elements of `y`
      real(r8)                           :: dinpr_intr
        !! Inner product of `x` and `y`
    end function dinpr_intr

    function dnorm_intr(n, x, sx)
      !! Interface for function which calculates vector norms. This
      !! has the same interface as the BLAS routine dnrm2.
      import :: r8
      integer, intent(in)                :: n
        !! The length of the array
      real(r8), dimension(*), intent(in) :: x
        !! The input vector
      integer, intent(in)                :: sx
        !! The stride in memory between consecutive elements of `x`
      real(r8) :: dnorm_intr
        !! The vector norm of `x`
    end function dnorm_intr

    function mat_mult(v, xcur, rhs, rpar, ipar, success)
      !! Interface for operations representing the multiplication of a
      !! vector by a matrix, such as that for a linear operator or a
      !! preconditioner.
      import :: r8
      real(r8), dimension(:), intent(in)    :: v
        !! The vector to be multiplied
      real(r8), dimension(:), intent(in)    :: xcur
        !! Array containing the current estimate of the independent
        !! variables in the linear system. This may not be needed, but
        !! is provided just in case.
      real(r8), dimension(:), intent(in)    :: rhs
        !! Array containing the right hand side of the linear
        !! system. This may not be needed, but is provided just in
        !! case.
      real(r8), dimension(*), intent(inout) :: rpar
        !! Parameter/work array 
      integer, dimension(*), intent(inout)  :: ipar
        !! Parameter/work array
      logical, intent(out)                  :: success
        !! Indicates whether operation was completed succesfully
      real(r8), dimension(size(xcur))       :: mat_mult
        !! Result of the operation
    end function mat_mult
  end interface


  interface
    subroutine nitsol(n, x, f, jacv, ftol, stptol, input, info, rwork, &
                      rpar, ipar, iterm, dinpr, dnorm)
      !* Author: Chris MacMackin
      !  Date: July 2016
      !
      ! An explicit interface to the 
      ! [nitsol](http://users.wpi.edu/~walker/Papers/nitsol,SISC_19,1998,302-318.pdf)
      ! Newton iterative nonlinear solver.
      !
      !####Input parameters
      !
      ! The `input` array argument allows the user to specify various
      ! options. It should be declared an integer vector of length 11
      ! [sic.] in the calling program. To specify an option, set the
      ! appropriate input component to the desired value according to
      ! the specifications below.
      !
      ! **Usage Note:** Setting a particular input component to zero gives the 
      ! default option for that component in all cases. 
      !
      ! The first five input components are things that every user might wish 
      ! to modify; the remainder will usually be of interest only to more 
      ! experienced users. 
      !
      ! Optional every-user input:
      !
      !    input(1) = nnimax = maximum number of nonlinear iterations (default 200).
      ! 
      !    input(2) = ijacv = flag for determining the method of J*v evaluation:
      !                 0 => finite-difference evaluation (default) 
      !                 1 => analytic evaluation
      !
      !    input(3) = ikrysl = flag for determining the Krylov solver: 
      !                 0 => GMRES (default)
      !                 1 => BiCGSTAB
      !                 2 => TFQMR
      !
      !               For brief descriptions of the solvers plus references, 
      !               see the subroutines nitgm, nitstb, and nittfq. 
      !
      !    input(4) = kdmax = maximum Krylov subspace dimension when GMRES is used 
      !               (default 20). 
      !
      !    input(5) = irpre = flag for right preconditioning: 
      !                 0 => no right preconditioning
      !                 1 => right preconditioning
      !
      ! Optional experienced user input:
      !
      !    input(6) = iksmax = maximum allowable number of iterations per call 
      !               to the Krylov solver routine (default 1000). 
      !
      !    input(7) = iresup = residual update flag when GMRES is used; on 
      !               restarts, the residual is updated as follows: 
      !                 0 => linear combination (default) 
      !                 1 => direct evaluation
      !               The first is cheap (one n-vector saxpy) but may lose 
      !               accuracy with extreme residual reduction; the second 
      !               retains accuracy better but costs one J*v product per 
      !               restart. 
      !
      !    input(8) = ifdord = order of the finite-difference formula (sometimes) 
      !               used when input(2) = ijacv = 0. When input(2) = ijacv = 0, 
      !               this must be 0, 1, 2, or 4 on input; otherwise, it is 
      !               irrelevant. With input(2) = ijacv = 0, the precise 
      !               meaning is as follows: 
      !
      !               If GMRES is used, then ifdord matters only if input(7) = 
      !               iresup = 1, in which case it determines the order of 
      !               the finite-difference formula used in evaluating the 
      !               initial residual at each GMRES restart (default 2); if 
      !               ifdord = 0 on input, then it is set to 2 below. NOTE: This 
      !               only affects initial residuals at restarts; first-order 
      !               differences are always used within each GMRES cycle. Using 
      !               higher-order differences at restarts only should give 
      !               the same accuracy as if higher-order differences were 
      !               used throughout; see K. Turner and H. F. Walker, "Efficient 
      !               high accuracy solutions with GMRES(m)," SIAM J. Sci. 
      !               Stat. Comput., 13 (1992), pp. 815--825. 
      !               
      !               If BiCGSTAB or TFQMR is used, then ifdord determines the 
      !               order of the finite-difference formula used at each 
      !               iteration (default 1); if ifdord = 0 on input, then it 
      !               is set to 1 below. 
      !
      !    input(9) = ibtmax = maximum allowable number of backtracks (step 
      !               reductions) per call to nitbt (default 10). 
      !
      !               USAGE NOTE: Backtracking can be turned off by setting 
      !	              ibtmax = -1. Other negative values of ibtmax are not 
      !               valid. 
      !
      !    input(10) = ieta = flag determining the forcing term eta as follows: 
      !                 0 => abs( ||fcur|| - ||fprev+Jprev*sprev|| )/||fprev||
      !                      (default) 
      !                 1 => (||fcur||/||fprev||)**2 
      !                 2 => gamma*(||fcur||/||fprev||)**alpha 
      !                      for user-supplied gamma in (0,1] and alpha in (1,2] 
      !                 3 => fixed (constant) eta in (0,1), either 0.1 (default) 
      !	                     or specified by the user (see USAGE NOTE below) 
      !               Here, fcur = current f, fprev = previous f, etc. The Krylov 
      !               iterations are terminated when an iterate s satisfies 
      !               an inexact Newton condition ||F + J*s|| .le. eta*||F||.
      !
      !               USAGE NOTE: If input(10) = ieta = 2, then alpha and gamma 
      !               must be set in common block nitparam.h as described below. 
      !	              If input(10) = ieta = 3, then the desired constant eta may 
      !	              be similarly set in nitparam.h if a value other than the 
      !	              default of 0.1 is desired. 
      !               
      !               The first three expressions above are from S. C. Eisenstat 
      !               and H. F. Walker, "Choosing the forcing terms in an inexact 
      !               Newton method", SIAM J. Scientific Computing, 17 (1996), 
      !               pp. 16--32. (They may be modified according to certain 
      !               safeguards in subroutine nitdrv.) The first gives convergence 
      !               that is q-superlinear and of r-order (1+sqrt(5))/2. The 
      !               second gives convergence that is r-quadratic and of q-order 
      !               p for every p in [1,2). The third gives convergence that is 
      !               of q-order alpha when gamma < 1 and, when gamma = 1, of 
      !               r-order alpha and q-order p for every p in [1,alpha). The 
      !               fourth gives q-linear convergence with asymptotic rate 
      !               constant eta in a certain norm; see R. S. Dembo, S. C. 
      !	              Eisenstat, and T. Steihaug, "Inexact Newton methods", 
      !               SIAM J. Numer. Anal., 18 (1982), pp. 400-408. 
      !
      !               Of these four choices, the 1st is usually satisfactory, 
      !               the 2nd or 3rd is sometimes preferred, and the 4th may be 
      !               useful in some situations, e.g., it may be desirable to 
      !               choose a fairly large fixed eta in (0,1), such as eta = .1, 
      !               when numerical inaccuracy prevents the Krylov solver 
      !               from obtaining much residual reduction. 
      !
      !
      !####Output parameters
      !
      ! On output, the components of the `info` argument are as follows:
      !
      !     info(1)   = nfe (number of function evaluations)
      !     info(2)   = njve (number of J*v evaluations)
      !     info(3)   = nrpre (number of P(inverse)*v evaluations)
      !     info(4)   = nli (number of linear iterations)
      !     info(5)   = nni (number of nonlinear iterations)
      !     info(6)   = nbt (number of backtracks)
      import :: r8
      implicit none
      integer, intent(in)                   :: n
        !! Dimension of the problem
      real(r8), dimension(n), intent(inout) :: x
        !! Vector of length n, initial guess on input and final
        !! approximate solution on output
      procedure(f_intr)                     :: f
        !! User-supplied subroutine for evaluating the function
        !! the zero of which is sought.
      procedure(jacv_intr)                  :: jacv
        !! User-supplied subroutine for optionally evaluating
        !! \(J\vec{v}\) or \(P^{-1}\vec{v}\), where \(J\) is the Jacobian
        !! of \(f\) and \(P\) is a right preconditioning operator. If
        !! neither analytic \(J\vec{v}\) evaluations nor right
        !! preconditioning is used, this can be a dummy subroutine;
        !! if right preconditioning is used but not analytic
        !! \(J\vec{v}\) evaluations, this need only evaluate
        !! \(P^{-1}\vec{v}\).
      real(r8), intent(in)                  :: ftol
        !! Stopping tolerance of the f-norm
      real(r8), intent(in)                  :: stptol
        !! Stopping tolerance of the step-length
      integer, dimension(10), intent(in)    :: input
        !! Array containing various user-specified inputs; see above
      integer, dimension(6), intent(out)    :: info
        !! Array containing various outputs; see above
      real(r8), dimension(*), intent(inout) :: rwork
        !! Work array with length depending on the solver used as follows:
        !!
        !!GMRES
        !!:    \(n\times(\text{kdmax}+5)+\text{kdmax}\times(\text{kdmax}+3)\),
        !!     where kdmax is the maximum Krylove subspace dimension, either
        !!     the default value of 20 or another value specified by the user
        !!
        !!BiCGSTAB
        !!:    \(11n\)
        !!
        !!TFQMR
        !!:    \(14n\)
      real(r8), dimension(*), intent(inout) :: rpar
        !! Parameter/work array passed to the `f` and `jacv` routines
      integer, dimension(*), intent(inout)  :: ipar
        !! Parameter/work array passed to the `f` and `jacv` routines
      integer, intent(out)                  :: iterm
        !! Termination flag. Values have the following meanings:
        !!
        !!-k
        !!:    illegal value in `input(k)`
        !!
        !!0
        !!:    normal termination: \(||F|| < \text{ftol}\) or \(||\text{step}||
        !!     < \text{stptol}\)
        !!
        !!1
        !!:    `nnimax` nonlinar iterations reached without success
        !!
        !!2
        !!:    failure to evaluate \(F\)
        !!
        !!3
        !!:    in `nitjv`, \(J\vec{v}\)  failure
        !!
        !!4
        !!:    in `nitjv`, \(P^{-1}\vec{v}\) failure
        !!
        !!5
        !!:    in `nitdrv`, insufficient initial model norm reduction for
        !!     adequate progress. **Note:** This can occur for several
        !!     reasons; examine `itrmks` on return from the Krylov solver
        !!     for further information. (This will be printed out if
        !!     \(\text{iplvl}\ge 3\); see the discussion of optional common
        !!     blocks below.)
        !!
        !!6
        !!:    in `nitbt`, failure to reach an acceptable step through
        !!     backtracking
      procedure(dinpr_intr)                 :: dinpr
        !! User-supplied function for calculating vector inner products.
        !! This has the same interace as the BLAS routine
        !! [ddot](http://www.netlib.org/lapack/explore-html/de/da4/group__double__blas__level1_ga75066c4825cb6ff1c8ec4403ef8c843a.html).
        !! If the Euclidean inner product is desired then user can link
        !! to a local BLAS library and provide the name `ddot` to `nitsol`.
        !! `dinpr` must be declared as an external function that returns
        !! a double precision in the calling program.
      procedure(dnorm_intr)                 :: dnorm
        !! User-supplied function for calculating vector norms. This
        !! has the same interface as the BLAS routine dnrm2; if the
        !! Euclidean norm is desired the user can link to a local
        !! BLAS library and provide the name dnrm2 to nitsol.  dnorm
        !! must be declared as an external function that returns a
        !! double precision value in the calling program.
    end subroutine nitsol

    subroutine nitgm2(n, xcur, fcur, step, eta, f, jacv, rpar, ipar,    &
                      ijacv, irpre, iksmax, iresup, ifdord, nfe, njve,  &
                      nrpre, nli, kdmax, kdmaxp1, vv, rr, svbig, svsml, &
                      w, rwork, rsnrm, dinpr, dnorm, itrmks)
      !* Author: Chris MacMackin
      !  Date: March 2017
      !
      ! An interface to my modified versino of the
      ! [nitsol](http://users.wpi.edu/~walker/Papers/nitsol,SISC_19,1998,302-318.pdf)
      ! implementation of the generalised minimal residual method
      ! ([GMRES](https://en.wikipedia.org/wiki/Generalized_minimal_residual_method))
      ! for iteratively solving linear systems. It has been modified
      ! so that the user provides a non-zero initial guess of the
      ! solution.
      !
      import :: r8
      implicit none
      integer, intent(in)                            :: n
        !! Dimension of the problem
      real(r8), dimension(n), intent(in)             :: xcur
        !! Array of length `n` containing the current \(x\) value
      real(r8), dimension(n), intent(in)             :: fcur
        !! Array of length `n` containing current approximate solution
      real(r8), dimension(n), intent(inout)          :: step
        !! Vector of of length `n` containing trial step
      real(r8), intent(in)                           :: eta
        !! Relative residual reduction factor
      procedure(f_intr)                              :: f
        !! User-supplied subroutine for evaluating the function
        !! the zero of which is sought.
      procedure(jacv_intr)                           :: jacv
        !! User-supplied subroutine for optionally evaluating
        !! \(J\vec{v}\) or \(P^{-1}\vec{v}\), where \(J\) is the Jacobian
        !! of \(f\) and \(P\) is a right preconditioning operator. If
        !! neither analytic \(J\vec{v}\) evaluations nor right
        !! preconditioning is used, this can be a dummy subroutine;
        !! if right preconditioning is used but not analytic
        !! \(J\vec{v}\) evaluations, this need only evaluate
        !! \(P^{-1}\vec{v}\).
      real(r8), dimension(*), intent(inout)          :: rpar
        !! Parameter/work array passed to the `f` and `jacv` routines
      integer, dimension(*), intent(inout)           :: ipar
        !! Parameter/work array passed to the `f` and `jacv` routines
      integer, intent(in)                            :: ijacv
        !! Flag for determining method of \(J\vec{v}\) evaluation. 0
        !! (default) indicates finite-difference evaluation, while 1
        !! indicates analytic evaluation.
      integer, intent(in)                            :: irpre
        !! Flag for right preconditioning. 0 indicates no
        !! preconditioning, while 1 inidcates right preconditioning.
      integer, intent(in)                            :: iksmax
        !! Maximum allowable number of GMRES iterations
      integer, intent(in)                            :: iresup
        !! Residual update flag. On GMRES restarts, the residual can
        !! be updated using a linear combination (`iresup == 0`) or by
        !! direct evaluation (`iresup == 1`). The first is cheap (one
        !! n-vector saxpy) but may lose accuracy with extreme residual
        !! reduction; the second retains accuracy better but costs one
        !! \(J\vec{v}\) product.
      integer, intent(in)                            :: ifdord
        !! Order of the finite-difference formula (sometimes) used on
        !! GMRES restarts when \(J\vec{v}\) products are evaluated
        !! using finite- differences. When ijacv = 0 on input to
        !! nitsol, ifdord is set to 1, 2, or 4 in nitsol; otherwise,
        !! it is irrelevant. When ijacv = 0 on input to this
        !! subroutine, the precise meaning is as follows:
        !!                                                             
        !! With GMRES, ifdord matters only if iresup = 1, in which case  
        !! it determines the order of the finite-difference formula used 
        !! in evaluating the initial residual at each GMRES restart      
        !! (default 2). If iresup = 1 and ijacv = 0 on input to this     
        !! subroutine, then ijacv is temporarily reset to -1 at each     
        !! restart below to force a finite-difference evaluation of order
        !! ifdord. NOTE: This only affects initial residuals at restarts;
        !! first-order differences are always used within each GMRES     
        !! cycle. Using higher-order differences at restarts only should 
        !! give the same accuracy as if higher-order differences were    
        !! used throughout; see K. Turner and H. F. Walker, "Efficient   
        !! high accuracy solutions with GMRES(m)," SIAM J. Sci. Stat.    
        !! Comput., 13 (1992), pp. 815--825.
      integer, intent(inout)                         :: nfe
        !! Number of function evaluations
      integer, intent(inout)                         :: njve
        !! Number of \(J\vec{v}\) evaluations
      integer, intent(inout)                         :: nrpre
        !! Number of \(P^{-1}\vec{v}\) evaluations
      integer, intent(inout)                         :: nli
        !! Number of linear iterations
      integer, intent(in)                            :: kdmax
        !! Maximum Krylov subspace dimension; default 10.
      integer, intent(in)                            :: kdmaxp1
        !! kdmax + 1
      real(r8), dimension(n, kdmaxp1), intent(out)   :: vv
        !! Matrix for storage of Krylov basis in GMRES; on return, the
        !! residual vector is contained in the first column.
      real(r8), dimension(kdmax, kdmax), intent(out) :: rr
        !! Matrix for storage of triangular matrix in GMRES.
      real(r8), dimension(kdmax), intent(out)        :: svbig
        !! Vector for storage of estimate of singular vector of `rr`
        !! with largest singular value.
      real(r8), dimension(kdmax), intent(out)        :: svsml
        !! Vector for storage of estimate of singular vector of `rr`
        !! with smallest singular value.
      real(r8), dimension(kdmax), intent(out)        :: w
        !! Vector containing right-hand side of triangular system and
        !! least-squares residual norm in GMRES.
      real(r8), dimension(n), intent(out)            :: rwork
        !! Work array
      real(r8), intent(out)                          :: rsnrm
        !! GMRES residual norm on return
      procedure(dinpr_intr)                          :: dinpr
        !! Inner-product routine, either user-supplied or BLAS `ddot`.
      procedure(dnorm_intr)                          :: dnorm
        !! Norm routine, either user supplied or BLAS dnrm2.
      integer, intent(out)                           :: itrmks
        !! Termination flag. Values have the following meanings:
        !!
        !!0
        !!:    normal termination: acceptable step found
        !!
        !!1
        !!:    \(J\vec{v}\)  failure in `nitjv`
        !!
        !!2
        !!:    \(P^{-1}\vec{v}\) failure in `nitjv`
        !!
        !!3
        !!:    Acceptable step not found in `iksmax` GMRES iterations
        !!
        !!4
        !!:    Insignificant residual norm reduction of a cycle of `kdmax` 
        !!     steps (stagnation) before an acceptable step has been found.
        !!
        !!5
        !!:    Dangerous ill-conditioning detected before an acceptable 
        !!     step has been found.
    end subroutine nitgm2
  end interface


  interface
    function ddot(n, x, sx, y, sy)
      !* Author: Chris MacMackin
      !  Date: November 2016
      !
      ! An interface to the BLAS routine for calculating Euclidean
      ! inner product. This can be passed to [[nitsol]] for the
      ! argument `dinpr`.
      !
      import :: r8
      implicit none
      integer, intent(in)                :: n
      real(r8), dimension(*), intent(in) :: x
      integer, intent(in)                :: sx
      real(r8), dimension(*), intent(in) :: y
      integer, intent(in)                :: sy
      real(r8)                           :: ddot
    end function ddot

    function dnrm2(n, x, sx)
      !* Author: Chris MacMackin
      !  Date: November 2016
      !
      ! An interface to the BLAS routine for calculating Euclidean
      ! norm. This can be passed to [[nitsol]] for the argument
      ! `dnorm`.
      !
      import :: r8
      implicit none
      integer, intent(in)                :: n
      real(r8), dimension(*), intent(in) :: x
      integer, intent(in)                :: sx
      real(r8)                           :: dnrm2
    end function dnrm2
  end interface 

contains

  subroutine dummy_jacv(n, xcur, fcur, ijob, v, z, rpar, ipar, itrmjv)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! A dummy subroutine which does not apply a preconditioner or
    ! calculate an analytic Jacobian. This can be passed to [[nitsol]]
    ! for the argument `jacv`.
    ! 
    integer, intent(in)                   :: n
      ! Dimension of the problem
    real(r8), dimension(n), intent(in)    :: xcur
      ! Array of length `n` containing the current \(x\) value
    real(r8), dimension(n), intent(in)    :: fcur
      ! Array of length `n` containing the current \(f(x)\) value
    integer, intent(in)                   :: ijob
      ! Integer flag indicating which product is desired. 0
      ! indicates \(z = J\vec{v}\). 1 indicates \(z = P^{-1}\vec{v}\).
    real(r8), dimension(n), intent(in)    :: v
      ! An array of length `n` to be multiplied
    real(r8), dimension(n), intent(out)   :: z
      ! An array of length n containing the desired product on
      ! output.
    real(r8), dimension(*), intent(inout) :: rpar
      ! Parameter/work array 
    integer, dimension(*), intent(inout)  :: ipar
      ! Parameter/work array
    integer, intent(out)                  :: itrmjv
      ! Termination flag. 0 indcates normal termination, 1
      ! indicatesfailure to prodce \(J\vec{v}\), and 2 indicates
      ! failure to produce \(P^{-1}\vec{v}\)
  end subroutine dummy_jacv


  subroutine dummy_f(n, xcur, fcur, rpar, ipar, itrmf)
    !* Author: Chris MacMackin
    !  Date: March 2017
    !
    ! A dummy subroutine which does not calculate the function.
    ! 
    integer, intent(in)                   :: n
      !! Dimension of the problem
    real(r8), dimension(n), intent(in)    :: xcur
      !! Array of length `n` containing the current \(x\) value
    real(r8), dimension(n), intent(out)   :: fcur
      !! Array of length `n` containing f(xcur) on output
    real(r8), dimension(*), intent(inout) :: rpar
      !! Parameter/work array
    integer, dimension(*), intent(inout)  :: ipar
      !! Parameter/work array
    integer, intent(out)                  :: itrmf
      !! Termination flag. 0 means normal termination, 1 means
      !! failure to produce f(xcur)
  end subroutine dummy_f


  subroutine gmres_solve(solution, lhs, rhs, resid_norm, flag, tol, precond, &
                         rpar, ipar, resid_update, iter_max, krylov_dim,     &
                         inner_prod, norm)
    !* Author: Chris MacMackin
    !  Date: March 2017
    !
    ! A wraper for the
    ! [nitsol](http://users.wpi.edu/~walker/Papers/nitsol,SISC_19,1998,302-318.pdf)
    ! implementation of the generalised minimal residual method
    ! ([GMRES](https://en.wikipedia.org/wiki/Generalized_minimal_residual_method))
    ! for iteratively solving linear systems. This provides a more
    ! general interface not specifically intended for use with Newton
    ! iteration. It also uses Fortran 90 features to provide a more
    ! convenient call signature.
    !
    real(r8), dimension(:), intent(inout)           :: solution
      !! On input, an initial guess of the solution to the linear
      !! system. On output, the iteratively determined solution.
    procedure(mat_mult)                             :: lhs
      !! The linear operator on the left hand side of the linear
      !! system.
    real(r8), dimension(:), intent(in)              :: rhs
      !! The right hand side of the linear system being solved
    real(r8), intent(out)                           :: resid_norm
      !! GMRES residual norm on return
    integer, intent(out)                            :: flag
      !! Termination flag. Values have the following meanings:
      !!
      !!0
      !!:    normal termination: acceptable solution found
      !!
      !!1
      !!:    \(J\vec{v}\)  failure
      !!
      !!2
      !!:    \(P^{-1}\vec{v}\) failure
      !!
      !!3
      !!:    Acceptable solution not found in `iter_max` GMRES iterations
      !!
      !!4
      !!:    Insignificant residual norm reduction of a cycle of `kdmax` 
      !!     steps (stagnation) before an acceptable step has been found.
      !!
      !!5
      !!:    Dangerous ill-conditioning detected before an acceptable 
      !!     step has been found.
    real(r8), intent(in), optional                  :: tol
      !! The tolerance for the solution. Default is `size(solution) * 1e-8`.
    procedure(mat_mult), optional                   :: precond
      !! A right-preconditioner which may be used to improve
      !! convergence of the solution.
    real(r8), dimension(*), intent(inout), optional :: rpar
      !! Parameter/work array passed to the `lhs` and `precond` routines.
    integer, dimension(*), intent(inout), optional  :: ipar
      !! Parameter/work array passed to the `lhs` and `precond` routines
    integer, intent(in), optional                   :: resid_update
      !! Residual update flag. On GMRES restarts, the residual can
      !! be updated using a linear combination (`iresup == 0`) or by
      !! direct evaluation (`iresup == 1`). The first is cheap (one
      !! n-vector saxpy) but may lose accuracy with extreme residual
      !! reduction; the second retains accuracy better but costs one
      !! \(J\vec{v}\) product. Default is 0.
    integer, intent(in), optional                   :: iter_max
      !! Maximum allowable number of GMRES iterations. Default is
      !! 1000.
    integer, intent(in), optional                   :: krylov_dim
      !! Maximum Krylov subspace dimension; default 10.    
    procedure(dinpr_intr), optional                 :: inner_prod
      !! Inner-product routine, either user-supplied or BLAS `ddot`.
    procedure(dnorm_intr), optional                 :: norm
      !! Norm routine, either user supplied or BLAS dnrm2.
    
    integer  :: npoints, preflag, resup, itmax, kdim, nfe, njve, nrpre, nli
    real(r8) :: eta
    procedure(dinpr_intr), pointer :: dinpr
    procedure(dnorm_intr), pointer :: dnorm
    real(r8), dimension(:), allocatable, save :: xcur, svbig, svsml, w, rwork
    real(r8), dimension(:,:), allocatable, save :: vv, rr

    npoints = size(solution)
    if (present(precond)) then
       preflag = 1
    else
       preflag = 0
    end if
    if (present(tol)) then
      eta = tol
    else
      eta = 1.e-8_r8 * npoints
    end if
    if (present(resid_update)) then
      resup = resid_update
    else
      resup = 0
    end if
    if (present(iter_max)) then
      itmax = iter_max
    else
      itmax = 1000
    end if
    if (present(krylov_dim)) then
      kdim = krylov_dim
    else
      kdim = 10
    end if
    if (present(inner_prod)) then
      dinpr => inner_prod
    else
      dinpr => ddot
    end if
    if (present(norm)) then
      dnorm => norm
    else
      dnorm => dnrm2
    end if

    xcur = solution

    if (allocated(vv)) then
      if (any(shape(vv) /= [npoints, kdim+1])) then
        deallocate(vv)
        allocate(vv(npoints, kdim+1))
      end if
      if (size(svbig) /= kdim) then
        deallocate(rr)
        deallocate(svbig)
        deallocate(svsml)
        deallocate(w)
        allocate(rr(kdim, kdim))
        allocate(svbig(kdim))
        allocate(svsml(kdim))
        allocate(w(kdim))
      end if
      if (size(rwork) /= npoints) then
        deallocate(rwork)
        allocate(rwork(npoints))
      end if
    else
      allocate(vv(npoints, kdim+1))
      allocate(rr(kdim, kdim))
      allocate(svbig(kdim))
      allocate(svsml(kdim))
      allocate(w(kdim))
      allocate(rwork(npoints))
    end if

    call nitgm2(npoints, xcur, -rhs, solution, eta, dummy_f, jacv,  &
                rpar, ipar, 1, preflag, itmax, resup, 1, nfe, njve, &
                nrpre, nli, kdim, kdim+1, vv, rr, svbig, svsml, w,  &
                rwork, resid_norm, dinpr, dnorm, flag)

  contains
    
    subroutine jacv(n, xcur, fcur, ijob, v, z, rpar, ipar, itrmjv)
      !! A wrapper on the user-provided routines for the linear
      !! operator and preconditioner, to put them in the form NITSOL
      !! expects.
      !! 
      integer, intent(in)                   :: n
        ! Dimension of the problem
      real(r8), dimension(n), intent(in)    :: xcur
        ! Array of length `n` containing the current \(x\) value
      real(r8), dimension(n), intent(in)    :: fcur
        ! Array of length `n` containing the current \(f(x)\) value
      integer, intent(in)                   :: ijob
        ! Integer flag indicating which product is desired. 0
        ! indicates \(z = J\vec{v}\). 1 indicates \(z = P^{-1}\vec{v}\).
      real(r8), dimension(n), intent(in)    :: v
        ! An array of length `n` to be multiplied
      real(r8), dimension(n), intent(out)   :: z
        ! An array of length n containing the desired product on
        ! output.
      real(r8), dimension(*), intent(inout) :: rpar
        ! Parameter/work array 
      integer, dimension(*), intent(inout)  :: ipar
        ! Parameter/work array
      integer, intent(out)                  :: itrmjv
        ! Termination flag. 0 indcates normal termination, 1
        ! indicatesfailure to prodce \(J\vec{v}\), and 2 indicates
        ! failure to produce \(P^{-1}\vec{v}\)
      logical :: success
      itrmjv = 0
      if (ijob == 0) then
        z = lhs(v, xcur, fcur, rpar, ipar, success)
        if (.not. success) itrmjv = 1
      else if (ijob == 1) then
        z = precond(v, xcur, fcur, rpar, ipar, success)
        if (.not. success) itrmjv = 2
      else
        error stop ('`ijob` not equal to 0 or 1.')
      end if
    end subroutine jacv

  end subroutine  gmres_solve

!  function scaled_norm(n, x, sx)
!    !* Author: Chris MacMackin
!    !  Date: December 2016
!    !
!    ! A wraper to the BLAS routine for calculating Euclidean norm
!    ! which then scales the norm by the size of the vector. This can
!    ! be passed to [[nitsol]] for the argument `dnorm`.
!    !
!    implicit none
!    integer, intent(in)                :: n
!    real(r8), dimension(*), intent(in) :: x
!    integer, intent(in)                :: sx
!    real(r8)                           :: scaled_norm
!    scaled_norm = dnrm2(n,x,sx)/real(n,r8)
!  end function scaled_norm
!
end module nitsol_mod

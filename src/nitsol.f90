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
  ! package. At some point I may produce a proper object oriented interface
  ! for it.
  !
  use iso_fortran_env, only: r8 => real64

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
      import :: r8
      integer, intent(in)                   :: n
        !! Dimension of the problem
      real(r8), dimension(n), intent(inout) :: x
        !! Vector of length n, initial guess on input and final
        !! approximaet solution on output

      interface
         subroutine f(n, xcur, fcur, rpar, ipar, itrmf)
           !! User-supplied subroutine for evaluating the function
           !! the zero of which is sought.
           import :: r8
           integer                                  :: n
             !! Dimension of the problem
           real(r8), dimension(*), intent(in)       :: xcur
             !! Array of length `n` containing the current $x$ value
           real(r8), dimension(*), intent(out)      :: fcur
             !! Array of length `n` containing f(xcur) on output
           real(r8), dimension(*), intent(inout)    :: rpar
             !! Parameter/work array
           integer(r8), dimension(*), intent(inout) :: ipar
             !! Parameter/work array
           integer(r8), intent(out)                 :: itrmf
             !! Termination flag. 0 means normal termination, 1 means
             !! failure to produce f(xcur)
         end subroutine f

         subroutine jacv(n, xcur, fcur, ijob, v, z, rpar, ipar, itrmjv)
           !! User-supplied subroutine for optionally evaluating
           !! $J\vec{v}$ or $P^{-1}\vec{v}$, where $J$ is the Jacobian
           !! of $f$ and $P$ is a right preconditioning operator. If
           !! neither analytic $J\vec{v}$ evaluations nor right
           !! preconditioning is used, this can be a dummy subroutine;
           !! if right preconditioning is used but not analytic
           !! $J\vec{v}$ evaluations, this need only evaluate
           !! $P^{-1}\vec{v}$.
           import :: r8
           integer, intent(in)                   :: n
             !! Dimension of the problem
           real, dimension(*), intent(in)        :: xcur
             !! Array of lenght `n` containing the current $x$ value
           real, dimension(*), intent(in)        :: fcur
             !! Array of lenght `n` containing the current $f(x)$ value
           integer, intent(in)                   :: ijob
             !! Integer flat indicating which product is desired. 0
             !! indicates $z = J\vec{v}). 1 indicates $z = P^{-1}\vec{v}$.
           real(r8), dimension(*), intent(in)    :: v
             !! An array of length `n` to be multiplied
           real(r8), dimension(*), intent(out)   :: z
             !! An array of length n containing the desired product on
             !! output.
           real(r8), dimension(*), intent(inout) :: rpar
             !! Parameter/work array 
           integer, dimension(*), intent(inout)  :: ipar
             !! Parameter/work array
           integer, intent(out)                  :: itrmjv
             !! Termination flag. 0 indcates normal termination, 1
             !! indicatesfailure to prodce $J\vec{v}$, and 2 indicates
             !! failure to produce $P^{-1}\vec{v}$/
         end subroutine jacv
      end interface

      real(r8), intent(in)                  :: ftol
        !! Stopping tolerance of the f-norm
      real(r8), intent(in)                  :: stptol
        !! Stopping tolerance of the step-length
      integer, dimension(10), intent(in)    :: input
        !! Array containing various user-specified inputs; see below
      integer, dimension(6), intent(out)    :: info
        !! Array containing various outputs; see below
      real(r8), dimension(*), intent(inout) :: rwork
        !! Work array with length depending on the solver used as follows:
        !!
        !!GMRES
      !!:    $n\times(\text{kdmax}+5)+\text{kdmax}\times(\text{kdmax}+3)$,
        !!     where kdmax is the maximum Krylove subspace dimension, either
        !!     the default value of 20 or another value specified by the user
        !!BiCGSTAB
        !!:    $11n$
        !!TFQMR
        !!:    $14n$
      real(r8), dimension(*), intent(inout) :: rpar
        !! Parameter/work array passed to the `f` and `jacv` routines
      integer, dimension(*), intent(inout)  :: ipar
        !! Parameter/work array passed to the `f` and `jacv` routines
      integer, intent(out)                  :: iterm
        !! Termination flag. Values have the following meanings:
        !!
        !!-k
        !!:    illegal value in `input(k)`
        !!0
        !!:    normal termination: $||F|| < \text{ftol}$ or $||\text{step}||
        !!     < \text{stptol}$
        !!1
        !!:    `nnimax` nonlinar iterations reached without success
        !!2
        !!:    failure to evaluate $F$
        !!3
        !!:    in `nitjv`, $J\vec{v}$  failure
        !!4
        !!:    in `nitjv`, $P^{-1}\vec{v}$ failure
        !!5
        !!:    in `nitdrv`, insufficient initial model norm reduction for
        !!     adequate progress. **Note:** This can occur for several
        !!     reasons; examine `itrmks` on return from the Krylov solver
        !!     for further information. (This will be printed out if
        !!     $\text{iplvl}\ge 3$; see the discussion of optional common
        !!     blocks below.)
        !!6
        !!:    in `nitbt`, failure to reach an acceptable step through
        !!     backtracking
      
      interface
        function dinpr(n, x, sx, y, sy)
          !! User-supplied function for calculating vector inner products.
          !! This has the same interace as the BLAS routine
          !! [ddot](http://www.netlib.org/lapack/explore-html/de/da4/group__double__blas__level1_ga75066c4825cb6ff1c8ec4403ef8c843a.html).
          !! If the Euclidean inner product is desired then user can link
          !! to a local BLAS library and provide the name `ddot` to `nitsol`.
          !! `dinpr` must be declared as an external function that returns
          !! a double precision in the calling program.
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
          real(r8)                           :: dinpr
          !! Inner product of `x` and `y`
        end function dinpr

        function dnorm(n, x, sx)
          !! User-supplied function for calculating vector norms. This
          !! has the same interface as the BLAS routine dnrm2; if the
          !! Euclidean norm is desired the user can link to a local
          !! BLAS library and provide the name dnrm2 to nitsol.  dnorm
          !! must be declared as an external function that returns a
          !! double precision value in the calling program.
          import :: r8
          integer, intent(in)                :: n
            !! The length of the array
          real(r8), dimension(*), intent(in) :: x
            !! The input vector
          integer, intent(in)                :: sx
            !! The stride in memory between consecutive elements of `x`
          real(r8) :: dnorm
            !! The vector norm of `x`
        end function dnorm
      end interface

    end subroutine nitsol
  end interface 

end module nitsol_mod

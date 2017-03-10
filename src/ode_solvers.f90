!
!  ode_solvers.f90
!  This file is part of ISOFT.
!  
!  Copyright 2017 Chris MacMackin <cmacmackin@gmail.com>
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

module ode_solvers_mod
  !* Author: Christopher MacMackin
  !  Date: March 2017
  !  License: GPLv3
  !
  ! Provides routines to solve systems of ODEs.
  !
  use iso_fortran_env, only: r8 => real64
  use nitsol_mod, only: gmres_solve
  implicit none

  abstract interface
    function L_intr(u)
      !! An interface for the (linear) left-hand-side of an ODE
      !! being solved by quasilinearisation.
      import :: r8
      implicit none
      real(r8), dimension(:), intent(in) :: u
        !! The state vector for the system of differential equations
      real(r8), dimension(size(u)) :: L_intr
    end function L_intr

    function f_intr(u)
      !! An interface for the (nonlinear) right-hand-side of an ODE
      !! being solved by quasilinearisation.
      import :: r8
      implicit none
      real(r8), dimension(:,:), intent(in) :: u
        !! The state vector for the system of differential equations,
        !! and its derivatives. Column \(i\) represents the \(i-1\)
        !! derivative.
      real(r8), dimension(size(u,1)) :: f_intr
    end function f_intr

    function diff_intr(u, n)
      !! An interface for a function evaluating the derivative of the
      !! state vector.
      import :: r8
      implicit none
      real(r8), dimension(:), intent(in) :: u
        !! The state vector for the system of differential equations
      integer, intent(in)                :: n
        !! The order of the derivative to take
      real(r8), dimension(size(u)) :: diff_intr
    end function diff_intr
    
    function pre_intr(v, u, L, f, Lcur, fcur)
      !! An interface for a preconditioner to be used with the
      !! quasilinearisation ODE solver.
      import :: r8
      implicit none
      real(r8), dimension(:), intent(in)   :: v
        !! The vector to be preconditioned.
      real(r8), dimension(:,:), intent(in) :: u
        !! The current state vector for the system of differential
        !! equations, and its derivatives. Column \(i\) represents the
        !! \(i-1\) derivative.
      procedure(L_intr)                    :: L
        !! The linear, left-hand-side of the ODE being solved.
      procedure(f_intr)                    :: f
        !! The nonlinear, right-hand-side of the ODE being solved.
      real(r8), dimension(:), intent(in)   :: Lcur
        !! The result of `L(u(:,1))`
      real(r8), dimension(:), intent(in)   :: fcur
        !! The result of `f(u)`
      real(r8), dimension(size(v)) :: pre_intr
        !! The result of applying the preconditioner.
    end function pre_intr

  end interface

contains

  subroutine quasilinear_solve(L, f, solution, order, resid_norm, &
                               flag, tol, precond, differentiate, &
                               iter_max, gmres_iter_max, krylov_dim)
    !* Author: Chris MacMackin
    !  Date: March 2017
    !
    ! This is an iterative method to solve nonlinear systems of
    ! ODEs. Consider $$ L^{(n)}\vec{u}(x) = \vec{f}(\vec{u}(x),
    ! \vec{u}^{(1)}(x), \ldots, \vec{u}^{n-1}(x), x), \qquad \vec{u}
    ! \in \mathbb{R}^{m}. $$ For a domain from \( (0,b) \), boundary
    ! conditions are specified by $$ g_k(\vec{u}(0), \vec{u}^{(1)}(0),
    ! \ldots, \vec{u}^{n-1}(0)), \qquad k = 1, \ldots, l $$ $$
    ! g_k(\vec{u}(b), \vec{u}^{(1)}(b), \ldots, \vec{u}^{n-1}(b)),
    ! \qquad k = l+1, \ldots, mn. $$ Here, \(L^{(n)}\) is an \(n\)th
    ! order ordinary differential operator, and \(\vec{f}, g_1, g_2,
    ! \ldots, g_{nm}\) are nonlinear functions of \(\vec{u}(x)\) and
    ! its first \(n-1\) derivatives.
    ! 
    ! With quasilinearisation method, we iterate according to $$
    ! L^{(n)}\vec{u}_{r+1}(x) = \vec{f}(\vec{u}_r(x),
    ! \vec{u}^{(1)}_r(x), \ldots, \vec{u}^{n-1}_r(x), x) + \\
    ! \qquad\sum_{s=0}^{n-1}\vec{f}_{\vec{u}^{(s)}}(\vec{u}_r(x),
    ! \vec{u}^{(1)}_r(x), \ldots, \vec{u}^{n-1}_r(x),
    ! x)(\vec{u}^{(s)}_{r+1}(x) - \vec{u}^{(s)}_r(x)). $$ The boundary
    ! conditions are iterated according to $$ \qquad \sum_{s=0}^{n-1}
    ! g_{k,\vec{u}^{(s)}}(\vec{u}_r(0), \vec{u}^{(1)}_r(0), \ldots,
    ! \vec{u}^{n-1}_r(0), x)\cdot(\vec{u}^{(s)}_{r+1}(0) -
    ! \vec{u}^{(s)}_r(0)) = 0, k=1,\ldots,l$$ and $$ \qquad
    ! \sum_{s=0}^{n-1} g_{k,\vec{u}^{(s)}}(\vec{u}_r(b),
    ! \vec{u}^{(1)}_r(b), \ldots, \vec{u}^{n-1}_r(b),
    ! x)\cdot(\vec{u}^{(s)}_{r+1}(b) - \vec{u}^{(s)}_r(b)) = 0,
    ! k=l+1,\ldots,mn.$$ Here, \(\vec{f}_{\vec{u}^{(s)}}\) is the
    ! derivative of \(\vec{f}\) with respect to \(\vec{u}^{(s)}\) and
    ! is a tensor, while \(g_{k,\vec{u}^{(s)}}\) is the derivative of
    ! \(g_k\) with respect to \(\vec{u}^{(s)}\) and is a vector.
    ! 
    ! The user must provide functions for \(L\) and
    ! \(\vec{f}\). Currently this implementation only handles linear
    ! boundary conditions, which remain the same between
    ! iterations. The boundary conditions should be handled in the
    ! functions passed as arguments. In most situations, a
    ! preconditioner will also be needed. A simple one is to
    ! approximate \(L^{-1}\). The derivatives of \(\vec{f}\) and
    ! \(g_k\) are estimated using a finite-difference.
    !
    procedure(L_intr)                     :: L
      !! A function providing the linear, left-hand-side of the ODE
      !! being solved.
    procedure(f_intr)                     :: f
      !! A function providing the nonlinear, right-hand-side of the
      !! ODE being solved.
    real(r8), dimension(:), intent(inout) :: solution
      !! On input, an estimate of the solution to the ODE. On output,
      !! the actual solution.
    integer, intent(in)                   :: order
      !! The order of the derivative taken by `L`
    real(r8), intent(out)                 :: resid_norm
      !! Norm of the residual of the final solution.
    integer, intent(out)                  :: flag
      !! Status flag indicating whether the iterations ended succesfully.
      !!
      !!< 0:
      !!:    `|flag|` is [[gmres_solve]] return code
      !!
      !!0
      !!:    Normal termination: acceptable solution found
      !!
      !!1
      !!:    Convergence, but residual greater than required tolerance
      !!
      !!2
      !!:    Solution did not converge within `iter_max` iterations
      !!
      !!3
      !!:    No `diff` procedure provided when `order > 1`
      !!
    real(r8), intent(in), optional        :: tol
      !! The tolerance for the solution. Default is `size(solution) * 1e-8`.
    procedure(pre_intr), optional         :: precond
      !! A right-preconditioner which may be used to improve
      !! convergence of the solution.
    procedure(diff_intr), optional        :: differentiate
      !! A procedure which will evaluate the `n`th derivative of the
      !! state vector, when `n` is less than `order`.
    integer, intent(in), optional         :: iter_max
      !! Maximum allowable number of quasilinearised
      !! iterations. Default is 15.
    integer, intent(in), optional         :: gmres_iter_max
      !! Maximum allowable number of GMRES iterations. Default is
      !! 1000.
    integer, intent(in), optional         :: krylov_dim
      !! Maximum Krylov subspace dimension; default 10. Larger values
      !! will allow for faster convergence (and in some cases be the
      !! difference between whether or not convergence is possible),
      !! but require more memory.

    integer :: npoints, itmax, gitmax, kdim
    real(r8) :: eta

    real(r8), dimension(size(solution),order) :: u, u_prev

    if (.not. present(differentiate) .and. order > 1) then
      flag = 3
      return
    end if

    npoints = size(solution)
    if (present(tol)) then
      eta = tol
    else
      eta = 1.e-8_r8 * npoints
    end if
    if (present(iter_max)) then
      itmax = iter_max
    else
      itmax = 15
    end if
    if (present(gmres_iter_max)) then
      gitmax = gmres_iter_max
    else
      gitmax = 1000
    end if
    if (present(krylov_dim)) then
      kdim = krylov_dim
    else
      kdim = 10
    end if
    
    

  end subroutine quasilinear_solve

end module ode_solvers_mod

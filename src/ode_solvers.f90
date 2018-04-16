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
  use nitsol_mod, only: gmres_solve, dnrm2, iplvl
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

    function jac_intr(u, du)
      !! An interface for the product of the Jacobian of the
      !! (nonlinear) right-hand-side of an ODE and another vector.
      import :: r8
      implicit none
      real(r8), dimension(:,:), intent(in) :: u
        !! The state vector for the system of differential equations,
        !! and its derivatives, for which the Jacobian should be
        !! evaluated. Column \(i\) represents the \(i-1\) derivative.
      real(r8), dimension(:,:), intent(in) :: du
        !! The state vector for the system of differential equations,
        !! and its derivatives, which the Jacobian operates on. Column
        !! \(i\) represents the \(i-1\) derivative.
      real(r8), dimension(size(u,1)) :: jac_intr
    end function jac_intr

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
    
    function pre_intr(v, u, L, f, fcur, rhs)
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
      real(r8), dimension(:), intent(in)   :: fcur
        !! The result of `f(u)`
      real(r8), dimension(:), intent(in)   :: rhs
        !! The right hand side of the linear system being
        !! preconditioned.
      real(r8), dimension(size(v)) :: pre_intr
        !! The result of applying the preconditioner.
    end function pre_intr

  end interface

contains

  subroutine quasilinear_solve(L, f, jac_prod, solution, order, resid_norm, &
                               flag, info, tol, precond, differentiate,     &
                               iter_max, gmres_iter_max, krylov_dim)
    !* Author: Chris MacMackin
    !  Date: March 2017
    !
    ! This is an iterative method to solve nonlinear systems of
    ! ODEs. Consider $$ L^{(n)}\vec{u}(x) = \vec{f}(\vec{u}(x),
    ! \vec{u}^{(1)}(x), \ldots, \vec{u}^{n-1}(x), x), \qquad \vec{u}
    ! \in \mathbb{R}^{m}. $$ For a domain from \( (0,b) \), boundary
    ! conditions are specified by $$ g_k(\vec{u}(0) = 0,
    ! \vec{u}^{(1)}(0), \ldots, \vec{u}^{n-1}(0)), \qquad k = 1,
    ! \ldots, l $$ $$ g_k(\vec{u}(b), \vec{u}^{(1)}(b), \ldots,
    ! \vec{u}^{n-1}(b)) = 0, \qquad k = l+1, \ldots, mn. $$ Here,
    ! \(L^{(n)}\) is an \(n\)th order ordinary differential operator,
    ! and \(\vec{f}, g_1, g_2, \ldots, g_{nm}\) are nonlinear
    ! functions of \(\vec{u}(x)\) and its first \(n-1\) derivatives.
    ! 
    ! With quasilinearisation method, we iterate according to $$
    ! L^{(n)}\vec{u}_{r+1}(x) = \vec{f}(\vec{u}_r(x),
    ! \vec{u}^{(1)}_r(x), \ldots, \vec{u}^{n-1}_r(x), x) + \\
    ! \qquad\sum_{s=0}^{n-1}\vec{f}_{\vec{u}^{(s)}}(\vec{u}_r(x),
    ! \vec{u}^{(1)}_r(x), \ldots, \vec{u}^{n-1}_r(x),
    ! x)(\vec{u}^{(s)}_{r+1}(x) - \vec{u}^{(s)}_r(x)). $$ The boundary
    ! conditions are iterated according to $$ \sum_{s=0}^{n-1}
    ! g_{k,\vec{u}^{(s)}}(\vec{u}_r(0), \vec{u}^{(1)}_r(0), \ldots,
    ! \vec{u}^{n-1}_r(0), x)\cdot(\vec{u}^{(s)}_{r+1}(0) -
    ! \vec{u}^{(s)}_r(0)) = 0, k=1,\ldots,l$$ and $$ \sum_{s=0}^{n-1}
    ! g_{k,\vec{u}^{(s)}}(\vec{u}_r(b), \vec{u}^{(1)}_r(b), \ldots,
    ! \vec{u}^{n-1}_r(b), x)\cdot(\vec{u}^{(s)}_{r+1}(b) -
    ! \vec{u}^{(s)}_r(b)) = 0, k=l+1,\ldots,mn.$$ Here,
    ! \(\vec{f}_{\vec{u}^{(s)}}\) is the derivative of \(\vec{f}\)
    ! with respect to \(\vec{u}^{(s)}\) and is a tensor, while
    ! \(g_{k,\vec{u}^{(s)}}\) is the derivative of \(g_k\) with
    ! respect to \(\vec{u}^{(s)}\) and is a vector.
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
    !####Output parameters
    !
    ! On output, the components of the `info` argument are as follows:
    !
    !     info(1)   = nLe   (number of evaluations of `L`)
    !     info(2)   = nfe   (number of evaluations of `f`)
    !     info(3)   = nrpre (number of preconditioner evaluations)
    !     info(4)   = nli   (number of linear iterations)
    !     info(5)   = nni   (number of nonlinear iterations)
    !
    procedure(L_intr)                            :: L
      !! A function providing the linear, left-hand-side of the ODE
      !! being solved.
    procedure(f_intr)                            :: f
      !! A function providing the nonlinear, right-hand-side of the
      !! ODE being solved.
    procedure(jac_intr)                          :: jac_prod
      !! A function providing the product of the Jacobian of the
      !! nonlinear, right-hand-side of the ODE being solved and
      !! another vector.
    real(r8), dimension(:), intent(inout)        :: solution
      !! On input, an estimate of the solution to the ODE. On output,
      !! the actual solution.
    integer, intent(in)                          :: order
      !! The order of the derivative taken by `L`
    real(r8), intent(out)                        :: resid_norm
      !! Norm of the residual of the final solution.
    integer, intent(out)                         :: flag
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
      !!:    Solution began to diverge
      !!
      !!4
      !!:    No `diff` procedure provided when `order > 1`
      !!
    integer, dimension(5), intent(out), optional :: info
      !! Array containing various outputs; see above
    real(r8), intent(in), optional               :: tol
      !! The required reduction in the solution residual. Default is
      !! `size(solution) * 1e-8`.
    procedure(pre_intr), optional                :: precond
      !! A right-preconditioner which may be used to improve
      !! convergence of the solution.
    procedure(diff_intr), optional               :: differentiate
      !! A procedure which will evaluate the `n`th derivative of the
      !! state vector, when `n` is less than `order`.
    integer, intent(in), optional                :: iter_max
      !! Maximum allowable number of quasilinearised
      !! iterations. Default is 15.
    integer, intent(in), optional                :: gmres_iter_max
      !! Maximum allowable number of GMRES iterations. Default is
      !! 1000.
    integer, intent(in), optional                :: krylov_dim
      !! Maximum Krylov subspace dimension; default 10. Larger values
      !! will allow for faster convergence (and in some cases be the
      !! difference between whether or not convergence is possible),
      !! but require more memory.

    integer :: npoints, itmax, gitmax, kdim
    real(r8) :: eta, gmres_eta
    real(r8) :: unorm, eps=1e-5_r8
    real(r8), parameter :: eps_m = epsilon(1._r8)

    integer :: i, stagnant_iters, gmres_flag
    integer :: nlhs, nrpre, nli, tnlhs, tnrpre, tnli
    real(r8) :: init_resid, old_resid, gmres_norm
    real(r8), dimension(size(solution),order) :: u_prev
    real(r8), dimension(size(solution))       :: f_prev, rhs

    tnlhs  = 0
    tnrpre = 0
    tnli   = 0

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

    i = 0
    stagnant_iters = 0

    u_prev = get_derivs(solution)
    unorm = dnrm2(npoints, u_prev, 1)
    f_prev = f(u_prev)
    resid_norm = dnrm2(npoints, L(solution) - f_prev, 1)
    init_resid = resid_norm
    old_resid = resid_norm * 1e3_r8
!    print*, L(solution) - f_prev

!iplvl=4
    do while(resid_norm > eta)
      print*, resid_norm, tnli
      !print*, L(solution) - f_prev
      i = i + 1
      if (abs(old_resid - resid_norm)/resid_norm < 1e-2_r8) then
        stagnant_iters = stagnant_iters + 1
      else
        stagnant_iters = 0
      end if
      if (stagnant_iters > 3) then
        flag = 1
        return
      end if
      if (i > itmax) then
        flag = 2
        return
      end if
      if (20*old_resid < resid_norm) then
        flag = 3
!    print*, L(solution) - f_prev
        return
      end if

      !eps = sqrt((1+unorm)*eps_m)*max(10._r8, 10._r8**(5-i))/unorm
      !print*,'RHS epsilon',eps
      !rhs = f_prev - (f(u_prev + eps*u_prev) - f_prev)/eps
      rhs = f_prev - (f(u_prev + epsilon*u_prev) - f_prev)/epsilon
      gmres_eta = max(min(eta*10._r8**min(i+2,6),1e-4_r8),1e-10_r8)
      gmres_eta = gmres_eta * 10._r8**(-2*stagnant_iters)
      call gmres_solve(solution, lin_op, rhs, gmres_norm, gmres_flag, &
                       nlhs, nrpre, nli, gmres_eta, preconditioner,   &
                       iter_max=gitmax, krylov_dim=kdim)

      tnlhs  = tnlhs  + nlhs
      tnrpre = tnrpre + nrpre
      tnli   = tnli   + nli
      if(gmres_flag > 0) print*, 'Warning, GMRES returned with flag', gmres_flag, tnli

      u_prev = get_derivs(solution)
      unorm = dnrm2(npoints, u_prev, 1)
      f_prev = f(u_prev)
    !print*,f_prev(size(f_prev)-10:)
      old_resid = resid_norm
      resid_norm = dnrm2(npoints, L(solution) - f_prev, 1)

      if (gmres_flag /= 0 .and. resid_norm > old_resid) then
        if (present(info)) then
          info(1) = i + tnlhs + tnrpre
          info(2) = 2*i + tnlhs
          info(3) = tnrpre
          info(4) = tnli
          info(5) = i
        end if
        if (init_resid < resid_norm) then
          flag = 3
        else
          flag = -gmres_flag
        end if
        return
      end if
    end do

    if (present(info)) then
      info(1) = 1 + i + tnlhs + tnrpre
      info(2) = 1 + 2*i + tnlhs
      info(3) = tnrpre
      info(4) = tnli
      info(5) = i
    end if
    flag = 0

  contains

    function get_derivs(v)
      !! Calculates the necessary number of derivatives and assembles
      !! them in 2-D array.
      real(r8), dimension(:), intent(in) :: v
      real(r8), dimension(size(v), order) :: get_derivs
      integer :: j
      get_derivs(:,1) = v
      do j = 2, order
        get_derivs(:,j) = differentiate(v, i-1)
      end do
    end function get_derivs

    function lin_op(v, xcur, rhs, rpar, ipar, success)
      !! The linear operator for the quasilinearised system.
      real(r8), dimension(:), intent(in)    :: v
        !! The vector to be operated upon
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
      real(r8), dimension(size(xcur))       :: lin_op
        !! Result of the operation
      real(r8), dimension(size(xcur),order) :: v_derivs
      real(r8) :: vnorm
      v_derivs = get_derivs(v)
      vnorm = dnrm2(npoints, v_derivs, 1)
      !eps = sqrt((1+unorm)*eps_m)*max(10._r8, 10._r8**(5-i))/vnorm
      !print*,'LHS epsilon',eps
      !lin_op = L(v) - (f(u_prev + eps*v_derivs) - f_prev)/eps
      lin_op = L(v) - (f(u_prev + epsilon*v_derivs) - f_prev)/epsilon
      success = .true.
    end function lin_op

    function preconditioner(v, xcur, rhs, rpar, ipar, success)
      !! The preconditioner for the quasilinearised system.
      real(r8), dimension(:), intent(in)    :: v
        !! The vector to be preconditioned
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
      real(r8), dimension(size(xcur))       :: preconditioner
        !! Result of the operation
      preconditioner = precond(v, reshape(xcur, [size(xcur), 1]), L, f, f_prev, rhs)
      success = .true.
    end function preconditioner

  end subroutine quasilinear_solve

end module ode_solvers_mod

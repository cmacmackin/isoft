!
!  pseudospectral_block.f90
!  This file is part of ISOFT.
!  
!  Copyright 2018 Chris MacMackin <cmacmackin@gmail.com>
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

module coriolis_block_mod
  !* Author: Christopher MacMackin
  !  Date: march 2017
  !  License: GPLv3
  !
  ! Provides a derived type which representes the operator acting on
  ! velocity and its derivative with the Coriolis force. This can be
  ! used for preconditioning in the [[plume]] solver.
  !
  use iso_fortran_env, only: r8 => real64
  use factual_mod!, only: abstract_field, scalar_field, vector_field
  use f95_lapack, only: la_gesvx
  use pseudospectral_block_mod, only: pseudospec_block
  use penf, only: str
  use logger_mod, only: logger => master_logger
  implicit none
  private

  character(len=1), parameter :: trans = 'N'
    !! The LAPACK parameter indicating not to operate on the transpose
    !! of the matrix when solving for boundary conditions.
  
  type(uniform_scalar_field) :: zero

  type, public :: coriolis_block
    !* Author: Chris MacMackin
    !  Date: January 2017
    !
    ! A data type representing a matrix operator for the momentum
    ! components of the linear parts plume equations, with the
    ! Coriolis force. It can be useful when preconditioning a the
    ! plume solver. It is inherently 1-D in its implementation, but
    ! has a transverse velocity component.
    !
    ! This type solves the linear, coupled, differential equation 
    ! \[ \frac{d}{dx} \bm{Y} = \bm{A}\bm{Y} + \bm{F} \]
    ! where \( \bm{Y} = [U, V, U', V']^{T} \),
    ! \[ \bm{A} = \begin{bmatrix}
    ! 0        & 0         & 1 & 0 \\
    ! 0        & 0         & 0 & 1 \\
    ! 0        & -\Phi/\nu & 0 & 0 \\
    ! \Phi/\nu & 0         & 0 & 0
    ! \end{bmatrix},
    ! \]
    ! and \( \bm{F} \) is the input value to which this block is 
    ! applied in the preconditioner. \(\Phi\) is the dimensionless
    ! coriolis parameter, while \(\nu\) is the eddy viscosity.
    !
    ! The matrix \(\bm{A}\) can be diagonalised by factoring it such
    ! that \(\bm{A} = \bm{V}\bm{D}\bm{V}^{-1}\), where \(\bm{V}\) is a
    ! change of basis matrix with columns made up of the eigenvectors 
    ! of \(\bm{A}\) and \(\bm{D}\) is a diagonal matrix made up of
    ! the corresponding eigenvalues. These have the values
    ! \[ \bm{V} = \begin{bmatrix}
    ! \tfrac{-1-i}{\alpha} & \tfrac{-1+i}{\alpha} & \tfrac{1-i}{\alpha} & \tfrac{1+i}{\alpha} \\
    ! \tfrac{-1+i}{\alpha} & \tfrac{-1-i}{\alpha} & \tfrac{1+i}{\alpha} & \tfrac{1-i}{\alpha} \\
    ! i                    & -i                   & -i                  & i                   \\
    ! 1                    & 1                    & 1                   & 1
    ! \end{bmatrix},
    ! \]
    ! \[ \bm{D} = \beta\begin{bmatrix}
    ! -1-i & 0    & 0   & 0   \\
    ! 0    & -1+i &     & 0   \\
    ! 0    & 0    & 1-i & 0   \\
    ! 0    & 0    & 0   & 1+i
    ! \end{bmatrix},
    ! \]
    ! \[ \bm{V}^{-1} = \frac{1}{4}\begin{bmatrix}
    ! (-1+i)\beta & -(1+i)\beta & -i & 1 \\
    ! (-1-i)\beta & -(1-i)\beta & i  & 1 \\
    ! (1+i)\beta  & (1-i)\beta  & i  & 1 \\
    ! (1-i)\beta  & (1+i)\beta  & -i & 1
    ! \end{bmatrix},
    ! \]
    ! \[ \alpha = \sqrt{\frac{2\Phi}{\nu}}, \quad \beta = \sqrt{\frac{\Phi}{2\nu}}. \]
    ! It can be shown that the solution to the differential equation is
    ! \[ \bm{Y} = \bm{V}\left(e^{\bm{D}x}\bm{B} + 
    !    e^{\bm{D}x}\int_0^{x}e^{\bm{D}x'}\bm{V}^{-1}\bm{F}dx' \right),\]
    ! where \(\bm{B}\in\mathbb{R}^{4}\) is a vector chosen to satisfy
    ! the boundary conditions on the system. It can be found by solving a
    ! 4×4 linear system.
    !
    ! This type inherits 
    !
    private
    real(r8), dimension(4)   :: D_r
      !! Real component of the diagonal matrix, \(\bm{D}\), with only
      !! diagonal values stored
    real(r8), dimension(4)   :: D_i
      !! Imaginary component of the diagonal matrix, \(\bm{D}\), with
      !! only diagonal values stored
    type(cheb1d_scalar_field), dimension(4,4) :: emDxVinv_r
      !! Real component of \(e^{-\bm{D}x}\bm{V}^{-1}\)
    type(cheb1d_scalar_field), dimension(4,4) :: emDxVinv_i
      !! Imaginary component of \(e^{-\bm{D}x}\bm{V}^{-1}\)
    type(cheb1d_scalar_field), dimension(4) :: eDx_r
      !! Real component of \(e^{\bm{D}x}\), with only diagonal values
      !! stored
    type(cheb1d_scalar_field), dimension(4) :: eDx_i
      !! Imaginary component of \(e^{\bm{D}x}\), with only diagonal
      !! values stored
    real(r8), dimension(4,4) :: V_r
      !! Real component of the change of basis matrix, \(\bm{V}\)
    real(r8), dimension(4,4) :: V_i
      !! Imaginary component of the change of basis matrix, \(\bm{V}\)
    type(pseudospec_block)   :: integrator
      !! A pseudospectral differentiation block which can be used to
      !! perform integration
    integer                  :: vel_bound_loc
      !! Location code for the velocity's boundary condition
    integer                  :: dvel_bound_loc
      !! Location code for the velocity derivative's boundary
      !! condition
    integer                  :: integrate_bound
      !! Location from which to perform the integration
    real(r8), dimension(4)   :: xbounds
      !! Boundary location for each component of the solution vector
    complex(r8), dimension(4,4) :: bound_matrix
      !! Matrix for the system to solve in order to satisfy the
      !! boundary conditions
    complex(r8), dimension(4,4) :: bound_matrix_scaled
      !! Matrix for the system to solve in order to satisfy the
      !! boundary conditions, which has been scaled by LAPACK95 
      !! to improve conditioning.
    complex(r8), dimension(4,4) :: factored_matrix
      !! Factored matrix for the system to solve in order to satisfy
      !! the boundary conditions
    integer, dimension(4)    :: pivots
      !! The pivots used in the factorisation of the matrix used to
      !! satisfy boundary conditions
    real(r8), dimension(4)   :: r_scales
      !! Row scale factors from equilibrating the bound_matrix
    real(r8), dimension(4)   :: c_scales
      !! Column scale factors from equilibrating the bound_matrix
    character(len=1)         :: equed
      !! The method used to equilibrate bound_matrix
    integer :: int
  contains
    private
    procedure, public :: solve_for
    procedure         :: assign
    generic, public   :: assignment(=) => assign
  end type coriolis_block

  interface coriolis_block
    module procedure constructor
  end interface coriolis_block

contains

  function constructor(phi, nu, velbound, dvelbound, integrate_bound, template) &
       result(this)
    !* Author: Chris MacMackin
    !  Date: January 2018
    !
    ! Builds a Coriolis block which can be used to solve the inverse
    ! problem for the linear components of the plume momentum
    ! equations. The result can only be used with fields having the
    ! same grid as the template.
    !
    real(r8), intent(in)              :: phi
      !! The dimensionless coriolis parameter
    real(r8), intent(in)              :: nu
      !! The dimensionless eddy diffusivity
    integer, intent(in)               :: velbound
      !! Location code for the velocity's boundary condition. 1
      !! indicates upper boundary, -1 indicates lower boundary.
    integer, intent(in)               :: dvelbound
      !! Location code for the velocity's boundary condition. 1
      !! indicates upper boundary, -1 indicates lower boundary.
    integer, intent(in)               :: integrate_bound
      !! Location code for the boundary to perform integrations
      !! from. This should be the opposite boundary from where
      !! boundary data is stored.
    class(abstract_field), intent(in) :: template
      !! A scalar field with the same grid as any fields passed as
      !! arguments to the [[pseudospec_block(type):solve_for]] method.
    type(coriolis_block)              :: this

    type(cheb1d_scalar_field) :: xvals
    integer :: i, j, info
    real(r8) :: alpha, beta, rcond, s
    real(r8), dimension(:,:), allocatable :: domain
    complex(r8), dimension(4) :: dummy_in
    complex(r8), dimension(4) :: dummy_out
    character(len=:), allocatable :: msg
    call template%guard_temp()

    zero = uniform_scalar_field(0._r8)
    this%integrator = pseudospec_block(template)
    domain = template%domain()
    this%vel_bound_loc = velbound
    this%dvel_bound_loc = dvelbound
    if (velbound == 1) then
      this%xbounds(1:2) = domain(1,2)
    else if (velbound == -1) then
      this%xbounds(1:2) = domain(1,1)
    else
      call logger%fatal('coriolis_block', &
                        'Only boundary location codes 1 or -1 are accepted.')
      error stop
    end if
    if (dvelbound == 1) then
      this%xbounds(3:4) = domain(1,2)
    else if (dvelbound == -1) then
      this%xbounds(3:4) = domain(1,1)
    else
      call logger%fatal('coriolis_block', &
                        'Only boundary location codes 1 or -1 are accepted.')
      error stop
    end if
    this%integrate_bound = integrate_bound
    xvals = cheb1d_scalar_field(template%elements(), linear, domain(1,1), &
                                domain(1,2))
    ! Compute the diagonal elements of \(\bm{B}\)
    alpha = sqrt(abs(2*phi/nu))
    beta = sqrt(abs(0.5_r8*phi/nu))
    s = sign(1._r8,phi)
    this%D_r = [-beta, -beta, beta, beta]
    this%D_i = [-beta, beta, -beta, beta]
    ! Use the associations as temporary work arrays, prior to
    ! computing their actual values
    associate (emDx_r => this%eDx_r, emDx_i => this%eDx_i, & 
               Vinv_r => this%V_r, Vinv_i => this%V_i)
      ! Compute \(e^{-\bm{D}x}\). Need to do this in two steps to avoid a compiler bug
      emDx_i = [(exp(-this%D_r(i)*xvals), i=1,4)]
      emDx_r = [(emDx_i(i)*cos(-this%D_i(i)*xvals), i=1,4)]
      emDx_i = [(emDx_i(i)*sin(-this%D_i(i)*xvals), i=1,4)]
      ! Compute \(\bm{V}^{-1}\)
      Vinv_r = 0.25_r8*reshape([-s*beta, -s*beta, s*beta, s*beta,  &
                                  -beta,   -beta,   beta,   beta,  &
                                  0._r8,   0._r8,  0._r8,  0._r8,  &
                                  1._r8,   1._r8,  1._r8,  1._r8], &
                               [4,4])
      Vinv_i = 0.25_r8*reshape([s*beta, -s*beta, s*beta, -s*beta,  &
                                 -beta,    beta,  -beta,    beta,  &
                                    -s,       s,      s,      -s,  &
                                 0._r8,   0._r8,  0._r8,   0._r8], &
                               [4,4])
      ! Compute \(e^{-\bm{D}x}\bm{V}^{-1}\)
      this%emDxVinv_r = reshape([((emDx_r(i)*Vinv_r(i,j) - emDx_i(i)*Vinv_i(i,j), &
                                   i=1,4), j=1,4)], [4, 4])
      this%emDxVinv_i = reshape([((emDx_r(i)*Vinv_i(i,j) + emDx_i(i)*Vinv_r(i,j), &
                                   i=1,4), j=1,4)], [4, 4])
    end associate
    ! Compute \(e^{\bm{D}x}\). Need to do this in two steps to avoid a compiler bug
    this%eDx_i = [(exp(this%D_r(i)*xvals), i=1,4)]
    this%eDx_r = [(this%eDx_i(i)*cos(this%D_i(i)*xvals), i=1,4)]
    this%eDx_i = [(this%eDx_i(i)*sin(this%D_i(i)*xvals), i=1,4)]
    ! Compute \(\bm{V}\)
    this%V_r = reshape([-s/alpha, -1._r8/alpha, 0._r8, 1._r8,  &
                        -s/alpha, -1._r8/alpha, 0._r8, 1._r8,  &
                         s/alpha,  1._r8/alpha, 0._r8, 1._r8,  &
                         s/alpha,  1._r8/alpha, 0._r8, 1._r8], &
                       [4,4])
    this%V_i = reshape([-s/alpha,  1._r8/alpha,  s, 0._r8,  &
                         s/alpha, -1._r8/alpha, -s, 0._r8,  &
                        -s/alpha,  1._r8/alpha, -s, 0._r8,  &
                         s/alpha, -1._r8/alpha,  s, 0._r8], &
                       [4,4])
    ! Construct and factor matrix used for satisfying boundary conditions
    dummy_in = [(1,0), (0,1), (0.5, 0.5), (-0.5, 0.5)]
    this%bound_matrix = reshape([((cmplx(this%V_r(i,j), this%V_i(i,j), r8)* &
                                   exp(cmplx(this%D_r(j), this%D_i(j), r8)* &
                                   this%xbounds(i)), i=1,4), j=1,4)], [4,4])
    this%bound_matrix_scaled = this%bound_matrix
    call la_gesvx(this%bound_matrix_scaled, dummy_in, dummy_out,  &
                   this%factored_matrix, this%pivots, 'E', trans, &
                   this%equed, this%r_scales, this%c_scales,      &
                   rcond=rcond, info=info)
    if (info /= 0) then
      msg = 'Tridiagonal matrix solver returned with flag '//str(info)
      call logger%error('coriolis_block',msg)
    end if
    msg = 'Boundary matrix factored with estimated condition '// &
          'number '//str(1._r8/rcond)
    call logger%trivia('coriolis_block',msg)
#ifdef DEBUG
    call logger%debug('coriolis_block', &
                      'Successfully constructed Coriolis preconditioner block.')
#endif

    call template%clean_temp()

  contains
    pure function linear(x) result(scalar)
      real(r8), dimension(:), intent(in) :: x 
        !! The position at which this function is evaluated
      real(r8) :: scalar
      scalar = x(1)
    end function linear
  end function constructor

  subroutine solve_for(this, velocity, velocity_dx)
    !* Author: Chris MacMackin
    !  Date: January 2018
    !
    ! Inverts the linear portions of the plume momentum equation with
    ! the provided data. This is done by solving the linear ODE
    ! described in the documentation for the [[coriolis_block]]
    ! type. The block object must first have been initialised using
    ! the constructor.
    !
    ! @Warning Currently this is only implemented for a 1-D field.
    !
    class(coriolis_block), intent(inout) :: this
    class(vector_field), intent(inout)   :: velocity
      !! On input, the velocity value being preconditioned. On output,
      !! the preconditioned velocity.
    class(vector_field), intent(inout)   :: velocity_dx
      !! On input, the velocity derivative being preconditioned. On
      !! output, the preconditioned velocity derivative.

    type(cheb1d_scalar_field), dimension(4) :: F, eDxC_r, eDxC_i, E_r, E_i
    integer :: i
    real(r8), dimension(1) :: rtmp, ctmp
    real(r8), dimension(4) :: bound_vals
    complex(r8), dimension(4) :: rhs, B, C_bounds
    type(cheb1d_scalar_field) :: tmp

    integer :: info
    real(r8) :: rcond

    call velocity%guard_temp(); call velocity_dx%guard_temp()

    ! Construct array of scalar matrices for easier manipulation
    F(1) = velocity%component(1)
    F(2) = velocity%component(2)
    F(3) = velocity_dx%component(1)
    F(4) = velocity_dx%component(2)
    ! Get boundary values
    tmp = F(1)%get_boundary(this%integrate_bound, 1)
    rtmp = tmp%raw()
    bound_vals(1) = rtmp(1)
    tmp = F(2)%get_boundary(this%integrate_bound, 1)
    rtmp = tmp%raw()
    bound_vals(2) = rtmp(1)
    tmp = F(3)%get_boundary(this%integrate_bound, 1)
    rtmp = tmp%raw()
    bound_vals(3) = rtmp(1)
    tmp = F(4)%get_boundary(this%integrate_bound, 1)
    rtmp = tmp%raw()
    bound_vals(4) = rtmp(1)
    ! If have boundary conditions at both boundaries, correct the
    ! input fields
!    if (this%vel_bound_loc /= this%integrate_bound) then
!      F(1) = this%integrator%solve_for(F(1), this%vel_bound_loc, zero)
!      F(1) = F(1)%d_dx(1)
!      F(2) = this%integrator%solve_for(F(2), this%vel_bound_loc, zero)
!      F(2) = F(2)%d_dx(1)
!    end if
!    if (this%dvel_bound_loc /= this%integrate_bound) then
!      F(3) = this%integrator%solve_for(F(3), this%dvel_bound_loc, zero)
!      F(3) = F(3)%d_dx(1)
!      F(4) = this%integrator%solve_for(F(4), this%dvel_bound_loc, zero)
!      F(4) = F(4)%d_dx(1)
!    end if
    ! Calculate \(e^{-\bm{D}x}\bm{V}^{-1}\bm{F}\), aliasing variables
    ! which are not needed yet
    associate (emDxVinvF_r => eDxC_r, emDxVinvF_i => eDxC_i, &
               M_r => this%emDxVinv_r, M_i => this%emDxVinv_i, &
               C_r => E_r, C_i => E_i)
      emDxVinvF_r = [(M_r(i,1)*F(1) + M_r(i,2)*F(2) + M_r(i,3)*F(3) + &
                      M_r(i,4)*F(4), i=1,4)]
      emDxVinvF_i = [(M_i(i,1)*F(1) + M_i(i,2)*F(2) + M_i(i,3)*F(3) + &
                      M_i(i,4)*F(4), i=1,4)]
      ! Integrate \(e^{-\bm{D}x}\bm{V}^{-1}\bm{F}\) to get \(\bm{C}\)
      C_r = [(this%integrator%solve_for(emDxVinvF_r(i), this%integrate_bound, &
              zero), i=1,4)]
      emDxVinvF_r(1) = C_r(4)%d_dx(1)
      C_i = [(this%integrator%solve_for(emDxVinvF_i(i), this%integrate_bound, &
              zero), i=1,4)]
      ! Get the values of \(e^{\bm{D}x}\bm{C}\) at the upper boundary
      do i=1,4
        tmp = C_r(i)%get_boundary(-this%integrate_bound, 1)
        rtmp = tmp%raw()
        tmp = C_i(i)%get_boundary(-this%integrate_bound, 1)
        ctmp = tmp%raw()
        C_bounds(i) = cmplx(rtmp(1), ctmp(1), r8)
      end do
      ! Calculate \(e^{\bm{D}x}\bm{C}\)
      eDxC_r = [(this%eDx_r(i)*C_r(i) - this%eDx_i(i)*C_i(i), i=1,4)]
      eDxC_i = [(this%eDx_r(i)*C_i(i) + this%eDx_i(i)*C_r(i), i=1,4)]
    end associate
    
    ! Compute RHS for the linear system satisfying the boundary conditions
    if (this%vel_bound_loc == this%integrate_bound) then
      rhs(1:2) = bound_vals(1:2)
    else
      rhs(1:2) = bound_vals(1:2) - matmul(this%bound_matrix(1:2,:), &
                                          C_bounds)
    end if
    if (this%dvel_bound_loc == this%integrate_bound) then
      rhs(3:4) = bound_vals(3:4)
    else
      rhs(3:4) = bound_vals(3:4) - matmul(this%bound_matrix(3:4,:), &
                                          C_bounds)
    end if

    ! Compute coefficients for inhomogeneous components of solution so
    ! that boundary conditions are satisfied
    call la_gesvx(this%bound_matrix_scaled, rhs, B, this%factored_matrix, &
                  this%pivots, 'F', trans, this%equed, this%r_scales,     &
                  this%c_scales, info=i)
!print*,this%D_r
!print*,this%D_i
!print*,B
    ! Calculate \(\bm{E} = e^{\bm{D}x}\bm{B} + e^{\bm{D}x}\bm{C}\)
    E_r = [(this%eDx_r(i)*real(B(i)) - this%eDx_i(i)*aimag(B(i)) + &
            eDxC_r(i), i=1,4)]
    E_i = [(this%eDx_i(i)*real(B(i)) + this%eDx_r(i)*aimag(B(i)) + &
            eDxC_i(i), i=1,4)]

    ! Transform back to proper basis to get proper solution
    F = [(this%V_r(i,1)*E_r(1) - this%V_i(i,1)*E_i(1) + &
          this%V_r(i,2)*E_r(2) - this%V_i(i,2)*E_i(2) + &
          this%V_r(i,3)*E_r(3) - this%V_i(i,3)*E_i(3) + &
          this%V_r(i,4)*E_r(4) - this%V_i(i,4)*E_i(4), i=1,4)]
    velocity = F(1:2)
    velocity_dx = F(3:4)
    call velocity%clean_temp(); call velocity_dx%clean_temp()
  end subroutine solve_for

  subroutine assign(this, rhs)
    !* Author: Chris MacMackin
    !  Date: January 2017
    !
    ! Safely assigns the value of one coriolis block to another.
    !
    class(coriolis_block), intent(inout) :: this
    class(coriolis_block), intent(in)  :: rhs
      !! The value being assigned
    this%D_r = rhs%D_r
    this%D_i = rhs%D_i
    this%emDxVinv_r = rhs%emDxVinv_r
    this%emDxVinv_i = rhs%emDxVinv_i
    this%eDx_r = rhs%eDx_r
    this%eDx_i = rhs%eDx_i
    this%V_r = rhs%V_r
    this%V_i = rhs%V_i
    this%integrator = rhs%integrator
    this%vel_bound_loc = rhs%vel_bound_loc
    this%dvel_bound_loc = rhs%dvel_bound_loc
    this%integrate_bound = rhs%integrate_bound
    this%xbounds= rhs%xbounds
    this%bound_matrix = rhs%bound_matrix
    this%bound_matrix_scaled = rhs%bound_matrix_scaled
    this%factored_matrix = rhs%factored_matrix
    this%pivots = rhs%pivots
    this%r_scales = rhs%r_scales
    this%c_scales = rhs%c_scales
    this%equed = rhs%equed
  end subroutine assign

end module coriolis_block_mod

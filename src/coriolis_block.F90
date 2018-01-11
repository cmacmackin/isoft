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
  use chebyshev_mod, only: collocation_points, integrate_1d
  use pseudospectral_block_mod, only: pseudospec_block
  use penf, only: str
  use logger_mod, only: logger => master_logger
  implicit none
  private

  integer, parameter :: no_extra_derivative = -1
  
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
    ! 0         & 0        & 1 & 0 \\
    ! 0         & 0        & 0 & 1 \\
    ! 0         & \Phi/\nu & 0 & 0 \\
    ! -\Phi/\nu & 0        & 0 & 0
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
    ! \tfrac{1+i}{\alpha}  & \tfrac{1-i}{\alpha}  & \tfrac{-1+i}{\alpha} & \tfrac{-1-i}{\alpha} \\
    ! \tfrac{-1+i}{\alpha} & \tfrac{-1-i}{\alpha} & \tfrac{1+i}{\alpha}  & \tfrac{1-i}{\alpha}  \\
    ! -i                   & i                    & i                    & -i                   \\
    ! 1                    & 1                    & 1                    & 1
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
    ! (1-i)\beta  & -(1+i)\beta & i  & 1 \\
    ! (1+i)\beta  & -(1-i)\beta & -i & 1 \\
    ! -(1+i)\beta & (1-i)\beta  & -i & 1 \\
    ! -(1-i)\beta & (1+i)\beta  & i  & 1
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
    class(scalar_field), allocatable :: xvals_field
      !! Coordinates at which to find the solution (same as xvals in
      !! [[pseudospec_block]], except stored as a field)
    real(r8), dimension(4)   :: D_r
      !! Real component of the diagonal matrix, \(\bm{D}\), with only
      !! diagonal values stored
    real(r8), dimension(4)   :: D_i
      !! Imaginary component of the diagonal matrix, \(\bm{D}\), with
      !! only diagonal values stored
    real(r8), dimension(4,4) :: V_r
      !! Real component of the change of basis matrix, \(\bm{V}\)
    real(r8), dimension(4,4) :: V_i
      !! Imaginary component of the change of basis matrix, \(\bm{V}\)
    real(r8), dimension(4,4) :: Vinv_r
      !! Real component of the inverse change of basis matrix,
      !! \(\bm{V}^{-1}\)
    real(r8), dimension(4,4) :: Vinv_i
      !! Imaginary component of the inverse change of basis matrix,
      !! \(\bm{V}^{-1}\)
    type(pseudospec_block)   :: integrator
      !! A pseudospectral differentiation block which can be used to
      !! perform integration.
  contains
    private
  end type coriolis_block

  interface coriolis_block
    module procedure constructor
  end interface coriolis_block

contains

  function constructor(template) result(this)
    !* Author: Chris MacMackin
    !  Date: January 2018
    !
    ! Builds a Coriolis block which can be used to solve the inverse
    ! problem for the linear components of the plume momentum
    ! equations. The result can only be used with fields having the
    ! same grid as the template.
    !
    class(abstract_field), intent(in)     :: template
      !! A scalar field with the same grid as any fields passed as
      !! arguments to the [[pseudospec_block(type):solve_for]] method.
    type(coriolis_block)                :: this
    this%integrator = pseudospec_block(template)
    ! Do more stuff...
  end function constructor


end module coriolis_block_mod

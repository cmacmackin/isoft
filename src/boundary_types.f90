!
!  boundary_types.f90
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

module boundary_types_mod
  !* Author: Christopher MacMackin
  !  Date: January 2017
  !  License: GPLv3
  !
  ! Provides parameters which can be used to represent different types
  ! of boundary conditions.
  !
  implicit none
  private

  integer, parameter, public :: free_boundary = -1
    !! Indicates that the solution may take any value at the boundary.
  integer, parameter, public :: dirichlet = 0
    !! Indicates that the value of the solution at this boundary is
    !! prescribed.
  integer, parameter, public :: neumann = 1
    !! Indicates that the first derivative of the solution at this
    !! boundary is prescribed.
  integer, parameter, public :: cauchy = 5
    !! Indicates that the value of the solution and its first
    !! derivative are both prescribed at this boundary.
  integer, parameter, public :: robin = 6
    !! Indicates that the linear combination of the solution's value
    !! and first derivative is prescribed at this boundary.

end module boundary_types_mod

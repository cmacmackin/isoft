!
!  melt_relationship.f90
!  This file is part of ISOFT.
!  
!  Copyright 2016 Chris MacMackin <cmacmackin@physics.ox.ac.uk>
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

module melt_relationship_mod
  !* Author: Christopher MacMackin
  !  Date: October 2016
  !  License: GPLv3
  !
  ! Provides an abstract data type to model melting of an ice shelf into a 
  ! vertically integrated plume. 
  !
  use iso_fortran_env, only: r8 => real64
  use factual_mod, only: scalar_field, vector_field
  implicit none
  private

  type, abstract, public :: abstract_melt_relationship
    !* Author: Christopher MacMackin
    !  Date: October 2016
    !
    ! An abstract data type for calculating melting of an ice shelf into
    ! a vertically integrated [[plume]]. The melt rate, as well as
    ! effect on termperature and salinity, are calculated by calling
    ! [[abstract_melt_relationship:solve_for_melt]] and then accessed
    ! using [[abstract_melt_relationship:melt_rate]],
    ! [[abstract_melt_relationship:heat_equation_terms]],
    ! [[abstract_melt_relationship:salt_equation_terms]]. 
    ! 
  contains
    procedure(solve), deferred      :: solve_for_melt
    procedure(get_scalar), deferred :: salt_equation_terms
      !! Returns the terms this melt formulation contributes to the
      !! salt equation, after they have been solved for using
      !! [[abstract_melt_relationship:solve_for_melt]].
    procedure(get_scalar), deferred :: heat_equation_terms
      !! Returns the terms this melt formulation contributes to the
      !! heat equation, after they have been solved for using
      !! [[abstract_melt_relationship:solve_for_melt]]. 
    procedure(get_scalar), deferred :: melt_rate
      !! Returns the melt rate calculated using this formulation,
      !! after it has been solved for using 
      !! [[abstract_melt_relationship:solve_for_melt]]. 
    procedure(has_terms), deferred  :: has_heat_terms
      !! Whether this formulation of melting contributes any terms to
      !! a plume's heat equation.
    procedure(has_terms), deferred  :: has_salt_terms
      !! Whether this formulation of melting contributes any terms to
      !! a plume's salinity equation.
 end type abstract_melt_relationship

  abstract interface
    subroutine solve(this, velocity, pressure, temperature, salinity, &
                     plume_thickness, time)
      import :: abstract_melt_relationship
      import :: scalar_field
      import :: vector_field
      import :: r8
      class(abstract_melt_relationship), intent(inout) :: this
      class(vector_field), intent(in)                  :: velocity
        !! The velocity field of the plume into which fluid is melting.
      class(scalar_field), intent(in)                  :: pressure
        !! The water pressure at the interface where the melting occurs.
      class(scalar_field), intent(in)                  :: temperature
        !! The temperature of the plume into which fluid is melting.
      class(scalar_field), intent(in)                  :: salinity
        !! The salinity of the plume into which fluid is melting.
      class(scalar_field), intent(in)                  :: plume_thickness
        !! The thickness of the plume into which fluid is melting.
      real(r8), intent(in), optional                   :: time
        !! The time at which the melting is being solved for. If not
        !! present then assumed to be same as previous value passed.
    end subroutine solve

    pure function get_scalar(this) result(property)
      import :: abstract_melt_relationship
      import :: scalar_field
      class(abstract_melt_relationship), intent(in) :: this
      class(scalar_field), allocatable :: property
        !! The value of whatever property is being returned.
    end function get_scalar

    pure function has_terms(this)
      import :: abstract_melt_relationship
      class(abstract_melt_relationship), intent(in) :: this
      logical :: has_terms
        !! Whether this formulation of melting contributes terms to
        !! the heat or salinity equations.
    end function has_terms
  end interface

end module melt_relationship_mod

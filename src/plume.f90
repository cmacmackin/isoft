!
!  plume.f90
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

module plume_mod
  !* Author: Christopher MacMackin
  !  Date: April 2016
  !  License: GPLv3
  !
  ! Provides a concrete implementation of the [[basal_surface(type)]] data type,
  ! representing a buoyant plume beneath an ice shelf.
  !
  use iso_fortran_env, only: r8 => real64
  use basal_surface_mod, only: basal_surface
  use factual_mod, only: scalar_field, cheb1d_scalar_field, cheb1d_vector_field
  use entrainment_mod, only: abstract_entrainment
  use melt_relationship_mod, only: abstract_melt_relationship
  use plume_boundary_mod, only: plume_boundary
  use ambient_mod, only: ambient_conditions
  implicit none
  private

  type, extends(basal_surface), public :: plume
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! A concrete implementation of the [[basal_surface(type)]]
    ! abstract data type, representing the buoyant plume beneath an
    ! ice shelf.
    !
    private
    type(cheb1d_scalar_field) :: thickness
      !! The thickness of the plume
    type(cheb1d_vector_field) :: velocity
      !! The velocity of the plume
    type(cheb1d_scalar_field) :: temperature
      !! The temperature of the plume
    type(cheb1d_scalar_field) :: salinity
      !! The salinity of the plume
    class(abstract_entrainment), allocatable :: entrainment_formulation
      !! An object which provides the parameterisation for entrainment
      !! of water into the plume.
    class(abstract_melt_relationship), allocatable :: melt_formulation
      !! An object which provides the parameterisation for melting,
      !! salt, and heat fluxes from the plume to the ice.
    class(ambient_conditions), allocatable :: ambient_conds
      !! An object specifying the temperature and salinity of the
      !! ambient ocean at its interface with the plume.
    class(plume_boundary), allocatable :: boundaries
      !! An object specifying the boundary conditions for the plume.
    real(r8)                  :: delta
      !! The dimensionless ratio $\delta \equiv \frac{D_0}{h_0}$
    real(r8)                  :: nu
      !! The dimensionless ratio $\nu \equiv \frac{\kappa_0}{x_0U_o}$
    real(r8)                  :: mu
      !! The dimensionless ratio $\mu \equiv \frac{\C_dx_0}{D_0}$
    real(r8)                  :: time
      !! The time at which the ice shelf is in this state
  contains
    procedure :: basal_melt => plume_melt
    procedure :: basal_drag_parameter => plume_drag_parameter
    procedure :: water_density => plume_water_density
    procedure :: residual => plume_residual
    procedure :: update => plume_update
    procedure :: set_time => plume_set_time
    procedure :: data_size => plume_data_size
    procedure :: state_vector => plume_state_vector
  end type plume

  abstract interface
    function scalar_func(location) result(thickness)
      !* Author: Chris MacMackin
      !  Date: April 2016
      !
      ! Abstract interface for function providing the initial values
      ! for the scalar properties of a [[plume(type)]] object when it
      ! is being instantiated.
      !
      import :: r8
      real(r8), dimension(:), intent(in) :: location
        !! The position $\vec{x}$ at which to compute the property
      real(r8)                           :: thickness
        !! The value of the scalar quantity at `location`
    end function scalar_func
      
    function velocity_func(location) result(velocity)
      !* Author: Chris MacMackin
      !  Date: April 2016
      !
      ! Abstract interface for function providing the [[plume(type)]] velocity
      ! when an object is being instantiated.
      !
      import :: r8
      real(r8), dimension(:), intent(in) :: location
        !! The position $\vec{x}$ at which to compute the thickness
      real(r8), dimension(2)             :: velocity
        !! The velocity vector of the water in the plume at `location`
    end function velocity_func
  end interface

contains

  function constructor(domain, resolution, thickness, velocity, temperature, &
                       salinity, entrainment_formulation, melt_formulation,  &
                       ambient_conds, boundaries, delta, nu, mu, sigma) result(this)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    ! 
    ! Instantiates a [[plume(type)]] object with initial coniditions
    ! provided by the arguments.At present only a 1D model is
    ! supported. If information is provided for higher dimensions then
    ! it will be ignored.
    !
    real(r8), dimension(:,:), intent(in) :: domain
      !! An array containing the upper and lower limits of the domain for
      !! the plume. The first index represents the dimension for which the
      !! boundaries apply. If the second index is 1 then it corresponds to
      !! the lower bound. If the second index is 2 then it corresponds to
      !! the upper bound.
    integer, dimension(:), intent(in)    :: resolution
      !! The number of data points in each dimension
    procedure(scalar_func)               :: thickness
      !! A function which calculates the initial value of the thickness of 
      !! the plume at a given location.
    procedure(velocity_func)             :: velocity
      !! A function which calculates the initial value of the velocity 
      !! (vector) of the water at a given location in a plume.
    procedure(scalar_func)               :: temperature
      !! A function which calculates the initial value of the temperature of 
      !! the plume at a given location.
    procedure(scalar_func)               :: salinity
      !! A function which calculates the initial value of the salinity of 
      !! the plume at a given location.
    class(abstract_entrainment), optional, intent(in) :: entrainment_formulation
      !! An object which calculates entrainment into the plume.
    class(abstract_melt_relationship), optional, intent(in) :: melt_formulation
      !! An object which calculates melting and the resulting thermal
      !! transfer into/out of the plume.
    class(ambient_conditions), optional, intent(in) :: ambient_conds
      !! An object specifying the salinity and temperature of the
      !! ambient ocean.
    class(plume_boundary), optional, intent(in) :: boundaries
      !! An object providing the boundary conditions for the plume.
    real(r8), optional, intent(in) :: delta
      !! The dimensionless ratio $\delta \equiv \frac{D_0}{h_0}$
    real(r8), optional, intent(in) :: nu
      !! The dimensionless ratio $\nu \equiv \frac{\kappa_0}{x_0U_o}$
    real(r8), optional, intent(in) :: mu
      !! The dimensionless ratio $\mu \equiv \frac{\C_dx_0}{D_0}$
    real(r8), optional, intent(in) :: sigma
      !! The dimensionless ratio 
      !! $\sigma \equiv \frac{U_0^2}{h_0g} = S_0\beta_S$
    type(plume) :: this
      !! A plume object with its domain and initial conditions set according
      !! to the arguments of the constructor function.
    
  end function constructor

  function plume_melt(this) result(melt)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Computes and returns the melt rate at the bottom of the ice
    ! shelf due to interaction with the plume.
    !
    class(plume), intent(in)         :: this
    class(scalar_field), allocatable :: melt
      !! The melt rate at the base of the ice shelf.
    allocate(melt, source=this%melt_formulation%melt_rate())
  end function plume_melt

  function plume_drag_parameter(this) result(drag)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Computes and returns a quantity which may be necessary to determine
    ! the frictional drag the plume exerts on the bottom of the ice
    ! shelf. The plume would actually tend to exert no drag on the bottom
    ! of the ice shelf, but this method is present so that there is a
    ! consistent interface with the [[ground(type)]] data type.
    !
    class(plume), intent(in)         :: this
    class(scalar_field), allocatable :: drag
      !! The melt rate at the base of the ice sheet.
  end function plume_drag_parameter

  function plume_water_density(this) result(density)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Computes and returns the density of the plume water beneath the ice
    ! shelf. The density of this water would vary depending on how much 
    ! saline ambient water has been entrained into the plume versus how
    ! much fresh water has been released due to melting. However, the
    ! Boussinesq approximation is used here and only a single reference 
    ! density is returned.
    !
    class(plume), intent(in) :: this
    real(r8)                 :: density
      !! The density of the water at the base of the ice sheet.
  end function plume_water_density
   
  function plume_residual(this, ice_thickness, ice_density, ice_temperature) &
                                                            result(residual)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Using the current state of the plume, this computes the residual
    ! of the system of equations which is used to describe the it.
    ! The residual takes the form of a 1D array, with each element 
    ! respresenting the residual for one of the equations in the system.
    !
    class(plume), intent(in)            :: this
    class(scalar_field), intent(in)     :: ice_thickness
      !! Thickness of the ice above the plume.
    real(r8), intent(in)                :: ice_density
      !! The density of the ice above the plume, assumed uniform.
    real(r8), intent(in)                :: ice_temperature
      !! The temperature of the ice above the plume, assumed uniform.
    real(r8), dimension(:), allocatable :: residual
      !! The residual of the system of equations describing the plume.
  end function plume_residual

  subroutine plume_update(this, state_vector)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Updates the state of the plume from its state vector. The state
    ! vector is a real array containing the value of each of the plume's
    ! properties at each of the locations on the grid used in discretization.
    !
    class(plume), intent(inout)        :: this
    real(r8), dimension(:), intent(in) :: state_vector
      !! A real array containing the data describing the state of the
      !! plume.
  end subroutine plume_update

  subroutine plume_set_time(this, time)
    !* Author: Christopher MacMackin
    !  Date: November 2016
    !
    ! Sets the time information held by the plume object. This is
    ! the time at which the plume is in its current state.
    !
    class(plume), intent(inout) :: this
    real(r8), intent(in)        :: time
      !! The time at which the plume is in the present state.
    this%time = time
  end subroutine plume_set_time

  pure function plume_data_size(this)
    !* Author: Christopher MacMackin
    !  Date: August 2016
    !
    ! Returns the number of elements in the plume's state vector.
    ! This is the size of the vector returned by [[plume(type):residual]]
    ! and [[plume(type):state_vector]] and taken as an argument by 
    ! [[plume(type):update]].
    !
    class(plume), intent(in) :: this
    integer                  :: plume_data_size
      !! The number of elements in the plume's state vector.
  end function plume_data_size

  pure function plume_state_vector(this) result(state_vector) 
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Returns the state vector for the current state of the plume. 
    ! This takes the form of a 1D array.
    !
    class(plume), intent(in)            :: this
    real(r8), dimension(:), allocatable :: state_vector
      !! The state vector describing the plume.
  end function plume_state_vector

end module plume_mod

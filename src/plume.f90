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
  use equation_of_state_mod, only: equation_of_state
  use jenkins1991_entrainment_mod, only: jenkins1991_entrainment
  use dallaston2015_melt_mod, only: dallaston2015_melt
  use uniform_ambient_mod, only: uniform_ambient_conditions
  use dallaston2015_plume_boundary_mod, only: dallaston2015_plume_boundary
  use linear_eos_mod, only: linear_eos
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
    class(equation_of_state), allocatable :: eos
      !! An object specifying the equation of state relating the plume
      !! water's density to its temperature and salinity.
    class(plume_boundary), allocatable :: boundaries
      !! An object specifying the boundary conditions for the plume.
    real(r8)                  :: delta
      !! The dimensionless ratio $\delta \equiv \frac{D_0}{h_0}$
    real(r8)                  :: nu
      !! The dimensionless ratio $\nu \equiv \frac{\kappa_0}{x_0U_o}$
    real(r8)                  :: mu
      !! The dimensionless ratio $\mu \equiv \frac{\C_dx_0}{D_0}$
    real(r8)                  :: epsilon
      !! A coefficient on the melt term in the continuity equation,
      !! which may be set to 0 to remove this term. This can be useful
      !! when testing with simple systems, such as that of Dallaston
      !! et al. (2015).
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
    pure function scalar_func(location) result(scalar)
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
      real(r8)                           :: scalar
        !! The value of the scalar quantity at `location`
    end function scalar_func
      
    pure function velocity_func(location) result(vector)
      !* Author: Chris MacMackin
      !  Date: April 2016
      !
      ! Abstract interface for function providing the [[plume(type)]] velocity
      ! when an object is being instantiated.
      !
      import :: r8
      real(r8), dimension(:), intent(in)  :: location
        !! The position $\vec{x}$ at which to compute the thickness
      real(r8), dimension(:), allocatable :: vector
        !! The velocity vector of the water in the plume at `location`
    end function velocity_func
  end interface

contains

  function constructor(domain, resolution, thickness, velocity, temperature, &
                       salinity, entrainment_formulation, melt_formulation,  &
                       ambient_conds, eos, boundaries, delta, nu, mu,        &
                       epsilon) result(this)
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
    class(abstract_entrainment), allocatable, optional, &
                           intent(inout) :: entrainment_formulation
      !! An object which calculates entrainment into the plume. Will
      !! be unallocated on exit. Defaults to that used by Jenkins
      !! (1991) with the coefficient $E_0 = 1$.
    class(abstract_melt_relationship), allocatable, optional, &
                           intent(inout) :: melt_formulation
      !! An object which calculates melting and the resulting thermal
      !! transfer into/out of the plume. Will be unallocated on
      !! exit. Defaults to that used by Dallaston et al. (2015),
      !! scaled to be consistent with the nondimensionalisation used
      !! here.
    class(ambient_conditions), allocatable, optional, &
                           intent(inout) :: ambient_conds
      !! An object specifying the salinity and temperature of the
      !! ambient ocean. Will be unallocated on exit. Defaults to
      !! uniform ambient salinity and temperature, both of which are
      !! set to 0 (as temperature and salinity are measured relative
      !! to some reference value).
    class(equation_of_state), allocatable, optional, &
                           intent(inout) :: eos
      !! An object specifying the equation of state for the water in
      !! the plume. Will be unallocated on exit. Defaults to
      !! linearised equation of state with no temperature dependence
      !! and a haline contraction coefficient of 1. The reference
      !! density is set to be 1 in the dimensionless units when
      !! salinity and temeprature are 0.
    class(plume_boundary), allocatable, optional, &
                           intent(inout) :: boundaries
      !! An object providing the boundary conditions for the
      !! plume. Will be unallocated on exit. Defaults to those used by
      !! Dallaston et al. (2015).
    real(r8), optional, intent(in)       :: delta
      !! The dimensionless ratio $\delta \equiv
      !! \frac{D_0}{h_0}$. Defaults to 0.036.
    real(r8), optional, intent(in)       :: nu
      !! The dimensionless ratio $\nu \equiv
      !! \frac{\kappa_0}{x_0U_o}$. Defaults to 0.
    real(r8), optional, intent(in)       :: mu
      !! The dimensionless ratio $\mu \equiv
      !! \frac{\C_dx_0}{D_0}$. Defaults to 0.
    real(r8), optional, intent(in)       :: epsilon
      !! A coefficient on the melting term in the continuity
      !! equation. Can be used to turn off this contribution, which
      !! may be useful depending on the scalings used. Defaults to
      !! 1.0.
    type(plume) :: this
      !! A plume object with its domain and initial conditions set according
      !! to the arguments of the constructor function.
    this%thickness = cheb1d_scalar_field(resolution(1),thickness,domain(1,1),domain(1,2))
    this%velocity = cheb1d_vector_field(resolution(1),velocity,domain(1,1),domain(1,2))
    this%thickness = cheb1d_scalar_field(resolution(1),temperature,domain(1,1),domain(1,2))
    this%thickness = cheb1d_scalar_field(resolution(1),salinity,domain(1,1),domain(1,2))
    if (present(entrainment_formulation)) then
      call move_alloc(entrainment_formulation, this%entrainment_formulation)
    else
      allocate(jenkins1991_entrainment :: this%entrainment_formulation)
    end if
    if (present(melt_formulation)) then
      call move_alloc(melt_formulation, this%melt_formulation)
    else
      allocate(dallaston2015_melt :: this%melt_formulation)
    end if
    if (present(ambient_conds)) then
      call move_alloc(ambient_conds, this%ambient_conds)
    else
      allocate(uniform_ambient_conditions :: this%ambient_conds)
    end if
    if (present(eos)) then
      call move_alloc(eos, this%eos)
    else
      allocate(linear_eos :: this%eos)
    end if
    if (present(boundaries)) then
      call move_alloc(boundaries, this%boundaries)
    else
      allocate(dallaston2015_plume_boundary :: this%boundaries)
    end if
    if (present(delta)) then
      this%delta = delta
    else
      this%delta = 0.036_r8
    end if
    if (present(nu)) then
      this%nu = nu
    else
      this%nu = 0.0_r8
    end if
    if (present(mu)) then
      this%mu = mu
    else
      this%mu = 0.0_r8
    end if
    if (present(epsilon)) then
      this%epsilon = epsilon
    else
      this%epsilon = 1.0_r8
    end if
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
    ! @NOTE Based on my approach to non-dimensionalisation, I'm pretty
    ! sure the density should always be 1, making this method
    ! unneccessary.
    !
    class(plume), intent(in) :: this
    real(r8)                 :: density
      !! The density of the water at the base of the ice sheet.
    density = 1.0_r8
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
    integer :: count, thick_count, vel_count, salt_count, temp_count
    
  end function plume_state_vector

end module plume_mod

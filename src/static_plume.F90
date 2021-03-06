!
!  plume.f90
!  This file is part of ISOFT.
!  
!  Copyright 2016 Chris MacMackin <cmacmackin@gmail.com>
!  
!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 3 of the License, or
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

#ifdef DEBUG
#define pure 
#define elemental 
#endif

module static_plume_mod
  !* Author: Christopher MacMackin
  !  Date: April 2016
  !  License: GPLv3
  !
  ! Provides a concrete implementation of the [[basal_surface(type)]]
  ! data type, representing a buoyant plume beneath an ice shelf. This
  ! implementation does not evolve the plume.
  !
  use iso_fortran_env, only: r8 => real64
  use basal_surface_mod, only: basal_surface, hdf_type_attr
  use factual_mod, only: scalar_field, vector_field, cheb1d_scalar_field, &
                         cheb1d_vector_field, uniform_scalar_field,       &
                         uniform_vector_field
  use ode_solvers_mod, only: quasilinear_solve
  use entrainment_mod, only: abstract_entrainment
  use melt_relationship_mod, only: abstract_melt_relationship
  use plume_boundary_mod, only: plume_boundary
  use upstream_plume_mod, only: upstream_plume_boundary
  use boundary_types_mod, only: free_boundary, dirichlet, neumann
  use ambient_mod, only: ambient_conditions
  use equation_of_state_mod, only: equation_of_state
  use jenkins1991_entrainment_mod, only: jenkins1991_entrainment
  use dallaston2015_melt_mod, only: dallaston2015_melt
  use uniform_ambient_mod, only: uniform_ambient_conditions
  use simple_plume_boundary_mod, only: simple_plume_boundary
  use pseudospectral_block_mod, only: pseudospec_block
  use linear_eos_mod, only: linear_eos
  use hdf5
  use h5lt
  use logger_mod, only: logger => master_logger
  use penf, only: str
  implicit none
  private

  character(len=9),  parameter, public :: hdf_type_name = 'plume'
  character(len=9),  parameter, public :: hdf_thickness = 'thickness'
  character(len=8),  parameter, public :: hdf_velocity = 'velocity'
  character(len=11), parameter, public :: hdf_temperature = 'temperature'
  character(len=8),  parameter, public :: hdf_salinity = 'salinity'
  character(len=5),  parameter, public :: hdf_delta = 'delta'
  character(len=2),  parameter, public :: hdf_nu = 'nu'
  character(len=2),  parameter, public :: hdf_mu = 'mu'
  character(len=5),  parameter, public :: hdf_r = 'r_val'
  character(len=3),  parameter, public :: hdf_phi = 'phi'

  type, extends(basal_surface), public :: static_plume
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! A concrete implementation of the [[basal_surface(type)]]
    ! abstract data type, representing the buoyant plume beneath an
    ! ice shelf, but unchanging in time. It keeps the values assigned
    ! at creation or with the [[static_plume:write_data]] method. It
    ! is useful if you want to uncouple the ice shelf from the plume.
    !
    private
    type(cheb1d_scalar_field) :: thickness
      !! The thickness of the plume
    type(cheb1d_vector_field) :: velocity
      !! The velocity of the plume
    type(cheb1d_vector_field) :: velocity_dx
      !! The derivative of the velocity field
    type(cheb1d_scalar_field) :: temperature
      !! The temperature of the plume
    type(cheb1d_scalar_field) :: temperature_dx
      !! The derivative of the temperature of the plume
    type(cheb1d_scalar_field) :: salinity
      !! The salinity of the plume
    type(cheb1d_scalar_field) :: salinity_dx
      !! The derivative of the salinity of the plume
    class(abstract_entrainment), allocatable :: entrainment_formulation
      !! An object which provides the parameterisation for entrainment
      !! of water into the plume.
    class(abstract_melt_relationship), allocatable :: melt_formulation
      !! An object which provides the parameterisation for melting,
      !! salt, and heat fluxes from the plume to the ice.
    class(ambient_conditions), allocatable :: ambient_conds
      !! An object specifying the temperature and salinity of the
      !! ambient ocean at its interface with the plume.
    class(equation_of_state), allocatable, public :: eos
      !! An object specifying the equation of state relating the plume
      !! water's density to its temperature and salinity.
    class(plume_boundary), allocatable :: boundaries
      !! An object specifying the boundary conditions for the plume.
    real(r8)                  :: delta
      !! The dimensionless ratio \(\delta \equiv \frac{D_0}{h_0}\)
    real(r8), public                  :: nu
      !! The dimensionless ratio \(\nu \equiv \frac{\kappa_0}{x_0U_o}\)
    real(r8)                  :: mu
      !! The dimensionless ratio \(\mu \equiv \frac{C_dx_0}{D_0}\)
    real(r8)                  :: r_val
      !! The dimensionless ratio of the ocean water density to the
      !! density of the overlying ice shelf.
    real(r8), public                  :: phi
      !! The inverse Rossby number, \(\Phi \equiv \frac{fx_0}{U_0}\)
    real(r8)                  :: time
      !! The time at which the ice shelf is in this state
    integer                   :: thickness_size
      !! The number of data values in the thickness field
    integer                   :: velocity_size
      !! The number of data values in the velocity field
    integer                   :: temperature_size
      !! The number of data values in the temperature field
    integer                   :: salinity_size
      !! the number of data values in the salinity field
    logical, dimension(7)     :: lower_bounds = .false.
      !! Which variables have boundary conditions at the grounding
      !! line.
    logical, dimension(7)     :: upper_bounds = .false.
      !! Which variables have boundary conditions at the calving
      !! front.
    type(pseudospec_block)    :: precond
      !! A pseudospectral differentiation block which can be used for
      !! preconditioning.
  contains
    procedure :: initialise => static_plume_initialise
    procedure :: basal_melt => static_plume_melt
    procedure :: basal_drag_parameter => static_plume_drag_parameter
    procedure :: water_density => static_plume_water_density
    procedure :: update => static_plume_update
    procedure :: data_size => static_plume_data_size
    procedure :: state_vector => static_plume_state_vector
    procedure :: read_data => static_plume_read_data
    procedure :: write_data => static_plume_write_data
    procedure :: solve => static_plume_solve
  end type static_plume


  abstract interface
     
#ifdef DEBUG
#undef pure
#undef elemental 
#endif

     pure function scalar_func(location) result(scalar)
      !* Author: Chris MacMackin
      !  Date: April 2016
      !
      ! Abstract interface for function providing the initial values
      ! for the scalar properties of a [[static_plume(type)]] object when it
      ! is being instantiated.
      !
      import :: r8
      real(r8), dimension(:), intent(in) :: location
        !! The position $\vec{x}$ at which to compute the property
      real(r8)                 :: scalar
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
    
#ifdef DEBUG
#define pure 
#define elemental 
#endif
    
  end interface

contains

  subroutine static_plume_initialise(this, domain, resolution, thickness, velocity,    &
                              temperature, salinity, entrainment_formulation,   &
                              melt_formulation, ambient_conds, eos, boundaries, &
                              delta, nu, mu, r_val, phi)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    ! 
    ! Instantiates a [[plume(type)]] object with initial coniditions
    ! provided by the arguments.At present only a 1D model is
    ! supported. If information is provided for higher dimensions then
    ! it will be ignored.
    !
    class(static_plume), intent(out)  :: this
      !! A plume object with its domain and initial conditions set according
      !! to the arguments of the constructor function.
    real(r8), dimension(:,:), intent(in) :: domain
      !! An array containing the upper and lower limits of the domain for
      !! the plume. The first index represents the dimension for which the
      !! boundaries apply. If the second index is 1 then it corresponds to
      !! the lower bound. If the second index is 2 then it corresponds to
      !! the upper bound.
    integer, dimension(:), intent(in)    :: resolution
      !! The number of data points in each dimension
    procedure(scalar_func)     :: thickness
      !! A function which calculates the initial value of the thickness of 
      !! the plume at a given location.
    procedure(velocity_func)   :: velocity
      !! A function which calculates the initial value of the velocity 
      !! (vector) of the water at a given location in a plume.
    procedure(scalar_func)     :: temperature
      !! A function which calculates the initial value of the temperature of 
      !! the plume at a given location.
    procedure(scalar_func)     :: salinity
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
      !! The dimensionless ratio \(\delta \equiv
      !! \frac{D_0}{h_0}\). Defaults to 0.036.
    real(r8), optional, intent(in)       :: nu
      !! The dimensionless ratio \(\nu \equiv
      !! \frac{\kappa_0}{x_0U_o}\). Defaults to 0.
    real(r8), optional, intent(in)       :: mu
      !! The dimensionless ratio \(\mu \equiv
      !! \frac{\C_dx_0}{D_0}\). Defaults to 0.
    real(r8), optional, intent(in)       :: r_val
      !! The dimensionless ratio of the water density to the ice shelf
      !! density, \( r = \rho_0/\rho_i. \) Defaults to 1.12.
    real(r8), optional, intent(in)       :: phi
      !! The inverse Rossby number, \(\Phi \equif
      !! \frac{fx_0}{U_0}\). Defaults to 0.

    integer :: i, btype_l, btype_u, bdepth_l, bdepth_u

    i = size(velocity([0._r8]))
    this%thickness = cheb1d_scalar_field(resolution(1),thickness,domain(1,1),domain(1,2))
    this%velocity = cheb1d_vector_field(resolution(1),velocity,domain(1,1),domain(1,2),i-1)
    this%temperature = cheb1d_scalar_field(resolution(1),temperature,domain(1,1),domain(1,2))
    this%salinity = cheb1d_scalar_field(resolution(1),salinity,domain(1,1),domain(1,2))
    this%thickness_size = this%thickness%raw_size()
    this%velocity_size = this%velocity%raw_size()
    this%temperature_size = this%temperature%raw_size()
    this%salinity_size = this%salinity%raw_size()
    this%velocity_dx = this%velocity%d_dx(1)
    this%salinity_dx = this%salinity%d_dx(1)
    this%temperature_dx = this%temperature%d_dx(1)
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
      allocate(simple_plume_boundary :: this%boundaries)
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
    if (present(r_val)) then
      this%r_val = r_val
    else
      this%r_val = 1.12_r8
    end if
    if (present(phi)) then
      this%phi = phi
    else
      this%phi = 0.0_r8
    end if
    this%time = 0.0_r8

#ifdef DEBUG
    call logger%debug('static_plume','Initialised new ice shelf object.')
#endif
  end subroutine static_plume_initialise


  function static_plume_melt(this) result(melt)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Computes and returns the melt rate at the bottom of the ice
    ! shelf due to interaction with the plume.
    !
    class(static_plume), intent(in)     :: this
    class(scalar_field), pointer :: melt
      !! The melt rate at the base of the ice shelf.
    melt => this%melt_formulation%melt_rate()
  end function static_plume_melt


  function static_plume_drag_parameter(this) result(drag)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Computes and returns a quantity which may be necessary to determine
    ! the frictional drag the plume exerts on the bottom of the ice
    ! shelf. The plume would actually tend to exert no drag on the bottom
    ! of the ice shelf, but this method is present so that there is a
    ! consistent interface with the [[ground(type)]] data type.
    !
    class(static_plume), intent(in)     :: this
    class(scalar_field), pointer :: drag
      !! The melt rate at the base of the ice sheet.
    type(uniform_scalar_field) :: dummy
    call dummy%allocate_scalar_field(drag)
    drag = uniform_scalar_field(0.0_r8)
    call drag%set_temp()
#ifdef DEBUG
    call logger%debug('static_plume%drag_parameter','Returned plume drag parameter.')
#endif
  end function static_plume_drag_parameter


  function static_plume_water_density(this) result(density)
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
    class(static_plume), intent(in) :: this
    real(r8)       :: density
      !! The density of the water at the base of the ice sheet.
    density = 1.0_r8
#ifdef DEBUG
    call logger%debug('static_plume%water_density','static_plume has average density '// &
                     trim(str(density))//'.')
#endif
  end function static_plume_water_density


  subroutine static_plume_update(this, state_vector, ice_thickness)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Updates the state of the plume from its state vector. The state
    ! vector is a real array containing the value of each of the plume's
    ! properties at each of the locations on the grid used in discretization.
    !
    class(static_plume), intent(inout)               :: this
    real(r8), dimension(:), intent(in)        :: state_vector
      !! A real array containing the data describing the state of the
      !! plume.
    class(scalar_field), optional, intent(in) :: ice_thickness
      !! The ice thickness which, if present, will be used to update
      !! the calculation of the melt rate.
    integer :: i
    !TODO: Add some assertion-like checks that the state vector is the right size
    call this%thickness%set_from_raw(state_vector(1:this%thickness_size))
    i = 1 + this%thickness_size
    call this%velocity%set_from_raw(state_vector(i:i + this%velocity_size - 1))
    i = i + this%velocity_size
    call this%velocity_dx%set_from_raw(state_vector(i:i + this%velocity_size - 1))
    i = i + this%velocity_size
    call this%temperature%set_from_raw(state_vector(i:i + this%temperature_size - 1))
    i = i + this%temperature_size
    call this%temperature_dx%set_from_raw(state_vector(i:i + this%temperature_size - 1))
    i = i + this%temperature_size
    call this%salinity%set_from_raw(state_vector(i:i + this%salinity_size - 1))
    i = i + this%salinity_size
    call this%salinity_dx%set_from_raw(state_vector(i:i + this%salinity_size - 1))
    if (present(ice_thickness)) then
      call this%melt_formulation%solve_for_melt(this%velocity, &
                                               -ice_thickness/this%r_val, &
                                               this%temperature, &
                                               this%salinity, &
                                               this%thickness, &
                                               this%time)      
    end if
#ifdef DEBUG
    call logger%debug('static_plume%update','Updated state of plume.')
#endif
  end subroutine static_plume_update


  function static_plume_data_size(this)
    !* Author: Christopher MacMackin
    !  Date: August 2016
    !
    ! Returns the number of elements in the plume's state vector.
    ! This is the size of the vector returned by
    ! [[static_plume(type):state_vector]] and taken as an argument by
    ! [[static_plume(type):update]].
    !
    class(static_plume), intent(in) :: this
    integer                  :: static_plume_data_size
      !! The number of elements in the plume's state vector.
    static_plume_data_size = this%thickness%raw_size() + this%velocity%raw_size() +      &
                      this%velocity_dx%raw_size() + this%temperature%raw_size() + &
                      this%temperature_dx%raw_size() + this%salinity%raw_size() + &
                      this%salinity_dx%raw_size()
#ifdef DEBUG
    call logger%debug('static_plume%data_size','static_plume shelf has '//     &
                      trim(str(static_plume_data_size))//' elements '// &
                      'in its state vector.')
#endif
  end function static_plume_data_size


  function static_plume_state_vector(this) result(state_vector) 
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Returns the state vector for the current state of the plume. 
    ! This takes the form of a 1D array.
    !
    class(static_plume), intent(in)  :: this
    real(r8), dimension(:), allocatable :: state_vector
      !! The state vector describing the plume.
    state_vector = [this%thickness%raw(), this%velocity%raw(),      &
                    this%velocity_dx%raw(), this%temperature%raw(), &
                    this%temperature_dx%raw(), this%salinity%raw(), &
                    this%salinity_dx%raw()]
#ifdef DEBUG
    call logger%debug('static_plume%state_vector','Returning state vector '// &
                      'for plume.')
#endif
  end function static_plume_state_vector


  subroutine static_plume_read_data(this,file_id,group_name,error)
    !* Author: Chris MacMackin
    !  Date: April 2017
    !
    ! Reads the state of the plume object from an HDF file in the
    ! specified group. This sets the thickness, velocity, temperature,
    ! salinity dataset, and parameter values.
    !
    class(static_plume), intent(inout)  :: this
    integer(hid_t), intent(in)   :: file_id
      !! The identifier for the HDF5 file/group in which this data is
      !! meant to be written.
    character(len=*), intent(in) :: group_name
      !! The name to give the group in the HDF5 file storing the
      !! ice shelf's data.
    integer, intent(out)         :: error
      !! Flag indicating whether routine ran without error. If no
      !! error occurs then has value 0.
    integer(hid_t) :: group_id
    integer :: ret_err
    real(r8), dimension(1) :: param
    character(len=50) :: base_type

    ret_err = 0
    call h5gopen_f(file_id, group_name, group_id, error)
    if (error /= 0) then
      call logger%error('static_plume%read_data','Could not open HDF group "'// &
                        group_name//'", so no IO performed.')
      return
    end if

    call h5ltget_attribute_string_f(file_id, group_name, hdf_type_attr, &
                                    base_type, error)
    if (trim(base_type) /= hdf_type_name) then
      call logger%error('static_plume%read_data','Trying to read data from '// &
                        'basal_surface of type other than plume.')
      error = -1
      return
    end if

   !call h5ltget_attribute_double_f(file_id, group_name, hdf_delta, &
   !                                param, error)
   !this%delta = param(1)
   !call h5ltget_attribute_double_f(file_id, group_name, hdf_nu, &
   !                                param, error)
   !this%nu = param(1)
   !call h5ltget_attribute_double_f(file_id, group_name, hdf_mu, &
   !                                param, error)
   !this%mu = param(1)
   !call h5ltget_attribute_double_f(file_id, group_name, hdf_r, &
   !                                param, error)
   !this%r_val = param(1)
   !call h5ltget_attribute_double_f(file_id, group_name, hdf_phi, &
   !                                param, error)
   !this%phi = param(1)
   !if (error /= 0) then
   !  call logger%warning('static_plume%read_data','Error code '//     &
   !                      trim(str(error))//' returned when '//  &
   !                      'reading attributes from HDF group '// &
   !                      group_name)
   !  ret_err = error
   !end if

    call this%thickness%read_hdf(group_id, hdf_thickness, error)
    if (error /= 0) then
      call logger%warning('static_plume%read_data','Error code '//     &
                          trim(str(error))//' returned when '// &
                          'reading plume thickness field '//   &
                          'from HDF file')
      if (ret_err == 0) ret_err = error
    end if

    call this%velocity%read_hdf(group_id, hdf_velocity, error)
    this%velocity_dx = this%velocity%d_dx(1)
    if (error /= 0) then
      call logger%warning('static_plume%read_data','Error code '//     &
                          trim(str(error))//' returned when '// &
                          'reading plume velocity field '//     &
                          'from HDF file')
      if (ret_err == 0) ret_err = error
    end if

    call this%temperature%read_hdf(group_id, hdf_temperature, error)
    this%temperature_dx = this%temperature%d_dx(1)
    if (error /= 0) then
      call logger%warning('static_plume%read_data','Error code '//     &
                          trim(str(error))//' returned when '// &
                          'reading plume temperature field '//  &
                          'from HDF file')
      if (ret_err == 0) ret_err = error
    end if

    call this%salinity%read_hdf(group_id, hdf_salinity, error)
    this%salinity_dx = this%salinity%d_dx(1)
    if (error /= 0) then
      call logger%warning('static_plume%read_data','Error code '//     &
                          trim(str(error))//' returned when '// &
                          'reading plume salinity field '//     &
                          'from HDF file')
      if (ret_err == 0) ret_err = error
    end if

    call h5gclose_f(group_id, error)
    if (error /= 0) then
      call logger%warning('static_plume%read_data','Error code '//     &
                          trim(str(error))//' returned when '// &
                          'closing HDF group '//group_name)
      if (ret_err == 0) ret_err = error
    end if
    error = ret_err
    call logger%trivia('static_plume%read_data','Read plume data from HDF group '// &
                       group_name)
  end subroutine static_plume_read_data


  subroutine static_plume_write_data(this,file_id,group_name,error)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Writes the state of the plume object to an HDF file in the
    ! specified group. This will consist of a thickness, a velocity, a
    ! temperature, and a salinity dataset.
    !
    class(static_plume), intent(in)     :: this
    integer(hid_t), intent(in)   :: file_id
      !! The identifier for the HDF5 file/group in which this data is
      !! meant to be written.
    character(len=*), intent(in) :: group_name
      !! The name to give the group in the HDF5 file storing the
      !! ice shelf's data.
    integer, intent(out)         :: error
      !! Flag indicating whether routine ran without error. If no
      !! error occurs then has value 0.
    integer(hid_t) :: group_id
    integer :: ret_err

    ret_err = 0
    call h5gcreate_f(file_id, group_name, group_id, error)
    if (error /= 0) then
      call logger%warning('static_plume%write_data','Error code '// &
                          trim(str(error))//' returned '//   &
                          'when creating HDF group "'//group_name//'"')
      call logger%error('static_plume%write_data','Data IO not performed for plume')
      return
    end if

    call h5ltset_attribute_string_f(file_id, group_name, hdf_type_attr, &
                                    hdf_type_name, error)
    call h5ltset_attribute_double_f(file_id, group_name, hdf_delta, &
                                    [this%delta], 1_size_t, error)
    call h5ltset_attribute_double_f(file_id, group_name, hdf_nu, &
                                    [this%nu], 1_size_t, error)
    call h5ltset_attribute_double_f(file_id, group_name, hdf_mu, &
                                    [this%mu], 1_size_t, error)
    call h5ltset_attribute_double_f(file_id, group_name, hdf_r, &
                                    [this%r_val], 1_size_t, error)
    call h5ltset_attribute_double_f(file_id, group_name, hdf_phi, &
                                    [this%phi], 1_size_t, error)
    if (error /= 0) then
      call logger%warning('static_plume%write_data','Error code '//    &
                          trim(str(error))//' returned when '// &
                          'writing attribute to HDF group '//   &
                          group_name)
      ret_err = error
    end if

    call this%thickness%write_hdf(group_id, hdf_thickness, error)
    if (error /= 0) then
      call logger%warning('static_plume%write_data','Error code '//    &
                          trim(str(error))//' returned when '// &
                          'writing plume thickness field to HDF')
      if (ret_err == 0) ret_err = error
    end if

    call this%velocity%write_hdf(group_id, hdf_velocity, error)
    if (error /= 0) then
      call logger%warning('static_plume%write_data','Error code '//    &
                          trim(str(error))//' returned when '// &
                          'writing plume velocity field to HDF')
      if (ret_err == 0) ret_err = error
    end if

    call this%temperature%write_hdf(group_id, hdf_temperature, error)
    if (error /= 0) then
      call logger%warning('static_plume%write_data','Error code '//       &
                          trim(str(error))//' returned when '//    &
                          'writing plume temperature field to HDF')
      if (ret_err == 0) ret_err = error
    end if

    call this%salinity%write_hdf(group_id, hdf_salinity, error)
    if (error /= 0) then
      call logger%warning('static_plume%write_data','Error code '//    &
                          trim(str(error))//' returned when '// &
                          'writing plume salinity field to HDF')
      if (ret_err == 0) ret_err = error
    end if

    call h5gclose_f(group_id, error)
    if (error /= 0) then
      call logger%warning('static_plume%write_data','Error code '//    &
                          trim(str(error))//' returned when '// &
                          'closing HDF group '//group_name)
      if (ret_err == 0) ret_err = error
    end if
    error = ret_err
#ifdef DEBUG
    call logger%debug('static_plume%write_data','Wrote plume data to HDF group '// &
                      group_name)
#endif
  end subroutine static_plume_write_data


  subroutine static_plume_solve(this, ice_thickness, ice_density, ice_temperature, &
                         time, success)
    !* Author: Chris MacMackin
    !  Date: March 2017
    !
    ! Would normally solve, but for this static implementation of the
    ! plume it does nothing.
    !
    class(static_plume), intent(inout)     :: this
    class(scalar_field), intent(in) :: ice_thickness
      !! Thickness of the ice above the basal surface
    real(r8), intent(in)            :: ice_density
      !! The density of the ice above the basal surface, assumed uniform
    real(r8), intent(in)            :: ice_temperature
      !! The temperature of the ice above the basal surface, assumed uniform
    real(r8), intent(in)            :: time
      !! The time to which the basal surface should be solved
    logical, intent(out)            :: success
      !! True if the solver is successful, false otherwise
    call ice_thickness%guard_temp()
    this%time = time
    call this%melt_formulation%solve_for_melt(this%velocity, &
                                             -ice_thickness/this%r_val, &
                                             this%temperature, &
                                             this%salinity, &
                                             this%thickness, &
                                             time)      
    success = .true.
    call ice_thickness%clean_temp()
  end subroutine static_plume_solve

end module static_plume_mod

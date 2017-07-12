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

module plume_mod
  !* Author: Christopher MacMackin
  !  Date: April 2016
  !  License: GPLv3
  !
  ! Provides a concrete implementation of the [[basal_surface(type)]] data type,
  ! representing a buoyant plume beneath an ice shelf.
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
  use boundary_types_mod, only: free_boundary, dirichlet, neumann
  use ambient_mod, only: ambient_conditions
  use equation_of_state_mod, only: equation_of_state
  use jenkins1991_entrainment_mod, only: jenkins1991_entrainment
  use dallaston2015_melt_mod, only: dallaston2015_melt
  use uniform_ambient_mod, only: uniform_ambient_conditions
  use simple_plume_boundary_mod, only: simple_plume_boundary
  use finite_difference_block_mod, only: fin_diff_block
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
      !! The dimensionless ratio $\delta \equiv \frac{D_0}{h_0}$
    real(r8), public                  :: nu
      !! The dimensionless ratio $\nu \equiv \frac{\kappa_0}{x_0U_o}$
    real(r8)                  :: mu
      !! The dimensionless ratio $\mu \equiv \frac{\C_dx_0}{D_0}$
    real(r8)                  :: r_val
      !! The dimensionless ratio of the ocean water density to the
      !! density of the overlying ice shelf.
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
    type(fin_diff_block), pointer :: precond_free_free => null()
      !! Preconditioner where upper and lower boundaries are both
      !! free.
    type(fin_diff_block), pointer :: precond_free_fixed => null()
      !! Preconditioner where lower boundary is free and upper
      !! boundary is fixed (Dirichlet).
    type(fin_diff_block), pointer :: precond_fixed_free => null()
      !! Preconditioner where lower boundary is fixed (Dirichlet) and
      !! upper boundary is free.
    type(fin_diff_block), pointer :: precond_fixed_fixed => null()
      !! Preconditioner where upper and lower boundaries are both
      !! fixed (Dirichlet).
    type(fin_diff_block), pointer :: thickness_precond => null()
      !! Preconditioner for plume thickness
    type(fin_diff_block), pointer :: velocity_precond => null()
      !! Preconditioner for plume velocity
    type(fin_diff_block), pointer :: velocity_dx_precond => null()
      !! Preconditioner for plume velocity derivative
    type(fin_diff_block), pointer :: temperature_precond => null()
      !! Preconditioner for plume temperature
    type(fin_diff_block), pointer :: temperature_dx_precond => null()
      !! Preconditioner for plume temperature derivative
    type(fin_diff_block), pointer :: salinity_precond => null()
      !! Preconditioner for plume salinity
    type(fin_diff_block), pointer :: salinity_dx_precond => null()
      !! Preconditioner for plume salinity derivative
  contains
    procedure :: initialise => plume_initialise
    procedure :: basal_melt => plume_melt
    procedure :: basal_drag_parameter => plume_drag_parameter
    procedure :: water_density => plume_water_density
    procedure :: update => plume_update
    procedure :: data_size => plume_data_size
    procedure :: state_vector => plume_state_vector
    procedure :: read_data => plume_read_data
    procedure :: write_data => plume_write_data
    procedure :: solve => plume_solve
    final :: plume_finalise
  end type plume


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
      ! for the scalar properties of a [[plume(type)]] object when it
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

  subroutine plume_initialise(this, domain, resolution, thickness, velocity,    &
                              temperature, salinity, entrainment_formulation,   &
                              melt_formulation, ambient_conds, eos, boundaries, &
                              delta, nu, mu, r_val)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    ! 
    ! Instantiates a [[plume(type)]] object with initial coniditions
    ! provided by the arguments.At present only a 1D model is
    ! supported. If information is provided for higher dimensions then
    ! it will be ignored.
    !
    class(plume), intent(out)  :: this
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

    integer :: btype_l, btype_u, bdepth_l, bdepth_u

    this%thickness = cheb1d_scalar_field(resolution(1),thickness,domain(1,1),domain(1,2))
    this%velocity = cheb1d_vector_field(resolution(1),velocity,domain(1,1),domain(1,2))
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
    this%time = 0.0_r8

    ! Initialise preconditioners
    call this%boundaries%thickness_bound_info(-1, btype_l, bdepth_l)
    call this%boundaries%thickness_bound_info(1, btype_u, bdepth_u)
#ifdef DEBUG
    if (abs(bdepth_l) > 1) then
      error stop ('Lower thickness boundary has depth greater than 1, '// &
                  'which is not supported by plume.')
    end if
    if (abs(bdepth_u) > 1) then
      error stop ('Upper thickness boundary has depth greater than 1, '// &
                  'which is not supported by plume.')
    end if
#endif
    select case(btype_l)
    case(free_boundary)
      select case(btype_u)
      case(free_boundary)
        call allocate_precond(this%precond_free_free)
        this%thickness_precond => this%precond_free_free
      case(dirichlet)
        call allocate_precond(this%precond_free_fixed, [1], [dirichlet])
        this%thickness_precond => this%precond_free_fixed
      case default
        error stop ('Only free, and Dirichlet boundary conditions '// &
                    'supported for plume thickness.')
      end select
    case(dirichlet)
      select case(btype_u)
      case(free_boundary)
        call allocate_precond(this%precond_fixed_free, resolution(1:1), &
                              [dirichlet])
        this%thickness_precond => this%precond_fixed_free
      case(dirichlet)
        call allocate_precond(this%precond_fixed_fixed, [1,resolution(1)], &
                              [dirichlet,dirichlet])
        this%thickness_precond => this%precond_fixed_fixed
      case default
        error stop ('Only free, Dirichlet, and Neumann boundary conditions '// &
                    'supported for plume.')
      end select
    case default
      error stop ('Only free, Dirichlet, and Neumann boundary conditions '// &
                  'supported for plume.')
    end select

    call this%boundaries%velocity_bound_info(-1, btype_l, bdepth_l)
    call this%boundaries%velocity_bound_info(1, btype_u, bdepth_u)
#ifdef DEBUG
    if (abs(bdepth_l) > 1) then
      error stop ('Lower velocity boundary has depth greater than 1, '// &
                  'which is not supported by plume.')
    end if
    if (abs(bdepth_u) > 1) then
      error stop ('Upper velocity boundary has depth greater than 1, '// &
                  'which is not supported by plume.')
    end if
#endif
    call set_preconditioners(btype_l, btype_u, this%velocity_precond, &
                             this%velocity_dx_precond)

    call this%boundaries%temperature_bound_info(-1, btype_l, bdepth_l)
    call this%boundaries%temperature_bound_info(1, btype_u, bdepth_u)
#ifdef DEBUG
    if (abs(bdepth_l) > 1) then
      error stop ('Lower temperature boundary has depth greater than 1, '// &
                  'which is not supported by plume.')
    end if
    if (abs(bdepth_u) > 1) then
      error stop ('Upper temperature boundary has depth greater than 1, '// &
                  'which is not supported by plume.')
    end if
#endif
    call set_preconditioners(btype_l, btype_u, this%temperature_precond, &
                             this%temperature_dx_precond)

    call this%boundaries%salinity_bound_info(-1, btype_l, bdepth_l)
    call this%boundaries%salinity_bound_info(1, btype_u, bdepth_u)
#ifdef DEBUG
    if (abs(bdepth_l) > 1) then
      error stop ('Lower salinity boundary has depth greater than 1, '// &
                  'which is not supported by plume.')
    end if
    if (abs(bdepth_u) > 1) then
      error stop ('Upper salinity boundary has depth greater than 1, '// &
                  'which is not supported by plume.')
    end if
#endif
    call set_preconditioners(btype_l, btype_u, this%salinity_precond, &
                             this%salinity_dx_precond)
    
#ifdef DEBUG
    call logger%debug('plume','Initialised new ice shelf object.')
#endif

  contains

    subroutine set_preconditioners(ltype, utype, pre, d_pre)
      integer, intent(in) :: ltype, utype
      type(fin_diff_block), pointer, intent(out)   :: pre, d_pre

      select case(ltype)
      case(free_boundary)
        select case(utype)
        case(free_boundary)
          call allocate_precond(this%precond_free_free)
          pre => this%precond_free_free
          d_pre => this%precond_free_free
        case(dirichlet)
          call allocate_precond(this%precond_free_free)
          call allocate_precond(this%precond_free_fixed, [1], [dirichlet])
          pre => this%precond_free_fixed
          d_pre => this%precond_free_free
        case(neumann)
          call allocate_precond(this%precond_free_free)
          call allocate_precond(this%precond_free_fixed, [1], [dirichlet])
          pre => this%precond_free_free
          d_pre => this%precond_free_fixed
        case default
          error stop ('Only free, Dirichlet, and Neumann boundary conditions '// &
                      'supported for plume.')
        end select
      case(dirichlet)
        select case(utype)
        case(free_boundary)
          call allocate_precond(this%precond_free_free)
          call allocate_precond(this%precond_fixed_free, resolution(1:1), &
                                [dirichlet])
          pre => this%precond_fixed_free
          d_pre => this%precond_free_free
        case(dirichlet)
          call allocate_precond(this%precond_free_free)
          call allocate_precond(this%precond_fixed_fixed, [1,resolution(1)], &
                                [dirichlet,dirichlet])
          pre => this%precond_fixed_fixed
          d_pre => this%precond_free_free
        case(neumann)
          call allocate_precond(this%precond_fixed_free, resolution(1:1), &
                                [dirichlet])
          call allocate_precond(this%precond_free_fixed, [1], [dirichlet])
          pre => this%precond_fixed_free
          d_pre => this%precond_free_fixed
        case default
          error stop ('Only free, Dirichlet, and Neumann boundary conditions '// &
                      'supported for plume.')
        end select
      case(neumann)
        select case(utype)
        case(free_boundary)
          call allocate_precond(this%precond_free_free)
          call allocate_precond(this%precond_fixed_free, resolution(1:1), &
                                [dirichlet])
          pre => this%precond_free_free
          d_pre => this%precond_fixed_free
        case(dirichlet)
          call allocate_precond(this%precond_fixed_free, resolution(1:1), &
                                [dirichlet])
          call allocate_precond(this%precond_free_fixed, [1], [dirichlet])
          pre => this%precond_free_fixed
          d_pre => this%precond_fixed_free
        case(neumann)
          call allocate_precond(this%precond_free_free)
          call allocate_precond(this%precond_fixed_fixed, [1,resolution(1)], &
                                [dirichlet,dirichlet])
          pre => this%precond_free_free
          d_pre => this%precond_fixed_fixed
        case default
          error stop ('Only free, Dirichlet, and Neumann boundary conditions '// &
                      'supported for plume.')
        end select
      case default
        error stop ('Only free, Dirichlet, and Neumann boundary conditions '// &
                    'supported for plume.')
      end select
    end subroutine set_preconditioners
    
    subroutine allocate_precond(precond, bound_locs, bound_types)
      type(fin_diff_block), pointer, intent(inout) :: precond
      integer, dimension(:), optional, intent(in)  :: bound_locs, bound_types
      if (.not. associated(precond)) then
        allocate(precond)
        precond = fin_diff_block(this%thickness, bound_locs, bound_types)
      end if
    end subroutine allocate_precond
  end subroutine plume_initialise


  function plume_melt(this) result(melt)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Computes and returns the melt rate at the bottom of the ice
    ! shelf due to interaction with the plume.
    !
    class(plume), intent(in)     :: this
    class(scalar_field), pointer :: melt
      !! The melt rate at the base of the ice shelf.
    melt => this%melt_formulation%melt_rate()
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
    class(plume), intent(in)     :: this
    class(scalar_field), pointer :: drag
      !! The melt rate at the base of the ice sheet.
    type(uniform_scalar_field) :: dummy
    call dummy%allocate_scalar_field(drag)
    drag = uniform_scalar_field(0.0_r8)
    call drag%set_temp()
#ifdef DEBUG
    call logger%debug('plume%drag_parameter','Returned plume drag parameter.')
#endif
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
    real(r8)       :: density
      !! The density of the water at the base of the ice sheet.
    density = 1.0_r8
#ifdef DEBUG
    call logger%debug('plume%water_density','Plume has average density '// &
                     trim(str(density))//'.')
#endif
  end function plume_water_density


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
#ifdef DEBUG
    call logger%debug('plume%update','Updated state of plume.')
#endif
  end subroutine plume_update


  function plume_data_size(this)
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
    plume_data_size = this%thickness%raw_size() + this%velocity%raw_size() +      &
                      this%velocity_dx%raw_size() + this%temperature%raw_size() + &
                      this%temperature_dx%raw_size() + this%salinity%raw_size() + &
                      this%salinity_dx%raw_size()
#ifdef DEBUG
    call logger%debug('plume%data_size','Plume shelf has '//     &
                      trim(str(plume_data_size))//' elements '// &
                      'in its state vector.')
#endif
  end function plume_data_size


  function plume_state_vector(this) result(state_vector) 
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Returns the state vector for the current state of the plume. 
    ! This takes the form of a 1D array.
    !
    class(plume), intent(in)  :: this
    real(r8), dimension(:), allocatable :: state_vector
      !! The state vector describing the plume.
    state_vector = [this%thickness%raw(), this%velocity%raw(),      &
                    this%velocity_dx%raw(), this%temperature%raw(), &
                    this%temperature_dx%raw(), this%salinity%raw(), &
                    this%salinity_dx%raw()]
#ifdef DEBUG
    call logger%debug('plume%state_vector','Returning state vector '// &
                      'for plume.')
#endif
  end function plume_state_vector


  subroutine plume_read_data(this,file_id,group_name,error)
    !* Author: Chris MacMackin
    !  Date: April 2017
    !
    ! Reads the state of the plume object from an HDF file in the
    ! specified group. This sets the thickness, velocity, temperature,
    ! salinity dataset, and parameter values.
    !
    class(plume), intent(inout)  :: this
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
      call logger%error('plume%read_data','Could not open HDF group "'// &
                        group_name//'", so no IO performed.')
      return
    end if

    call h5ltget_attribute_string_f(file_id, group_name, hdf_type_attr, &
                                    base_type, error)
    if (trim(base_type) /= hdf_type_name) then
      call logger%error('plume%read_data','Trying to read data from '// &
                        'basal_surface of type other than plume.')
      error = -1
      return
    end if

    call h5ltget_attribute_double_f(file_id, group_name, hdf_delta, &
                                    param, error)
    this%delta = param(1)
    call h5ltget_attribute_double_f(file_id, group_name, hdf_nu, &
                                    param, error)
    this%nu = param(1)
    call h5ltget_attribute_double_f(file_id, group_name, hdf_mu, &
                                    param, error)
    this%mu = param(1)
    call h5ltget_attribute_double_f(file_id, group_name, hdf_r, &
                                    param, error)
    this%r_val = param(1)
    if (error /= 0) then
      call logger%warning('plume%read_data','Error code '//     &
                          trim(str(error))//' returned when '//  &
                          'reading attributes from HDF group '// &
                          group_name)
      ret_err = error
    end if

    call this%thickness%read_hdf(group_id, hdf_thickness, error)
    if (error /= 0) then
      call logger%warning('plume%read_data','Error code '//     &
                          trim(str(error))//' returned when '// &
                          'reading plume thickness field '//   &
                          'from HDF file')
      if (ret_err == 0) ret_err = error
    end if

    call this%velocity%read_hdf(group_id, hdf_velocity, error)
    this%velocity_dx = this%velocity%d_dx(1)
    if (error /= 0) then
      call logger%warning('plume%read_data','Error code '//     &
                          trim(str(error))//' returned when '// &
                          'reading plume velocity field '//     &
                          'from HDF file')
      if (ret_err == 0) ret_err = error
    end if

    call this%temperature%read_hdf(group_id, hdf_temperature, error)
    this%temperature_dx = this%temperature%d_dx(1)
    if (error /= 0) then
      call logger%warning('plume%read_data','Error code '//     &
                          trim(str(error))//' returned when '// &
                          'reading plume temperature field '//  &
                          'from HDF file')
      if (ret_err == 0) ret_err = error
    end if

    call this%salinity%read_hdf(group_id, hdf_salinity, error)
    this%salinity_dx = this%salinity%d_dx(1)
    if (error /= 0) then
      call logger%warning('plume%read_data','Error code '//     &
                          trim(str(error))//' returned when '// &
                          'reading plume salinity field '//     &
                          'from HDF file')
      if (ret_err == 0) ret_err = error
    end if

    call h5gclose_f(group_id, error)
    if (error /= 0) then
      call logger%warning('plume%read_data','Error code '//     &
                          trim(str(error))//' returned when '// &
                          'closing HDF group '//group_name)
      if (ret_err == 0) ret_err = error
    end if
    error = ret_err
    call logger%trivia('plume%read_data','Read plume data from HDF group '// &
                       group_name)
  end subroutine plume_read_data


  subroutine plume_write_data(this,file_id,group_name,error)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Writes the state of the plume object to an HDF file in the
    ! specified group. This will consist of a thickness, a velocity, a
    ! temperature, and a salinity dataset.
    !
    class(plume), intent(in)     :: this
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
      call logger%warning('plume%write_data','Error code '// &
                          trim(str(error))//' returned '//   &
                          'when creating HDF group "'//group_name//'"')
      call logger%error('plume%write_data','Data IO not performed for plume')
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
    if (error /= 0) then
      call logger%warning('plume%write_data','Error code '//    &
                          trim(str(error))//' returned when '// &
                          'writing attribute to HDF group '//   &
                          group_name)
      ret_err = error
    end if

    call this%thickness%write_hdf(group_id, hdf_thickness, error)
    if (error /= 0) then
      call logger%warning('plume%write_data','Error code '//    &
                          trim(str(error))//' returned when '// &
                          'writing plume thickness field to HDF')
      if (ret_err == 0) ret_err = error
    end if

    call this%velocity%write_hdf(group_id, hdf_velocity, error)
    if (error /= 0) then
      call logger%warning('plume%write_data','Error code '//    &
                          trim(str(error))//' returned when '// &
                          'writing plume velocity field to HDF')
      if (ret_err == 0) ret_err = error
    end if

    call this%temperature%write_hdf(group_id, hdf_temperature, error)
    if (error /= 0) then
      call logger%warning('plume%write_data','Error code '//       &
                          trim(str(error))//' returned when '//    &
                          'writing plume temperature field to HDF')
      if (ret_err == 0) ret_err = error
    end if

    call this%salinity%write_hdf(group_id, hdf_salinity, error)
    if (error /= 0) then
      call logger%warning('plume%write_data','Error code '//    &
                          trim(str(error))//' returned when '// &
                          'writing plume salinity field to HDF')
      if (ret_err == 0) ret_err = error
    end if

    call h5gclose_f(group_id, error)
    if (error /= 0) then
      call logger%warning('plume%write_data','Error code '//    &
                          trim(str(error))//' returned when '// &
                          'closing HDF group '//group_name)
      if (ret_err == 0) ret_err = error
    end if
    error = ret_err
#ifdef DEBUG
    call logger%debug('plume%write_data','Wrote plume data to HDF group '// &
                      group_name)
#endif
  end subroutine plume_write_data


  subroutine plume_solve(this, ice_thickness, ice_density, ice_temperature, &
                         time, success)
    !* Author: Chris MacMackin
    !  Date: March 2017
    !
    ! Solves the state of the plume for the specified ice properties,
    ! at the specified time. This is done using the a
    ! quasilinearisation method and a GMRES iterative linear solver.
    !
    class(plume), intent(inout)     :: this
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
   
    real(r8), dimension(:), allocatable :: solution
    real(r8) :: residual
    integer, dimension(5) :: info
    integer :: flag
    class(scalar_field), pointer :: b

    call ice_thickness%guard_temp()
    b => -ice_thickness/this%r_val
    call b%guard_temp()
    this%time = time

    solution = this%state_vector()
#ifdef DEBUG
    call logger%debug('plume%solve','Calling QLM ODE solver')
#endif
    call quasilinear_solve(L, f, solution, 1, residual, flag, info,         &
                           1.e-9_r8*size(solution), precond=preconditioner, &
                           iter_max=100, krylov_dim=100, gmres_iter_max=5000)
    call this%update(solution)
#ifdef DEBUG
    call logger%debug('plume%solve','QLM solver required '//         &
                      trim(str(info(5)))//' nonlinear iterations '// &
                      'and '//trim(str(info(1)+info(2)))//           &
                      ' function calls.')
#endif

    select case(flag)
    case(0)
      call logger%trivia('plume%solver','Solved plume at time '//trim(str(time)))
      success = .true.
    case(1)
      call logger%warning('plume%solver','Plume solver stagnated with '// &
                       'residual of '//trim(str(residual)))
      success = .true.
    case(2)
      call logger%error('plume%solver','Reached maximum number of '// &
                        'iterations solving plume')
      success = .false.
    case(3)
      call logger%error('plume%solver','Plume solution began to diverge.')
      success = .false.
    case default
      call logger%error('plume%solve','QLM solver failed for plume with '// &
                        'error code '//trim(str(flag)))
      success = .false.
    end select
    
    call ice_thickness%clean_temp(); call b%clean_temp()

  contains

    function L(v)
      !! The linear differentiation operator
      real(r8), dimension(:), intent(in) :: v
        !! The state vector for the system of differential equations
      real(r8), dimension(size(v))       :: L
      
      integer :: st, en, btype_l, btype_u, bdepth_l, bdepth_u
      type(cheb1d_scalar_field) :: scalar_tmp
      type(cheb1d_vector_field) :: vector_tmp

      call this%update(v)
      
      ! Thickness
      call this%boundaries%thickness_bound_info(-1, btype_l, bdepth_l)
      call this%boundaries%thickness_bound_info(1, btype_u, bdepth_u)
      scalar_tmp = this%thickness%d_dx(1)
      select case(btype_l)
      case(dirichlet)
        call scalar_tmp%set_boundary(-1, bdepth_l,  &
                                     this%thickness%get_boundary(-1, bdepth_l))
      case(free_boundary)
        continue
      case default
        call logger%fatal('plume%solve','Plume can only handle Dirichlet '// &
                         'or free thickness boundaries.')
        error stop
      end select
      select case(btype_u)
      case(dirichlet)
        call scalar_tmp%set_boundary(1, bdepth_u, &
                                     this%thickness%get_boundary(1, bdepth_u))
      case(free_boundary)
        continue
      case default
        call logger%fatal('plume%solve','Plume can only handle Dirichlet '// &
                          'or free thickness boundaries.')
        error stop
      end select
      st = 1
      en = st + this%thickness_size - 1
      L(st:en) = scalar_tmp%raw()
      
      ! Velocity
      call this%boundaries%velocity_bound_info(-1, btype_l, bdepth_l)
      call this%boundaries%velocity_bound_info(1, btype_u, bdepth_u)
      vector_tmp = this%velocity%d_dx(1)
      if (btype_l == dirichlet) then
        call vector_tmp%set_boundary(-1, bdepth_l, &
                                     this%velocity%get_boundary(-1, bdepth_l))
      end if
      if (btype_u == dirichlet) then
        call vector_tmp%set_boundary(1, bdepth_u, &
                                     this%velocity%get_boundary(1, bdepth_u))
      end if
      st = en + 1
      en = st + this%velocity_size - 1
      L(st:en) = vector_tmp%raw()

      vector_tmp = this%velocity_dx%d_dx(1)
      if (btype_l == neumann) then
        call vector_tmp%set_boundary(-1, bdepth_l, &
                                     this%velocity_dx%get_boundary(-1, bdepth_l))
      else if (btype_l /= dirichlet .and. btype_l /= free_boundary) then
        call logger%fatal('plume%solve','Plume can only handle Dirichlet, '// &
                          'Neumann or free velocity boundaries.')
        error stop
      end if
      if (btype_u == neumann) then
        call vector_tmp%set_boundary(1, bdepth_u, &
                                     this%velocity_dx%get_boundary(1, bdepth_u))
      else if (btype_u /= dirichlet .and. btype_u /= free_boundary) then
        call logger%fatal('plume%solve','Plume can only handle Dirichlet, '// &
                          'Neumann or free velocity boundaries.')
        error stop
      end if
      st = en + 1
      en = st + this%velocity_size - 1
      L(st:en) = vector_tmp%raw()

      ! Temperature
      call this%boundaries%temperature_bound_info(-1, btype_l, bdepth_l)
      call this%boundaries%temperature_bound_info(1, btype_u, bdepth_u)
      scalar_tmp = this%temperature%d_dx(1)
      if (btype_l == dirichlet) then
        call scalar_tmp%set_boundary(-1, bdepth_l, &
                                     this%temperature%get_boundary(-1, bdepth_l))
      end if
      if (btype_u == dirichlet) then
        call scalar_tmp%set_boundary(1, bdepth_u, &
                                     this%temperature%get_boundary(1, bdepth_u))
      end if
      st = en + 1
      en = st + this%temperature_size - 1
      L(st:en) = scalar_tmp%raw()

      scalar_tmp = this%temperature_dx%d_dx(1)
      if (btype_l == neumann) then
        call scalar_tmp%set_boundary(-1, bdepth_l, &
                                     this%temperature_dx%get_boundary(-1, bdepth_l))
      else if (btype_l /= dirichlet .and. btype_l /= free_boundary) then
        call logger%fatal('plume%solve','Plume can only handle Dirichlet, '// &
                          'Neumann or free temperature boundaries.')
        error stop
      end if
      if (btype_u == neumann) then
        call scalar_tmp%set_boundary(1, bdepth_u, &
                                     this%temperature_dx%get_boundary(1, bdepth_u))
      else if (btype_u /= dirichlet .and. btype_u /= free_boundary) then
        call logger%fatal('plume%solve','Plume can only handle Dirichlet, '// &
                          'Neumann or free temperature boundaries.')
        error stop
      end if
      st = en + 1
      en = st + this%temperature_size - 1
      L(st:en) = scalar_tmp%raw()

      ! Salinity
      call this%boundaries%salinity_bound_info(-1, btype_l, bdepth_l)
      call this%boundaries%salinity_bound_info(1, btype_u, bdepth_u)
      scalar_tmp = this%salinity%d_dx(1)
      if (btype_l == dirichlet) then
        call scalar_tmp%set_boundary(-1, bdepth_l, &
                                     this%salinity%get_boundary(-1, bdepth_l))
      end if
      if (btype_u == dirichlet) then
        call scalar_tmp%set_boundary(1, bdepth_u, &
                                     this%salinity%get_boundary(1, bdepth_u))
      end if
      st = en + 1
      en = st + this%salinity_size - 1
      L(st:en) = scalar_tmp%raw()

      scalar_tmp = this%salinity_dx%d_dx(1)
      if (btype_l == neumann) then
        call scalar_tmp%set_boundary(-1, bdepth_l, &
                                     this%salinity_dx%get_boundary(-1, bdepth_l))
      else if (btype_l /= dirichlet .and. btype_l /= free_boundary) then
        call logger%fatal('plume%solve','Plume can only handle Dirichlet, '// &
                          'Neumann or free salinity boundaries.')
        error stop
      end if
      if (btype_u == neumann) then
        call scalar_tmp%set_boundary(1, bdepth_u, &
                                     this%salinity_dx%get_boundary(1, bdepth_u))
      else if (btype_u /= dirichlet .and. btype_u /= free_boundary) then
        call logger%fatal('plume%solve','Plume can only handle Dirichlet, '// &
                          'Neumann or free salinity boundaries.')
        error stop
      end if
      st = en + 1
      en = st + this%salinity_size - 1
      L(st:en) = scalar_tmp%raw()
    end function L

    subroutine non_diff_terms(D, Uvec, T, S, b, DU_x, DUU_x, DUT_x, DUS_x)
      !! Computes the values of \((DU)_x\), \((DU\vec{U})_x\),
      !! \((DUT)_x\), \((DUS)_x\), when diffusion is not
      !! included. This should be able to handle uniform field types,
      !! for use in an ODE solver when integrating near the boundary.
      class(scalar_field), intent(in)  :: D
        !! The plume thickness
      class(vector_field), intent(in)  :: Uvec
        !! The plume velocity
      class(scalar_field), intent(in)  :: T
        !! The plume temperature
      class(scalar_field), intent(in)  :: S
        !! The plume salinity
      class(scalar_field), intent(in)  :: b
        !! The debth of the base of the ice shelf
      class(scalar_field), intent(out) :: DU_x
        !! The derivative of the product DU
      class(vector_field), intent(out) :: DUU_x
        !! The derivative of the product DUU
      class(scalar_field), intent(out) :: DUT_x
        !! The derivative of the product DUT
      class(scalar_field), intent(out) :: DUS_x
        !! The derivative of the product DUS
     
      integer :: i, dims
      class(scalar_field), pointer :: m, rho, e, S_a, U, V, &
                                      T_a, rho_a, rho_x, Unorm
      class(scalar_field), allocatable, dimension(:) :: tmp
      call D%guard_temp(); call Uvec%guard_temp(); call T%guard_temp()
      call S%guard_temp(); call b%guard_temp()

      e => this%entrainment_formulation%entrainment_rate(Uvec, D, b, this%time)
      S_a => this%ambient_conds%ambient_salinity(b,this%time)
      T_a => this%ambient_conds%ambient_temperature(b,this%time)
      rho => this%eos%water_density(T, S)
      U => Uvec%component(1)
      V => Uvec%component(2)
      call e%guard_temp(); call S_a%guard_temp(); call T_a%guard_temp()
      call rho%guard_temp(); call U%guard_temp(); call V%guard_temp()
      rho_a => this%eos%water_density(T_a, S_a)
      call rho_a%guard_temp()

      call this%melt_formulation%solve_for_melt(Uvec, b, T, S, D, time)
      m => this%melt_formulation%melt_rate()
      call m%guard_temp()

      DU_x = e + m
      if (this%melt_formulation%has_heat_terms()) then
        DUT_x = e*T_a - this%melt_formulation%heat_equation_terms()
      else
        DUT_x = e*T_a
      end if
      if (this%melt_formulation%has_salt_terms()) then
        DUS_x = e*S_a - this%melt_formulation%salt_equation_terms()
      else
        DUS_x = e*S_a
      end if
      
      Unorm => Uvec%norm()
      rho_x => this%eos%water_density_derivative(T, T%d_dx(1), S, S%d_dx(1), 1)
      call Unorm%guard_temp()
      call rho_x%guard_temp()
      dims = Uvec%raw_size()/Uvec%elements()
      allocate(tmp(dims), mold=D)
      tmp(1) = (D*(rho_a - rho)*(b%d_dx(1) - 2*this%delta*DU_x/U) &
               + 0.5*this%delta*D**2*rho_x - this%mu*Unorm*U)/ &
               (1._r8 - this%delta*(rho_a - rho)/U**2)
      if (dims > 1) then
        tmp(2) = -this%mu*Unorm*V
      end if
      DUU_x = tmp

      call e%clean_temp(); call S_a%clean_temp(); call T_a%clean_temp()
      call rho%clean_temp(); call m%clean_temp(); call rho_a%clean_temp()
      call U%clean_temp(); call V%clean_temp(); call Unorm%clean_temp()
      call rho_x%clean_temp()
      call D%clean_temp(); call Uvec%clean_temp(); call T%clean_temp()
      call S%clean_temp(); call b%clean_temp()
    end subroutine non_diff_terms

    function f(v)
      !! The nonlinear operator
      real(r8), dimension(:,:), intent(in) :: v
        !! The state vector for the system of differential equations,
        !! and its derivatives. Column \(i\) represents the \(i-1\)
        !! derivative.
      real(r8), dimension(size(v,1)) :: f

      integer :: st, en, btype_l, btype_u, bdepth_l, bdepth_u
      type(cheb1d_scalar_field) :: scalar_tmp, D_x, D_nd, S_nd, T_nd
      type(cheb1d_vector_field) :: vector_tmp, U_nd
      class(scalar_field), pointer :: U, U_x

      call this%update(v(:,1))
      call this%boundaries%set_time(this%time)

      ! Use same or similar notation for variables as used in equations
      associate(D => this%thickness, Uvec => this%velocity,     &
                Uvec_x => this%velocity_dx, S => this%salinity, &
                S_x => this%salinity_dx, T => this%temperature, &
                T_x => this%temperature_dx, mf => this%melt_formulation, &
                h => ice_thickness, delta => this%delta, nu => this%nu,  &
                mu => this%mu, r => this%r_val, bounds => this%boundaries)

        call non_diff_terms(D, Uvec, T, S, b, D_nd, U_nd, T_nd, S_nd)

        ! FIXME: Alter this so that can take advantage of
        ! parameterisations returning uniform fields.
        U => this%velocity%component(1)
        U_x => this%velocity_dx%component(1)
        call U%guard_temp(); call U_x%guard_temp()

        ! Thickness
        call bounds%thickness_bound_info(-1, btype_l, bdepth_l)
        call bounds%thickness_bound_info(1, btype_u, bdepth_u)
        scalar_tmp = (D_nd - D*U_x)/U
        D_x = scalar_tmp
        select case(btype_l)
        case(dirichlet)
          call scalar_tmp%set_boundary(-1, bdepth_l, bounds%thickness_bound(-1))
        case(free_boundary)
          continue
        case default
          call logger%fatal('plume%solve','Plume can only handle Dirichlet '// &
                            'or free thickness boundaries.')
          error stop
        end select
        select case(btype_u)
        case(dirichlet)
          call scalar_tmp%set_boundary(1, bdepth_u, bounds%thickness_bound(1))
        case(free_boundary)
          continue
        case default
          call logger%fatal('plume%solve','Plume can only handle Dirichlet '// &
                            'or free thickness boundaries.')
          error stop
        end select
        st = 1
        en = st + this%thickness_size - 1
        f(st:en) = scalar_tmp%raw()
      
        ! Velocity
        call bounds%velocity_bound_info(-1, btype_l, bdepth_l)
        call bounds%velocity_bound_info(1, btype_u, bdepth_u)
        vector_tmp = Uvec_x
        if (btype_l == dirichlet) then
          call vector_tmp%set_boundary(-1, bdepth_l, bounds%velocity_bound(-1))
        end if
        if (btype_u == dirichlet) then
          call vector_tmp%set_boundary(1, bdepth_u, bounds%velocity_bound(1))
        end if
        st = en + 1
        en = st + this%velocity_size - 1
        f(st:en) = vector_tmp%raw()

        vector_tmp = D*U*Uvec_x !Needed due to compiler bug
        vector_tmp = (vector_tmp + D*U_x*Uvec + D_x*U*Uvec - &
                      U_nd - nu*D_x*Uvec_x)/(nu*D)
        if (btype_l == neumann) then
          call vector_tmp%set_boundary(-1, bdepth_l, bounds%velocity_bound(-1))
        else if (btype_l /= dirichlet .and. btype_l /= free_boundary) then
          call logger%fatal('plume%solve','Plume can only handle Dirichlet, '// &
                            'Neumann or free velocity boundaries.')
          error stop
        end if
        if (btype_u == neumann) then
          call vector_tmp%set_boundary(1, bdepth_u, bounds%velocity_bound(1))
        else if (btype_l /= dirichlet .and. btype_l /= free_boundary) then
          call logger%fatal('plume%solve','Plume can only handle Dirichlet, '// &
                            'Neumann or free velocity boundaries.')
          error stop
        end if
        st = en + 1
        en = st + this%velocity_size - 1
        f(st:en) = vector_tmp%raw()

        ! Temperature
        call bounds%temperature_bound_info(-1, btype_l, bdepth_l)
        call bounds%temperature_bound_info(1, btype_u, bdepth_u)
        scalar_tmp = T_x
        if (btype_l == dirichlet) then
          call scalar_tmp%set_boundary(-1, bdepth_l, bounds%temperature_bound(-1))
        end if
        if (btype_u == dirichlet) then
          call scalar_tmp%set_boundary(1, bdepth_u, bounds%temperature_bound(1))
        end if
        st = en + 1
        en = st + this%temperature_size - 1
        f(st:en) = scalar_tmp%raw()

        scalar_tmp = (D*U*T_x + D*U_x*T + D_x*U*T - T_nd - nu*D_x*T_x)/(nu*D)
        if (btype_l == neumann) then
          call scalar_tmp%set_boundary(-1, bdepth_l, bounds%temperature_bound(-1))
        else if (btype_l /= dirichlet .and. btype_l /= free_boundary) then
          call logger%fatal('plume%solve','Plume can only handle Dirichlet, '// &
                            'Neumann or free temperature boundaries.')
          error stop
        end if
        if (btype_u == neumann) then
          call scalar_tmp%set_boundary(1, bdepth_u, bounds%temperature_bound(1))
        else if (btype_l /= dirichlet .and. btype_l /= free_boundary) then
          call logger%fatal('plume%solve','Plume can only handle Dirichlet, '// &
                            'Neumann or free temperature boundaries.')
          error stop
        end if
        st = en + 1
        en = st + this%salinity_size - 1
        f(st:en) = scalar_tmp%raw()

        ! Salinity
        call bounds%salinity_bound_info(-1, btype_l, bdepth_l)
        call bounds%salinity_bound_info(1, btype_u, bdepth_u)
        scalar_tmp = S_x
        if (btype_l == dirichlet) then
          call scalar_tmp%set_boundary(-1, bdepth_l, bounds%salinity_bound(-1))
        end if
        if (btype_u == dirichlet) then
          call scalar_tmp%set_boundary(1, bdepth_u, bounds%salinity_bound(1))
        end if
        st = en + 1
        en = st + this%salinity_size - 1
        f(st:en) = scalar_tmp%raw()

        scalar_tmp = (D*U*S_x + D*U_x*S + D_x*U*S - S_nd - nu*D_x*S_x)/(nu*D)
        if (btype_l == neumann) then
          call scalar_tmp%set_boundary(-1, bdepth_l, bounds%salinity_bound(-1))
        else if (btype_l /= dirichlet .and. btype_l /= free_boundary) then
          call logger%fatal('plume%solve','Plume can only handle Dirichlet, '// &
                            'Neumann or free salinity boundaries.')
          error stop
        end if
        if (btype_u == neumann) then
          call scalar_tmp%set_boundary(1, bdepth_u, bounds%salinity_bound(1))
        else if (btype_l /= dirichlet .and. btype_l /= free_boundary) then
          call logger%fatal('plume%solve','Plume can only handle Dirichlet, '// &
                            'Neumann or free salinity boundaries.')
          error stop
        end if
        st = en + 1
        en = st + this%salinity_size - 1
        f(st:en) = scalar_tmp%raw()

        call U%clean_temp(); call U_x%clean_temp()
      end associate
    end function f

    function preconditioner(v, state, L_op, f_op, Lcur, fcur)
      !! The preconditioner, which approximates an inverse of `L`.
      real(r8), dimension(:), intent(in)   :: v
        !! The vector to be preconditioned.
      real(r8), dimension(:,:), intent(in) :: state
        !! The current state vector for the system of differential
        !! equations, and its derivatives. Column \(i\) represents the
        !! \(i-1\) derivative.
      procedure(L)                         :: L_op
        !! The linear, left-hand-side of the ODE being solved.
      procedure(f)                         :: f_op
        !! The nonlinear, right-hand-side of the ODE being solved.
      real(r8), dimension(:), intent(in)   :: Lcur
        !! The result of `L(u(:,1))`
      real(r8), dimension(:), intent(in)   :: fcur
        !! The result of `f(u)`
      real(r8), dimension(size(v)) :: preconditioner
        !! The result of applying the preconditioner.

      integer :: st, en, ust, uen
      real(r8) :: nu
      type(plume) :: v_plume
      type(cheb1d_scalar_field) :: scalar_tmp
      type(cheb1d_vector_field) :: vector_tmp
      class(scalar_field), pointer :: U, U_x

      v_plume%thickness_size = this%thickness_size
      call v_plume%thickness%assign_meta_data(this%thickness)
      v_plume%velocity_size = this%velocity_size
      call v_plume%velocity%assign_meta_data(this%velocity)
      call v_plume%velocity_dx%assign_meta_data(this%velocity_dx)
      v_plume%temperature_size = this%temperature_size
      call v_plume%temperature%assign_meta_data(this%temperature)
      call v_plume%temperature_dx%assign_meta_data(this%temperature_dx)
      v_plume%salinity_size = this%salinity_size
      call v_plume%salinity%assign_meta_data(this%salinity)
      call v_plume%salinity_dx%assign_meta_data(this%salinity_dx)
      call v_plume%update(v)

      call this%update(state(:,1))
      U => this%velocity%component(1)
      U_x => this%velocity_dx%component(1)
      call U%guard_temp(); call U_x%guard_temp()
      nu = this%nu

      v_plume%thickness = this%thickness_precond%solve_for(v_plume%thickness)!, U_x/U)
      st = 1
      en = st + this%thickness_size - 1
      preconditioner(st:en) = v_plume%thickness%raw()
  
      v_plume%velocity = this%velocity_precond%solve_for(v_plume%velocity)
      st = en + 1
      en = st + this%velocity_size - 1
      preconditioner(st:en) = v_plume%velocity%raw()
  
      ! Precondition the U_x term after have preconditioned values for S and T
      st = en + 1
      en = st + this%velocity_size - 1
      ust = st
      uen = en
  
      v_plume%temperature = this%temperature_precond%solve_for(v_plume%temperature)
      st = en + 1
      en = st + this%temperature_size - 1
      preconditioner(st:en) = v_plume%temperature%raw()
  
      ! Store diagonal offset in unused field
      v_plume%thickness = U/(-nu)
      v_plume%temperature_dx = &
           this%temperature_dx_precond%solve_for(v_plume%temperature_dx)!, &
!                                                v_plume%thickness)
      st = en + 1
      en = st + this%temperature_size - 1
      preconditioner(st:en) = v_plume%temperature_dx%raw()
  
      v_plume%salinity = this%salinity_precond%solve_for(v_plume%salinity)
      st = en + 1
      en = st + this%salinity_size - 1
      preconditioner(st:en) = v_plume%salinity%raw()
  
      v_plume%salinity_dx = &
           this%salinity_dx_precond%solve_for(v_plume%salinity_dx)!, &
!                                                v_plume%thickness)
      st = en + 1
      en = st + this%salinity_size - 1
      preconditioner(st:en) = v_plume%salinity_dx%raw()

      ! Precondition the U_x term using the values of S and T,
      ! allowing buoyance to be included.
      !v_plume%velocity = -uniform_vector_field([2._r8,1._r8])/nu * U
      !v_plume%velocity = (.grad. this%thickness)*this%eos%haline_contraction(this%temperature, this%salinity)/(this%r_val*nu)
      !v_plume%velocity = this%eos%haline_contraction(this%temperature, this%salinity)/(this%r_val*nu) * v_plume%velocity
      v_plume%velocity_dx = &
           this%velocity_dx_precond%solve_for(v_plume%velocity_dx)!, v_plume%velocity)
      preconditioner(ust:uen) = v_plume%velocity_dx%raw()

      call U%clean_temp(); call U_x%clean_temp()
    end function preconditioner
  end subroutine plume_solve


  elemental subroutine plume_finalise(this)
    type(plume), intent(inout) :: this
    if (associated(this%precond_free_free)) then
      deallocate(this%precond_free_free)
    end if
    if (associated(this%precond_fixed_free)) then
      deallocate(this%precond_fixed_free)
    end if
    if (associated(this%precond_free_fixed)) then
      deallocate(this%precond_free_fixed)
    end if
    if (associated(this%precond_fixed_fixed)) then
      deallocate(this%precond_fixed_fixed)
    end if
  end subroutine plume_finalise

end module plume_mod

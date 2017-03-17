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
  use factual_mod, only: scalar_field, cheb1d_scalar_field, cheb1d_vector_field, &
                         uniform_scalar_field
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
  implicit none
  private

  character(len=9),  parameter, public :: hdf_type_name = 'ice_shelf'
  character(len=9),  parameter, public :: hdf_thickness = 'thickness'
  character(len=8),  parameter, public :: hdf_velocity = 'velocity'
  character(len=11), parameter, public :: hdf_temperature = 'temperature'
  character(len=8),  parameter, public :: hdf_salinity = 'salinity'

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
    class(equation_of_state), allocatable :: eos
      !! An object specifying the equation of state relating the plume
      !! water's density to its temperature and salinity.
    class(plume_boundary), allocatable :: boundaries
      !! An object specifying the boundary conditions for the plume.
    real(r8)        :: delta
      !! The dimensionless ratio $\delta \equiv \frac{D_0}{h_0}$
    real(r8)        :: nu
      !! The dimensionless ratio $\nu \equiv \frac{\kappa_0}{x_0U_o}$
    real(r8)        :: mu
      !! The dimensionless ratio $\mu \equiv \frac{\C_dx_0}{D_0}$
    real(r8)        :: epsilon
      !! A coefficient on the melt term in the continuity equation,
      !! which may be set to 0 to remove this term. This can be useful
      !! when testing with simple systems, such as that of Dallaston
      !! et al. (2015).
    real(r8)        :: r_val
      !! The dimensionless ratio of the ocean water density to the
      !! density of the overlying ice shelf.
    real(r8)        :: time
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
    procedure :: residual => plume_residual
    procedure :: update => plume_update
    procedure :: set_time => plume_set_time
    procedure :: data_size => plume_data_size
    procedure :: state_vector => plume_state_vector
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
                              delta, nu, mu, epsilon, r_val)
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
    real(r8), optional, intent(in)       :: epsilon
      !! A coefficient on the melting term in the continuity
      !! equation. Can be used to turn off this contribution, which
      !! may be useful depending on the scalings used. Defaults to
      !! 1.0.
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
    if (present(epsilon)) then
      this%epsilon = epsilon
    else
      this%epsilon = 1.0_r8
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
    class(plume), intent(in)         :: this
    class(scalar_field), allocatable :: melt
      !! The melt rate at the base of the ice shelf.
    allocate(melt, source=this%melt_formulation%melt_rate())
    call melt%set_temp()
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
    allocate(uniform_scalar_field :: drag)
    drag = uniform_scalar_field(0.0_r8)
    call drag%set_temp()
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
  end function plume_water_density

   
  function plume_residual(this, ice_thickness, ice_density, ice_temperature) &
                                                            result(residual)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !  Deprecated: True
    !
    ! Using the current state of the plume, this computes the residual
    ! of the system of equations which is used to describe the it.
    ! The residual takes the form of a 1D array, with each element 
    ! respresenting the residual for one of the equations in the system.
    !
    ! @Note This is no longer used, as the plume is instead solved
    ! using a quaslinearisation method.
    !
    class(plume), intent(inout)         :: this
    class(scalar_field), intent(in)     :: ice_thickness
      !! Thickness of the ice above the plume.
    real(r8), intent(in)      :: ice_density
      !! The density of the ice above the plume, assumed uniform.
    real(r8), intent(in)      :: ice_temperature
      !! The temperature of the ice above the plume, assumed uniform.
    real(r8), dimension(:), allocatable :: residual
      !! The residual of the system of equations describing the plume.
    type(cheb1d_scalar_field) :: scalar_tmp
    type(cheb1d_vector_field) :: vector_tmp
    type(cheb1d_scalar_field) :: b, U, m, rho, e, S_a, T_a, rho_a
    integer :: start, finish
    integer, dimension(:), allocatable :: lower, upper
    
    call ice_thickness%guard_temp()

!    allocate(residual(this%data_size()))
!    residual = 0.0_r8
!    call this%melt_formulation%solve_for_melt(this%velocity, ice_thickness/this%r_val, &
!                                              this%temperature, this%salinity,        &
!                                              this%thickness, this%time)
!    start = 1
!
!    ! Use same or similar notation for variables as used in equations
!    associate(D => this%thickness, Uvec => this%velocity, &
!              S => this%salinity, &
!              T => this%temperature, &
!              mf => this%melt_formulation, h => ice_thickness, &              
!              delta => this%delta, nu => this%nu, mu => this%mu, &
!              epsilon => this%epsilon, r => this%r_val)
!
!      b = h/r
!      U = this%velocity%component(1)
!      m = this%melt_formulation%melt_rate()
!      rho = this%eos%water_density(this%temperature, this%salinity)
!      e = this%entrainment_formulation%entrainment_rate(this%velocity, this%thickness, &
!                                                        b,this%time)
!      call S_a%assign_meta_data(m)
!      call T_a%assign_meta_data(m)
!      call rho_a%assign_meta_data(m)
!      S_a = this%ambient_conds%ambient_salinity(b,this%time)
!      T_a = this%ambient_conds%ambient_temperature(b,this%time)
!      rho_a = this%eos%water_density(T_a, S_a)
!
!      ! Continuity equation
!      scalar_tmp = e + epsilon*m - .div.(D*Uvec)
!
!      lower = this%boundaries%thickness_lower_bound()
!      upper = this%boundaries%thickness_upper_bound()
!      finish = start + scalar_tmp%raw_size(lower,upper) - 1
!      residual(start:finish) = scalar_tmp%raw(lower,upper)
!      start = finish + 1
!
!      ! Momentum equation, x-component
!      scalar_tmp = .div.(D*.grad.U) ! Needed due to stupid
!                                    ! compiler issue. Works when
!                                    ! tested with trunk.
!      scalar_tmp = scalar_tmp*nu + &
!                   D*(rho_a - rho)*(this%delta*D%d_dx(1) - b%d_dx(1)) &
!                   - mu*Uvec%norm()*U + 0.5_r8*delta*(D**2)*rho%d_dx(1) &
!                   - .div.(D*Uvec*U)
!
!      lower = this%boundaries%velocity_lower_bound()
!      upper = this%boundaries%velocity_upper_bound()
!      finish = start + scalar_tmp%raw_size(lower,upper) - 1
!      residual(start:finish) = scalar_tmp%raw(lower,upper)
!      start = finish + 1
!
!      ! Salinity equation
!      scalar_tmp = .div.(D*.grad.S)
!      if (mf%has_salt_terms()) then
!        scalar_tmp = scalar_tmp*nu + e*S_a + mf%salt_equation_terms() &
!                   - .div.(D*Uvec*S)
!      else
!        scalar_tmp = scalar_tmp*nu + e*S_a - .div.(D*Uvec*S)
!      end if
!
!      lower = this%boundaries%salinity_lower_bound()
!      upper = this%boundaries%salinity_upper_bound()
!      finish = start + scalar_tmp%raw_size(lower,upper) - 1
!      residual(start:finish) = scalar_tmp%raw(lower,upper)
!      start = finish + 1
!
!      ! Temperature equation
!      scalar_tmp = D*T
!      scalar_tmp = .div.(D*.grad.T)
!      if (mf%has_heat_terms()) then
!        scalar_tmp = scalar_tmp*nu + e*T_a + mf%heat_equation_terms() &
!                   - .div.(D*Uvec*T)
!      else
!        scalar_tmp = scalar_tmp*nu + e*T_a - .div.(D*Uvec*T)
!      end if
!
!      lower = this%boundaries%temperature_lower_bound()
!      upper = this%boundaries%temperature_upper_bound()
!      finish = start + scalar_tmp%raw_size(lower,upper) - 1
!      residual(start:finish) = scalar_tmp%raw(lower,upper)
!
!      ! Boundary conditions
!      residual(finish+1:) = this%boundaries%boundary_residuals(D,Uvec,T,S,this%time)
!    end associate

    call ice_thickness%clean_temp()
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
    plume_data_size = this%thickness%raw_size() + this%velocity%raw_size() +      &
                      this%velocity_dx%raw_size() + this%temperature%raw_size() + &
                      this%temperature_dx%raw_size() + this%salinity%raw_size() + &
                      this%salinity_dx%raw_size()
  end function plume_data_size


  pure function plume_state_vector(this) result(state_vector) 
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
  end function plume_state_vector


  subroutine plume_write_data(this,file_id,group_name,error)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Writes the state of the plume object to an HDF file in the
    ! specified group. This will consist of a thickness, a velocity, a
    ! temperature, and a salinity dataset.
    !
    class(plume), intent(in) :: this
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
      write(*,*) 'WARNING: Error code',error,' returned when creating HDF '// &
                 'group', group_name
      write(*,*) '         Data IO not performed for ice shelf'
      return
    end if

    call h5ltset_attribute_string_f(file_id, group_name, hdf_type_attr, &
                                    hdf_type_name, error)
    if (error /= 0) then
      write(*,*) 'WARNING: Error code', error,' returned when writing '// &
                 'attribute to HDF group', group_name
      write(*,*) '         Output file will have missing information'
      ret_err = error
    end if

    call this%thickness%write_hdf(group_id, hdf_thickness, error)
    if (error /= 0) then
      write(*,*) 'WARNING: Error code',error,' returned when writing plume '// &
                 'thickness field to HDF'
      write(*,*) '         Data likely missing'
      if (ret_err == 0) ret_err = error
    end if

    call this%velocity%write_hdf(group_id, hdf_velocity, error)
    if (error /= 0) then
      write(*,*) 'WARNING: Error code',error,' returned when writing plume '// &
                 'velocity field to HDF'
      write(*,*) '         Data likely missing'
      if (ret_err == 0) ret_err = error
    end if

    call this%temperature%write_hdf(group_id, hdf_temperature, error)
    if (error /= 0) then
      write(*,*) 'WARNING: Error code',error,' returned when writing plume '// &
                 'temperature field to HDF'
      write(*,*) '         Data likely missing'
      if (ret_err == 0) ret_err = error
    end if

    call this%salinity%write_hdf(group_id, hdf_salinity, error)
    if (error /= 0) then
      write(*,*) 'WARNING: Error code',error,' returned when writing plume '// &
                 'salinity field to HDF'
      write(*,*) '         Data likely missing'
      if (ret_err == 0) ret_err = error
    end if

    call h5gclose_f(group_id, error)
    if (error /= 0) then
      write(*,*) 'WARNING: Error code',error,' returned when closing HDF '// &
                 'group', group_name
      write(*,*) '         Possible bad IO'
      if (ret_err == 0) ret_err = error
    end if
    error = ret_err
  end subroutine plume_write_data


  subroutine plume_solve(this, ice_thickness, ice_density, ice_temperature, time)
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
    real(r8), intent(in)  :: ice_density
      !! The density of the ice above the basal surface, assumed uniform
    real(r8), intent(in)  :: ice_temperature
      !! The temperature of the ice above the basal surface, assumed uniform
    real(r8), intent(in)  :: time
      !! The time to which the basal surface should be solved
    
    real(r8), dimension(:), allocatable :: solution
    real(r8) :: residual
    integer :: flag

    call ice_thickness%guard_temp()
    this%time = time

    solution = this%state_vector()
    call quasilinear_solve(L, f, solution, 1, residual, flag, &
                           precond=preconditioner)
    call this%update(solution)
    
    call ice_thickness%clean_temp()

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
        error stop ('Plume can only handle Dirichlet or free thickness boundaries.')
      end select
      select case(btype_u)
      case(dirichlet)
        call scalar_tmp%set_boundary(1, bdepth_u, &
                                     this%thickness%get_boundary(1, bdepth_u))
      case(free_boundary)
        continue
      case default
        error stop ('Plume can only handle Dirichlet or free thickness '// &
                    'boundaries.')
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
      st = en
      en = st + this%velocity_size - 1
      L(st:en) = vector_tmp%raw()

      vector_tmp = this%velocity_dx%d_dx(1)
      if (btype_l == neumann) then
        call vector_tmp%set_boundary(-1, bdepth_l, &
                                     this%velocity_dx%get_boundary(-1, bdepth_l))
      else if (btype_l /= dirichlet .and. btype_l /= free_boundary) then
        error stop ('Plume can only handle Dirichlet, Neumann or free '// &
                    'velocity boundaries.')
      end if
      if (btype_u == neumann) then
        call vector_tmp%set_boundary(1, bdepth_u, &
                                     this%velocity_dx%get_boundary(1, bdepth_u))
      else if (btype_l /= dirichlet .and. btype_l /= free_boundary) then
        error stop ('Plume can only handle Dirichlet, Neumann or free '// &
                    'velocity boundaries.')
      end if
      st = en
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
      st = en
      en = st + this%temperature_size - 1
      L(st:en) = scalar_tmp%raw()

      scalar_tmp = this%temperature_dx%d_dx(1)
      if (btype_l == neumann) then
        call scalar_tmp%set_boundary(-1, bdepth_l, &
                                     this%temperature_dx%get_boundary(-1, bdepth_l))
      else if (btype_l /= dirichlet .and. btype_l /= free_boundary) then
        error stop ('Plume can only handle Dirichlet, Neumann or free '// &
                    'temperature boundaries.')
      end if
      if (btype_u == neumann) then
        call scalar_tmp%set_boundary(1, bdepth_u, &
                                     this%temperature_dx%get_boundary(1, bdepth_u))
      else if (btype_l /= dirichlet .and. btype_l /= free_boundary) then
        error stop ('Plume can only handle Dirichlet, Neumann or free '// &
                    'temperature boundaries.')
      end if
      st = en
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
      st = en
      en = st + this%salinity_size - 1
      L(st:en) = scalar_tmp%raw()

      scalar_tmp = this%salinity_dx%d_dx(1)
      if (btype_l == neumann) then
        call scalar_tmp%set_boundary(-1, bdepth_l, &
                                     this%salinity_dx%get_boundary(-1, bdepth_l))
      else if (btype_l /= dirichlet .and. btype_l /= free_boundary) then
        error stop ('Plume can only handle Dirichlet, Neumann or free '// &
                    'salinity boundaries.')
      end if
      if (btype_u == neumann) then
        call scalar_tmp%set_boundary(1, bdepth_u, &
                                     this%salinity_dx%get_boundary(1, bdepth_u))
      else if (btype_l /= dirichlet .and. btype_l /= free_boundary) then
        error stop ('Plume can only handle Dirichlet, Neumann or free '// &
                    'salinity boundaries.')
      end if
      st = en
      en = st + this%salinity_size - 1
      L(st:en) = scalar_tmp%raw()
    end function L


    function f(v)
      !! The nonlinear operator
      real(r8), dimension(:,:), intent(in) :: v
        !! The state vector for the system of differential equations,
        !! and its derivatives. Column \(i\) represents the \(i-1\)
        !! derivative.
      real(r8), dimension(size(v,1)) :: f

      integer :: st, en, btype_l, btype_u, bdepth_l, bdepth_u
      type(cheb1d_scalar_field) :: scalar_tmp
      type(cheb1d_vector_field) :: vector_tmp
      type(cheb1d_scalar_field) :: b, D_x, U, U_x, m, rho, rho_x, e, S_a, &
                                   T_a, rho_a

      call this%update(v(:,1))

      ! Use same or similar notation for variables as used in equations
      associate(D => this%thickness, Uvec => this%velocity,      &
                Uvec_x => this%velocity_dx, S => this%salinity, &
                S_x => this%salinity_dx, T => this%temperature,  &
                T_x => this%temperature_dx, mf => this%melt_formulation, &
                h => ice_thickness, delta => this%delta, nu => this%nu,  &
                mu => this%mu, epsilon => this%epsilon, r => this%r_val, &
                bounds => this%boundaries)
  
        b = h/r
        call mf%solve_for_melt(Uvec, b, T, S, D, time)

        ! FIXME: Alter this so that can take advantage of
        ! parameterisations returning uniform fields.
        U = this%velocity%component(1)
        U_x = this%velocity_dx%component(1)
        m = mf%melt_rate()
        rho = this%eos%water_density(T, S)
        rho_x = this%eos%water_density_derivative(T,T_x,S,S_x,1)
        e = this%entrainment_formulation%entrainment_rate(Uvec, D, b, this%time)
        call S_a%assign_meta_data(m)
        call T_a%assign_meta_data(m)
        call rho_a%assign_meta_data(m)
        S_a = this%ambient_conds%ambient_salinity(b,this%time)
        T_a = this%ambient_conds%ambient_temperature(b,this%time)
        rho_a = this%eos%water_density(T_a, S_a)
        
        ! Thickness
        call bounds%thickness_bound_info(-1, btype_l, bdepth_l)
        call bounds%thickness_bound_info(1, btype_u, bdepth_u)
        D_x = (e + m - D*U_x)/U
        select case(btype_l)
        case(dirichlet)
          call D_x%set_boundary(-1, bdepth_l, bounds%thickness_bound(-1))
        case(free_boundary)
          continue
        case default
          error stop ('Plume can only handle Dirichlet or free thickness boundaries.')
        end select
        select case(btype_u)
        case(dirichlet)
          call D_x%set_boundary(1, bdepth_u, bounds%thickness_bound(1))
        case(free_boundary)
          continue
        case default
          error stop ('Plume can only handle Dirichlet or free thickness '// &
                      'boundaries.')
        end select
        st = 1
        en = st + this%thickness_size - 1
        f(st:en) = D_x%raw()
      
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
        st = en
        en = st + this%velocity_size - 1
        f(st:en) = vector_tmp%raw()

        vector_tmp = D*U*Uvec_x !Needed due to compiler bug
        vector_tmp = (vector_tmp + D*U_x*Uvec + D_x*U*Uvec + &
                      mu*Uvec%norm()*Uvec - nu*D_x*Uvec_x -  &
                      D*(rho_a - rho)*(.grad.(b - delta*D))  &
                      - 0.5_r8*delta*D**2*(.grad. rho))/(nu*D)
        if (btype_l == neumann) then
          call vector_tmp%set_boundary(-1, bdepth_l, bounds%velocity_bound(-1))
        else if (btype_l /= dirichlet .and. btype_l /= free_boundary) then
          error stop ('Plume can only handle Dirichlet, Neumann or free '// &
                      'velocity boundaries.')
        end if
        if (btype_u == neumann) then
          call vector_tmp%set_boundary(1, bdepth_u, bounds%velocity_bound(1))
        else if (btype_l /= dirichlet .and. btype_l /= free_boundary) then
          error stop ('Plume can only handle Dirichlet, Neumann or free '// &
                      'velocity boundaries.')
        end if
        st = en
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
        st = en
        en = st + this%temperature_size - 1
        f(st:en) = scalar_tmp%raw()

        if (mf%has_heat_terms()) then
          scalar_tmp = (D*U*T_x + D*U_x*T + D_x*U*T - e*T_a &
                       - nu*D_x*T_x)/(nu*D)
        else
          scalar_tmp = (D*U*T_x + D*U_x*T + D_x*U*T - e*T_a &
                       - nu*D_x*T_x - mf%heat_equation_terms())/(nu*D)
        end if
        if (btype_l == neumann) then
          call scalar_tmp%set_boundary(-1, bdepth_l, bounds%temperature_bound(-1))
        else if (btype_l /= dirichlet .and. btype_l /= free_boundary) then
          error stop ('Plume can only handle Dirichlet, Neumann or free '// &
                      'temperature boundaries.')
        end if
        if (btype_u == neumann) then
          call scalar_tmp%set_boundary(1, bdepth_u, bounds%temperature_bound(1))
        else if (btype_l /= dirichlet .and. btype_l /= free_boundary) then
          error stop ('Plume can only handle Dirichlet, Neumann or free '// &
                      'temperature boundaries.')
        end if
        st = en
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
        st = en
        en = st + this%salinity_size - 1
        f(st:en) = scalar_tmp%raw()

        if (mf%has_salt_terms()) then
          scalar_tmp = (D*U*S_x + D*U_x*S + D_x*U*S - e*S_a &
                       - nu*D_x*S_x)/(nu*D)
        else
          scalar_tmp = (D*U*S_x + D*U_x*S + D_x*U*S - e*S_a &
                       - nu*D_x*S_x - mf%salt_equation_terms())/(nu*D)
        end if
        if (btype_l == neumann) then
          call scalar_tmp%set_boundary(-1, bdepth_l, bounds%salinity_bound(-1))
        else if (btype_l /= dirichlet .and. btype_l /= free_boundary) then
          error stop ('Plume can only handle Dirichlet, Neumann or free '// &
                      'salinity boundaries.')
        end if
        if (btype_u == neumann) then
          call scalar_tmp%set_boundary(1, bdepth_u, bounds%salinity_bound(1))
        else if (btype_l /= dirichlet .and. btype_l /= free_boundary) then
          error stop ('Plume can only handle Dirichlet, Neumann or free '// &
                      'salinity boundaries.')
        end if
        st = en
        en = st + this%salinity_size - 1
        f(st:en) = scalar_tmp%raw()
      end associate
    end function f

    function preconditioner(v, u, L_op, f_op, Lcur, fcur)
      !! The preconditioner, which approximates an inverse of `L`.
      real(r8), dimension(:), intent(in)   :: v
        !! The vector to be preconditioned.
      real(r8), dimension(:,:), intent(in) :: u
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

      integer :: st, en
      type(plume) :: v_plume
      type(cheb1d_scalar_field) :: scalar_tmp
      type(cheb1d_vector_field) :: vector_tmp

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

      v_plume%thickness = this%thickness_precond%solve_for(v_plume%thickness)
      st = 1
      en = st + this%thickness_size - 1
      preconditioner(st:en) = v_plume%thickness%raw()

      v_plume%velocity = this%velocity_precond%solve_for(v_plume%velocity)
      st = en
      en = st + this%velocity_size - 1
      preconditioner(st:en) = v_plume%velocity%raw()

      v_plume%velocity_dx = &
           this%velocity_dx_precond%solve_for(v_plume%velocity_dx)
      st = en
      en = st + this%velocity_size - 1
      preconditioner(st:en) = v_plume%velocity_dx%raw()

      v_plume%temperature = this%temperature_precond%solve_for(v_plume%temperature)
      st = en
      en = st + this%temperature_size - 1
      preconditioner(st:en) = v_plume%temperature%raw()

      v_plume%temperature_dx = &
           this%temperature_dx_precond%solve_for(v_plume%temperature_dx)
      st = en
      en = st + this%temperature_size - 1
      preconditioner(st:en) = v_plume%temperature_dx%raw()

      v_plume%salinity = this%salinity_precond%solve_for(v_plume%salinity)
      st = en
      en = st + this%salinity_size - 1
      preconditioner(st:en) = v_plume%salinity%raw()

      v_plume%salinity_dx = &
           this%salinity_dx_precond%solve_for(v_plume%salinity_dx)
      st = en
      en = st + this%salinity_size - 1
      preconditioner(st:en) = v_plume%salinity_dx%raw()
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

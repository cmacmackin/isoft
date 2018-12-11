!
!  asymmetric_plume.f90
!  This file is part of ISOFT.
!  
!  Copyright 2017 Chris MacMackin <cmacmackin@gmail.com>
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

module asymmetric_plume_mod
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
  use upstream_plume_mod, only: upstream_plume_boundary
  use boundary_types_mod, only: free_boundary, dirichlet, neumann
  use ambient_mod, only: ambient_conditions
  use equation_of_state_mod, only: equation_of_state
  use jenkins1991_entrainment_mod, only: jenkins1991_entrainment
  use dallaston2015_melt_mod, only: dallaston2015_melt
  use uniform_ambient_mod, only: uniform_ambient_conditions
  use simple_plume_boundary_mod, only: simple_plume_boundary
  use pseudospectral_block_mod, only: pseudospec_block
  use coriolis_block_mod, only: coriolis_block
  use linear_eos_mod, only: linear_eos
  use ave_linear_eos_mod, only: ave_linear_eos
  use hdf5
  use h5lt
  use logger_mod, only: logger => master_logger
  use penf, only: str
  implicit none
  private

  character(len=10), parameter, public :: hdf_type_name = 'asym_plume'
  character(len=9),  parameter, public :: hdf_thickness = 'thickness'
  character(len=8),  parameter, public :: hdf_velocity = 'velocity'
  character(len=11), parameter, public :: hdf_temperature = 'temperature'
  character(len=8),  parameter, public :: hdf_salinity = 'salinity'
  character(len=5),  parameter, public :: hdf_delta = 'delta'
  character(len=2),  parameter, public :: hdf_nu = 'nu'
  character(len=2),  parameter, public :: hdf_mu = 'mu'
  character(len=5),  parameter, public :: hdf_r = 'r_val'
  character(len=3),  parameter, public :: hdf_phi = 'phi'

  type, public :: plume_shape
    !* Author: Christopher MacMackin
    !  Date: August 2018
    !
    ! A type containing the data necessary to specify the transverse
    ! shape of the plume. All variables in the plume are assumed to be
    ! seperable with the form \(D(x,y) = f_D2(y)\hat{D}(x)\), \(U(x,y)
    ! = f_U2(y)\hat{U}(x)\), etc. The magnitude of the velocity does
    ! not necessary take the form \(|\vec{U}| = \sqrt{U^2 + V^2}\) and
    ! is treated as an independent seperable variable \(|\vec{U}| =
    ! f_{|\vec{U}|}(y)\widehat{|\vec{U}|}(x)\). The transverse
    ! functions \(f_{D}(y)\), \(f_U2(y)\), etc. have all been
    ! normalised over the integration width \(y_1\) to \(y_2\).
    !
    ! @Note that a number of variables are definied using the
    ! averaging operator $$ \overline{A} = \frac{1}{y_2 - y_1}
    ! \int^{y_2}_{y_1} A(y) dy. $$
    ! 
    real(r8) :: f_D1 = 1.0_r8
      !! \(f_D2(y_1)\)
    real(r8) :: f_D2 = 1.0_r8
      !! \(f_D2(y_2)\)
    real(r8) :: f_U1 = 1.0_r8
      !! \(f_U2(y_1)\)
    real(r8) :: f_U2 = 1.0_r8
      !! \(f_U2(y_2)\)
    real(r8) :: f_V1 = 1.0_r8
      !! \(f_V2(y_1)\)
    real(r8) :: f_V2 = 1.0_r8
      !! \(f_V2(y_2)\)
    real(r8) :: f_S1 = 1.0_r8
      !! \(f_S2(y_1)\) 
    real(r8) :: f_S2 = 1.0_r8
      !! \(f_S2(y_2)\)
    real(r8) :: f_T1 = 1.0_r8
      !! \(f_T2(y_1)\)
    real(r8) :: f_T2 = 1.0_r8
      !! \(f_T2(y_2)\)
    real(r8) :: f_Up = 0.0_r8
      !! \(df_U2(y_2)/dy\)
    real(r8) :: f_Vp = 0.0_r8
      !! \(df_V2(y_2)/dy\)
    real(r8) :: f_Sp = 0.0_r8
      !! \(df_S2(y_2)/dy\)
    real(r8) :: f_Tp = 0.0_r8
      !! \(df_T2(y_2)/dy\)
    real(r8) :: a_DU = 1.0_r8
      !! \(\alpha_{DU} = \overline{f_D2 f_U2} \)
    real(r8) :: a_DV = 1.0_r8
      !! \(\alpha_{DV} = \overline{f_D2 f_V2} \)
    real(r8) :: a_DU2 = 1.0_r8
      !! \(\alpha_{DU^2} = \overline{f_D2 f_U2^2} \)
    real(r8) :: a_DUV = 1.0_r8
      !! \(\alpha_{DUV} = \overline{f_D2 f_U2 f_V2} \)
    real(r8) :: a_D2 = 1.0_r8
      !! \(\alpha_{D^2} = \overline{f_D2^2}\)
    real(r8) :: a_UabsU = 1.0_r8
      !! \(\alpha_{|\vec{U}|U} = \overline{f_{|\vec{U}|} f_U2\)
    real(r8) :: a_UabsV = 1.0_r8
      !! \(\alpha_{|\vec{U}|V} = \overline{f_{|\vec{U}|} f_V2\)
    real(r8) :: a_DUS = 1.0_r8
      !! \(\alpha_{DUS} = \overline{f_D2 f_U2 f_S2} \)
    real(r8) :: a_DUT = 1.0_r8
      !! \(\alpha_{DUT} = \overline{f_D2 f_U2 f_T2} \)
    real(r8) :: a_UabsT = 1.0_r8
      !! \(\alpha_{|\vec{U}|T} = \overline{f_{|\vec{U}|} f_T2\)
    real(r8) :: a_DS = 1.0_r8
      !! \(\alpha_{DS} = \overline{f_D2 f_S2} \)
    real(r8) :: a_DT = 1.0_r8
      !! \(\alpha_{DT} = \overline{f_D2 f_T2} \)
    real(r8) :: a_DS_t = 1.0_r8
      !! \(\tilde{\alpha}_{DS} = \overline{f_D2 f_S2}/\alpha_{D^2} \)
    real(r8) :: a_DT_t = 1.0_r8
      !! \(\tilde{\alpha}_{DT} = \overline{f_D2 f_T2}/\alpha_{D^2} \)
  end type plume_shape

  type, extends(basal_surface), public :: asym_plume
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! A concrete implementation of the [[basal_surface(type)]]
    ! abstract data type, representing the buoyant plume beneath an
    ! ice shelf, which has been horizontally integrated over width
    ! \(\Delta y = y_2 - y_1\) in the y-direction. Transverse
    ! variation is assumed to be seperable, with variables having the
    ! form \(D(x,y) = f_D2(y)\hat{D}(x)\), \(U(x,y) =
    ! f_U2(y)\hat{U}(x)\). The shape of the transverse functions are
    ! described using the [[plume_shape(type)]] type.
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
    real(r8), public          :: nu
      !! The dimensionless ratio \(\nu \equiv \frac{\kappa_0}{x_0U_o}\)
    real(r8)                  :: mu
      !! The dimensionless ratio \(\mu \equiv \frac{\C_dx_0}{D_0}\)
    real(r8)                  :: r_val
      !! The dimensionless ratio of the ocean water density to the
      !! density of the overlying ice shelf.
    real(r8), public          :: phi
      !! The inverse Rossby number, \(\Phi \equiv \frac{fx_0}{U_0}\)
    type(plume_shape), public         :: shape
      !! The coefficients describing the transverse shape of the
      !! plume.
    real(r8)                  :: dy
      !! The width \(\Delta y = y_2-y_1\) over which the plume has
      !! been horizontally averaged.
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
    integer                   :: vel_dims
      !! The number of vector dimensions for the velocity
    logical, dimension(7)     :: lower_bounds = .false.
      !! Which variables have boundary conditions at the grounding
      !! line.
    logical, dimension(7)     :: upper_bounds = .false.
      !! Which variables have boundary conditions at the calving
      !! front.
    type(pseudospec_block)    :: precond
      !! A pseudospectral differentiation block which can be used for
      !! preconditioning.
    type(coriolis_block)      :: vel_precond
      !! A representation of the operation on the plume velocity terms
      !! which can be used for preconditioning.
  contains
    procedure :: initialise => asym_plume_initialise
    procedure :: basal_melt => asym_plume_melt
    procedure :: basal_drag_parameter => asym_plume_drag_parameter
    procedure :: water_density => asym_plume_water_density
    procedure :: update => asym_plume_update
    procedure :: data_size => asym_plume_data_size
    procedure :: state_vector => asym_plume_state_vector
    procedure :: read_data => asym_plume_read_data
    procedure :: write_data => asym_plume_write_data
    procedure :: solve => asym_plume_solve
  end type asym_plume


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
      ! for the scalar properties of an [[asym_plume(type)]] object
      ! when it is being instantiated.
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
      ! Abstract interface for function providing the
      ! [[asym_plume(type)]] velocity when an object is being
      ! instantiated.
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

  subroutine asym_plume_initialise(this, domain, resolution,            &
                              thickness, velocity, temperature,         &
                              salinity, shape, entrainment_formulation, &
                              melt_formulation, ambient_conds, eos,     &
                              boundaries, delta, nu, mu, r_val, phi, dy)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    ! 
    ! Instantiates an [[asym_plume(type)]] object with initial
    ! coniditions provided by the arguments.At present only a 1D model
    ! is supported. If information is provided for higher dimensions
    ! then it will be ignored.
    !
    class(asym_plume), intent(out)  :: this
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
    type(plume_shape), optional, intent(in) :: shape
      !! An object with components describing the transverse profile
      !! of the plume. Default is a plume that is uniform in the
      !! transverse direction.
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
      !! \frac{C_dx_0}{D_0}\). Defaults to 0.
    real(r8), optional, intent(in)       :: r_val
      !! The dimensionless ratio of the water density to the ice shelf
      !! density, \( r = \rho_0/\rho_i. \) Defaults to 1.12.
    real(r8), optional, intent(in)       :: phi
      !! The inverse Rossby number, \(\Phi \equiv
      !! \frac{fx_0}{U_0}\). Defaults to 0.
    real(r8), optional, intent(in)       :: dy
      !! The width over which the plume has been averaged. Defualts to 1.

    integer :: i, btype_l, btype_u, bdepth_l, bdepth_u

    this%vel_dims = size(velocity([0._r8]))
    this%thickness = cheb1d_scalar_field(resolution(1),thickness,domain(1,1),domain(1,2))
    this%velocity = cheb1d_vector_field(resolution(1),velocity,domain(1,1),domain(1,2),this%vel_dims-1)
    this%temperature = cheb1d_scalar_field(resolution(1),temperature,domain(1,1),domain(1,2))
    this%salinity = cheb1d_scalar_field(resolution(1),salinity,domain(1,1),domain(1,2))
    this%thickness_size = this%thickness%raw_size()
    this%velocity_size = this%velocity%raw_size()
    this%temperature_size = this%temperature%raw_size()
    this%salinity_size = this%salinity%raw_size()
    this%velocity_dx = this%velocity%d_dx(1)
    this%salinity_dx = this%salinity%d_dx(1)
    this%temperature_dx = this%temperature%d_dx(1)
    if (present(shape)) this%shape = shape
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
    if (present(dy)) then
      this%dy = dy
    else
      this%dy = 1.0_r8
    end if
    this%time = 0.0_r8

    ! Initialise preconditioner
    this%precond = pseudospec_block(this%thickness)

    ! Store information on boundary conditions
    call this%boundaries%thickness_bound_info(-1, btype_l, bdepth_l)
    call this%boundaries%thickness_bound_info(1, btype_u, bdepth_u)
#ifdef DEBUG
    if (abs(bdepth_l) > 1) then
      error stop ('Lower thickness boundary has depth greater than 1, '// &
                  'which is not supported by asym_plume.')
    end if
    if (abs(bdepth_u) > 1) then
      error stop ('Upper thickness boundary has depth greater than 1, '// &
                  'which is not supported by asym_plume.')
    end if
#endif
    select case(btype_l)
    case(free_boundary)
    case(dirichlet)
      this%lower_bounds(1) = .true.
    case default
      error stop ('Only free, and Dirichlet boundary conditions '// &
                  'supported for asym_plume thickness.')
    end select

    select case(btype_u)
    case(free_boundary)
    case(dirichlet)
      this%upper_bounds(1) = .true.
    case default
      error stop ('Only free, Dirichlet, and Neumann boundary conditions '// &
                  'supported for asym_plume.')
    end select


    call this%boundaries%velocity_bound_info(-1, btype_l, bdepth_l)
    call this%boundaries%velocity_bound_info(1, btype_u, bdepth_u)
#ifdef DEBUG
    if (abs(bdepth_l) > 1) then
      error stop ('Lower velocity boundary has depth greater than 1, '// &
                  'which is not supported by asym_plume.')
    end if
    if (abs(bdepth_u) > 1) then
      error stop ('Upper velocity boundary has depth greater than 1, '// &
                  'which is not supported by asym_plume.')
    end if
#endif
    call set_preconditioners(btype_l, btype_u, 2)

    call this%boundaries%temperature_bound_info(-1, btype_l, bdepth_l)
    call this%boundaries%temperature_bound_info(1, btype_u, bdepth_u)
#ifdef DEBUG
    if (abs(bdepth_l) > 1) then
      error stop ('Lower temperature boundary has depth greater than 1, '// &
                  'which is not supported by asym_plume.')
    end if
    if (abs(bdepth_u) > 1) then
      error stop ('Upper temperature boundary has depth greater than 1, '// &
                  'which is not supported by asym_plume.')
    end if
#endif
    call set_preconditioners(btype_l, btype_u, 4)

    call this%boundaries%salinity_bound_info(-1, btype_l, bdepth_l)
    call this%boundaries%salinity_bound_info(1, btype_u, bdepth_u)
#ifdef DEBUG
    if (abs(bdepth_l) > 1) then
      error stop ('Lower salinity boundary has depth greater than 1, '// &
                  'which is not supported by asym_plume.')
    end if
    if (abs(bdepth_u) > 1) then
      error stop ('Upper salinity boundary has depth greater than 1, '// &
                  'which is not supported by asym_plume.')
    end if
#endif
    call set_preconditioners(btype_l, btype_u, 6)

    if (this%phi /= 0._r8) then
      associate (velbound => btype_l, dvelbound => btype_u)
        if (this%lower_bounds(2)) then
          velbound = -1
        else if (this%upper_bounds(2)) then
          velbound = 1
        end if
        if (this%lower_bounds(3)) then
          dvelbound = -1
        else if (this%upper_bounds(3)) then
          dvelbound = 1
        end if

        this%vel_precond = coriolis_block(this%phi, this%nu, velbound, &
                                          dvelbound, 1, this%thickness)
      end associate
    end if
    
#ifdef DEBUG
    call logger%debug('asym_plume','Initialised new ice shelf object.')
#endif

  contains

    subroutine set_preconditioners(ltype, utype, comp_id)
      integer, intent(in) :: ltype, utype, comp_id

      select case(ltype)
      case(free_boundary)
      case(dirichlet)
        this%lower_bounds(comp_id) = .true.
      case(neumann)
        this%lower_bounds(comp_id + 1) = .true.
      case default
        error stop ('Only free, Dirichlet, and Neumann boundary conditions '// &
                    'supported for asym_plume.')
      end select

      select case(utype)
      case(free_boundary)
      case(dirichlet)
        this%upper_bounds(comp_id) = .true.
      case(neumann)
        this%upper_bounds(comp_id + 1) = .true.
      case default
        error stop ('Only free, Dirichlet, and Neumann boundary conditions '// &
                    'supported for asym_plume.')
      end select
    end subroutine set_preconditioners
    
  end subroutine asym_plume_initialise


  function asym_plume_melt(this) result(melt)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Computes and returns the melt rate at the bottom of the ice
    ! shelf due to interaction with the plume.
    !
    class(asym_plume), intent(in)     :: this
    class(scalar_field), pointer :: melt
      !! The melt rate at the base of the ice shelf.
    melt => this%melt_formulation%melt_rate()
#ifdef DEBUG
    call logger%debug('asym_plume%melt','Returned melt rate.')
#endif
  end function asym_plume_melt


  function asym_plume_drag_parameter(this) result(drag)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Computes and returns a quantity which may be necessary to determine
    ! the frictional drag the plume exerts on the bottom of the ice
    ! shelf. The plume would actually tend to exert no drag on the bottom
    ! of the ice shelf, but this method is present so that there is a
    ! consistent interface with the [[ground(type)]] data type.
    !
    class(asym_plume), intent(in)     :: this
    class(scalar_field), pointer :: drag
      !! The melt rate at the base of the ice sheet.
    type(uniform_scalar_field) :: dummy
    call dummy%allocate_scalar_field(drag)
    drag = uniform_scalar_field(0.0_r8)
    call drag%set_temp()
#ifdef DEBUG
    call logger%debug('asym_plume%drag_parameter','Returned plume drag parameter.')
#endif
  end function asym_plume_drag_parameter


  function asym_plume_water_density(this) result(density)
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
    class(asym_plume), intent(in) :: this
    real(r8)       :: density
      !! The density of the water at the base of the ice sheet.
    density = 1.0_r8
#ifdef DEBUG
    call logger%debug('asym_plume%water_density','Plume has average density '// &
                     trim(str(density))//'.')
#endif
  end function asym_plume_water_density


  subroutine asym_plume_update(this, state_vector, ice_thickness)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Updates the state of the plume from its state vector. The state
    ! vector is a real array containing the value of each of the plume's
    ! properties at each of the locations on the grid used in discretization.
    !
    class(asym_plume), intent(inout)               :: this
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
    call logger%debug('asym_plume%update','Updated state of plume.')
#endif
  end subroutine asym_plume_update


  function asym_plume_data_size(this)
    !* Author: Christopher MacMackin
    !  Date: August 2016
    !
    ! Returns the number of elements in the plume's state vector.
    ! This is the size of the vector returned by
    ! [[asym_plume(type):state_vector]] and taken as an argument by
    ! [[plume(type):update]].
    !
    class(asym_plume), intent(in) :: this
    integer                       :: asym_plume_data_size
      !! The number of elements in the plume's state vector.
    asym_plume_data_size = this%thickness%raw_size() + this%velocity%raw_size() +      &
                           this%velocity_dx%raw_size() + this%temperature%raw_size() + &
                           this%temperature_dx%raw_size() + this%salinity%raw_size() + &
                           this%salinity_dx%raw_size()
#ifdef DEBUG
    call logger%debug('asym_plume%data_size','Plume shelf has '//     &
                      trim(str(asym_plume_data_size))//' elements '// &
                      'in its state vector.')
#endif
  end function asym_plume_data_size


  function asym_plume_state_vector(this) result(state_vector) 
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Returns the state vector for the current state of the plume. 
    ! This takes the form of a 1D array.
    !
    class(asym_plume), intent(in)  :: this
    real(r8), dimension(:), allocatable :: state_vector
      !! The state vector describing the plume.
    state_vector = [this%thickness%raw(), this%velocity%raw(),      &
                    this%velocity_dx%raw(), this%temperature%raw(), &
                    this%temperature_dx%raw(), this%salinity%raw(), &
                    this%salinity_dx%raw()]
#ifdef DEBUG
    call logger%debug('asym_plume%state_vector','Returning state vector '// &
                      'for plume.')
#endif
  end function asym_plume_state_vector


  subroutine asym_plume_read_data(this,file_id,group_name,error)
    !* Author: Chris MacMackin
    !  Date: April 2017
    !
    ! Reads the state of the plume object from an HDF file in the
    ! specified group. This sets the thickness, velocity, temperature,
    ! salinity dataset, and parameter values.
    !
    class(asym_plume), intent(inout)  :: this
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
      call logger%error('asym_plume%read_data','Could not open HDF group "'// &
                        group_name//'", so no IO performed.')
      return
    end if

    call h5ltget_attribute_string_f(file_id, group_name, hdf_type_attr, &
                                    base_type, error)
    if (trim(base_type) /= hdf_type_name) then
      call logger%error('asym_plume%read_data','Trying to read data from '// &
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
   !  call logger%warning('plume%read_data','Error code '//     &
   !                      trim(str(error))//' returned when '//  &
   !                      'reading attributes from HDF group '// &
   !                      group_name)
   !  ret_err = error
   !end if

    call this%thickness%read_hdf(group_id, hdf_thickness, error)
    if (error /= 0) then
      call logger%warning('asym_plume%read_data','Error code '//     &
                          trim(str(error))//' returned when '// &
                          'reading plume thickness field '//   &
                          'from HDF file')
      if (ret_err == 0) ret_err = error
    end if

    call this%velocity%read_hdf(group_id, hdf_velocity, error)
    this%velocity_dx = this%velocity%d_dx(1)
    if (error /= 0) then
      call logger%warning('asym_plume%read_data','Error code '//     &
                          trim(str(error))//' returned when '// &
                          'reading plume velocity field '//     &
                          'from HDF file')
      if (ret_err == 0) ret_err = error
    end if

    call this%temperature%read_hdf(group_id, hdf_temperature, error)
    this%temperature_dx = this%temperature%d_dx(1)
    if (error /= 0) then
      call logger%warning('asym_plume%read_data','Error code '//     &
                          trim(str(error))//' returned when '// &
                          'reading plume temperature field '//  &
                          'from HDF file')
      if (ret_err == 0) ret_err = error
    end if

    call this%salinity%read_hdf(group_id, hdf_salinity, error)
    this%salinity_dx = this%salinity%d_dx(1)
    if (error /= 0) then
      call logger%warning('asym_plume%read_data','Error code '//     &
                          trim(str(error))//' returned when '// &
                          'reading plume salinity field '//     &
                          'from HDF file')
      if (ret_err == 0) ret_err = error
    end if

    call h5gclose_f(group_id, error)
    if (error /= 0) then
      call logger%warning('asym_plume%read_data','Error code '//     &
                          trim(str(error))//' returned when '// &
                          'closing HDF group '//group_name)
      if (ret_err == 0) ret_err = error
    end if
    error = ret_err
    call logger%trivia('asym_plume%read_data','Read plume data from HDF group '// &
                       group_name)
  end subroutine asym_plume_read_data


  subroutine asym_plume_write_data(this,file_id,group_name,error)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Writes the state of the plume object to an HDF file in the
    ! specified group. This will consist of a thickness, a velocity, a
    ! temperature, and a salinity dataset.
    !
    class(asym_plume), intent(in)     :: this
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
      call logger%warning('asym_plume%write_data','Error code '// &
                          trim(str(error))//' returned '//   &
                          'when creating HDF group "'//group_name//'"')
      call logger%error('asym_plume%write_data','Data IO not performed for plume')
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
      call logger%warning('asym_plume%write_data','Error code '//    &
                          trim(str(error))//' returned when '// &
                          'writing attribute to HDF group '//   &
                          group_name)
      ret_err = error
    end if

    call this%thickness%write_hdf(group_id, hdf_thickness, error)
    if (error /= 0) then
      call logger%warning('asym_plume%write_data','Error code '//    &
                          trim(str(error))//' returned when '// &
                          'writing plume thickness field to HDF')
      if (ret_err == 0) ret_err = error
    end if

    call this%velocity%write_hdf(group_id, hdf_velocity, error)
    if (error /= 0) then
      call logger%warning('asym_plume%write_data','Error code '//    &
                          trim(str(error))//' returned when '// &
                          'writing plume velocity field to HDF')
      if (ret_err == 0) ret_err = error
    end if

    call this%temperature%write_hdf(group_id, hdf_temperature, error)
    if (error /= 0) then
      call logger%warning('asym_plume%write_data','Error code '//       &
                          trim(str(error))//' returned when '//    &
                          'writing plume temperature field to HDF')
      if (ret_err == 0) ret_err = error
    end if

    call this%salinity%write_hdf(group_id, hdf_salinity, error)
    if (error /= 0) then
      call logger%warning('asym_plume%write_data','Error code '//    &
                          trim(str(error))//' returned when '// &
                          'writing plume salinity field to HDF')
      if (ret_err == 0) ret_err = error
    end if

    call h5gclose_f(group_id, error)
    if (error /= 0) then
      call logger%warning('asym_plume%write_data','Error code '//    &
                          trim(str(error))//' returned when '// &
                          'closing HDF group '//group_name)
      if (ret_err == 0) ret_err = error
    end if
    error = ret_err
#ifdef DEBUG
    call logger%debug('asym_plume%write_data','Wrote plume data to HDF group '// &
                      group_name)
#endif
  end subroutine asym_plume_write_data


  subroutine asym_plume_solve(this, ice_thickness, ice_density, ice_temperature, &
                         time, success)
    !* Author: Chris MacMackin
    !  Date: March 2017
    !
    ! Solves the state of the plume for the specified ice properties,
    ! at the specified time. This is done using the a
    ! quasilinearisation method and a GMRES iterative linear solver.
    !
    class(asym_plume), intent(inout)     :: this
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
    select type(bound => this%boundaries)
    class is(upstream_plume_boundary)
      call bound%calculate(time, non_diff_terms, b)
    class default
      call bound%set_time(time)
    end select
    solution = this%state_vector()
#ifdef DEBUG
    call logger%debug('asym_plume%solve','Calling QLM ODE solver')
#endif
    call quasilinear_solve(L, f, jac_prod, solution, 1, residual, flag, info, &
                           1.e-9_r8*size(solution), precond=preconditioner, &
                           iter_max=100, krylov_dim=85, gmres_iter_max=5000)
    call this%update(solution)
#ifdef DEBUG
    call logger%debug('plume%solve','QLM solver required '//         &
                      trim(str(info(5)))//' nonlinear iterations '// &
                      'and '//trim(str(info(1)+info(2)))//           &
                      ' function calls.')
#endif

    select case(flag)
    case(0)
      call logger%trivia('asym_plume%solver','Solved plume at time '//trim(str(time)))
      success = .true.
      this%time = time
    case(1)
      call logger%warning('asym_plume%solver','Plume solver stagnated with '// &
                       'residual of '//trim(str(residual)))
      success = .false.
    case(2)
      call logger%error('asym_plume%solver','Reached maximum number of '// &
                        'iterations solving plume')
      success = .false.
    case(3)
      call logger%error('asym_plume%solver','Plume solution began to diverge.')
      success = .false.
    case default
      call logger%error('asym_plume%solve','QLM solver failed for plume with '// &
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
      type(cheb1d_scalar_field) :: scalar_tmp, ctmp(2)
      type(cheb1d_vector_field) :: vector_tmp
      type(cheb1d_vector_field) :: coriolis

      call this%update(v)
      
      ! Thickness
      scalar_tmp = this%thickness%d_dx(1)
      call this%boundaries%thickness_bound_info(-1, btype_l, bdepth_l)
      call this%boundaries%thickness_bound_info(1, btype_u, bdepth_u)
      if (this%lower_bounds(1)) then
        call scalar_tmp%set_boundary(1, 1, this%thickness%get_boundary(-1, 1))
      end if
      if (this%upper_bounds(1)) then
        call scalar_tmp%set_boundary(1, 1, this%thickness%get_boundary(1, 1))
      end if
      st = 1
      en = st + this%thickness_size - 1
      L(st:en) = scalar_tmp%raw()
      
      ! Velocity
      vector_tmp = this%velocity%d_dx(1) - this%velocity_dx
      if (this%lower_bounds(2)) then
        call vector_tmp%set_boundary(1, 1, this%velocity%get_boundary(-1, 1))
      end if
      if (this%upper_bounds(2)) then
        call vector_tmp%set_boundary(1, 1, this%velocity%get_boundary(1, 1))
      end if
      st = en + 1
      en = st + this%velocity_size - 1
      L(st:en) = vector_tmp%raw()

      if (this%phi /= 0._r8) then
        ctmp(1) = -this%phi*this%shape%a_DV/(this%nu*this%shape%a_DU) &
                  *this%velocity%component(2)
        ctmp(2) = this%phi*this%shape%a_DU/(this%nu*this%shape%a_DV) &
                  *this%velocity%component(1)
        coriolis = ctmp
        vector_tmp = this%velocity_dx%d_dx(1) - coriolis
      else
        vector_tmp = this%velocity_dx%d_dx(1)
      end if
      if (this%lower_bounds(3)) then
        call vector_tmp%set_boundary(1, 1, this%velocity_dx%get_boundary(-1, 1))
      end if
      if (this%upper_bounds(3)) then
        call vector_tmp%set_boundary(1, 1, this%velocity_dx%get_boundary(1, 1))
      end if
      st = en + 1
      en = st + this%velocity_size - 1
      L(st:en) = vector_tmp%raw()
 !     print*,vector_tmp%raw()

      ! Temperature
      scalar_tmp = this%temperature%d_dx(1) - this%temperature_dx
      if (this%lower_bounds(4)) then
        call scalar_tmp%set_boundary(1, 1, this%temperature%get_boundary(-1, 1))
      end if
      if (this%upper_bounds(4)) then
        call scalar_tmp%set_boundary(1, 1, this%temperature%get_boundary(1, 1))
      end if
      st = en + 1
      en = st + this%temperature_size - 1
      L(st:en) = scalar_tmp%raw()

      scalar_tmp = this%temperature_dx%d_dx(1)
      if (this%lower_bounds(5)) then
        call scalar_tmp%set_boundary(1, 1, this%temperature_dx%get_boundary(-1, 1))
      end if
      if (this%upper_bounds(5)) then
        call scalar_tmp%set_boundary(1, 1, this%temperature_dx%get_boundary(1, 1))
      end if
      st = en + 1
      en = st + this%temperature_size - 1
      L(st:en) = scalar_tmp%raw()

      ! Salinity
      scalar_tmp = this%salinity%d_dx(1) - this%salinity_dx
      if (this%lower_bounds(6)) then
        call scalar_tmp%set_boundary(1, 1, this%salinity%get_boundary(-1, 1))
      end if
      if (this%upper_bounds(6)) then
        call scalar_tmp%set_boundary(1, 1, this%salinity%get_boundary(1, 1))
      end if
      st = en + 1
      en = st + this%salinity_size - 1
      L(st:en) = scalar_tmp%raw()

      scalar_tmp = this%salinity_dx%d_dx(1)
      if (this%lower_bounds(7)) then
        call scalar_tmp%set_boundary(1, 1, this%salinity_dx%get_boundary(-1, 1))
      end if
      if (this%upper_bounds(7)) then
        call scalar_tmp%set_boundary(1, 1, this%salinity_dx%get_boundary(1, 1))
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
     
      class(scalar_field), pointer :: m, rho_b, rho_t, e, S_a, U, V, &
                                      T_a, rho_a, rho_x, Unorm
      class(scalar_field), allocatable, dimension(:) :: tmp
      type(cheb1d_vector_field) :: coriolis

      call D%guard_temp(); call Uvec%guard_temp(); call T%guard_temp()
      call S%guard_temp(); call b%guard_temp()
      S_a => this%ambient_conds%ambient_salinity(b,time)
      T_a => this%ambient_conds%ambient_temperature(b,time)
      select type(eos => this%eos)
      class is(ave_linear_eos)
        rho_b => eos%water_density_ave1(T, S)
        rho_t => eos%water_density_ave2(T, S)
      class default
        rho_b => eos%water_density(T, S)
        rho_t => rho_b
      end select
      U => Uvec%component(1)
      V => Uvec%component(2)
      call S_a%guard_temp(); call T_a%guard_temp(); call rho_b%guard_temp()
      call rho_t%guard_temp(); call U%guard_temp(); call V%guard_temp()
      rho_a => this%eos%water_density(T_a, S_a)
      call rho_a%guard_temp()
      e => this%entrainment_formulation%entrainment_rate(Uvec, D, b, rho_a - rho_b, time)
      call e%guard_temp()

      call this%melt_formulation%solve_for_melt(Uvec, b, T, S, D, time)
      m => this%melt_formulation%melt_rate()
      call m%guard_temp()

      associate(f_D1 => this%shape%f_D1, f_D2 => this%shape%f_D2,       &
                f_U2 => this%shape%f_U2, f_V2 => this%shape%f_V2,       &
                f_S1 => this%shape%f_S1, f_S2 => this%shape%f_S2,       &
                f_T1 => this%shape%f_T1, f_T2 => this%shape%f_T2,       &
                a_DU => this%shape%a_DU, a_DU2 => this%shape%a_DU2,     &
                a_DUV => this%shape%a_DUV, a_DUS => this%shape%a_DUS,   &
                a_DUT => this%shape%a_DUT, a_D2 => this%shape%a_D2,     &
                a_DV => this%shape%a_DV, a_UabsU => this%shape%a_UabsU, &
                a_UabsV => this%shape%a_UabsV, dy => this%dy)
        DU_x = (e + m - f_D2*D * f_V2*V/dy)/a_DU
        if (this%melt_formulation%has_heat_terms()) then
          DUT_x = (e*T_a - f_D2*D * f_V2*V * f_T2*T/dy - &
                   this%melt_formulation%heat_equation_terms())/a_DUT
        else
          DUT_x = (e*T_a - f_D2*D * f_V2*V * f_T2*T/dy)/a_DUT
        end if
        if (this%melt_formulation%has_salt_terms()) then
          DUS_x = (e*S_a - f_D2*D * f_V2*V * f_S2*S/dy - &
                   this%melt_formulation%salt_equation_terms())/a_DUS
        else
          DUS_x = (e*S_a - f_D2*D * f_V2*V * f_S2*S/dy)/a_DUS
        end if
  
        Unorm => Uvec%norm()
        call Unorm%guard_temp()
        select type(Uvec)
        class is(uniform_vector_field)
          rho_x => this%eos%water_density_derivative(T, (DUT_x - DU_x*T)/(D*U), &
                                                    S, (DUS_x - DU_x*S)/(D*U), 1)
          call rho_x%guard_temp()
          allocate(tmp(this%vel_dims), mold=D)
          tmp(1) = (D*(rho_a - rho_b)*b%d_dx(1) &
                 - 2._r8*D*this%delta*a_D2*(rho_a - rho_t)*DU_x/U &
                 + 0.5*this%delta*a_D2*D**2*rho_x &
                 - this%mu*a_UabsU*Unorm*U &
                 + this%phi*a_DV*D*V &
                 - f_D2*D * f_U2*U * f_V2*V/dy) / &
                   (a_DU2 - this%delta*a_D2*D*(rho_a - rho_t)/U**2)
          if (this%vel_dims > 1) then
             ! Use entrainment and melt as work-arrays to hold the
             ! upper and lower boundary density values. Should be able
             ! to just multiply the function results directly, but
             ! there's a compiler bug.
             call e%clean_temp(); call m%clean_temp()
             e => this%eos%water_density(f_T2*T, f_S2*S)
             m => this%eos%water_density(f_T1*T, f_S1*S)
             call e%guard_temp(); call m%guard_temp()
             tmp(2) = ((f_D2**2*(e - rho_a) - f_D1**2*(m - rho_a))*  &
                      0.5_r8*this%delta/dy*D**2 - this%mu*a_UabsV*Unorm*V       &
                      - a_DU*this%phi*D*U - f_D2*D * f_V2**2*V**2/dy)/a_DUV
          end if
          DUU_x = tmp
          call rho_x%clean_temp()
        class default
           allocate(tmp(this%vel_dims), mold=D)
           tmp(1) = b%d_dx(1)
           tmp(1) = (D*(rho_a - rho_b)*b%d_dx(1) &
                  - a_D2*D*(rho_a - rho_t)*this%delta*D%d_dx(1) &
                  + 0.5_r8*a_D2*this%delta*D**2*rho_t%d_dx(1) &
                  - this%mu*a_UabsU*U*Unorm &
                  - f_D2*D * f_V2*V * f_U2*U/dy)/a_DU2
           if (this%vel_dims > 1) then
             ! Use entrainment and melt as work-arrays to hold the
             ! upper and lower boundary density values. Should be able
             ! to just multiply the function results directly, but
             ! there's a compiler bug.
             call e%clean_temp(); call m%clean_temp()
             e => this%eos%water_density(f_T2*T, f_S2*S)
             m => this%eos%water_density(f_T1*T, f_S1*S)
             call e%guard_temp(); call m%guard_temp()
             tmp(2) = ((f_D2**2*(e - rho_a) - f_D1**2*(m - rho_a))* &
                      0.5_r8*this%delta/dy*D**2 - this%mu*a_UabsV*V*Unorm &
                      - f_D2*D * f_V2**2 * V**2/dy)/a_DUV
           end if
           DUU_x = tmp
        end select
        call Unorm%clean_temp()
      end associate
      call e%clean_temp(); call S_a%clean_temp(); call T_a%clean_temp()
      call rho_b%clean_temp(); call m%clean_temp(); call rho_a%clean_temp()
      call rho_t%clean_temp(); call U%clean_temp(); call V%clean_temp()
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

      call this%update(v(:,1))
      call nonlinear(f, .false.)
    end function f

    
    function jac_prod(v, dv)
      !! The product of the Jacobian of the nonlienar operator at v,
      !! multiplying dv.
      real(r8), dimension(:,:), intent(in) :: v
        !! The state vector for the system of differential equations,
        !! and its derivatives. Column \(i\) represents the \(i-1\)
        !! derivative.
      real(r8), dimension(:,:), intent(in) :: dv
        !! The state vector for the system of differential equations,
        !! and its derivatives, to be multiplied by the
        !! Jacobian. Column \(i\) represents the \(i-1\) derivative.
      real(r8), dimension(size(v,1)) :: jac_prod
      type(cheb1d_scalar_field) :: stmp
      type(cheb1d_vector_field) :: vtmp
      integer :: i

      call this%update(v(:,1))

      call stmp%assign_meta_data(this%thickness)
      call vtmp%assign_meta_data(this%velocity)
      
      call stmp%set_from_raw(dv(1:this%thickness_size, 1))
      call this%thickness%set_derivative(stmp)
      i = 1 + this%thickness_size
      call vtmp%set_from_raw(dv(i:i + this%velocity_size - 1, 1))
      call this%velocity%set_derivative(vtmp)
      i = i + this%velocity_size
      call vtmp%set_from_raw(dv(i:i + this%velocity_size - 1, 1))
      call this%velocity_dx%set_derivative(vtmp)
      i = i + this%velocity_size
      call stmp%set_from_raw(dv(i:i + this%temperature_size - 1, 1))
      call this%temperature%set_derivative(stmp)
      i = i + this%temperature_size
      call stmp%set_from_raw(dv(i:i + this%temperature_size - 1, 1))
      call this%temperature_dx%set_derivative(stmp)
      i = i + this%temperature_size
      call stmp%set_from_raw(dv(i:i + this%salinity_size - 1, 1))
      call this%salinity%set_derivative(stmp)
      i = i + this%salinity_size
      call stmp%set_from_raw(dv(i:i + this%salinity_size - 1, 1))
      call this%salinity_dx%set_derivative(stmp)
      
      call nonlinear(jac_prod, .true.)

      call this%thickness%unset_derivative()
      call this%velocity%unset_derivative()
      call this%velocity_dx%unset_derivative()
      call this%temperature%unset_derivative()
      call this%temperature_dx%unset_derivative()
      call this%salinity%unset_derivative()
      call this%salinity_dx%unset_derivative()
    end function jac_prod

    subroutine nonlinear(f, deriv)
      real(r8), dimension(:), intent(out) :: f
      logical, intent(in) :: deriv
        !! If true, return Jacobian product, otherwise return result
        !! of nonlienar operator.

      integer :: st, en
      type(cheb1d_scalar_field) :: scalar_tmp, D_x, D_nd, S_nd, T_nd
      type(cheb1d_scalar_field), allocatable, dimension(:) :: vtmp
      type(cheb1d_vector_field) :: vector_tmp, U_nd
      class(scalar_field), pointer :: U, V, U_x, V_x

      ! Use same or similar notation for variables as used in equations
      associate(D => this%thickness, Uvec => this%velocity,     &
                Uvec_x => this%velocity_dx, S => this%salinity, &
                S_x => this%salinity_dx, T => this%temperature, &
                T_x => this%temperature_dx, mf => this%melt_formulation, &
                h => ice_thickness, delta => this%delta, nu => this%nu,  &
                mu => this%mu, r => this%r_val, bounds => this%boundaries,&
                a_DU2 => this%shape%a_DU2, a_DUV => this%shape%a_DUV,     &
                a_DUS => this%shape%a_DUS, a_DUT => this%shape%a_DUT,     &
                a_DU => this%shape%a_DU, a_DV => this%shape%a_DV,         &
                a_DS => this%shape%a_DS, a_DT => this%shape%a_DT,         &
                f_D2 => this%shape%f_D2, f_Up => this%shape%f_Up,           &
                f_Vp => this%shape%f_Vp, f_Sp => this%shape%f_Sp,         &
                f_Tp => this%shape%f_Tp, dy => this%dy )

        call non_diff_terms(D, Uvec, T, S, b, D_nd, U_nd, T_nd, S_nd)

        U => this%velocity%component(1)
        U_x => this%velocity_dx%component(1)
        call U%guard_temp(); call U_x%guard_temp()
        if (this%vel_dims > 1) then
          V => this%velocity%component(2)
          V_x => this%velocity_dx%component(2)
          call V%guard_temp(); call V_x%guard_temp()
        end if

        ! Thickness
        scalar_tmp = (D_nd - D*U_x)/U
        D_x = scalar_tmp
        if (this%lower_bounds(1)) then
          call scalar_tmp%set_boundary(1, 1, bounds%thickness_bound(-1))
        end if
        if (this%upper_bounds(1)) then
          call scalar_tmp%set_boundary(1, 1, bounds%thickness_bound(1))
        end if
        st = 1
        en = st + this%thickness_size - 1
        if (deriv) then
          scalar_tmp = scalar_tmp%get_derivative()
        end if
        f(st:en) = scalar_tmp%raw()

        ! Velocity
        vector_tmp = 0._r8 * Uvec
        if (this%lower_bounds(2)) then
          call vector_tmp%set_boundary(1, 1, bounds%velocity_bound(-1))
        end if
        if (this%upper_bounds(2)) then
          call vector_tmp%set_boundary(1, 1, bounds%velocity_bound(1))
        end if
        st = en + 1
        en = st + this%velocity_size - 1
        if (deriv) then
          vector_tmp = vector_tmp%get_derivative()
        end if
        f(st:en) = vector_tmp%raw()
        
        allocate (vtmp(this%vel_dims), mold=D)
        vtmp(1) = (a_DU2*(2._r8*D*U*U_x + D_x*U**2 - U_nd%component(1)) &
                - nu*a_DU*D_x*U_x - nu*f_D2*D*f_Up*U/dy)/(nu*a_DU*D)
        if (this%vel_dims > 1) then
          vtmp(2) = (a_DUV*(D*U*V_x + D*U_x*V + D_x*U*V - U_nd%component(2)) &
                  - nu*a_DV*D_x*V_x - nu*f_D2*D*f_Vp*V/dy)/(nu*a_DV*D)
        end if
        vector_tmp = vtmp
        if (this%lower_bounds(3)) then
          call vector_tmp%set_boundary(1, 1, bounds%velocity_bound(-1))
        end if
        if (this%upper_bounds(3)) then
          call vector_tmp%set_boundary(1, 1, bounds%velocity_bound(1))
        end if
        st = en + 1
        en = st + this%velocity_size - 1
        if (deriv) then
          vector_tmp = vector_tmp%get_derivative()
        end if
        f(st:en) = vector_tmp%raw()

        ! Temperature
        scalar_tmp = uniform_scalar_field(0._r8)
        if (this%lower_bounds(4)) then
          call scalar_tmp%set_boundary(1, 1, bounds%temperature_bound(-1))
        end if
        if (this%upper_bounds(4)) then
          call scalar_tmp%set_boundary(1, 1, bounds%temperature_bound(1))
        end if
        st = en + 1
        en = st + this%temperature_size - 1
        if (deriv) then
          scalar_tmp = scalar_tmp%get_derivative()
        end if
        f(st:en) = scalar_tmp%raw()

        scalar_tmp = (a_DUT*(D*U*T_x + D*U_x*T + D_x*U*T - T_nd) - &
                      nu*a_DT*D_x*T_x - nu*f_D2*D*f_Tp*T/dy)/(nu*a_DT*D)
        if (this%lower_bounds(5)) then
          call scalar_tmp%set_boundary(1, 1, bounds%temperature_bound(-1))
        end if
        if (this%upper_bounds(5)) then
          call scalar_tmp%set_boundary(1, 1, bounds%temperature_bound(1))
        end if
        st = en + 1
        en = st + this%salinity_size - 1
        if (deriv) then
          scalar_tmp = scalar_tmp%get_derivative()
        end if
        f(st:en) = scalar_tmp%raw()

        ! Salinity
        scalar_tmp = uniform_scalar_field(0._r8)
        if (this%lower_bounds(6)) then
          call scalar_tmp%set_boundary(1, 1, bounds%salinity_bound(-1))
        end if
        if (this%upper_bounds(6)) then
          call scalar_tmp%set_boundary(1, 1, bounds%salinity_bound(1))
        end if
        st = en + 1
        en = st + this%salinity_size - 1
        if (deriv) then
          scalar_tmp = scalar_tmp%get_derivative()
        end if
        f(st:en) = scalar_tmp%raw()

        scalar_tmp = (a_DUS*(D*U*S_x + D*U_x*S + D_x*U*S - S_nd) - &
                      nu*a_DS*D_x*S_x - nu*f_D2*D*f_Sp*S/dy)/(nu*a_DS*D)
        if (this%lower_bounds(7)) then
          call scalar_tmp%set_boundary(1, 1, bounds%salinity_bound(-1))
        end if
        if (this%upper_bounds(7)) then
          call scalar_tmp%set_boundary(1, 1, bounds%salinity_bound(1))
        end if
        st = en + 1
        en = st + this%salinity_size - 1
        if (deriv) then
          scalar_tmp = scalar_tmp%get_derivative()
        end if
        f(st:en) = scalar_tmp%raw()

        call U%clean_temp(); call U_x%clean_temp()
        if (this%vel_dims > 0) then
          call V%clean_temp(); call V_x%clean_temp()
        end if
      end associate
    end subroutine nonlinear

    function preconditioner(v, state, L_op, f_op, fcur, rhs)
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
      real(r8), dimension(:), intent(in)   :: fcur
        !! The result of `f(u)`
      real(r8), dimension(:), intent(in)   :: rhs
        !! The right hand side of the linear system being
        !! preconditioned.
      real(r8), dimension(size(v)) :: preconditioner
        !! The result of applying the preconditioner.

      integer :: st, en, ust, uen, pst, pen
      integer :: bloc, i
      
      real(r8) :: nu
      type(asym_plume) :: v_plume
      type(cheb1d_scalar_field) :: scalar_tmp
      type(cheb1d_vector_field) :: vector_tmp, tmp2
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

      nu = this%nu
      bloc = get_bound_loc(1)  
      v_plume%thickness = this%precond%solve_for(v_plume%thickness, bloc, &
           v_plume%thickness%get_boundary(1, 1), -1)
      st = 1
      en = st + this%thickness_size - 1
      preconditioner(st:en) = v_plume%thickness%raw()
  
      bloc = get_bound_loc(3)
      v_plume%velocity_dx = this%precond%solve_for(v_plume%velocity_dx, bloc, &
           v_plume%velocity_dx%get_boundary(1, 1), -1)

      bloc = get_bound_loc(2)
      vector_tmp = v_plume%velocity + v_plume%velocity_dx
      v_plume%velocity = this%precond%solve_for(vector_tmp, bloc, &
           v_plume%velocity%get_boundary(1, 1), -1)
      st = en + 1
      en = st + this%velocity_size - 1
      preconditioner(st:en) = v_plume%velocity%raw()
      st = en + 1
      en = st + this%velocity_size - 1
      preconditioner(st:en) = v_plume%velocity_dx%raw()        

      ! Precondition T_x terms before T
      st = en + 1
      en = st + this%temperature_size - 1
      pst = st
      pen = en
  
      bloc = get_bound_loc(5)
      v_plume%temperature_dx = this%precond%solve_for(v_plume%temperature_dx, bloc, &
           v_plume%temperature_dx%get_boundary(1, 1), -1)
      st = en + 1
      en = st + this%temperature_size - 1
      preconditioner(st:en) = v_plume%temperature_dx%raw()
  
      bloc = get_bound_loc(4)
      v_plume%temperature = this%precond%solve_for(v_plume%temperature + &
           v_plume%temperature_dx, bloc, v_plume%temperature%get_boundary(1, 1), -1)
      preconditioner(pst:pen) = v_plume%temperature%raw()

      ! Precondition S_x terms before S
      st = en + 1
      en = st + this%salinity_size - 1
      pst = st
      pen = en
  
      bloc = get_bound_loc(7)
      v_plume%salinity_dx = this%precond%solve_for(v_plume%salinity_dx, bloc, &
           v_plume%salinity_dx%get_boundary(1, 1), -1)
      st = en + 1
      en = st + this%temperature_size - 1
      preconditioner(st:en) = v_plume%salinity_dx%raw()

      bloc = get_bound_loc(6)
      v_plume%salinity = this%precond%solve_for(v_plume%salinity+v_plume%salinity_dx, bloc, &
           v_plume%salinity%get_boundary(1, 1), -1)
      preconditioner(pst:pen) = v_plume%salinity%raw()
    end function preconditioner

    integer function get_bound_loc(component_id)
      integer :: component_id
      if (this%lower_bounds(component_id)) then
        get_bound_loc = -1
      else if (this%upper_bounds(component_id)) then
        get_bound_loc = 1
      else
        get_bound_loc = 0
      end if
    end function get_bound_loc
  end subroutine asym_plume_solve

end module asymmetric_plume_mod

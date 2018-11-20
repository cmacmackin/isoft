!
!  ice_shelf.f90
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

#ifdef DEBUG
#define pure 
#define elemental 
#endif

module ice_shelf_mod
  !* Author: Christopher MacMackin
  !  Date: April 2016
  !  License: GPLv3
  !
  ! Provides a concrete implementation of the [[glacier]] type, using
  ! a vertically integrated model of an ice shelf.
  !
  use iso_fortran_env, only: r8 => real64
  !use foodie, only: integrand
  use glacier_mod, only: glacier, thickness_func, velocity_func, hdf_type_attr
  use factual_mod, only: scalar_field, vector_field, cheb1d_scalar_field, &
                         cheb1d_vector_field, maxval, minval, abs
  use viscosity_mod, only: abstract_viscosity
  use newtonian_viscosity_mod, only: newtonian_viscosity
  use glacier_boundary_mod, only: glacier_boundary
  use dallaston2015_glacier_boundary_mod, only: dallaston2015_glacier_boundary
  use jacobian_block_mod, only: jacobian_block
  use boundary_types_mod
  use nitsol_mod
  use hdf5
  use h5lt
  use logger_mod, only: logger => master_logger
  use penf, only: str
  implicit none
  private

  character(len=9), parameter, public  :: hdf_type_name = 'ice_shelf'
  character(len=9), parameter, public  :: hdf_thickness = 'thickness'
  character(len=8), parameter, public  :: hdf_velocity = 'velocity'
  character(len=6), parameter, public  :: hdf_lambda = 'lambda'
  character(len=4), parameter, public  :: hdf_zeta = 'zeta'
  character(len=3), parameter, public  :: hdf_chi = 'chi'
  character(len=10), parameter, public :: hdf_n_kappa = 'num_kappas'
  character(len=15), parameter, public :: hdf_kappa = '("kappa_",i0.4)'

  type, extends(glacier), public :: ice_shelf
    !* Author: Chris MacMackin
    !  Date: April 2016
    !
    ! A concrete implementation of the [[glacier]] type, using a vertically
    ! integrated model of an ice shelf. This model is 1-dimensional only.
    !
    private
    type(cheb1d_scalar_field) :: thickness
      !! Thickness of ice shelf, \(h\)
    type(cheb1d_vector_field) :: velocity
      !! Flow velocity of ice shelf, \(\vec{u}\)
    type(cheb1d_scalar_field) :: eta
      !! Viscosity of the ice, \(\eta\)
    type(cheb1d_scalar_field), dimension(:), allocatable :: kappa
      !! Taylor coefficients for the vertical structure of a
      !! Lagrangian tracer representing englacial layers/internal
      !! reflectors.
    real(r8)                  :: lambda
      !! The dimensionless ratio 
      !! $\lambda \equiv \frac{\rho_0m_0x_0}{\rho_iH-0u_0}$
    real(r8)                  :: chi
      !! The dimensionless ratio \(\chi \equiv
      !! \frac{\rho_igh_0x_0}{2\eta_0u_0} \left(1 -
      !! \frac{\rho_i}{\rho_o}\right)\)
    real(r8)                  :: zeta
      !! The dimensionless ratio \(\zeta \equiv
      !! \frac{\rho_iu_0x_0}{\eta_0}\), corresponding to the Reynolds
      !! number. Currently unused.
    real(r8)                  :: courant
      !! The Courant number to use when calculating the time step.
    class(abstract_viscosity), allocatable :: viscosity_law
      !! An object representing the model used for ice viscosity.
    class(glacier_boundary), allocatable   :: boundaries
      !! An object specifying the boundary conditions for the ice
      !! shelf.
    real(r8)                  :: max_dt
      !! The maximu  allowable time step
    real(r8)                  :: time
      !! The time at which the ice shelf is in this state.
    integer                   :: thickness_size
      !! The number of data values in the thickness field.
    integer                   :: velocity_size
      !! The number of data values in the velocity field.
    integer                   :: boundary_start
      !! The number of data values needed to represent the boundary
      !! conditions.
    integer                   :: thickness_lower_bound_size
      !! The number of data values needed to represent the lower
      !! boundary conditions for thickness.
    integer                   :: thickness_upper_bound_size
      !! The number of data values needed to represent the upper
      !! boundary conditions for thickness.
    integer                   :: velocity_lower_bound_size
      !! The number of data values needed to represent the lower
      !! boundary conditions for velocity.
    integer                   :: velocity_upper_bound_size
      !! The number of data values needed to represent the upper
      !! boundary conditions for thickness.
    type(jacobian_block)      :: thickness_jacobian
      !! A representation of the Jacobian for the ice shelf thickness.
    type(jacobian_block)      :: velocity_jacobian
      !! A representation of the Jacobian for the ice shelf velocity.
    logical                   :: stale_eta
      !! Indicates whether the viscosity needs updating.
    logical                   :: stale_jacobian
      !! Indicates if the Jacobians are stale and in need of updating.
  contains
    procedure :: initialise => shelf_initialise
    procedure :: ice_thickness => shelf_thickness
!$    procedure :: ice_velocity => shelf_velocity
    procedure :: ice_density => shelf_density
    procedure :: ice_temperature => shelf_temperature
    procedure :: residual => shelf_residual
    procedure :: update => shelf_update
    procedure :: precondition =>  shelf_precondition
    procedure :: set_time => shelf_set_time
    procedure :: data_size => shelf_data_size
    procedure :: state_vector => shelf_state_vector
    procedure :: kappa_vector => shelf_kappa_vector
    procedure :: read_data => shelf_read_data
    procedure :: write_data => shelf_write_data
    procedure :: time_step => shelf_time_step
    procedure :: solve_velocity => shelf_solve_velocity
!    procedure :: integrate => shelf_integrate
    procedure, private :: assign => shelf_assign
    procedure :: integrate_layers => ice_shelf_integrate_layers
  end type ice_shelf

!  interface ice_shelf
!    module procedure constructor
!  end interface ice_shelf

  abstract interface
    pure function kappa_init_func(n, location) result(kappa)
      !* Author: Chris MacMackin
      !  Date: September 2018
      !
      ! Abstract interface for function providing the Taylor
      ! coefficients describing the distribution of internal
      ! reflectors within an ice shelf.
      !
      import :: r8
      integer, intent(in)                :: n
        !! The index of the Taylor coefficient being calculated
      real(r8), dimension(:), intent(in) :: location
        !! The position $\vec{x}$ at which to compute the coefficient
      real(r8) :: kappa
        !! The value of coefficient `n` at `location`
    end function kappa_init_func
  end interface

contains
  
  subroutine shelf_initialise(this, domain, resolution, thickness, velocity, &
                              temperature, viscosity_law, boundaries, lambda, &
                              chi, zeta, courant, max_dt, kappa, n_kappa)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Initialises an [[ice_shelf]] object with initial conditions provided
    ! by the arguments. At present only a 1D model is supported. If
    ! information is provided for higher dimensions then it will be ignored.
    !
    class(ice_shelf), intent(out)        :: this
    real(r8), dimension(:,:), intent(in) :: domain
      !! An array containing the upper and lower limits of the domain for
      !! the ice shelf. The first index represents the dimension for which
      !! the boundaries apply. If the second index is 1 then it corresponds
      !! to the lower bound. If the second index is 2 then it corresponds to
      !! the upper bound.
    integer, dimension(:), intent(in)    :: resolution
      !! The number of data points in each dimension.
    procedure(thickness_func)            :: thickness
      !! A function which calculates the initial value of the thickness of 
      !! the ice shelf at a given location.
    procedure(velocity_func)             :: velocity
      !! A function which calculates the initial value of the velocity 
      !! (vector) of the ice at a given location in an ice shelf.
    real(r8), intent(in), optional       :: temperature
      !! The temperature of the ice in the ice shelf.
    class(abstract_viscosity), allocatable, optional, &
                           intent(inout) :: viscosity_law
      !! An object which calculates the viscosity of the ice. If not
      !! specified, then Glen's law will be used with $n=3$. Will be
      !! unallocated on return.
    class(glacier_boundary), allocatable, optional, &
                         intent(inout)   :: boundaries
      !! An object specifying the boundary conditions for the ice
      !! shelf. Will be unallocated on return.
    real(r8), intent(in), optional       :: lambda
      !! The dimensionless ratio 
      !! $\lambda \equiv \frac{\rho_0m_0x_0}{\rho_ih_0u_0}$.
    real(r8), intent(in), optional       :: chi
      !! The dimensionless ratio
      !! $\chi \equiv \frac{\rho_igh_0x_x}{2\eta_0u_0}\left(1 -
      !! \frac{\rho_i}{\rho_0}\right)$.
    real(r8), intent(in), optional       :: zeta
      !! The dimensionless ratio $\zeta \equiv
      !! \frac{\rho_iu_0x_0}{\eta_0}$, corresponding to the Reynolds
      !! number. Currently this is unused and always treated as 0.
    real(r8), intent(in), optional       :: courant
      !! The Courant number to use when calculating the time
      !! step. Defaults to 100. Too large a value will pose
      !! difficulties for the nonlinear solver, while too small a
      !! value can be numerically unstable. Typically, smaller values
      !! are needed for lower resolution.
    real(r8), intent(in), optional       :: max_dt
      !! The maximum allowable time step. This defaults to \(1\times
      !! 10^{99}\) (effectively no maximum).
    procedure(kappa_init_func), optional :: kappa
      !! A function which specifies the initial values of the Taylor
      !! coefficients describing the vertical distribution of internal
      !! reflectors within the ice. The initial conditions at the
      !! grounding line will provide the boundary conditions there
      !! throughout the simulation. If this parameter is not provided
      !! then these layers will not be included in the
      !! integration. Both this parameter and `n_kappa` must be
      !! specified for the calculation to take place.
    integer, optional, intent(in)        :: n_kappa
      !! The number of Taylor coefficients used to describe internal
      !! reflectors. If not provided then these reflectors will not be
      !! included in the integration. Both this parameter and `kappa`
      !! must be specified for the calculation to take place.

    integer :: n
    integer, dimension(:), allocatable :: lower, upper

    this%thickness = cheb1d_scalar_field(resolution(1),thickness,domain(1,1), &
                                         domain(1,2))
    this%velocity = cheb1d_vector_field(resolution(1),velocity,domain(1,1), &
                                        domain(1,2))
    this%eta = cheb1d_scalar_field(resolution(1), lower_bound=domain(1,1), &
                                   upper_bound=domain(1,2))
    this%thickness_size = this%thickness%raw_size()
    this%velocity_size = this%velocity%raw_size()

    if (present(kappa) .and. present(n_kappa)) then
      allocate(this%kappa(n_kappa))
      do n = 1, n_kappa
        this%kappa(n) = cheb1d_scalar_field(resolution(1),kappa_func,domain(1,1), &
                                            domain(1,2))
      end do
    end if

    if (present(temperature)) then
      continue ! This doesn't do anything at the moment
    else
      continue ! Again, doesn't do anything at the moment
    end if

    if (present(viscosity_law)) then
      call move_alloc(viscosity_law, this%viscosity_law)
    else
      allocate(newtonian_viscosity:: this%viscosity_law)
    end if

    if (present(boundaries)) then
      call move_alloc(boundaries, this%boundaries)
    else
      allocate(dallaston2015_glacier_boundary :: this%boundaries)
    end if

    if (present(lambda)) then
      this%lambda = lambda
    else
      this%lambda = 0.37_r8
    end if

    if (present(chi)) then
      this%chi = chi
    else
      this%chi = 4.0_r8
    end if

    if (present(zeta)) then
      this%zeta = zeta
    else
      this%zeta = 1.3e-11_r8
    end if

    if (present(courant)) then
      this%courant = courant
    else
      this%courant = 1e2_r8
    end if

    if (present(max_dt)) then
      this%max_dt = max_dt
    else
      this%max_dt = 1e99_r8
    end if

    lower = this%boundaries%thickness_lower_bound()
    upper = this%boundaries%thickness_upper_bound()
    this%thickness_lower_bound_size = this%thickness_size &
                                    - this%thickness%raw_size(lower)
    this%thickness_upper_bound_size = this%thickness_size &
                                    - this%thickness%raw_size(upper)

    lower = this%boundaries%velocity_lower_bound()
    upper = this%boundaries%velocity_upper_bound()
    this%velocity_lower_bound_size = this%velocity_size &
                                   - this%velocity%raw_size(lower)
    this%velocity_upper_bound_size = this%velocity_size &
                                   - this%velocity%raw_size(upper)

    this%boundary_start = this%thickness_size + this%velocity_size + 1 &
                        - this%thickness_lower_bound_size &
                        - this%thickness_upper_bound_size &
                        - this%velocity_lower_bound_size &
                        - this%velocity_upper_bound_size

    this%time = 0.0_r8
    this%stale_eta = .true.
    this%stale_jacobian = .true.
#ifdef DEBUG
    call logger%debug('ice_shelf','Initialised new ice shelf object')
#endif

  contains

    pure function kappa_func(location) result(thickness)
      !* Author: Chris MacMackin
      !  Date: April 2016
      !
      ! Wrapper around user-provided `kappa` routine, converting it to
      ! the form needed to initialise field-types.
      !
      real(r8), dimension(:), intent(in) :: location
        !! The position $\vec{x}$ at which to compute the thickness
      real(r8) :: thickness
        !! The thickness of the glacier at `location`
      thickness = kappa(n, location)
    end function kappa_func

  end subroutine shelf_initialise


  function shelf_thickness(this) result(thickness)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Returns the thickness of the ice shelf across its domain.
    !
    class(ice_shelf), intent(in) :: this
    class(scalar_field), pointer :: thickness !! The ice thickness.
    call this%thickness%allocate_scalar_field(thickness)
    thickness = this%thickness
#ifdef DEBUG
    call logger%debug('ice_shelf%thickness','Returned ice shelf thickness')
#endif    
  end function shelf_thickness


  function shelf_velocity(this) result(velocity)
    !* Author: Christopher MacMackin
    !  Date: July 2016
    !
    ! Returns the velocity of the ice shelf across its domain.
    !
    class(ice_shelf), intent(in) :: this
    class(vector_field), pointer :: velocity !! The ice velocity.
    call this%velocity%allocate_vector_field(velocity)
    velocity = this%velocity
#ifdef DEBUG
    call logger%debug('ice_shelf%velocity','Returned ice shelf velocity')
#endif
  end function shelf_velocity


  pure function shelf_density(this) result(density)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Returns the density of the ice in the shelf, which is assumed to be
    ! uniform across its domain.
    !
    ! @NOTE Based on my approach to non-dimensionalisation, I'm pretty
    ! sure the density should always be 1, making this method
    ! unneccessary.
    !
    class(ice_shelf), intent(in) :: this
    real(r8)                     :: density !! The ice density.
    density = 1.0_r8/1.12_r8 !TODO: Will probably want to change this at some point
#ifdef DEBUG
    call logger%debug('ice_shelf%density','Ice shelf has density '// &
                      trim(str(density))//'.')
#endif
  end function shelf_density


  pure function shelf_temperature(this) result(temperature)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Returns the density of the ice in the shelf, which is assumed to be
    ! uniform across its domain.
    !
    class(ice_shelf), intent(in) :: this
    real(r8)                     :: temperature !! The ice density.
    temperature = -15.0_r8 !TODO: Will probably want to change this at some point.
#ifdef DEBUG
    call logger%debug('ice_shelf%temperature','Ice shelf has temperature '// &
                      trim(str(temperature)))
#endif
  end function shelf_temperature


  function shelf_residual(this, previous_states, melt_rate, &
                          basal_drag_parameter, water_density) result(residual) 
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Returns the residual when the current state of the glacier is run
    ! through the system of equations describing it. The residual takes the
    ! form of a 1D array, with each element respresenting the residual for
    ! one of the equations in the system.
    !
    class(ice_shelf), intent(in)             :: this
    class(glacier), dimension(:), intent(in) :: previous_states
      !! The states of the glacier in the previous time steps. The
      !! first element of the array should be the most recent. The
      !! default implementation will only make use of the most recent
      !! state, but the fact that this is an array allows potential
      !! other implementations to use older states for higher-order
      !! integration methods.
    class(scalar_field), intent(in)          :: melt_rate
      !! Thickness of the ice above the glacier.
    class(scalar_field), intent(in)          :: basal_drag_parameter
      !! A paramter, e.g. coefficient of friction, needed to calculate the
      !! drag on basal surface of the glacier.
    real(r8), intent(in)                     :: water_density
      !! The density of the water below the glacier.
    real(r8), dimension(:), allocatable      :: residual
      !! The residual of the system of equations describing the glacier.
    type(cheb1d_scalar_field) :: scalar_tmp
    integer :: start, finish, bounds_start, bounds_finish
    integer, dimension(:), allocatable :: lower, upper
    real(r8), dimension(:), allocatable :: bounds
    logical :: success

    call melt_rate%guard_temp(); call basal_drag_parameter%guard_temp()
    allocate(residual(this%data_size()))
    start = 1

    ! Use same or similar notation for variables as in equations
    select type(previous_states)
    class is(ice_shelf)
      associate(h => this%thickness, h_old => previous_states(1)%thickness, &
                uvec => this%velocity, m => melt_rate, eta => this%eta,     &
                lambda => this%lambda, t_old => previous_states(1)%time)
        ! Boundary conditions
        if (this%stale_eta) then
          bounds = this%boundaries%boundary_residuals(h, uvec, &
               this%viscosity_law%ice_viscosity(uvec, this%ice_temperature(), &
                                               this%time), this%time)
        else
          bounds = this%boundaries%boundary_residuals(h, uvec, eta, this%time)
        end if

        ! Continuity equation
        scalar_tmp = (h - h_old)/(this%time - t_old) + .div.(h*uvec) + lambda*m
        
        lower = this%boundaries%thickness_lower_bound()
        upper = this%boundaries%thickness_upper_bound()
        ! TODO: Figure out how to make this independent of order which
        ! values are stored in the field
        finish = start + this%thickness_upper_bound_size - 1
        bounds_start = this%thickness_lower_bound_size + 1
        bounds_finish = this%thickness_lower_bound_size &
                      + this%thickness_upper_bound_size
        residual(start:finish) = bounds(bounds_start:bounds_finish)
        start = finish + 1
        finish = start + scalar_tmp%raw_size(lower,upper) - 1
        residual(start:finish) = scalar_tmp%raw(lower,upper)
        start = finish + 1
        finish = start + this%thickness_lower_bound_size - 1
        bounds_start = 1
        bounds_finish = this%thickness_lower_bound_size
        residual(start:finish) = bounds(bounds_start:bounds_finish)
        start = finish + 1

      end associate
    class default
      call logger%fatal('ice_shelf%residual','Type other than `ice_shelf` '// &
                        'passed to `ice_shelf` object as a previous state.')
      error stop
    end select

    call melt_rate%clean_temp(); call basal_drag_parameter%clean_temp()
#ifdef DEBUG
    call logger%debug('ice_shelf%residual','Calculated residual of ice shelf.')
#endif
  end function shelf_residual


  subroutine shelf_update(this, state_vector)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Updates the state of the ice shelf from its state vector. The state
    ! vector is a real array containing the value of each of the ice shelf's
    ! properties at each of the locations on the grid used in discretization.
    !
    class(ice_shelf), intent(inout)     :: this
    real(r8), dimension(:), intent(in)  :: state_vector
      !! A real array containing the data describing the state of the
      !! glacier.
    !TODO: Add some assertion-like checks that the state vector is the right size
    call this%thickness%set_from_raw(state_vector(1:this%thickness_size))
    this%stale_jacobian = .true.
#ifdef DEBUG
    call logger%debug('ice_shelf%update','Updated state of ice shelf.')
#endif
  end subroutine shelf_update


  function shelf_precondition(this, previous_states, melt_rate, &
                              basal_drag_parameter, water_density, &
                              delta_state) result(preconditioned)
    !* Author: Chris MacMackin
    !  Date: January 2016
    !
    ! Provides a preconditioner for the nonlinear solver trying to
    ! bring the residual to zero. The Jacobian is approximated as a
    ! block matrix, where each block is a tridiagonal matrix using a
    ! finite difference method for differentiation.
    !
    class(ice_shelf), intent(inout)          :: this
    class(glacier), dimension(:), intent(in) :: previous_states
      !! The states of the glacier in the previous time steps. The
      !! first element of the array should be the most recent. The
      !! default implementation will only make use of the most
      !! recent state, but the fact that this is an array allows
      !! overriding methods to use older states for higher-order
      !! integration methods.
    class(scalar_field), intent(in)          :: melt_rate
      !! Thickness of the ice above the glacier
    class(scalar_field), intent(in)          :: basal_drag_parameter
      !! A paramter, e.g. coefficient of friction, needed to calculate
      !! the drag on basal surface of the glacier.
    real(r8), intent(in)                     :: water_density
      !! The density of the water below the glacier
    real(r8), dimension(:), intent(in)       :: delta_state
      !! The change to the state vector which is being preconditioned.
    real(r8), dimension(:), allocatable      :: preconditioned
      !! The result of applying the preconditioner to `delta_state`.

    type(cheb1d_scalar_field) :: delta_h
    integer :: i, sl, el, su, eu
    integer, dimension(:), allocatable :: boundary_locations
    real(r8) :: delta_t

    call melt_rate%guard_temp(); call basal_drag_parameter%guard_temp()
    allocate(preconditioned(size(delta_state)))
    select type(previous_states)
    class is(ice_shelf)
      delta_t = this%time - previous_states(1)%time
    class default
      call logger%fatal('ice_shelf%precondition','Type other than `ice_shelf` '// &
                        'passed to `ice_shelf` object as a previous state.')
      error stop
    end select
    call delta_h%assign_meta_data(this%thickness)
    call delta_h%set_from_raw(delta_state)

    associate(jac => this%thickness_jacobian, uvec => this%velocity)
      if (this%stale_jacobian) then
        sl = this%thickness_size - this%thickness_lower_bound_size + 1
        el = this%thickness_size
        su = 1
        eu = this%thickness_upper_bound_size
        boundary_locations = [(i, i=sl,el), (i, i=su,eu)]

        jac = jacobian_block(uvec%component(1), 1, boundary_locs=boundary_locations) &
            + 1._r8/delta_t

        this%stale_jacobian = .false.
      end if

      delta_h = jac%solve_for(delta_h)
    end associate
    
    preconditioned = delta_h%raw()
    call melt_rate%clean_temp(); call basal_drag_parameter%clean_temp()
#ifdef DEBUG
    call logger%debug('ice_shelf%precondition','Applied preconditioner for ice shelf.')
#endif
  end function shelf_precondition


  subroutine shelf_set_time(this, time)
    !* Author: Christopher MacMackin
    !  Date: November 2016
    !
    ! Sets the time information held by the ice shelf object. This is
    ! the time at which the ice sheet is in its current state.
    !
    class(ice_shelf), intent(inout) :: this
    real(r8), intent(in)            :: time
      !! The time at which the glacier is in the present state.
    this%time = time
#ifdef DEBUG
    call logger%debug('ice_shelf%set_time','Updating time for ice shelf to '// &
                      trim(str(time)))
#endif
  end subroutine shelf_set_time


  pure function shelf_data_size(this)
    !* Author: Christopher MacMackin
    !  Date: August 2016
    !
    ! Returns the number of elements in the ice shelf's state vector.
    ! This is the size of the vector returned by [[ice_shelf:residual]]
    ! and [[ice_shelf:state_vector]] and taken as an argument by 
    ! [[ice_shelf:update]].
    !
    class(ice_shelf), intent(in) :: this
    integer                      :: shelf_data_size
      !! The number of elements in the ice shelf's state vector.
    shelf_data_size = this%thickness_size
#ifdef DEBUG
    call logger%debug('ice_shelf%data_size','Ice shelf has '//   &
                      trim(str(shelf_data_size))//' elements '// &
                      'in its state vector.')
#endif
  end function shelf_data_size


  function shelf_state_vector(this) result(state_vector) 
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Returns the state vector for the current state of the ice shelf. 
    ! This takes the form of a 1D array.
    !
    class(ice_shelf), intent(in)        :: this
    real(r8), dimension(:), allocatable :: state_vector
      !! The state vector describing the ice shelf.
    state_vector = this%thickness%raw()
#ifdef DEBUG
    call logger%debug('ice_shelf%state_vector','Returning state vector '// &
                      'for ice shelf.')
#endif
  end function shelf_state_vector


  function shelf_kappa_vector(this) result(kappa_vector) 
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Returns the a vector representing the current state of the
    ! internal reflectors in the ice shelf.  This takes the form of a
    ! 1D array. The routien is only used for debugging purposes.
    !
    class(ice_shelf), intent(in)        :: this
    real(r8), dimension(:), allocatable :: kappa_vector
      !! The state vector describing the ice shelf.
    integer :: i
    if (allocated(this%kappa)) then
      allocate(kappa_vector(this%thickness_size*size(this%kappa)))
      do i = 1, size(this%kappa)
        kappa_vector((i-1)*this%thickness_size + 1:this%thickness_size*i) = this%kappa(i)%raw()
      end do
    end if
#ifdef DEBUG
    call logger%debug('ice_shelf%state_vector','Returning state vector '// &
                      'for ice shelf.')
#endif
  end function shelf_kappa_vector


  subroutine shelf_read_data(this,file_id,group_name,error)
    !* Author: Chris MacMackin
    !  Date: April 2017
    !
    ! Reads the state of the ice shelf object from the specified group
    ! in an HDF5 file. This sets the thickness, the velocity, and
    ! parameter values.
    !
    class(ice_shelf), intent(inout) :: this
    integer(hid_t), intent(in)      :: file_id
      !! The identifier for the HDF5 file/group in which this data is
      !! meant to be written.
    character(len=*), intent(in)    :: group_name
      !! The name to give the group in the HDF5 file storing the
      !! ice shelf's data.
    integer, intent(out)            :: error
      !! Flag indicating whether routine ran without error. If no
      !! error occurs then has value 0.
    integer(hid_t) :: group_id
    integer :: ret_err, i, nkap
    real(r8), dimension(1) :: param
    integer, dimension(1) :: iparam
    character(len=20) :: fieldname
    character(len=50) :: ice_type

    ret_err = 0
    call h5gopen_f(file_id, group_name, group_id, error)
    if (error /= 0) then
      call logger%error('ice_shelf%read_data','Could not open HDF group "'// &
                        group_name//'", so no IO performed.')
      return
    end if

    call h5ltget_attribute_string_f(file_id, group_name, hdf_type_attr, &
                                    ice_type, error)
    if (trim(ice_type) /= hdf_type_name) then
      call logger%error('ice_shelf%read_data','Trying to read data from '// &
                        'glacier of type other than ice_shelf.')
      error = -1
      return
    end if
    call h5ltget_attribute_double_f(file_id, group_name, hdf_lambda, &
                                    param, error)
    this%lambda = param(1)
    call h5ltget_attribute_double_f(file_id, group_name, hdf_chi, &
                                    param, error)
    this%chi = param(1)
    call h5ltget_attribute_double_f(file_id, group_name, hdf_zeta, &
                                    param, error)
    this%zeta = param(1)
    if (error /= 0) then
      call logger%warning('ice_shelf%read_data','Error code '//  &
                          trim(str(error))//' returned when '//  &
                          'reading attributes from HDF group '// &
                          group_name)
      ret_err = error
    end if
    call h5ltget_attribute_int_f(file_id, group_name, hdf_n_kappa, &
                                 iparam, error)
    if (error /= 0) then
      ! For backwards compatibility, don't crash if this attribute is
      ! not present. Instead just realise there is no internal layer
      ! data in the HDF file.
      nkap = 0
    else
      nkap = iparam(1)
    end if

    call this%thickness%read_hdf(group_id, hdf_thickness, error)
    if (error /= 0) then
      call logger%warning('ice_shelf%read_data','Error code '//        &
                          trim(str(error))//' returned when reading '// &
                          'ice shelf thickness field from HDF file')
      if (ret_err == 0) ret_err = error
    end if

    call this%velocity%read_hdf(group_id, hdf_velocity, error)
    if (error /= 0) then
      call logger%warning('ice_shelf%read_data','Error code '// &
                          trim(str(error))//' returned when '//  &
                          'reading ice shelf velocity field '//  &
                          'from HDF file')
      if (ret_err == 0) ret_err = error
    end if

    if (nkap > 0) allocate(this%kappa(nkap))
    do i = 1, nkap
      write(fieldname , hdf_kappa), i
      call this%kappa(i)%read_hdf(group_id, fieldname, error)
      if (error /= 0) then
        call logger%warning('ice_shelf%read_data','Error code '//        &
                            trim(str(error))//' returned when reading '// &
                            'ice shelf kappa field from HDF file')
        if (ret_err == 0) ret_err = error
      end if
    end do

    call h5gclose_f(group_id, error)
    if (error /= 0) then
      call logger%warning('ice_shelf%read_data','Error code '// &
                          trim(str(error))//' returned when '//  &
                          'closing HDF group '//group_name)
      if (ret_err == 0) ret_err = error
    end if
    error = ret_err
#ifdef DEBUG
    call logger%debug('ice_shelf%read_data','Read ice shelf data from '// &
                      'HDF group '//group_name)
#endif
  end subroutine shelf_read_data


  subroutine shelf_write_data(this,file_id,group_name,error)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Writes the state of the ice shelf object to an HDF file in the
    ! specified group. This will consist of a thickness and a velocity
    ! dataset.
    !
    class(ice_shelf), intent(in) :: this
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
    integer :: ret_err, i, nkap
    character(len=20) :: fieldname

    ret_err = 0
    call h5gcreate_f(file_id, group_name, group_id, error)
    if (error /= 0) then
      call logger%warning('ice_shelf%write_data','Error code '// &
                          trim(str(error))//' returned when '//  &
                          'creating HDF group "'//group_name//'"')
      call logger%error('ice_shelf%write_data','Data IO not performed for '// &
                        'ice shelf')
      return
    end if

    if (allocated(this%kappa)) then
      nkap = size(this%kappa)
    else
      nkap = 0
    end if

    call h5ltset_attribute_string_f(file_id, group_name, hdf_type_attr, &
                                    hdf_type_name, error)
    call h5ltset_attribute_double_f(file_id, group_name, hdf_lambda, &
                                    [this%lambda], 1_size_t, error)
    call h5ltset_attribute_double_f(file_id, group_name, hdf_chi, &
                                    [this%chi], 1_size_t, error)
    call h5ltset_attribute_double_f(file_id, group_name, hdf_zeta, &
                                    [this%zeta], 1_size_t, error)
    call h5ltset_attribute_int_f(file_id, group_name, hdf_n_kappa, &
                                 [nkap], 1_size_t, error)
    if (error /= 0) then
      call logger%warning('ice_shelf%write_data','Error code '// &
                          trim(str(error))//' returned when '//  &
                          'writing attribute to HDF group '//    &
                          group_name)
      ret_err = error
    end if

    call this%thickness%write_hdf(group_id, hdf_thickness, error)
    if (error /= 0) then
      call logger%warning('ice_shelf%write_data','Error code '//        &
                          trim(str(error))//' returned when writing '// &
                          'ice shelf thickness field to HDF file')
      if (ret_err == 0) ret_err = error
    end if

    call this%velocity%write_hdf(group_id, hdf_velocity, error)
    if (error /= 0) then
      call logger%warning('ice_shelf%write_data','Error code '// &
                          trim(str(error))//' returned when '//  &
                          'writing ice shelf velocity field '//  &
                          'to HDF file')
      if (ret_err == 0) ret_err = error
    end if

    do i = 1, nkap
      write(fieldname , hdf_kappa), i
      call this%kappa(i)%write_hdf(group_id, trim(fieldname), error)
      if (error /= 0) then
        call logger%warning('ice_shelf%write_data','Error code '// &
                            trim(str(error))//' returned when '//  &
                            'writing ice shelf kappa field '//  &
                            'to HDF file')
        if (ret_err == 0) ret_err = error
      end if
    end do

    call h5gclose_f(group_id, error)
    if (error /= 0) then
      call logger%warning('ice_shelf%write_data','Error code '// &
                          trim(str(error))//' returned when '//  &
                          'closing HDF group '//group_name)
      if (ret_err == 0) ret_err = error
    end if
    error = ret_err
    call logger%trivia('ice_shelf%write_data','Wrote ice shelf data to '// &
                       'HDF group '//group_name)
  end subroutine shelf_write_data


  function shelf_time_step(this) result(dt)
    !* Author: Chris MacMackin
    !  Date: December 2016
    !
    ! Calculates the time step for integrating the ice shelf, using
    ! the CFL condition.
    !
    class(ice_shelf), intent(in) :: this
    class(scalar_field), pointer :: u, dx1
    class(vector_field), pointer :: dx
    real(r8) :: dt
      !! The time-step to use
    u => this%velocity%component(1)
    dx => this%velocity%grid_spacing()
    dx1 => dx%component(1)
    call dx1%guard_temp()
    dt = min(minval(abs(this%courant*dx1/u)), this%max_dt)
    call dx1%clean_temp()
    call logger%trivia('ice_shelf%time_step','Calculated time step of '// &
                        trim(str(dt))//' using Courant number of '//       &
                        trim(str(this%courant)))
  end function shelf_time_step


  subroutine shelf_solve_velocity(this, basal_drag, success)
    !* Author: Chris MacMackin
    !  Date: May 2017
    !
    ! Computes the ice shelf velocity at the current time with the
    ! current ice thickness.
    !
    class(ice_shelf), intent(inout)          :: this
    class(scalar_field), intent(in)          :: basal_drag
      !! A paramter, e.g. coefficient of friction, needed to calculate
      !! the drag on basal surface of the glacier.
    logical, intent(out)                     :: success
      !! True if the integration is successful, false otherwise
    
    integer, save                             :: nval, kdmax = 20
    real(r8), dimension(:), allocatable       :: state
    integer, dimension(10)                    :: input
    integer, dimension(6)                     :: info
    real(r8), dimension(:), allocatable, save :: work
    real(r8), dimension(1)                    :: real_param
    integer, dimension(1)                     :: int_param
    integer                                   :: flag
 
    call basal_drag%guard_temp()
    nval = this%velocity%raw_size()
    if (allocated(work)) then
      if (size(work) < nval*(kdmax+5) + kdmax*(kdmax+3)) then
        deallocate(work)
        allocate(work(nval*(kdmax+5) + kdmax*(kdmax+3)))
      end if
    else
      allocate(work(nval*(kdmax+5) + kdmax*(kdmax+3)))
    end if

    input = 0
    input(4) = kdmax
    input(5) = 1
    input(9) = -1
    input(10) = 3

    state = this%velocity%raw()
    call nitsol(nval, state, nitsol_residual, nitsol_precondition, &
                1.e-10_r8*nval, 1.e-10_r8*nval, input, info, work, &
                real_param, int_param, flag, ddot, dnrm2)
    call this%velocity%set_from_raw(state)
    this%eta = this%viscosity_law%ice_viscosity(this%velocity, this%ice_temperature(), &
                                               this%time)
    this%stale_jacobian = .true.

    select case(flag)
    case(0)
      call logger%trivia('ice_shelf%solve_velocity','Found ice velocity at time '// &
                         trim(str(this%time)))
      success = .true.
    case(1)
      call logger%error('ice_shelf%solve_velocity','Reached maximum number of'// &
                        ' iterations finding ice velocity')
      success = .false.
    case default
      call logger%error('ice_shelf%solve_velocity','NITSOL failed when finding'// &
                        ' ice velocity with error code '//trim(str(flag)))
      success = .false.
    end select
    
    call basal_drag%clean_temp()

  contains
    
    subroutine nitsol_residual(n, xcur, fcur, rpar, ipar, itrmf)
      !! A routine matching the interface expected by NITSOL which
      !! returns the residual for the glacier.
      integer, intent(in)                   :: n
        !! Dimension of the problem
      real(r8), dimension(n), intent(in)    :: xcur
        !! Array of length `n` containing the current \(x\) value
      real(r8), dimension(n), intent(out)   :: fcur
        !! Array of length `n` containing f(xcur) on output
      real(r8), dimension(*), intent(inout) :: rpar
        !! Parameter/work array
      integer, dimension(*), intent(inout)  :: ipar
        !! Parameter/work array
      integer, intent(out)                  :: itrmf
        !! Termination flag. 0 means normal termination, 1 means
        !! failure to produce f(xcur)

      type(cheb1d_scalar_field) :: scalar_tmp
      integer :: start, finish, bounds_start, bounds_finish
      integer, dimension(:), allocatable :: lower, upper
      real(r8), dimension(:), allocatable :: bounds

      call this%velocity%set_from_raw(xcur)
      this%stale_jacobian = .true.
      associate(h => this%thickness, uvec => this%velocity, chi => this%chi, &
                zeta => this%zeta, eta => this%eta)
        eta = this%viscosity_law%ice_viscosity(uvec, this%ice_temperature(), &
                                               this%time)
        bounds = this%boundaries%boundary_residuals(h, uvec, eta, this%time)

        scalar_tmp = uvec%component(1)
        scalar_tmp = 4.0_r8*eta*h*scalar_tmp%d_dx(1)
        scalar_tmp = -2.0_r8*chi*h*h%d_dx(1) + scalar_tmp%d_dx(1)

        lower = this%boundaries%velocity_lower_bound()
        upper = this%boundaries%velocity_upper_bound()
        ! TODO: Figure out how to make this independent of order which
        ! values are stored in the field
        start = 1
        finish = start + this%velocity_upper_bound_size - 1
        bounds_start = this%thickness_lower_bound_size &
                     + this%thickness_upper_bound_size &
                     + this%velocity_lower_bound_size + 1
        bounds_finish = bounds_start + this%velocity_upper_bound_size - 1
        fcur(start:finish) = bounds(bounds_start:bounds_finish)
        start = finish + 1
        finish = start + scalar_tmp%raw_size(lower,upper) - 1
        fcur(start:finish) = scalar_tmp%raw(lower,upper)
        start = finish + 1
        finish = start + this%velocity_lower_bound_size - 1
        bounds_start = this%thickness_lower_bound_size &
                     + this%thickness_upper_bound_size + 1
        bounds_finish = bounds_start + this%velocity_lower_bound_size - 1
        fcur(start:finish) = bounds(bounds_start:bounds_finish)
      end associate
      itrmf = 0
    end subroutine nitsol_residual

    subroutine nitsol_precondition(n, xcur, fcur, ijob, v, z, rpar, ipar, itrmjv)
      !! A subroutine matching the interface expected by NITSOL, which
      !! acts as a preconditioner.
      integer, intent(in)                   :: n
        ! Dimension of the problem
      real(r8), dimension(n), intent(in)    :: xcur
        ! Array of lenght `n` containing the current $x$ value
      real(r8), dimension(n), intent(in)    :: fcur
        ! Array of lenght `n` containing the current \(f(x)\) value
      integer, intent(in)                   :: ijob
        ! Integer flat indicating which product is desired. 0
        ! indicates \(z = J\vec{v}\). 1 indicates \(z = P^{-1}\vec{v}\).
      real(r8), dimension(n), intent(in)    :: v
        ! An array of length `n` to be multiplied
      real(r8), dimension(n), intent(out)   :: z
        ! An array of length n containing the desired product on
        ! output.
      real(r8), dimension(*), intent(inout) :: rpar
        ! Parameter/work array 
      integer, dimension(*), intent(inout)  :: ipar
        ! Parameter/work array
      integer, intent(out)                  :: itrmjv
        ! Termination flag. 0 indcates normal termination, 1
        ! indicatesfailure to prodce $J\vec{v}$, and 2 indicates
        ! failure to produce \(P^{-1}\vec{v}\)

      type(cheb1d_scalar_field) :: delta_u, tmp
      integer :: i, sl, el, su, eu
      integer, dimension(2) :: upper_type, lower_type
      integer, dimension(:), allocatable :: boundary_types, boundary_locations
      real(r8) :: eta_val

      if (ijob /= 1) then
        itrmjv = 0
        return
      end if

      call delta_u%assign_meta_data(this%velocity)
      call delta_u%set_from_raw(v)
      
      sl = this%velocity_size - this%velocity_lower_bound_size + 1
      el = this%velocity_size
      su = 1
      eu = this%velocity_upper_bound_size
      boundary_locations = [(i, i=sl,el), (i, i=su,eu)]
      upper_type = this%boundaries%velocity_upper_type()
      lower_type = this%boundaries%velocity_lower_type()
      boundary_types = [(lower_type(1), i=sl,el), (upper_type(1), i=su,eu)]
      associate(h => this%thickness, chi => this%chi, zeta => this%zeta, &
                jac => this%velocity_jacobian, eta => this%eta)
        if (this%stale_jacobian) then
          ! Mathematically, it should be 4._r8*eta*h which I pass, but
          ! this sometimes results in an ill-conditioned Jacobian (for
          ! reasons I'm not clear on). If eta ~ 1 then turns out I get
          ! good results just ignoring it.
          select type(visc => this%viscosity_law)
          class is(newtonian_viscosity)
            eta_val = eta%get_element(1)
            if (eta_val <= 0._r8) eta_val = 1._r8
            jac = jacobian_block(4._r8*eta*h, 1, 1, boundary_locs=boundary_locations, &
                                 boundary_types=boundary_types)
          class default
            jac = jacobian_block(4._r8*h, 1, 1, boundary_locs=boundary_locations, &
                                 boundary_types=boundary_types)
          end select
          this%stale_jacobian = .false.
        end if
        delta_u = jac%solve_for(delta_u)
        z(1:n) = delta_u%raw()
      end associate
      itrmjv = 0
    end subroutine nitsol_precondition

  end subroutine shelf_solve_velocity


  subroutine shelf_integrate(this, old_states, basal_melt, basal_drag, &
                             water_density, time, success)
    !* Author: Chris MacMackin
    !  Date: May 2017
    !
    ! Integrates the glacier's state forward to `time`. This is done
    ! using an explicit method for the thickness and a Newton's solver
    ! for velocity.
    !
    class(ice_shelf), intent(inout)          :: this
    class(glacier), dimension(:), intent(in) :: old_states
      !! Previous states of the glacier, with the most recent one
      !! first.
    class(scalar_field), intent(in)          :: basal_melt
      !! The melt rate that the bottom of the glacier experiences
      !! during this time step.
    class(scalar_field), intent(in)          :: basal_drag
      !! A paramter, e.g. coefficient of friction, needed to calculate
      !! the drag on basal surface of the glacier.
    real(r8), intent(in)                     :: water_density
      !! The density of the water below the glacier.
    real(r8), intent(in)                     :: time
      !! The time to which the glacier should be integrated
    logical, intent(out)                     :: success
      !! True if the integration is successful, false otherwise
    
    associate(h => this%thickness, uvec => this%velocity, m => basal_melt, &
              lambda => this%lambda, chi => this%chi, zeta => this%zeta,   &
              t_old => this%time)
      this%thickness = this%thickness - (time - t_old)*(lambda*m + .div.(h*uvec))
    end associate
  end subroutine shelf_integrate

  subroutine shelf_assign(this, rhs)
    !* Author: Chris MacMackin
    !  Date: February 2017
    !
    ! Copies the data from one ice shelf into another. This is only
    ! needed due to a bug in gfortran which means that the intrinsic
    ! assignment for glacier types is not using the appropriate
    ! defined assignment for the field components.
    !
    ! It does not assign the Jacobian object as it would take up quite
    ! a bit of extra space and it is unlikely that it would ever be
    ! needed without first having to be recalculated.
    !
    class(ice_shelf), intent(out) :: this
    class(glacier), intent(in)    :: rhs
      !! The ice shelf to be assigned to this one.
    select type(rhs)
    class is(ice_shelf)
      this%thickness = rhs%thickness
      this%velocity = rhs%velocity
      this%eta = rhs%eta
      this%lambda = rhs%lambda
      this%chi = rhs%chi
      this%zeta = rhs%zeta
      this%courant = rhs%courant
      this%max_dt = rhs%max_dt
      allocate(this%viscosity_law, source=rhs%viscosity_law)
      allocate(this%boundaries, source=rhs%boundaries)
      this%time = rhs%time
      this%thickness_size = rhs%thickness_size
      this%velocity_size = rhs%velocity_size
      this%boundary_start = rhs%boundary_start
      this%thickness_lower_bound_size = rhs%thickness_lower_bound_size
      this%thickness_upper_bound_size = rhs%thickness_upper_bound_size
      this%velocity_lower_bound_size = rhs%velocity_lower_bound_size
      this%velocity_upper_bound_size = rhs%velocity_upper_bound_size
      this%stale_jacobian = .true.
      this%stale_eta = .true.
      if (allocated(rhs%kappa)) then
        allocate(this%kappa(size(rhs%kappa)))
        this%kappa = rhs%kappa
      end if
    class default
      call logger%fatal('ice_shelf%assign','Type other than `ice_shelf` '// &
                        'requested to be assigned.')
      error stop
    end select
#ifdef DEBUG
    call logger%debug('ice_shelf%assign','Copied ice shelf data.')
#endif
  end subroutine shelf_assign


  subroutine ice_shelf_integrate_layers(this, old_states, time, success)
    !* Author: Chris MacMackin
    !  Date: September 2018
    !
    ! Integrate the Taylor coefficients representing the vertical
    ! structure of internal reflectors forward to the specified
    ! time. This is done using an implicit method, with the resulting
    ! linear system solved using GMRES.
    !
    class(ice_shelf), intent(inout)          :: this
    class(glacier), dimension(:), intent(in) :: old_states
      !! Previous states of the ice_shelf, with the most recent one
      !! first.
    real(r8), intent(in)                     :: time
      !! The time to which the ice_shelf should be integrated
    logical, intent(out)                     :: success
      !! True if the integration is successful, false otherwise

    integer :: n, flag
    real(r8), dimension(:), allocatable :: solution
    real(r8) :: dt, resid
    type(jacobian_block) :: precond_block

    if (.not. allocated(this%kappa)) return
    select type(old_states)
    class is(ice_shelf)
      dt = time - old_states(1)%time
      do n = 1, size(this%kappa)
        solution = this%kappa(n)%raw()
        precond_block = jacobian_block(dt*this%velocity%component(1), 1, &
                                       boundary_locs=[this%thickness_size], &
                                       boundary_types=[dirichlet], &
                                       coef=real(-n, r8)) + 1._r8
        call gmres_solve(solution, operator, old_states(1)%kappa(n)%raw(), resid, &
                         flag, precond=preconditioner, krylov_dim=40)
        call this%kappa(n)%set_from_raw(solution)
        if (flag /= 0) then
          call logger%error('ice_shelf%integrate_layers','GMRES solver '// &
                            'returned with error code '//str(flag))
          success = .false.
          return
        end if
        call this%kappa(n)%set_from_raw(solution)
      end do
    class default
      call logger%fatal('ice_shelf%integrate_layers','Type other than `ice_shelf` '// &
                        'passed to `ice_shelf` object as a previous state.')
      error stop
    end select

  contains

    function operator(v, xcur, rhs, rpar, ipar, success)
      real(r8), dimension(:), intent(in)    :: v
        !! The vector to be multiplied
      real(r8), dimension(:), intent(in)    :: xcur
        !! Array containing the current estimate of the independent
        !! variables in the linear system. This may not be needed, but
        !! is provided just in case.
      real(r8), dimension(:), intent(in)    :: rhs
        !! Array containing the right hand side of the linear
        !! system. This may not be needed, but is provided just in
        !! case.
      real(r8), dimension(*), intent(inout) :: rpar
        !! Parameter/work array 
      integer, dimension(*), intent(inout)  :: ipar
        !! Parameter/work array
      logical, intent(out)                  :: success
        !! Indicates whether operation was completed succesfully
      real(r8), dimension(size(xcur))       :: operator
        !! Result of the operation

      type(cheb1d_scalar_field) :: kappa
      class(scalar_field), pointer :: tmp
      
      call kappa%assign_meta_data(this%kappa(n))
      call kappa%set_from_raw(v)
      tmp => dt*this%velocity .dot. (.grad. kappa)
      kappa = (1._r8 - n*dt*(.div.this%velocity))*kappa + tmp
      operator = kappa%raw()
      operator(this%thickness_size) = v(this%thickness_size)
      success = .true.
    end function operator


    function preconditioner(v, xcur, rhs, rpar, ipar, success)
      real(r8), dimension(:), intent(in)    :: v
        !! The vector to be multiplied
      real(r8), dimension(:), intent(in)    :: xcur
        !! Array containing the current estimate of the independent
        !! variables in the linear system. This may not be needed, but
        !! is provided just in case.
      real(r8), dimension(:), intent(in)    :: rhs
        !! Array containing the right hand side of the linear
        !! system. This may not be needed, but is provided just in
        !! case.
      real(r8), dimension(*), intent(inout) :: rpar
        !! Parameter/work array 
      integer, dimension(*), intent(inout)  :: ipar
        !! Parameter/work array
      logical, intent(out)                  :: success
        !! Indicates whether operation was completed succesfully
      real(r8), dimension(size(xcur))       :: preconditioner
        !! Result of the operation

      type(cheb1d_scalar_field) :: tmp
      
      call tmp%assign_meta_data(this%kappa(n))
      call tmp%set_from_raw(v)
      tmp = precond_block%solve_for(tmp)
      preconditioner = tmp%raw()
      success = .true.
    end function preconditioner
  end subroutine ice_shelf_integrate_layers


end module ice_shelf_mod

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
  use preconditioner_mod, only: preconditioner
  use hdf5
  use h5lt
  use logger_mod, only: logger => master_logger
  use penf, only: str
  implicit none
  private

  character(len=9), parameter, public :: hdf_type_name = 'ice_shelf'
  character(len=9), parameter, public :: hdf_thickness = 'thickness'
  character(len=8), parameter, public :: hdf_velocity = 'velocity'

  type, extends(glacier), public :: ice_shelf
    !* Author: Chris MacMackin
    !  Date: April 2016
    !
    ! A concrete implementation of the [[glacier]] type, using a vertically
    ! integrated model of an ice shelf. This model is 1-dimensional only.
    !
    private
    type(cheb1d_scalar_field) :: thickness 
      !! Thickness of ice shelf, $h$
    type(cheb1d_vector_field) :: velocity  
      !! Flow velocity of ice shelf, $\vec{u}$
    real(r8)                  :: lambda
      !! The dimensionless ratio 
      !! $\lambda \equiv \frac{\rho_0m_0x_0}{\rho_iH-0u_0}$
    real(r8)                  :: chi
      !! The dimensionless ratio
      !! $\chi \equiv \frac{\rho_igh_0x_x}{2\eta_0u_0}$
    class(abstract_viscosity), allocatable :: viscosity_law
      !! An object representing the model used for ice viscosity.
    class(glacier_boundary), allocatable   :: boundaries
      !! An object specifying the boundary conditions for the ice
      !! shelf.
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
    real(r8)                  :: jacobian_time
      !! The time at which the Jacobian was last updated.
    type(jacobian_block), dimension(2,2) :: jacobian
      !! A representation of the Jacobian for this ice shelf.
    type(preconditioner)      :: precondition_obj
      !! An object with a method to apply a block-Jacobian
      !! preconditioner, with the specified convergence properties.
  contains
!$    procedure            :: t => shelf_dt
!$    procedure            :: local_error => shelf_local_error
!$    procedure            :: integrand_multiply_integrand => shelf_m_shelf
!$    procedure            :: integrand_multiply_real => shelf_m_real
!$    procedure, pass(rhs) :: real_multiply_integrand => real_m_shelf
!$    procedure            :: add => shelf_add
!$    procedure            :: sub => shelf_sub
!$    procedure            :: assign_integrand => shelf_assign
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
    procedure :: write_data => shelf_write_data
    procedure :: time_step => shelf_time_step
    procedure, private :: assign => shelf_assign
  end type ice_shelf

!  interface ice_shelf
!    module procedure constructor
!  end interface ice_shelf

contains
  
  subroutine shelf_initialise(this, domain, resolution, thickness, velocity, &
                              temperature, viscosity_law, boundaries,  &
                              lambda, chi)
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
      !! $\lambda \equiv \frac{\rho_0m_0x_0}{\rho_iH_0u_0}$.
    real(r8), intent(in), optional       :: chi
      !! The dimensionless ratio
      !! $\chi \equiv \frac{\rho_igh_0x_x}{2\eta_0u_0}\left(1 -
      !! \frac{\rho_i}{\rho_0}\right)$.

    integer, dimension(:), allocatable :: lower, upper

    this%thickness = cheb1d_scalar_field(resolution(1),thickness,domain(1,1),domain(1,2))
    this%velocity = cheb1d_vector_field(resolution(1),velocity,domain(1,1),domain(1,2))
    this%thickness_size = this%thickness%raw_size()
    this%velocity_size = this%velocity%raw_size()

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
    this%jacobian_time = -1._r8
    this%precondition_obj = preconditioner(1.e-3_r8, 10)
#ifdef DEBUG
    call logger%debug('ice_shelf','Instantiated new ice shelf object')
#endif
  end subroutine shelf_initialise

!$  function shelf_dt(self,t)
!$    !* Author: Christopher MacMackin
!$    !  Date: April 2016
!$    !
!$    ! Computes the derivative of the ice shelf with respect to time. As
!$    ! the only property of the ice which explicitely changes with time is
!$    ! the ice thickness, that will be the only portion of the returned type
!$    ! which actually corresponds to the derivative.
!$    !
!$    class(ice_shelf), intent(in)   :: self
!$    real(r8), intent(in), optional :: t
!$      !! Time at which to evaluate the derivative
  !$    class(integrand), allocatable  :: shelf_dt
!$      !! The time rate of change of the ice shelf. Has dynamic type
!$      !! [[ice_shelf]].
!$  end function shelf_dt
!$
!$  function shelf_local_error(lhs, rhs) result(error)
!$    !* Author: Christopher MacMackin
!$    !  Date: April 2016
!$    !
!$    ! Calculates a real scalar to represent the absolute difference between
!$    ! two ice_shelf objects. `rhs` must be a a [[ice_shelf]] object, or a
!$    ! runtime error will occur.
!$    !
!$    class(ice_shelf), intent(in) :: lhs
!$      !! Self
!$    class(integrand), intent(in) :: rhs
!$      !! The ice shelf object which is being compared against.
!$    real(r8) :: error
!$      !! The scalar representation of the absolute difference between these
!$      !! two ice shelves.
!$  end function shelf_local_error
!$
!$  function shelf_m_shelf(lhs, rhs) result(product)
!$    !* Author: Christopher MacMackin
!$    !  Date: April 2016
!$    !
!$    ! Multiplies one ice shelf object by another. That is to say, it 
!$    ! performs element-wise multiplication of the state vectors 
!$    ! representing the two arguments. `rhs` must be an [[ice_shelf]]
!$    ! object, or a runtime error will occur.
!$    !
!$    class(ice_shelf), intent(in)  :: lhs
!$      !! Self
!$    class(integrand), intent(in)  :: rhs
!$      !! The ice shelf object being multiplied by.
!$    class(integrand), allocatable :: product
!$      !! The product of the two arguments. Has dynamic type [[ice_shelf]].
!$  end function shelf_m_shelf
!$
!$  function shelf_m_real(lhs, rhs) result(product)
!$    !* Author: Christopher MacMackin
!$    !  Date: April 2016
!$    !
!$    ! Multiplies one ice shelf object by a scalar. That is to say, it 
!$    ! performs element-wise multiplication of the state vector 
!$    ! representing the ice shelf.
!$    !
!$    class(ice_shelf), intent(in)  :: lhs
!$      !! Self
!$    real(r8), intent(in)          :: rhs
!$      !! The scalar being multiplied by.
!$    class(integrand), allocatable :: product
!$      !! The product of the two arguments. Has dynamic type [[ice_shelf]].
!$  end function shelf_m_real
!$
!$  function real_m_shelf(lhs, rhs) result(product)
!$    !* Author: Christopher MacMackin
!$    !  Date: April 2016
!$    !
!$    ! Multiplies one ice shelf object by a scalar. That is to say, it 
!$    ! performs element-wise multiplication of the state vector 
!$    ! representing the ice shelf.
!$    !
!$    real(r8), intent(in)          :: lhs
!$      !! The scalar being multiplied by.
!$    class(ice_shelf), intent(in)  :: rhs
!$      !! Self
!$    class(integrand), allocatable :: product
!$      !! The product of the two arguments. Has dynamic type [[ice_shelf]].
!$  end function real_m_shelf
!$
!$  function shelf_add(lhs, rhs) result(sum)
!$    !* Author: Christopher MacMackin
!$    !  Date: April 2016
!$    !
!$    ! Adds one ice shelf object to another. That is to say, it performs
!$    ! element-wise addition of the state vectors representing the two
!$    ! arguments. `rhs` must be an [[ice_shelf]] object, or a runtime
!$    ! error will occur.
!$    !
!$    class(ice_shelf), intent(in)  :: lhs
!$      !! Self
!$    class(integrand), intent(in)  :: rhs
!$      !! The ice shelf object being added.
!$    class(integrand), allocatable :: sum
!$      !! The sum of the two arguments. Has dynamic type [[ice_shelf]].
!$  end function shelf_add
!$
!$  function shelf_sub(lhs, rhs) result(difference)
!$    !* Author: Christopher MacMackin
!$    !  Date: April 2016
!$    !
!$    ! Subtracts one ice shelf object from another. That is to say, it 
!$    ! performs element-wise addition of the state vectors representing 
!$    ! the two arguments. `rhs` must be a a [[ice_shelf]] object, or a
!$    ! runtime error will occur.
!$    !
!$    class(ice_shelf), intent(in)  :: lhs
!$      !! Self
!$    class(integrand), intent(in)  :: rhs
!$      !! The ice shelf object being subtracted.
!$    class(integrand), allocatable :: difference
!$      !! The difference of the two arguments. Has dynamic type [[ice_shelf]].
!$  end function shelf_sub
!$
!$  subroutine shelf_assign(lhs, rhs)
!$    !* Author: Christopher MacMackin
!$    !  Date: April 2016
!$    !
!$    ! Assigns the `rhs` ice shelf to this, `lhs`, one. All components
!$    ! will be the same following the assignment.
!$    !
!$    class(ice_shelf), intent(inout) :: lhs
!$      !! Self
!$    class(integrand), intent(in)    :: rhs
!$      !! The object to be assigned. Must have dynamic type [[ice_shelf]],
!$      !! or a runtime error will occur.
!$  end subroutine shelf_assign

  pure function shelf_thickness(this) result(thickness)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Returns the thickness of the ice shelf across its domain.
    !
    class(ice_shelf), intent(in)     :: this
    class(scalar_field), allocatable :: thickness !! The ice thickness.
    allocate(thickness, source=this%thickness)
#ifdef DEBUG
    call logger%debug('ice_shelf%thickness','Returned ice shelf thickness')
#endif    
  end function shelf_thickness


  pure function shelf_velocity(this) result(velocity)
    !* Author: Christopher MacMackin
    !  Date: July 2016
    !
    ! Returns the velocity of the ice shelf across its domain.
    !
    class(ice_shelf), intent(in)     :: this
    class(vector_field), allocatable :: velocity !! The ice velocity.
    allocate(velocity, source=this%velocity)
    call velocity%set_temp()
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
    call logger%debug('ice_shelf%density','Ice shelf has '// &
                      str(density))
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
                      str(temperature))
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

    call melt_rate%guard_temp(); call basal_drag_parameter%guard_temp()
    allocate(residual(this%data_size()))
    start = 1

    ! Use same or similar notation for variables as in equations
    select type(previous_states)
    class is(ice_shelf)
      ! TODO: Either move the function results over to proper
      ! assignment or figure out some way to have temporary results in
      ! associate constructs.
      associate(h => this%thickness, h_old => previous_states(1)%thickness, &
                uvec => this%velocity, u => this%velocity%component(1), &
                eta => this%viscosity_law%ice_viscosity(this%velocity, &
                                   this%ice_temperature(), this%time), &
                m => melt_rate, lambda => this%lambda, chi => this%chi, &
                t_old => previous_states(1)%time)
        ! Boundary conditions
        bounds = this%boundaries%boundary_residuals(h, uvec, eta, this%time)

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
  
        ! Momentum equation, x-component
        scalar_tmp = 4.0_r8*eta*h*u%d_dx(1)
        scalar_tmp = scalar_tmp%d_dx(1) - 2.0_r8*chi*h*h%d_dx(1)
        
        lower = this%boundaries%velocity_lower_bound()
        upper = this%boundaries%velocity_upper_bound()
        ! TODO: Figure out how to make this independent of order which
        ! values are stored in the field
        finish = start + this%velocity_upper_bound_size - 1
        bounds_start = this%thickness_lower_bound_size &
                     + this%thickness_upper_bound_size &
                     + this%velocity_lower_bound_size + 1
        bounds_finish = bounds_start + this%velocity_upper_bound_size - 1
        residual(start:finish) = bounds(bounds_start:bounds_finish)
        start = finish + 1
        finish = start + scalar_tmp%raw_size(lower,upper) - 1
        residual(start:finish) = scalar_tmp%raw(lower,upper)
        start = finish + 1
        finish = start + this%velocity_lower_bound_size - 1
        bounds_start = this%thickness_lower_bound_size &
                     + this%thickness_upper_bound_size + 1
        bounds_finish = bounds_start + this%velocity_lower_bound_size - 1
        residual(start:finish) = bounds(bounds_start:bounds_finish)
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
    integer :: i
    !TODO: Add some assertion-like checks that the state vector is the right size
    call this%thickness%set_from_raw(state_vector(1:this%thickness_size))
    i = 1 + this%thickness_size
    call this%velocity%set_from_raw(state_vector(i:i + this%velocity_size - 1))
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

    type(ice_shelf) :: delta_shelf
    type(cheb1d_scalar_field), dimension(2) :: vector, estimate
    class(scalar_field), allocatable :: eta, u
    integer :: sl, el, su, eu
    integer :: i, elem
    integer, dimension(2) :: upper_type, lower_type
    integer, dimension(:), allocatable :: boundary_types, boundary_locations
    real(r8) :: delta_t

    call melt_rate%guard_temp(); call basal_drag_parameter%guard_temp()
    allocate(preconditioned(size(delta_state)))
    select type(previous_states)
    class is(ice_shelf)
      delta_t = this%time - previous_states(1)%time
    class default
      error stop ('Type other than `ice_shelf` passed to `ice_shelf` '// &
                  'object as a previous state.')
    end select
    call delta_shelf%thickness%assign_meta_data(this%thickness)
    call delta_shelf%velocity%assign_meta_data(this%velocity)
    delta_shelf%thickness_size = this%thickness_size
    delta_shelf%velocity_size = this%velocity_size
    call delta_shelf%update(delta_state)
    allocate(eta, source=this%viscosity_law%ice_viscosity(this%velocity, &
                                       this%ice_temperature(), this%time))
    call eta%unset_temp()
    vector(1) = delta_shelf%thickness
    vector(2) = delta_shelf%velocity%component(1)
    call this%velocity%allocate_scalar_field(u)
    u = this%velocity%component(1)

    associate(h => this%thickness, chi => this%chi)
      if (this%jacobian_time < this%time) then
        sl = this%thickness_size - this%thickness_lower_bound_size + 1
        el = this%thickness_size
        su = 1
        eu = this%thickness_upper_bound_size
        boundary_locations = [(i, i=sl,el), (i, i=su,eu)]
        upper_type = this%boundaries%thickness_upper_type()
        lower_type = this%boundaries%thickness_lower_type()
        boundary_types = [(lower_type(1), i=sl,el), (upper_type(1), i=su,eu)]

        this%jacobian(1,1) = jacobian_block(u, 1, boundary_locs=boundary_locations) &
                           + 1._r8/delta_t
        this%jacobian(1,2) = jacobian_block(h, 1, boundary_locs=boundary_locations, &
                                            boundary_types=boundary_types)

        sl = this%velocity_size - this%velocity_lower_bound_size + 1
        el = this%velocity_size
        su = 1
        eu = this%velocity_upper_bound_size
        boundary_locations = [(i, i=sl,el), (i, i=su,eu)]
        upper_type = this%boundaries%velocity_upper_type()
        lower_type = this%boundaries%velocity_lower_type()
        boundary_types = [(lower_type(1), i=sl,el), (upper_type(1), i=su,eu)]

        this%jacobian(2,1) = jacobian_block(4._r8*eta*u%d_dx(1) - 2._r8*chi*h, 1, &
                                            boundary_locs=boundary_locations,     &
                                            boundary_operations=jacobian_bounds_2_1)
        this%jacobian(2,2) = jacobian_block(4._r8*eta*h, 1, 1,                &
                                            boundary_locs=boundary_locations, &
                                            boundary_types=boundary_types)
        this%jacobian_time = this%time
      end if

      elem = h%elements()
      do i = 1, size(estimate)
        estimate(i) = cheb1d_scalar_field(elem)
        call estimate(i)%assign_meta_data(vector(i))
      end do
      call this%precondition_obj%apply(this%jacobian, vector, estimate)
    end associate
    
    preconditioned = [(estimate(i)%raw(), i=1,size(estimate))]

    call melt_rate%clean_temp(); call basal_drag_parameter%clean_temp()
#ifdef DEBUG
    call logger%debug('ice_shelf%precondition','Applied preconditioner for ice shelf.')
#endif

  contains
    
    subroutine jacobian_bounds_2_1(contents, derivative, rhs,     &
                                   boundary_locs, boundary_types, &
                                   boundary_values)
      !* Author: Chris MacMackin
      !  Date: January 2016
      !
      ! Specifies the boundary values for the result multiplying the
      ! thickness field by the bottom-left Jacobian block.
      !
      class(scalar_field), intent(in)                  :: contents
        !! The field used to construct the Jacobian block
      class(scalar_field), intent(in)                  :: derivative
        !! The first spatial derivative of the field used to construct
        !! the Jacobian block, in the direction specified
      class(scalar_field), intent(in)                  :: rhs
        !! The scalar field representing the vector being multiplied
        !! by Jacobian
      integer, dimension(:), allocatable, intent(in)   :: boundary_locs
        !! The locations in the raw representation of `rhs` containing
        !! the boundaries.
      integer, dimension(:), allocatable, intent(in)   :: boundary_types
        !! Integers specifying the type of boundary condition. The type
        !! of boundary condition corresponding to a given integer is
        !! specified in [[boundary_types_mod]]. Only Dirichlet and
        !! Neumann conditions are supported. The storage order must
        !! correspond to that of `boundary_locs`.
      real(r8), dimension(:), allocatable, intent(out) :: boundary_values
        !! The values to go at the boundaries when multiplying a field
        !! by the Jacobian block. The storage order must be the same as
        !! for `boundary_locs`.
      integer :: n
      n = size(boundary_locs)
      allocate(boundary_values(n))
      !TODO: Make this general
      boundary_values = [0._r8,-this%thickness%get_element(boundary_locs(n))/ &
                         eta%get_element(boundary_locs(n))*.25_r8*this%chi]
    end subroutine jacobian_bounds_2_1
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
                      str(time))
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
    shelf_data_size = this%thickness_size + this%velocity_size
#ifdef DEBUG
    call logger%debug('ice_shelf%data_size','Ice shelf has '//                &
                      str(shelf_data_size)//' elements in its state vector.')
#endif
  end function shelf_data_size


  pure function shelf_state_vector(this) result(state_vector) 
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Returns the state vector for the current state of the ice shelf. 
    ! This takes the form of a 1D array.
    !
    class(ice_shelf), intent(in)        :: this
    real(r8), dimension(:), allocatable :: state_vector
      !! The state vector describing the ice shelf.
    state_vector = [this%thickness%raw(),this%velocity%raw()]
#ifdef DEBUG
    call logger%debug('ice_shelf%state_vector','Returning state vector '// &
                      'for ice shelf.')
#endif
  end function shelf_state_vector


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
    integer :: ret_err

    ret_err = 0
    call h5gcreate_f(file_id, group_name, group_id, error)
    if (error /= 0) then
      call logger%warning('ice_shelf%write_data', 'Error code '//str(error)// &
                          ' returned when creating HDF group "'//group_name//'"')
      call logger%error('ice_shelf%write_data','Data IO not performed for '// &
                        'ice shelf')
      return
    end if

    call h5ltset_attribute_string_f(file_id, group_name, hdf_type_attr, &
                                    hdf_type_name, error)
    if (error /= 0) then
      call logger%warning('ice_shelf%write_data','Error code'//str(error)// &
                          ' returned when writing attribute to HDF group'// &
                          group_name)
      ret_err = error
    end if

    call this%thickness%write_hdf(group_id, hdf_thickness, error)
    if (error /= 0) then
      call logger%warning('ice_shelf%write_data','Error code '//str(error)// &
                          ' returned when writing ice shelf thickness '//    &
                          'field to HDF file')
      if (ret_err == 0) ret_err = error
    end if

    call this%velocity%write_hdf(group_id, hdf_velocity, error)
    if (error /= 0) then
      call logger%warning('ice_shelf%write_data','Error code '//str(error)// &
                          ' returned when writing ice shelf velocity '//     &
                          'field to HDF file')
      if (ret_err == 0) ret_err = error
    end if

    call h5gclose_f(group_id, error)
    if (error /= 0) then
      call logger%warning('ice_shelf%write_data','Error code'//str(error)// &
                          ' returned when closing HDF group'//group_name)
      if (ret_err == 0) ret_err = error
    end if
    error = ret_err
#ifdef DEBUG
    call logger%debug('ice_shelf%write_data','Wrote ice shelf data to '// &
                      'HDF group'//group_name)
#endif
  end subroutine shelf_write_data

  function shelf_time_step(this) result(dt)
    !* Author: Chris MacMackin
    !  Date: December 2016
    !
    ! Calculates the time step for integrating the ice shelf, using
    ! the CFL condition.
    !
    class(ice_shelf), intent(in) :: this
    real(r8) :: dt
      !! The time-step to use
    associate(u => this%velocity%component(1), dx => this%velocity%grid_spacing(), &
              C => 10.0_r8)
      associate(dx1 => dx%component(1))
        dt = minval(abs(C*dx1/u))
        call logger%trivia('ice_shelf%time_step','Calculated time step of '// &
                           str(dt)//' using Courant number of '//str(C))
      end associate
    end associate
  end function shelf_time_step

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
      this%lambda = rhs%lambda
      this%chi = rhs%chi
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
      this%precondition_obj = rhs%precondition_obj
    class default
      call logger%fatal('ice_shelf%assign','Type other than `ice_shelf` '// &
                        'requested to be assigned.')
      error stop
    end select
#ifdef DEBUG
    call logger%debug('ice_shelf%assign','Copied ice shelf data.')
#endif
  end subroutine shelf_assign

end module ice_shelf_mod

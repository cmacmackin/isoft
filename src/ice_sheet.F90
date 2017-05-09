!
!  ice_sheet.f90
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

module ice_sheet_mod
  !* Author: Christopher MacMackin
  !  Date: April 2016
  !  License: GPLv3
  !
  ! Provides a concrete implementation of the [[glacier]] type, using
  ! a vertically integrated model of an ice sheet.
  !
  use iso_fortran_env, only: r8 => real64
  !use foodie, only: integrand
  use glacier_mod, only: glacier, thickness_func, velocity_func
  use factual_mod, only: scalar_field, vector_field, cheb1d_scalar_field, &
                         cheb1d_vector_field, maxval
  use viscosity_mod, only: abstract_viscosity
  use jacobian_block_mod, only: jacobian_block
  use hdf5
  use logger_mod, only: logger => master_logger
  implicit none
  private

  type, extends(glacier), public :: ice_sheet
    !* Author: Chris MacMackin
    !  Date: April 2016
    !
    ! A concrete implementation of the [[glacier]] type, using a vertically
    ! integrated model of an ice sheet. This model is 1-dimensional only.
    !
    private
    type(cheb1d_scalar_field) :: thickness
      !! Thickness of ice sheet, $h$
    type(cheb1d_vector_field) :: velocity
      !! Flow velocity of ice sheet, $\vec{u}$
    real(r8)                  :: lambda
      !! The dimensionless ratio 
      !! $\lambda \equiv \frac{\rho_0m_0x_0}{\rho_iH-0u_0}$
    real(r8)                  :: chi
      !! The dimensionless ratio
      !! $\chi \equiv \frac{\rho_igh_0x_x}{2\eta_0u_0}$
    class(abstract_viscosity), allocatable :: viscosity_law
      !! An object representing the model used for ice viscosity.
    real(r8)                  :: time
      !! The time at which the ice sheet is in this state
  contains
!$    procedure            :: t => sheet_dt
!$    procedure            :: local_error => sheet_local_error
!$    procedure            :: integrand_multiply_integrand => sheet_m_sheet
!$    procedure            :: integrand_multiply_real => sheet_m_real
!$    procedure, pass(rhs) :: real_multiply_integrand => real_m_sheet
!$    procedure            :: add => sheet_add
!$    procedure            :: sub => sheet_sub
!$    procedure            :: assign_integrand => sheet_assign
    procedure :: ice_thickness => sheet_thickness
!$    procedure :: ice_velocity => sheet_velocity
    procedure :: ice_density => sheet_density
    procedure :: ice_temperature => sheet_temperature
    procedure :: residual => sheet_residual
    procedure :: update => sheet_update
    procedure :: precondition =>  sheet_precondition
    procedure :: set_time => sheet_set_time
    procedure :: data_size => sheet_data_size
    procedure :: state_vector => sheet_state_vector
    procedure :: read_data => sheet_read_data
    procedure :: write_data => sheet_write_data
    procedure :: time_step => sheet_time_step
    procedure, private :: assign => sheet_assign
  end type ice_sheet

  interface ice_sheet
    module procedure constructor
  end interface ice_sheet

contains
  
  function constructor(domain, resolution, thickness, velocity, &
                       viscosity_law, lambda, chi) result(this)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Creates a new [[ice_sheet]] object with initial conditions provided
    ! by the arguments. At present only a 1D model is supported. If
    ! information is provided for higher dimensions then it will be ignored.
    !
    real(r8), dimension(:,:), intent(in)            :: domain
      !! An array containing the upper and lower limits of the domain for
      !! the ice sheet. The first index represents the dimension for which
      !! the boundaries apply. If the second index is 1 then it corresponds
      !! to the lower bound. If the second index is 2 then it corresponds to
      !! the upper bound.
    integer, dimension(:), intent(in)               :: resolution
      !! The number of data points in each dimension
    procedure(thickness_func)                       :: thickness
      !! A function which calculates the initial value of the thickness of 
      !! the ice sheet at a given location.
    procedure(velocity_func)                        :: velocity
      !! A function which calculates the initial value of the velocity 
      !! (vector) of the ice at a given location in an ice sheet.
    class(abstract_viscosity), intent(in), optional :: viscosity_law
      !! An object which calculates the viscosity of the ice.
    real(r8), intent(in), optional                  :: lambda
      !! The dimensionless ratio 
      !! $\lambda \equiv \frac{\rho_0m_0x_0}{\rho_iH-0u_0}$
    real(r8), intent(in), optional                  :: chi
      !! The dimensionless ratio
      !! $\chi \equiv \frac{\rho_igh_0x_x}{2\eta_0u_0}$
    type(ice_sheet)                                 :: this
      !! An ice sheet object with its domain and initial conditions set
      !! according to the arguments of the constructor function.
  end function constructor
  
!$  function sheet_dt(self,t)
!$    !* Author: Christopher MacMackin
!$    !  Date: April 2016
!$    !
!$    ! Computes the derivative of the ice sheet with respect to time. As
!$    ! the only property of the ice which explicitely changes with time is
!$    ! the ice thickness, that will be the only portion of the returned type
!$    ! which actually corresponds to the derivative.
!$    !
!$    class(ice_sheet), intent(in)   :: self
!$    real(r8), intent(in), optional :: t
!$      !! Time at which to evaluate the derivative
!$    class(integrand), allocatable  :: sheet_dt
!$      !! The time rate of change of the ice sheet. Has dynamic type
!$      !! [[ice_sheet]].
!$  end function sheet_dt
!$
!$  function sheet_local_error(lhs, rhs) result(error)
!$    !* Author: Christopher MacMackin
!$    !  Date: April 2016
!$    !
!$    ! Calculates a real scalar to represent the absolute difference between
!$    ! two ice_sheet objects. `rhs` must be a a [[ice_sheet]] object, or a
!$    ! runtime error will occur.
!$    !
!$    class(ice_sheet), intent(in) :: lhs
!$      !! Self
!$    class(integrand), intent(in) :: rhs
!$      !! The ice sheet object which is being compared against.
!$    real(r8) :: error
!$      !! The scalar representation of the absolute difference between these
!$      !! two ice shelves.
!$  end function sheet_local_error
!$
!$  function sheet_m_sheet(lhs, rhs) result(product)
!$    !* Author: Christopher MacMackin
!$    !  Date: April 2016
!$    !
!$    ! Multiplies one ice sheet object by another. That is to say, it 
!$    ! performs element-wise multiplication of the state vectors 
!$    ! representing the two arguments. `rhs` must be an [[ice_sheet]]
!$    ! object, or a runtime error will occur.
!$    !
!$    class(ice_sheet), intent(in)  :: lhs
!$      !! Self
!$    class(integrand), intent(in)  :: rhs
!$      !! The ice sheet object being multiplied by.
!$    class(integrand), allocatable :: product
!$      !! The product of the two arguments. Has dynamic type [[ice_sheet]].
!$  end function sheet_m_sheet
!$
!$  function sheet_m_real(lhs, rhs) result(product)
!$    !* Author: Christopher MacMackin
!$    !  Date: April 2016
!$    !
!$    ! Multiplies one ice sheet object by a scalar. That is to say, it 
!$    ! performs element-wise multiplication of the state vector 
!$    ! representing the ice sheet.
!$    !
!$    class(ice_sheet), intent(in)  :: lhs
!$      !! Self
!$    real(r8), intent(in)          :: rhs
!$      !! The scalar being multiplied by.
!$    class(integrand), allocatable :: product
!$      !! The product of the two arguments. Has dynamic type [[ice_sheet]].
!$  end function sheet_m_real
!$
!$  function real_m_sheet(lhs, rhs) result(product)
!$    !* Author: Christopher MacMackin
!$    !  Date: April 2016
!$    !
!$    ! Multiplies one ice sheet object by a scalar. That is to say, it 
!$    ! performs element-wise multiplication of the state vector 
!$    ! representing the ice sheet.
!$    !
!$    real(r8), intent(in)          :: lhs
!$      !! The scalar being multiplied by.
!$    class(ice_sheet), intent(in)  :: rhs
!$      !! Self
!$    class(integrand), allocatable :: product
!$      !! The product of the two arguments. Has dynamic type [[ice_sheet]].
!$  end function real_m_sheet
!$
!$  function sheet_add(lhs, rhs) result(sum)
!$    !* Author: Christopher MacMackin
!$    !  Date: April 2016
!$    !
!$    ! Adds one ice sheet object to another. That is to say, it performs
!$    ! element-wise addition of the state vectors representing the two
!$    ! arguments. `rhs` must be an [[ice_sheet]] object, or a runtime
!$    ! error will occur.
!$    !
!$    class(ice_sheet), intent(in)  :: lhs
!$      !! Self
!$    class(integrand), intent(in)  :: rhs
!$      !! The ice sheet object being added.
!$    class(integrand), allocatable :: sum
!$      !! The sum of the two arguments. Has dynamic type [[ice_sheet]].
!$  end function sheet_add
!$
!$  function sheet_sub(lhs, rhs) result(difference)
!$    !* Author: Christopher MacMackin
!$    !  Date: April 2016
!$    !
!$    ! Subtracts one ice sheet object from another. That is to say, it 
!$    ! performs element-wise addition of the state vectors representing 
!$    ! the two arguments. `rhs` must be a a [[ice_sheet]] object, or a
!$    ! runtime error will occur.
!$    !
!$    class(ice_sheet), intent(in)  :: lhs
!$      !! Self
!$    class(integrand), intent(in)  :: rhs
!$      !! The ice sheet object being subtracted.
!$    class(integrand), allocatable :: difference
!$      !! The difference of the two arguments. Has dynamic type [[ice_sheet]].
!$  end function sheet_sub
!$
!$  subroutine sheet_assign(lhs, rhs)
!$    !* Author: Christopher MacMackin
!$    !  Date: April 2016
!$    !
!$    ! Assigns the `rhs` ice sheet to this, `lhs`, one. All components
!$    ! will be the same following the assignment.
!$    !
!$    class(ice_sheet), intent(inout) :: lhs
!$      !! Self
!$    class(integrand), intent(in)    :: rhs
!$      !! The object to be assigned. Must have dynamic type [[ice_sheet]],
!$      !! or a runtime error will occur.
!$  end subroutine sheet_assign

  pure function sheet_thickness(this) result(thickness)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Returns the thickness of the ice sheet across its domain.
    !
    class(ice_sheet), intent(in) :: this
    class(scalar_field), pointer :: thickness !! The ice thickness.
  end function sheet_thickness

  function sheet_velocity(this) result(velocity)
    !* Author: Christopher MacMackin
    !  Date: July 2016
    !
    ! Returns the velocity of the ice sheet across its domain.
    !
    class(ice_sheet), intent(in)     :: this
    class(vector_field), allocatable :: velocity !! The ice velocity.
    allocate(velocity, source=this%velocity)
  end function sheet_velocity

  pure function sheet_density(this) result(density)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Returns the density of the ice in the sheet, which is assumed to be
    ! uniform across its domain.
    !
    class(ice_sheet), intent(in) :: this
    real(r8)                     :: density !! The ice density.
  end function sheet_density

  pure function sheet_temperature(this) result(temperature)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Returns the density of the ice in the sheet, which is assumed to be
    ! uniform across its domain.
    !
    class(ice_sheet), intent(in) :: this
    real(r8)                     :: temperature !! The ice density.
  end function sheet_temperature

  function sheet_residual(this, previous_states, melt_rate, basal_drag_parameter, &
                          water_density) result(residual)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Returns the residual when the current state of the glacier is run
    ! through the system of equations describing it. The residual takes the
    ! form of a 1D array, with each element respresenting the residual for
    ! one of the equations in the system.
    !
    class(ice_sheet), intent(in)             :: this
    class(glacier), dimension(:), intent(in) :: previous_states
      !! The states of the glacier in the previous time steps. The
      !! first element of the array should be the most recent. The
      !! default implementation will only make use of the most recent
      !! state, but the fact that this is an array allows overriding
      !! methods to use older states for higher-order integration
      !! methods.
    class(scalar_field), intent(in)          :: melt_rate
      !! Thickness of the ice above the glacier.
    class(scalar_field), intent(in)          :: basal_drag_parameter
      !! A paramter, e.g. coefficient of friction, needed to calculate the
      !! drag on basal surface of the glacier.
    real(r8), intent(in)                     :: water_density
      !! The density of the water below the glacier.
    real(r8), dimension(:), allocatable      :: residual
      !! The residual of the system of equations describing the glacier.
  end function sheet_residual

  subroutine sheet_update(this, state_vector)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Updates the state of the ice sheet from its state vector. The state
    ! vector is a real array containing the value of each of the ice sheet's
    ! properties at each of the locations on the grid used in descretization.
    !
    class(ice_sheet), intent(inout)    :: this
    real(r8), dimension(:), intent(in) :: state_vector
      !! A real array containing the data describing the state of the
      !! glacier.
  end subroutine sheet_update


  function sheet_precondition(this, previous_states, melt_rate, &
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
    class(ice_sheet), intent(inout)     :: this
    class(glacier), dimension(:), intent(in) :: previous_states
      !! The states of the glacier in the previous time steps. The
      !! first element of the array should be the most recent. The
      !! default implementation will only make use of the most
      !! recent state, but the fact that this is an array allows
      !! overriding methods to use older states for higher-order
      !! integration methods.
    class(scalar_field), intent(in)          :: melt_rate
      !! Thickness of the ice above the glacier
    class(scalar_field), intent(in)     :: basal_drag_parameter
      !! A paramter, e.g. coefficient of friction, needed to calculate
      !! the drag on basal surface of the glacier.
    real(r8), intent(in)                :: water_density
      !! The density of the water below the glacier
    real(r8), dimension(:), intent(in)  :: delta_state
      !! The change to the state vector which is being preconditioned.
    real(r8), dimension(:), allocatable :: preconditioned
      !! The result of applying the preconditioner to `delta_state`.
  end function sheet_precondition


  subroutine sheet_set_time(this, time)
    !* Author: Christopher MacMackin
    !  Date: November 2016
    !
    ! Sets the time information held by the ice sheet object. This is
    ! the time at which the ice sheet is in its current state.
    !
    class(ice_sheet), intent(inout) :: this
    real(r8), intent(in)            :: time
      !! The time at which the glacier is in the present state.
    this%time = time
  end subroutine sheet_set_time

  pure function sheet_data_size(this)
    !* Author: Christopher MacMackin
    !  Date: August 2016
    !
    ! Returns the number of elements in the ice sheet's state vector.
    ! This is the size of the vector returned by [[ice_sheet:residual]]
    ! and [[ice_sheet:state_vector]] and taken as an argument by 
    ! [[ice_sheet:update]].
    !
    class(ice_sheet), intent(in) :: this
    integer :: sheet_data_size
      !! The number of elements in the ice sheet's state vector.
  end function sheet_data_size

  pure function sheet_state_vector(this) result(state_vector) 
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Returns the state vector for the current state of the ice sheet. 
    ! This takes the form of a 1D array.
    !
    class(ice_sheet), intent(in)        :: this
    real(r8), dimension(:), allocatable :: state_vector
      !! The state vector describing the ice sheet.
  end function sheet_state_vector

  subroutine sheet_read_data(this,file_id,group_name,error)
    !* Author: Chris MacMackin
    !  Date: April 2017
    !
    ! Reads the state of the ice shelf object from the specified group
    ! in an HDF5 file. This sets the thickness, the velocity, and
    ! parameter values.
    !
    class(ice_sheet), intent(inout) :: this
    integer(hid_t), intent(in)      :: file_id
      !! The identifier for the HDF5 file/group in which this data is
      !! meant to be written.
    character(len=*), intent(in)    :: group_name
      !! The name to give the group in the HDF5 file storing the
      !! ice shelf's data.
    integer, intent(out)            :: error
      !! Flag indicating whether routine ran without error. If no
      !! error occurs then has value 0.
  end subroutine sheet_read_data

  subroutine sheet_write_data(this,file_id,group_name,error)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Writes the state of the ice sheet object to an HDF file in the
    ! specified group. This will consist of a thickness and a velocity
    ! dataset.
    !
    class(ice_sheet), intent(in) :: this
    integer(hid_t), intent(in)   :: file_id
      !! The identifier for the HDF5 file/group in which this data is
      !! meant to be written.
    character(len=*), intent(in) :: group_name
      !! The name to give the group in the HDF5 file storing the
      !! ice sheet's data.
    integer, intent(out)         :: error
      !! Flag indicating whether routine ran without error. If no
      !! error occurs then has value 0.
  end subroutine sheet_write_data

  function sheet_time_step(this) result(dt)
    !* Author: Chris MacMackin
    !  Date: December 2016
    !
    ! Calculates the time step for integrating the ice sheet, using
    ! the CFL condition.
    !
    class(ice_sheet), intent(in) :: this
    real(r8) :: dt
      !! The time-step to use
  end function sheet_time_step

  subroutine sheet_assign(this, rhs)
    !* Author: Chris MacMackin
    !  Date: February 2017
    !
    ! Copies the data from one ice sheet into another. This is only
    ! needed due to a bug in gfortran which means that the intrinsic
    ! assignment for glacier types is not using the appropriate
    ! defined assignment for the field components.
    !
    ! It does not assign the Jacobian object as it would take up quite
    ! a bit of extra space and it is unlikely that it would ever be
    ! needed without first having to be recalculated.
    !
    class(ice_sheet), intent(out) :: this
    class(glacier), intent(in)    :: rhs
      !! The ice sheet to be assigned to this one.
    select type(rhs)
    class is(ice_sheet)
      this%thickness = rhs%thickness
      this%velocity = rhs%velocity
      this%lambda = rhs%lambda
      this%chi = rhs%chi
      allocate(this%viscosity_law, source=rhs%viscosity_law)
      this%time = rhs%time
    class default
      call logger%fatal('ice_sheet%assign','Type other than `ice_sheet` '// &
                        'requested to be assigned.')
      error stop
    end select
#ifdef DEBUG
    call logger%debug('ice_sheet%assign','Copied ice sheet data.')
#endif
  end subroutine sheet_assign

end module ice_sheet_mod

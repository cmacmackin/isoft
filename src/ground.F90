!
!  ground.f90
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

module ground_mod
  !* Author: Christopher MacMackin
  !  Date: April 2016
  !  License: GPLv3
  !
  ! Provides a concrete implementation of the [[basal_surface]] data type,
  ! representing solid ground.
  !
  use iso_fortran_env, only: r8 => real64
  use basal_surface_mod, only: basal_surface
  use factual_mod, only: scalar_field
  use hdf5
  implicit none
  private

  type, extends(basal_surface), public :: ground
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! A concrete implementation of the [[basal_surface]] abstract data type,
    ! representing the ground beneath an ice sheet. At the moment this
    ! doesn't actually do anything.
    !
  contains
    procedure :: basal_melt => ground_melt
    procedure :: basal_drag_parameter => ground_drag_parameter
    procedure :: water_density => ground_water_density
    procedure :: update => ground_update
    procedure :: data_size => ground_data_size
    procedure :: state_vector => ground_state_vector
    procedure :: read_data => ground_read_data
    procedure :: write_data => ground_write_data
    procedure :: solve => ground_solve
  end type ground

  interface ground
    module procedure constructor
  end interface ground

contains

  function constructor() result(this)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    ! 
    ! Instantiates a [[ground]] object.
    !
    type(ground) :: this
  end function constructor


  function ground_melt(this) result(melt)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Computes and returns the melt rate at the bottom of the ice
    ! sheet due to interaction with the ground.
    !
    class(ground), intent(in)    :: this
    class(scalar_field), pointer :: melt
      !! The melt rate at the base of the ice sheet.
  end function ground_melt


  function ground_drag_parameter(this) result(drag)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Computes and returns a quantity which may be necessary to determine
    ! the frictional drag the ground exerts on the bottom of the ice
    ! sheet. An example would be the coefficient of friction. The 
    ! description of this method is left deliberately vague so that as not
    ! to constrain how the drag is parameterized.
    !
    class(ground), intent(in)    :: this
    class(scalar_field), pointer :: drag
      !! The value of a paramter describing the drag of the ground on the
      !! ice sheet.
  end function ground_drag_parameter


  function ground_water_density(this) result(density)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Computes and returns the density of the water beneath the ice sheet.
    ! This water would be subglacial discharge and would tend to lubricate
    ! the motion of the ice sheet. The density probably won't be important
    ! in the case of an ice sheet, but is included so that the ground data
    ! type can have the same interface as the [[plume]] data type.
    !
    class(ground), intent(in) :: this
    real(r8)                  :: density
      !! The density of any water at the base of the ice sheet.
  end function ground_water_density


  subroutine ground_update(this, state_vector, ice_thickness)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Updates the state of the ground from its state vector. The state
    ! vector is a real array containing the value of each of the ground's
    ! properties at each of the locations on the grid used in discretization.
    !
    class(ground), intent(inout)       :: this
    real(r8), dimension(:), intent(in) :: state_vector
      !! A real array containing the data describing the state of the
      !! ground.
    class(scalar_field), optional, intent(in) :: ice_thickness
      !! The ice thickness which, if present, will be used to update
      !! the calculation of the melt rate and/or drag parameter.
  end subroutine ground_update


  pure function ground_data_size(this)
    !* Author: Christopher MacMackin
    !  Date: August 2016
    !
    ! Returns the number of elements in the ground's state vector.
    ! This is the size of the vector returned by
    ! [[ground(type):state_vector]] and taken as an argument by
    ! [[ground(type):update]].
    !
    class(ground), intent(in) :: this
    integer                   :: ground_data_size
      !! The number of elements in the ground's state vector.
  end function ground_data_size


  pure function ground_state_vector(this) result(state_vector) 
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Returns the state vector for the current state of the ground. 
    ! This takes the form of a 1D array.
    !
    class(ground), intent(in)           :: this
    real(r8), dimension(:), allocatable :: state_vector
      !! The state vector describing the ground.
  end function ground_state_vector


  subroutine ground_read_data(this,file_id,group_name,error)
    !* Author: Chris MacMackin
    !  Date: April 2017
    !
    ! Reads the state of the ground object from the specified group in
    ! an HDF file.
    !
    class(ground), intent(inout) :: this
    integer(hid_t), intent(in)   :: file_id
      !! The identifier for the HDF5 file/group from which this data is
      !! meant to be read.
    character(len=*), intent(in) :: group_name
      !! The name of the group in the HDF5 file storing the
      !! ground's data.
    integer, intent(out)         :: error
      !! Flag indicating whether routine ran without error. If no
      !! error occurs then has value 0.
    integer(hid_t) :: group_id
    error = 0
  end subroutine ground_read_data


  subroutine ground_write_data(this,file_id,group_name,error)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Writes the state of the ground object to an HDF file in the
    ! specified group.
    !
    class(ground), intent(in)    :: this
    integer(hid_t), intent(in)   :: file_id
      !! The identifier for the HDF5 file/group in which this data is
      !! meant to be written.
    character(len=*), intent(in) :: group_name
      !! The name to give the group in the HDF5 file storing the
      !! ground's data.
    integer, intent(out)         :: error
      !! Flag indicating whether routine ran without error. If no
      !! error occurs then has value 0.
    integer(hid_t) :: group_id
    call h5gcreate_f(file_id, group_name, group_id, error)
    if (error /= 0) then
      write(*,*) 'WARNING: Error code',error,' returned when creating HDF '// &
                 'group', group_name
      write(*,*) '         Data IO not performed for ice shelf'
      return
    end if

    call h5gclose_f(group_id, error)
    if (error /= 0) then
      write(*,*) 'WARNING: Error code',error,' returned when closing HDF '// &
                 'group', group_name
      write(*,*) '         Possible bad IO'
    end if
  end subroutine ground_write_data


  subroutine ground_solve(this, ice_thickness, ice_density, &
                         ice_temperature, time, success)
    !* Author: Chris MacMackin
    !  Date: April 2017
    !
    ! Solves the state of the ground for the specified ice properties,
    ! at the specified time.
    !
    class(ground), intent(inout)    :: this
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
  end subroutine ground_solve

end module ground_mod

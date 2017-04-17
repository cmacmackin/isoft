!
!  cryosphere.f90
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

module cryosphere_mod
  !* Author: Christopher MacMackin
  !  Date: April 2016
  !  License: GPLv3
  !
  ! Provides a data structure representing a system of ice sheets and/or ice
  ! shelves, as well as the ground and/or ocean which they interact with. This
  ! is the fundamental data type of the ISOFT software suite.
  !
  use iso_fortran_env, only: r8 => real64
  !use foodie, only: integrand
  use basal_surface_mod, only: basal_surface
  use glacier_mod, only: glacier
  use meta_mod
  use hdf5
  use h5lt
  use logger_mod, only: logger => master_logger
  use penf, only: str
  implicit none
  private

  ! Names for objects and attributes in HDF output
  character(len=12), parameter, public :: hdf_glacier = 'glacier'
  character(len=18), parameter, public :: hdf_basal = 'basal_surface'
  character(len=13), parameter, public :: hdf_version = 'isoft_version'
  character(len=23), parameter, public :: hdf_comp_time = 'binary_compilation_time'
  character(len=16), parameter, public :: hdf_write_time = 'data_output_time'

  type, public :: cryosphere
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! A data structure representing glaciers, either as ice shelves or
    ! (eventually) ice sheets. It will allow coupled systems of
    ! glaciers as well as different basal couplings with the ocean or
    ! ground. This type is a subclass of the
    ! [FOODIE](https://github.com/Fortran-FOSS-Programmers/FOODIE)
    ! [integrand](http://fortran-foss-programmers.github.io/FOODIE/type/integrand.html),
    ! allowing it to take advantage of that set of integration
    ! libraries for evolution in time.
    !
    private
    class(glacier), allocatable               :: ice
      !! A model for the ice shelf or ice sheet
    class(basal_surface), allocatable         :: sub_ice
      !! A model for the ground or ocean underneath the ice
    real(r8)                                  :: time 
      !! The time in the simulation
  contains
    procedure :: integrate
    procedure :: write_data
    procedure :: time_step
  end type cryosphere

  interface cryosphere
    module procedure constructor
  end interface cryosphere

contains

  function constructor(ice, sub_ice) result(this)
    !* Author: Christopher MacMackin
    !  Date: November 2016
    !
    ! Create a new [[cryosphere]] object from the provided
    ! components. This will model the evolution of a glacier/ice
    ! shelf/ice sheet and its surroundings.
    !
    class(glacier), allocatable, intent(inout)       :: ice
      !! An object modelling the ice sheet or shelf component of this
      !! system. Will be deallocated on return.
    class(basal_surface), allocatable, intent(inout) :: sub_ice
      !! An object modelling the component of this system beneath the
      !! ice. Will be deallocated on return.
    type(cryosphere) :: this
    call move_alloc(ice, this%ice)
    call move_alloc(sub_ice, this%sub_ice)
    this%time = 0.0_r8
#ifdef DEBUG
    call logger%debug('cryosphere','Instantiated new cryosphere object.')
#endif
  end function constructor

  function time_step(this)
    !* Author: Chris MacMackin
    !  Date: December 2016
    !
    ! Calculates an appropriate time step with which to integrate the
    ! cryosphere so as not to cause numerical instability.
    !
    class(cryosphere), intent(in) :: this
    real(r8) :: time_step
    time_step = this%ice%time_step()
  end function time_step

  subroutine integrate(this,time)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Integrates the cryosphere forward until the specified `time` is
    ! reached.
    !
    class(cryosphere), intent(inout) :: this
    real(r8), intent(in)             :: time
      !! The time to which to integrate the cryosphere
    class(glacier), dimension(:), allocatable :: old_glaciers
    logical, save :: first_call = .true.
    logical :: success

    ! Normally the plume should be solved at the end of the previous
    ! iteration, but in the first iteration obviously there hasn't
    ! been a chance for this to happen yet.
    if (first_call) then
      ! As I am integrating only semi-implicitly and solving the plume
      ! for the current (rather than future) state, I think I should
      ! pass the current time. I only *think* that this is correct,
      ! however.
      call this%sub_ice%solve(this%ice%ice_thickness(), this%ice%ice_density(), &
                              this%ice%ice_temperature(), this%time, success)
      first_call = .false.
    end if

    allocate(old_glaciers(1), mold=this%ice)
    old_glaciers(1) = this%ice
    call this%ice%integrate(old_glaciers, this%sub_ice%basal_melt(), &
                            this%sub_ice%basal_drag_parameter(),  &
                            this%sub_ice%water_density(), time, success)

    ! Solve the plume so that it is ready for use in the next step of
    ! the time integration.
    call this%sub_ice%solve(this%ice%ice_thickness(), this%ice%ice_density(), &
                            this%ice%ice_temperature(), time, success)
    this%time = time
  end subroutine integrate

    
  subroutine write_data(this,outfile)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Writes the data describing the cryosphere to the disc as an HDF5
    ! file. `h5open_f` must have been once called prior to using this
    ! method. After the method has been used, `h5close_f` must be
    ! called once before the end of the program.
    !
    class(cryosphere), intent(in) :: this
    character(len=*), intent(in) :: outfile
      !! The file to which to write the data describing the state of the 
      !! cryosphere
    integer(hid_t) :: file_id, error_code
    call h5fopen_f(outfile, H5F_ACC_TRUNC_F, file_id, error_code)
    if (error_code /= 0) then
      call logger%error('cryosphere%write_data','Error code '//       &
                        str(error_code)//' returned when creating '// &
                        'HDF5 file '//outfile)
      return
    end if

    ! Write any whole-system data...
    call h5ltset_attribute_string_f(file_id,'/',hdf_version,version(),error_code)
    call h5ltset_attribute_string_f(file_id,'/',hdf_comp_time,compile_time(), &
                                    error_code)
    call h5ltset_attribute_string_f(file_id,'/',hdf_write_time,current_time(), &
                                    error_code)
    if (error_code /= 0) then
      call logger%warning('cryosphere%write_data','Error code '//      &
                          str(error_code)//' returned when writing '// &
                          'attributes to HDF5 file '//outfile)
    end if
    
    ! Call for subobjects
    call this%ice%write_data(file_id, hdf_glacier, error_code)
    call this%sub_ice%write_data(file_id, hdf_basal, error_code)

    call h5fclose_f(file_id, error_code)
    if (error_code /= 0) then
      call logger%warning('cryosphere%write_data','Error code '//      &
                          str(error_code)//' returned when closing '// &
                          'HDF5 file '//outfile)
    end if

#ifdef DEBUG
    call logger%debug('cryosphere%write_data','Wrote cryosphere data to '// &
                      'HDF file '//outfile)
#endif
  end subroutine write_data
  
end module cryosphere_mod

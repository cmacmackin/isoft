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
  use nitsol_mod, only: iplvl
  implicit none
  private

  ! Names for objects and attributes in HDF output
  character(len=12), parameter, public :: hdf_glacier = 'glacier'
  character(len=18), parameter, public :: hdf_basal = 'basal_surface'
  character(len=13), parameter, public :: hdf_version = 'isoft_version'
  character(len=23), parameter, public :: hdf_comp_time = 'binary_compilation_time'
  character(len=16), parameter, public :: hdf_write_time = 'data_output_time'
  character(len=15), parameter, public :: hdf_simulation_time = 'simulation_time'

  character(len=25), parameter, public :: hdf_crash_file = &
                                            'isoft_termination_dump.h5'

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
    class(glacier), allocatable       :: ice
      !! A model for the ice shelf or ice sheet
    class(basal_surface), allocatable :: sub_ice
      !! A model for the ground or ocean underneath the ice
    real(r8)                          :: time 
      !! The time in the simulation
    logical                           :: first_integration
      !! Indicates whether the cryosphere has been integrated before
      !! or not.
    real(r8)                          :: dt_factor = 1.0_r8
      !! A factor by which to reduce the time step
    real(r8)                          :: min_dt_factor = 1e-3_r8
      !! The smallest time step reduction to allow
    logical                           :: performing_time_step
      !! True if in the process of trying to get a time-step to
      !! successfully integrate.
  contains
    procedure :: initialise
    procedure :: integrate
    procedure :: write_data
    procedure :: time_step
    procedure :: reduce_time_step
    procedure :: increase_time_step
  end type cryosphere

contains

  subroutine initialise(this, ice, sub_ice)
    !* Author: Christopher MacMackin
    !  Date: November 2016
    !
    ! Initialise a cryosphere object from the provided
    ! components. This object will model the evolution of a
    ! glacier/ice shelf/ice sheet and its surroundings.
    !
    class(cryosphere), intent(out) :: this
    class(glacier), allocatable, intent(inout)       :: ice
      !! An object modelling the ice sheet or shelf component of this
      !! system. Will be deallocated on return.
    class(basal_surface), allocatable, intent(inout) :: sub_ice
      !! An object modelling the component of this system beneath the
      !! ice. Will be deallocated on return.
    call move_alloc(ice, this%ice)
    call move_alloc(sub_ice, this%sub_ice)
    this%time = 0.0_r8
    this%first_integration = .true.
#ifdef DEBUG
    call logger%debug('cryosphere','Instantiated new cryosphere object.')
#endif
  end subroutine initialise


  function time_step(this)
    !* Author: Chris MacMackin
    !  Date: December 2016
    !
    ! Calculates an appropriate time step with which to integrate the
    ! cryosphere so as not to cause numerical instability.
    !
    class(cryosphere), intent(inout) :: this
    real(r8) :: time_step
    time_step = this%dt_factor*this%ice%time_step()
!  contains
!    function factor()
!      real(r8) :: factor
!      factor = this%dt_factor
!      this%dt_factor = min(this%dt_factor+0.1_r8,1._r8)
!    end function factor
  end function time_step


  subroutine reduce_time_step(this)
    !* Author: Christopher MacMackin
    !  Date: April 2017
    !
    ! Reuces the time step by a factor of 2, unless doing so would
    ! take it below the minimum value.
    !
    class(cryosphere), intent(inout) :: this
    real(r8) :: new_factor
    new_factor = 0.7_r8 * this%dt_factor
    if (new_factor < this%min_dt_factor) then
      this%dt_factor = this%min_dt_factor
      call logger%warning('cryosphere%reduce_time_step','Attempting to '// &
                          'reduce time step factor below minimum value '// &
                          'of '//trim(str(this%min_dt_factor)))
    else
      this%dt_factor = new_factor
    end if
#ifdef DEBUG
    call logger%debug('cryosphere%reduce_time_step','Reducing time '// &
                      'step by factor of '//trim(str(this%dt_factor)))
#endif
  end subroutine reduce_time_step


  subroutine increase_time_step(this)
    !* Author: Christopher MacMackin
    !  Date: April 2017
    !
    ! Increases the time step by a factor of 2, unless doing so would
    ! take it above the maximum.
    !
    class(cryosphere), intent(inout) :: this
    this%dt_factor = min(1.2_r8 * this%dt_factor, 1._r8)
#ifdef DEBUG
    call logger%debug('cryosphere%increase_time_step','Reducing time '// &
                      'step by factor of '//trim(str(this%dt_factor)))
#endif
  end subroutine increase_time_step


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
    logical :: success
    real(r8) :: t, old_t

    if (time < this%time) then
      call logger%warning('cryosphere%integrate','Request made to '// &
                          'integrate cryosphere to earlier time '//  &
                          'than present state. No action taken.')
      return
    end if

    ! Normally the plume should be solved at the end of the previous
    ! iteration, but in the first iteration obviously there hasn't
    ! been a chance for this to happen yet.
    if (this%first_integration) then
      ! As I am integrating only semi-implicitly and solving the plume
      ! for the current (rather than future) state, I think I should
      ! pass the current time. I only *think* that this is correct,
      ! however.
      call this%sub_ice%solve(this%ice%ice_thickness(), this%ice%ice_density(), &
                              this%ice%ice_temperature(), this%time, success)
      this%first_integration = .false.
    end if
    allocate(old_glaciers(1), mold=this%ice)

    old_t = this%time
    t = min(this%time + this%time_step(), time)
    do while (t <= time)
      old_glaciers(1) = this%ice
      call this%ice%integrate(old_glaciers, this%sub_ice%basal_melt(), &
                              this%sub_ice%basal_drag_parameter(),     &
                              this%sub_ice%water_density(), t, success)
      if (.not. success) then
        this%ice = old_glaciers(1)
        if (this%dt_factor > this%min_dt_factor) then
          call logger%warning('cryosphere%integrate','Failure in nonlinear '// &
                              'solver. Reducing time step and trying again.')
          call this%reduce_time_step()
          t = min(old_t + this%time_step(), time) ! Not quite right...
          cycle
        else
          call logger%fatal('cryosphere%integrate','Failed to integrate '//  &
                            'glacier to time '//trim(str(t))//'! Writing '// &
                            'cryosphere state to file "'//hdf_crash_file//'".')
          call this%write_data(hdf_crash_file)
          iplvl = 4
          call this%ice%integrate(old_glaciers, this%sub_ice%basal_melt(), &
                                  this%sub_ice%basal_drag_parameter(),     &
                                  this%sub_ice%water_density(), t, success)
          error stop
        end if
      end if

      ! Solve the plume so that it is ready for use in the next step of
      ! the time integration.
      call this%sub_ice%solve(this%ice%ice_thickness(), this%ice%ice_density(), &
                              this%ice%ice_temperature(), t, success)

      if (success) then
        call logger%trivia('cryosphere%integrate','Successfully integrated '// &
                           'cryosphere to time '//trim(str(t)))
      else
        call logger%fatal('cryosphere%integrate','Failed to solve plume '//   &
                          'at time '//trim(str(t))//'! Writing cryosphere '// &
                          'state to file "'//hdf_crash_file//'".')
        call this%write_data(hdf_crash_file)
        error stop
      end if

      call this%increase_time_step()
      if (t >= time) exit
      old_t = t
      t = min(t + this%time_step(), time)
    end do

    this%time = time
    call logger%info('cryosphere%integrate','Successfully integrated '// &
                     'cryosphere to time '//trim(str(t)))
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
    call h5fcreate_f(outfile, H5F_ACC_TRUNC_F, file_id, error_code)
    if (error_code /= 0) then
      call logger%error('cryosphere%write_data','Error code '//    &
                        trim(str(error_code))//' returned when '// &
                        'creating HDF5 file '//outfile)
      return
    end if

    ! Write any whole-system data...
    call h5ltset_attribute_string_f(file_id,'/',hdf_version,version(),error_code)
    call h5ltset_attribute_string_f(file_id,'/',hdf_comp_time,compile_time(), &
                                    error_code)
    call h5ltset_attribute_string_f(file_id,'/',hdf_write_time,current_time(), &
                                    error_code)
    call h5ltset_attribute_double_f(file_id,'/',hdf_simulation_time,[this%time], &
                                    1_size_t,error_code)
    if (error_code /= 0) then
      call logger%warning('cryosphere%write_data','Error code '//    &
                          trim(str(error_code))//' returned when '// &
                          'writing attributes to HDF5 file '//outfile)
    end if
    
    ! Call for subobjects
    call this%ice%write_data(file_id, hdf_glacier, error_code)
    call this%sub_ice%write_data(file_id, hdf_basal, error_code)

    call h5fclose_f(file_id, error_code)
    if (error_code /= 0) then
      call logger%warning('cryosphere%write_data','Error code '//    &
                          trim(str(error_code))//' returned when '// &
                          'closing HDF5 file '//outfile)
    end if

#ifdef DEBUG
    call logger%debug('cryosphere%write_data','Wrote cryosphere data to '// &
                      'HDF file '//outfile)
#endif
  end subroutine write_data
  
end module cryosphere_mod

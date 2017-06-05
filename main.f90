!
!  main.f90
!  This file is part of ISOFT.
!  
!  Copyright 2017 Chris MacMackin <cmacmackin@gmail.com>
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

program isoft
  !* Author: Chris MacMackin
  !  Date: April 2017
  !  License: GPLv3
  !
  ! This is the driver program for ISOFT.
  !
  use iso_fortran_env, only: r8 => real64
  use logger_mod, only: debug, trivia, info, warning, error, fatal, &
                        logger => master_logger, logger_init
  use hdf5, only: h5open_f, h5close_f
  use penf, only: str
  use meta_mod
  use cryosphere_mod, only: cryosphere

  use glacier_mod, only: glacier
  use ice_shelf_mod, only: ice_shelf
  use viscosity_mod, only: abstract_viscosity
  use newtonian_viscosity_mod, only: newtonian_viscosity
  use glacier_boundary_mod, only: glacier_boundary
  use dallaston2015_glacier_boundary_mod, only: dallaston2015_glacier_boundary

  use basal_surface_mod, only: basal_surface
  use plume_mod, only: plume
  use entrainment_mod, only: abstract_entrainment
  use jenkins1991_entrainment_mod, only: jenkins1991_entrainment
  use melt_relationship_mod, only: abstract_melt_relationship
  use one_equation_melt_mod, only: one_equation_melt
  use ambient_mod, only: ambient_conditions
  use uniform_ambient_mod, only: uniform_ambient_conditions
  use equation_of_state_mod, only: equation_of_state
  use linear_eos_mod, only: linear_eos
  use plume_boundary_mod, only: plume_boundary
  use simple_plume_boundary_mod, only: simple_plume_boundary
  use dallaston2015_seasonal_mod, only: dallaston2015_seasonal_boundary
  use specfun_mod, only: ei

  implicit none

  ! Simulation parameters
  integer  :: grid_points
  real(r8), dimension(2,2) :: domain
  real(r8) :: end_time
  character(len=:), allocatable :: restart_file
  logical  :: end_on_steady, restart_from_file, restart_at_0
  
  ! Output parameters
  real(r8) :: output_interval
  character(len=:), allocatable :: hdf_base_name, logfile
  integer  :: stdout_lim, stderr_lim, logfile_lim, output_start

  ! Ice shelf parameters
  real(r8) :: chi, lambda, zeta, ice_temperature, courant

  ! Viscosity parameters
  real(r8) :: visc_coefficient

  ! Ice shelf boundary parameters
  real(r8) :: ice_thickness_lower
  real(r8), dimension(2) :: ice_velocity_lower

  ! Entrainment parameters
  real(r8) :: ent_coefficient

  ! Plume parameters
  real(r8) :: delta, nu, mu, r_val, nu_init
  logical  :: initialise_iteratively
  integer  :: initial_steps

  ! Melt parameters
  real(r8) :: alpha1, alpha2

  ! Ambient condition parameters
  real(r8) :: ambient_salinity, ambient_temperature

  ! Equation of state parameters
  real(r8) :: ref_density, ref_temperature, ref_salinity, &
              beta_s, beta_t

  ! Plume boundary parameters
  real(r8) :: discharge, offset
  real(r8) :: plume_thickness_lower, plume_temperature_lower, &
              plume_salinity_lower
  real(r8), dimension(2) :: plume_velocity_lower

  ! Variables for use in the program
  integer                      :: i, hdf_err
  real(r8)                     :: time
  real                         :: cpu_start, cpu_setup, cpu_end
  character(len=18), parameter :: hdf_file_format = '(a,"-",i0.4,".h5")'
  character(len=50)            :: hdf_filename

  type(cryosphere) :: system
    !! The ice-ocean system which is being simulated

  class(glacier), allocatable            :: ice
  type(ice_shelf), allocatable           :: shelf
  class(abstract_viscosity), allocatable :: viscosity
  class(glacier_boundary), allocatable   :: ice_bound

  class(basal_surface), allocatable              :: sub_ice
  type(plume), allocatable                       :: water
  class(abstract_entrainment), allocatable       :: entrainment
  class(abstract_melt_relationship), allocatable :: melt_relationship
  class(equation_of_state), allocatable          :: eos
  class(ambient_conditions), allocatable         :: ambient
  class(plume_boundary), allocatable             :: water_bound
  real(r8)                                       :: alpha, nu_tmp, nu_fact
  logical                                        :: success

  call cpu_time(cpu_start)

  ! Initialise variables to defaults
  grid_points = 50
  domain(1,:) = [0._r8, 6._r8]
  domain(2,:) = [-1._r8, 1._r8]
  end_time = 10.0_r8
  restart_file = "isoft-0020.h5"
  end_on_steady = .true.
  restart_from_file = .false.
  restart_at_0 = .false.

  output_interval = 0.5_r8
  hdf_base_name = "isoft"
  logfile = "isoft.log"
  stdout_lim = info
  stderr_lim = error
  logfile_lim = trivia
  output_start = 0

  chi = 4._r8
  lambda = 1e2_r8
  zeta = 1e-11_r8 ! Not actually used
  ice_temperature = -7._r8
  courant = 50.0_r8

  visc_coefficient = 1._r8

  ice_thickness_lower = 1._r8
  ice_velocity_lower = [1._r8, 0._r8]

  ent_coefficient = 1._r8

  delta = 3.6e-2_r8
  nu = 3.69e-2_r8
  mu = 0.799_r8
  r_val = 1.12_r8
  nu_init = 6.71e2_r8
  initialise_iteratively = .false.
  initial_steps = 217
  nu = 369

  alpha1 = 0.0182_r8
  alpha2 = 0.0238

  ambient_salinity = 1._r8
  ambient_temperature = 1._r8

  ref_density = 1.0_r8
  ref_temperature = 1._r8
  ref_salinity = 1._r8
  beta_s = 0.0271_r8
  beta_t = 0._r8

  discharge = 1.178e-3_r8
  offset = 0.001_r8
  plume_thickness_lower = offset
  plume_temperature_lower = 0._r8
  plume_salinity_lower = 0._r8
  plume_velocity_lower = [discharge/plume_thickness_lower, 0._r8]
  alpha = discharge**(1._r8/3._r8)/nu

  ! Set up IO
  call logger_init(logfile, stderr_lim, stdout_lim, logfile_lim)
  call h5open_f(hdf_err)
  if (hdf_err < 0) error stop

  ! Print welcome message
  call logger%info('isoft', 'Welcome to ISOFT: Ice Shelf/Ocean '// &
                   'Fluid- and Thermodynamics')
  call logger%info('isoft','ISOFT v'//version()//', compiled on '// &
                   compile_time())
  call logger%trivia('isoft', trim(compile_info()))
  
  ! Initialise cryosphere
  allocate(viscosity, source=newtonian_viscosity(visc_coefficient))
  allocate(ice_bound,                                                 &
           source=dallaston2015_glacier_boundary(ice_thickness_lower, &
           ice_velocity_lower(1), chi))
  allocate(shelf)
  call shelf%initialise(domain, [grid_points], h, u_ice, ice_temperature, &
                        viscosity, ice_bound, lambda, chi, zeta, courant)

  allocate(entrainment, source=jenkins1991_entrainment(ent_coefficient))
  allocate(melt_relationship, source=one_equation_melt(alpha1, alpha2))
  allocate(ambient,                                               &
           source=uniform_ambient_conditions(ambient_temperature, &
           ambient_salinity))
  allocate(eos,                                                          &
           source=linear_eos(ref_density, ref_temperature, ref_salinity, &
           beta_t, beta_s))
  allocate(water_bound, &
           source=simple_plume_boundary(plume_thickness_lower, &
           plume_velocity_lower, plume_temperature_lower, &
           plume_salinity_lower))
!  allocate(water_bound, &
!           source=dallaston2015_seasonal_boundary(plume_thickness_lower, &
!           13.823_r8, 0.9_r8, 1.0_r8, plume_temperature_lower))
  allocate(water)
  call water%initialise(domain, [grid_points], D, U_plume, T, S,      &
                        entrainment, melt_relationship, ambient, eos, &
                        water_bound, delta, nu, mu, r_val)

  if (initialise_iteratively .and. .not. restart_from_file) then
    ! Start by solving for plume at high diffusivity and reduce to
    ! desired value.
    nu_tmp = nu_init
    nu_fact = (nu/nu_init)**(1._r8/real(initial_steps, r8))
    do i = 1, initial_steps
      water%nu = nu_tmp
      print*, 'Nu = ', water%nu
      call water%solve(shelf%ice_thickness(), shelf%ice_density(), &
                       shelf%ice_temperature(), time, success)
      if (.not. success) then
        call logger%fatal('isoft', 'Failed to succesfully solve for plume '// &
                          'with diffusivity '//str(nu_tmp))
        error stop
      end if
      nu_tmp = nu_fact*nu_tmp
    end do
    water%nu = nu
    print*, 'Nu = ', water%nu
    call water%solve(shelf%ice_thickness(), shelf%ice_density(), &
                     shelf%ice_temperature(), time, success)
    if (.not. success) then
      call logger%fatal('isoft', 'Failed to succesfully solve for plume '// &
                        'with diffusivity '//str(nu))
      error stop
    end if    
  end if

  call move_alloc(shelf, ice)
  call move_alloc(water, sub_ice)
  call system%initialise(ice, sub_ice)

  if (restart_from_file) then
    call system%read_data(restart_file, .not. restart_at_0)
  end if
  
  call cpu_time(cpu_setup)
  call logger%info('isoft', 'Finished setting up simulation. Took '// &
                   trim(str(cpu_setup-cpu_start))//'s.')

  ! Run simulation
  i = output_start
  time = system%get_time()
  do while (end_time - time > courant*1.e-5_r8/grid_points)
    write(hdf_filename, hdf_file_format) hdf_base_name, i
    time = min(time + output_interval, end_time)
    call system%write_data(trim(hdf_filename))
    call system%integrate(time)
    i = i + 1
  end do

  ! Clean up
  call logger%info('isoft', 'Successfully simulated ice shelf '// &
                   'evolution to time '//trim(str(time)))
  write(hdf_filename, hdf_file_format) hdf_base_name, i
  call system%write_data(trim(hdf_filename))
  call h5close_f(hdf_err)
  if (hdf_err < 0) error stop
  call logger%info('isoft', 'Wrote final ice shelf state to file '// &
                   trim(hdf_filename))
  
  ! Print goodbye message
  call cpu_time(cpu_end)
  call logger%info('isoft','Finished running isoft. This took '// &
                   trim(str(cpu_end-cpu_start))//'s, of which '// &
                   trim(str(cpu_end-cpu_setup))//'s were '//      &
                   'needed to run the simulation.')

contains

  pure function h(x)
    !! Initial ice thickness.
    real(r8), dimension(:), intent(in) :: x
    real(r8)                           :: h
    h = ice_thickness_lower*(1._r8 - x(1)/(1.1_r8*domain(1,2)))
  end function h

  pure function u_ice(x)
    !! Initial ice velocity
    real(r8), dimension(:), intent(in)  :: x
    real(r8), dimension(:), allocatable :: u_ice
    real(r8) :: big_x
    big_x = 1._r8/lambda
    allocate(u_ice(1))
    u_ice(1) = ice_velocity_lower(1) + 0.25_r8*ice_thickness_lower*chi*(x(1) &
             - x(1)**2/(2.2*domain(1,2)))
  end function u_ice

  pure function D(x)
    !! Initial guess for plume thickness
    real(r8), dimension(:), intent(in) :: x
    real(r8)                           :: D
    D = offset + 1_r8*(h([0._r8]) - h(x))/r_val
  end function D

  pure function U_plume(x)
    !! Initial guess for plume velocity
    real(r8), dimension(:), intent(in)  :: x
    real(r8), dimension(:), allocatable :: U_plume
    allocate(U_plume(1))
    U_plume = plume_velocity_lower(1)!/(0.5*x(1) + 1._r8)
  end function U_plume

  pure function S(x)
    !! Initial guess of the plume salinity
    real(r8), dimension(:), intent(in) :: x
    real(r8)                           :: S
    S = plume_salinity_lower - 0.1*x(1)*(x(1) - 2*domain(1,2))!*exp(x(1)/domain(1,2))! - 0.0001_r8*x(1)*(x(1) - 2*domain(1,2))
  end function S

  pure function T(x)
    !! Initial guess of the plume temperature
    real(r8), dimension(:), intent(in) :: x
    real(r8)                           :: T
    T = plume_temperature_lower - 0.1*x(1)*(x(1) - 2*domain(1,2))!*exp(x(1)/domain(1,2))! - 0.0001_r8*x(1)*(x(1) - 2*domain(1,2))
  end function T

  
end program isoft

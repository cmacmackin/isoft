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
  use dallaston2015_melt_mod, only: dallaston2015_melt
  use ambient_mod, only: ambient_conditions
  use uniform_ambient_mod, only: uniform_ambient_conditions
  use equation_of_state_mod, only: equation_of_state
  use linear_eos_mod, only: linear_eos
  use plume_boundary_mod, only: plume_boundary
  use simple_plume_boundary_mod, only: simple_plume_boundary
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
  integer  :: stdout_lim, stderr_lim, logfile_lim

  ! Ice shelf parameters
  real(r8) :: chi, lambda, ice_temperature, courant

  ! Viscosity parameters
  real(r8) :: visc_coefficient

  ! Ice shelf boundary parameters
  real(r8) :: ice_thickness_lower
  real(r8), dimension(2) :: ice_velocity_lower

  ! Entrainment parameters
  real(r8) :: ent_coefficient

  ! Plume parameters
  real(r8) :: delta, nu, mu, r_val

  ! Melt parameters
  real(r8) :: beta, epsilon_g, epsilon_m

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
  real(r8)                                       :: alpha

  call cpu_time(cpu_start)

  ! Initialise variables to defaults
  grid_points = 150
  domain(1,:) = [0._r8, 2.5_r8]
  domain(2,:) = [-1._r8, 1._r8]
  end_time = 5._r8
  restart_file = "steady_state.h5"
  end_on_steady = .true.
  restart_from_file = .false.
  restart_at_0 = .true.

  output_interval = 0.5_r8
  hdf_base_name = "isoft"
  logfile = "isoft.log"
  stdout_lim = info
  stderr_lim = error
  logfile_lim = trivia

  chi = 4._r8
  lambda = 0.37_r8
  ice_temperature = -1._r8
  courant = 20._r8

  visc_coefficient = 1._r8

  ice_thickness_lower = 1._r8
  ice_velocity_lower = [1._r8, 0._r8]

  ent_coefficient = 1._r8

  delta = 3.6e-2_r8
  nu = 1e-2_r8
  mu = 1.27_r8
  r_val = 1.12_r8

  beta = 2.4e-2_r8
  epsilon_m = 6.9e-4_r8
  epsilon_g = 1.1e-3_r8

  ambient_salinity = 0._r8
  ambient_temperature = 0._r8

  ref_density = 1.0_r8
  ref_temperature = 0._r8
  ref_salinity = 0._r8
  beta_s = -1.0_r8
  beta_t = 0._r8

  discharge = 1._r8
  offset = 0.1_r8
  plume_thickness_lower = offset
  plume_temperature_lower = -1._r8
  plume_salinity_lower = discharge**(2._r8/3._r8)/plume_thickness_lower
  plume_velocity_lower = [discharge**(1._r8/3._r8), 0._r8]
  alpha = discharge**(1._r8/3._r8)/nu

  ! Set up IO
  call logger_init(logfile, stderr_lim, stdout_lim, logfile_lim)
  call h5open_f(hdf_err)
  if (hdf_err < 0) error stop

  ! Print welcome message
  call logger%info('isoft', 'Welcome to ISOFT: Ice Shelf/Ocean '// &
                   'Fluid- and Thermodynamics')
  call logger%info('isoft','ISOFT v'//version()//' compiled on '// &
                   compile_time())
  call logger%info('isoft', trim(compile_info()))
  
  ! Initialise cryosphere
  allocate(viscosity, source=newtonian_viscosity(visc_coefficient))
  allocate(ice_bound,                                                 &
           source=dallaston2015_glacier_boundary(ice_thickness_lower, &
           ice_velocity_lower(1), chi))
  allocate(shelf)
  call shelf%initialise(domain, [grid_points], h, u_ice, ice_temperature, &
                        viscosity, ice_bound, lambda, chi, courant)

  allocate(entrainment, source=jenkins1991_entrainment(ent_coefficient))
  allocate(melt_relationship,                                     &
           source=dallaston2015_melt(beta, epsilon_m, epsilon_g))
  allocate(ambient,                                               &
           source=uniform_ambient_conditions(ambient_temperature, &
           ambient_salinity))
  allocate(eos,                                                          &
           source=linear_eos(ref_density, ref_temperature, ref_salinity, &
           beta_t, beta_s))
  allocate(water_bound,                                        &
           source=simple_plume_boundary(plume_thickness_lower, &
           plume_velocity_lower, plume_temperature_lower,      &
           plume_salinity_lower))
  allocate(water)
  call water%initialise(domain, [grid_points], D, U_plume, T, S,      &
                        entrainment, melt_relationship, ambient, eos, &
                        water_bound, delta, nu, mu, r_val)

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
  i = 0
  time = system%get_time()
  do while (time < end_time)
    write(hdf_filename, hdf_file_format) hdf_base_name, i
    call system%write_data(trim(hdf_filename))
    time = min(time + output_interval, end_time)
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
    real(r8), dimension(1) :: vel
    real(r8)               :: big_x
    big_x = 1._r8/lambda
    vel = u_ice(x)
    h = (1._r8 - x(1)/big_x)/vel(1)
  end function h

  pure function u_ice(x)
    !! Initial ice velocity
    real(r8), dimension(:), intent(in)  :: x
    real(r8), dimension(:), allocatable :: u_ice
    real(r8) :: big_x
    big_x = 1._r8/lambda
    allocate(u_ice(1))
    u_ice(1) = sqrt(1._r8 + 0.25_r8*chi*big_x - &
                    0.25_r8*chi*big_x*(1._r8 - x(1)/big_x)**2)
  end function u_ice

  pure function D(x)
    !! Initial guess for plume thickness
    real(r8), dimension(:), intent(in) :: x
    real(r8)                           :: D
    D = offset + (h([0._r8]) - h(x))/r_val
  end function D

  pure function U_plume(x)
    !! Initial guess for plume velocity
    real(r8), dimension(:), intent(in)  :: x
    real(r8), dimension(:), allocatable :: U_plume
    allocate(U_plume(1))
    U_plume = discharge**(1._r8/3._r8)
  end function U_plume

  pure function s1(x)
    !! First component of analytic solution for plume salinity
    real(r8), intent(in) :: x
    real(r8)             :: s1
    s1 = exp(alpha*x)
  end function s1

  pure function s2(x)
    !! Second component of analytic solution for plume salinity
    real(r8), intent(in) :: x
    real(r8)             :: s2
    s2 = s1(x)*ei(-alpha*(x+offset))
  end function s2

  pure function S(x)
    !! Initial guess of the plume salinity
    real(r8), dimension(:), intent(in) :: x
    real(r8)                           :: S
    real(r8) :: a12, a21, a22, phi, theta
    a12 = ei(-alpha*offset)
    a21 = alpha*s1(domain(1,2))
    a22 = alpha*s2(domain(1,2)) + exp(-alpha*offset)/(domain(1,2) + offset)
    phi = -a21 * plume_salinity_lower/(a22 - a21*a12)
    theta = plume_salinity_lower - a12*phi
    S = theta*s1(x(1)) + phi*s2(x(1))
  end function S

  pure function T(x)
    !! Initial guess of the plume temperature
    real(r8), dimension(:), intent(in) :: x
    real(r8)                           :: T
    real(r8) :: zeta, a12, a21, a22, phi, theta
    zeta = r_val*(beta+1._r8)
    a12 = ei(-alpha*offset)
    a21 = alpha*s1(domain(1,2))
    a22 = alpha*s2(domain(1,2)) + exp(-alpha*offset)/(domain(1,2) + offset)
    phi = -a21 * (plume_temperature_lower - zeta)/(a22 - a21*a12)
    theta = plume_temperature_lower - zeta - a12*phi
    T = theta*s1(x(1)) + phi*s2(x(1)) + zeta
  end function T

  
end program isoft

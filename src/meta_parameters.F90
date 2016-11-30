!
!  meta_parameters.f90
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

module meta_parameters
  !* Author: Chris MacMackin
  !  Date: November 2016
  !  License: GPLv3
  !
  ! Provides parameters specifying the version of ISOFT, time of compilation, etc.
  !
  use iso_fortran_env, only: i8 => int64
  implicit none
  
  character(len=20), parameter :: compile_time = __DATE__//' '//__TIME__
  character(len=5),  parameter :: version = '0.1.0'

  character(len=3), dimension(12), parameter :: months = ['Jan', &
                                                          'Feb', &
                                                          'Mar', &
                                                          'Apr', &
                                                          'May', &
                                                          'Jun', &
                                                          'Jul', &
                                                          'Aug', &
                                                          'Sep', &
                                                          'Oct', &
                                                          'Nov', &
                                                          'Dec']
  character(len=42), parameter :: time_format = '(a3,1x,i2,1x,i4,1x,i2.2,":",'// &
                                                'i2.2,":",i2.2)'

contains

  function get_time()
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Returns the current time in the same format as the `compile_time` parameter.
    !
    character(len=20) :: get_time
    integer(i8), dimension(8) :: time_vals
    call date_and_time(values=time_vals)
    write(get_time,time_format) months(time_vals(2)), time_vals(3), time_vals(1), &
                                time_vals(5), time_vals(6), time_vals(7)
  end function get_time

end module meta_parameters

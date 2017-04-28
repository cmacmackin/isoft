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

submodule(meta_mod) meta_implementation_mod
  !* Author: Chris MacMackin
  !  Date: November 2016
  !  License: GPLv3
  !
  ! Implements functions specifying the version of ISOFT, time of
  ! compilation, etc. This is done in a submodule because it should be
  ! recompiled with every build (so that compilation time is accurate)
  ! but without a submodule that would require every module using it
  ! to be recompiled as well.
  !
  implicit none
  
  character(len=20), parameter :: compile_time_val = __DATE__//' '//__TIME__
  character(len=5),  parameter :: version_num = '0.1.0'

contains

  module function version()
    version = version_num
  end function version

  module function compile_time()
    compile_time = compile_time_val
  end function compile_time

  module function compile_info()
    character(len=40), parameter :: result_format = &
                                    '("Compiled with ",a," using options ",a)'
    write(compile_info,result_format) compiler_version(), compiler_options()
    compile_info = 'Compiled with "'//compiler_version()//   &
                   '" using options "'//compiler_options()// &
                   '"'
  end function compile_info

end submodule meta_implementation_mod

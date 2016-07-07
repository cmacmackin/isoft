!
!  nitsol.f90
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

module nitsol_mod
  !* Author: Christopher MacMackin
  !  Date: July 2016
  !  License: GPLv3
  !
  ! Provides an explicit interface to the 
  ! [NITSOL](http://users.wpi.edu/~walker/Papers/nitsol,SISC_19,1998,302-318.pdf) 
  ! package. At some point I may produce a proper object oriented interface
  ! for it.
  !
  use iso_fortran_env, only: r8 => real64

end module nitsol_mod

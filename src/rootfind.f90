!
!  rootfind.f90
!  This file is part of ISOFT.
!  
!  Copyright 2017 Chris MacMackin <cmacmackin@gmail.com>
!  
!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 3 of the License, or
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

MODULE rootfind
!! author: Chris MacMackin
!!
!! Provides subroutines implementing various root-finding algorithms, 
!! as well as a global bracket finder.

CONTAINS

    PURE SUBROUTINE bis_secant   ( error   , fnctn   , maxerr  , maxsteps,  &
                                   steps   , xcur    , xleft   , xright     )
      !! A root finder using the hybrid bisection-secant algorithm.
      !! Returns a computed value of the root within xleft and xright
      !! once error < maxerr or has run maximum number of iterations.
        IMPLICIT NONE

        interface
            pure function fnctn ( x        )
              !! The Fortran function for which the root will be
              !! found. Must take only one real argument, return a
              !! real argument.
                real(8), intent(in) ::  x
                real(8)             ::  fnctn
            end function fnctn
        end interface
        
        integer, intent(in)     ::  maxsteps
          !! An integer value for the maximum number of iterations
          !! applied to the algorithm before terminating.
        real(8), intent(in)     ::  maxerr
          !! A real value which specifies the the maximum allowable
          !! error in the computed root. The subroutine terminates
          !! once the error passes below this level.
        real(8), intent(inout)  ::  xleft
          !! A real value specifying lower bound within which to
          !! search for root.
        real(8), intent(inout)  ::  xright
          !! A real value specifying upper bound within which to
          !! search for root.
        integer, intent(out)    ::  steps
          !! Returns the number of iterations needed before error
          !! falls below maxerr (or returns maxsteps if that many
          !! iterations occur first).
        real(8), intent(out)    ::  error
          !! A real variable in which an estimate of the error in the
          !! computed root will be stored and returned.
        real(8), intent(out)    ::  xcur
          !! A real variable in which the computed value of the root
          !! will be returned.
        
        ! Other variables:
        REAL(8) ::  dfxcur,                                     &
                    fxcur,                                      &
                    fxleft,                                     &
                    fxprev1,                                    &
                    fxright,                                    &
                    xnext,                                      &
                    xprev1
    !------------------------------------------------------------------!

        ! Initialize variables
        fxleft = fnctn(xleft)
        fxright = fnctn(xright)
        dfxcur = 0.d0
        fxcur = 0.d0
        fxprev1 = 0.d0
        xnext = 0.d0
        xprev1 = 0.d0
        
        ! Check that there is a root bracketed by 'xleft' and 'xright'
        IF ( fxleft * fxright > 0.0d0 ) THEN
           error stop ('No root in specified interval.')
        END IF
        
        ! Initialize variables)
        xcur = xleft
        fxcur = fxleft
        xprev1 = xright
        fxprev1 = fxright
        
        ! Iterate until root converges or reach 'maxsteps'
        DO steps = 1, maxsteps
            ! If not 1st iteration, update variables for next iteration
            IF ( steps /= 1 ) THEN
                fxprev1 = fxcur
                fxcur = fnctn(xcur)
                IF ( fxleft * fxcur < 0.d0 ) THEN
                    xright = xcur
                    fxright = fxcur
                ELSE
                    xleft = xcur
                    fxleft = fxcur
                END IF
            END IF
            
            ! Compute a secant and use to find 'xnext'
            dfxcur = (fxcur - fxprev1) / (xcur - xprev1)
            xnext = xcur - fxcur / dfxcur
            ! If this 'xnext' value is outside of left or right bounds, 
            ! instead use a bisection to find 'xnext'
            IF ( ( xnext < xleft ) .OR. ( xnext > xright ) ) THEN
                xnext = 5.d-1 * (xleft + xright)
            END IF
            
            ! Estimate error in approximating root and update variables
            error = ABS(xnext - xcur)
            xprev1 = xcur
            xcur = xnext
            
            ! If error less than tolerance, return
            IF ( error < maxerr ) RETURN
        END DO
        
        ! If reached maximum iterations, print a warning message
        error stop ('No solution found in specified number of steps.')

    END SUBROUTINE bis_secant



    SUBROUTINE global_bis_sec   ( dx      , error   , fnctn   , maxerr  ,      &
                                  maxsteps, numroots, roots   , steps   ,      &
                                  verbose , xmax    , xmin     )
      !! A subroutine which finds the values of all roots of a
      !! function (or as many as will fit into the provided arrays)
      !! within a given range of values. This subroutine uses the
      !! hybrid bisection-secant root-finding algorithm.
        IMPLICIT NONE

        interface
            pure function fnctn ( x        )
              !! The Fortran function for which the root will be
              !! found. Must take only one real argument, return a
              !! real argument.
                real(8), intent(in) ::  x
                real(8)             ::  fnctn
            end function fnctn
        end interface
        
        ! Input and output variables:
        integer, intent(in)                 ::  maxsteps
          !! An integer value for the maximum number of iterations
          !! applied to the algorithm before terminating.
        logical, intent(in)                 ::  verbose
          !! A logical variable which specifies whether to print
          !! progress to stdout as brackets found and at each
          !! iteration as root found. Also says whether to print a
          !! warning if 'dx' set to 'dxmin' during 'globrack' routine
          !! and if maximum number of iterations reached while finding
          !! root.
        real(8), intent(in)                 ::  maxerr
          !! A real value which specifies the the maximum allowable
          !! error in the computed root. The subroutine terminates
          !! once the error passes below this level.
        real(8), intent(in)                 ::  xmax
          !! The upper limit of the range on which the subroutine will
          !! search for roots and brackets.
        real(8), intent(in)                 ::  xmin
          !! The lower limit of the range on which the subroutine will
          !! search for roots and brackets.
        real(8), intent(inout)              ::  dx
          !! The initial size of increment to use when examining
          !! function.  Minimum interval will be 0.01 of this.
        integer, intent(out)                ::  numroots
          !! An integer value which will return the number of roots
          !! for which brackets were found. A negative number
          !! indicates that an error occurred.
        integer, intent(out), dimension(:)  ::  steps
          !! An integer array in which the number of iterations needed
          !! before error falls below maxerr for each root is stored
          !! and returned.
        real(8), intent(out), dimension(:)  ::  error
          !! A real array in which an estimate of the error in each
          !! computed root will be stored and returned.
        real(8), intent(out), dimension(:)  ::  roots
          !! A real array in which the computed values of each root
          !! will be returned.
        
        ! Other variables:
        INTEGER                                 ::  counter = 0,       &
                                                    maxroots = 0
        REAL(8), DIMENSION(:,:), ALLOCATABLE    ::  brackets
    !------------------------------------------------------------------!
    
        ! Maximum number of roots is given by the size of the smallest
        ! of the arrays 'error', 'roots', and 'steps' (of course, if
        ! it is written properly then the calling program should have
        ! made these all the same size)
        maxroots = MIN(SIZE(error),SIZE(roots),SIZE(steps))
        ALLOCATE(brackets(2,maxroots))
        brackets = 0.d0
        
        ! Get the sets of bracket pairs for the roots of 'fnctn', using
        ! the 'globrack' subroutine 
        CALL globrack(brackets, dx, fnctn, numroots, verbose, xmax,    &
                      xmin)
        
        ! For each bracket pair, find the enclosed root using the 
        ! 'biscnt' subroutine
        DO counter = 1,numroots
            CALL biscnt(error(counter), fnctn, maxerr, maxsteps,       &
                        steps(counter), verbose, roots(counter),       &
                        brackets(1,counter), brackets(2,counter))
        END DO

        RETURN
    END SUBROUTINE global_bis_sec



    SUBROUTINE global_brackets  ( brackets, dx      , fnctn   , numroots,      &
                                  verbose , xmax    , xmin     )
      !! A global bracket finder. For a given function it finds values
      !! on each side of each of the function's roots within a given
      !! range.
        IMPLICIT NONE

        interface
            pure real(8) function fnctn ( x        )
              !!The Fortran function for which the brackets will be
              !! found. Must take only one real argument, return a
              !! real argument.
                real(8), intent(in) ::  x
            end function fnctn
        end interface
        
        ! Input and output variables:
        logical, intent(in)                     ::  verbose
          !! A logical variable which specifies whether to print to
          !! stdout any bracket values which are found and warning
          !! messages when 'dx' set to 'dxmin'.
        real(8), intent(in)                     ::  xmax
          !! The upper limit of the range on which the subroutine will
          !! search for roots and brackets.
        real(8), intent(in)                     ::  xmin
          !! The lower limit of the range on which the subroutine will
          !! search for roots and brackets.
        real(8), intent(inout)                  ::  dx
          !! The initial size of increment to use when examining
          !! function.  Minimum interval will be 0.01 of this.
        integer, intent(out)                    ::  numroots
          !! An integer value which will return the number of roots
          !! for which brackets were found. A negative number
          !! indicates that an error occurred.
        real(8), intent(out), dimension(:,:)    ::  brackets
          !! A 2 by n real array in which the left and right brackets
          !! will be stored and returned. Will find up to n sets of
          !! brackets.
        
        ! Other variables:
        REAL(8) ::  dfdx = 0.d0,                                       &
                    dfdxh = 0.d0,                                      &
                    dldx = 0.d0,                                       &
                    dldxh = 0.d0,                                      &
                    dxmin = 0.d0,                                      &
                    dxh = 0.d0,                                        &
                    fleft = 0.d0,                                      &
                    fmid = 0.d0,                                       &
                    fright = 0.d0,                                     &
                    scaleval = 0.d0,                                   &
                    xleft = 0.d0,                                      &
                    xright = 0.d0,                                     &
                    xmid = 0.d0
    !------------------------------------------------------------------!
   
        ! Check if passed array is big enough to hold returned data
        IF ( SIZE(brackets, 1) < 2 ) THEN
            numroots = -1
            WRITE( 0,2000)
            RETURN
        END IF
        
        ! Initialize variables
        dxmin = 1.d-2 * dx
        numroots = 0
        xleft = xmin
        fleft = fnctn(xmin)
        
        ! Repeat this loop until have traversed range [xmin:xmax], or
        ! until brackets array is filled
        DO WHILE ( xleft <= xmax )
            ! Make sure that xright <= xmax
            dx = MIN(dx, xmax - xmin)
            
            ! Get values and derivatives for next step along the 
            ! function
            xright = xleft + dx
            fright = fnctn(xright)
            scaleval = MAX(1.0d0, ABS(xleft)) / MAX(ABS(fleft),        &
                       ABS(fright)) 
            dfdx = scaleval * (fright - fleft) / dx
            dldx = SIGN(1.0d0, dfdx) * SQRT(1.0d0 + dfdx**2.d0)
            
            ! Find derivative of secant length (called 'dldx' when found
            ! above) at midpoint between 'xleft' and 'xright' ('dldxh'). 
            ! If difference between 'dldx' and 'dldxh' is too great, 
            ! reduce 'dx' and repeat.
            DO
                dxh = 5d-1 * dx
                xmid = xleft + dxh
                fmid = fnctn(xmid)
                dfdxh = scaleval * (fmid - fleft) / dxh
                dldxh = SIGN(1.0d0, dfdxh) * SQRT(1.0d0 + dfdxh**2)
                IF ( ( ABS(dldx - dldxh) > 5.d-1*ABS(dldx) )           &
                     .AND. ( dxh > dxmin ) ) THEN
                    dx = dxh
                    xright = xmid
                    fright = fmid
                    dldx = dldxh
                ELSE
                    EXIT
                END IF
            END DO
            
            ! If difference between 'dldx' and 'dldxh' is too small,
            ! make step-size larger for next iteration
            IF ( ABS(dldx - dldxh) < 1.d-1 * ABS(dldx) ) THEN
                dx = 1.5d0 * dx
            ! If diference between 'dldx' and 'dldxh' is too small, but
            ! can't go smaller without falling below 'dxmin', set 'dx =
            ! dxmin' and print warning message
            ELSE IF ( ( ABS(dldx - dldxh) > 5.d-1 * ABS(dldx) )   &
                      .AND. ( dxh < dxmin ) ) THEN
                dx = dxmin
                xright = xleft + dx
                fright = fnctn(xright)
                IF ( verbose .EQV. .TRUE. ) WRITE( 6,2010) dxmin
            END IF
            ! Otherwise keep 'dx' the same (no action needs to be taken)
            
            ! If the signs of 'fleft' and 'fright' are opposite, and if
            ! signs of 'fright' and the derivative of 'fnctn(xleft)' are
            ! the same, then 'xleft' and 'xright' bracket a root
            IF ( ( fleft * fright < 0 ) .AND. ( fright * (fleft -      &
                 fnctn(xleft - dxmin) ) > 0 ) ) THEN
                numroots = numroots + 1 
                brackets(1:2,numroots) = (/ xleft, xright /)
                IF ( verbose .EQV. .TRUE. ) WRITE( 6,2020) xleft, xright
                ! If have filled 'brackets', then leave loop
                IF ( numroots >= SIZE(brackets,2) ) EXIT
            END IF
            
            ! Update variables for next iteration
            xleft = xright
            fleft = fright
        END DO
        
        !--------------------------------------------------------------!
        !                    Write format statements                   !
        !--------------------------------------------------------------!
        2000 FORMAT('GLOBRACK: ERROR: Passed array too small to hold ',&
                    'returned values. Must have at least 2 columns.')
        2010 FORMAT('GLOBRACK: WARNING: Value of dx set to minimum '   &
                   ,'value: ',1PG22.15)
        2020 FORMAT('GLOBRACK: Found root bracketed by ',1PG22.15,     &
                    ' and',/,                                          &
                    'GLOBRACK: ',1PG22.15)
        !--------------------------------------------------------------!

        RETURN
    END SUBROUTINE global_brackets
    !==================================================================!
    !                  E N D    S U B R O U T I N E :                  !
    !                  G L O B A L _ B R A C K E T S                   !
    !==================================================================!

END MODULE rootfind

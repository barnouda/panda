subroutine intersfacex(xf,y1,y2,z1,z2,x0,y0,z0,a,b,c,inters,xi,yi,zi)

    implicit none

    real(8) :: xf, y1, y2, z1, z2    ! rectangle definition
    real(8) :: x0, y0, z0, a, b, c   ! line definition
    logical :: inters                ! true if intersection
    real(8) :: xi, yi, zi            ! intersection coordintes
      
    real(8) :: t

    ! parametric equation of a straight line
    ! x = a*t + x0
    ! y = b*t + y0
    ! z = c*t + z0

    inters = .false.
    xi = xf

    if (a /= 0) then

        ! intersection of line with plane x=xf
        ! xf = a*t + x0
        t = (xi-x0) / a

        ! check if intersect plane within face
        yi = b*t + y0
        zi = c*t + z0

        ! last condition to ensure we are not looking behind the detector
        if (yi>=y1 .and. yi<=y2 .and. zi>=z1 .and. zi<=z2) then!TMP .and. zi<=z0) then
            inters = .true.
        end if

!    else
!        if (xi == x0) then
!            yi = yf
!            if (b /= 0) then
!                
!            end if
!            inters = .true.
!        end if
    end if

end subroutine intersfacex

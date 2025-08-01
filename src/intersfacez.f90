subroutine intersfacez(x1,x2,y1,y2,zf,x0,y0,z0,a,b,c,inters,xi,yi,zi)

    implicit none

    real(8) :: x1, x2, y1, y2, zf    ! rectangle definition
    real(8) :: x0, y0, z0, a, b, c   ! line definition
    logical :: inters                ! true if intersection
    real(8) :: xi, yi, zi            ! intersection coordintes
      
    real(8) :: t

    ! parametric equation of a straight line
    ! x = a*t + x0
    ! y = b*t + y0
    ! z = c*t + z0

    inters = .false.
    zi = zf

    if (c /= 0) then

        ! intersection of line with plane x=xf
        ! xf = a*t + x0
        t = (zi-z0) / c

        ! check if intersect plane within face
        xi = a*t + x0
        yi = b*t + y0
    
        if (xi>=x1 .and. xi<=x2 .and. yi>=y1 .and. yi<=y2) then!TMP .and. zi<=z0) then
            inters = .true.
        end if

    end if

end subroutine intersfacez

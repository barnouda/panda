subroutine intersfacey(x1,x2,yf,z1,z2,x0,y0,z0,a,b,c,inters,xi,yi,zi)

    implicit none

    real(8) :: x1, x2, yf, z1, z2    ! rectangle definition
    real(8) :: x0, y0, z0, a, b, c   ! line definition
    logical :: inters                ! true if intersection
    real(8) :: xi, yi, zi            ! intersection coordintes
      
    real(8) :: t

    ! parametric equation of a straight line
    ! x = a*t + x0
    ! y = b*t + y0
    ! z = c*t + z0

    inters = .false.
    yi = yf

    if (b /= 0) then

        ! intersection of line with plane x=xf
        ! xf = a*t + x0
        t = (yi-y0) / b

        ! check if intersect plane within face
        xi = a*t + x0
        zi = c*t + z0
    
        if (xi>=x1 .and. xi<=x2 .and. zi>=z1 .and. zi<=z2) then!TMP .and. zi<=z0) then
            inters = .true.
        end if

    end if

end subroutine intersfacey

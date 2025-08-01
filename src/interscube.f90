subroutine interscube(x1,x2,y1,y2,z1,z2,x0,y0,z0,a,b,c,inters,xi1,yi1,zi1,xi2,yi2,zi2)

    implicit none

    real(8) :: x1, x2, y1, y2, z1, z2    ! rectangle definition
    real(8) :: x0, y0, z0, a, b, c       ! line definition
    logical :: inters                    ! true if intersection with cube
    !real(8), dimension(2) :: xi, yi, zi
    ! intersections with 6 faces
    real(8), dimension(6) :: xi, yi, zi  ! intersection with the 6 cube faces
    real(8), dimension(6) :: si          ! distance from intersection to telescope
    logical, dimension(6) :: iface
    integer :: i1, i2
    real(8) :: xi1, xi2, yi1, yi2, zi1, zi2
      
    ! parametric equation of a straight line
    ! x = a*t + x0
    ! y = b*t + y0
    ! z = c*t + z0

    ! 6 faces of a cube
    !
    !    7------8
    !  5------6 |
    !  | |    | |
    !  | 3 - -|-4
    !  1------2
    !
    ! face 1243 zmin
    ! face 5687 zmax
    ! face 1265 ymin
    ! face 3487 ymax
    ! face 1375 xmin
    ! face 2486 xmax

    inters = .false.
    xi = x0
    yi = y0
    zi = z0

!write(*,*) 'xi = ', xi

    call intersfacex(x1,y1,y2,z1,z2,x0,y0,z0,a,b,c,iface(1),xi(1),yi(1),zi(1)) ! xmin
    call intersfacex(x2,y1,y2,z1,z2,x0,y0,z0,a,b,c,iface(2),xi(2),yi(2),zi(2)) ! xmax
    call intersfacey(x1,x2,y1,z1,z2,x0,y0,z0,a,b,c,iface(3),xi(3),yi(3),zi(3)) ! ymin
    call intersfacey(x1,x2,y2,z1,z2,x0,y0,z0,a,b,c,iface(4),xi(4),yi(4),zi(4)) ! ymax
    call intersfacez(x1,x2,y1,y2,z1,x0,y0,z0,a,b,c,iface(5),xi(5),yi(5),zi(5)) ! zmin
    call intersfacez(x1,x2,y1,y2,z2,x0,y0,z0,a,b,c,iface(6),xi(6),yi(6),zi(6)) ! zmax

    si = sqrt((xi-x0)**2+(yi-y0)**2+(zi-z0)**2)
    i1 = minloc(si,1,iface)
    i2 = maxloc(si,1,iface)

!write(*,*) 'FACE INTER =', iface, i1, i2

!write(*,*) 'INTERSCUBE', si, iface, i1, i2

    ! line crosses cube only if 2 distinct intersections
    !if (i1/=0 .and. i2/=0 .and. i1/=i2) then
    ! if starts inside cube
    if (x0>=x1 .and. x0<=x2 .and. y0>=y1 .and. y0<=y2 .and. z0>=z1 .and. z0<=z2) then
        if (i1/=0 .and. i2/=0 .and. i1/=i2) then
            inters = .true.
            xi1 = xi(i1)
            yi1 = yi(i1)
            zi1 = zi(i1)
            xi2 = xi(i2)
            yi2 = yi(i2)
            zi2 = zi(i2)
        elseif (i1/=0 .and. i2/=0 .and. i1==i2) then
            inters = .true.
            xi1 = x0
            yi1 = y0
            zi1 = z0
            xi2 = xi(i1)
            yi2 = yi(i1)
            zi2 = zi(i1)
        !elseif (i1/=0 .and. i2==0) then
        !    inters = .true.
        !    xi1 = x0
        !    yi1 = y0
        !    zi1 = z0
        !    xi2 = xi(i1)
        !    yi2 = yi(i1)
        !    zi2 = zi(i1)
        !elseif (i1==0 .and. i2/=0) then
        !    inters = .true.
        !    xi1 = x0
        !    yi1 = y0
        !    zi1 = z0
        !    xi2 = xi(i2)
        !    yi2 = yi(i2)
        !    zi2 = zi(i2)
        end if
    ! if line starts outside cube, crosses cube only if two distinct intersections
    elseif (x0<x1 .or. x0>x2 .or. y0<y1 .or. y0>y2 .or. z0<z1 .or. z0>z2) then
        if (i1/=0 .and. i2/=0 .and. i1/=i2) then
            inters = .true.
            xi1 = xi(i1)
            yi1 = yi(i1)
            zi1 = zi(i1)
            xi2 = xi(i2)
            yi2 = yi(i2)
            zi2 = zi(i2)
        end if
    end if

end subroutine interscube

subroutine bilinint(nx,ny,x,y,v,ni,xi,yi,vi,nanvalue)

    use omp_lib

    implicit none

    integer :: nx,ny
    real(8), dimension(nx) :: x
    real(8), dimension(ny) :: y
    real(8), dimension(nx*ny) :: v

    integer :: ni
    real(8), dimension(ni) :: xi,yi
    real(8), dimension(ni) :: vi

    real(8) :: nanvalue

    integer :: i

    real(8) :: a1,a2,b1,b2
    real(8) :: x1,x2,y1,y2
    integer, dimension(1) :: vtest
    integer :: test
    integer :: i1,i2,j1,j2

    integer :: l11,l12,l21,l22
    integer :: v11,v12,v21,v22

    real(8) :: xmin,xmax,ymin,ymax

    xmin = minval(x)
    xmax = maxval(x)
    ymin = minval(y)
    ymax = maxval(y)

    ! parallel default(none), &
        ! shared(ni,xi,xmin,xmax,yi,ymin,ymax,vi,x,y,nanvalue,nx,ny,v), &
        ! private(i,vtest,test,i1,i2,j1,j2,x1,x2,y1,y2,l11,l12,l21,l22,v11,v12,v21,v22, &
        ! a1,a2,b1,b2)
    ! do schedule(static,ni/omp_get_num_threads())
    do i = 1,ni
!write(*,*) 'OMP LOOP', i, '/', ni

        ! NaN if outside grid coordinates
        if (xi(i) < xmin .or. xi(i) > xmax .or. yi(i) < ymin .or. yi(i) > ymax) then
            vi(i) = nanvalue
            cycle
        end if

        ! Find indices i1 and i2 of neighbours along x axis
        vtest = minloc(abs(x-xi(i)))
        test = vtest(1)
        if (xi(i) >= x(test) .and. test < nx) then
            i1 = test
            i2 = test+1
        else if (xi(i) < x(test) .or. test == nx) then
            i1 = test-1
            i2 = test
        end if

        ! Find indices y1 and y2 of neighbours along y axis
        vtest = minloc(abs(y-yi(i)))
        test = vtest(1)
        if (yi(i) >= y(test) .and. test < ny) then
            j1 = test
            j2 = test+1
        else if (yi(i) < y(test) .or. test == ny) then
            j1 = test-1
            j2 = test
        end if

        ! Coordinates of the 4 nodes
        x1 = x(i1)
        x2 = x(i2)
        y1 = y(j1)
        y2 = y(j2)

        ! Linear indices of the 4 nodes
        !l11 = i1+(j1-1)*nx
        !l12 = i1+(j2-1)*nx
        !l21 = i2+(j1-1)*nx
        !l22 = i2+(j2-1)*nx
        l11 = j1+(i1-1)*ny
        l12 = j1+(i2-1)*ny
        l21 = j2+(i1-1)*ny
        l22 = j2+(i2-1)*ny

        ! Extract values at the 4 nodes
        v11 = v(l11)
        v12 = v(l12)
        v21 = v(l21)
        v22 = v(l22)

        ! Compute coefficients a1,a2,b1,b2
        a1 = (x2-xi(i))/(x2-x1)
        a2 = (xi(i)-x1)/(x2-x1)
        b1 = (y2-yi(i))/(y2-y1)
        b2 = (yi(i)-y1)/(y2-y1)

        ! Compute interpolated value
        vi(i) = a1*b1*v11 + a1*b2*v12 + a2*b1*v21 + a2*b2*v22

    end do
    ! end do
    ! end parallel

end subroutine bilinint


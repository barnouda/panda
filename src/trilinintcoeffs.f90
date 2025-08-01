subroutine trilinintcoeffs(nx,ny,nz,x,y,z,v,ni,xi,yi,zi,vi,ci,nanvalue)

    use omp_lib

    implicit none

    integer :: nx,ny,nz
    real(8), dimension(nx) :: x
    real(8), dimension(ny) :: y
    real(8), dimension(nz) :: z
    real(8), dimension(nx*ny*nz) :: v

    integer :: ni
    real(8), dimension(ni) :: xi,yi,zi
    real(8), dimension(ni) :: vi
    real(8), dimension(ni,nx*ny*nz) :: ci

    real(8) :: nanvalue

    integer :: i

    real(8) :: a1,a2,b1,b2,c1,c2
    real(8) :: x1,x2,y1,y2,z1,z2
    integer, dimension(1) :: vtest
    integer :: test
    integer :: i1,i2,j1,j2,k1,k2

    integer :: l111,l121,l211,l221,l112,l122,l212,l222
    integer :: v111,v121,v211,v221,v112,v122,v212,v222

    real(8) :: xmin,xmax,ymin,ymax,zmin,zmax

    xmin = minval(x)
    xmax = maxval(x)
    ymin = minval(y)
    ymax = maxval(y)
    zmin = minval(z)
    zmax = maxval(z)

    ci = 0
    vi = 0

    ! parallel default(none), &
        ! shared(nx,ny,nz,x,y,z,v,ni,xi,yi,zi,vi,xmin,xmax,ymin,ymax,zmin,zmax,nanvalue,ci), &
        ! private(i,a1,a2,b1,b2,c1,c2,x1,x2,y1,y2,z1,z2,vtest,test,i1,i2,j1,j2,k1,k2, &
        ! l111,l121,l211,l221,l112,l122,l212,l222,v111,v121,v211,v221,v112,v122,v212,v222)
    ! do schedule(static,ni/omp_get_num_threads())
    do i = 1,ni

        ! NaN if outside grid coordinates
        if (xi(i) < xmin .or. xi(i) > xmax .or. yi(i) < ymin .or. yi(i) > ymax .or. zi(i) < zmin .or. zi(i) > zmax) then
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

        ! Find indices z1 and z2 of neighbours along z axis
        vtest = minloc(abs(z-zi(i)))
        test = vtest(1)
        if (zi(i) >= z(test) .and. test < nz) then
            k1 = test
            k2 = test+1
        else if (zi(i) < z(test) .or. test == nz) then
            k1 = test-1
            k2 = test
        end if

        ! Coordinates of the 8 nodes
        x1 = x(i1)
        x2 = x(i2)
        y1 = y(j1)
        y2 = y(j2)
        z1 = z(k1)
        z2 = z(k2)

        ! Linear indices of the 8 nodes
        l111 = i1+(j1-1)*nx+(k1-1)*nx*ny
        l121 = i1+(j2-1)*nx+(k1-1)*nx*ny
        l211 = i2+(j1-1)*nx+(k1-1)*nx*ny
        l221 = i2+(j2-1)*nx+(k1-1)*nx*ny
        l112 = i1+(j1-1)*nx+(k2-1)*nx*ny
        l122 = i1+(j2-1)*nx+(k2-1)*nx*ny
        l212 = i2+(j1-1)*nx+(k2-1)*nx*ny
        l222 = i2+(j2-1)*nx+(k2-1)*nx*ny

        ! Extract values at the 8 nodes
        v111 = v(l111)
        v121 = v(l121)
        v211 = v(l211)
        v221 = v(l221)
        v112 = v(l112)
        v122 = v(l122)
        v212 = v(l212)
        v222 = v(l222)

        ! Compute coefficients a1,a2,b1,b2,c1,c2
        a1 = (x2-xi(i))/(x2-x1)
        a2 = (xi(i)-x1)/(x2-x1)
        b1 = (y2-yi(i))/(y2-y1)
        b2 = (yi(i)-y1)/(y2-y1)
        c1 = (z2-zi(i))/(z2-z1)
        c2 = (zi(i)-z1)/(z2-z1)

        ! Fill up grid of coefficients
        ci(i,l111) = a1*b1*c1
        ci(i,l121) = a1*b2*c1
        ci(i,l211) = a2*b1*c1
        ci(i,l221) = a2*b2*c1
        ci(i,l112) = a1*b1*c2
        ci(i,l122) = a1*b2*c2
        ci(i,l212) = a2*b1*c2
        ci(i,l222) = a2*b2*c2

        ! Compute interpolated value
        vi(i) = a1*b1*c1*v111 + a1*b2*c1*v121 + a2*b1*c1*v211 + a2*b2*c1*v221 &
              + a1*b1*c2*v112 + a1*b2*c2*v122 + a2*b1*c2*v212 + a2*b2*c2*v222

    end do
    ! end do
    ! end parallel

end subroutine trilinintcoeffs

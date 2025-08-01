subroutine trilinintcoeffs(x1,x2,y1,y2,z1,z2,xa,ya,za,xb,yb,zb,x0,y0,z0,azim,elev,ci)
! Analytical computation of the contribution of the 8 surrounding nodes to
! a segment within these 8 nodes and belonging to a line defined by
! (x0,y0,z0,azim,elev)
! The contribution of the 8 nodes is stored in ci

    implicit none

    real(8) :: x1, x2, y1, y2, z1, z2           ! cube limits
    real(8) :: xa, ya, za, xb, yb, zb           ! start and end segment points
    real(8) :: x0, y0, z0                       ! origin of the los
    real(8) :: azim, elev                       ! azimuth and elevation of los
    real(8), dimension(8) :: ci                 ! node contribution

!    real(8) :: nanvalue

    real(16) :: sa, sb
    real(16) :: ds1, ds2, ds3, ds4
    real(16) :: alpha, beta
    !real(8), parameter :: pi=3.14159265
    real(16), parameter :: pi = 4*atan(1.0_16)

    real(16) :: cosa, sina, cosb, sinb
    real(16) :: a, b, c
    real(16) :: int_1_ds
    real(16) :: int_x_ds, int_y_ds, int_z_ds
    real(16) :: int_xy_ds, int_xz_ds, int_yz_ds
    real(16) :: int_xyz_ds

    integer :: i, j, k, kk, kkk

    real(8) :: ddd
    integer :: di1, di2, dj1, dj2, dk1, dk2
    real(16) :: ddx, ddy, ddz, ddi, ddj, ddk

    real(16) :: c000, c100, c010, c001, c110, c101, c011, c111
    real(16), dimension(8) :: ccc
    real(16) :: ctmp

    !real(8) :: ztmp

    ! TEST TRANSFORM Z POSITIF VERS LE HAUT...
    !ztmp = z1
    !z1 = -z2
    !z2 = -ztmp
    !za = -za
    !zb = -zb
    !z0 = -z0

!write(*,*) 'pi = ', pi

    ci = 0

    sa = sqrt((xa-x0)**2+(ya-y0)**2+(za-z0)**2)
    sb = sqrt((xb-x0)**2+(yb-y0)**2+(zb-z0)**2)

    !if (sb < sa) then
    !    write(*,*) 'DEBUGTLIC (sa,sb) =', sa, sb
    !end if

    ds1 = sb-sa
    ds2 = (sb**2-sa**2)/2
    ds3 = (sb**3-sa**3)/3
    ds4 = (sb**4-sa**4)/4

    ddd = (x2-x1)*(y2-y1)*(z2-z1)

    !if (x2-x1 < 0 .or. y2-y1 < 0 .or. z2-z1 < 0) then
    !    write(*,*) 'DEBUGTLIC (x1,x2,y1,y2,z1,z2) =', x1, x2, y1, y2, z1, z2
    !end if

    alpha = elev
    beta = azim

    sina = sin(alpha)
    cosa = cos(alpha)
    sinb = sin(beta)
    cosb = cos(beta)

    ! coefficient directeur de la droite
    a = sinb*cosa
    b = cosb*cosa
    c = -sina ! car z positif vers le bas
    !c = sina

!write(*,*) 'DEBUG TLIC (x1 x2 y1 y2 z1 z2)', x1, x2, y1, y2, z1, z2
!write(*,*) 'DEBUG TLIC (xa ya za xb yb zb)', xa, ya, za, xb, yb, zb
!write(*,*) 'DEBUG TLIC (x0 y0 z0 azim elev)', x0, y0, z0, azim, elev
!write(*,*) 'DEBUG TLIC (sa sb ds1 ds2 ds3 ds4)', sa, sb, ds1, ds2, ds3, ds4
!write(*,*) 'DEBUG TLIC (ddd a b c)', ddd, a, b, c

    !write(*,*) '---'
    !write(*,*) a, b, c
    !write(*,*) ds1, ds2, ds3, ds4
    !write(*,*) ddd

    ! Integrals
    int_1_ds = ds1
    int_x_ds = x0*ds1 + a*ds2
    int_y_ds = y0*ds1 + b*ds2
    int_z_ds = z0*ds1 + c*ds2
    int_xy_ds = x0*y0*ds1 + (x0*b+y0*a)*ds2 + a*b*ds3
    int_xz_ds = x0*z0*ds1 + (x0*c+z0*a)*ds2 + a*c*ds3
    int_yz_ds = y0*z0*ds1 + (y0*c+z0*b)*ds2 + b*c*ds3
    int_xyz_ds = x0*y0*z0*ds1 + (x0*y0*c+x0*z0*b+y0*z0*a)*ds2 &
       + (x0*b*c+y0*a*c+z0*a*b)*ds3 &
       + a*b*c*ds4

    if (b==0 .and. c==0) then
        write(*,*) 'TRILIN', int_1_ds, int_x_ds, int_y_ds, int_z_ds
        write(*,*) 'TRILIN', int_xy_ds, int_xz_ds, int_yz_ds, int_xyz_ds
    end if

    ! Contribution of the 8 nodes

    do i = 1,2
        if (i == 1) then
            di1 = 1
            di2 = 0
        elseif (i == 2) then
            di1 = 0
            di2 = 1
        end if
        ddi = di2-di1
        ddx = di1*x2-di2*x1
        do j = 1,2
            if (j == 1) then
                dj1 = 1
                dj2 = 0
            elseif (j == 2) then
                dj1 = 0
                dj2 = 1
            end if
            ddj = dj2-dj1
            ddy = dj1*y2-dj2*y1
            do k = 1,2
                if (k == 1) then
                    dk1 = 1
                    dk2 = 0
                elseif (k == 2) then
                    dk1 = 0
                    dk2 = 1
                end if
                ddk = dk2-dk1
                ddz = dk1*z2-dk2*z1
                c000 = ddx*ddy*ddz
                c100 = ddy*ddz*ddi
                c010 = ddx*ddz*ddj
                c001 = ddx*ddy*ddk
                c110 = ddz*ddi*ddj
                c101 = ddy*ddi*ddk
                c011 = ddx*ddj*ddk
                c111 = ddi*ddj*ddk
                ! pour eviter les erreurs numeriques lors de l'addition dans ci
                ! on range par ordre croissant (NB : 4 positifs et 4 negatifs)
                ! et on additionne d'abord 2 par 2 (1 positif + 1 negatif)
                ccc(1) = c000*int_1_ds
                ccc(2) = c100*int_x_ds
                ccc(3) = c010*int_y_ds
                ccc(4) = c001*int_z_ds
                ccc(5) = c110*int_xy_ds
                ccc(6) = c101*int_xz_ds
                ccc(7) = c011*int_yz_ds
                ccc(8) = c111*int_xyz_ds
                do kk = 2,8
                    do kkk = kk,2,-1
                    if (ccc(kkk-1) > ccc(kkk)) then
                        ctmp = ccc(kkk)
                        ccc(kkk) = ccc(kkk-1)
                        ccc(kkk-1) = ctmp
                    end if
                    end do
                end do

    if (b==0 .and. c==0) then
        write(*,*) 'TRILIN', ddi, ddj, ddk
        write(*,*) 'TRILIN', ddx, ddy, ddz, ddd
        write(*,*) 'TRILIN', c000, c100, c010, c001
        write(*,*) 'TRILIN', c110, c101, c011, c111
    end if

                           !l111 = 1
                           !l211 = 2
                           !l121 = 1+2 = 3
                           !l221 = 2+2 = 4
                           !l112 = 1+4 = 5
                           !l212 = 2+4 = 6
                           !l122 = 1+2+4 = 7
                           !l222 = 2+2+4 = 8
                ! maxval: si contribution negative (erreurs numeriques), on clip a zero
                ci(i+(j-1)*2+(k-1)*4) = maxval((/0.0_8,dble((ccc(1)+ccc(8))+(ccc(2)+ccc(7))+(ccc(3)+ccc(6))+(ccc(4)+ccc(5)))/))
                !ci(i+(j-1)*2+(k-1)*4) = c000*int_1_ds &
                !    + c100*int_x_ds + c010*int_y_ds + c001*int_z_ds &
                !    + c110*int_xy_ds + c101*int_xz_ds + c011*int_yz_ds &
                !    + c111*int_xyz_ds

                !if (ci(i+(j-1)*2+(k-1)*4) < 0) then
                !    write(*,*) 'DEBUGTLIC C000', c000, int_1_ds, c000*int_1_ds
                !    write(*,*) 'DEBUGTLIC C100', c100, int_x_ds, c100*int_x_ds
                !    write(*,*) 'DEBUGTLIC C010', c010, int_y_ds, c010*int_y_ds
                !    write(*,*) 'DEBUGTLIC C001', c001, int_z_ds, c001*int_z_ds
                !    write(*,*) 'DEBUGTLIC C110', c110, int_xy_ds, c110*int_xy_ds
                !    write(*,*) 'DEBUGTLIC C101', c101, int_xz_ds, c101*int_xz_ds
                !    write(*,*) 'DEBUGTLIC C011', c011, int_yz_ds, c011*int_yz_ds
                !    write(*,*) 'DEBUGTLIC C111', c111, int_xyz_ds, c111*int_xyz_ds
                !    !write(*,*) 'DEBUGTLIC', x1,x2,y1,y2,z1,z2,xa,ya,za,xb,yb,zb,x0,y0,z0,azim,elev
                !end if
                !write(*,*) 'DEBUGTLIC', (ci(i+(j-1)*2+(k-1)*4)<0), &
                !    (ccc(1) > 0), &
                !    (ccc(2) > 0), &
                !    (ccc(3) > 0), &
                !    (ccc(4) > 0), &
                !    (ccc(5) > 0), &
                !    (ccc(6) > 0), &
                !    (ccc(7) > 0), &
                !    (ccc(8) > 0), &
                    !c000*int_1_ds, &
                    !c000*int_1_ds + c100*int_x_ds, &
                    !c000*int_1_ds + c100*int_x_ds + c010*int_y_ds, &
                    !c000*int_1_ds + c100*int_x_ds + c010*int_y_ds + c001*int_z_ds, &
                    !c000*int_1_ds + c100*int_x_ds + c010*int_y_ds + c001*int_z_ds + c110*int_xy_ds, &
                    !c000*int_1_ds + c100*int_x_ds + c010*int_y_ds + c001*int_z_ds + c110*int_xy_ds + &
                    !    c101*int_xz_ds, &
                    !c000*int_1_ds + c100*int_x_ds + c010*int_y_ds + c001*int_z_ds + c110*int_xy_ds + &
                    !    c101*int_xz_ds + c011*int_yz_ds, &
                !    ci(i+(j-1)*2+(k-1)*4), &
                !    c000*int_1_ds + c100*int_x_ds + c010*int_y_ds + c001*int_z_ds + c110*int_xy_ds + &
                 !       c101*int_xz_ds + c011*int_yz_ds + c111*int_xyz_ds
            end do
        end do
    end do

    ci = ci/ddd

end subroutine trilinintcoeffs

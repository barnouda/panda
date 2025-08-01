program panda

    use omp_lib
!    use vtr

    implicit none

!    logical, parameter :: debug = .true.
!    logical :: test

    character(len=100) :: inputfile
    character(len=100) :: outputfile, topofile, gridparamfile, densitygridfile, televiewfile
    character(len=1) :: telestr
    character(len=100) :: muogdatacorrfile
    logical :: muogdatacorr

    ! density node grid
    real(8) :: xmin,xmax,ymin,ymax,zmin,zmax
    real(8) :: dx,dy,dz, dxy
    integer :: nx,ny,nz
    real(8), dimension(:), allocatable :: x,y,z

    real(8), dimension(:), allocatable :: rho

    ! telescope view rays
    integer :: nmuog
    real(8), dimension(:), allocatable :: xtele, ytele, ztele
    real(8), dimension(:), allocatable :: dazim, delev
    real(8), dimension(:), allocatable :: azim, elev
!    real(8), dimension(:), allocatable :: dmax
    real(8) :: az, el
    integer :: iaz, iel

    ! topography
    integer :: nxt, nyt
    real(8) :: xmint, xmaxt, ymint, ymaxt, dxyt
    real(8), dimension(:), allocatable :: xt, yt
    real(8), dimension(:,:), allocatable :: topo

    ! blocs
    integer :: nxb, nyb, nzb, nb
    real(8), dimension(:), allocatable :: x1b, x2b, y1b, y2b, z1b, z2b
    integer, dimension(:), allocatable :: l111, l211, l121, l112, l221, l212, l122, l222 ! corresponding indices
    integer :: ib ! to iterate over blocs

    ! lines of sight
    real(8), dimension(8) :: dist = 0 ! to fine maximal distance
    real(8) :: smax ! maximal distance for a line of sight
    real(8) :: azmoy, elmoy, az1, az2, el1, el2, azinf, azsup, elinf, elsup
    real(8) :: alostmp, elostmp
    real(8) :: daz, del
    real(8) :: dalos, delos
    integer :: nalos, nelos
    real(8) :: nlos ! nb of effective los per bin weighted by cos(elev)
    real(8), dimension(:), allocatable :: alos, elos
    integer :: ia, ie ! to iterate over lines of sight
    logical, dimension(:), allocatable :: btl
    integer :: nin, nout
    integer, dimension(:), allocatable :: iin, iout
    real(8), dimension(:), allocatable :: xin, yin, zin, xout, yout, zout
    integer :: i1, i2
    real(8) :: t1, t2
    real(8) :: alphalos, alphatopo
    real(8) :: tin, tout

    real(8) :: a, b, c ! coefficients directeurs des los
    logical :: imesh ! intersection of los with mesh or not
    integer, dimension(:), allocatable :: indb ! indexes of blocs to investigate
    real(8) :: xi1, yi1, zi1, xi2, yi2, zi2 ! coordinates of intersections with mesh
    logical :: ibloc
    real(8) :: xi1b, yi1b, zi1b, xi2b, yi2b, zi2b ! coordinates of intersections with blocs

    real(8), dimension(:), allocatable :: slos, xslos, yslos, zslos, ztopo
    integer :: nslos
    real(8) :: dslos

    real(8), dimension(:), allocatable :: clos  ! node contribution to one LoS
    real(8) :: llos                             ! length of LoS
    real(8), dimension(8) :: ci8                ! node contribution in a block
    real(8) :: wci8                             ! weight of ci8 if multiplication (shared faces/vertices)
    !real(8), dimension(:), allocatable :: ci    ! node contribution in the grid
    !real(8) :: li                               ! line of sight length
    integer :: ni                               ! number of lines of sight per bin
!    real(8), dimension(:), allocatable :: cosi  ! cos(elevation) factor

    ! outputs
    real(8), dimension(:,:), allocatable :: m   ! matrix of node contributions
    real(8), dimension(:), allocatable :: l     ! averaged rock depth per bin
    real(8), dimension(:), allocatable :: r     ! resulting integrated densities
!    real(8), dimension(:,:,:), allocatable :: mm    ! for vtr files

    !integer, dimension(:), allocatable :: inds
    !real(8) :: ds
    !real(8), dimension(:), allocatable :: s, rs
    !integer, dimension(:), allocatable :: is
    !real(8), dimension(:), allocatable :: rho_int
    !real(8), dimension(:,:), allocatable :: cs, cs2

    integer :: i, j, k

    !real(8), parameter :: pi=3.14159265
    real(8), parameter :: pi = 4*atan(1.0_8)
    real(8), parameter :: nodatavalue=0

!    type(vtr_file_handle) :: fd

!    call omp_set_num_threads(8)

    call get_command_argument(1,inputfile)
    call readparam(inputfile, outputfile, topofile, gridparamfile, densitygridfile, televiewfile, muogdatacorrfile)
write(*,*) gridparamfile, televiewfile, topofile
call flush()

    muogdatacorr = .not. (muogdatacorrfile(1:4) == 'none' &
        .or. muogdatacorrfile(1:4) == 'None' &
        .or. muogdatacorrfile(1:4) == 'NONE')

    ! Read density grid parameter file
!    if (debug) then
!        xmin = 0
!        ymin = 0
!        zmin = -1
!        nx = 2
!        ny = 2
!        nz = 2
!        dxy = 1
!        dz = 1
!    else
    call readgridparam(gridparamfile, xmin, ymin, zmin, nx, ny, nz, dxy, dz)
!    end if
    allocate(x(nx))
    allocate(y(ny))
    allocate(z(nz))
    x = (/(xmin+(i-1)*dxy,i=1,nx)/)
    y = (/(ymin+(i-1)*dxy,i=1,ny)/)
    !zmin = -zmin-(nz-1)*dz
    z = (/(zmin+(i-1)*dz,i=1,nz)/)
    xmax = maxval(x)
    ymax = maxval(y)
    zmax = maxval(z)
write(*,*) xmin, ymin, zmin
write(*,*) nx,ny,nz
write(*,*) dxy, dz
call flush()

    ! Read density grid binary file
    allocate(rho(nx*ny*nz))
!    if (debug) then
!        rho = 100
!    else
    call readdensitygrid(densitygridfile, nx, ny, nz, rho)
!    end if
write(*,*) 'rho', minval(rho), maxval(rho)
call flush()

    ! Read telescope view parameter file
!    if (debug) then
!        nmuog = 1
!    else
    call numberofmuogdata(televiewfile, nmuog)
!    end if
    allocate(xtele(nmuog))
    allocate(ytele(nmuog))
    allocate(ztele(nmuog))
    allocate(dazim(nmuog))
    allocate(delev(nmuog))
    allocate(azim(nmuog))
    allocate(elev(nmuog))
!    allocate(dmax(nmuog))
!    if (debug) then
!        xtele = 0
!        ytele = 0
!        ztele = 0
!        dazim = 5
!        delev = 5
!        azim = 45
!        elev = 35.26!45
!        dmax = 5
!    else
    call readmuogdata(televiewfile, nmuog, xtele, ytele, ztele, dazim, delev, azim, elev)
!    end if
    ztele = -ztele ! car z positif vers le bas

    ! Read topography
!    if (debug) then
!        nxt = 2
!        nyt = 2
!        xmint = 0
!        xmaxt = 1
!        ymint = 0
!        ymaxt = 1
!        dxyt = 1
!    else
    call readtopoparam(topofile, nxt, nyt, xmint, xmaxt, ymint, ymaxt, dxyt)
!    end if
    allocate(xt(nxt))
    allocate(yt(nyt))
    allocate(topo(nyt,nxt))
!    if (debug) then
!        topo = 1
!        xt = (/(xmint+(k-1)*dxyt,k=1,nxt)/)
!        yt = (/(ymint+(k-1)*dxyt,k=1,nyt)/)
!        !topo(2,2) = 0.9
!    else
    call readtopo(topofile, nxt, nyt, xmint, xmaxt, ymint, ymaxt, dxyt, xt, yt, topo)
!    end if
    topo = -topo ! car z positif vers le bas
!    if (debug) then
!        write(*,*) 'XT = ', xt
!        write(*,*) 'YT = ', yt
!        write(*,*) 'TOPO = ', topo
!    end if

    ! Identify blocs limited by 8 nodes

    nxb = nx-1
    nyb = ny-1
    nzb = nz-1
    nb = nxb*nyb*nzb
    allocate(x1b(nb))
    allocate(x2b(nb))
    allocate(y1b(nb))
    allocate(y2b(nb))
    allocate(z1b(nb))
    allocate(z2b(nb))
    allocate(l111(nb))
    allocate(l121(nb))
    allocate(l211(nb))
    allocate(l221(nb))
    allocate(l112(nb))
    allocate(l122(nb))
    allocate(l212(nb))
    allocate(l222(nb))
    do i = 1,nxb
        do j = 1,nyb
            do k = 1,nzb
                x1b(i+(j-1)*nxb+(k-1)*nxb*nyb) = x(i)
                x2b(i+(j-1)*nxb+(k-1)*nxb*nyb) = x(i+1)
                y1b(i+(j-1)*nxb+(k-1)*nxb*nyb) = y(j)
                y2b(i+(j-1)*nxb+(k-1)*nxb*nyb) = y(j+1)
                z1b(i+(j-1)*nxb+(k-1)*nxb*nyb) = z(k)
                z2b(i+(j-1)*nxb+(k-1)*nxb*nyb) = z(k+1)
                ! Corresponding linear indices of the 8 nodes
                l111(i+(j-1)*nxb+(k-1)*nxb*nyb) = i+(j-1)*nx+(k-1)*nx*ny
                l121(i+(j-1)*nxb+(k-1)*nxb*nyb) = i+j*nx+(k-1)*nx*ny
                l211(i+(j-1)*nxb+(k-1)*nxb*nyb) = i+1+(j-1)*nx+(k-1)*nx*ny
                l221(i+(j-1)*nxb+(k-1)*nxb*nyb) = i+1+j*nx+(k-1)*nx*ny
                l112(i+(j-1)*nxb+(k-1)*nxb*nyb) = i+(j-1)*nx+k*nx*ny
                l122(i+(j-1)*nxb+(k-1)*nxb*nyb) = i+j*nx+k*nx*ny
                l212(i+(j-1)*nxb+(k-1)*nxb*nyb) = i+1+(j-1)*nx+k*nx*ny
                l222(i+(j-1)*nxb+(k-1)*nxb*nyb) = i+1+j*nx+k*nx*ny
            end do
        end do
    end do

    !CHECK => OK
    !open(60, file = trim(outputfile)//'_blocs.txt')
    !do k = 1,nxb*nyb*nzb
    !    write(60,*) l111(k), l121(k), l211(k), l221(k), l112(k), l122(k), l212(k), l222(k)
    !end do
    !close(60)
    !open(60, file = trim(outputfile)//'_blxyz.txt')
    !do k = 1,nxb*nyb*nzb
    !    write(60,*) x1b(k), x2b(k), y1b(k), y2b(k), z1b(k), z2b(k)
    !end do
    !close(60)

    ! Determine maximum length for line of sight smax

    dist(1) = maxval(sqrt((xtele-xmin)**2+(ytele-ymin)**2+(ztele-zmin)**2))
    dist(2) = maxval(sqrt((xtele-xmax)**2+(ytele-ymin)**2+(ztele-zmin)**2))
    dist(3) = maxval(sqrt((xtele-xmin)**2+(ytele-ymax)**2+(ztele-zmin)**2))
    dist(4) = maxval(sqrt((xtele-xmax)**2+(ytele-ymax)**2+(ztele-zmin)**2))
    dist(5) = maxval(sqrt((xtele-xmin)**2+(ytele-ymin)**2+(ztele-zmax)**2))
    dist(6) = maxval(sqrt((xtele-xmax)**2+(ytele-ymin)**2+(ztele-zmax)**2))
    dist(7) = maxval(sqrt((xtele-xmin)**2+(ytele-ymax)**2+(ztele-zmax)**2))
    dist(8) = maxval(sqrt((xtele-xmax)**2+(ytele-ymax)**2+(ztele-zmax)**2))
    smax = maxval(dist)+dxy
write(*,*) 'Maximum line of sight length smax = ', smax
call flush()

    ! Initialize output variables

    allocate(l(nmuog))
    allocate(m(nmuog,nx*ny*nz))
    !allocate(mm(nx,ny,nz))
    !allocate(ci(nx*ny*nz))
    allocate(clos(nx*ny*nz))
    l = 0
    m = 0

    ! Loop over solid angles

    !$OMP PARALLEL DEFAULT(NONE), &
        !$OMP SHARED(nmuog, azim, elev, dazim, delev), &
        !$OMP SHARED(dxy, dz, smax, xmin, xmax, ymin, ymax, zmin, zmax), &
        !$OMP SHARED(xtele, ytele, ztele, dxyt, nxt, nyt, xt, yt, topo), &
        !$OMP SHARED(nb, x1b, x2b, y1b, y2b, z1b, z2b), &
        !$OMP SHARED(l111, l211, l121, l221, l112, l212, l122, l222), &
        !$OMP SHARED(l, m), &
        !$OMP PRIVATE(k, daz, del, azmoy, elmoy, azinf, azsup, elinf, elsup), &
        !$OMP PRIVATE(alostmp, elostmp, nalos, dalos, nelos, delos, alos, elos), &
        !$OMP PRIVATE(i, nlos, ia, ie, az, el, a, b, c), &
        !$OMP PRIVATE(imesh, xi1, yi1, xi2, yi2, zi1, zi2), &
        !$OMP PRIVATE(dslos, nslos, slos, xslos, yslos, zslos, ztopo), &
        !$OMP PRIVATE(btl, iin, iout, nin, nout, xin, yin, zin, i1, i2, t1, t2), &
        !$OMP PRIVATE(alphalos, alphatopo, tin, xout, yout, zout, tout), &
        !$OMP PRIVATE(clos, llos, ib, ibloc, xi1b, yi1b, zi1b, xi2b, yi2b, zi2b), &
        !$OMP PRIVATE(ci8)
    !$OMP DO SCHEDULE(DYNAMIC)
    do k = 1,nmuog

        write(*,*) 'Bin (ie data)', k, '/', nmuog, azim(k), elev(k)
        call flush()

        daz = dazim(k)
        del = delev(k)
        azmoy = azim(k)
        elmoy = elev(k)
        azinf = azmoy - daz/2
        azsup = azmoy + daz/2
        elinf = elmoy - del/2
        elsup = elmoy + del/2

!        if (debug) then
            !write(*,*) 'daz, del = ', daz, del
            !write(*,*) 'azmoy, elmoy = ', azmoy, elmoy
!            write(*,*) 'azinf, azsup = ', azinf, azsup
!            write(*,*) 'elinf, elsup = ', elinf, elsup
!        end if
        
!        ! Find closest neighbors (BRICOLAGE => daz del a ajouter dans les donnees ???)
!        azmoy = azim(k)
!        elmoy = elev(k)
!        az1 = azim(minloc(abs(azmoy-azim), 1, &
!            xtele==xtele(k) .and. ytele==ytele(k) .and. ztele==ztele(k) &
!            .and. azim/=azim(k) .and. elev/=elev(k)))
!        el1 = elev(minloc(abs(elmoy-elev), 1, &
!            xtele==xtele(k) .and. ytele==ytele(k) .and. ztele==ztele(k) &
!            .and. azim/=azim(k) .and. elev/=elev(k)))
!        az2 = azim(minloc(abs(azmoy-azim), 1, &
!            xtele==xtele(k) .and. ytele==ytele(k) .and. ztele==ztele(k) &
!            .and. azim/=azim(k) .and. elev/=elev(k) .and. azim/=az1))
!        el2 = elev(minloc(abs(elmoy-elev), 1, &
!            xtele==xtele(k) .and. ytele==ytele(k) .and. ztele==ztele(k) &
!            .and. azim/=azim(k) .and. elev/=elev(k) .and. elev/=el1))
!
!        if (az1 > azmoy) then
!            azsup = az1
!            if (az2 < azmoy) then
!                azinf = az2
!            else
!                azsup = min(az1,az2)
!                azinf = azmoy - (azsup-azmoy)
!            end if
!        else
!            azinf = az1
!            if (az2 > azmoy) then
!                azsup = az2
!            else
!                azinf = max(az1,az2)
!                azsup = azmoy + (azmoy-azinf)
!            end if
!        end if
!
!        daz = azsup - azinf
!        del = elsup - elinf
!
!        write(*,*) 'VERIF angle steps', daz, del 
!        write(*,*) azinf, azsup
!        write(*,*) elinf, elsup
!        call flush()

        ! Determine optimal line of sight (los) spacings thaz thel
        alostmp = 2*atan(dxy/25/smax)*180/pi
        elostmp = 2*atan(dz/25/smax)*180/pi
!        alostmp = 0.2 ! in degrees
!        elostmp = 0.2 ! in degrees
        nalos = ceiling(daz/alostmp)
        dalos = daz/nalos
        nelos = ceiling(del/elostmp)
        delos = del/nelos
        !write(*,*) 'SUBLOS', alostmp, elostmp, nalos, dalos, nelos, delos
        allocate(alos(nalos))
        allocate(elos(nelos))
        alos = (/(azinf+dalos/2+(i-1)*dalos,i=1,nalos)/)
        elos = (/(elinf+delos/2+(i-1)*delos,i=1,nelos)/)

!        if (debug) then
!            write(*,*) 'nalos, nelos = ', nalos, nelos
!        end if

!        allocate(cosi(nalos*nelos))
!        cosi = 0
        nlos = 0 ! number of contributing line of sight per bin

!        write(*,*) 'Line of sight spacing', dalos, delos
!        call flush()
    
!        ci = 0

        ! Loop over lines of sight
        do ia = 1,nalos
            do ie = 1,nelos

               !write(*,*) 'LOS (azim,elev) =', alos(ia), elos(ie)

               ! from here: angles in rad
               az = alos(ia)*pi/180
               el = elos(ie)*pi/180

               ! Coefficients directeurs de la droite
               a = sin(az)*cos(el)
               b = cos(az)*cos(el)
               c = -sin(el) ! car z positif vers le bas

!               if (debug) then
!                   write(*,*) 'LoS: ', alos(ia), elos(ia), a, b, c
 !              end if

               ! Submesh limits
               !write(*,*) xmin, xmax, ymin, ymax, zmin, zmax
               !write(*,*) xtele(k), ytele(k), ztele(k), a, b, c
               !call flush()
               ! minval((/zmax,ztele(k)/)): we ignore volume below detector (hence topography behind it)
               call interscube(xmin,xmax,ymin,ymax,zmin,minval((/zmax,ztele(k)/)), &
                   xtele(k),ytele(k),ztele(k),a,b,c,imesh,xi1,yi1,zi1,xi2,yi2,zi2)
               !if (debug) then
               !    write(*,*) 'MESH INTER = ', imesh
               !end if
               if (.not. imesh) then ! if los outside mesh => go to the next los
                   cycle
               end if
!write(*,*) 'imesh = ', imesh
!write(*,*) xi1, yi1, zi1, xi2, yi2, zi2
!call flush()
        
               ! Find topography crossings
               dslos = dxyt ! integration step = topography resolution
               !nslos = smax/dslos + 1 ! in meters, nb of samples along path
               nslos = sqrt((xi2-xi1)**2+(yi2-yi1)**2+(zi2-zi1)**2)/dslos + 1 ! in meters, nb of samples along path
               allocate(slos(nslos))
               allocate(xslos(nslos))
               allocate(yslos(nslos))
               allocate(zslos(nslos))
               allocate(ztopo(nslos))
               slos = (/((i-1)*dslos,i=1,nslos)/)
               xslos = xi1+slos*sin(az)*cos(el)
               yslos = yi1+slos*cos(az)*cos(el)
               zslos = +(zi1-slos*sin(el)) ! car z positif vers le bas
               !do i = 1,nslos
               !    write(*,*) 'LOS COORD', xslos(i), yslos(i), zslos(i)
               !end do
               call bilinint(nxt,nyt,xt,yt,topo,nslos,xslos,yslos,ztopo,nodatavalue)
!               if (debug) then
!                   write(*,*) 'xslos = ', xslos
!                   write(*,*) 'yslos = ', yslos
!                   write(*,*) 'ztopo = ', ztopo
!               end if

               allocate(btl(nslos)) ! below topography level
               btl = (ztopo < zslos) ! car z positif vers le bas
!               if (debug) then
!                   write(*,*) 'ztopo = ', ztopo
!                   write(*,*) 'zslos = ', zslos
!                   write(*,*) 'btl = ', btl
!               end if
!write(*,*) ztopo
!call flush()
               allocate(iin(floor(nslos/2.0)))
               allocate(iout(floor(nslos/2.0)))
               nin = 0
               nout = 0
               if (btl(1)) then
                   nin = nin+1
                   iin(nin) = 1
                   write(*,*) 'WARNING: line of sight initiates inside topography'
                   write(*,*) 'line of sight: ', alos(ia), elos(ie)
                   write(*,*) 'coordinates: ', xslos(1), yslos(1), zslos(1)
                   write(*,*) 'topography: ', ztopo(1)
               end if
               do i = 2,nslos
                   if (btl(i) .and. (.not. btl(i-1))) then
                      nin = nin+1
                      iin(nin) = i
!                      if (debug) then
!                          write(*,*) 'IN', xslos(i), yslos(i), zslos(i)
!                      end if
                   elseif ((.not. btl(i)) .and. btl(i-1)) then
                      nout = nout+1
                      iout(nout) = i
!                      if (debug) then
!                          write(*,*) 'OUT', xslos(i), yslos(i), zslos(i)
!                      end if
                   end if
               end do
               if (btl(nslos)) then
                   nout = nout+1
                   iout(nout) = nslos
                   write(*,*) 'WARNING: line of sight finishes inside topography'
                   write(*,*) 'line of sight: ', alos(ia), elos(ie)
                   write(*,*) 'coordinates: ', xslos(nslos), yslos(nslos), zslos(nslos)
                   write(*,*) 'topography: ', ztopo(nslos)
               end if

!               if (debug) then
!                   write(*,*) 'Nin, Nout = ', nin, nout
!                   call flush()
!               end if

               allocate(xin(nin))
               allocate(yin(nin))
               allocate(zin(nin))
               !if (btl(1)) then
               !    xin(1) = xslos(1)
               !end if
               do i = 1,nin
                   i2 = iin(i)
                   if (i2 == 1) then
                       xin(1) = xslos(1)
                       yin(1) = yslos(1)
                       zin(1) = zslos(1)
                   else
                       i1 = i2-1
                       t1 = sqrt((xslos(i1)-xi1)**2+(yslos(i1)-yi1)**2)
                       t2 = sqrt((xslos(i2)-xi1)**2+(yslos(i2)-yi1)**2)
                       alphalos = -(zslos(i2)-zslos(i1))/(t2-t1)
                       alphatopo = -(ztopo(i2)-ztopo(i1))/(t2-t1)
                       tin = -(ztopo(i1)-zslos(i1))/(alphalos-alphatopo)+t1
                       zin(i) = (zslos(i2)-zslos(i1))*(ztopo(i1)-zslos(i1)) &
                           /(zslos(i2)-zslos(i1)-ztopo(i2)+ztopo(i1))+zslos(i1)
                       xin(i) = xi1 + tin*sin(az)
                       yin(i) = yi1 + tin*cos(az)
                   end if
               end do
               allocate(xout(nout))
               allocate(yout(nout))
               allocate(zout(nout))
               do i = 1,nout
!                   if (debug) then
!                       write(*,*) 'CHECK OUT', i, nout, iout(i), nslos
!                   end if
                   i2 = iout(i)
                   if (i2 == nslos) then
                       xout(nout) = xslos(nslos)
                       yout(nout) = yslos(nslos)
                       zout(nout) = zslos(nslos)
                   else
                       i1 = i2-1
                       t1 = sqrt((xslos(i1)-xi1)**2+(yslos(i1)-yi1)**2)
                       t2 = sqrt((xslos(i2)-xi1)**2+(yslos(i2)-yi1)**2)
                       alphalos = (zslos(i2)-zslos(i1))/(t2-t1)
                       alphatopo = (ztopo(i2)-ztopo(i1))/(t2-t1)
                       tout = (ztopo(i1)-zslos(i1))/(alphalos-alphatopo)+t1
                       zout(i) = (zslos(i2)-zslos(i1))*(ztopo(i1)-zslos(i1)) &
                           /(zslos(i2)-zslos(i1)-ztopo(i2)+ztopo(i1))+zslos(i1)
                       xout(i) = xi1 + tout*sin(az)
                       yout(i) = yi1 + tout*cos(az)
                   end if
               end do

               ! Loop over topography crossings (in/out pairs)
               if (nin /= nout) then
                   write(*,*) 'WARNING: unconsistent number of in and out topography crossings'
                   write(*,*) 'nin = ', nin, 'nout = ', nout
               !else
                   !open(62, file = trim(outputfile)//'_losio.txt')
                   !do i = 1,nin
                   !    write(*,*) 'LOS/TOPO X', alos(ia), elos(ie), xin(i), yin(i), zin(i), xout(i), yout(i), zout(i)
                   !end do
                   !close(62)
               end if

               if (nin > 0) then

                   !li = 0 ! line of sight length
                   !ci = 0 ! contribution
                   !ni = 0
                   clos = 0 ! contribution of nodes to LoS
                   llos = 0 ! length of LoS

                   do i = 1,nin

!                       if (debug) then
!                           write(*,*) '   IN', xin(i), yin(i), zin(i)
!                           write(*,*) '   OUT', xout(i), yout(i), zout(i)
!                       end if

                       ! Loop over blocs
                       do ib = 1,nb
!                           if (debug) then
!                               write(*,*) 'bloc,', ib, '/', nb
!                               write(*,*) 'xin, yin, zin = ', xin(ib), yin(ib), zin(ib)
!                               write(*,*) x1b(ib)<=minval((/xin(i),xout(i)/))-dxy
 !                              write(*,*) x2b(ib)>=maxval((/xin(i),xout(i)/))+dxy
 !                              write(*,*) y1b(ib)<=minval((/yin(i),yout(i)/))-dxy
 !                              write(*,*) y2b(ib)>=maxval((/yin(i),yout(i)/))+dxy
 !                              write(*,*) z1b(ib)<=minval((/zin(i),zout(i)/))-dz
 !                              write(*,*) z2b(ib)>=maxval((/zin(i),zout(i)/))+dz
 !                          end if
!                           test = (alos(ia)==78.75 .and. elos(ie)==3.75 &
!                               .and. x1b(ib)==-100 .and. x2b(ib)==0 &
!                               .and. y1b(ib)==200 .and. y2b(ib)==300 &
!                               .and. z1b(ib)==-100 .and. z2b(ib)==0)
!                           if (test) then
!                               write(*,*) 'BLOC TEST', &
!                                   .not. (x1b(ib)<=minval((/xin(i),xout(i)/))-dxy &
!                                   .or. x2b(ib)>=maxval((/xin(i),xout(i)/))+dxy &
 !                                  .or. y1b(ib)<=minval((/yin(i),yout(i)/))-dxy &
 !                                  .or. y2b(ib)>=maxval((/yin(i),yout(i)/))+dxy &
 !                                  .or. z1b(ib)<=minval((/zin(i),zout(i)/))-dz &
 !                                  .or. z2b(ib)>=maxval((/zin(i),zout(i)/))+dz)
 !                              write(*,*) xin(i), yin(i), zin(i)
 !                              write(*,*) xout(i), yout(i), zout(i)
 !                          end if
                           if (x1b(ib)<minval((/xin(i),xout(i)/))-dxy .or. x2b(ib)>maxval((/xin(i),xout(i)/))+dxy &
                               .or. y1b(ib)<minval((/yin(i),yout(i)/))-dxy .or. y2b(ib)>maxval((/yin(i),yout(i)/))+dxy &
                               .or. z1b(ib)<minval((/zin(i),zout(i)/))-dz .or. z2b(ib)>maxval((/zin(i),zout(i)/))+dz) then
                               cycle ! if bloc outside submesh => go to the next bloc
                           end if

                           call interscube(x1b(ib),x2b(ib),y1b(ib),y2b(ib),z1b(ib),z2b(ib), &
                               xtele(k),ytele(k),ztele(k),a,b,c,ibloc,xi1b,yi1b,zi1b,xi2b,yi2b,zi2b)
!write(*,*) 'INTERSCUBE', x1b(ib), x2b(ib), y1b(ib), y2b(ib), z1b(ib), z2b(ib), ibloc
                           if (.not. ibloc) then
                               cycle ! if los does not cross bloc => go to the next bloc
                           end if
!ICI
!                           ! if along vertice, counted 4 times
!                           if (abs(a)==1 .and. b==0 .and. c==0 &
!                               .and. (minval((/yi1b, yi2b/))==y1b(ib) .or. maxval((/yi1b, yi2b/))==y2b(ib)) &
!                               .and. (minval((/zi1b, zi2b/))==z1b(ib) .or. maxval((/zi1b, zi2b/))==z2b(ib)) &
!                               ) then
!                               wci8 = 0.25
!                           elseif (a==0 .and. abs(b)==1 .and. c==0 &
!                               .and. (minval((/xi1b, xi2b/))==x1b(ib) .or. maxval((/xi1b, xi2b/))==x2b(ib)) &
!                               .and. (minval((/zi1b, zi2b/))==z1b(ib) .or. maxval((/zi1b, zi2b/))==z2b(ib)) &
!                               ) then
!                               wci8 = 0.25
!                           elseif (a==0 .and. b==0 .and. abs(c)==1 &
!                               .and. (minval((/xi1b, xi2b/))==x1b(ib) .or. maxval((/xi1b, xi2b/))==x2b(ib)) &
!                               .and. (minval((/yi1b, yi2b/))==y1b(ib) .or. maxval((/yi1b, yi2b/))==y2b(ib)) &
!                               ) then
!                               wci8 = 0.25
!                           ! if along face, counted twice
!                           elseif (a==0 &
!                               .and. (minval((/xi1b, xi2b/))==x1b(ib) .or. maxval((/xi1b, xi2b/))==x2b(ib)) &
!                               ) then
!                               wci8 = 0.5
!                           end if
                           !write(*,*) 'BLOC', floor(x1b(ib)), floor(x2b(ib)), &
                           !    floor(y1b(ib)), floor(y2b(ib)), &
                           !    floor(z1b(ib)), floor(z2b(ib)), &
                           !    sqrt((xi1b-xi2b)**2+(yi1b-yi2b)**2+(zi1b-zi2b)**2)

                           ! analytic integral within 8 nodes and contributions
                           llos = llos + sqrt((xi1b-xi2b)**2+(yi1b-yi2b)**2+(zi1b-zi2b)**2)
                           !ni = ni+1
                           ci8 = 0
                           call trilinintcoeffs(x1b(ib),x2b(ib),y1b(ib),y2b(ib),z1b(ib),z2b(ib), &
                               xi1b,yi1b,zi1b,xi2b,yi2b,zi2b,xtele(k),ytele(k),ztele(k),az,el,ci8)
                               ! -el car z positif vers le bas
                           ! construction ci(i+(j-1)*2+(k-1)*4) 
                           clos(l111(ib)) = clos(l111(ib)) + ci8(1)
                           clos(l211(ib)) = clos(l211(ib)) + ci8(2)
                           clos(l121(ib)) = clos(l121(ib)) + ci8(3)
                           clos(l221(ib)) = clos(l221(ib)) + ci8(4)
                           clos(l112(ib)) = clos(l112(ib)) + ci8(5)
                           clos(l212(ib)) = clos(l212(ib)) + ci8(6)
                           clos(l122(ib)) = clos(l122(ib)) + ci8(7)
                           clos(l222(ib)) = clos(l222(ib)) + ci8(8)
                           !write(*,*) ci8
                           !call flush()
if (minval(ci8) < 0) then
    write(*,*) 'DEBUG TLIC', minval(ci8), maxval(ci8), sum(ci8), sum(ci8)/sqrt((xi1b-xi2b)**2+(yi1b-yi2b)**2+(zi1b-zi2b)**2)
end if
                  
                       end do ! end loop over blocs 
                   end do ! end loop over segments inside topography along los

!write(*,*) 'li = ', li
!write(*,*) 'az = ', az, 'el = ', el
!write(*,*) xi1, yi1, zi1, xi2, yi2, zi2
!call flush()

                   if (llos == 0) then
                       write(*,*) 'WARNING: llos null'
!                   if (llos > 0) then
                   else
                       nlos = nlos + 1*cos(el)
                       l(k) = l(k) + llos*cos(el)
!write(*,*) nlos, llos, cos(el), sum(clos)
!                   cosi(ia+(ie-1)*nalos) = cos(el)
                       m(k,:) = m(k,:) + clos(:)*cos(el)/llos
!                   end if

                       write(*,*) '    line of sight (az, el, l) =', az*180/pi, el*180/pi, llos*cos(el)
                       !write(*,*) '    DEBUG', k, cos(el), llos, minval(clos(:)), maxval(clos(:))
                   end if

               end if ! if nin > 0

               deallocate(slos)
               deallocate(xslos, yslos, zslos, ztopo)
               deallocate(btl)
               deallocate(iin, iout)
               deallocate(xin, yin, zin, xout, yout, zout)

            end do ! end loop over elev of lines of sight
        end do ! end loop over azim of lines of sight

!write(*,*) nlos, sum(m(k,:)), sum(m(k,:)/nlos)
        if (nlos > 0) then
!write(*,*) '--', sum(cosi), '--'
            !m(k,:) = m(k,:)/sum(cosi)
            !m(k,:) = m(k,:)/l(k) ! we divide by the total length of all LoS in bin
            m(k,:) = m(k,:)/nlos
            write(*,*) sum(m(k,:)), sum(m(k,:)/sum(m(k,:)))
            m(k,:) = m(k,:)/sum(m(k,:)) ! to counteract the very small numerical errors
!        else
!            write(*,*) 'WARNING: no contribution to data ', k, '/', nmuog
            l(k) = l(k)/nlos !(nalos*nelos)
        end if
write(*,*) l(k)

        deallocate(alos, elos)
!        deallocate(cosi)


!    ! write sensivitiy in vtr file
!    call vtr_open_file(prefix=outputfile, fd=fd)
!    call vtr_write_mesh(fd=fd, x=x, y=y, z=-z)
!    mm = reshape(m(k,:),(/nx, ny, nz/))
!    call vtr_write_var(fd=fd, name="sensitivity", field=mm)
!    call vtr_close_file(fd=fd)

    end do ! end loop over data nmuog (bins)
    !$OMP END DO
    !$OMP END PARALLEL

!    deallocate(mm)
!    call vtr_collect_file(fd)

!m(k,:) = m(k,:) + sum(cs2,1)/nsin!*ds!*sum(is)

    deallocate(x,y,z)
    deallocate(xtele, ytele, ztele)
    deallocate(dazim, delev)
!    deallocate(dmax)
    deallocate(xt, yt, topo)
    deallocate(x1b, x2b, y1b, y2b, z1b, z2b)
    deallocate(l111, l121, l221, l112, l122, l212, l211, l222)
!    deallocate(ci)
    deallocate(clos)

    !
    ! Write output files
    !

write(*,*) minval(m), maxval(m), sum(m)
call flush()
    open(90, file = trim(outputfile)//'_muog.bin', form = 'unformatted', access = 'stream')
    write(90) nx, ny, nz, nmuog
    do k=1,nmuog
        write(90) m(k,:)
    end do
    close(90)

    allocate(r(nmuog))
    r = 0
    do k = 1,nmuog
        do i = 1,nx*ny*nz
            r(k) = r(k) + m(k,i)*rho(i)
        end do
        !write(*,*) maxval(rho), sum(abs(m(k,:))), r(k)
    end do

    open(80, file = trim(outputfile)//'_dcal.txt')
    do k = 1,nmuog
        write(80,*) azim(k), elev(k), l(k), r(k)
    end do
    close(80)

    deallocate(azim, elev)
    deallocate(m, l, r)
    deallocate(rho)

end program panda

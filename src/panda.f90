program panda

    use omp_lib

    implicit none

    character(len=100) :: inputfile
    character(len=100) :: outputfile, topofile, gridparamfile, densitygridfile, televiewfile
    character(len=1) :: telestr

    ! density node grid
    real(8) :: xmin,xmax,ymin,ymax,zmin,zmax
    real(8) :: dx,dy,dz, dxy
    integer :: nx,ny,nz
    real(8), dimension(:), allocatable :: x,y,z

    real(8), dimension(:), allocatable :: rho

    ! telescope view rays
    integer :: nmuog
    real(8), dimension(:), allocatable :: xtele, ytele, ztele
    real(8), dimension(:), allocatable :: azim, elev
    real(8) :: az, el
    integer :: iaz, iel

    integer :: nxt, nyt
    real(8) :: xmint, xmaxt, ymint, ymaxt, dxyt
    real(8), dimension(:), allocatable :: xt, yt
    real(8), dimension(:,:), allocatable :: topo

    real(8), dimension(:), allocatable :: xs, ys, zs, zstopo
    integer :: ns, nsin
    integer, dimension(:), allocatable :: inds
    real(8) :: ds
    real(8), dimension(:), allocatable :: s, rs
    integer, dimension(:), allocatable :: is
    real(8), dimension(:), allocatable :: rho_int
    real(8), dimension(:,:), allocatable :: m, cs, cs2
    integer :: km

    integer :: k, i

    real(8), parameter :: pi=3.1415
    real(8), parameter :: nodatavalue=0

    call get_command_argument(1,inputfile)
    call readparam(inputfile, outputfile, topofile, gridparamfile, densitygridfile, televiewfile)
write(*,*) gridparamfile, televiewfile, topofile
call flush()

    ! Read density grid parameter file
    call readgridparam(gridparamfile, xmin, ymin, zmin, nx, ny, nz, dxy, dz)
    allocate(x(nx))
    allocate(y(ny))
    allocate(z(nz))
    x = (/(xmin+(i-1)*dxy,i=1,nx)/)
    y = (/(ymin+(i-1)*dxy,i=1,ny)/)
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
    call readdensitygrid(densitygridfile, nx, ny, nz, rho)
write(*,*) rho(1:10)
call flush()

    ! Read telescope view parameter file
    call numberofmuogdata(televiewfile, nmuog)
!write(*,*) 'number of telescopes = ', ntele
    allocate(xtele(nmuog))
    allocate(ytele(nmuog))
    allocate(ztele(nmuog))
    allocate(azim(nmuog))
    allocate(elev(nmuog))
    call readmuogdata(televiewfile, nmuog, xtele, ytele, ztele, azim, elev)
!    call readteleparam(televiewfile, ntele, xtele, ytele, ztele, azmin, azmax, daz, elmin, elmax, del)
!write(*,*) xtele
!write(*,*) ytele
!write(*,*) ztele
!write(*,*) azmin
!write(*,*) azmax
!write(*,*) daz
!write(*,*) elmin
!write(*,*) elmax
!write(*,*) del

    ! Read topography
    call readtopoparam(topofile, nxt, nyt, xmint, xmaxt, ymint, ymaxt, dxyt)
    allocate(xt(nxt))
    allocate(yt(nyt))
    allocate(topo(nyt,nxt))
    call readtopo(topofile, nxt, nyt, xmint, xmaxt, ymint, ymaxt, dxyt, xt, yt, topo)
    

    ! Compute contribution to average density

    ds = 1 ! in meters, integration step along path
    ns = 10000 ! in meters, number of samples along path
    allocate(s(ns))
    s = (/((i-1)*ds,i=1,ns)/)

    allocate(xs(ns))
    allocate(ys(ns))
    allocate(zs(ns))
    allocate(zstopo(ns))
    allocate(is(ns))
    allocate(cs(ns,nx*ny*nz))
    allocate(rs(ns))

    allocate(m(nmuog,nx*ny*nz))
    m = 0

    !OMP PARALLEL DEFAULT(NONE), &
        !OMP SHARED(nmuog,azim,elev,xtele,ytele,ztele,s,nxt,nyt,xt,yt,topo), &
        !OMP SHARED(ns,xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz,x,y,z,rho,m), &
        !OMP PRIVATE(k,az,el,xs,ys,zs,is,zstopo,nsin,i,inds,rs,cs,cs2)
    !OMP DO SCHEDULE(DYNAMIC)
! DO SCHEDULE(STATIC,nmuog/omp_get_num_threads())
    do k = 1,nmuog
write(*,*) k, '/', nmuog, ' starts'
call flush()

!
!           muon trajectory
!   o  o  o/ o  o
!   o  o  o  o  o
!   o--o-/o--o--o topogrpahy
!   o  o/ o  o  o
!   o  o  o  o  o
!   o /o  o  o  o
!    /
!   /
!  /
! /
! 

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

        az = azim(k)*pi/180
        el = elev(k)*pi/180
                
        xs = xtele(k)+s*sin(az)*cos(el)
        ys = ytele(k)+s*cos(az)*cos(el)
        zs = -(ztele(k)+s*sin(el)) ! car z positif vers le bas
        is = 0

        ! Get altitude from bilinear interpolation of topography
write(*,*) 'bilinint'
call flush()
        call bilinint(nxt,nyt,xt,yt,topo,ns,xs,ys,zstopo,nodatavalue)
write(*,*) 'bilinint ends'
call flush()

        nsin = 0
        !omp parallel default(none), &
        !omp shared(ns,xs,xmin,xmax,ys,ymin,ymax,zs,zmin,zmax,zstopo,is,nsin), &
        !omp private(i)
        !omp do reduction(+:nsin)
        do i = 1,ns
                ! a priori is=0 (ie ignored)
                ! test inside/outside density grid : if inside grid => is=1
                if (xs(i)>=xmin .and. xs(i)<=xmax .and. ys(i)>=ymin .and. ys(i)<=ymax .and. zs(i)>=zmin .and. zs(i)<=zmax) then
                    ! test above/below topography : if above topography => is=0
                    if (-zs(i)<=zstopo(i)) then
                        is(i) = 1
                        nsin = nsin+1
                    end if
                end if
        end do
        !omp end do
        !omp end parallel

        if (nsin /= 0) then
            allocate(inds(nsin))
            inds = 0
            nsin = 0
            !omp parallel default(none), &
            !omp shared(ns,xs,xmin,xmax,ys,ymin,ymax,zs,zmin,zmax,zstopo,is,nsin,inds), &
            !omp private(i)
            !omp do reduction(+:nsin)
            do i = 1,ns
                ! test inside/outside density grid : if inside grid => is=1
                if (xs(i)>=xmin .and. xs(i)<=xmax .and. ys(i)>=ymin .and. ys(i)<=ymax .and. zs(i)>=zmin .and. zs(i)<=zmax) then
                    ! test above/below topography : if above topography => is=0
                    if (-zs(i)<=zstopo(i)) then
                        nsin = nsin+1
                        inds(nsin) = i
                    end if
                end if
            end do
            !omp end do
            !omp end parallel

            !open(90, file = 'indices.bin', form = 'unformatted', access = 'stream')
            !write(90) is
            !close(90)

            rs = 0
            cs = 0
            !call trilinint(nx,ny,nz,x,y,z,rho,count(is==1), &
            !    xs(where(is==1)),ys(where(is==1)),zs(where(is==1)), &
            !    rs(where(is==1)),cs(where(is==1),:),0)
write(*,*) 'trilinintcoeffs'
call flush()
            call trilinintcoeffs(nx,ny,nz,x,y,z,rho,ns,xs,ys,zs,rs,cs,nodatavalue)
write(*,*) 'trilinintcoeffs ends'
call flush()
            allocate(cs2(nsin,nx*ny*nz))
            cs2 = 0
            cs2 = cs(inds,:)
            deallocate(inds)

            m(k,:) = m(k,:) + sum(cs2,1)/nsin!*ds!*sum(is)
            deallocate(cs2)
        end if

        write(*,*) k, '/', nmuog, ' ends'
        call flush()

    end do
    !OMP END DO
    !OMP END PARALLEL

    deallocate(xtele,ytele,ztele,azim,elev)
    deallocate(s,xs,ys,zs,is,cs)

    ! Write output files

write(*,*) minval(m), maxval(m), sum(m)
call flush()
    open(90, file = trim(outputfile)//'_muogs.bin', form = 'unformatted', access = 'stream')
    write(90) nx, ny, nz, nmuog
    do k=1,nmuog
        write(90) m(k,:)
    end do
    close(90)
       
    deallocate(m)

end program panda

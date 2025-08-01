subroutine readtopo(filename, nx, ny, xmin, xmax, ymin, ymax, dxy, x, y, topo)

    implicit none

    character(len=100) :: filename
    real(8) :: xmin, xmax, ymin, ymax, dxy
    integer :: nx, ny
    real(8), dimension(ny,nx) :: topo
    real(4), dimension(nx) :: topo4
    real(8), dimension(nx) :: x
    real(8), dimension(ny) :: y
    real(4) :: tmp
      
    integer :: k

    !
    ! LECTURE TOPOGRAPHY
    !

    open(10, file = filename, form = 'unformatted', access='stream')
    read(10) tmp
    read(10) tmp
    read(10) tmp
    read(10) tmp
    read(10) tmp
    read(10) nx
    read(10) ny
    do k = 1,ny
        read(10) topo4
        topo(k,:) = real(topo4,8)
    end do
    close(10)

    x = (/(xmin+(k-1)*dxy,k=1,nx)/)
    y = (/(ymin+(k-1)*dxy,k=1,ny)/)

!    open(80, file = 'topo_profil_eo.txt')
!    write(80,*) topo(floor(real(nx)/2),:)
!    close(80)
!    open(81, file = 'topo_profil_ns.txt')
!    write(81,*) topo(:,floor(real(ny)/2))
!    close(81)

end subroutine readtopo

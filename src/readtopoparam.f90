subroutine readtopoparam(filename, nx, ny, xmin, xmax, ymin, ymax, dxy)

    implicit none

    character(len=100) :: filename
    real(4) :: xmin4, xmax4, ymin4, ymax4, dxy4
    real(8) :: xmin, xmax, ymin, ymax, dxy
    integer :: nx, ny
      
    !
    ! LECTURE TOPOGRAPHY
    !

    open(10, file = filename, form = 'unformatted', access='stream')
    read(10) ymin4
    read(10) ymax4
    read(10) xmin4
    read(10) xmax4
    read(10) dxy4
    read(10) nx
    read(10) ny
    close(10)

    ymin = real(ymin4,8)
    ymax = real(ymax4,8)
    xmin = real(xmin4,8)
    xmax = real(xmax4,8)
    dxy = real(dxy4,8)

end subroutine readtopoparam

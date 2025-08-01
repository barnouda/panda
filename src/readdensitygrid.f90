subroutine readdensitygrid(filename, nx, ny, nz, grid)

    implicit none

    character(len=100) :: filename
    integer :: nx, ny, nz
    real(8), dimension(nx*ny*nz) :: grid
    real(4), dimension(nx*ny*nz) :: grid4

    integer :: k
      
    !
    ! LECTURE GRILLE
    !

    open(10, file = filename, form = 'unformatted', access='stream')
    read(10) grid4
    close(10)

    grid = real(grid4,8)

    !do k = 1,nx*ny*nz
    !    if (grid(k) /= 0) then
    !       grid(k) = 100
    !    end if
    !end do
    !grid = 100

end subroutine readdensitygrid

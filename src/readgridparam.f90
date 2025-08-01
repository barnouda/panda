subroutine readgridparam(filename, xmin, ymin, zmin, nx, ny, nz, dxy, dz)

    implicit none

    character(len=100) :: filename
    real(8) :: xmin, ymin, zmin
    integer :: nx, ny, nz
    real(8) :: dxy, dz
        
    integer :: tmp
    real(8) :: dx, dy
        
    !
    ! LECTURE DES PARAMETRES DE LA GRILLE
    !

    open(10,file = trim(filename))
    read(10,*) tmp
    read(10,*) xmin, ymin, zmin
    read(10,*) nx, ny, nz
    read(10,*) dx, dy, dz
    close(10)

    if (dx /= dy) then
        write(*,*) 'ERROR: dx/=dy IN ', trim(filename), ' INPUT FILE!'
        stop
    end if

    dxy = dx

end subroutine readgridparam

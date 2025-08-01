subroutine readmuogdata2(filename, nmuog, xtele, ytele, ztele, dazim, delev, azim, elev)

    implicit none

    character(len=100) :: filename
    integer :: nmuog
    real(8), dimension(nmuog) :: xtele, ytele, ztele
    real(8), dimension(nmuog) :: dazim, delev
    real(8), dimension(nmuog) :: azim, elev!, dmax ! dmax=bin to search for

    integer :: k
        
    !
    ! LECTURE DES PARAMETRES DE VISEE DU TELESCOPE A MUONS
    !

    open(10,file = trim(filename))
    do k = 1,nmuog
        read(10,*) xtele(k), ytele(k), ztele(k), azim(k), elev(k), dazim(k), delev(k)
        ! CONE
        !dazim(k) = 5
        !delev(k) = 5
        !dmax(k) = 12 !TODO test
        ! CUBE
        !dazim(k) = 45
        !delev(k) = 45
        ! PDD
        !dazim(k) = 1
        !delev(k) = 1
    end do
    close(10)

end subroutine readmuogdata2

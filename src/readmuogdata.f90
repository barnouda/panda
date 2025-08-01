subroutine readmuogdata(filename, nmuog, xtele, ytele, ztele, azim, elev)

    implicit none

    character(len=100) :: filename
    integer :: nmuog
    real(8), dimension(nmuog) :: xtele, ytele, ztele
    real(8), dimension(nmuog) :: azim, elev

    integer :: k
        
    !
    ! LECTURE DES PARAMETRES DE VISEE DU TELESCOPE A MUONS
    !

    open(10,file = trim(filename))
    do k = 1,nmuog
        read(10,*) xtele(k), ytele(k), ztele(k), azim(k), elev(k)
    end do
    close(10)

end subroutine readmuogdata

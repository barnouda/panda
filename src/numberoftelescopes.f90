subroutine numberoftelescopes(filename, ntele)

    implicit none

    character(len=100) :: filename
    integer :: ntele

    !
    ! LECTURE DES PARAMETRES DE VISEE DU TELESCOPE A MUONS
    !

    open(10,file = trim(filename))
    read(10,*) ntele
    close(10)

end subroutine numberoftelescopes

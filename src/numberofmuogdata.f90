subroutine numberofmuogdata(filename, nmuog)

    implicit none

    character(len=100) :: filename
    integer :: nmuog

    character(len=200) :: tmp
    integer :: io
        
    !
    ! LECTURE FICHIER POUR COUNTER LE NOMBRE DE LIGNES = DE DONNEES
    !

    open(10,file = filename)
    nmuog = 0
    do
        read(10,*,iostat = io) tmp
        if (io < 0) exit
        nmuog = nmuog+1
    end do
    close(10)

end subroutine numberofmuogdata

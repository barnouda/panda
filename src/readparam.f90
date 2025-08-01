subroutine readparam(filename, outputfile, topofile, gridparamfile, densitygridfile, televiewfile, muogdatacorrfile)
!> Read parameter file
!! \param[in]   filename        Name of the input file containing parameters
!! \param[out]  outputfile      Prefix name for the output files
!! \param[out]  gridparamfile   
!! \param[out]  densitygridfile 
!! \param[out]  televiewfile  
!! \param[out]  muogdatacorrfile

    implicit none

    character(len=100) :: filename
    character(len=100) :: outputfile, topofile, gridparamfile, densitygridfile, televiewfile
    character(len=100) :: muogdatacorrfile

    !
    ! LECTURE DES PARAMETRES
    !

    open(10,file = trim(filename))
    read(10,*) outputfile
    read(10,*) topofile
    read(10,*) gridparamfile
    read(10,*) densitygridfile
    read(10,*) televiewfile
    read(10,*) muogdatacorrfile
    close(10)

end subroutine readparam

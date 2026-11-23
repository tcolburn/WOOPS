!-----------------------------------------------------------------------------------------
SUBROUTINE writeStuff1
USE files
IMPLICIT NONE
   WRITE(*,*) 'INPUT PARAMETERS READ FROM INPUT FILE :'
   WRITE(*,*) 'PSF FILE USED                  .... ', psffile
   WRITE(*,*) 'PDB FILE USED                  .... ', pdbfile
   WRITE(*,*) 'DCD FILE USED                  .... ', dcdfile
END SUBROUTINE writeStuff1
!-----------------------------------------------------------------------------------------

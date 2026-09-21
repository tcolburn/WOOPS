!-----------------------------------------------------------------------------------------
! READS PDBS GENERATED FROM PSFGEN
! THE HEADER DOES NOT CONTAIN BOX INFORMATION
SUBROUTINE pdbread(filenm,fileInt,num)
USE coords
USE atoms
IMPLICIT NONE
DOUBLE PRECISION :: x, y, z, cor, bta
INTEGER :: itmp1, itmp2, itmp3
INTEGER :: i, num, atomind
INTEGER :: fileInt
CHARACTER(LEN=60)  :: filenm
CHARACTER(LEN=6)   :: typ1, typ2, typ3
CHARACTER(LEN=60)  :: format1, format2
CHARACTER(LEN=200) :: dummyheader

   format1 = "(A4,1X,A6,1X,A4,1X,A4,A1,I4,4X,3(F8.3),2(F6.2),5X,A5,A2)"
   format2 = "(A4,1X,I7,1X,A4,1X,A4,A1,I4,4X,3(F8.3),2(F6.2),5X,A5,A2)"

   OPEN(UNIT=fileInt,FILE=filenm,FORM='FORMATTED')
   READ(fileInt,*) dummyheader
   DO i=1, num
!      READ(fileInt,format1) typ1, pdbaindex(i), pdbatomname(i), pdbresname(i), typ2, pdbresid(i), &
!                            xini(i), yini(i), zini(i), coor(i), beta(i), segname(i), atomname2(i)
      READ(fileInt,format1) typ1, pdbaindexC(i), pdbatomname(i), pdbresname(i), typ2, pdbresid(i), &
                            xini(i), yini(i), zini(i), coor(i), beta(i), segname(i), atomname2(i)
      pdbaindex(i) = i
!      WRITE(*,format2) typ1, pdbaindex(i), pdbatomname(i), pdbresname(i), typ2, pdbresid(i), &
!                       xini(i), yini(i), zini(i), coor(i), beta(i), segname(i), atomname2(i)
   END DO
   CLOSE(fileInt)

RETURN
END SUBROUTINE pdbread
!-----------------------------------------------------------------------------------------

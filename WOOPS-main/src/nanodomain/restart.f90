!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
SUBROUTINE restartIn(rsInt,rsFileInt,rsFileName,onconf)
USE p1p2
IMPLICIT NONE
INTEGER :: rsInt, rsFileInt, onconf
CHARACTER(LEN=60) :: rsFileName
INTEGER :: i,j,k
REAL*8 :: x,y,z
   ! restart information contained in noName file in main, here rsFileName
   OPEN(UNIT=rsFileInt, FILE=rsFileName, FORM='FORMATTED')
      READ(rsFileInt,916) rsInt, onconf
      ! assign starting values for restart information here
      ! this is designed for trajectories over multiple files
      ! noInt=1 is for initial trajectory file while noInt!=1
      ! is for the following files.
      ! The read and write statements here and below should have the same format
   CLOSE(rsFileInt)
RETURN
916   FORMAT(2(I8))
917   FORMAT(3(F20.10),I20,2(E20.10))
918   FORMAT(3(E20.10))
END SUBROUTINE restartIn
!-------------------------------------------------------------------------------
SUBROUTINE restartOut(rsInt,rsFileInt,rsFileName,inconf)
USE p1p2
IMPLICIT NONE
INTEGER :: rsInt, rsFileInt, inconf
CHARACTER(LEN=60) :: rsFileName
INTEGER :: i,j,k
REAL*8 :: x,y,z
   OPEN(UNIT=rsFileInt, FILE=rsFileName, FORM='FORMATTED')
      WRITE(rsFileInt,916) rsInt, inconf
      ! write starting values for restart information here
      ! this is designed for trajectories over multiple files
      ! noInt=1 is for initial trajectory file while noInt!=1
      ! is for the following files.
      ! The read and write statements here and above should have the same format
   CLOSE(rsFileInt)
RETURN
916   FORMAT(2(I8))
917   FORMAT(3(F20.10),I20,2(E20.10))
918   FORMAT(3(E20.10))
END SUBROUTINE restartOut
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

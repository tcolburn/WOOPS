!-------------------------------------------------------------------------------
! The Kirkwood-only executable analyzes each DCD leg independently.  restartOut
! is retained only to keep the established accumRESTART sidecar behavior used by
! existing setup directories.  restartIn is retained for source compatibility,
! although main.f90 explicitly requires noInt=1.
!-------------------------------------------------------------------------------
SUBROUTINE restartIn(rsInt,rsFileInt,rsFileName,onconf)
IMPLICIT NONE
INTEGER :: rsInt, rsFileInt, onconf
CHARACTER(LEN=60) :: rsFileName
   OPEN(UNIT=rsFileInt, FILE=rsFileName, FORM='FORMATTED')
      READ(rsFileInt,916) rsInt, onconf
   CLOSE(rsFileInt)
RETURN
916 FORMAT(2(I8))
END SUBROUTINE restartIn
!-------------------------------------------------------------------------------
SUBROUTINE restartOut(rsInt,rsFileInt,rsFileName,inconf)
IMPLICIT NONE
INTEGER :: rsInt, rsFileInt, inconf
CHARACTER(LEN=60) :: rsFileName
   OPEN(UNIT=rsFileInt, FILE=rsFileName, FORM='FORMATTED')
      WRITE(rsFileInt,916) rsInt, inconf
   CLOSE(rsFileInt)
RETURN
916 FORMAT(2(I8))
END SUBROUTINE restartOut
!-------------------------------------------------------------------------------

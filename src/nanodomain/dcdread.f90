!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
SUBROUTINE getnframes(filenm,fileInt,nframesdcd)
USE atoms
IMPLICIT NONE
INTEGER :: fileInt
CHARACTER(LEN=60) :: filenm
INTEGER :: i, natomdcd, nframesdcd

REAL             :: dummyr
INTEGER          :: dummyi
CHARACTER(LEN=4) :: dummyc

   OPEN(UNIT=fileInt,FILE=filenm,STATUS='OLD',FORM='UNFORMATTED')
   READ(fileInt) dummyc, nframesdcd, (dummyi,i=1,8), dummyr, (dummyi,i=1,9)
   READ(fileInt) dummyi, dummyr
   READ(fileInt) natomdcd
   CLOSE(fileInt)

   IF (natom.NE.natomdcd) THEN
      WRITE(*,*) 'SOMETHING IS WRONG ...'
      WRITE(*,*) 'THE NUMBER OF ATOMS READ FROM THE'
      WRITE(*,*) 'PSF FILE DOES NOT MATCH THE DCD FILE.'
      STOP
   END IF

RETURN
END SUBROUTINE getnframes
!-----------------------------------------------------------------------------------------
SUBROUTINE dcdread(filenm,fileInt)
USE atoms
USE coords
USE framesNum
IMPLICIT NONE
INTEGER :: fileInt
CHARACTER(LEN=60) :: filenm

REAL             :: dummyr
INTEGER          :: dummyi
CHARACTER(LEN=4) :: dummyc
DOUBLE PRECISION, DIMENSION(1:6) :: d
!DOUBLE PRECISION, DIMENSION(1:6) :: ddum
DOUBLE PRECISION :: ddum
REAL, DIMENSION(1:natom) :: rdum

INTEGER :: i, j, natomdcd, nframesdcd, nmem
LOGICAL, PARAMETER :: DCDUnitCell=.TRUE.

   OPEN(UNIT=fileInt,FILE=filenm,STATUS='OLD',FORM='UNFORMATTED')
   READ(fileInt) dummyc, nframesdcd, (dummyi,i=1,8), dummyr, (dummyi,i=1,9)
   READ(fileInt) dummyi, dummyr
   READ(fileInt) natomdcd

   nmem=1
   DO i=1,nframes
!   DO i=1,10
      IF (i.GE.startFrame.AND.i.LE.endFrame) THEN
         IF (DCDUnitCell) THEN
            READ(fileInt) (d(j), j=1, 6)
            box(nmem,1) = d(1); box(nmem,2) = d(3); box(nmem,3) = d(6)
         ENDIF
         READ(fileInt) (rxf(nmem,j),j=1,natom)
         READ(fileInt) (ryf(nmem,j),j=1,natom)
         READ(fileInt) (rzf(nmem,j),j=1,natom)
         nmem=nmem+1
      ELSE
!         WRITE(*,*) 'SKIPPING FRAME i = ', i
         IF (DCDUnitCell) THEN
            READ(fileInt) (ddum, j=1, 6)
         ENDIF
         READ(fileInt) (rdum(j),j=1,natom)
         READ(fileInt) (rdum(j),j=1,natom)
         READ(fileInt) (rdum(j),j=1,natom)
      END IF
   END DO

   CLOSE(fileInt)
 
END SUBROUTINE dcdread
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

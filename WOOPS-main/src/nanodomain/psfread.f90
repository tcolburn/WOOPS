!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
SUBROUTINE getnum(filenm,fileInt,istart,num)
IMPLICIT NONE
INTEGER :: fileInt, num, i, iwhere, istart
CHARACTER(LEN=60) :: filenm
CHARACTER(LEN=20) :: tmpChar1, tmpChar2

   iwhere = 0
   OPEN(UNIT=fileInt,FILE=filenm,FORM='FORMATTED',STATUS='OLD')

   ! Search for !NATOM.  The atom records begin immediately after that line.
   DO i=1,1000
      READ(fileInt,*,END=900) tmpChar1, tmpChar2
      IF (TRIM(tmpChar2).EQ.'!NATOM') THEN
         READ(tmpChar1,*) num
         iwhere = i
         EXIT
      END IF
   END DO

900 CONTINUE
   CLOSE(fileInt)

   IF (iwhere.EQ.0) THEN
      WRITE(*,*) 'ERROR: COULD NOT FIND !NATOM IN PSF FILE: ', TRIM(filenm)
      STOP
   END IF

   istart = iwhere
RETURN
END SUBROUTINE getnum
!-----------------------------------------------------------------------------------------
SUBROUTINE psfread(filenm,fileInt,istart,num)
USE atoms
USE topology
IMPLICIT NONE
INTEGER :: i, k, num, istart, itmp1
INTEGER :: fileInt, ios
INTEGER, ALLOCATABLE, DIMENSION(:) :: bondValues
CHARACTER(LEN=60) :: filenm
CHARACTER(LEN=20) :: tmpChar1, tmpChar2
CHARACTER(LEN=5) :: typ1
CHARACTER(LEN=512) :: line
LOGICAL :: foundBonds

   OPEN(UNIT=fileInt,FILE=filenm,FORM='FORMATTED',STATUS='OLD')

   ! Advance through the !NATOM line.
   DO i=1,istart
      READ(fileInt,*) tmpChar1, tmpChar2
   END DO

   ! Read PSF atom metadata.
   DO i=1,num
      READ(fileInt,*) aindex(i), typ1, resid(i), resname(i), &
                      atomname(i), atomtype(i), charge(i), mass(i), itmp1
   END DO

   ! Locate and read the PSF bond list.  Bond atom numbers use the same
   ! one-based PSF atom IDs stored in aindex(:).
   foundBonds = .FALSE.
   nBonds = 0

   DO
      READ(fileInt,'(A)',IOSTAT=ios) line
      IF (ios.NE.0) EXIT

      IF (INDEX(line,'!NBOND').GT.0) THEN
         READ(line,*,IOSTAT=ios) nBonds
         IF (ios.NE.0 .OR. nBonds.LT.0) THEN
            WRITE(*,*) 'ERROR: COULD NOT PARSE !NBOND LINE IN PSF FILE.'
            WRITE(*,*) TRIM(line)
            STOP
         END IF
         foundBonds = .TRUE.
         EXIT
      END IF
   END DO

   IF (.NOT.foundBonds) THEN
      WRITE(*,*) 'ERROR: COULD NOT FIND !NBOND SECTION IN PSF FILE: ', TRIM(filenm)
      STOP
   END IF

   IF (ALLOCATED(bondAtom1)) DEALLOCATE(bondAtom1)
   IF (ALLOCATED(bondAtom2)) DEALLOCATE(bondAtom2)
   ALLOCATE(bondAtom1(1:nBonds),bondAtom2(1:nBonds))

   IF (nBonds.GT.0) THEN
      ALLOCATE(bondValues(1:2*nBonds))
      READ(fileInt,*,IOSTAT=ios) bondValues

      IF (ios.NE.0) THEN
         WRITE(*,*) 'ERROR: FAILED WHILE READING PSF BOND ATOM NUMBERS.'
         STOP
      END IF

      DO k=1,nBonds
         bondAtom1(k) = bondValues(2*k-1)
         bondAtom2(k) = bondValues(2*k)
      END DO

      DEALLOCATE(bondValues)
   END IF

   CLOSE(fileInt)

   WRITE(*,*) 'NUMBER OF PSF BONDS READ: ', nBonds
RETURN
END SUBROUTINE psfread
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------
! Streaming trajectory loop: reads and analyzes one DCD frame at a time.
! Includes a tqdm-like progress bar with elapsed time, frame rate, and ETA.
!--------------------------------------------------------------------------------------------
SUBROUTINE mainLoop(inconf,oldnconf,frames)

USE files
USE coords
USE atoms
USE input, ONLY : betaFlg
USE framesNum
USE beta1
USE p1p2

IMPLICIT NONE

INTEGER :: inconf, oldnconf, frames
INTEGER :: i, j, nframesdcd, natomdcd, nselected, nfFinal
INTEGER :: noInt, nmem
INTEGER :: dummyi
REAL :: dummyr
CHARACTER(LEN=4) :: dummyc
DOUBLE PRECISION, DIMENSION(1:6) :: d
REAL, ALLOCATABLE, DIMENSION(:) :: rdum
LOGICAL, PARAMETER :: DCDUnitCell=.TRUE.

! Progress-bar variables
INTEGER :: clockStart, clockNow, clockRate, clockMax
INTEGER :: updateEvery, barWidth, filled
INTEGER :: elapsedSec, etaSec
INTEGER :: elapsedHour, elapsedMin, elapsedRemSec
INTEGER :: etaHour, etaMin, etaRemSec
DOUBLE PRECISION :: elapsedTime, frameRate, etaTime, fraction
CHARACTER(LEN=40) :: progressBar
CHARACTER(LEN=200) :: progressLine


noInt=0
IsaveInt=0
ftrack=0
nselected=0
nfFinal=oldnconf

! Update the progress display once every 1,000 analyzed frames.
updateEvery=1000
barWidth=40

ALLOCATE(rdum(1:natom))

OPEN(UNIT=dcdInt,FILE=dcdfile,STATUS='OLD',FORM='UNFORMATTED')

READ(dcdInt) dummyc, nframesdcd, (dummyi,j=1,8), &
             dummyr, (dummyi,j=1,9)

READ(dcdInt) dummyi, dummyr
READ(dcdInt) natomdcd

IF (natomdcd.NE.natom) THEN
   WRITE(*,*) 'SOMETHING IS WRONG ...'
   WRITE(*,*) 'THE NUMBER OF ATOMS READ FROM THE PSF FILE DOES NOT MATCH THE DCD FILE.'
   STOP
END IF

OPEN(UNIT=fIntp1p2,FILE=p1p2File)

WRITE(*,*)
WRITE(*,*) 'STREAMING TRAJECTORY FROM DCD FILE ...'
WRITE(*,*) 'ANALYZING DCD ONE FRAME AT A TIME ...'
WRITE(*,*)

CALL SYSTEM_CLOCK(clockStart,clockRate,clockMax)


DO i=1,nframesdcd

   IF (DCDUnitCell) THEN
      READ(dcdInt) (d(j),j=1,6)
   END IF

   READ(dcdInt) (rdum(j),j=1,natom)

   IF (i.GE.startFrame.AND.i.LE.endFrame) THEN
      rx(:)=DBLE(rdum(:))
   END IF

   READ(dcdInt) (rdum(j),j=1,natom)

   IF (i.GE.startFrame.AND.i.LE.endFrame) THEN
      ry(:)=DBLE(rdum(:))
   END IF

   READ(dcdInt) (rdum(j),j=1,natom)

   IF (i.GE.startFrame.AND.i.LE.endFrame) THEN
      rz(:)=DBLE(rdum(:))
   END IF


   IF (i.GE.startFrame.AND.i.LE.endFrame) THEN

      nselected=nselected+1
      inconf=oldnconf+nselected

      boxx=d(1)
      boxy=d(3)
      boxz=d(6)

      boxxSq=boxx*boxx
      boxySq=boxy*boxy
      boxzSq=boxz*boxz

      IF (betaFlg.EQ.1) THEN
         CALL popBetaArrays(natom)
      END IF

      IF (inconf.EQ.1) THEN
         IsaveInt=1
      ELSE
         IsaveInt=0
      END IF

      ftrack=inconf

      CALL p1p2Main

      noInt=noInt+1

      IF (nselected.EQ.frames) THEN
         CALL restartOut(noInt,noFileInt,noFileName,inconf)
      END IF


      !-----------------------------------------------------------------
      ! Progress display
      !-----------------------------------------------------------------
      IF (MOD(nselected,updateEvery).EQ.0 .OR. &
          nselected.EQ.frames) THEN

         CALL SYSTEM_CLOCK(clockNow)

         elapsedTime=DBLE(clockNow-clockStart)/DBLE(clockRate)

         fraction=DBLE(nselected)/DBLE(frames)
         fraction=MAX(0.0D0,MIN(1.0D0,fraction))

         IF (elapsedTime.GT.0.0D0) THEN
            frameRate=DBLE(nselected)/elapsedTime
         ELSE
            frameRate=0.0D0
         END IF

         IF (nselected.GT.0) THEN
            etaTime=elapsedTime*DBLE(frames-nselected) &
                    /DBLE(nselected)
         ELSE
            etaTime=0.0D0
         END IF

         elapsedSec=MAX(0,NINT(elapsedTime))
         etaSec=MAX(0,NINT(etaTime))

         elapsedHour=elapsedSec/3600
         elapsedMin=MOD(elapsedSec,3600)/60
         elapsedRemSec=MOD(elapsedSec,60)

         etaHour=etaSec/3600
         etaMin=MOD(etaSec,3600)/60
         etaRemSec=MOD(etaSec,60)

         filled=INT(fraction*DBLE(barWidth))
         filled=MAX(0,MIN(barWidth,filled))

         progressBar=REPEAT('=',filled) // &
                     REPEAT('-',barWidth-filled)

         WRITE(progressLine,1001) progressBar,              &
                                  100.0D0*fraction,          &
                                  nselected,                 &
                                  frames,                    &
                                  frameRate,                 &
                                  elapsedHour,               &
                                  elapsedMin,                &
                                  elapsedRemSec,             &
                                  etaHour,                   &
                                  etaMin,                    &
                                  etaRemSec

         WRITE(*,'(A)',ADVANCE='NO') CHAR(13) // &
                                      TRIM(progressLine)

         CALL FLUSH(6)

      END IF
      !-----------------------------------------------------------------


      frameNum=frameNum+1
      nfFinal=inconf

   END IF

END DO


! Move subsequent terminal output onto a new line.
WRITE(*,*)

CLOSE(dcdInt)


IF (nselected.NE.frames) THEN

   WRITE(*,*) 'ERROR: NUMBER OF FRAMES ANALYZED DOES NOT MATCH INPUT.'
   WRITE(*,*) 'EXPECTED: ',frames
   WRITE(*,*) 'ANALYZED: ',nselected

   STOP

END IF


p1Avg(:)=p1Acc(:)/DBLE(nselected)
p2Avg(:)=p2Acc(:)/DBLE(nselected)
p21Avg(:)=p21Acc(:)/DBLE(nselected)


DO j=1,nbeta1

   WRITE(fIntp1p2,1000) 'DATA DI : ',                &
                         nfFinal,                    &
                         aidbt1(j),                  &
                         ridbt1(j),                  &
                         anamebt1(j),                &
                         rnamebt1(j),                &
                         chrgbt1(j),                 &
                         p1Avg(j),                   &
                         p2Avg(j),                   &
                         p21Avg(j)

END DO


CLOSE(fIntp1p2)

DEALLOCATE(rdum)

RETURN


1000 FORMAT(A10,I10,I5,3X,3A10,4F16.10)

1001 FORMAT('[',A40,'] ',F6.2,'%  ',I0,'/',I0,       &
            '  ',F10.2,' frames/s  elapsed ',        &
            I2.2,':',I2.2,':',I2.2,                 &
            '  ETA ',                                &
            I2.2,':',I2.2,':',I2.2)


END SUBROUTINE mainLoop
!--------------------------------------------------------------------------------------------

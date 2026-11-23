!----------------------------------------------------------------------------------------
! Dedicated scalar Kirkwood analysis derived from the WOOPS streaming reader.
!----------------------------------------------------------------------------------------
PROGRAM main
USE input
USE files
USE water
USE atoms
USE framesNum
IMPLICIT NONE
INTEGER :: inconf, oldnconf
INTEGER :: noInt

   CALL copyright

   ! Preserve the established 11-line WOOPS input format so existing prepared
   ! Kirkwood leg directories do not need new input files.
   READ(*,'(A)') psffile
   READ(*,'(A)') pdbfile
   READ(*,'(A)') dcdfile
   READ(*,'(A)') noFileName
   READ(*,*) nresInput
   READ(*,*) nPAtomInput
   READ(*,*) wattype
   READ(*,*) ntrajFrame
   READ(*,*) startFrame
   READ(*,*) endFrame
   READ(*,*) noInt

   betaFlg = 1

   IF (endFrame-startFrame+1.NE.ntrajFrame) THEN
      WRITE(*,*) 'THE START AND END FRAME NUMBERS DO NOT MATCH THE TOTAL'
      WRITE(*,*) 'NUMBER OF FRAMES REQUESTED FOR ANALYSIS.'
      STOP
   END IF

   ! This dedicated program treats each DCD leg as an independent statistical
   ! sample.  Multi-leg combination belongs in post-processing, so no trajectory
   ! accumulators are silently carried through the legacy restart mechanism.
   IF (noInt.NE.1) THEN
      WRITE(*,*) 'ERROR: DEDICATED KIRKWOOD ANALYSIS REQUIRES INPUT LINE 11 = 1.'
      WRITE(*,*) 'ANALYZE EACH DCD LEG INDEPENDENTLY AND COMBINE LEGS AFTERWARD.'
      STOP
   END IF

   frameNum = 0
   oldnconf = 0

   CALL writeStuff1
   CALL setup

   WRITE(*,*) '----------------------------------------------------------------------'
   WRITE(*,*) 'STREAMING TRAJECTORY FROM DCD FILE ...'
   WRITE(*,*) '----------------------------------------------------------------------'
   CALL mainLoop(inconf,oldnconf,ntrajFrame)

   CALL exitProg
   WRITE(*,*) 'COMPLETED SUCCESSFULLY ...'
END PROGRAM main
!-----------------------------------------------------------------------------------------

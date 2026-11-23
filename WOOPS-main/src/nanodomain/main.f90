!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!  (c) 2015, Arizona Board of Regents on behalf of ASU
!----------------------------------------------------------------------------------------
!  THIS PROGRAM WAS WRITTEN BY DANIEL MARTIN AT ASU'S CENTER FOR BIOLOGICAL PHYSICS
!  IN DMITRY MATYUSHOV'S GROUP.
!  - TIP4P WATER IS NOT SUPPORTED
!  - ONLY 3 POINT WATER MODELS ARE SUPPORTED
!  - TESTED ONLY FOR TIP3P AND SPCE
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
PROGRAM main
USE input
USE files
USE water
USE coords
USE atoms
USE framesNum
USE beta1
USE beta2
USE p1p2
IMPLICIT NONE
INTEGER :: inconf, oldnconf
INTEGER :: noInt

INTEGER, PARAMETER :: res1Int    = 81
CHARACTER(LEN=60)  :: res1Name

   CALL copyright

   READ(*,'(A)') psffile
   READ(*,'(A)') pdbfile
   READ(*,'(A)') dcdfile
   READ(*,'(A)') noFileName
   READ(*,*) nres    ! number of residues in the solute
   READ(*,*) nPAtom  ! number of atoms in the solute
   READ(*,*) wattype ! 0=TIP3P, 1=TIP4P (NOT SUPPORTED), 2=SPCE
   READ(*,*) ntrajFrame
   READ(*,*) startFrame
   READ(*,*) endFrame
   READ(*,*) noInt
   ALLOCATE(p1Acc(1:nPAtom),p2Acc(1:nPAtom),p1Avg(1:nPAtom),p2Avg(1:nPAtom))
   ALLOCATE(p21Acc(1:nPAtom),p21Avg(1:nPAtom))
   ALLOCATE(rDifSave(1:nPAtom))
   ALLOCATE(nwSave(1:nPAtom),nwSaveAcc(1:nPAtom),aidbt2Save(1:nPAtom))
   ALLOCATE(nWatProtAtom(1:nPAtom))
   betaFlg=1
   IF (endFrame-startFrame+1.NE.ntrajFrame) THEN
      WRITE(*,*) 'THE START AND END FRAME NUMBERS DO NOT '
      WRITE(*,*) 'MATCH THE TOTAL FRAMES BEING ANALYZED ...'
      WRITE(*,*) 'EXITING THE PROGRAM ...'
      STOP
   END IF
   frameNum = 0
   CALL writeStuff1
   IF (noInt.EQ.1) THEN
      WRITE(*,*) 'THIS IS THE STARTING POINT OF THE ANALYSIS!'
   ELSE
      WRITE(*,*) 'CONTINUING FROM PREVIOUS ANALYSIS DIRECTORY ...'
   END IF
   CALL setup
   p1Acc(:)=0.0d0; p2Acc(:)=0.0d0; p21Acc(:)=0.0d0
   p1Avg(:)=0.0d0; p2Avg(:)=0.0d0; p21Avg(:)=0.0d0
   nwSave(:)=0
   nwSaveAcc(:)=0
   aidbt2Save(:)=0
   nWatProtAtom(:)=0
   oldnconf = 0
   ! restart information contained in noName file
   IF (noInt.NE.1) THEN
      CALL restartIn(noInt,noFileInt,noFileName,oldnconf)
   END IF
   WRITE(*,*) '----------------------------------------------------------------------'
   WRITE(*,*) 'STREAMING TRAJECTORY FROM DCD FILE ...'
   WRITE(*,*) '----------------------------------------------------------------------'
   !######################
   !## loop over frames ##
   !######################
   ! loop over frames in the dcd file.
   ! the loop is from '1+oldnconf' to 'nframes+oldnconf'
   ! oldconf is the number of frames already analysed.
   WRITE(*,*) 'ANALYZING DCD ONE FRAME AT A TIME ... '
   CALL mainLoop(inconf,oldnconf,ntrajFrame)
   !##########################
   !## end loop over frames ##
   !##########################
   DEALLOCATE(p1Acc,p2Acc,p1Avg,p2Avg)
   DEALLOCATE(p21Acc,p21Avg)
   DEALLOCATE(rDifSave)
   DEALLOCATE(nwSave,nwSaveAcc,aidbt2Save)
   DEALLOCATE(nWatProtAtom)
   CALL exitProg
   WRITE(*,*) 'COMPLETED SUCCESSFULLY ...'
END PROGRAM main
!-----------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------
! Streaming trajectory loop for the dedicated scalar Kirkwood analysis.
!--------------------------------------------------------------------------------------------
SUBROUTINE mainLoop(inconf,oldnconf,frames)
USE files
USE coords
USE atoms
USE input, ONLY : betaFlg
USE framesNum
USE kirkwood
IMPLICIT NONE
INTEGER :: inconf, oldnconf, frames
INTEGER :: i, j, region, nframesdcd, natomdcd, nselected, nfFinal
INTEGER :: dummyi, ios, validInt
REAL :: dummyr
CHARACTER(LEN=4) :: dummyc
CHARACTER(LEN=256) :: ioMsg
DOUBLE PRECISION, DIMENSION(1:6) :: d
REAL, ALLOCATABLE, DIMENSION(:) :: rdum
LOGICAL, PARAMETER :: DCDUnitCell=.TRUE.

   ftrack = 0
   nselected = 0
   nfFinal = oldnconf

   IF (writeTimeSeries .AND. timeSeriesStride.LT.1) THEN
      WRITE(*,*) 'ERROR: timeSeriesStride MUST BE >= 1.'
      STOP
   END IF

   IF (writeDiagnosticSums .AND. diagnosticStride.LT.1) THEN
      WRITE(*,*) 'ERROR: diagnosticStride MUST BE >= 1.'
      STOP
   END IF

   ALLOCATE(rdum(1:natom))

   OPEN(UNIT=dcdInt,FILE=dcdfile,STATUS='OLD',FORM='UNFORMATTED', &
        IOSTAT=ios,IOMSG=ioMsg)
   CALL checkDcdIO(ios,ioMsg,'opening DCD',0)

   READ(dcdInt,IOSTAT=ios,IOMSG=ioMsg) dummyc, nframesdcd, &
        (dummyi,j=1,8), dummyr, (dummyi,j=1,9)
   CALL checkDcdIO(ios,ioMsg,'reading DCD header record 1',0)

   READ(dcdInt,IOSTAT=ios,IOMSG=ioMsg) dummyi, dummyr
   CALL checkDcdIO(ios,ioMsg,'reading DCD header record 2',0)

   READ(dcdInt,IOSTAT=ios,IOMSG=ioMsg) natomdcd
   CALL checkDcdIO(ios,ioMsg,'reading DCD atom-count record',0)

   IF (natomdcd.NE.natom) THEN
      WRITE(*,*) 'ERROR: PSF AND DCD ATOM COUNTS DO NOT MATCH.'
      WRITE(*,*) 'PSF / DCD: ', natom, natomdcd
      STOP
   END IF

   IF (writeTimeSeries) THEN
      DO region=1,nKirkwoodRegions
         OPEN(UNIT=timeSeriesUnit(region),FILE=TRIM(timeSeriesFile(region)), &
              STATUS='REPLACE',ACTION='WRITE')
         WRITE(timeSeriesUnit(region),'(A)') &
            '# Compact full Kirkwood time series; valid=0 frames are excluded from averages.'
         WRITE(timeSeriesUnit(region),'(A,A)') '# region: ',TRIM(regionLabel(region))
         WRITE(timeSeriesUnit(region),'(A)') &
            '# CSV columns: frame,valid,N,normal_coherence,g_k,g_perp,g_parallel,decomposition_residual'
         WRITE(timeSeriesUnit(region),'(A)') &
            'frame,valid,N,normal_coherence,g_k,g_perp,g_parallel,decomposition_residual'
      END DO
   END IF

   IF (writeDiagnosticSums) THEN
      DO region=1,nKirkwoodRegions
         OPEN(UNIT=diagnosticUnit(region),FILE=TRIM(diagnosticFile(region)), &
              STATUS='REPLACE',ACTION='WRITE')
         WRITE(diagnosticUnit(region),'(A)') &
            '# Optional Kirkwood audit sums.  These are not required for routine post-processing.'
         WRITE(diagnosticUnit(region),'(A,A)') '# region: ',TRIM(regionLabel(region))
         WRITE(diagnosticUnit(region),'(A)') &
            '# CSV columns: frame,valid,N,sum_ux,sum_uy,sum_uz,sum_p,' // &
            'sum_upar_x,sum_upar_y,sum_upar_z,sum_nx,sum_ny,sum_nz'
         WRITE(diagnosticUnit(region),'(A)') &
            'frame,valid,N,sum_ux,sum_uy,sum_uz,sum_p,' // &
            'sum_upar_x,sum_upar_y,sum_upar_z,sum_nx,sum_ny,sum_nz'
      END DO
   END IF

   DO i=1,nframesdcd
      IF (DCDUnitCell) THEN
         READ(dcdInt,IOSTAT=ios,IOMSG=ioMsg) (d(j),j=1,6)
         CALL checkDcdIO(ios,ioMsg,'reading DCD unit-cell record',i)
      END IF

      READ(dcdInt,IOSTAT=ios,IOMSG=ioMsg) (rdum(j),j=1,natom)
      CALL checkDcdIO(ios,ioMsg,'reading DCD X coordinates',i)
      IF (i.GE.startFrame.AND.i.LE.endFrame) rx(:)=DBLE(rdum(:))

      READ(dcdInt,IOSTAT=ios,IOMSG=ioMsg) (rdum(j),j=1,natom)
      CALL checkDcdIO(ios,ioMsg,'reading DCD Y coordinates',i)
      IF (i.GE.startFrame.AND.i.LE.endFrame) ry(:)=DBLE(rdum(:))

      READ(dcdInt,IOSTAT=ios,IOMSG=ioMsg) (rdum(j),j=1,natom)
      CALL checkDcdIO(ios,ioMsg,'reading DCD Z coordinates',i)
      IF (i.GE.startFrame.AND.i.LE.endFrame) rz(:)=DBLE(rdum(:))

      IF (i.GE.startFrame.AND.i.LE.endFrame) THEN
         nselected = nselected + 1
         inconf = oldnconf + nselected

         boxx = d(1)
         boxy = d(3)
         boxz = d(6)
         boxxSq = boxx*boxx
         boxySq = boxy*boxy
         boxzSq = boxz*boxz

         IF (betaFlg.EQ.1) CALL popBetaArrays(natom)

         ftrack = inconf
         CALL kirkwoodMain

         IF (writeTimeSeries) THEN
            IF (MOD(nselected-1,timeSeriesStride).EQ.0) THEN
               DO region=1,nKirkwoodRegions
                  validInt = MERGE(1,0,regionFrameValid(region))
                  WRITE(timeSeriesUnit(region),1001) &
                     inconf, validInt, regionFrameN(region), &
                     regionFrameNormalCoherence(region), regionFrameGK(region), &
                     regionFrameGPerp(region), regionFrameGParallel(region), &
                     regionFrameResidual(region)
               END DO
            END IF
         END IF

         IF (writeDiagnosticSums) THEN
            IF (MOD(nselected-1,diagnosticStride).EQ.0) THEN
               DO region=1,nKirkwoodRegions
                  validInt = MERGE(1,0,regionFrameValid(region))
                  WRITE(diagnosticUnit(region),1003) &
                     inconf, validInt, regionFrameN(region), &
                     regionFrameSumU(1,region), regionFrameSumU(2,region), &
                     regionFrameSumU(3,region), regionFrameSumP(region), &
                     regionFrameSumUPar(1,region), regionFrameSumUPar(2,region), &
                     regionFrameSumUPar(3,region), regionFrameSumN(1,region), &
                     regionFrameSumN(2,region), regionFrameSumN(3,region)
               END DO
            END IF
         END IF

         IF (MOD(nselected,1000).EQ.0) THEN
            WRITE(*,*) 'FRAME) ', inconf, nselected, ' DCD FRAME ', i
         END IF

         frameNum = frameNum + 1
         nfFinal = inconf
      END IF
   END DO

   CLOSE(dcdInt)
   IF (writeTimeSeries) THEN
      DO region=1,nKirkwoodRegions
         CLOSE(timeSeriesUnit(region))
      END DO
   END IF

   IF (writeDiagnosticSums) THEN
      DO region=1,nKirkwoodRegions
         CLOSE(diagnosticUnit(region))
      END DO
   END IF

   IF (nselected.NE.frames) THEN
      WRITE(*,*) 'ERROR: NUMBER OF FRAMES ANALYZED DOES NOT MATCH INPUT.'
      WRITE(*,*) 'EXPECTED: ', frames, ' ANALYZED: ', nselected
      STOP
   END IF

   CALL restartOut(1,noFileInt,noFileName,nfFinal)

   DO region=1,nKirkwoodRegions
      CALL writeKirkwoodSummary(region,nselected)
   END DO

   DEALLOCATE(rdum)
RETURN
1001 FORMAT(I0,',',I0,',',I0,5(',',ES22.14))
1003 FORMAT(I0,',',I0,',',I0,10(',',ES22.14))
END SUBROUTINE mainLoop
!--------------------------------------------------------------------------------------------
SUBROUTINE checkDcdIO(ios,ioMsg,stage,dcdFrame)
USE files, ONLY : dcdfile
IMPLICIT NONE
INTEGER, INTENT(IN) :: ios, dcdFrame
CHARACTER(LEN=*), INTENT(IN) :: ioMsg, stage

   IF (ios.EQ.0) RETURN

   WRITE(*,*) '----------------------------------------------------------------------'
   WRITE(*,*) 'ERROR: DCD READ FAILED.'
   WRITE(*,*) 'FILE:  ', TRIM(dcdfile)
   WRITE(*,*) 'STAGE: ', TRIM(stage)
   IF (dcdFrame.GT.0) WRITE(*,*) 'DCD FRAME: ', dcdFrame
   WRITE(*,*) 'IOSTAT: ', ios
   WRITE(*,*) 'IOMSG:  ', TRIM(ioMsg)
   WRITE(*,*) 'THIS IS A FILE/STREAM READ FAILURE; TEST THE DCD WITH A RAW READ'
   WRITE(*,*) '(FOR EXAMPLE dd if=<DCD> of=/dev/null) BEFORE BLAMING THE ESTIMATOR.'
   WRITE(*,*) '----------------------------------------------------------------------'
   STOP 2
END SUBROUTINE checkDcdIO
!--------------------------------------------------------------------------------------------
SUBROUTINE writeKirkwoodSummary(region,totalFrames)
USE kirkwood
IMPLICIT NONE
INTEGER, INTENT(IN) :: region, totalFrames
INTEGER :: validFrames, zeroFrames, ratioValid
DOUBLE PRECISION :: meanN, validFraction, meanCoherence
DOUBLE PRECISION :: gK, gPerp, gParallel, residual
DOUBLE PRECISION :: ratioPerp, ratioParallel, devPerp, devParallel
DOUBLE PRECISION, PARAMETER :: tiny = 1.0d-15

   validFrames = regionValidFrames(region)
   zeroFrames = totalFrames - validFrames

   meanN = 0.0d0
   validFraction = 0.0d0
   meanCoherence = 0.0d0
   gK = 0.0d0
   gPerp = 0.0d0
   gParallel = 0.0d0
   residual = 0.0d0
   ratioPerp = 0.0d0
   ratioParallel = 0.0d0
   devPerp = 0.0d0
   devParallel = 0.0d0
   ratioValid = 0

   IF (validFrames.GT.0) THEN
      meanN = regionNAcc(region)/DBLE(validFrames)
      validFraction = DBLE(validFrames)/DBLE(totalFrames)
      meanCoherence = regionNormalCoherenceAcc(region)/DBLE(validFrames)
      gK = regionGKAcc(region)/DBLE(validFrames)
      gPerp = regionGPerpAcc(region)/DBLE(validFrames)
      gParallel = regionGParallelAcc(region)/DBLE(validFrames)
      residual = regionResidualAcc(region)/DBLE(validFrames)

      IF (ABS(gK).GT.tiny) THEN
         ratioValid = 1
         ratioPerp = gPerp/gK
         ratioParallel = gParallel/gK
         devPerp = ratioPerp - 1.0d0/3.0d0
         devParallel = ratioParallel - 2.0d0/3.0d0
      END IF
   END IF

   OPEN(UNIT=summaryUnit(region),FILE=TRIM(summaryFile(region)), &
        STATUS='REPLACE',ACTION='WRITE')
   WRITE(summaryUnit(region),'(A,A)') '# region: ',TRIM(regionLabel(region))
   WRITE(summaryUnit(region),'(A)') &
      '# Frames with N=0 are excluded from all bracket averages.'
   WRITE(summaryUnit(region),'(A)') &
      '# decomposition_residual = <g_k - g_perp - g_parallel>.'
   WRITE(summaryUnit(region),'(A)') &
      '# For a common normal, decomposition_residual should vanish exactly.'
   WRITE(summaryUnit(region),'(A)') &
      '# For an isotropic ensemble with a common normal, g_perp/g_k -> 1/3 and g_parallel/g_k -> 2/3.'
   WRITE(summaryUnit(region),'(A)') &
      '# These common-normal identities are diagnostics, not constraints for the curved full interface.'
   WRITE(summaryUnit(region),'(A)') &
      '# total_frames valid_frames zero_occupancy_frames ratio_valid mean_N ' // &
      'valid_fraction mean_normal_coherence g_k g_perp g_parallel ' // &
      'decomposition_residual g_perp_over_g_k g_parallel_over_g_k ' // &
      'g_perp_ratio_minus_one_third g_parallel_ratio_minus_two_thirds'
   WRITE(summaryUnit(region),1002) totalFrames, validFrames, zeroFrames, ratioValid, &
        meanN, validFraction, meanCoherence, gK, gPerp, gParallel, residual, &
        ratioPerp, ratioParallel, devPerp, devParallel
   CLOSE(summaryUnit(region))

   WRITE(*,*) '----------------------------------------------------------------------'
   WRITE(*,*) 'KIRKWOOD SUMMARY: ', TRIM(regionLabel(region))
   WRITE(*,*) 'TOTAL / VALID / ZERO-OCCUPANCY FRAMES: ', &
              totalFrames, validFrames, zeroFrames
   WRITE(*,*) 'MEAN N:                    ', meanN
   WRITE(*,*) 'MEAN NORMAL COHERENCE:     ', meanCoherence
   WRITE(*,*) '<g_K>:                     ', gK
   WRITE(*,*) '<g_perp>:                  ', gPerp
   WRITE(*,*) '<g_parallel>:              ', gParallel
   WRITE(*,*) '<g_K-g_perp-g_parallel>:   ', residual
   IF (ratioValid.EQ.1) THEN
      WRITE(*,*) 'g_perp/g_K:                ', ratioPerp
      WRITE(*,*) 'g_parallel/g_K:            ', ratioParallel
   END IF
RETURN
1002 FORMAT(4(I12,1X),11(ES22.14,1X))
END SUBROUTINE writeKirkwoodSummary
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------

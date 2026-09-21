!-----------------------------------------------------------------------------------------
! Cell-list nanodomain WOOPS estimator.
!
! Per frame:
!   1. Build a periodic linked-cell list of all beta=1 protein atoms.
!   2. For each water oxygen, search only neighboring protein cells.
!   3. In the same candidate pass:
!        a. test whether the oxygen is within 7 A of any protein atom;
!        b. locate its nearest protein heavy atom.
!   4. Calculate P1, P2, and P21 relative to that heavy-atom-to-water vector.
!   5. Average waters assigned to the same heavy atom in the frame.
!   6. Propagate the heavy-atom frame mean to its bonded protein hydrogens.
!   7. Add one atom-level value per frame to the original trajectory accumulators.
!
! The cell width is at least assignmentSearchRadius = 7.0 + 2.5 A.  A water
! selected by proximity to a protein hydrogen must have that hydrogen's parent
! heavy atom inside this radius.  Protein H-parent distances are checked every
! frame, making the cell result identical to the global brute-force definition.
!
! Frames in which an atom's parent heavy atom receives no water contribute zero,
! preserving the original WOOPS normalization by the total number of frames.
!-----------------------------------------------------------------------------------------
SUBROUTINE p1p2Main
USE atoms, ONLY : nbeta1, nbeta2
USE coords, ONLY : boxx, boxy, boxz
USE beta1
USE beta2
USE water, ONLY : qhw
USE p1p2, ONLY : p1Acc, p2Acc, p21Acc, nWatProtAtom
USE nanodomain
USE omp_lib, ONLY : omp_get_thread_num, omp_get_num_threads
IMPLICIT NONE

INTEGER :: nw, j, h, parentHeavy
INTEGER :: nearestHeavy, atomInCell
INTEGER :: ix, iy, iz, jx, jy, jz
INTEGER :: ox, oy, oz, cellId, q
INTEGER :: nCandidateCells
INTEGER :: tid, actualThreads, touchIndex
INTEGER :: errorCode, errorWaterId, errorHeavyId
INTEGER, DIMENSION(27) :: candidateCells
LOGICAL :: cellAlreadyPresent
INTEGER :: newNx, newNy, newNz, newCellCount
REAL*8 :: rxiOx, ryiOx, rziOx
REAL*8 :: rxiH1, ryiH1, rziH1
REAL*8 :: rxiH2, ryiH2, rziH2
REAL*8 :: rDifX, rDifY, rDifZ, rDifSq, rDif
REAL*8 :: bestDifX, bestDifY, bestDifZ
REAL*8 :: minAnySq, minHeavySq, cutoffSq, assignmentRadiusSq
REAL*8 :: rxOH1, ryOH1, rzOH1, rOH1
REAL*8 :: rxOH2, ryOH2, rzOH2, rOH2
REAL*8 :: exOH1, eyOH1, ezOH1
REAL*8 :: exOH2, eyOH2, ezOH2
REAL*8 :: rxiDip, ryiDip, rziDip, magDip
REAL*8 :: exiDip, eyiDip, eziDip
REAL*8 :: dotDipi, roh1n, roh2n
REAL*8 :: rABCrossx, rABCrossy, rABCrossz, rABCrossMag
REAL*8 :: nmCrossx, nmCrossy, nmCrossz, nmCrossMag
REAL*8 :: cChiw, p1, p2, p21
REAL*8 :: lx, ly, lz, invLx, invLy, invLz
REAL*8 :: cellWidthX, cellWidthY, cellWidthZ
REAL*8 :: wrappedX, wrappedY, wrappedZ
REAL*8 :: parentBondSq, maxParentBondSq
REAL*8 :: errorValue1, errorValue2
REAL*8, PARAMETER :: tinyLength = 1.0d-12

   cutoffSq = hydrationCutoff*hydrationCutoff
   assignmentRadiusSq = assignmentSearchRadius*assignmentSearchRadius

   lx = ABS(boxx)
   ly = ABS(boxy)
   lz = ABS(boxz)

   IF (lx.LE.tinyLength .OR. ly.LE.tinyLength .OR. lz.LE.tinyLength) THEN
      WRITE(*,*) 'ERROR: INVALID OR MISSING ORTHORHOMBIC BOX LENGTHS.'
      WRITE(*,*) 'BOX LENGTHS: ', lx, ly, lz
      STOP
   END IF

   invLx = 1.0d0/lx
   invLy = 1.0d0/ly
   invLz = 1.0d0/lz

   ! Verify the geometric bound used to make nearest-heavy assignment local.
   ! This is only O(Nprotein) and is negligible compared with the water search.
   maxParentBondSq = maxProteinHydrogenBond*maxProteinHydrogenBond
   DO j=1,nbeta1
      IF (.NOT.isProteinHeavy(j)) THEN
         parentHeavy = parentHeavyIndex(j)

         rDifX = rxbt1(j) - rxbt1(parentHeavy)
         rDifY = rybt1(j) - rybt1(parentHeavy)
         rDifZ = rzbt1(j) - rzbt1(parentHeavy)

         rDifX = rDifX - lx*DNINT(rDifX*invLx)
         rDifY = rDifY - ly*DNINT(rDifY*invLy)
         rDifZ = rDifZ - lz*DNINT(rDifZ*invLz)

         parentBondSq = rDifX*rDifX + rDifY*rDifY + rDifZ*rDifZ

         IF (parentBondSq.GT.maxParentBondSq) THEN
            WRITE(*,*) 'ERROR: PROTEIN H-PARENT DISTANCE EXCEEDS CELL-SEARCH BOUND.'
            WRITE(*,*) 'HYDROGEN ID / PARENT ID: ', aidbt1(j), aidbt1(parentHeavy)
            WRITE(*,*) 'DISTANCE / BOUND: ', SQRT(parentBondSq), maxProteinHydrogenBond
            STOP
         END IF
      END IF
   END DO

   ! Choose as many cells as possible while keeping every cell width at least
   ! assignmentSearchRadius.  Therefore every relevant atom lies in the water's
   ! own cell or one of the 26 periodic neighboring cells.
   newNx = MAX(1,INT(lx/assignmentSearchRadius))
   newNy = MAX(1,INT(ly/assignmentSearchRadius))
   newNz = MAX(1,INT(lz/assignmentSearchRadius))
   newCellCount = newNx*newNy*newNz

   IF (.NOT.ALLOCATED(proteinCellHead) .OR. &
       newCellCount.NE.proteinCellCount) THEN
      IF (ALLOCATED(proteinCellHead)) DEALLOCATE(proteinCellHead)
      ALLOCATE(proteinCellHead(1:newCellCount))
   END IF

   proteinCellNx = newNx
   proteinCellNy = newNy
   proteinCellNz = newNz
   proteinCellCount = newCellCount

   cellWidthX = lx/DBLE(proteinCellNx)
   cellWidthY = ly/DBLE(proteinCellNy)
   cellWidthZ = lz/DBLE(proteinCellNz)

   proteinCellHead(:) = 0
   proteinCellNext(:) = 0

   ! Rebuild the small protein cell list for the current frame.
   DO j=1,nbeta1
      wrappedX = MODULO(rxbt1(j),lx)
      wrappedY = MODULO(rybt1(j),ly)
      wrappedZ = MODULO(rzbt1(j),lz)

      ix = MIN(proteinCellNx,INT(wrappedX/cellWidthX)+1)
      iy = MIN(proteinCellNy,INT(wrappedY/cellWidthY)+1)
      iz = MIN(proteinCellNz,INT(wrappedZ/cellWidthZ)+1)

      cellId = ix + proteinCellNx*((iy-1) + proteinCellNy*(iz-1))
      proteinCellNext(j) = proteinCellHead(cellId)
      proteinCellHead(cellId) = j
   END DO

   waterCountFrame(:) = 0
   p1FrameSum(:) = 0.0d0
   p2FrameSum(:) = 0.0d0
   p21FrameSum(:) = 0.0d0
   touchedCountThread(:) = 0

   actualThreads = 1
   errorCode = 0
   errorWaterId = 0
   errorHeavyId = 0
   errorValue1 = 0.0d0
   errorValue2 = 0.0d0

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP& PRIVATE(tid,nw,j,nearestHeavy,atomInCell,ix,iy,iz,jx,jy,jz, &
!$OMP&         ox,oy,oz,cellId,q,nCandidateCells,candidateCells, &
!$OMP&         cellAlreadyPresent,touchIndex,rxiOx,ryiOx,rziOx, &
!$OMP&         rxiH1,ryiH1,rziH1,rxiH2,ryiH2,rziH2, &
!$OMP&         rDifX,rDifY,rDifZ,rDifSq,rDif,bestDifX,bestDifY,bestDifZ, &
!$OMP&         minAnySq,minHeavySq,rxOH1,ryOH1,rzOH1,rOH1, &
!$OMP&         rxOH2,ryOH2,rzOH2,rOH2,exOH1,eyOH1,ezOH1, &
!$OMP&         exOH2,eyOH2,ezOH2,rxiDip,ryiDip,rziDip,magDip, &
!$OMP&         exiDip,eyiDip,eziDip,dotDipi,roh1n,roh2n, &
!$OMP&         rABCrossx,rABCrossy,rABCrossz,rABCrossMag, &
!$OMP&         nmCrossx,nmCrossy,nmCrossz,nmCrossMag,cChiw,p1,p2,p21, &
!$OMP&         wrappedX,wrappedY,wrappedZ)
   tid = omp_get_thread_num() + 1

!$OMP SINGLE
   actualThreads = omp_get_num_threads()
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
   DO nw=1,nbeta2,3
      rxiOx = rxbt2(nw)
      ryiOx = rybt2(nw)
      rziOx = rzbt2(nw)

      wrappedX = MODULO(rxiOx,lx)
      wrappedY = MODULO(ryiOx,ly)
      wrappedZ = MODULO(rziOx,lz)

      ix = MIN(proteinCellNx,INT(wrappedX/cellWidthX)+1)
      iy = MIN(proteinCellNy,INT(wrappedY/cellWidthY)+1)
      iz = MIN(proteinCellNz,INT(wrappedZ/cellWidthZ)+1)

      ! Collect unique periodic neighbor cells.  Deduplication matters when a
      ! box dimension contains only one or two cells.
      nCandidateCells = 0
      DO oz=-1,1
         jz = MODULO(iz-1+oz,proteinCellNz) + 1
         DO oy=-1,1
            jy = MODULO(iy-1+oy,proteinCellNy) + 1
            DO ox=-1,1
               jx = MODULO(ix-1+ox,proteinCellNx) + 1
               cellId = jx + proteinCellNx*((jy-1) + proteinCellNy*(jz-1))

               cellAlreadyPresent = .FALSE.
               DO q=1,nCandidateCells
                  IF (candidateCells(q).EQ.cellId) THEN
                     cellAlreadyPresent = .TRUE.
                     EXIT
                  END IF
               END DO

               IF (.NOT.cellAlreadyPresent) THEN
                  nCandidateCells = nCandidateCells + 1
                  candidateCells(nCandidateCells) = cellId
               END IF
            END DO
         END DO
      END DO

      ! One candidate pass performs both the any-protein cutoff test and the
      ! nearest-heavy assignment.  No second heavy-atom loop is needed.
      minAnySq = HUGE(1.0d0)
      minHeavySq = HUGE(1.0d0)
      nearestHeavy = 0
      bestDifX = 0.0d0
      bestDifY = 0.0d0
      bestDifZ = 0.0d0

      DO q=1,nCandidateCells
         atomInCell = proteinCellHead(candidateCells(q))

         DO WHILE (atomInCell.NE.0)
            j = atomInCell

            rDifX = rxiOx - rxbt1(j)
            rDifY = ryiOx - rybt1(j)
            rDifZ = rziOx - rzbt1(j)

            rDifX = rDifX - lx*DNINT(rDifX*invLx)
            rDifY = rDifY - ly*DNINT(rDifY*invLy)
            rDifZ = rDifZ - lz*DNINT(rDifZ*invLz)

            rDifSq = rDifX*rDifX + rDifY*rDifY + rDifZ*rDifZ

            IF (rDifSq.LT.minAnySq) minAnySq = rDifSq

            IF (isProteinHeavy(j)) THEN
               IF (rDifSq.LT.minHeavySq) THEN
                  minHeavySq = rDifSq
                  nearestHeavy = j
                  bestDifX = rDifX
                  bestDifY = rDifY
                  bestDifZ = rDifZ
               END IF
            END IF

            atomInCell = proteinCellNext(atomInCell)
         END DO
      END DO

      IF (minAnySq.GT.cutoffSq) CYCLE

      IF (nearestHeavy.EQ.0 .OR. minHeavySq.GT.assignmentRadiusSq) THEN
!$OMP CRITICAL(nanodomain_error)
         IF (errorCode.EQ.0) THEN
            errorCode = 1
            errorWaterId = aidbt2(nw)
            errorValue1 = SQRT(minAnySq)
            IF (nearestHeavy.NE.0) errorValue2 = SQRT(minHeavySq)
         END IF
!$OMP END CRITICAL(nanodomain_error)
         CYCLE
      END IF

      rDif = SQRT(minHeavySq)
      IF (rDif.LE.tinyLength) THEN
!$OMP CRITICAL(nanodomain_error)
         IF (errorCode.EQ.0) THEN
            errorCode = 2
            errorWaterId = aidbt2(nw)
            errorHeavyId = aidbt1(nearestHeavy)
         END IF
!$OMP END CRITICAL(nanodomain_error)
         CYCLE
      END IF

      ! Minimum-image O-H vectors keep the water geometry continuous even
      ! when a molecule straddles a periodic boundary.
      rxiH1 = rxbt2(nw+1)
      ryiH1 = rybt2(nw+1)
      rziH1 = rzbt2(nw+1)
      rxiH2 = rxbt2(nw+2)
      ryiH2 = rybt2(nw+2)
      rziH2 = rzbt2(nw+2)

      rxOH1 = rxiH1 - rxiOx
      ryOH1 = ryiH1 - ryiOx
      rzOH1 = rziH1 - rziOx
      rxOH1 = rxOH1 - lx*DNINT(rxOH1*invLx)
      ryOH1 = ryOH1 - ly*DNINT(ryOH1*invLy)
      rzOH1 = rzOH1 - lz*DNINT(rzOH1*invLz)

      rxOH2 = rxiH2 - rxiOx
      ryOH2 = ryiH2 - ryiOx
      rzOH2 = rziH2 - rziOx
      rxOH2 = rxOH2 - lx*DNINT(rxOH2*invLx)
      ryOH2 = ryOH2 - ly*DNINT(ryOH2*invLy)
      rzOH2 = rzOH2 - lz*DNINT(rzOH2*invLz)

      rOH1 = SQRT(rxOH1*rxOH1 + ryOH1*ryOH1 + rzOH1*rzOH1)
      rOH2 = SQRT(rxOH2*rxOH2 + ryOH2*ryOH2 + rzOH2*rzOH2)

      IF (rOH1.LE.tinyLength .OR. rOH2.LE.tinyLength) THEN
!$OMP CRITICAL(nanodomain_error)
         IF (errorCode.EQ.0) THEN
            errorCode = 3
            errorWaterId = aidbt2(nw)
         END IF
!$OMP END CRITICAL(nanodomain_error)
         CYCLE
      END IF

      exOH1 = rxOH1/rOH1
      eyOH1 = ryOH1/rOH1
      ezOH1 = rzOH1/rOH1
      exOH2 = rxOH2/rOH2
      eyOH2 = ryOH2/rOH2
      ezOH2 = rzOH2/rOH2

      ! For a neutral three-site model, qO + 2*qH = 0, so the molecular
      ! dipole is obtained directly from the two minimum-image O-H vectors.
      rxiDip = qhw*(rxOH1 + rxOH2)
      ryiDip = qhw*(ryOH1 + ryOH2)
      rziDip = qhw*(rzOH1 + rzOH2)
      magDip = SQRT(rxiDip*rxiDip + ryiDip*ryiDip + rziDip*rziDip)

      IF (magDip.LE.tinyLength) THEN
!$OMP CRITICAL(nanodomain_error)
         IF (errorCode.EQ.0) THEN
            errorCode = 4
            errorWaterId = aidbt2(nw)
         END IF
!$OMP END CRITICAL(nanodomain_error)
         CYCLE
      END IF

      exiDip = rxiDip/magDip
      eyiDip = ryiDip/magDip
      eziDip = rziDip/magDip

      dotDipi = (exiDip*bestDifX + eyiDip*bestDifY + eziDip*bestDifZ)/rDif
      dotDipi = MAX(-1.0d0,MIN(1.0d0,dotDipi))

      roh1n = (exOH1*bestDifX + eyOH1*bestDifY + ezOH1*bestDifZ)/rDif
      roh2n = (exOH2*bestDifX + eyOH2*bestDifY + ezOH2*bestDifZ)/rDif

      ! Retain the original WOOPS water-plane convention.  Reversing the
      ! cross product changes only its sign, while P21 depends on cChiw^2.
      IF (roh1n.GE.0.0d0 .AND. roh2n.LE.0.0d0) THEN
         rABCrossx = eyOH1*ezOH2 - ezOH1*eyOH2
         rABCrossy = ezOH1*exOH2 - exOH1*ezOH2
         rABCrossz = exOH1*eyOH2 - eyOH1*exOH2
      ELSE
         rABCrossx = eyOH2*ezOH1 - ezOH2*eyOH1
         rABCrossy = ezOH2*exOH1 - exOH2*ezOH1
         rABCrossz = exOH2*eyOH1 - eyOH2*exOH1
      END IF

      nmCrossx = (bestDifY*eziDip - bestDifZ*eyiDip)/rDif
      nmCrossy = (bestDifZ*exiDip - bestDifX*eziDip)/rDif
      nmCrossz = (bestDifX*eyiDip - bestDifY*exiDip)/rDif

      nmCrossMag = SQRT(nmCrossx*nmCrossx + nmCrossy*nmCrossy + nmCrossz*nmCrossz)
      rABCrossMag = SQRT(rABCrossx*rABCrossx + rABCrossy*rABCrossy + &
                         rABCrossz*rABCrossz)

      IF (nmCrossMag.LE.tinyLength .OR. rABCrossMag.LE.tinyLength) THEN
         cChiw = 0.0d0
      ELSE
         cChiw = (rABCrossx*nmCrossx + rABCrossy*nmCrossy + &
                  rABCrossz*nmCrossz)/(nmCrossMag*rABCrossMag)
         cChiw = MAX(-1.0d0,MIN(1.0d0,cChiw))
      END IF

      p1 = dotDipi
      p2 = 0.50d0*(3.0d0*dotDipi*dotDipi - 1.0d0)
      p21 = 0.50d0*((1.0d0-dotDipi*dotDipi)*(2.0d0*cChiw*cChiw-1.0d0))

      ! Each thread writes only to its own accumulator column.  The first
      ! observation of a heavy atom records it in a compact touched list, so
      ! the serial reduction does not scan every atom-thread pair each frame.
      IF (waterCountThread(nearestHeavy,tid).EQ.0) THEN
         touchIndex = touchedCountThread(tid) + 1
         touchedCountThread(tid) = touchIndex
         touchedHeavyThread(touchIndex,tid) = nearestHeavy
      END IF

      p1ThreadSum(nearestHeavy,tid) = p1ThreadSum(nearestHeavy,tid) + p1
      p2ThreadSum(nearestHeavy,tid) = p2ThreadSum(nearestHeavy,tid) + p2
      p21ThreadSum(nearestHeavy,tid) = p21ThreadSum(nearestHeavy,tid) + p21
      waterCountThread(nearestHeavy,tid) = waterCountThread(nearestHeavy,tid) + 1
   END DO
!$OMP END DO
!$OMP END PARALLEL

   IF (actualThreads.GT.nOpenMPThreads) THEN
      WRITE(*,*) 'ERROR: OPENMP USED MORE THREADS THAN WERE ALLOCATED.'
      WRITE(*,*) 'USED / ALLOCATED: ', actualThreads, nOpenMPThreads
      STOP
   END IF

   IF (errorCode.NE.0) THEN
      SELECT CASE (errorCode)
      CASE (1)
         WRITE(*,*) 'ERROR: SELECTED WATER HAS NO HEAVY ATOM INSIDE THE'
         WRITE(*,*) 'VALIDATED CELL-SEARCH RADIUS.'
         WRITE(*,*) 'WATER OXYGEN PSF ID: ', errorWaterId
         WRITE(*,*) 'MINIMUM ANY-ATOM DISTANCE: ', errorValue1
         IF (errorValue2.GT.0.0d0) THEN
            WRITE(*,*) 'MINIMUM CANDIDATE HEAVY DISTANCE: ', errorValue2
         END IF
      CASE (2)
         WRITE(*,*) 'ERROR: WATER OXYGEN OVERLAPS ITS NEAREST PROTEIN HEAVY ATOM.'
         WRITE(*,*) 'WATER OXYGEN PSF ID: ', errorWaterId
         WRITE(*,*) 'PROTEIN HEAVY PSF ID: ', errorHeavyId
      CASE (3)
         WRITE(*,*) 'ERROR: ZERO-LENGTH WATER O-H BOND ENCOUNTERED.'
         WRITE(*,*) 'WATER OXYGEN PSF ID: ', errorWaterId
      CASE (4)
         WRITE(*,*) 'ERROR: ZERO WATER DIPOLE MAGNITUDE ENCOUNTERED.'
         WRITE(*,*) 'WATER OXYGEN PSF ID: ', errorWaterId
      END SELECT
      STOP
   END IF

   ! Reduce only the heavy atoms touched by each thread, then clear those
   ! thread-local entries for reuse in the next frame.
   DO tid=1,actualThreads
      DO touchIndex=1,touchedCountThread(tid)
         j = touchedHeavyThread(touchIndex,tid)
         waterCountFrame(j) = waterCountFrame(j) + waterCountThread(j,tid)
         p1FrameSum(j) = p1FrameSum(j) + p1ThreadSum(j,tid)
         p2FrameSum(j) = p2FrameSum(j) + p2ThreadSum(j,tid)
         p21FrameSum(j) = p21FrameSum(j) + p21ThreadSum(j,tid)

         waterCountThread(j,tid) = 0
         p1ThreadSum(j,tid) = 0.0d0
         p2ThreadSum(j,tid) = 0.0d0
         p21ThreadSum(j,tid) = 0.0d0
         touchedHeavyThread(touchIndex,tid) = 0
      END DO
      touchedCountThread(tid) = 0
   END DO

   ! Convert each occupied heavy atom's water sums into one frame value.
   DO h=1,nProteinHeavy
      j = heavyBetaIndex(h)
      IF (waterCountFrame(j).GT.0) THEN
         p1FrameSum(j) = p1FrameSum(j)/DBLE(waterCountFrame(j))
         p2FrameSum(j) = p2FrameSum(j)/DBLE(waterCountFrame(j))
         p21FrameSum(j) = p21FrameSum(j)/DBLE(waterCountFrame(j))
      END IF
   END DO

   ! Heavy atoms are their own parents; bonded hydrogens inherit their parent's
   ! frame value.  Unoccupied parents add zero, preserving total-frame averaging.
   DO j=1,nbeta1
      parentHeavy = parentHeavyIndex(j)

      IF (parentHeavy.LE.0 .OR. parentHeavy.GT.nbeta1) THEN
         WRITE(*,*) 'ERROR: INVALID HEAVY-ATOM PARENT FOR BETA=1 INDEX ', j
         STOP
      END IF

      IF (waterCountFrame(parentHeavy).GT.0) THEN
         p1Acc(j) = p1Acc(j) + p1FrameSum(parentHeavy)
         p2Acc(j) = p2Acc(j) + p2FrameSum(parentHeavy)
         p21Acc(j) = p21Acc(j) + p21FrameSum(parentHeavy)
         nWatProtAtom(j) = nWatProtAtom(j) + 1
      END IF
   END DO
RETURN
END SUBROUTINE p1p2Main
!-----------------------------------------------------------------------------------------

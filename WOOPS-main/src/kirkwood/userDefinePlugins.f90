!-----------------------------------------------------------------------------------------
! Dedicated scalar Kirkwood estimator.
!
! For each analyzed frame, three water ensembles are evaluated simultaneously:
!   1. interface_7A : water O within 7 A of ANY protein atom, including H;
!   2. PRO16_7A    : water O within 7 A of the PRO16 carbonyl O;
!   3. PRO16_10A   : water O within 10 A of the PRO16 carbonyl O.
!
! For every selected water i, the local normal is ALWAYS water-specific:
!
!   n_i = (nearest protein heavy atom -> water O) / distance.
!
! With u_i the unit water dipole,
!
!   p_i              = u_i . n_i
!   u_i,parallel      = u_i - n_i (u_i . n_i)
!
! For a region with instantaneous occupancy N,
!
!   g_K        = |sum_i u_i|^2 / N
!   g_perp     = (sum_i p_i)^2 / N
!   g_parallel = |sum_i u_i,parallel|^2 / N
!
! These are algebraically identical to the corresponding ij double sums but
! require only O(N) work.  Trajectory brackets are taken only over frames with
! N>0.  No atom-resolved Kirkwood map and no P1/P2/P21 calculation are performed.
!-----------------------------------------------------------------------------------------
SUBROUTINE kirkwoodMain
USE atoms, ONLY : nbeta1, nbeta2
USE coords, ONLY : boxx, boxy, boxz
USE beta1
USE beta2
USE water, ONLY : qhw
USE surface_geometry
USE kirkwood
USE omp_lib, ONLY : omp_get_thread_num, omp_get_num_threads
IMPLICIT NONE

INTEGER :: nw, j, parentHeavy
INTEGER :: nearestHeavy, atomInCell
INTEGER :: ix, iy, iz, jx, jy, jz
INTEGER :: ox, oy, oz, cellId, q
INTEGER :: nCandidateCells
INTEGER :: tid, actualThreads, region
INTEGER :: errorCode, errorWaterId, errorHeavyId
INTEGER :: newNx, newNy, newNz, newCellCount
INTEGER, DIMENSION(27) :: candidateCells
LOGICAL :: cellAlreadyPresent
LOGICAL, DIMENSION(nKirkwoodRegions) :: selectedRegion
REAL*8 :: rxiOx, ryiOx, rziOx
REAL*8 :: rxiH1, ryiH1, rziH1
REAL*8 :: rxiH2, ryiH2, rziH2
REAL*8 :: rDifX, rDifY, rDifZ, rDifSq, rDif
REAL*8 :: bestDifX, bestDifY, bestDifZ
REAL*8 :: minAnySq, minHeavySq
REAL*8 :: interfaceCutoffSq, pro16Cutoff7Sq, pro16Cutoff10Sq
REAL*8 :: searchRadiusSq
REAL*8 :: rxOH1, ryOH1, rzOH1, rOH1Sq
REAL*8 :: rxOH2, ryOH2, rzOH2, rOH2Sq
REAL*8 :: rxiDip, ryiDip, rziDip, magDip
REAL*8 :: ux, uy, uz
REAL*8 :: nx, ny, nz, p
REAL*8 :: uparx, upary, uparz
REAL*8 :: lx, ly, lz, invLx, invLy, invLz
REAL*8 :: cellWidthX, cellWidthY, cellWidthZ
REAL*8 :: wrappedX, wrappedY, wrappedZ
REAL*8 :: parentBondSq, maxParentBondSq
REAL*8 :: errorValue1, errorValue2
REAL*8 :: dxPro16, dyPro16, dzPro16, r2Pro16
REAL*8, PARAMETER :: tinyLength = 1.0d-12

   interfaceCutoffSq = interfaceCutoff*interfaceCutoff
   pro16Cutoff7Sq = pro16Cutoff7*pro16Cutoff7
   pro16Cutoff10Sq = pro16Cutoff10*pro16Cutoff10
   searchRadiusSq = proteinSearchRadius*proteinSearchRadius

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

   ! The interface selection includes protein H atoms.  Verify that every H is
   ! close enough to its heavy parent for the 10 A cell-search bound to remain
   ! rigorous when a water is selected by proximity to H.
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
            WRITE(*,*) 'ERROR: PROTEIN H-PARENT DISTANCE EXCEEDS SEARCH BOUND.'
            WRITE(*,*) 'HYDROGEN ID / PARENT ID: ', aidbt1(j), aidbt1(parentHeavy)
            WRITE(*,*) 'DISTANCE / BOUND: ', SQRT(parentBondSq), maxProteinHydrogenBond
            STOP
         END IF
      END IF
   END DO

   ! Each cell width is at least proteinSearchRadius, so every relevant protein
   ! atom is in the water's own cell or one of the 26 periodic neighbor cells.
   newNx = MAX(1,INT(lx/proteinSearchRadius))
   newNy = MAX(1,INT(ly/proteinSearchRadius))
   newNz = MAX(1,INT(lz/proteinSearchRadius))
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

   regionCountThread(:,:) = 0
   regionSumUThread(:,:,:) = 0.0d0
   regionSumPThread(:,:) = 0.0d0
   regionSumUParThread(:,:,:) = 0.0d0
   regionSumNThread(:,:,:) = 0.0d0

   regionFrameN(:) = 0
   regionFrameValid(:) = .FALSE.
   regionFrameSumU(:,:) = 0.0d0
   regionFrameSumP(:) = 0.0d0
   regionFrameSumUPar(:,:) = 0.0d0
   regionFrameSumN(:,:) = 0.0d0
   regionFrameNormalCoherence(:) = 0.0d0
   regionFrameGK(:) = 0.0d0
   regionFrameGPerp(:) = 0.0d0
   regionFrameGParallel(:) = 0.0d0
   regionFrameResidual(:) = 0.0d0

   actualThreads = 1
   errorCode = 0
   errorWaterId = 0
   errorHeavyId = 0
   errorValue1 = 0.0d0
   errorValue2 = 0.0d0

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP& PRIVATE(tid,nw,j,nearestHeavy,atomInCell,ix,iy,iz,jx,jy,jz, &
!$OMP&         ox,oy,oz,cellId,q,nCandidateCells,candidateCells, &
!$OMP&         cellAlreadyPresent,selectedRegion,region, &
!$OMP&         rxiOx,ryiOx,rziOx,rxiH1,ryiH1,rziH1,rxiH2,ryiH2,rziH2, &
!$OMP&         rDifX,rDifY,rDifZ,rDifSq,rDif,bestDifX,bestDifY,bestDifZ, &
!$OMP&         minAnySq,minHeavySq,rxOH1,ryOH1,rzOH1,rOH1Sq, &
!$OMP&         rxOH2,ryOH2,rzOH2,rOH2Sq,rxiDip,ryiDip,rziDip,magDip, &
!$OMP&         ux,uy,uz,nx,ny,nz,p,uparx,upary,uparz, &
!$OMP&         wrappedX,wrappedY,wrappedZ,dxPro16,dyPro16,dzPro16,r2Pro16)
   tid = omp_get_thread_num() + 1

!$OMP SINGLE
   actualThreads = omp_get_num_threads()
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
   DO nw=1,nbeta2,3
      rxiOx = rxbt2(nw)
      ryiOx = rybt2(nw)
      rziOx = rzbt2(nw)

      ! PRO16 membership is independent of the full-interface selection.
      dxPro16 = rxiOx - rxbt1(pro16OIndex)
      dyPro16 = ryiOx - rybt1(pro16OIndex)
      dzPro16 = rziOx - rzbt1(pro16OIndex)
      dxPro16 = dxPro16 - lx*DNINT(dxPro16*invLx)
      dyPro16 = dyPro16 - ly*DNINT(dyPro16*invLy)
      dzPro16 = dzPro16 - lz*DNINT(dzPro16*invLz)
      r2Pro16 = dxPro16*dxPro16 + dyPro16*dyPro16 + dzPro16*dzPro16

      wrappedX = MODULO(rxiOx,lx)
      wrappedY = MODULO(ryiOx,ly)
      wrappedZ = MODULO(rziOx,lz)

      ix = MIN(proteinCellNx,INT(wrappedX/cellWidthX)+1)
      iy = MIN(proteinCellNy,INT(wrappedY/cellWidthY)+1)
      iz = MIN(proteinCellNz,INT(wrappedZ/cellWidthZ)+1)

      ! Deduplicate periodic neighbor cells for boxes with only one or two cells
      ! along a dimension.
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

      ! One protein pass determines both the nearest ANY protein atom (for the
      ! 7 A interface selection) and the nearest protein HEAVY atom (for n_i).
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

      selectedRegion(:) = .FALSE.
      selectedRegion(regionInterface) = (minAnySq.LE.interfaceCutoffSq)
      selectedRegion(regionPro16_7) = (r2Pro16.LE.pro16Cutoff7Sq)
      selectedRegion(regionPro16_10) = (r2Pro16.LE.pro16Cutoff10Sq)

      IF (.NOT.ANY(selectedRegion)) CYCLE

      IF (nearestHeavy.EQ.0 .OR. minHeavySq.GT.searchRadiusSq) THEN
!$OMP CRITICAL(kirkwood_error)
         IF (errorCode.EQ.0) THEN
            errorCode = 1
            errorWaterId = aidbt2(nw)
            IF (minAnySq.LT.HUGE(1.0d0)/2.0d0) errorValue1 = SQRT(minAnySq)
            IF (nearestHeavy.NE.0) errorValue2 = SQRT(minHeavySq)
         END IF
!$OMP END CRITICAL(kirkwood_error)
         CYCLE
      END IF

      rDif = SQRT(minHeavySq)
      IF (rDif.LE.tinyLength) THEN
!$OMP CRITICAL(kirkwood_error)
         IF (errorCode.EQ.0) THEN
            errorCode = 2
            errorWaterId = aidbt2(nw)
            errorHeavyId = aidbt1(nearestHeavy)
         END IF
!$OMP END CRITICAL(kirkwood_error)
         CYCLE
      END IF

      ! Minimum-image O-H vectors keep water geometry continuous across PBC.
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

      rOH1Sq = rxOH1*rxOH1 + ryOH1*ryOH1 + rzOH1*rzOH1
      rOH2Sq = rxOH2*rxOH2 + ryOH2*ryOH2 + rzOH2*rzOH2
      IF (rOH1Sq.LE.tinyLength*tinyLength .OR. &
          rOH2Sq.LE.tinyLength*tinyLength) THEN
!$OMP CRITICAL(kirkwood_error)
         IF (errorCode.EQ.0) THEN
            errorCode = 3
            errorWaterId = aidbt2(nw)
         END IF
!$OMP END CRITICAL(kirkwood_error)
         CYCLE
      END IF

      ! For the supported neutral three-site models, qO+2qH=0 and the dipole
      ! direction is proportional to the sum of the two O-H vectors.
      rxiDip = qhw*(rxOH1 + rxOH2)
      ryiDip = qhw*(ryOH1 + ryOH2)
      rziDip = qhw*(rzOH1 + rzOH2)
      magDip = SQRT(rxiDip*rxiDip + ryiDip*ryiDip + rziDip*rziDip)

      IF (magDip.LE.tinyLength) THEN
!$OMP CRITICAL(kirkwood_error)
         IF (errorCode.EQ.0) THEN
            errorCode = 4
            errorWaterId = aidbt2(nw)
         END IF
!$OMP END CRITICAL(kirkwood_error)
         CYCLE
      END IF

      ux = rxiDip/magDip
      uy = ryiDip/magDip
      uz = rziDip/magDip

      nx = bestDifX/rDif
      ny = bestDifY/rDif
      nz = bestDifZ/rDif

      p = ux*nx + uy*ny + uz*nz
      p = MAX(-1.0d0,MIN(1.0d0,p))

      uparx = ux - nx*p
      upary = uy - ny*p
      uparz = uz - nz*p

      DO region=1,nKirkwoodRegions
         IF (.NOT.selectedRegion(region)) CYCLE

         regionCountThread(region,tid) = regionCountThread(region,tid) + 1

         regionSumUThread(1,region,tid) = regionSumUThread(1,region,tid) + ux
         regionSumUThread(2,region,tid) = regionSumUThread(2,region,tid) + uy
         regionSumUThread(3,region,tid) = regionSumUThread(3,region,tid) + uz

         regionSumPThread(region,tid) = regionSumPThread(region,tid) + p

         regionSumUParThread(1,region,tid) = &
            regionSumUParThread(1,region,tid) + uparx
         regionSumUParThread(2,region,tid) = &
            regionSumUParThread(2,region,tid) + upary
         regionSumUParThread(3,region,tid) = &
            regionSumUParThread(3,region,tid) + uparz

         regionSumNThread(1,region,tid) = regionSumNThread(1,region,tid) + nx
         regionSumNThread(2,region,tid) = regionSumNThread(2,region,tid) + ny
         regionSumNThread(3,region,tid) = regionSumNThread(3,region,tid) + nz
      END DO
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
         WRITE(*,*) 'VALIDATED 10 A CELL-SEARCH RADIUS.'
         WRITE(*,*) 'WATER OXYGEN PSF ID: ', errorWaterId
         IF (errorValue1.GT.0.0d0) THEN
            WRITE(*,*) 'MINIMUM ANY-PROTEIN DISTANCE: ', errorValue1
         END IF
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

   ! O(nthreads) reduction for each of the three scalar regions.
   DO region=1,nKirkwoodRegions
      DO tid=1,actualThreads
         regionFrameN(region) = regionFrameN(region) + &
                                regionCountThread(region,tid)
         regionFrameSumU(:,region) = regionFrameSumU(:,region) + &
                                     regionSumUThread(:,region,tid)
         regionFrameSumP(region) = regionFrameSumP(region) + &
                                   regionSumPThread(region,tid)
         regionFrameSumUPar(:,region) = regionFrameSumUPar(:,region) + &
                                        regionSumUParThread(:,region,tid)
         regionFrameSumN(:,region) = regionFrameSumN(:,region) + &
                                     regionSumNThread(:,region,tid)
      END DO

      CALL computeKirkwoodObservables( &
           regionFrameN(region), regionFrameSumU(:,region), &
           regionFrameSumP(region), regionFrameSumUPar(:,region), &
           regionFrameSumN(:,region), regionFrameValid(region), &
           regionFrameNormalCoherence(region), regionFrameGK(region), &
           regionFrameGPerp(region), regionFrameGParallel(region), &
           regionFrameResidual(region))

      IF (regionFrameValid(region)) THEN
         regionValidFrames(region) = regionValidFrames(region) + 1
         regionNAcc(region) = regionNAcc(region) + DBLE(regionFrameN(region))
         regionNormalCoherenceAcc(region) = regionNormalCoherenceAcc(region) + &
                                            regionFrameNormalCoherence(region)
         regionGKAcc(region) = regionGKAcc(region) + regionFrameGK(region)
         regionGPerpAcc(region) = regionGPerpAcc(region) + &
                                  regionFrameGPerp(region)
         regionGParallelAcc(region) = regionGParallelAcc(region) + &
                                      regionFrameGParallel(region)
         regionResidualAcc(region) = regionResidualAcc(region) + &
                                     regionFrameResidual(region)
      END IF
   END DO
RETURN
END SUBROUTINE kirkwoodMain
!-----------------------------------------------------------------------------------------
SUBROUTINE computeKirkwoodObservables(n,sumU,sumP,sumUPar,sumN,valid, &
                                      normalCoherence,gK,gPerp,gParallel,residual)
IMPLICIT NONE
INTEGER, INTENT(IN) :: n
DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: sumU, sumUPar, sumN
DOUBLE PRECISION, INTENT(IN) :: sumP
LOGICAL, INTENT(OUT) :: valid
DOUBLE PRECISION, INTENT(OUT) :: normalCoherence, gK, gPerp, gParallel, residual
DOUBLE PRECISION :: invN

   valid = (n.GT.0)
   normalCoherence = 0.0d0
   gK = 0.0d0
   gPerp = 0.0d0
   gParallel = 0.0d0
   residual = 0.0d0

   IF (.NOT.valid) RETURN

   invN = 1.0d0/DBLE(n)

   gK = DOT_PRODUCT(sumU,sumU)*invN
   gPerp = sumP*sumP*invN
   gParallel = DOT_PRODUCT(sumUPar,sumUPar)*invN
   residual = gK - gPerp - gParallel
   normalCoherence = SQRT(DOT_PRODUCT(sumN,sumN))*invN
RETURN
END SUBROUTINE computeKirkwoodObservables
!-----------------------------------------------------------------------------------------
SUBROUTINE kirkwoodSanityCheck
IMPLICIT NONE
INTEGER :: i, j, k, n
DOUBLE PRECISION, PARAMETER :: tol = 5.0d-13
DOUBLE PRECISION :: u(3,3), nvec(3,3), upar(3,3), p(3)
DOUBLE PRECISION :: sumU(3), sumUPar(3), sumN(3), sumP
DOUBLE PRECISION :: gK, gPerp, gParallel, residual, coherence
DOUBLE PRECISION :: gKExp, gPerpExp, gParallelExp
DOUBLE PRECISION :: gKMean, gPerpMean, gParallelMean, residualMean
DOUBLE PRECISION :: basis(3,3), oneUPar(3), oneN(3), oneP
LOGICAL :: valid

   WRITE(*,*) '----------------------------------------------------------------------'
   WRITE(*,*) 'RUNNING INTERNAL KIRKWOOD SANITY CHECKS ...'

   ! Check 1: O(N) reduced sums must reproduce the explicit ij definitions.
   n = 3
   u(:,1) = [1.0d0,0.0d0,0.0d0]
   u(:,2) = [0.0d0,1.0d0,0.0d0]
   u(:,3) = [SQRT(0.5d0),0.0d0,SQRT(0.5d0)]
   nvec(:,1) = [0.0d0,0.0d0,1.0d0]
   nvec(:,2) = [1.0d0,0.0d0,0.0d0]
   nvec(:,3) = [0.0d0,1.0d0,0.0d0]
   DO i=1,n
      p(i) = DOT_PRODUCT(u(:,i),nvec(:,i))
      upar(:,i) = u(:,i) - nvec(:,i)*p(i)
   END DO

   sumU = 0.0d0
   sumUPar = 0.0d0
   sumN = 0.0d0
   sumP = 0.0d0
   DO i=1,n
      sumU = sumU + u(:,i)
      sumUPar = sumUPar + upar(:,i)
      sumN = sumN + nvec(:,i)
      sumP = sumP + p(i)
   END DO

   CALL computeKirkwoodObservables(n,sumU,sumP,sumUPar,sumN,valid, &
        coherence,gK,gPerp,gParallel,residual)

   gKExp = 0.0d0
   gPerpExp = 0.0d0
   gParallelExp = 0.0d0
   DO i=1,n
      DO j=1,n
         gKExp = gKExp + DOT_PRODUCT(u(:,i),u(:,j))/DBLE(n)
         gPerpExp = gPerpExp + p(i)*p(j)/DBLE(n)
         gParallelExp = gParallelExp + &
                        DOT_PRODUCT(upar(:,i),upar(:,j))/DBLE(n)
      END DO
   END DO

   IF (.NOT.valid .OR. ABS(gK-gKExp).GT.tol .OR. &
       ABS(gPerp-gPerpExp).GT.tol .OR. &
       ABS(gParallel-gParallelExp).GT.tol) THEN
      WRITE(*,*) 'ERROR: REDUCED KIRKWOOD SUMS DO NOT MATCH EXPLICIT IJ SUMS.'
      WRITE(*,*) 'gK reduced/explicit:       ', gK, gKExp
      WRITE(*,*) 'gPerp reduced/explicit:    ', gPerp, gPerpExp
      WRITE(*,*) 'gParallel reduced/explicit:', gParallel, gParallelExp
      STOP
   END IF
   WRITE(*,*) 'SANITY CHECK 1 PASS: REDUCED SUMS = EXPLICIT IJ SUMS'
   WRITE(*,*) '                     FOR WATER-SPECIFIC LOCAL NORMALS.'

   ! Check 2: common-normal decomposition and isotropic 1/3 : 2/3 limit.
   ! Three equally weighted synthetic N=1 frames use u=x,y,z and n=z.
   basis(:,1) = [1.0d0,0.0d0,0.0d0]
   basis(:,2) = [0.0d0,1.0d0,0.0d0]
   basis(:,3) = [0.0d0,0.0d0,1.0d0]
   oneN = [0.0d0,0.0d0,1.0d0]

   gKMean = 0.0d0
   gPerpMean = 0.0d0
   gParallelMean = 0.0d0
   residualMean = 0.0d0

   DO k=1,3
      oneP = DOT_PRODUCT(basis(:,k),oneN)
      oneUPar = basis(:,k) - oneN*oneP
      CALL computeKirkwoodObservables(1,basis(:,k),oneP,oneUPar,oneN, &
           valid,coherence,gK,gPerp,gParallel,residual)
      gKMean = gKMean + gK/3.0d0
      gPerpMean = gPerpMean + gPerp/3.0d0
      gParallelMean = gParallelMean + gParallel/3.0d0
      residualMean = residualMean + residual/3.0d0
   END DO

   IF (ABS(gKMean-1.0d0).GT.tol .OR. &
       ABS(gPerpMean-1.0d0/3.0d0).GT.tol .OR. &
       ABS(gParallelMean-2.0d0/3.0d0).GT.tol .OR. &
       ABS(residualMean).GT.tol) THEN
      WRITE(*,*) 'ERROR: COMMON-NORMAL / ISOTROPIC SANITY CHECK FAILED.'
      WRITE(*,*) 'gK, gPerp, gParallel, residual: ', &
                 gKMean, gPerpMean, gParallelMean, residualMean
      STOP
   END IF
   WRITE(*,*) 'SANITY CHECK 2 PASS: COMMON NORMAL gPerp + gParallel = gK.'
   WRITE(*,*) 'SANITY CHECK 3 PASS: ISOTROPIC LIMIT gPerp/gK = 1/3,'
   WRITE(*,*) '                     gParallel/gK = 2/3.'
   WRITE(*,*) 'NOTE: THESE IDENTITIES ARE NOT REQUIRED FOR A CURVED INTERFACE'
   WRITE(*,*) '      WHEN DIFFERENT WATERS HAVE DIFFERENT LOCAL NORMALS.'
   WRITE(*,*) '----------------------------------------------------------------------'
RETURN
END SUBROUTINE kirkwoodSanityCheck
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
SUBROUTINE popBetaStructArrays(num)
USE atoms
USE beta1
USE beta2
IMPLICIT NONE
INTEGER :: i, num, nbt0, nbt1, nbt2, nbt

   nbt0 = 0
   nbt1 = 0
   nbt2 = 0

   DO i=1,num
      IF (INT(beta(i)).EQ.0.AND.nbeta0.GT.0) THEN
         nbt0 = nbt0 + 1
      ELSE IF (INT(beta(i)).EQ.1.AND.nbeta1.GT.0) THEN
         nbt1 = nbt1 + 1
         aidbt1(nbt1) = aindex(i)
         anamebt1(nbt1) = atomname(i)
         atypbt1(nbt1) = atomtype(i)
         rnamebt1(nbt1) = resname(i)
         ridbt1(nbt1) = resid(i)
         chrgbt1(nbt1) = charge(i)
         massbt1(nbt1) = mass(i)
         segbt1(nbt1) = segname(i)
         atmnm2bt1(nbt1) = atomname2(i)
      ELSE IF (INT(beta(i)).EQ.2.AND.nbeta2.GT.0) THEN
         nbt2 = nbt2 + 1
         aidbt2(nbt2) = aindex(i)
         anamebt2(nbt2) = atomname(i)
         atypbt2(nbt2) = atomtype(i)
         rnamebt2(nbt2) = resname(i)
         ridbt2(nbt2) = resid(i)
         chrgbt2(nbt2) = charge(i)
         massbt2(nbt2) = mass(i)
         segbt2(nbt2) = segname(i)
         atmnm2bt2(nbt2) = atomname2(i)
      ELSE
         WRITE(*,*) 'THERE ARE BETA VALUES IN YOUR PDB FILE OUTSIDE 0,1,2.'
         STOP
      END IF
   END DO

   nbt = nbt0 + nbt1 + nbt2
   IF (nbt.NE.natom) THEN
      WRITE(*,*) 'ERROR: BETA COUNTS DO NOT SUM TO THE PSF ATOM COUNT.'
      STOP
   END IF
RETURN
END SUBROUTINE popBetaStructArrays
!-----------------------------------------------------------------------------------------
SUBROUTINE popBetaArrays(num)
USE atoms
USE coords
USE beta1
USE beta2
IMPLICIT NONE
INTEGER :: i, num, nbt0, nbt1, nbt2, nbt

   nbt0 = 0
   nbt1 = 0
   nbt2 = 0

   DO i=1,num
      IF (INT(beta(i)).EQ.0.AND.nbeta0.GT.0) THEN
         nbt0 = nbt0 + 1
      ELSE IF (INT(beta(i)).EQ.1.AND.nbeta1.GT.0) THEN
         nbt1 = nbt1 + 1
         rxbt1(nbt1) = rx(i)
         rybt1(nbt1) = ry(i)
         rzbt1(nbt1) = rz(i)
      ELSE IF (INT(beta(i)).EQ.2.AND.nbeta2.GT.0) THEN
         nbt2 = nbt2 + 1
         rxbt2(nbt2) = rx(i)
         rybt2(nbt2) = ry(i)
         rzbt2(nbt2) = rz(i)
      ELSE
         WRITE(*,*) 'THERE ARE BETA VALUES IN YOUR PDB FILE OUTSIDE 0,1,2.'
         STOP
      END IF
   END DO

   nbt = nbt0 + nbt1 + nbt2
   IF (nbt.NE.natom) THEN
      WRITE(*,*) 'ERROR: FRAME BETA COUNTS DO NOT SUM TO THE PSF ATOM COUNT.'
      STOP
   END IF
RETURN
END SUBROUTINE popBetaArrays
!-----------------------------------------------------------------------------------------
SUBROUTINE buildKirkwoodTopology
USE atoms, ONLY : natom, nbeta1
USE beta1
USE topology
USE surface_geometry
USE kirkwood
USE omp_lib, ONLY : omp_get_max_threads
IMPLICIT NONE
INTEGER :: j, b, atomId1, atomId2, betaIndex1, betaIndex2
INTEGER :: nProteinHydrogen, nMappedHydrogen
INTEGER, ALLOCATABLE, DIMENSION(:) :: globalToBeta1, parentCount

   IF (.NOT.ALLOCATED(bondAtom1) .OR. .NOT.ALLOCATED(bondAtom2)) THEN
      WRITE(*,*) 'ERROR: PSF BOND ARRAYS WERE NOT INITIALIZED.'
      STOP
   END IF

   ALLOCATE(isProteinHeavy(1:nbeta1))
   ALLOCATE(parentHeavyIndex(1:nbeta1))
   ALLOCATE(parentCount(1:nbeta1))
   ALLOCATE(globalToBeta1(1:natom))
   ALLOCATE(proteinCellNext(1:nbeta1))

   isProteinHeavy(:) = .FALSE.
   parentHeavyIndex(:) = 0
   parentCount(:) = 0
   globalToBeta1(:) = 0
   proteinCellNext(:) = 0

   nProteinHeavy = 0
   nProteinHydrogen = 0

   DO j=1,nbeta1
      IF (massbt1(j).GT.hydrogenMassMax) THEN
         isProteinHeavy(j) = .TRUE.
         parentHeavyIndex(j) = j
         nProteinHeavy = nProteinHeavy + 1
      ELSE
         nProteinHydrogen = nProteinHydrogen + 1
      END IF

      atomId1 = aidbt1(j)
      IF (atomId1.LT.1 .OR. atomId1.GT.natom) THEN
         WRITE(*,*) 'ERROR: BETA=1 PSF ATOM ID IS OUT OF RANGE: ', atomId1
         STOP
      END IF
      IF (globalToBeta1(atomId1).NE.0) THEN
         WRITE(*,*) 'ERROR: DUPLICATE BETA=1 PSF ATOM ID: ', atomId1
         STOP
      END IF
      globalToBeta1(atomId1) = j
   END DO

   pro16OIndex = 0
   DO j=1,nbeta1
      IF (TRIM(ridbt1(j)).EQ.'16' .AND. TRIM(anamebt1(j)).EQ.'O') THEN
         IF (pro16OIndex.NE.0) THEN
            WRITE(*,*) 'ERROR: MULTIPLE BETA=1 PRO16 O ATOMS FOUND.'
            STOP
         END IF
         pro16OIndex = j
      END IF
   END DO

   IF (pro16OIndex.EQ.0) THEN
      WRITE(*,*) 'ERROR: PRO16 CARBONYL O NOT FOUND IN BETA=1 PROTEIN.'
      STOP
   END IF

   ! Map protein H atoms to their heavy-atom parents.  This is used only to
   ! validate the geometric bound underlying the 7 A any-protein selection.
   DO b=1,nBonds
      atomId1 = bondAtom1(b)
      atomId2 = bondAtom2(b)

      IF (atomId1.LT.1 .OR. atomId1.GT.natom) CYCLE
      IF (atomId2.LT.1 .OR. atomId2.GT.natom) CYCLE

      betaIndex1 = globalToBeta1(atomId1)
      betaIndex2 = globalToBeta1(atomId2)
      IF (betaIndex1.EQ.0 .OR. betaIndex2.EQ.0) CYCLE

      IF (.NOT.isProteinHeavy(betaIndex1) .AND. isProteinHeavy(betaIndex2)) THEN
         parentHeavyIndex(betaIndex1) = betaIndex2
         parentCount(betaIndex1) = parentCount(betaIndex1) + 1
      ELSE IF (.NOT.isProteinHeavy(betaIndex2) .AND. isProteinHeavy(betaIndex1)) THEN
         parentHeavyIndex(betaIndex2) = betaIndex1
         parentCount(betaIndex2) = parentCount(betaIndex2) + 1
      END IF
   END DO

   nMappedHydrogen = 0
   DO j=1,nbeta1
      IF (.NOT.isProteinHeavy(j)) THEN
         IF (parentCount(j).NE.1 .OR. parentHeavyIndex(j).EQ.0) THEN
            WRITE(*,*) 'ERROR: PROTEIN HYDROGEN DOES NOT HAVE EXACTLY ONE'
            WRITE(*,*) 'BETA=1 HEAVY-ATOM PARENT IN THE PSF TOPOLOGY.'
            WRITE(*,*) 'ATOM ID / RESID / NAME: ', aidbt1(j), &
                       TRIM(ridbt1(j)), TRIM(anamebt1(j))
            WRITE(*,*) 'NUMBER OF HEAVY PARENTS FOUND: ', parentCount(j)
            STOP
         END IF
         nMappedHydrogen = nMappedHydrogen + 1
      END IF
   END DO

   nOpenMPThreads = MAX(1,omp_get_max_threads())
   ALLOCATE(regionCountThread(1:nKirkwoodRegions,1:nOpenMPThreads))
   ALLOCATE(regionSumUThread(1:3,1:nKirkwoodRegions,1:nOpenMPThreads))
   ALLOCATE(regionSumPThread(1:nKirkwoodRegions,1:nOpenMPThreads))
   ALLOCATE(regionSumUParThread(1:3,1:nKirkwoodRegions,1:nOpenMPThreads))
   ALLOCATE(regionSumNThread(1:3,1:nKirkwoodRegions,1:nOpenMPThreads))

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

   regionValidFrames(:) = 0
   regionNAcc(:) = 0.0d0
   regionNormalCoherenceAcc(:) = 0.0d0
   regionGKAcc(:) = 0.0d0
   regionGPerpAcc(:) = 0.0d0
   regionGParallelAcc(:) = 0.0d0
   regionResidualAcc(:) = 0.0d0

   WRITE(*,*) '----------------------------------------------------------------------'
   WRITE(*,*) 'KIRKWOOD GEOMETRY INITIALIZED'
   WRITE(*,*) 'PROTEIN ATOMS (BETA=1):         ', nbeta1
   WRITE(*,*) 'PROTEIN HEAVY ATOMS:            ', nProteinHeavy
   WRITE(*,*) 'PROTEIN HYDROGENS:              ', nProteinHydrogen
   WRITE(*,*) 'HYDROGENS MAPPED TO A PARENT:   ', nMappedHydrogen
   WRITE(*,*) 'INTERFACE CUTOFF (ANGSTROM):     ', interfaceCutoff
   WRITE(*,*) 'PRO16 CUTOFFS (ANGSTROM):        ', pro16Cutoff7, pro16Cutoff10
   WRITE(*,*) 'HEAVY SEARCH RADIUS (ANGSTROM):  ', proteinSearchRadius
   WRITE(*,*) 'OPENMP MAXIMUM THREADS:           ', nOpenMPThreads
   WRITE(*,*) 'PRO16 CARBONYL O PSF ID:          ', aidbt1(pro16OIndex)
   WRITE(*,*) '----------------------------------------------------------------------'

   DEALLOCATE(globalToBeta1,parentCount)
RETURN
END SUBROUTINE buildKirkwoodTopology
!-----------------------------------------------------------------------------------------
SUBROUTINE deKirkwoodArrays
USE topology
USE surface_geometry
USE kirkwood
IMPLICIT NONE

   IF (ALLOCATED(bondAtom1)) DEALLOCATE(bondAtom1)
   IF (ALLOCATED(bondAtom2)) DEALLOCATE(bondAtom2)
   IF (ALLOCATED(isProteinHeavy)) DEALLOCATE(isProteinHeavy)
   IF (ALLOCATED(parentHeavyIndex)) DEALLOCATE(parentHeavyIndex)
   IF (ALLOCATED(proteinCellHead)) DEALLOCATE(proteinCellHead)
   IF (ALLOCATED(proteinCellNext)) DEALLOCATE(proteinCellNext)

   IF (ALLOCATED(regionCountThread)) DEALLOCATE(regionCountThread)
   IF (ALLOCATED(regionSumUThread)) DEALLOCATE(regionSumUThread)
   IF (ALLOCATED(regionSumPThread)) DEALLOCATE(regionSumPThread)
   IF (ALLOCATED(regionSumUParThread)) DEALLOCATE(regionSumUParThread)
   IF (ALLOCATED(regionSumNThread)) DEALLOCATE(regionSumNThread)

   nBonds = 0
   nProteinHeavy = 0
   proteinCellNx = 0
   proteinCellNy = 0
   proteinCellNz = 0
   proteinCellCount = 0
   nOpenMPThreads = 1
   pro16OIndex = 0
RETURN
END SUBROUTINE deKirkwoodArrays
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

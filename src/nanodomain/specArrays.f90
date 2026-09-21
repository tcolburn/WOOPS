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
         WRITE(*,*) 'THERE ARE BETA VALUES IN YOUR PDB FILE'
         WRITE(*,*) 'THAT ARE OUT OF THE RANGE FOR THIS ANALYSIS.'
         WRITE(*,*) 'THE ALLOWED VALUES ARE 0, 1, AND 2.'
         WRITE(*,*) 'THE PROGRAM WILL NOW EXIT.'
         STOP
      END IF
   END DO

   nbt = nbt0 + nbt1 + nbt2

   IF (nbt.NE.natom) THEN
      WRITE(*,*) 'THE TOTAL NUMBER OF ATOMS CALCULATED IN popBetaStructArrays'
      WRITE(*,*) 'DOES NOT EQUAL THE NUMBER OF ATOMS READ FROM THE PSF FILE.'
      WRITE(*,*) 'THE PROGRAM WILL NOW EXIT.'
      STOP
   END IF
RETURN
END SUBROUTINE popBetaStructArrays
!-----------------------------------------------------------------------------------------
SUBROUTINE popBetaArrays(num)
USE input
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
         WRITE(*,*) 'THERE ARE BETA VALUES IN YOUR PDB FILE'
         WRITE(*,*) 'THAT ARE OUT OF THE RANGE FOR THIS ANALYSIS.'
         WRITE(*,*) 'THE ALLOWED VALUES ARE 0, 1, AND 2.'
         WRITE(*,*) 'THE PROGRAM WILL NOW EXIT.'
         STOP
      END IF
   END DO

   nbt = nbt0 + nbt1 + nbt2

   IF (nbt.NE.natom) THEN
      WRITE(*,*) 'THE TOTAL NUMBER OF ATOMS CALCULATED IN popBetaArrays'
      WRITE(*,*) 'DOES NOT EQUAL THE ORIGINAL NUMBER OF ATOMS READ FROM THE PSF FILE.'
      WRITE(*,*) 'THE PROGRAM WILL NOW EXIT.'
      STOP
   END IF
RETURN
END SUBROUTINE popBetaArrays
!-----------------------------------------------------------------------------------------
SUBROUTINE buildProteinTopology
USE atoms, ONLY : natom, nbeta1
USE beta1
USE topology
USE nanodomain
USE p1p2, ONLY : nPAtom
USE omp_lib, ONLY : omp_get_max_threads
IMPLICIT NONE
INTEGER :: j, h, b, atomId1, atomId2, betaIndex1, betaIndex2
INTEGER :: nProteinHydrogen, nMappedHydrogen
INTEGER, ALLOCATABLE, DIMENSION(:) :: globalToBeta1, parentCount

   IF (nbeta1.NE.nPAtom) THEN
      WRITE(*,*) 'ERROR: INPUT NUMBER OF PROTEIN ATOMS DOES NOT MATCH BETA=1 COUNT.'
      WRITE(*,*) 'INPUT nPAtom: ', nPAtom
      WRITE(*,*) 'BETA=1 COUNT: ', nbeta1
      STOP
   END IF

   IF (.NOT.ALLOCATED(bondAtom1) .OR. .NOT.ALLOCATED(bondAtom2)) THEN
      WRITE(*,*) 'ERROR: PSF BOND ARRAYS WERE NOT INITIALIZED.'
      STOP
   END IF

   ALLOCATE(isProteinHeavy(1:nbeta1))
   ALLOCATE(parentHeavyIndex(1:nbeta1))
   ALLOCATE(parentCount(1:nbeta1))
   ALLOCATE(globalToBeta1(1:natom))

   isProteinHeavy(:) = .FALSE.
   parentHeavyIndex(:) = 0
   parentCount(:) = 0
   globalToBeta1(:) = 0

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

   ALLOCATE(heavyBetaIndex(1:nProteinHeavy))
   h = 0
   DO j=1,nbeta1
      IF (isProteinHeavy(j)) THEN
         h = h + 1
         heavyBetaIndex(h) = j
      END IF
   END DO

   ! Read hydrogen-parent relationships directly from the PSF bond topology.
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
            WRITE(*,*) 'ATOM ID / RESID / NAME: ', aidbt1(j), TRIM(ridbt1(j)), TRIM(anamebt1(j))
            WRITE(*,*) 'NUMBER OF HEAVY PARENTS FOUND: ', parentCount(j)
            STOP
         END IF
         nMappedHydrogen = nMappedHydrogen + 1
      END IF
   END DO

   ALLOCATE(waterCountFrame(1:nbeta1))
   ALLOCATE(p1FrameSum(1:nbeta1))
   ALLOCATE(p2FrameSum(1:nbeta1))
   ALLOCATE(p21FrameSum(1:nbeta1))
   ALLOCATE(proteinCellNext(1:nbeta1))

   nOpenMPThreads = MAX(1,omp_get_max_threads())
   ALLOCATE(waterCountThread(1:nbeta1,1:nOpenMPThreads))
   ALLOCATE(p1ThreadSum(1:nbeta1,1:nOpenMPThreads))
   ALLOCATE(p2ThreadSum(1:nbeta1,1:nOpenMPThreads))
   ALLOCATE(p21ThreadSum(1:nbeta1,1:nOpenMPThreads))
   ALLOCATE(touchedCountThread(1:nOpenMPThreads))
   ALLOCATE(touchedHeavyThread(1:nProteinHeavy,1:nOpenMPThreads))

   waterCountFrame(:) = 0
   p1FrameSum(:) = 0.0d0
   p2FrameSum(:) = 0.0d0
   p21FrameSum(:) = 0.0d0
   waterCountThread(:,:) = 0
   p1ThreadSum(:,:) = 0.0d0
   p2ThreadSum(:,:) = 0.0d0
   p21ThreadSum(:,:) = 0.0d0
   touchedCountThread(:) = 0
   touchedHeavyThread(:,:) = 0

   WRITE(*,*) '----------------------------------------------------------------------'
   WRITE(*,*) 'NANODOMAIN TOPOLOGY INITIALIZED'
   WRITE(*,*) 'PROTEIN ATOMS (BETA=1):       ', nbeta1
   WRITE(*,*) 'PROTEIN HEAVY ATOMS:          ', nProteinHeavy
   WRITE(*,*) 'PROTEIN HYDROGENS:            ', nProteinHydrogen
   WRITE(*,*) 'HYDROGENS MAPPED TO A PARENT: ', nMappedHydrogen
   WRITE(*,*) 'HYDRATION CUTOFF (ANGSTROM):  ', hydrationCutoff
   WRITE(*,*) 'HEAVY SEARCH RADIUS (ANGSTROM): ', assignmentSearchRadius
   WRITE(*,*) 'OPENMP MAXIMUM THREADS:         ', nOpenMPThreads
   WRITE(*,*) '----------------------------------------------------------------------'

   DEALLOCATE(globalToBeta1,parentCount)
RETURN
END SUBROUTINE buildProteinTopology
!-----------------------------------------------------------------------------------------
SUBROUTINE deNanodomainArrays
USE topology
USE nanodomain
IMPLICIT NONE

   IF (ALLOCATED(bondAtom1)) DEALLOCATE(bondAtom1)
   IF (ALLOCATED(bondAtom2)) DEALLOCATE(bondAtom2)

   IF (ALLOCATED(isProteinHeavy)) DEALLOCATE(isProteinHeavy)
   IF (ALLOCATED(heavyBetaIndex)) DEALLOCATE(heavyBetaIndex)
   IF (ALLOCATED(parentHeavyIndex)) DEALLOCATE(parentHeavyIndex)
   IF (ALLOCATED(waterCountFrame)) DEALLOCATE(waterCountFrame)
   IF (ALLOCATED(p1FrameSum)) DEALLOCATE(p1FrameSum)
   IF (ALLOCATED(p2FrameSum)) DEALLOCATE(p2FrameSum)
   IF (ALLOCATED(p21FrameSum)) DEALLOCATE(p21FrameSum)
   IF (ALLOCATED(proteinCellHead)) DEALLOCATE(proteinCellHead)
   IF (ALLOCATED(proteinCellNext)) DEALLOCATE(proteinCellNext)
   IF (ALLOCATED(waterCountThread)) DEALLOCATE(waterCountThread)
   IF (ALLOCATED(p1ThreadSum)) DEALLOCATE(p1ThreadSum)
   IF (ALLOCATED(p2ThreadSum)) DEALLOCATE(p2ThreadSum)
   IF (ALLOCATED(p21ThreadSum)) DEALLOCATE(p21ThreadSum)
   IF (ALLOCATED(touchedCountThread)) DEALLOCATE(touchedCountThread)
   IF (ALLOCATED(touchedHeavyThread)) DEALLOCATE(touchedHeavyThread)

   nBonds = 0
   nProteinHeavy = 0
   proteinCellNx = 0
   proteinCellNy = 0
   proteinCellNz = 0
   proteinCellCount = 0
   nOpenMPThreads = 1
RETURN
END SUBROUTINE deNanodomainArrays
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
SUBROUTINE assignWatType(wtype)
USE water
IMPLICIT NONE
INTEGER :: wtype

   WRITE(*,*) '----------------------------------------------------------------------'

   IF (wattype.EQ.0) THEN
      watertype = 'TIP3'
      oxname = 'OH2'
      hy1name = 'H1'
      hy2name = 'H2'
      omname = 'XX'
      qow = -0.834d0
      qhw = 0.417d0
      WRITE(*,*) 'USING TIP3P WATER TYPE ...'
   ELSE IF (wattype.EQ.1) THEN
      WRITE(*,*) 'TIP4P WATER TYPE IS NOT YET SUPPORTED FOR THIS PROGRAM!'
      WRITE(*,*) 'EXITING PROGRAM'
      STOP
   ELSE IF (wattype.EQ.2) THEN
      watertype = 'SPCE'
      oxname = 'OH2'
      hy1name = 'H1'
      hy2name = 'H2'
      omname = 'XX'
      qow = -0.8476d0
      qhw = 0.4238d0
      WRITE(*,*) 'USING SPCE WATER TYPE ...'
   ELSE
      WRITE(*,*) 'UNRECOGNIZED WATER TYPE:'
      WRITE(*,*) 'TIP3P -> wattype = 0 in input file'
      WRITE(*,*) 'TIP4P -> wattype = 1 in input file - NOT SUPPORTED'
      WRITE(*,*) 'SPCE  -> wattype = 2 in input file'
      STOP
   END IF

   WRITE(*,*) '----------------------------------------------------------------------'
RETURN
END SUBROUTINE assignWatType
!-------------------------------------------------------------------------------
SUBROUTINE setup
USE files
USE water
USE atoms
USE input
USE framesNum
IMPLICIT NONE
INTEGER :: nbetaTot

   CALL assignWatType(wattype)

   WRITE(*,*) '----------------------------------------------------------------------'
   WRITE(*,*) 'GETTING THE NUMBER OF ATOMS FROM THE PSF FILE ...'
   CALL getnum(psffile,psfInt,istart,natom)
   WRITE(*,*) 'THE NUMBER OF ATOMS REFERENCED IN THE PSF FILE IS ', natom
   WRITE(*,*) '----------------------------------------------------------------------'

   WRITE(*,*) 'ALLOCATING APPROPRIATE STRUCTURAL ARRAYS ...'
   CALL alStructArrays
   WRITE(*,*) 'DONE ALLOCATING APPROPRIATE STRUCTURAL ARRAYS ...'
   WRITE(*,*) '----------------------------------------------------------------------'

   WRITE(*,*) 'GETTING THE NUMBER OF FRAMES FROM THE DCD FILE ...'
   CALL getnframes(dcdfile,dcdInt,nframes)
   WRITE(*,*) 'THE NUMBER OF FRAMES REFERENCED IN THE DCD FILE IS ', nframes
   WRITE(*,*) '----------------------------------------------------------------------'

   WRITE(*,*) 'ALLOCATING APPROPRIATE DYNAMIC ARRAYS ...'
   CALL alCoordArrays(ntrajFrame)
   WRITE(*,*) 'DONE ALLOCATING APPROPRIATE DYNAMIC ARRAYS ...'
   WRITE(*,*) '----------------------------------------------------------------------'

   CALL psfread(psffile,psfInt,istart,natom)
   CALL pdbread(pdbfile,pdbInt,natom)
   CALL configinfo(natom,wattype,nwater,nbeta0,nbeta1,nbeta2,nnwater)

   WRITE(*,*) '     # WATERS  |  # BETA 0  |  # BETA 1  |  # BETA 2  |  # NON-WATERS '
   WRITE(*,*) '    ------------------------------------------------------------------'
   WRITE(*,*) nwater, nbeta0, nbeta1, nbeta2, nnwater
   WRITE(*,*) '    ------------------------------------------------------------------'

   IF (nbeta1.GT.0.AND.betaFlg.NE.0) CALL alBeta1Arrays(nbeta1)
   IF (nbeta2.GT.0.AND.betaFlg.NE.0) CALL alBeta2Arrays(nbeta2)

   nbetaTot = nbeta1 + nbeta2
   IF (nbetaTot.GT.0.AND.betaFlg.NE.0) CALL popBetaStructArrays(natom)

   IF (nbeta1.LE.0) THEN
      WRITE(*,*) 'ERROR: KIRKWOOD ANALYSIS REQUIRES BETA=1 PROTEIN ATOMS.'
      STOP
   END IF

   IF (nbeta2.LE.0 .OR. MOD(nbeta2,3).NE.0) THEN
      WRITE(*,*) 'ERROR: KIRKWOOD ANALYSIS REQUIRES BETA=2 THREE-SITE WATERS.'
      STOP
   END IF

   IF (nbeta1.NE.nPAtomInput) THEN
      WRITE(*,*) 'ERROR: INPUT PROTEIN-ATOM COUNT DOES NOT MATCH BETA=1 COUNT.'
      WRITE(*,*) 'INPUT COUNT / BETA=1 COUNT: ', nPAtomInput, nbeta1
      STOP
   END IF

   CALL buildKirkwoodTopology
   CALL kirkwoodSanityCheck
RETURN
END SUBROUTINE setup
!-------------------------------------------------------------------------------
SUBROUTINE exitProg
USE atoms, ONLY : nbeta1, nbeta2
IMPLICIT NONE

   CALL deKirkwoodArrays
   CALL deCoordArrays
   CALL deStructArrays

   IF (nbeta1.GT.0) CALL deBeta1Arrays
   IF (nbeta2.GT.0) CALL deBeta2Arrays
RETURN
END SUBROUTINE exitProg
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

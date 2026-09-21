!------------------------------------------------------------------------------
! Core modules for the dedicated scalar Kirkwood analysis.
!------------------------------------------------------------------------------
MODULE input
   IMPLICIT NONE
   INTEGER :: betaFlg
   INTEGER :: nresInput
   INTEGER :: nPAtomInput
END MODULE input
!------------------------------------------------------------------------------
MODULE framesNum
   IMPLICIT NONE
   INTEGER :: frameNum, ntrajFrame, startFrame, endFrame, ftrack
END MODULE framesNum
!------------------------------------------------------------------------------
MODULE files
   IMPLICIT NONE
   CHARACTER(LEN=60) :: psffile, pdbfile, dcdfile
   INTEGER, PARAMETER :: psfInt  = 100
   INTEGER, PARAMETER :: pdbInt  = 101
   INTEGER, PARAMETER :: dcdInt  = 102
   INTEGER :: istart
   INTEGER, PARAMETER :: noFileInt  = 201
   CHARACTER(LEN=60) :: noFileName
END MODULE files
!------------------------------------------------------------------------------
MODULE coords
   IMPLICIT NONE
   REAL*4, ALLOCATABLE, DIMENSION(:,:) :: rxf, ryf, rzf
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: rx, ry, rz
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: box
   DOUBLE PRECISION :: boxx, boxy, boxz, boxxSq, boxySq, boxzSq
END MODULE coords
!------------------------------------------------------------------------------
MODULE water
   IMPLICIT NONE
   INTEGER :: wattype
   CHARACTER(LEN=4) :: watertype
   CHARACTER(LEN=3) :: oxname
   CHARACTER(LEN=2) :: hy1name
   CHARACTER(LEN=2) :: hy2name
   CHARACTER(LEN=2) :: omname
   DOUBLE PRECISION :: qow, qhw, qmw
END MODULE water
!------------------------------------------------------------------------------
MODULE atoms
   IMPLICIT NONE
   INTEGER :: natom, nframes, nwater, nbeta0, nbeta1, nbeta2, nnwater
   INTEGER, ALLOCATABLE, DIMENSION(:) :: aindex
   CHARACTER(LEN=20), ALLOCATABLE, DIMENSION(:) :: resid, resname
   CHARACTER(LEN=20), ALLOCATABLE, DIMENSION(:) :: atomname, atomtype
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: charge, mass
   INTEGER, ALLOCATABLE, DIMENSION(:) :: pdbaindex, pdbresid
   CHARACTER(LEN=6), ALLOCATABLE, DIMENSION(:) :: pdbaindexC
   CHARACTER(LEN=20), ALLOCATABLE, DIMENSION(:) :: pdbresname
   CHARACTER(LEN=20), ALLOCATABLE, DIMENSION(:) :: pdbatomname
   CHARACTER(LEN=20), ALLOCATABLE, DIMENSION(:) :: segname
   CHARACTER(LEN=2), ALLOCATABLE, DIMENSION(:) :: atomname2
   REAL, ALLOCATABLE, DIMENSION(:) :: coor, beta
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: xini, yini, zini
END MODULE atoms
!------------------------------------------------------------------------------
MODULE beta1
   IMPLICIT NONE
   INTEGER, ALLOCATABLE, DIMENSION(:) :: aidbt1
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: rxbt1, rybt1, rzbt1
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: chrgbt1, massbt1
   CHARACTER(LEN=20), ALLOCATABLE, DIMENSION(:) :: rnamebt1, ridbt1
   CHARACTER(LEN=20), ALLOCATABLE, DIMENSION(:) :: anamebt1, atypbt1
   CHARACTER(LEN=20), ALLOCATABLE, DIMENSION(:) :: segbt1
   CHARACTER(LEN=2), ALLOCATABLE, DIMENSION(:) :: atmnm2bt1
END MODULE beta1
!------------------------------------------------------------------------------
MODULE beta2
   IMPLICIT NONE
   INTEGER, ALLOCATABLE, DIMENSION(:) :: aidbt2
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: rxbt2, rybt2, rzbt2
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: chrgbt2, massbt2
   CHARACTER(LEN=20), ALLOCATABLE, DIMENSION(:) :: rnamebt2, ridbt2
   CHARACTER(LEN=20), ALLOCATABLE, DIMENSION(:) :: anamebt2, atypbt2
   CHARACTER(LEN=20), ALLOCATABLE, DIMENSION(:) :: segbt2
   CHARACTER(LEN=2), ALLOCATABLE, DIMENSION(:) :: atmnm2bt2
END MODULE beta2
!------------------------------------------------------------------------------
MODULE constants
   IMPLICIT NONE
   DOUBLE PRECISION, PARAMETER :: mOx = 15.9994d0
   DOUBLE PRECISION, PARAMETER :: mHy = 1.0080d0
END MODULE constants
!------------------------------------------------------------------------------
MODULE topology
   IMPLICIT NONE
   INTEGER :: nBonds = 0
   INTEGER, ALLOCATABLE, DIMENSION(:) :: bondAtom1, bondAtom2
END MODULE topology
!------------------------------------------------------------------------------
MODULE surface_geometry
   IMPLICIT NONE

   ! Interface selection: water O within 7 A of ANY beta=1 protein atom,
   ! including protein hydrogens.
   DOUBLE PRECISION, PARAMETER :: interfaceCutoff = 7.0d0

   ! PRO16 selections are evaluated simultaneously in the same trajectory pass.
   DOUBLE PRECISION, PARAMETER :: pro16Cutoff7 = 7.0d0
   DOUBLE PRECISION, PARAMETER :: pro16Cutoff10 = 10.0d0

   DOUBLE PRECISION, PARAMETER :: hydrogenMassMax = 1.5d0
   DOUBLE PRECISION, PARAMETER :: maxProteinHydrogenBond = 2.5d0

   ! A 7 A interface water selected by a protein H has its H-parent heavy atom
   ! within 7 + 2.5 = 9.5 A.  A PRO16-10 A water always has PRO16 O itself
   ! within 10 A.  Therefore a 10 A heavy-atom search is sufficient for every
   ! water contributing to any of the three scalar selections.
   DOUBLE PRECISION, PARAMETER :: proteinSearchRadius = 10.0d0

   INTEGER :: nProteinHeavy = 0
   INTEGER :: proteinCellNx = 0
   INTEGER :: proteinCellNy = 0
   INTEGER :: proteinCellNz = 0
   INTEGER :: proteinCellCount = 0
   INTEGER, ALLOCATABLE, DIMENSION(:) :: proteinCellHead
   INTEGER, ALLOCATABLE, DIMENSION(:) :: proteinCellNext

   LOGICAL, ALLOCATABLE, DIMENSION(:) :: isProteinHeavy
   INTEGER, ALLOCATABLE, DIMENSION(:) :: parentHeavyIndex
END MODULE surface_geometry
!------------------------------------------------------------------------------
MODULE kirkwood
   IMPLICIT NONE

   ! Region numbering is fixed and shared by the analysis and output routines.
   INTEGER, PARAMETER :: nKirkwoodRegions = 3
   INTEGER, PARAMETER :: regionInterface = 1
   INTEGER, PARAMETER :: regionPro16_7 = 2
   INTEGER, PARAMETER :: regionPro16_10 = 3

   CHARACTER(LEN=32), PARAMETER, DIMENSION(nKirkwoodRegions) :: regionLabel = [ &
      CHARACTER(LEN=32) :: 'interface_7A', 'PRO16_7A', 'PRO16_10A' ]

   CHARACTER(LEN=48), PARAMETER, DIMENSION(nKirkwoodRegions) :: timeSeriesFile = [ &
      CHARACTER(LEN=48) :: 'kirkwood_interface.dat', &
                           'kirkwood_PRO16_7A.dat', &
                           'kirkwood_PRO16_10A.dat' ]

   CHARACTER(LEN=56), PARAMETER, DIMENSION(nKirkwoodRegions) :: summaryFile = [ &
      CHARACTER(LEN=56) :: 'kirkwood_interface_summary.dat', &
                           'kirkwood_PRO16_7A_summary.dat', &
                           'kirkwood_PRO16_10A_summary.dat' ]

   CHARACTER(LEN=64), PARAMETER, DIMENSION(nKirkwoodRegions) :: diagnosticFile = [ &
      CHARACTER(LEN=64) :: 'kirkwood_interface_diagnostics.dat', &
                           'kirkwood_PRO16_7A_diagnostics.dat', &
                           'kirkwood_PRO16_10A_diagnostics.dat' ]

   INTEGER, PARAMETER, DIMENSION(nKirkwoodRegions) :: timeSeriesUnit = [301,302,303]
   INTEGER, PARAMETER, DIMENSION(nKirkwoodRegions) :: summaryUnit = [311,312,313]
   INTEGER, PARAMETER, DIMENSION(nKirkwoodRegions) :: diagnosticUnit = [321,322,323]

   ! Compact full time series are written for every analyzed frame by default.
   ! They contain the observables needed for trimming, block analysis, plotting,
   ! and publication-level post-processing.  Summaries always use every valid
   ! analyzed frame, independent of the output stride.
   LOGICAL, PARAMETER :: writeTimeSeries = .TRUE.
   INTEGER, PARAMETER :: timeSeriesStride = 1

   ! Optional audit output containing the underlying vector/scalar sums.  This is
   ! useful while validating the estimator but is unnecessary for routine runs.
   ! Keeping it off avoids roughly doubling the formatted output volume.
   LOGICAL, PARAMETER :: writeDiagnosticSums = .FALSE.
   INTEGER, PARAMETER :: diagnosticStride = 1

   INTEGER :: nOpenMPThreads = 1
   INTEGER :: pro16OIndex = 0

   ! Thread-local scalar/vector sums.  Dimensions are (region,thread) or
   ! (xyz,region,thread).  These are reduced after each frame.
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: regionCountThread
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: regionSumUThread
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: regionSumPThread
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: regionSumUParThread
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: regionSumNThread

   ! Current-frame reduced sums and observables.
   INTEGER, DIMENSION(nKirkwoodRegions) :: regionFrameN = 0
   LOGICAL, DIMENSION(nKirkwoodRegions) :: regionFrameValid = .FALSE.
   DOUBLE PRECISION, DIMENSION(3,nKirkwoodRegions) :: regionFrameSumU = 0.0d0
   DOUBLE PRECISION, DIMENSION(nKirkwoodRegions) :: regionFrameSumP = 0.0d0
   DOUBLE PRECISION, DIMENSION(3,nKirkwoodRegions) :: regionFrameSumUPar = 0.0d0
   DOUBLE PRECISION, DIMENSION(3,nKirkwoodRegions) :: regionFrameSumN = 0.0d0
   DOUBLE PRECISION, DIMENSION(nKirkwoodRegions) :: regionFrameNormalCoherence = 0.0d0
   DOUBLE PRECISION, DIMENSION(nKirkwoodRegions) :: regionFrameGK = 0.0d0
   DOUBLE PRECISION, DIMENSION(nKirkwoodRegions) :: regionFrameGPerp = 0.0d0
   DOUBLE PRECISION, DIMENSION(nKirkwoodRegions) :: regionFrameGParallel = 0.0d0
   DOUBLE PRECISION, DIMENSION(nKirkwoodRegions) :: regionFrameResidual = 0.0d0

   ! Trajectory accumulators.  Frames with N=0 are excluded.
   INTEGER, DIMENSION(nKirkwoodRegions) :: regionValidFrames = 0
   DOUBLE PRECISION, DIMENSION(nKirkwoodRegions) :: regionNAcc = 0.0d0
   DOUBLE PRECISION, DIMENSION(nKirkwoodRegions) :: regionNormalCoherenceAcc = 0.0d0
   DOUBLE PRECISION, DIMENSION(nKirkwoodRegions) :: regionGKAcc = 0.0d0
   DOUBLE PRECISION, DIMENSION(nKirkwoodRegions) :: regionGPerpAcc = 0.0d0
   DOUBLE PRECISION, DIMENSION(nKirkwoodRegions) :: regionGParallelAcc = 0.0d0
   DOUBLE PRECISION, DIMENSION(nKirkwoodRegions) :: regionResidualAcc = 0.0d0
END MODULE kirkwood
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

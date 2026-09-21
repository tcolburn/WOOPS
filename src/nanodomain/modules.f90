!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
MODULE input
   IMPLICIT NONE
   INTEGER :: betaFlg
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
MODULE nanodomain
   IMPLICIT NONE

   ! A water is selected when its oxygen is within this distance of any
   ! beta=1 protein atom.  Assignment is then made to the nearest beta=1
   ! protein heavy atom.
   DOUBLE PRECISION, PARAMETER :: hydrationCutoff = 7.0d0
   DOUBLE PRECISION, PARAMETER :: hydrogenMassMax = 1.5d0

   ! A selected water is guaranteed to have a protein heavy atom within
   ! hydrationCutoff + its selecting hydrogen's parent-bond length.  The
   ! bound below is checked explicitly in every frame before the cell search.
   DOUBLE PRECISION, PARAMETER :: maxProteinHydrogenBond = 2.5d0
   DOUBLE PRECISION, PARAMETER :: assignmentSearchRadius = &
      hydrationCutoff + maxProteinHydrogenBond

   INTEGER :: nProteinHeavy = 0

   ! Periodic linked-cell list for beta=1 protein atoms.  The grid is rebuilt
   ! each frame, while its storage is reused whenever the box dimensions keep
   ! the same cell count.
   INTEGER :: proteinCellNx = 0
   INTEGER :: proteinCellNy = 0
   INTEGER :: proteinCellNz = 0
   INTEGER :: proteinCellCount = 0
   INTEGER, ALLOCATABLE, DIMENSION(:) :: proteinCellHead
   INTEGER, ALLOCATABLE, DIMENSION(:) :: proteinCellNext

   ! beta=1 indexing is used throughout these arrays.
   LOGICAL, ALLOCATABLE, DIMENSION(:) :: isProteinHeavy
   INTEGER, ALLOCATABLE, DIMENSION(:) :: heavyBetaIndex
   INTEGER, ALLOCATABLE, DIMENSION(:) :: parentHeavyIndex

   ! Per-frame temporary accumulators.  Only heavy-atom entries receive
   ! water observations directly.  Their means are later propagated to
   ! the heavy atom and all of its bonded protein hydrogens.
   INTEGER, ALLOCATABLE, DIMENSION(:) :: waterCountFrame
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: p1FrameSum
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: p2FrameSum
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: p21FrameSum

   ! OpenMP thread-local accumulation.  The touched-index lists let the
   ! serial reduction visit only heavy atoms actually observed by each thread.
   INTEGER :: nOpenMPThreads = 1
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: waterCountThread
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: p1ThreadSum
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: p2ThreadSum
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: p21ThreadSum
   INTEGER, ALLOCATABLE, DIMENSION(:) :: touchedCountThread
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: touchedHeavyThread
END MODULE nanodomain
!------------------------------------------------------------------------------
MODULE p1p2
   IMPLICIT NONE
   INTEGER :: nres
   INTEGER :: nPAtom
   CHARACTER(LEN=20), PARAMETER :: p1p2File = 'p1p2.dat'
   INTEGER, PARAMETER :: fIntp1p2 = 301
   INTEGER :: IsaveInt
   REAL*8, ALLOCATABLE, DIMENSION(:) :: p1Acc, p2Acc, p1Avg, p2Avg
   REAL*8, ALLOCATABLE, DIMENSION(:) :: p21Acc, p21Avg
   REAL*8, ALLOCATABLE, DIMENSION(:) :: rDifSave
   INTEGER, ALLOCATABLE, DIMENSION(:) :: nwSave, nwSaveAcc, aidbt2Save
   INTEGER, ALLOCATABLE, DIMENSION(:) :: nWatProtAtom
END MODULE p1p2
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

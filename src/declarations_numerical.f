!               *****************************
                MODULE DECLARATIONS_NUMERICAL
!               *****************************
!
!
!***********************************************************************
! 2D-DIFFUSION SOLVER - FINITE DIFFERENCES
!***********************************************************************
!
!brief    DECLARATION OF NUMERICAL VARIABLES FOR DIFFUSION SOLVER
!
!history  Sergio Castiblanco
!+        16/01/2021
!+        Translation for original Matlab implementation
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE DECLARATIONS_CSC
      IMPLICIT NONE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     LOGICALS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! DEBUGGER OPTION
!
      LOGICAL :: DEBUG
!
! ADI SOLVER OPTION
!
      LOGICAL :: ISADI
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     INTEGERS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! LOGICAL UNIT FOR WRITING OUPUTS
!      INTEGER :: LU=10102
!
! 
! NUMBER OF NODES IN DIRECTIONS X AND Y
!
      INTEGER :: NX
      INTEGER :: NY
!
! GLOBAL INDEX
!
      INTEGER, ALLOCATABLE :: POS(:)
      INTEGER, ALLOCATABLE :: MPOS(:,:)
      INTEGER, ALLOCATABLE :: MPOS2(:,:)
!
! BOUNDARY INDEXES
!
      INTEGER, ALLOCATABLE :: BOUND(:), UPBOUND(:), DOBOUND(:),
     &    RIBOUND(:), LEBOUND(:), UPBOUNDI(:), DOBOUNDI(:),
     &    RIBOUNDI(:), LEBOUNDI(:)
!
! NUMBER OF TIME STEPS
!
      INTEGER :: NT
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     DOUBLE PRECISION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! DIMENSIONS OF THE DOMAIN
!
      DOUBLE PRECISION :: XMIN
      DOUBLE PRECISION :: XMAX
      DOUBLE PRECISION :: YMIN
      DOUBLE PRECISION :: YMAX
!
! DIFFERENTIALS
!
      DOUBLE PRECISION :: DX
      DOUBLE PRECISION :: DY
!
! TIME DISCRETIZATION
!
!     INITIAL TIME
      DOUBLE PRECISION :: TO
!     FINAL TIME
      DOUBLE PRECISION :: TF
!
! GRID
!
      DOUBLE PRECISION, ALLOCATABLE :: X(:), Y(:)
!
! STIFFNESS MATRIX
!
      TYPE(CSC_OBJ) :: KM
!
! TRACER VECTOR
!
      !INITIAL CONDITION, SOLUTION AT T-1
      DOUBLE PRECISION, ALLOCATABLE :: CO(:)
      !ANALYTICA
      DOUBLE PRECISION, ALLOCATABLE :: CA(:)
      !ERROR VECTOR
      DOUBLE PRECISION, ALLOCATABLE :: ERRC(:)
      !CS, C0 FOR SOLVER 
      DOUBLE PRECISION, ALLOCATABLE :: CS(:)
      !CCS, CC FOR SOLVER
      DOUBLE PRECISION, ALLOCATABLE :: CCS(:)
      !CB, C ON BOUNDARIES
      DOUBLE PRECISION, ALLOCATABLE :: CB(:)
!
! TIME AND TIME STEP
!
      DOUBLE PRECISION, ALLOCATABLE  :: T(:)
      DOUBLE PRECISION :: DT
!
! DIFFUSION STABILITY PARAMETERS
!
      DOUBLE PRECISION :: SX
      DOUBLE PRECISION :: SY
!
! CG SOLVER OPTIONS AND PARAMETERS
      INTEGER :: MAXNITER, NITER
      DOUBLE PRECISION :: TOL
!
! TIME INTERVAL FOR WRITING
      DOUBLE PRECISION :: WTIME
      DOUBLE PRECISION :: NTIME = 0D0
!
! OPTIONS FOR ANALYTICAL
      DOUBLE PRECISION :: XZERO, YZERO, TZERO, AM, AL
!
!     ============================================
!
      END MODULE DECLARATIONS_NUMERICAL


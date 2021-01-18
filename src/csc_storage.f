!                     ******************
                      MODULE CSC_STORAGE
!                     ******************
!
!
!***********************************************************************
! STRUCTURE DECLARATION FOR CSC STRUCTURES
!***********************************************************************
!
!brief    STRUCTURE CSC FOR STORAGE VALUES, ROWS INDICES AND COLUMN
!         STARTS
!
!history  Sergio Castiblanco
!+        16/01/2021
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE DECLARATIONS_CSC
      INTERFACE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      SUBROUTINE ALL_CSC(OBJ, N, M, NZ, NAM)
        USE DECLARATIONS_CSC
        TYPE(CSC_OBJ), INTENT(INOUT) :: OBJ
        INTEGER, INTENT(IN)          :: N, M, NZ
        CHARACTER(LEN=6), INTENT(IN) :: NAM
      END SUBROUTINE
!
      SUBROUTINE LOG2INT(DIMI,LOGI,INTE)
        INTEGER, INTENT(IN) :: DIMI
        INTEGER, DIMENSION(DIMI), INTENT(OUT) :: INTE
        LOGICAL, DIMENSION(DIMI), INTENT(IN) :: LOGI
      END SUBROUTINE
!
      SUBROUTINE LOG2IND(DIMI,LOGI,INTE)
        INTEGER, INTENT(IN) :: DIMI
        INTEGER, ALLOCATABLE, INTENT(OUT) :: INTE(:)
        LOGICAL, DIMENSION(DIMI), INTENT(IN) :: LOGI
      END SUBROUTINE
! 
      SUBROUTINE CSC_DIAG(NCH,NRH,H2,D,MA,NAMMA)
        USE DECLARATIONS_CSC
        INTEGER, INTENT(IN)          :: NCH
        INTEGER, INTENT(IN)          :: NRH
        DOUBLE PRECISION, INTENT(IN), DIMENSION(NRH,NCH) :: H2
        INTEGER, INTENT(IN), DIMENSION(NCH)              :: D
        TYPE(CSC_OBJ), INTENT(INOUT)                     :: MA
        CHARACTER(LEN=6), INTENT(IN)                     :: NAMMA
      END SUBROUTINE CSC_DIAG
!
      SUBROUTINE CSC_KRON(MA,MB,MC,NAMC)
        USE DECLARATIONS_CSC
        TYPE(CSC_OBJ), INTENT(IN) :: MA, MB
        TYPE(CSC_OBJ), INTENT(INOUT) :: MC
        CHARACTER(LEN=6), INTENT(IN) :: NAMC
      END SUBROUTINE CSC_KRON
!
      SUBROUTINE CSC_SUM(MA,MB,MC,NAMC)
        USE DECLARATIONS_CSC
        TYPE(CSC_OBJ), INTENT(IN) :: MA, MB
        TYPE(CSC_OBJ), INTENT(INOUT) :: MC
        CHARACTER(LEN=6), INTENT(IN) :: NAMC
      END SUBROUTINE CSC_SUM
!
      SUBROUTINE UNION(VA,VB,VC)
        INTEGER, INTENT(IN) :: VA(:), VB(:)
        INTEGER, ALLOCATABLE, INTENT(OUT) :: VC(:)
      END SUBROUTINE UNION
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END INTERFACE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END MODULE CSC_STORAGE


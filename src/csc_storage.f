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
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END INTERFACE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END MODULE CSC_STORAGE


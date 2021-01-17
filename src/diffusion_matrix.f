!                 ***************************
                  SUBROUTINE DIFFUSION_MATRIX
!                 ***************************
     & (DX,DY,DT,VX,VY,BOUND,UPBOUND,DOBOUND,RIBOUND,LEBOUND,KM,
     &    SX,SY)
!
!***********************************************************************
! 2D-DIFFUSION SOLVER SOLVER - FINITE DIFFERENCES
!***********************************************************************
!
!brief    1) COMPUTE DIFFUSION COEFFICIENTS AND MATRIX.
!
!history  Sergio Castiblanco
!+        16/01/2021
!+        Translation for original Matlab implementation
!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| DX,DY     |-->| SPACE DIFFERENTIALS                                 |
!| DT        |-->| TIME DIFFERENTIAL                                   |
!| VX,VY     |-->| DIFFUSION COEFFICIENTS                              |
!| BOUND     |-->| BOUNDARY INDICES                                    |
!| UPBOUND   |-->| TOP BOUNDARY INDICES                                |
!| DOBOUND   |-->| BOTTOM BOUNDARY INDICES                             |
!| RIBOUND   |-->| RIGHT BOUNDARY INDICES                              |
!| LEBOUND   |-->| LEFT BOUNDARY INDICES                               |
!| KM        |<--| DIFFUSION MATRIX                                    |
!| SX,SY     |<->| DIFF COEFFICIENTS                                   |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE CSC_STORAGE
      USE DECLARATIONS_NUMERICAL, ONLY:NX,NY,DEBUG,LU 
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      DOUBLE PRECISION, INTENT(IN)                  :: DX,DY,DT,VX,VY
      INTEGER, DIMENSION((2*NX+2*NY)-4), INTENT(IN) :: BOUND
      INTEGER, DIMENSION(NX), INTENT(IN)            :: UPBOUND, DOBOUND
      INTEGER, DIMENSION(NY-2), INTENT(IN)          :: RIBOUND,LEBOUND
      DOUBLE PRECISION, INTENT(INOUT)               :: SX,SY
      TYPE(CSC_OBJ), INTENT(OUT)                    :: KM
!
! IN SUBROUTINE VARIABLES
!
      INTEGER :: I,J
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! SETTING DIFFUSSION STABILITY PARAMETERS
!
      IF(DEBUG) WRITE(LU,*) 'COMPUTING SX AND SY'
      SX = ((DT*VX)/(DX**2))
      SY = ((DT*VY)/(DY**2))
      IF(DEBUG) WRITE(LU,*) 'END COMPUTING SX AND SY'
!
! COMPUTING STTIFFNESS DIFFUSION MATRIX
!
!
! DISPLAYING COMPUTED SX AND SY
!
      WRITE(LU,*) REPEAT('~',72)
      WRITE(LU,*) 'SX VALUE: ', SX
      WRITE(LU,*) 'SY VALUE: ', SY
      WRITE(LU,*) REPEAT('~',72)
!
      END SUBROUTINE DIFFUSION_MATRIX
 

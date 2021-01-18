!                       *********************
                        SUBROUTINE POINT_DIFF
!                       *********************
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!***********************************************************************
! 2D DIFUSSION SOLVER - FINITE DIFFERENCES
!***********************************************************************
!
!brief    1) ALLOCATING AND COMPUTING INITIAL CONSTANTS AND SETUP.
!
!history  Sergio Castiblanco
!+        16/01/2021
!+        Translation for original Matlab implementation
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      USE DECLARATIONS_PHYSIC
      USE DECLARATIONS_NUMERICAL
!
      IMPLICIT NONE
!    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! IN SUBROUTINE VARIABLES
      INTEGER :: I,J,M
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ALLOCATING GRID
!
      IF(DEBUG) WRITE(LU,*) 'ALLOCATING GRID, X AND Y'
      ALLOCATE(X(NX))
      X = 0D0
      ALLOCATE(Y(NY))
      Y = 0D0
!
! ALLOCATE BOUNDARY POSITION VECTORS
!
      IF(DEBUG) WRITE(LU,*) 'ALLOCATING BOUNDARY POSITIONAL VECTORS'
      ALLOCATE(POS(NX*NY))
      POS = 1
      ALLOCATE(MPOS(NY,NX))
      MPOS = 1
      ALLOCATE(BOUND((2*NX + 2*NY)-4))
      BOUND = 1
      ALLOCATE(UPBOUND(NX))
      UPBOUND = 1
      ALLOCATE(DOBOUND(NX))
      DOBOUND = 1
      ALLOCATE(RIBOUND(NY-2))
      RIBOUND = 1
      ALLOCATE(LEBOUND(NY-2))
      LEBOUND = 1
!
! ALLOCATE TRACER VECTOR
!
      IF(DEBUG) WRITE(LU,*) 'ALLOCATING TRACER VECTOR'
      ALLOCATE(CO(NX*NY))
      CO = 0D0
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  GRID - GRID - GRID - GRID
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! COMPUTING DIFFERENTIALS
!
      DX = (XMAX - XMIN)/(NX-1)
      DY = (YMAX - YMIN)/(NY-1)
!
! COMPUTING X AND Y
!
      IF(DEBUG) WRITE(LU,*) 'COMPUTING X'
      X(1) = XMIN
      DO I=2,NX
        X(I) = X(I-1) + DX
      ENDDO
!
      IF(DEBUG) WRITE(LU,*) 'COMPUTING Y'
      Y(1) = YMIN
      DO I=2,NY
        Y(I) = Y(I-1) + DY
      ENDDO
!
! ENUMERATION AND INDEXING
!
      IF(DEBUG) WRITE(LU,*) 'COMPUTING POS (ENUMERATION)'
      DO I=1,NX*NY
        POS(I) = I
      ENDDO
      M=1
      DO J=1,NX
        DO I=1,NY
          MPOS(I,J) = POS(M)
          M = M+1
        ENDDO
      ENDDO
!
! BOUNDARIES INDEXING
!
      IF(DEBUG) WRITE(LU,*) 'COMPUTING BOUNDARIES INDEXING'
!     UP BOUNDARY
      UPBOUND = MPOS(1,:)
!
!     BOTTOM BOUNDARY
      DOBOUND = MPOS(NY,:)
!
!     RIGHT BOUNDARY
      RIBOUND = MPOS(2:NY-1,NX)
!
!     LEFT BOUNDARY
      LEBOUND = MPOS(2:NY-1,1)
!
!     ALL BOUNDARIES TOGETHER
      BOUND = [UPBOUND,DOBOUND,RIBOUND,LEBOUND]
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! END GRID - END GRID - END GRID
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TIME - TIME - TIME - TIME
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! SETTING TIME PARAMETERS
!
      IF(DEBUG) WRITE(LU,*) 'COMPUTING DT AND NT'
      DT = (0.5*MIN(DX,DY)**2)/MAX(VX,VY)
      DT = 0.025
      NT = INT(FLOOR((TF-TO)/DT))
!
! ALLOCATE TIME VECTOR
!
      IF(DEBUG) WRITE(LU,*) 'ALLOCATING AND COMPUTING TIME VECTOR'
      ALLOCATE(T(NT))
!
! FILLING TIME VECTOR
!
      T(1) = TO
      DO I=2,NT
        T(I) = T(I-1) + DT
      ENDDO
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! END TIME - END TIME - END TIME
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE POINT_DIFF






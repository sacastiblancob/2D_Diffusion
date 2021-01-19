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
      USE CSC_STORAGE
      USE DECLARATIONS_PHYSIC
      USE DECLARATIONS_NUMERICAL
!
      IMPLICIT NONE
!    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! IN SUBROUTINE VARIABLES
      INTEGER :: I,J,M
      !DT EXPLICIT TYPE LIMITATION
      DOUBLE PRECISION :: DTL
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ALLOCATING GRID
!
      IF(DEBUG) WRITE(*,*) 'ALLOCATING GRID, X AND Y'
      ALLOCATE(X(NX))
      X = 0D0
      ALLOCATE(Y(NY))
      Y = 0D0
!
! ALLOCATE BOUNDARY POSITION VECTORS
!
      IF(DEBUG) WRITE(*,*) 'ALLOCATING BOUNDARY POSITIONAL VECTORS'
      ALLOCATE(POS(NX*NY))
      POS = 1
      ALLOCATE(MPOS(NY,NX))
      MPOS = 1
!
!  EXTERNAL BOUNDARIES
      ALLOCATE(BOUND(2*NX + 2*(NY-2)))
      BOUND = 1
      ALLOCATE(UPBOUND(NX))
      UPBOUND = 1
      ALLOCATE(DOBOUND(NX))
      DOBOUND = 1
      ALLOCATE(LEBOUND(NY-2))
      LEBOUND = 1
      ALLOCATE(RIBOUND(NY-2))
      RIBOUND = 1
!  INTERNAL BOUNDARIES
!
      ALLOCATE(UPBOUNDI(NX-4))
      UPBOUNDI = 1
      ALLOCATE(DOBOUNDI(NX-4))
      DOBOUNDI = 1
      ALLOCATE(RIBOUNDI(NY-4))
      RIBOUNDI = 1
      ALLOCATE(LEBOUNDI(NY-4))
      LEBOUNDI = 1
!
! ALLOCATE TRACER VECTOR
!
      IF(DEBUG) WRITE(*,*) 'ALLOCATING TRACER VECTOR'
      ALLOCATE(CO(NX*NY))
      CO = 0D0
      ALLOCATE(CA(NX*NY))
      CA = 0D0
      ALLOCATE(ERRC(NX*NY))
      ERRC = 0D0
      ALLOCATE(CS((NX-2)*(NY-2)))
      CS = 0D0
      ALLOCATE(CCS((NX-2)*(NY-2)))
      CCS = 0D0
      ALLOCATE(CB(2*NX + 2*(NY-2)))
      CB = 0D0
!
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
      IF(DEBUG) WRITE(*,*) 'COMPUTING X'
      X(1) = XMIN
      DO I=2,NX
        X(I) = X(I-1) + DX
      ENDDO
!
      IF(DEBUG) WRITE(*,*) 'COMPUTING Y'
      Y(1) = YMIN
      DO I=2,NY
        Y(I) = Y(I-1) + DY
      ENDDO
!
! ENUMERATION AND INDEXING
!
      IF(DEBUG) WRITE(*,*) 'COMPUTING POS (ENUMERATION)'
      DO I=1,NX*NY
        POS(I) = I
      ENDDO
      M=1
      DO I = NY,1,-1
        DO J = 1,NX
          MPOS(I,J) = POS(M)
          M = M+1
        ENDDO
      ENDDO
      ! ! !DO I = 1,NY
      ! ! !  WRITE(*,*) 'MPOS',I,':',MPOS(I,:)
      ! ! !ENDDO
!
! BOUNDARIES INDEXING
!
      IF(DEBUG) WRITE(*,*) 'COMPUTING BOUNDARIES INDEXING'
!     INTERNAL BOUNDARIES
      UPBOUNDI = MPOS(2,3:NX-2)
      DOBOUNDI = MPOS(NY-1,3:NX-2)
      LEBOUNDI = MPOS(3:NY-2,2)
      RIBOUNDI = MPOS(3:NY-2,NX-1)
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
!
!  SORTING BOUNDARY INDICES VECTORS
      CALL QUICKSORT(BOUND,1,SIZE(BOUND))
      CALL QUICKSORT(LEBOUND,1,SIZE(LEBOUND))
      CALL QUICKSORT(RIBOUND,1,SIZE(RIBOUND))
      CALL QUICKSORT(UPBOUND,1,SIZE(UPBOUND))
      CALL QUICKSORT(DOBOUND,1,SIZE(DOBOUND))
      CALL QUICKSORT(LEBOUNDI,1,SIZE(LEBOUNDI))
      CALL QUICKSORT(RIBOUNDI,1,SIZE(RIBOUNDI))
      CALL QUICKSORT(UPBOUNDI,1,SIZE(UPBOUNDI))
      CALL QUICKSORT(DOBOUNDI,1,SIZE(DOBOUNDI))
!      
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
      IF(DEBUG) WRITE(*,*) 'COMPUTING DT AND NT'
      DTL = (0.5*MIN(DX,DY)**2)/MAX(VX,VY)
      IF(DT.GT.DTL) THEN
        DT = DTL
      ENDIF
      NT = INT(FLOOR((TF-TO)/DT))
!
! ALLOCATE TIME VECTOR
!
      IF(DEBUG) WRITE(*,*) 'ALLOCATING AND COMPUTING TIME VECTOR'
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






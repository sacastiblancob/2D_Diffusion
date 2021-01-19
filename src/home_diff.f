!                    *****************
                     PROGRAM HOME_DIFF
!                    *****************
!
!
!***********************************************************************
! 2D DIFUSSION SOLVER - FINITE DIFFERENCES
!***********************************************************************
!
!brief    1) MAIN PROGRAM 2D-DIFFUSION SOLVER.
!
!history  Sergio Castiblanco
!+        12/01/2021
!+        Translation for original Matlab implementation
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE DECLARATIONS_PHYSIC
      USE DECLARATIONS_NUMERICAL
!
      IMPLICIT NONE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  IN SUBROUTINE VARIABLES
!
      INTEGER :: I,J,K,TI
      INTEGER :: WTI = 0
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  READING USER INPUT FILE
!
      NAMELIST /DIFFCONF/ VX,VY,NX,NY,XMIN,XMAX,YMIN,YMAX,XZERO,YZERO,
     &    TZERO,AM,AL,TO,TF,DT,ISADI,MAXNITER,TOL,WTIME,DEBUG
!
      OPEN(10101, FILE = "diffconf.nml", STATUS = 'OLD')
      READ(10101, NML = DIFFCONF)
      CLOSE(10101)
      IF(DEBUG) WRITE(*,*) 'EXIT READING USER ENTRIES'
      IF(DEBUG) WRITE(*,*) 'DIFFCONF ', VX,VY,NX,NY,XMIN,XMAX,YMIN,
     &   ISADI, MAXNITER, TOL, DEBUG
!
!  WRITING HEADERS
!
      IF(DEBUG) WRITE(*,*) 'GOING INTO WRITE_HEADERS'
      CALL WRITE_HEADERS
      IF(DEBUG) WRITE(*,*) 'EXIT WRITE_HEADERS'
!
!  ALLOCATING MEMORY AND SETTING INITIAL CONDITION
!
      IF(DEBUG) WRITE(*,*) 'GOING INTO POINT_DIFF'
      CALL POINT_DIFF
      IF(DEBUG) WRITE(*,*) 'EXIT POINT_DIFF'
!
!  COMPUTING STIFFNESS DIFUSSION MATRIX
!
      IF(DEBUG) WRITE(*,*) 'GOING INTO DIFFUSION_MATRIX'
      CALL DIFFUSION_MATRIX(DX,DY,DT,VX,VY,KM,SX,SY) 
      IF(DEBUG) WRITE(*,*) 'EXIT DIFFUSION_MATRIX'
!
!  INITIAL CONDITION
!
      IF(DEBUG) WRITE(*,*) 'GOING INTO ANALYTICAL SOLUTION'
      K = 1
      DO J = 1,NY
        DO I = 1,NX
        CALL ANALYTICAL(X(I),Y(J),TO,0D0,0D0,VX,VY,XZERO,YZERO,
     &    TZERO,AM,AL,CO(K))
        K = K+1
        ENDDO
      ENDDO
      IF(DEBUG) WRITE(*,*) 'EXIT ANALYTICAL'
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  TIME LOOP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      DO TI = 2,SIZE(T)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF(DEBUG) WRITE(*,*) 'TIME LOOP INIT ON TIME = ',T(TI),':',TI
!
!  COMPUTING ANALYTICAL FOR T(TI)
!
      IF(DEBUG) WRITE(*,*) 'COMPUTING ANALYTICAL'
      K = 1
      DO J=1,NY
        DO I = 1,NX
        CALL ANALYTICAL(X(I),Y(J),T(TI),0D0,0D0,VX,VY,XZERO,YZERO,
     &      TZERO,AM,AL,CA(K))
        K = K + 1
        ENDDO
      ENDDO
      IF(DEBUG) WRITE(*,*) 'END ANALTICAL'
!
!  COMPUTING RIGHT HAND SIDE VECTOR
!
      IF(DEBUG) WRITE(*,*) 'GOING INTO POINT_RHS'
      CALL POINT_RHS(CO,CA,CS,CCS,CB,BOUND,UPBOUNDI,DOBOUNDI,LEBOUNDI,
     &    RIBOUNDI,VX,VY)
      IF(DEBUG) WRITE(*,*) 'EXIT FROM POINT_RHS'
!
!  SOLVING THE SYSTEM OF EQUATIONS WITH CONJUGATE GRADIENT
!
      IF(DEBUG) WRITE(*,*) 'CALLING SOLVER'
      CALL CSC_CG(KM,CS,CCS,MAXNITER,NITER,TOL,SIZE(CS))
      IF(DEBUG) WRITE(*,*) 'EXIT SOLVER'
!
!  UPDATING VARIABLES AND WRITING OUTPUTS
!
      IF(DEBUG) WRITE(*,*) 'CALLING UPDATE AND WRITING'
      CALL UPDATE_AND_WRITE(X,Y,CA,ERRC,CO,CCS,CB,BOUND,T,TI,NTIME,
     &     WTIME,NITER,TOL,WTI)
      IF(DEBUG) WRITE(*,*) 'EXIT UPDATE AND WRITING'
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  END TIME LOOP
      ENDDO
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      STOP 0
      END

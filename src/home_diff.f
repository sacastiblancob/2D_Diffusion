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
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  READING USER INPUT FILE
!
      NAMELIST /DIFFCONF/ VX,VY,NX,NY,XMIN,XMAX,YMIN,YMAX,TO,TF,ISADI,
     &     DEBUG
!
      OPEN(1010, FILE = "diffconf.nml", STATUS = 'OLD')
      READ(1010, NML = DIFFCONF)
      CLOSE(1010)
      IF(DEBUG) WRITE(LU,*) 'EXIT READING USER ENTRIES'
!
!  WRITING HEADERS
!
      IF(DEBUG) WRITE(LU,*) 'GOING INTO WRITE_HEADERS'
      CALL WRITE_HEADERS(LU)
      IF(DEBUG) WRITE(LU,*) 'EXIT WRITE_HEADERS'
!
!  ALLOCATING MEMORY AND SETTING INITIAL CONDITION
!
      IF(DEBUG) WRITE(LU,*) 'GOING INTO POINT_DIFF'
      CALL POINT_DIFF
      IF(DEBUG) WRITE(LU,*) 'EXIT POINT_DIFF'
!
!  COMPUTING STIFFNESS DIFUSSION MATRIX
!
      IF(DEBUG) WRITE(LU,*) 'GOING INTO DIFFUSION_MATRIX'
      CALL DIFFUSION_MATRIX(DX,DY,DT,VX,VY,BOUND,UPBOUND,DOBOUND,
     &     RIBOUND,LEBOUND,KM,SX,SY) 
      IF(DEBUG) WRITE(LU,*) 'EXIT DIFFUSION_MATRIX'
!
!  INITIAL CONDITION
!
      IF(DEBUG) WRITE(LU,*) 'GOING INTO ANALYTICAL SOLUTION'
      CALL ANALYTICAL(X,Y,1.5D0,DUMMYU,DUMMYU,VX,VY,CO,0D0,0D0,1D0,
     &    1D0,1D0)
      IF(DEBUG) WRITE(LU,*) 'EXIT ANALYTICAL'
!
      STOP 0
      END

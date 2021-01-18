!                    *********************
                     SUBROUTINE ANALYTICAL
!                    *********************
     & (X,Y,TIME,U,V,VX,VY,VC,XO,YO,TO,M,L)
!
!***********************************************************************
! 2D-DIFFUSION SOLVER SOLVER - FINITE DIFFERENCES
!***********************************************************************
!
!brief    1) COMPUTE DIFFUSION-ADVECTION ANALYTICAL SOLUTION.
!
!history  Sergio Castiblanco
!+        17/01/2021
!+        Translation for original Matlab implementation
!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| X,Y       |-->| COORDINATES VECTORS                                 |
!| TIME      |-->| TIME                                                |
!| U,V       |-->| VELOCITIES IN X(U) AND Y(V)                         |
!| VX,VY     |-->| DIFFUSION COEFFICIENTS                              |
!| VC        |<--| SOLUTION OF ANALYTICAL                              |
!| XO,YO,TO  |<--| CENTER OF THE ANALYTICAL SOLUTION                   |
!| M,L       |<--| CONSTANTS VALUES OF SLOPE M/L                       |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE DECLARATIONS_NUMERICAL, ONLY:NX,NY,DEBUG,LU 
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      DOUBLE PRECISION, DIMENSION(NX*NY), INTENT(IN) :: X,Y,U,V
      DOUBLE PRECISION, INTENT(IN) :: TIME, VX, VY, XO, YO, TO, M, L
      DOUBLE PRECISION, INTENT(INOUT) :: VC
!
!  IN SUBROUTINE VARIABLES
!
      INTEGER :: I
      DOUBLE PRECISION :: PI = 3.14159265358979323846D0
      DOUBLE PRECISION :: D
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      DO I = 1,(NX*NY)
      D = -1D0*(((X(I)-XO-U(I)*TIME)**2)*((4*VX*(TIME-TO))**(-1)))  
      D = D - (((Y(I)-YO-V(I)*TIME)**2)*((4*VY*(TIME-TO))**(-1)))
      D = EXP(D)
      D = D*((M/L)/((4*PI*SQRT(VX*VY))*(TIME-TO)))
      !VC(I) = D
      WRITE(*,*) 'D ', D
      ENDDO
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE ANALYTICAL





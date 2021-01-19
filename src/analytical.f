!                    *********************
                     SUBROUTINE ANALYTICAL
!                    *********************
     & (X,Y,TIME,U,V,VX,VY,XO,YO,TO,M,L,SOL)
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
!| X,Y       |-->| COORDINATES                                         |
!| TIME      |-->| TIME                                                |
!| U,V       |-->| VELOCITIES IN X(U) AND Y(V)                         |
!| VX,VY     |-->| DIFFUSION COEFFICIENTS                              |
!| SOL       |<--| SOLUTION OF ANALYTICAL                              |
!| XO,YO,TO  |-->| CENTER OF THE ANALYTICAL SOLUTION                   |
!| M,L       |-->| CONSTANTS VALUES OF SLOPE M/L                       |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      DOUBLE PRECISION, INTENT(IN) :: X,Y,U,V
      DOUBLE PRECISION, INTENT(IN) :: TIME, VX, VY, XO, YO, TO, M, L
      DOUBLE PRECISION, INTENT(OUT) :: SOL
!
!  IN SUBROUTINE VARIABLES
!
      DOUBLE PRECISION :: PI = 3.14159265358979323846D0
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      SOL = -1D0*(((X-XO-U*TIME)**2)*((4*VX*(TIME-TO))**(-1)))  
      SOL = SOL - (((Y-YO-V*TIME)**2)*((4*VY*(TIME-TO))**(-1)))
      SOL = EXP(SOL)
      SOL = SOL*((M/L)/((4*PI*SQRT(VX*VY))*(TIME-TO)))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE ANALYTICAL





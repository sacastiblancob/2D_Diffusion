!                    ********************
                     SUBROUTINE POINT_RHS
!                    ********************
     & (CO,CA,CS,CCS,CB,BOUND,UPBOUNDI,DOBOUNDI,LEBOUNDI,RIBOUNDI,VX,VY)
!
!***********************************************************************
! 2D-DIFFUSION SOLVER SOLVER - FINITE DIFFERENCES
!***********************************************************************
!
!brief    1) COMPUTE RHS VECTOR FOR SOLVING THE EQUATION.
!
!history  Sergio Castiblanco
!+        18/01/2021
!+        Translation for original Matlab implementation
!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| CO        |-->| INITIAL CONDITION / RESULT LAST TIME STEP           |
!| CA        |-->| ANALYTICAL SOLUTION                                 |
!| CS        |-->| RHS FOR THE SYSTEM OF EQUATIONS                     |
!| CCS       |-->| GUES FOR CG SOLVER, FILLED WITH CO                  |
!| CB        |-->| VALUES IN THE BOUNDARY = CA(BOUND)                  |
!| BOUND     |-->| BOUNDARY INDICES                                    |
!| UPBOUNDI  |-->| TOP INTERNAL BOUNDARY INDICES                       |
!| DOBOUNDI  |-->| BOTTOM INTERNAL BOUNDARY INDICES                    |
!| RIBOUNDI  |-->| RIGHT INTERNAL BOUNDARY INDICES                     |
!| LEBOUNDI  |-->| LEFT INTERNAL BOUNDARY INDICES                      |
!| VX,VY     |<--| DIFFUSION COEFFICIENTS                              |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE DECLARATIONS_NUMERICAL, ONLY:NX,NY,DX,DY,DT 
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      DOUBLE PRECISION, INTENT(IN), DIMENSION(NX*NY) :: CO, CA
      DOUBLE PRECISION, INTENT(OUT), DIMENSION((NX-2)*(NY-2)) :: CS,CCS
      DOUBLE PRECISION, INTENT(OUT), DIMENSION(2*NX + 2*(NY-2)) :: CB
      INTEGER, INTENT(IN), DIMENSION(2*NX + 2*(NY-2)) :: BOUND
      INTEGER, INTENT(IN), DIMENSION(NX-4) :: UPBOUNDI,DOBOUNDI
      INTEGER, INTENT(IN), DIMENSION(NY-4) :: LEBOUNDI,RIBOUNDI
      DOUBLE PRECISION, INTENT(IN) :: VX, VY
!
! IN SUBROUTINE VARIABLES
!
      INTEGER :: I,J,L
      DOUBLE PRECISION :: AX, AY
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!
!  COMPUTING AX AND AY
      AX = -(VX*DT)/(DX**2)
      AY = -(VY*DT)/(DY**2)
!
!  COMPUTING CS
      I = 1
      L = 1
      DO J = 1,NX*NY
        IF(ANY(J.EQ.BOUND)) THEN
          CB(I) = CA(J)
          I = I + 1
        ELSEIF(ANY(J.EQ.UPBOUNDI)) THEN
          CS(L) = CO(J) - AY*CA(J-NX)
          CCS(L) = CO(J)
          L = L + 1
        ELSEIF(ANY(J.EQ.RIBOUNDI)) THEN
          CS(L) = CO(J) - AX*CA(J+1)
          CCS(L) = CO(J)
          L = L + 1
        ELSEIF(ANY(J.EQ.DOBOUNDI)) THEN
          CS(L) = CO(J) - AY*CA(J+NX)
          CCS(L) = CO(J)
          L = L + 1
        ELSEIF(ANY(J.EQ.LEBOUNDI)) THEN
          CS(L) = CO(J) - AX*CA(J-1)
          CCS(L) = CO(J)
          L = L + 1
        ELSEIF(J.EQ.(NX+2)) THEN
          CS(L) = CO(J) - AY*CA(J-NX) - AX*CA(J-1)
          CCS(L) = CO(J)
          L = L + 1
        ELSEIF(J.EQ.(2*NX - 1)) THEN
          CS(L) = CO(J) - AY*CA(J-NX) - AX*CA(J+1)
          CCS(L) = CO(J)
          L = L + 1
        ELSEIF(J.EQ.(NX*NY - 2*NX + 2)) THEN
          CS(L) = CO(J) - AY*CA(J+NX) - AX*CA(J-1)
          CCS(L) = CO(J)
          L = L + 1
        ELSEIF(J.EQ.(NX*NY - NX - 1)) THEN
          CS(L) = CO(J) - AY*CA(J+NX) - AX*CA(J+1)
          CCS(L) = CO(J)
          L = L + 1
        ELSE
          CS(L) = CO(J)
          CCS(L) = CO(J)
          L = L + 1
        ENDIF
      ENDDO
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE POINT_RHS


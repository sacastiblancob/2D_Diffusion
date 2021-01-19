!                 ***************************
                  SUBROUTINE UPDATE_AND_WRITE
!                 ***************************
     & (X,Y,CA,CERR,CO,CCS,CB,BOUND,T,TI,NTIME,WTIME,NITER,TOL,WTI)
!
!***********************************************************************
! 2D-DIFFUSION SOLVER SOLVER - FINITE DIFFERENCES
!***********************************************************************
!
!brief    1) UPDATE VARIABLES, AND WRITE IF WTIME.
!
!history  Sergio Castiblanco
!+        18/01/2021
!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| X,Y       |-->| COORDINATES                                         |
!| CA        |-->| ANALYTICAL SOLUTION                                 |
!| CERR      |<->| ABSOLUTE ERROR OF SOLUTION VS ANALYTICAL            |
!| CO        |<->| SOLUTION TO BE UPDATED                              |
!| CCS       |-->| SOLUTION IN INTERNAL POINTS                         |
!| CB        |-->| SOLUTION IN BOUNDARY POINTS CB = CA(BOUND)          |
!| BOUND     |-->| BOUNDARY INDICES                                    |
!| T         |-->| VECTOR WITH TIMES                                   |
!| TI        |-->| TIME STEP AT WE ARE                                 |
!| NTIME     |<->| NEXT WRITING TIME                                   |
!| WTIME     |-->| TIME INTERVAL FOR WRITING                           |
!| NITER     |-->| NUMBER OF ITERATIONS IN SOLVER                      |
!| TOL       |-->| SOLVER TOLERANCE                                    |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE DECLARATIONS_NUMERICAL, ONLY:NX,NY,NT
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      DOUBLE PRECISION, INTENT(IN), DIMENSION(NX) :: X
      DOUBLE PRECISION, INTENT(IN), DIMENSION(NY) :: Y
      DOUBLE PRECISION, INTENT(IN), DIMENSION(NX*NY) :: CA
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(NX*NY) :: CERR, CO
      DOUBLE PRECISION, INTENT(IN), DIMENSION((NX-2)*(NY-2)) :: CCS
      DOUBLE PRECISION, INTENT(IN), DIMENSION(2*NX + 2*(NY-2)) :: CB
      INTEGER, INTENT(IN), DIMENSION(2*NX + 2*(NY-2)) :: BOUND
      DOUBLE PRECISION, INTENT(IN), DIMENSION(NT) :: T
      INTEGER, INTENT(IN) :: TI
      DOUBLE PRECISION, INTENT(INOUT) :: NTIME
      DOUBLE PRECISION, INTENT(IN) :: WTIME
      INTEGER, INTENT(IN) :: NITER
      DOUBLE PRECISION, INTENT(IN) :: TOL
      INTEGER, INTENT(INOUT) :: WTI
!
! IN SUBROUTINE VARIABLES
!
      INTEGER :: I,J,L,K
      DOUBLE PRECISION :: ERR,C1,C2,C3,C4
      CHARACTER(LEN=80) :: FOUT
      CHARACTER(LEN=9)  :: FMT2
      CHARACTER(LEN=30) :: FSPEC
      DOUBLE PRECISION :: EPSI = 1E-6
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  WRITING INITIAL CONDITION IF TI == 2
!
      IF(TI.EQ.2) THEN
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! WRITING IN TERMINAL
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        WRITE(*,*) REPEAT('~',72)
        WRITE(*,*) 'INIT DONE'
        WRITE(*,*) 'WRITING INTIAL CONDITION'
        WRITE(*,*) REPEAT('~',72)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        FMT2 = '(i4.4,a4)'
        FSPEC = './res/2D_Diffusion_'
        FOUT = FSPEC
!
        WRITE(UNIT=FOUT(20:30),FMT=FMT2) WTI,'.dat'
!
        OPEN(WTI,FILE=FOUT)
!
!  HEADER
        WRITE(WTI,*) 'TITLE = "DIFFUSION 2D - FD'
        WRITE(WTI,*) 'VARIABLES = "X" "Y" "C" "ERR"'
        WRITE(WTI,*) ' ZONE F=POINT, I=',NX,', J= ',NY
!
!  WRITING RESULTS
        K = 1
        C4 = 0
        DO J = 1,NY
          DO I = 1,NX
            C1 = X(I)
            C2 = Y(J)
            C3 = CO(K)
            WRITE(WTI,'(*(F14.8))') C1, C2, C3, C4
            K = K + 1
          ENDDO
        ENDDO
        CLOSE(WTI)
        WTI = WTI + 10
        NTIME = T(1) + WTIME
      ENDIF
!
!  UPDATING CO
!
      I = 1
      L = 1
      DO J = 1,NX*NY
        IF(ANY(J.EQ.BOUND)) THEN
          CO(J) = CB(I)
          I = I + 1
        ELSE
          CO(J) = CCS(L)
          L = L + 1
        ENDIF
      ENDDO
!
!  COMPUTING CR
!
      DO J = 1,NX*NY
        CERR = ABS(CA(J) - CO(J))
      ENDDO
      ERR = NORM2(CERR)
!
!  WRITING IF IS TIME TO WRITE
!
      IF((T(TI).GT.(NTIME-EPSI)).OR.(T(TI).EQ.T(NT))) THEN
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! WRITING IN TERMINAL
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        WRITE(*,*) REPEAT('~',72)
        WRITE(*,*) 'TIME: ',T(TI),'(s)'
        WRITE(*,*) 'TIME STEP: ',TI,' OF ',SIZE(T)
        WRITE(*,*) 'LAST SOLVER TOL: ', TOL
        WRITE(*,*) 'LAST SOLVER NUMBER OF ITERATIONS: ',NITER
        WRITE(*,*) 'ERROR (NORM2) WITH ANALYTICAL: ',ERR
        WRITE(*,*) 'WRITING FILE NUMBER: ', WTI
        WRITE(*,*) REPEAT('~',72)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  WRITING FILES
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        FMT2 = '(i4.4,a4)'
        FSPEC = './res/2D_Diffusion_'
        FOUT = FSPEC
!
        WRITE(UNIT=FOUT(20:30),FMT=FMT2) WTI,'.dat'
!
! HERE IS CALLING STOP 1 WITHOUT REASON
        OPEN(WTI,FILE=FOUT)
!
!  HEADER
        WRITE(WTI,*) 'TITLE = "DIFFUSION 2D - FD'
        WRITE(WTI,*) 'VARIABLES = "X" "Y" "C" "ERR"'
        WRITE(WTI,*) ' ZONE F=POINT, I=',NX,', J= ',NY
!
!  WRITING RESULTS
        K = 1
        DO J = 1,NY
          DO I = 1,NX
            C1 = X(I)
            C2 = Y(J)
            C3 = CO(K)
            C4 = CERR(K)
            WRITE(WTI,'(*(F14.8))') C1, C2, C3, C4
            K = K + 1
          ENDDO
        ENDDO
        CLOSE(WTI)
        WTI = WTI + 10
        NTIME = NTIME + WTIME
      ENDIF
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE UPDATE_AND_WRITE


!                     *******************
                      SUBROUTINE CSC_DIAG
!                     *******************
     & (H2,D)
!
!***********************************************************************
! 2D-DIFFUSION SOLVER SOLVER - FINITE DIFFERENCES
!***********************************************************************
!
!brief    1) CREATES A MATRIX WITH CSC STORAGE BY DIAGONALS
!
!history  Sergio Castiblanco
!+        16/01/2021
!+        Translation for original Matlab implementation
!
!+
!!
!!This function make a matrix in CSC storage by diagonals, the entries are:
!!      h2 --> matrix whose columns are gonna be the diagonals of A
!!      d ---> row-vector with the numbers of the diagonals in h2 cols
!!             0 for the main diagonal, >0 for lower part and >0 for upper
!!             none of the diagonals are mandatory
!!
!!      warning: if you have zeros in h2 you should remove them with another
!!          function.
!!      NOTE: This function computes a square matrix with the given
!!      diagonals.
!!
!!      Example:
!!                    h2              d
!!               | 6  1  6 |      [-1 0 1]
!!               | 7  2  7 |
!!               | 8  3  8 |
!!               | 9  4  9 |
!!               | 10 5 10 |
!!      Returns:
!!                       A
!!              | 1  6  0  0  0 |
!!              | 6  2  7  0  0 |
!!              | 0  7  3  8  0 |
!!              | 0  0  8  4  9 |
!!              | 0  0  0  9  5 |
!!
!!      Note that 10's of h2(5,1) and h2(5,3) have dissapeared, those could
!!      be zero anyway
!!
!!      Of course A is stored in CSC, so the function returns
!!      Av, Ar, Ac (values, rows, columns)
!!
!!      Sergio A. Castiblanco B. - Avanced Numerical Methods
!!      Pontificia Universidad Javeriana - Bogota
!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| HD        |-->| MATRIX WITH DIAGONALS                               |
!| D         |-->| VECTOR WITH DIAGONALS POSITIONS                     |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE CSC_STORAGE
      USE DECLARATIONS_NUMERICAL, ONLY:NX,NY,DEBUG,LU 
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!

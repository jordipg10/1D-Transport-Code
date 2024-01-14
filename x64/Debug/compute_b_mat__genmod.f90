        !COMPILER-GENERATED INTERFACE MODULE: Sun Jan 14 13:10:51 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_B_MAT__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_B_MAT(THIS,THETA,B_MAT,K)
              USE PDE_TRANSIENT_M
              CLASS (PDE_1D_TRANSIENT_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: THETA
              TYPE (TRIDIAG_MATRIX_C), INTENT(OUT) :: B_MAT
              INTEGER(KIND=4) ,OPTIONAL, INTENT(IN) :: K
            END SUBROUTINE COMPUTE_B_MAT
          END INTERFACE 
        END MODULE COMPUTE_B_MAT__genmod

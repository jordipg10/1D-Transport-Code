        !COMPILER-GENERATED INTERFACE MODULE: Tue Jan 30 15:09:03 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_A_MAT_LIN_SYST__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_A_MAT_LIN_SYST(THIS,THETA,A_MAT,K)
              USE PDE_TRANSIENT_M
              CLASS (PDE_1D_TRANSIENT_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: THETA
              TYPE (TRIDIAG_MATRIX_C), INTENT(OUT) :: A_MAT
              INTEGER(KIND=4) ,OPTIONAL, INTENT(IN) :: K
            END SUBROUTINE COMPUTE_A_MAT_LIN_SYST
          END INTERFACE 
        END MODULE COMPUTE_A_MAT_LIN_SYST__genmod

        !COMPILER-GENERATED INTERFACE MODULE: Sun Jan 14 13:10:56 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_B_LIN_SYST__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_B_LIN_SYST(THIS,THETA,CONC_OLD,B,K)
              USE BCS_SUBROUTINES_M
              CLASS (PDE_1D_TRANSIENT_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: THETA
              REAL(KIND=8), INTENT(IN) :: CONC_OLD(:)
              REAL(KIND=8), INTENT(OUT) :: B(:)
              INTEGER(KIND=4) ,OPTIONAL, INTENT(IN) :: K
            END SUBROUTINE COMPUTE_B_LIN_SYST
          END INTERFACE 
        END MODULE COMPUTE_B_LIN_SYST__genmod

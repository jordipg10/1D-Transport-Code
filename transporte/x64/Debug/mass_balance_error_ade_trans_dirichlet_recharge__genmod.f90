        !COMPILER-GENERATED INTERFACE MODULE: Mon Mar 11 10:51:53 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MASS_BALANCE_ERROR_ADE_TRANS_DIRICHLET_RECHARGE__genmod
          INTERFACE 
            FUNCTION MASS_BALANCE_ERROR_ADE_TRANS_DIRICHLET_RECHARGE(   &
     &THIS,CONC_OLD,CONC_NEW,DELTA_T,DELTA_X) RESULT(MASS_BAL_ERR)
              USE TRANSPORT_TRANSIENT_M
              CLASS (TRANSPORT_1D_TRANSIENT_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: CONC_OLD(:)
              REAL(KIND=8), INTENT(IN) :: CONC_NEW(:)
              REAL(KIND=8), INTENT(IN) :: DELTA_T
              REAL(KIND=8), INTENT(IN) :: DELTA_X
              REAL(KIND=8) :: MASS_BAL_ERR
            END FUNCTION MASS_BALANCE_ERROR_ADE_TRANS_DIRICHLET_RECHARGE
          END INTERFACE 
        END MODULE                                                      &
     &MASS_BALANCE_ERROR_ADE_TRANS_DIRICHLET_RECHARGE__genmod
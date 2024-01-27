        !COMPILER-GENERATED INTERFACE MODULE: Sat Jan 27 13:27:19 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MASS_BALANCE_ERROR_ADE_STAT_DIRICHLET_DISCHARGE__genmod
          INTERFACE 
            FUNCTION MASS_BALANCE_ERROR_ADE_STAT_DIRICHLET_DISCHARGE(   &
     &THIS) RESULT(MASS_BAL_ERR)
              USE TRANSPORT_M
              CLASS (TRANSPORT_1D_C), INTENT(IN) :: THIS
              REAL(KIND=8) :: MASS_BAL_ERR
            END FUNCTION MASS_BALANCE_ERROR_ADE_STAT_DIRICHLET_DISCHARGE
          END INTERFACE 
        END MODULE                                                      &
     &MASS_BALANCE_ERROR_ADE_STAT_DIRICHLET_DISCHARGE__genmod

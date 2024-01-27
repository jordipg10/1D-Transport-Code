        !COMPILER-GENERATED INTERFACE MODULE: Sat Jan 27 13:27:20 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_MIXING_RATIOS_DELTA_T_HOMOG__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_MIXING_RATIOS_DELTA_T_HOMOG(THIS,THETA)
              USE BCS_SUBROUTINES_M
              CLASS (PDE_1D_TRANSIENT_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: THETA
            END SUBROUTINE COMPUTE_MIXING_RATIOS_DELTA_T_HOMOG
          END INTERFACE 
        END MODULE COMPUTE_MIXING_RATIOS_DELTA_T_HOMOG__genmod

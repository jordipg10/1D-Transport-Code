        !COMPILER-GENERATED INTERFACE MODULE: Sun Jan 14 13:10:49 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_MIXING_RATIOS__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_MIXING_RATIOS(THIS,THETA,K)
              USE BCS_SUBROUTINES_M
              CLASS (DIFFUSION_1D_TRANSIENT_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: THETA
              INTEGER(KIND=4) ,OPTIONAL, INTENT(IN) :: K
            END SUBROUTINE COMPUTE_MIXING_RATIOS
          END INTERFACE 
        END MODULE COMPUTE_MIXING_RATIOS__genmod

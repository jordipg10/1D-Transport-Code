        !COMPILER-GENERATED INTERFACE MODULE: Thu Aug 14 18:29:14 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE UPDATE_MIXING_RATIOS_DELTA_T_HOMOG__genmod
          INTERFACE 
            SUBROUTINE UPDATE_MIXING_RATIOS_DELTA_T_HOMOG(THIS)
              USE TRANSPORT_STAB_PARAMS_M
              USE TRANSPORT_PROPERTIES_HETEROG_M
              USE STABILITY_PARAMETERS_M
              USE DIFF_STAB_PARAMS_M
              USE PROPERTIES_M
              USE DIFF_PROPS_HETEROG_M
              USE CONC_M
              USE TIME_DISCR_M
              USE VECTORS_M
              USE MATRICES_M
              USE TIME_FCT_M
              USE BCS_M
              USE TARGET_M
              USE SPATIAL_DISCR_M
              USE PDE_M
              USE PDE_TRANSIENT_M
              USE DIFFUSION_TRANSIENT_M
              USE TRANSPORT_TRANSIENT_M, ONLY :                         &
     &          TRANSPORT_1D_TRANSIENT_C,                               &
     &          PDE_1D_TRANSIENT_C,                                     &
     &          DIAG_MATRIX_C,                                          &
     &          TRIDIAG_MATRIX_C,                                       &
     &          TIME_DISCR_HOMOG_C,                                     &
     &          PROD_TRIDIAG_MAT_MAT,                                   &
     &          COMPUTE_INVERSE_TRIDIAG_MATRIX
              CLASS (TRANSPORT_1D_TRANSIENT_C) :: THIS
            END SUBROUTINE UPDATE_MIXING_RATIOS_DELTA_T_HOMOG
          END INTERFACE 
        END MODULE UPDATE_MIXING_RATIOS_DELTA_T_HOMOG__genmod

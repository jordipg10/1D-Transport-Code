        !COMPILER-GENERATED INTERFACE MODULE: Mon Feb  5 13:25:34 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE REFINE_MESH_HETEROG__genmod
          INTERFACE 
            SUBROUTINE REFINE_MESH_HETEROG(THIS,CONC,CONC_EXT,REL_TOL)
              USE SPATIAL_DISCR_1D_M
              CLASS (MESH_1D_EULER_HETEROG_C) :: THIS
              REAL(KIND=8), INTENT(INOUT) :: CONC(:,:)
              REAL(KIND=8), INTENT(INOUT) :: CONC_EXT(:,:)
              REAL(KIND=8), INTENT(IN) :: REL_TOL
            END SUBROUTINE REFINE_MESH_HETEROG
          END INTERFACE 
        END MODULE REFINE_MESH_HETEROG__genmod

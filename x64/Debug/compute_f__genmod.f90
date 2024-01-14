        !COMPILER-GENERATED INTERFACE MODULE: Sun Jan 14 13:10:54 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_F__genmod
          INTERFACE 
            FUNCTION COMPUTE_F(THIS,K) RESULT(F)
              USE TRANSPORT_TRANSIENT_M
              CLASS (PDE_1D_TRANSIENT_C) :: THIS
              INTEGER(KIND=4) ,OPTIONAL, INTENT(IN) :: K
              REAL(KIND=8) ,ALLOCATABLE :: F(:)
            END FUNCTION COMPUTE_F
          END INTERFACE 
        END MODULE COMPUTE_F__genmod

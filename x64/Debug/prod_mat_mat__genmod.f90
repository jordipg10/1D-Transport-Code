        !COMPILER-GENERATED INTERFACE MODULE: Tue Jan 30 15:12:55 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PROD_MAT_MAT__genmod
          INTERFACE 
            FUNCTION PROD_MAT_MAT(THIS,B_MAT) RESULT(C_MAT)
              USE MATRICES_M
              CLASS (MATRIX_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: B_MAT(:,:)
              REAL(KIND=8) ,ALLOCATABLE :: C_MAT(:,:)
            END FUNCTION PROD_MAT_MAT
          END INTERFACE 
        END MODULE PROD_MAT_MAT__genmod

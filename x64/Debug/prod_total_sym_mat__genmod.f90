        !COMPILER-GENERATED INTERFACE MODULE: Tue Jan 30 15:12:53 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PROD_TOTAL_SYM_MAT__genmod
          INTERFACE 
            FUNCTION PROD_TOTAL_SYM_MAT(A,Y0,B,TIME) RESULT(Y)
              USE MATRICES_M
              CLASS (SQ_MATRIX_C), INTENT(IN) :: A
              REAL(KIND=8), INTENT(IN) :: Y0(:)
              REAL(KIND=8), INTENT(IN) :: B(:)
              REAL(KIND=8), INTENT(IN) :: TIME
              REAL(KIND=8) ,ALLOCATABLE :: Y(:)
            END FUNCTION PROD_TOTAL_SYM_MAT
          END INTERFACE 
        END MODULE PROD_TOTAL_SYM_MAT__genmod

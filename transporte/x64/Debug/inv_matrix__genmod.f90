        !COMPILER-GENERATED INTERFACE MODULE: Mon Feb  5 13:25:31 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INV_MATRIX__genmod
          INTERFACE 
            SUBROUTINE INV_MATRIX(A,TOL,INV)
              REAL(KIND=8), INTENT(IN) :: A(:,:)
              REAL(KIND=8), INTENT(IN) :: TOL
              REAL(KIND=8), INTENT(OUT) :: INV(:,:)
            END SUBROUTINE INV_MATRIX
          END INTERFACE 
        END MODULE INV_MATRIX__genmod
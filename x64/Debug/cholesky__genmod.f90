        !COMPILER-GENERATED INTERFACE MODULE: Tue Jan 30 15:12:24 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CHOLESKY__genmod
          INTERFACE 
            SUBROUTINE CHOLESKY(A,L)
              REAL(KIND=8), INTENT(IN) :: A(:,:)
              REAL(KIND=8), INTENT(OUT) :: L(:,:)
            END SUBROUTINE CHOLESKY
          END INTERFACE 
        END MODULE CHOLESKY__genmod
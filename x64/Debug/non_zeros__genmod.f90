        !COMPILER-GENERATED INTERFACE MODULE: Tue Jan 30 15:12:48 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE NON_ZEROS__genmod
          INTERFACE 
            SUBROUTINE NON_ZEROS(A,N,IND)
              REAL(KIND=8), INTENT(IN) :: A(:,:)
              INTEGER(KIND=4), INTENT(OUT) :: N
              INTEGER(KIND=4), INTENT(OUT) :: IND(:,:)
            END SUBROUTINE NON_ZEROS
          END INTERFACE 
        END MODULE NON_ZEROS__genmod

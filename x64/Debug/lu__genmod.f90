        !COMPILER-GENERATED INTERFACE MODULE: Sat Jan 27 13:27:01 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LU__genmod
          INTERFACE 
            SUBROUTINE LU(A,L,U)
              REAL(KIND=8), INTENT(IN) :: A(:,:)
              REAL(KIND=8), INTENT(OUT) :: L(:,:)
              REAL(KIND=8), INTENT(OUT) :: U(:,:)
            END SUBROUTINE LU
          END INTERFACE 
        END MODULE LU__genmod
        !COMPILER-GENERATED INTERFACE MODULE: Mon Feb  5 13:25:15 2024
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

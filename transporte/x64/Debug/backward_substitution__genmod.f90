        !COMPILER-GENERATED INTERFACE MODULE: Mon Feb  5 13:25:15 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE BACKWARD_SUBSTITUTION__genmod
          INTERFACE 
            SUBROUTINE BACKWARD_SUBSTITUTION(U,B,X)
              REAL(KIND=8), INTENT(IN) :: U(:,:)
              REAL(KIND=8), INTENT(IN) :: B(:)
              REAL(KIND=8), INTENT(OUT) :: X(:)
            END SUBROUTINE BACKWARD_SUBSTITUTION
          END INTERFACE 
        END MODULE BACKWARD_SUBSTITUTION__genmod
        !COMPILER-GENERATED INTERFACE MODULE: Mon Feb  5 13:25:33 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE QR_HOUSEHOLDER__genmod
          INTERFACE 
            SUBROUTINE QR_HOUSEHOLDER(A,Q,R)
              REAL(KIND=8), INTENT(IN) :: A(:,:)
              REAL(KIND=8), INTENT(OUT) :: Q(:,:)
              REAL(KIND=8), INTENT(OUT) :: R(:,:)
            END SUBROUTINE QR_HOUSEHOLDER
          END INTERFACE 
        END MODULE QR_HOUSEHOLDER__genmod

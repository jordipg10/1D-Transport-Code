        !COMPILER-GENERATED INTERFACE MODULE: Tue Jan 30 15:09:06 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LEGENDRE_POLY__genmod
          INTERFACE 
            FUNCTION LEGENDRE_POLY(N,X) RESULT(P)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: X(:)
              REAL(KIND=8) ,ALLOCATABLE :: P(:)
            END FUNCTION LEGENDRE_POLY
          END INTERFACE 
        END MODULE LEGENDRE_POLY__genmod

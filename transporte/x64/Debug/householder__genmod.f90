        !COMPILER-GENERATED INTERFACE MODULE: Mon Feb  5 13:25:38 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE HOUSEHOLDER__genmod
          INTERFACE 
            FUNCTION HOUSEHOLDER(X)
              REAL(KIND=8), INTENT(IN) :: X(:)
              REAL(KIND=8) ,ALLOCATABLE :: HOUSEHOLDER(:,:)
            END FUNCTION HOUSEHOLDER
          END INTERFACE 
        END MODULE HOUSEHOLDER__genmod

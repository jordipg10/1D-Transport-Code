        !COMPILER-GENERATED INTERFACE MODULE: Tue Jan 30 15:08:36 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE STR_BELONGS_TO__genmod
          INTERFACE 
            SUBROUTINE STR_BELONGS_TO(STRING,ARRAY,FLAG,INDEX)
              CHARACTER(*), INTENT(IN) :: STRING
              CHARACTER(*), INTENT(IN) :: ARRAY(:)
              LOGICAL(KIND=4), INTENT(OUT) :: FLAG
              INTEGER(KIND=4) ,OPTIONAL, INTENT(OUT) :: INDEX
            END SUBROUTINE STR_BELONGS_TO
          END INTERFACE 
        END MODULE STR_BELONGS_TO__genmod

  é4  i   k820309    Ö          24.0        Ųše                                                                                                          
       C:\Users\user2319\source\repos\jordipg10\1D-Transport-Code\metodos_numericos\prod_diag_tridiag_mat.f90 PROD_DIAG_TRIDIAG_MAT &         @                                                                 #TRIDIAG_MATRIX_C    #DIAG_MATRIX_C    #A N   #B O   #TRIDIAG_MATRIX_C                      @              Ā                       '                   #TRIDIAG_SYM_MATRIX_C    #SUPER 9   #SET_TRIDIAG_MATRIX :   #COMPUTE_TRANSPOSE_TRIDIAG_MATRIX A   #COMPUTE_INVERSE_TRIDIAG_MATRIX E   #PROD_TRIDIAG_MAT_MAT J                 $                                           @                      #TRIDIAG_SYM_MATRIX_C                       @              Ā                       '@                   #DIAG_MATRIX_C    #SUB 3   #CHECK_EIGENVECTORS_TRIDIAG_SYM_MATRIX 4                 $                                           ø                       #DIAG_MATRIX_C                      @              Ā                       'ø                    #SQ_MATRIX_C    #DIAG -   #SET_DIAG_MATRIX .                 $                                           °                       #SQ_MATRIX_C                       @              @                       '°                    #MATRIX_C 	   #DIM #   #EIGENVALUES $   #EIGENVECTORS %   #COMPUTE_EIGENVALUES &   #COMPUTE_EIGENVECTORS *                 $                                      	                             #MATRIX_C 
                      @                                  
     '                      #ALLOCATE_MATRIX    #PROD_MAT_VEC    #GET_DIAG    #GET_SUB    #GET_SUPER    #COMPUTE_NORM_INF    #COMPUTE_NORM_1     1         Ā    $                                                       #ALLOCATE_MATRIX    #         @                                                                #THIS    #N                                                                            #PROD_DIAG_TRIDIAG_MAT%MATRIX_C              
                                                     1         Ā    $                                                      #PROD_MAT_VEC    (        D                                                                             
    #THIS    #B              &                                                     
                                                             #PROD_DIAG_TRIDIAG_MAT%MATRIX_C              
                                                            
              &                                           1         Ā    $                                                      #GET_DIAG    (        D                                                             %                
    #THIS              &                                                     
                                                              #PROD_DIAG_TRIDIAG_MAT%MATRIX_C    1         Ā    $                                                      #GET_SUB    (        D                                                             &                
    #THIS              &                                                     
                                                              #PROD_DIAG_TRIDIAG_MAT%MATRIX_C    1         Ā    $                                                      #GET_SUPER    (        D                                                             '                
    #THIS              &                                                     
                                                              #PROD_DIAG_TRIDIAG_MAT%MATRIX_C    1         Ā    $                                                      #COMPUTE_NORM_INF    %         @                                                              
       #THIS              
                                                              #PROD_DIAG_TRIDIAG_MAT%MATRIX_C    1         Ā    $                                                       #COMPUTE_NORM_1 !   %         @                                          !                    
       #THIS "             
                                          "                    #PROD_DIAG_TRIDIAG_MAT%MATRIX_C                                                       #                                                                     $                              
            &                                                                                            %            P                  
            &                   &                                           1         Ā    $                                    &                   #COMPUTE_EIGENVALUES '   #         @                                            '     	               #THIS (                                                     (     °               #PROD_DIAG_TRIDIAG_MAT%SQ_MATRIX_C )   1         Ā    $                                    *               	    #COMPUTE_EIGENVECTORS +   #         @                                            +     	               #THIS ,                                                     ,     °               #PROD_DIAG_TRIDIAG_MAT%SQ_MATRIX_C )                                                   -            °                  
            &                                           1         Ā    $                                    .               
    #SET_DIAG_MATRIX /   #         @                                            /                    #THIS 0   #DIAG 2                                                      0     ø               #PROD_DIAG_TRIDIAG_MAT%DIAG_MATRIX_C 1             
                                          2                   
 !             &                                                                                            3            ø                  
            &                                           1         Ā    $                                    4                   #CHECK_EIGENVECTORS_TRIDIAG_SYM_MATRIX 5   #         @                                            5     	               #THIS 6   #TOLERANCE 8             
                                         6     @             #PROD_DIAG_TRIDIAG_MAT%TRIDIAG_SYM_MATRIX_C 7             
                                         8     
                                                       9            @                 
            &                                           1         Ā    $                                    :                   #SET_TRIDIAG_MATRIX ;   #         @                                            ;                    #THIS <   #SUB >   #DIAG ?   #SUPER @                                                      <                   #PROD_DIAG_TRIDIAG_MAT%TRIDIAG_MATRIX_C =             
 @                                       >                   
 "             &                                                     
 @                                       ?                   
 #             &                                                     
 @                                       @                   
 $             &                                           1         Ā    $                                    A                   #COMPUTE_TRANSPOSE_TRIDIAG_MATRIX B   #         @                                            B                    #THIS C   #TRANSPOSE D             
                                          C                  #PROD_DIAG_TRIDIAG_MAT%TRIDIAG_MATRIX_C =                                                      D                   #PROD_DIAG_TRIDIAG_MAT%TRIDIAG_MATRIX_C =   1         Ā    $                                    E                   #COMPUTE_INVERSE_TRIDIAG_MATRIX F   #         @                                            F     	               #THIS G   #TOL H   #INV_MAT I             
                                         G                  #PROD_DIAG_TRIDIAG_MAT%TRIDIAG_MATRIX_C =             
                                         H     
                                                        I                   
                &                   &                                           1         Ā    $                                   J                   #PROD_TRIDIAG_MAT_MAT K   (        D                                          K                                   
    #THIS L   #B_MAT M             &                   &                                                     
                                         L                  #PROD_DIAG_TRIDIAG_MAT%TRIDIAG_MATRIX_C =             
                                         M                   
              &                   &                                                     
                                        N     ø              #DIAG_MATRIX_C              
                                        O                  #TRIDIAG_MATRIX_C                       @               Ā                  =     '                   #TRIDIAG_SYM_MATRIX_C P   #SUPER d   #SET_TRIDIAG_MATRIX e   #COMPUTE_TRANSPOSE_TRIDIAG_MATRIX f   #COMPUTE_INVERSE_TRIDIAG_MATRIX g   #PROD_TRIDIAG_MAT_MAT h                 $                                       P     @                     #PROD_DIAG_TRIDIAG_MAT%TRIDIAG_SYM_MATRIX_C 7                      @               Ā                  7     '@                   #DIAG_MATRIX_C Q   #SUB b   #CHECK_EIGENVECTORS_TRIDIAG_SYM_MATRIX c                 $                                       Q     ø                      #PROD_DIAG_TRIDIAG_MAT%DIAG_MATRIX_C 1                      @               Ā                  1     'ø                    #SQ_MATRIX_C R   #DIAG `   #SET_DIAG_MATRIX a                 $                                       R     °                      #PROD_DIAG_TRIDIAG_MAT%SQ_MATRIX_C )                      @               @                  )     '°                    #MATRIX_C S   #DIM [   #EIGENVALUES \   #EIGENVECTORS ]   #COMPUTE_EIGENVALUES ^   #COMPUTE_EIGENVECTORS _                 $                                       S                            #PROD_DIAG_TRIDIAG_MAT%MATRIX_C                       @                                       '                      #ALLOCATE_MATRIX T   #PROD_MAT_VEC U   #GET_DIAG V   #GET_SUB W   #GET_SUPER X   #COMPUTE_NORM_INF Y   #COMPUTE_NORM_1 Z   1         Ā    $                                     T                  #ALLOCATE_MATRIX    1         Ā    $                                    U                  #PROD_MAT_VEC    1         Ā    $                                    V                  #GET_DIAG    1         Ā    $                                    W                  #GET_SUB    1         Ā    $                                    X                  #GET_SUPER    1         Ā    $                                    Y                  #COMPUTE_NORM_INF    1         Ā    $                                    Z                  #COMPUTE_NORM_1 !                                                       [                                                                     \                             
            &                                                                                             ]            P                 
            &                   &                                           1         Ā    $                                     ^                  #COMPUTE_EIGENVALUES '   1         Ā    $                                     _              	    #COMPUTE_EIGENVECTORS +                                                    `            °                 
            &                                           1         Ā    $                                     a              
    #SET_DIAG_MATRIX /                                                     b            ø                 
            &                                           1         Ā    $                                     c                  #CHECK_EIGENVECTORS_TRIDIAG_SYM_MATRIX 5                                                     d            @                
            &                                           1         Ā    $                                     e                  #SET_TRIDIAG_MATRIX ;   1         Ā    $                                     f                  #COMPUTE_TRANSPOSE_TRIDIAG_MATRIX B   1         Ā    $                                     g                  #COMPUTE_INVERSE_TRIDIAG_MATRIX F   1         Ā    $                                    h                  #PROD_TRIDIAG_MAT_MAT K                fn#fn &   %  Ĩ   `   PROD_DIAG_TRIDIAG_MAT ,   Ę  ų   a   TRIDIAG_MATRIX_C+MATRICES_M A   Ã  r   a   TRIDIAG_MATRIX_C%TRIDIAG_SYM_MATRIX_C+MATRICES_M 0   5     a   TRIDIAG_SYM_MATRIX_C+MATRICES_M >   Ô  k   a   TRIDIAG_SYM_MATRIX_C%DIAG_MATRIX_C+MATRICES_M )   ?     a   DIAG_MATRIX_C+MATRICES_M 5   Į  i   a   DIAG_MATRIX_C%SQ_MATRIX_C+MATRICES_M '   0  Å   a   SQ_MATRIX_C+MATRICES_M 0   õ  f   a   SQ_MATRIX_C%MATRIX_C+MATRICES_M $   [  Ķ   a   MATRIX_C+MATRICES_M 4   .  e   a   MATRIX_C%ALLOCATE_MATRIX+MATRICES_M A     a       PROD_DIAG_TRIDIAG_MAT%ALLOCATE_MATRIX+MATRICES_M F   ô  t   a   PROD_DIAG_TRIDIAG_MAT%ALLOCATE_MATRIX%THIS+MATRICES_M C   h  H   a   PROD_DIAG_TRIDIAG_MAT%ALLOCATE_MATRIX%N+MATRICES_M 1   °  b   a   MATRIX_C%PROD_MAT_VEC+MATRICES_M >   	  ĩ       PROD_DIAG_TRIDIAG_MAT%PROD_MAT_VEC+MATRICES_M C   Į	  t   a   PROD_DIAG_TRIDIAG_MAT%PROD_MAT_VEC%THIS+MATRICES_M @   ;
     a   PROD_DIAG_TRIDIAG_MAT%PROD_MAT_VEC%B+MATRICES_M -   Ī
  ^   a   MATRIX_C%GET_DIAG+MATRICES_M :   -  Ž       PROD_DIAG_TRIDIAG_MAT%GET_DIAG+MATRICES_M ?   Û  t   a   PROD_DIAG_TRIDIAG_MAT%GET_DIAG%THIS+MATRICES_M ,   O  ]   a   MATRIX_C%GET_SUB+MATRICES_M 9   Ŧ  Ž       PROD_DIAG_TRIDIAG_MAT%GET_SUB+MATRICES_M >   Z  t   a   PROD_DIAG_TRIDIAG_MAT%GET_SUB%THIS+MATRICES_M .   Î  _   a   MATRIX_C%GET_SUPER+MATRICES_M ;   -  Ž       PROD_DIAG_TRIDIAG_MAT%GET_SUPER+MATRICES_M @   Û  t   a   PROD_DIAG_TRIDIAG_MAT%GET_SUPER%THIS+MATRICES_M 5   O  f   a   MATRIX_C%COMPUTE_NORM_INF+MATRICES_M B   ĩ  b       PROD_DIAG_TRIDIAG_MAT%COMPUTE_NORM_INF+MATRICES_M G     t   a   PROD_DIAG_TRIDIAG_MAT%COMPUTE_NORM_INF%THIS+MATRICES_M 3     d   a   MATRIX_C%COMPUTE_NORM_1+MATRICES_M @   ī  b       PROD_DIAG_TRIDIAG_MAT%COMPUTE_NORM_1+MATRICES_M E   Q  t   a   PROD_DIAG_TRIDIAG_MAT%COMPUTE_NORM_1%THIS+MATRICES_M +   Å  P   a   SQ_MATRIX_C%DIM+MATRICES_M 3        a   SQ_MATRIX_C%EIGENVALUES+MATRICES_M 4   ą  ´   a   SQ_MATRIX_C%EIGENVECTORS+MATRICES_M ;   e  i   a   SQ_MATRIX_C%COMPUTE_EIGENVALUES+MATRICES_M E   Î  Z       PROD_DIAG_TRIDIAG_MAT%COMPUTE_EIGENVALUES+MATRICES_M J   (  w   a   PROD_DIAG_TRIDIAG_MAT%COMPUTE_EIGENVALUES%THIS+MATRICES_M <     j   a   SQ_MATRIX_C%COMPUTE_EIGENVECTORS+MATRICES_M F   	  Z       PROD_DIAG_TRIDIAG_MAT%COMPUTE_EIGENVECTORS+MATRICES_M K   c  w   a   PROD_DIAG_TRIDIAG_MAT%COMPUTE_EIGENVECTORS%THIS+MATRICES_M .   Ú     a   DIAG_MATRIX_C%DIAG+MATRICES_M 9   v  e   a   DIAG_MATRIX_C%SET_DIAG_MATRIX+MATRICES_M A   Û  d       PROD_DIAG_TRIDIAG_MAT%SET_DIAG_MATRIX+MATRICES_M F   ?  y   a   PROD_DIAG_TRIDIAG_MAT%SET_DIAG_MATRIX%THIS+MATRICES_M F   ¸     a   PROD_DIAG_TRIDIAG_MAT%SET_DIAG_MATRIX%DIAG+MATRICES_M 4   L     a   TRIDIAG_SYM_MATRIX_C%SUB+MATRICES_M V   č  {   a   TRIDIAG_SYM_MATRIX_C%CHECK_EIGENVECTORS_TRIDIAG_SYM_MATRIX+MATRICES_M W   c  i       PROD_DIAG_TRIDIAG_MAT%CHECK_EIGENVECTORS_TRIDIAG_SYM_MATRIX+MATRICES_M \   Ė     a   PROD_DIAG_TRIDIAG_MAT%CHECK_EIGENVECTORS_TRIDIAG_SYM_MATRIX%THIS+MATRICES_M a   L  H   a   PROD_DIAG_TRIDIAG_MAT%CHECK_EIGENVECTORS_TRIDIAG_SYM_MATRIX%TOLERANCE+MATRICES_M 2        a   TRIDIAG_MATRIX_C%SUPER+MATRICES_M ?   0  h   a   TRIDIAG_MATRIX_C%SET_TRIDIAG_MATRIX+MATRICES_M D     x       PROD_DIAG_TRIDIAG_MAT%SET_TRIDIAG_MATRIX+MATRICES_M I     |   a   PROD_DIAG_TRIDIAG_MAT%SET_TRIDIAG_MATRIX%THIS+MATRICES_M H        a   PROD_DIAG_TRIDIAG_MAT%SET_TRIDIAG_MATRIX%SUB+MATRICES_M I         a   PROD_DIAG_TRIDIAG_MAT%SET_TRIDIAG_MATRIX%DIAG+MATRICES_M J   ´     a   PROD_DIAG_TRIDIAG_MAT%SET_TRIDIAG_MATRIX%SUPER+MATRICES_M M   H  v   a   TRIDIAG_MATRIX_C%COMPUTE_TRANSPOSE_TRIDIAG_MATRIX+MATRICES_M R   ž  i       PROD_DIAG_TRIDIAG_MAT%COMPUTE_TRANSPOSE_TRIDIAG_MATRIX+MATRICES_M W   '  |   a   PROD_DIAG_TRIDIAG_MAT%COMPUTE_TRANSPOSE_TRIDIAG_MATRIX%THIS+MATRICES_M \   Ŗ  |   a   PROD_DIAG_TRIDIAG_MAT%COMPUTE_TRANSPOSE_TRIDIAG_MATRIX%TRANSPOSE+MATRICES_M K      t   a   TRIDIAG_MATRIX_C%COMPUTE_INVERSE_TRIDIAG_MATRIX+MATRICES_M P      p       PROD_DIAG_TRIDIAG_MAT%COMPUTE_INVERSE_TRIDIAG_MATRIX+MATRICES_M U   !  |   a   PROD_DIAG_TRIDIAG_MAT%COMPUTE_INVERSE_TRIDIAG_MATRIX%THIS+MATRICES_M T   !  H   a   PROD_DIAG_TRIDIAG_MAT%COMPUTE_INVERSE_TRIDIAG_MATRIX%TOL+MATRICES_M X   Į!  Ŧ   a   PROD_DIAG_TRIDIAG_MAT%COMPUTE_INVERSE_TRIDIAG_MATRIX%INV_MAT+MATRICES_M A   s"  j   a   TRIDIAG_MATRIX_C%PROD_TRIDIAG_MAT_MAT+MATRICES_M F   Ũ"  Ņ       PROD_DIAG_TRIDIAG_MAT%PROD_TRIDIAG_MAT_MAT+MATRICES_M K   Ž#  |   a   PROD_DIAG_TRIDIAG_MAT%PROD_TRIDIAG_MAT_MAT%THIS+MATRICES_M L   *$  Ŧ   a   PROD_DIAG_TRIDIAG_MAT%PROD_TRIDIAG_MAT_MAT%B_MAT+MATRICES_M    Ö$  c   a   A    9%  f   a   B B   %  ų       PROD_DIAG_TRIDIAG_MAT%TRIDIAG_MATRIX_C+MATRICES_M W   &     a   PROD_DIAG_TRIDIAG_MAT%TRIDIAG_MATRIX_C%TRIDIAG_SYM_MATRIX_C+MATRICES_M F    '         PROD_DIAG_TRIDIAG_MAT%TRIDIAG_SYM_MATRIX_C+MATRICES_M T   ŋ'     a   PROD_DIAG_TRIDIAG_MAT%TRIDIAG_SYM_MATRIX_C%DIAG_MATRIX_C+MATRICES_M ?   @(         PROD_DIAG_TRIDIAG_MAT%DIAG_MATRIX_C+MATRICES_M K   Č(     a   PROD_DIAG_TRIDIAG_MAT%DIAG_MATRIX_C%SQ_MATRIX_C+MATRICES_M =   G)  Å       PROD_DIAG_TRIDIAG_MAT%SQ_MATRIX_C+MATRICES_M F   *  |   a   PROD_DIAG_TRIDIAG_MAT%SQ_MATRIX_C%MATRIX_C+MATRICES_M :   *  Ķ       PROD_DIAG_TRIDIAG_MAT%MATRIX_C+MATRICES_M J   [+  e   a   PROD_DIAG_TRIDIAG_MAT%MATRIX_C%ALLOCATE_MATRIX+MATRICES_M G   Ā+  b   a   PROD_DIAG_TRIDIAG_MAT%MATRIX_C%PROD_MAT_VEC+MATRICES_M C   ",  ^   a   PROD_DIAG_TRIDIAG_MAT%MATRIX_C%GET_DIAG+MATRICES_M B   ,  ]   a   PROD_DIAG_TRIDIAG_MAT%MATRIX_C%GET_SUB+MATRICES_M D   Ũ,  _   a   PROD_DIAG_TRIDIAG_MAT%MATRIX_C%GET_SUPER+MATRICES_M K   <-  f   a   PROD_DIAG_TRIDIAG_MAT%MATRIX_C%COMPUTE_NORM_INF+MATRICES_M I   ĸ-  d   a   PROD_DIAG_TRIDIAG_MAT%MATRIX_C%COMPUTE_NORM_1+MATRICES_M A   .  P   a   PROD_DIAG_TRIDIAG_MAT%SQ_MATRIX_C%DIM+MATRICES_M I   V.     a   PROD_DIAG_TRIDIAG_MAT%SQ_MATRIX_C%EIGENVALUES+MATRICES_M J   ō.  ´   a   PROD_DIAG_TRIDIAG_MAT%SQ_MATRIX_C%EIGENVECTORS+MATRICES_M Q   Ļ/  i   a   PROD_DIAG_TRIDIAG_MAT%SQ_MATRIX_C%COMPUTE_EIGENVALUES+MATRICES_M R   0  j   a   PROD_DIAG_TRIDIAG_MAT%SQ_MATRIX_C%COMPUTE_EIGENVECTORS+MATRICES_M D   y0     a   PROD_DIAG_TRIDIAG_MAT%DIAG_MATRIX_C%DIAG+MATRICES_M O   1  e   a   PROD_DIAG_TRIDIAG_MAT%DIAG_MATRIX_C%SET_DIAG_MATRIX+MATRICES_M J   z1     a   PROD_DIAG_TRIDIAG_MAT%TRIDIAG_SYM_MATRIX_C%SUB+MATRICES_M l   2  {   a   PROD_DIAG_TRIDIAG_MAT%TRIDIAG_SYM_MATRIX_C%CHECK_EIGENVECTORS_TRIDIAG_SYM_MATRIX+MATRICES_M H   2     a   PROD_DIAG_TRIDIAG_MAT%TRIDIAG_MATRIX_C%SUPER+MATRICES_M U   -3  h   a   PROD_DIAG_TRIDIAG_MAT%TRIDIAG_MATRIX_C%SET_TRIDIAG_MATRIX+MATRICES_M c   3  v   a   PROD_DIAG_TRIDIAG_MAT%TRIDIAG_MATRIX_C%COMPUTE_TRANSPOSE_TRIDIAG_MATRIX+MATRICES_M a   4  t   a   PROD_DIAG_TRIDIAG_MAT%TRIDIAG_MATRIX_C%COMPUTE_INVERSE_TRIDIAG_MATRIX+MATRICES_M W   4  j   a   PROD_DIAG_TRIDIAG_MAT%TRIDIAG_MATRIX_C%PROD_TRIDIAG_MAT_MAT+MATRICES_M 
  ń3  |   k820309    Ö          24.0        öŽe                                                                                                          
       C:\Users\Jordi\source\repos\jordipg10\1D-Transport-Code\diff_props_heterog_m.f90 DIFF_PROPS_HETEROG_M                                                             
                          @                                       '                     #NUM_TARGETS    #NUM_TARGETS_DEFINED    #TARGETS_FLAG    #MEASURE    #SCHEME    #ADAPT_REF    #SET_TARGETS 	   #SET_MEASURE    #SET_SCHEME    #READ_MESH    #GET_MESH_SIZE    #GET_DIM    #COMPUTE_MEASURE !   #REFINE_MESH $                                                                                                                                                                                                                                                                                          
                                                                                                                                               1         À    $                                     	                  #SET_TARGETS 
   #         @                                            
                    #THIS    #NUM_TARGETS    #FLAG                                                                            #SPATIAL_DISCR_C              
                                                               
                                                     1         À    $                                                       #SET_MEASURE    #         @                                                                #THIS    #MEASURE                                                                            #SPATIAL_DISCR_C              
                                               
      1         À    $                                                  	     #SET_SCHEME    #         @                                                                #THIS    #SCHEME                                                                            #SPATIAL_DISCR_C              
                                                     1         À    $                                                 
     #READ_MESH    #         @                                                	               #THIS    #FILENAME                                                                           #SPATIAL_DISCR_C              
                                                            1 1         À    $                                                     #GET_MESH_SIZE    %         @                                                             
       #THIS    #I                                                                           #SPATIAL_DISCR_C              
                                                   1         À    $                                                     #GET_DIM    %         @                                                                    #THIS                                                                             #SPATIAL_DISCR_C    1         À    $                                    !                  #COMPUTE_MEASURE "   #         @                                           "     	               #THIS #                                                     #                     #SPATIAL_DISCR_C    1         À    $                                    $                  #REFINE_MESH %   #         @                                           %     	               #THIS &   #CONC '   #CONC_EXT (   #REL_TOL )                                                     &                     #SPATIAL_DISCR_C            
                                        '                   
               &                   &                                                   
                                        (                   
               &                   &                                                     
                                         )     
      #         @                                           *     	               #THIS +   #FILENAME -   #SPATIAL_DISCR .                                                     +                     #PROPS_C ,             
                                        -                    1           
                                        .                    #SPATIAL_DISCR_C    #         @                                           /     	               #THIS 0                                                     0                     #PROPS_C ,   #         @                                            1                    #THIS 2   #SOURCE_TERM 3                                                      2                     #PROPS_C ,             
                                          3                   
              &                                           #         @                                            4                    #THIS 5   #SOURCE_TERM_ORDER 6                                                      5                     #PROPS_C ,             
                                          6           #         @                                            7                    #THIS 8   #BCS 9                                                      8                     #PROPS_C ,             
                                          9     0              #BCS_T :   #         @                                            ;                    #THIS <   #BCS =                                                      <     0               #BCS_T :             
                                          =                       p          p            p                          #         @                                            >                    #THIS ?   #EVAP @                                                      ?     0               #BCS_T :             
                                           @           #         @                                            A                    #THIS B   #FILENAME C                                                      B     0               #BCS_T :             
                                         C                    1 #         @                                            D                    #THIS E   #FILENAME F                                                      E     0               #BCS_T :             
                                         F                    1 #         @                                            G                    #THIS H   #FILENAME I                                                      H     0               #BCS_T :             
                                         I                    1 #         @                                            J                    #THIS K   #FILENAME L                                                      K     0               #BCS_T :             
                                         L                    1 #         @                                            M                    #THIS N   #CONC_INF O   #CONC_OUT P                                                      N     0               #BCS_T :             
                                          O     
                
                                          P     
      #         @                                            Q                    #THIS R   #FLUX S                                                      R     0               #BCS_T :             
                                          S     
                         @               À                  T     '0                   #PROPS_C U   #POROSITY _   #DISPERSION `   #SET_PROPS_DIFF_HETEROG a   #READ_PROPS f   #ARE_PROPS_HOMOG k                 $                                       U                            #PROPS_C ,                      @               @                  ,     '               	      #SOURCE_TERM_ORDER V   #SOURCE_TERM W   #SOURCE_TERM_FLAG X   #HOMOG_FLAG Y   #SET_SOURCE_TERM Z   #SET_SOURCE_TERM_ORDER [   #SET_SOURCE_TERM_FLAG \   #READ_PROPS ]   #ARE_PROPS_HOMOG ^                                                       V                                                                    W                             
            &                                                                                             X            P                             &                                                                                                Y                  1         À    $                                     Z                  #SET_SOURCE_TERM 1   1         À    $                                     [                  #SET_SOURCE_TERM_ORDER 4   1         À    $                                     \                  #SET_SOURCE_TERM_FLAG 7   1         À    $                                    ]                  #READ_PROPS *   1         À    $                                    ^             	     #ARE_PROPS_HOMOG /                                                    _                              
            &                                                                                             `            è                 
            &                                           1         À    $                                     a                  #SET_PROPS_DIFF_HETEROG b   #         @                                            b                    #THIS c   #POROSITY d   #DISPERSION e             D                                         c     0              #DIFF_PROPS_HETEROG_C T             
 @                                       d                   
              &                                                     
 @                                       e                   
              &                                           1         À    $                                    f                  #READ_PROPS_DIFF_HETEROG g   #         @                                            g                    #THIS h   #FILENAME i   #SPATIAL_DISCR j             D                                         h     0              #DIFF_PROPS_HETEROG_C T             
                                         i                    1           
                                         j                    #SPATIAL_DISCR_C    1         À    $                                    k                  #ARE_DIFF_PROPS_HOMOG l   #         @                                            l                    #THIS m             D                                         m     0              #DIFF_PROPS_HETEROG_C T                      @                                   :     '0                    #BCS_LABEL n   #EVAP o   #CONC_INF p   #CONC_OUT q   #FLUX_INF r   #FLUX_OUT s   #SET_BCS_LABEL t   #SET_EVAP u   #READ_BCS v   #READ_DIRICHLET_BCS w   #READ_ROBIN_BC_INFLOW x   #READ_FLUX_INF y   #SET_CONC_BOUNDARY z   #SET_CST_FLUX_BOUNDARY {                                                       n                                p          p            p                                                                               o                                                                      p               
                                                       q               
                                                       r                
                                                       s     (          
   1         À    $                                     t                  #SET_BCS_LABEL ;   1         À    $                                     u                  #SET_EVAP >   1         À    $                                     v             	     #READ_BCS A   1         À    $                                     w             
     #READ_DIRICHLET_BCS D   1         À    $                                     x                  #READ_ROBIN_BC_INFLOW G   1         À    $                                     y                  #READ_FLUX_INF J   1         À    $                                     z                  #SET_CONC_BOUNDARY M   1         À    $                                     {                  #SET_CST_FLUX_BOUNDARY Q          n      fn#fn      H   J   PROPERTIES_M 0   V  C      SPATIAL_DISCR_C+SPATIAL_DISCR_M <     P   a   SPATIAL_DISCR_C%NUM_TARGETS+SPATIAL_DISCR_M D   é  P   a   SPATIAL_DISCR_C%NUM_TARGETS_DEFINED+SPATIAL_DISCR_M =   9  P   a   SPATIAL_DISCR_C%TARGETS_FLAG+SPATIAL_DISCR_M 8     P   a   SPATIAL_DISCR_C%MEASURE+SPATIAL_DISCR_M 7   Ù  P   a   SPATIAL_DISCR_C%SCHEME+SPATIAL_DISCR_M :   )  P   a   SPATIAL_DISCR_C%ADAPT_REF+SPATIAL_DISCR_M <   y  a   a   SPATIAL_DISCR_C%SET_TARGETS+SPATIAL_DISCR_M ,   Ú  u       SET_TARGETS+SPATIAL_DISCR_M 1   O  e   a   SET_TARGETS%THIS+SPATIAL_DISCR_M 8   Ž  H   a   SET_TARGETS%NUM_TARGETS+SPATIAL_DISCR_M 1   ü  H   a   SET_TARGETS%FLAG+SPATIAL_DISCR_M <   D  a   a   SPATIAL_DISCR_C%SET_MEASURE+SPATIAL_DISCR_M ,   „  g       SET_MEASURE+SPATIAL_DISCR_M 1     e   a   SET_MEASURE%THIS+SPATIAL_DISCR_M 4   q  H   a   SET_MEASURE%MEASURE+SPATIAL_DISCR_M ;   č  `   a   SPATIAL_DISCR_C%SET_SCHEME+SPATIAL_DISCR_M +     f       SET_SCHEME+SPATIAL_DISCR_M 0     e   a   SET_SCHEME%THIS+SPATIAL_DISCR_M 2   ä  H   a   SET_SCHEME%SCHEME+SPATIAL_DISCR_M :   ,	  _   a   SPATIAL_DISCR_C%READ_MESH+SPATIAL_DISCR_M *   	  h       READ_MESH+SPATIAL_DISCR_M /   ó	  e   a   READ_MESH%THIS+SPATIAL_DISCR_M 3   X
  T   a   READ_MESH%FILENAME+SPATIAL_DISCR_M >   Ź
  c   a   SPATIAL_DISCR_C%GET_MESH_SIZE+SPATIAL_DISCR_M .     i       GET_MESH_SIZE+SPATIAL_DISCR_M 3   x  e   a   GET_MESH_SIZE%THIS+SPATIAL_DISCR_M 0   Ę  H   a   GET_MESH_SIZE%I+SPATIAL_DISCR_M 8   %  ]   a   SPATIAL_DISCR_C%GET_DIM+SPATIAL_DISCR_M (     b       GET_DIM+SPATIAL_DISCR_M -   ä  e   a   GET_DIM%THIS+SPATIAL_DISCR_M @   I  e   a   SPATIAL_DISCR_C%COMPUTE_MEASURE+SPATIAL_DISCR_M 0   ź  Z       COMPUTE_MEASURE+SPATIAL_DISCR_M 5     e   a   COMPUTE_MEASURE%THIS+SPATIAL_DISCR_M <   m  a   a   SPATIAL_DISCR_C%REFINE_MESH+SPATIAL_DISCR_M ,   Î         REFINE_MESH+SPATIAL_DISCR_M 1   M  e   a   REFINE_MESH%THIS+SPATIAL_DISCR_M 1   Č  Ź   a   REFINE_MESH%CONC+SPATIAL_DISCR_M 5   ^  Ź   a   REFINE_MESH%CONC_EXT+SPATIAL_DISCR_M 4   
  H   a   REFINE_MESH%REL_TOL+SPATIAL_DISCR_M (   R  {       READ_PROPS+PROPERTIES_M -   Í  ]   a   READ_PROPS%THIS+PROPERTIES_M 1   *  T   a   READ_PROPS%FILENAME+PROPERTIES_M 6   ~  e   a   READ_PROPS%SPATIAL_DISCR+PROPERTIES_M -   ă  Z       ARE_PROPS_HOMOG+PROPERTIES_M 2   =  ]   a   ARE_PROPS_HOMOG%THIS+PROPERTIES_M -     k       SET_SOURCE_TERM+PROPERTIES_M 2     ]   a   SET_SOURCE_TERM%THIS+PROPERTIES_M 9   b     a   SET_SOURCE_TERM%SOURCE_TERM+PROPERTIES_M 3   ö  q       SET_SOURCE_TERM_ORDER+PROPERTIES_M 8   g  ]   a   SET_SOURCE_TERM_ORDER%THIS+PROPERTIES_M E   Ä  H   a   SET_SOURCE_TERM_ORDER%SOURCE_TERM_ORDER+PROPERTIES_M 2     c       SET_SOURCE_TERM_FLAG+PROPERTIES_M 7   o  ]   a   SET_SOURCE_TERM_FLAG%THIS+PROPERTIES_M 6   Ì  [   a   SET_SOURCE_TERM_FLAG%BCS+PROPERTIES_M $   '  c       SET_BCS_LABEL+BCS_M )     [   a   SET_BCS_LABEL%THIS+BCS_M (   ć     a   SET_BCS_LABEL%BCS+BCS_M      d       SET_EVAP+BCS_M $   ć  [   a   SET_EVAP%THIS+BCS_M $   @  H   a   SET_EVAP%EVAP+BCS_M      h       READ_BCS+BCS_M $   đ  [   a   READ_BCS%THIS+BCS_M (   K  T   a   READ_BCS%FILENAME+BCS_M )     h       READ_DIRICHLET_BCS+BCS_M .     [   a   READ_DIRICHLET_BCS%THIS+BCS_M 2   b  T   a   READ_DIRICHLET_BCS%FILENAME+BCS_M +   ¶  h       READ_ROBIN_BC_INFLOW+BCS_M 0     [   a   READ_ROBIN_BC_INFLOW%THIS+BCS_M 4   y  T   a   READ_ROBIN_BC_INFLOW%FILENAME+BCS_M $   Í  h       READ_FLUX_INF+BCS_M )   5  [   a   READ_FLUX_INF%THIS+BCS_M -     T   a   READ_FLUX_INF%FILENAME+BCS_M (   ä  v       SET_CONC_BOUNDARY+BCS_M -   Z  [   a   SET_CONC_BOUNDARY%THIS+BCS_M 1   ”  H   a   SET_CONC_BOUNDARY%CONC_INF+BCS_M 1   ę  H   a   SET_CONC_BOUNDARY%CONC_OUT+BCS_M ,   E  d       SET_CST_FLUX_BOUNDARY+BCS_M 1   ©  [   a   SET_CST_FLUX_BOUNDARY%THIS+BCS_M 1      H   a   SET_CST_FLUX_BOUNDARY%FLUX+BCS_M %   L   Ä       DIFF_PROPS_HETEROG_C -   !  e   a   DIFF_PROPS_HETEROG_C%PROPS_C %   u!        PROPS_C+PROPERTIES_M 7   "  P   a   PROPS_C%SOURCE_TERM_ORDER+PROPERTIES_M 1   Ú"     a   PROPS_C%SOURCE_TERM+PROPERTIES_M 6   v#     a   PROPS_C%SOURCE_TERM_FLAG+PROPERTIES_M 0   $  P   a   PROPS_C%HOMOG_FLAG+PROPERTIES_M 5   b$  e   a   PROPS_C%SET_SOURCE_TERM+PROPERTIES_M ;   Ç$  k   a   PROPS_C%SET_SOURCE_TERM_ORDER+PROPERTIES_M :   2%  j   a   PROPS_C%SET_SOURCE_TERM_FLAG+PROPERTIES_M 0   %  `   a   PROPS_C%READ_PROPS+PROPERTIES_M 5   ü%  e   a   PROPS_C%ARE_PROPS_HOMOG+PROPERTIES_M .   a&     a   DIFF_PROPS_HETEROG_C%POROSITY 0   ę&     a   DIFF_PROPS_HETEROG_C%DISPERSION <   '  l   a   DIFF_PROPS_HETEROG_C%SET_PROPS_DIFF_HETEROG '   (  x       SET_PROPS_DIFF_HETEROG ,   }(  j   a   SET_PROPS_DIFF_HETEROG%THIS 0   ç(     a   SET_PROPS_DIFF_HETEROG%POROSITY 2   {)     a   SET_PROPS_DIFF_HETEROG%DISPERSION 0   *  m   a   DIFF_PROPS_HETEROG_C%READ_PROPS (   |*  {       READ_PROPS_DIFF_HETEROG -   ś*  j   a   READ_PROPS_DIFF_HETEROG%THIS 1   a+  T   a   READ_PROPS_DIFF_HETEROG%FILENAME 6   ”+  e   a   READ_PROPS_DIFF_HETEROG%SPATIAL_DISCR 5   ,  j   a   DIFF_PROPS_HETEROG_C%ARE_PROPS_HOMOG %   ,  Z       ARE_DIFF_PROPS_HOMOG *   Ț,  j   a   ARE_DIFF_PROPS_HOMOG%THIS    H-  O      BCS_T+BCS_M &   .  €   a   BCS_T%BCS_LABEL+BCS_M !   ;/  P   a   BCS_T%EVAP+BCS_M %   /  P   a   BCS_T%CONC_INF+BCS_M %   Û/  P   a   BCS_T%CONC_OUT+BCS_M %   +0  P   a   BCS_T%FLUX_INF+BCS_M %   {0  P   a   BCS_T%FLUX_OUT+BCS_M *   Ë0  c   a   BCS_T%SET_BCS_LABEL+BCS_M %   .1  ^   a   BCS_T%SET_EVAP+BCS_M %   1  ^   a   BCS_T%READ_BCS+BCS_M /   ê1  h   a   BCS_T%READ_DIRICHLET_BCS+BCS_M 1   R2  j   a   BCS_T%READ_ROBIN_BC_INFLOW+BCS_M *   Œ2  c   a   BCS_T%READ_FLUX_INF+BCS_M .   3  g   a   BCS_T%SET_CONC_BOUNDARY+BCS_M 2   3  k   a   BCS_T%SET_CST_FLUX_BOUNDARY+BCS_M 
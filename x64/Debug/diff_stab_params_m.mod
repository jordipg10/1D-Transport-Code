  ģ;     k820309    Ö          24.0        CÏĢe                                                                                                          
       C:\Users\Jordi\source\repos\jordipg10\1D-Transport-Code\diff_stab_params_m.f90 DIFF_STAB_PARAMS_M                                                             
                                                                   
                                                                   
                                                                   
                                                                   
                          @               @                       '               	      #SOURCE_TERM_ORDER    #SOURCE_TERM    #SOURCE_TERM_FLAG 	   #HOMOG_FLAG 
   #SET_SOURCE_TERM    #SET_SOURCE_TERM_ORDER    #SET_SOURCE_TERM_FLAG    #READ_PROPS    #ARE_PROPS_HOMOG                                                                                                                                                         
            &                                                                                             	            P                             &                                                                                                
                  1         Ā    $                                                       #SET_SOURCE_TERM    #         @                                                                #THIS    #SOURCE_TERM                                                                            #PROPS_C              
                                                             
              &                                           1         Ā    $                                                       #SET_SOURCE_TERM_ORDER    #         @                                                                #THIS    #SOURCE_TERM_ORDER                                                                            #PROPS_C              
                                                     1         Ā    $                                                       #SET_SOURCE_TERM_FLAG    #         @                                                                #THIS    #BCS                                                                            #PROPS_C              
                                               0              #BCS_T    1         Ā    $                                                      #READ_PROPS    #         @                                                	               #THIS    #FILENAME    #SPATIAL_DISCR                                                                           #PROPS_C              
                                                            1           
                                                            #SPATIAL_DISCR_C    1         Ā    $                                                 	     #ARE_PROPS_HOMOG    #         @                                                	               #THIS                                                                             #PROPS_C                       @               Ā                  !     '0                   #PROPS_C "   #POROSITY #   #DISPERSION $   #SET_PROPS_DIFF_HETEROG %   #READ_PROPS *   #ARE_PROPS_HOMOG /                 $                                       "                            #PROPS_C                                                     #                              
            &                                                                                             $            č                 
            &                                           1         Ā    $                                     %                  #SET_PROPS_DIFF_HETEROG &   #         @                                            &                    #THIS '   #POROSITY (   #DISPERSION )                                                      '     0              #DIFF_PROPS_HETEROG_C !             
 @                                       (                   
              &                                                     
 @                                       )                   
              &                                           1         Ā    $                                    *                  #READ_PROPS_DIFF_HETEROG +   #         @                                            +                    #THIS ,   #FILENAME -   #SPATIAL_DISCR .                                                      ,     0              #DIFF_PROPS_HETEROG_C !             
                                         -                    1           
                                         .                    #SPATIAL_DISCR_C    1         Ā    $                                    /                  #ARE_DIFF_PROPS_HOMOG 0   #         @                                            0                    #THIS 1                                                      1     0              #DIFF_PROPS_HETEROG_C !   #         @                                           2     	               #THIS 3   #PROPS_OBJ 5   #MESH_SIZE 6   #TIME_STEP 7                                                     3                    #STAB_PARAMS_C 4             
                                         5                    #PROPS_C              
                                         6     
                
                                         7     
      #         @                                            8                    #THIS 9   #BCS :                                                      9     0               #BCS_T              
                                          :                       p          p            p                          #         @                                            ;                    #THIS <   #EVAP =                                                      <     0               #BCS_T              
                                           =           #         @                                            >                    #THIS ?   #FILENAME @                                                      ?     0               #BCS_T              
                                         @                    1 #         @                                            A                    #THIS B   #FILENAME C                                                      B     0               #BCS_T              
                                         C                    1 #         @                                            D                    #THIS E   #FILENAME F                                                      E     0               #BCS_T              
                                         F                    1 #         @                                            G                    #THIS H   #FILENAME I                                                      H     0               #BCS_T              
                                         I                    1 #         @                                            J                    #THIS K   #CONC_INF L   #CONC_OUT M                                                      K     0               #BCS_T              
                                          L     
                
                                          M     
      #         @                                            N                    #THIS O   #FLUX P                                                      O     0               #BCS_T              
                                          P     
      #         @                                            Q                    #THIS R   #NUM_TARGETS S   #FLAG T                                                      R                     #SPATIAL_DISCR_C              
                                          S                     
                                          T           #         @                                            U                    #THIS V   #MEASURE W                                                      V                     #SPATIAL_DISCR_C              
                                          W     
      #         @                                            X                    #THIS Y   #SCHEME Z                                                      Y                     #SPATIAL_DISCR_C              
                                          Z           #         @                                           [     	               #THIS \   #FILENAME ]                                                     \                     #SPATIAL_DISCR_C              
                                        ]                    1 %         @                                         ^                    
       #THIS _   #I `                                                     _                     #SPATIAL_DISCR_C              
                                        `           %         @                                         a                           #THIS b                                                     b                     #SPATIAL_DISCR_C    #         @                                           c     	               #THIS d                                                     d                     #SPATIAL_DISCR_C    #         @                                           e     	               #THIS f   #CONC g   #CONC_EXT h   #REL_TOL i                                                     f                     #SPATIAL_DISCR_C            
                                        g                   
               &                   &                                                   
                                        h                   
               &                   &                                                     
                                         i     
                         @                                  j     '                    #STAB_PARAMS_C k   #BETA n   #COMPUTE_STAB_PARAMS o                 $                                       k                           #STAB_PARAMS_C 4                      @                                  4     '                    #DELTA_T_CRIT l   #COMPUTE_STAB_PARAMS m                                                       l                
   1         Ā    $                                    m                  #COMPUTE_STAB_PARAMS 2                                                       n               
   1         Ā    $                                    o                  #COMPUTE_STAB_PARAMS_DIFF p   #         @                                            p                    #THIS q   #PROPS_OBJ r   #MESH_SIZE s   #TIME_STEP t             D                                         q                    #STAB_PARAMS_DIFF_C j             
                                          r                    #PROPS_C              
                                          s     
                
                                          t     
                         @                                       '                     #NUM_TARGETS u   #NUM_TARGETS_DEFINED v   #TARGETS_FLAG w   #MEASURE x   #SCHEME y   #ADAPT_REF z   #SET_TARGETS {   #SET_MEASURE |   #SET_SCHEME }   #READ_MESH ~   #GET_MESH_SIZE    #GET_DIM    #COMPUTE_MEASURE    #REFINE_MESH                                                        u                                                                        v                                                                      w                                                                      x               
                                                       y                                                                      z                  1         Ā    $                                     {                  #SET_TARGETS Q   1         Ā    $                                     |                  #SET_MEASURE U   1         Ā    $                                     }             	     #SET_SCHEME X   1         Ā    $                                    ~             
     #READ_MESH [   1         Ā    $                                                     #GET_MESH_SIZE ^   1         Ā    $                                                     #GET_DIM a   1         Ā    $                                                      #COMPUTE_MEASURE c   1         Ā    $                                                      #REFINE_MESH e                      @                                        '0                    #BCS_LABEL    #EVAP    #CONC_INF    #CONC_OUT    #FLUX_INF    #FLUX_OUT    #SET_BCS_LABEL    #SET_EVAP    #READ_BCS    #READ_DIRICHLET_BCS    #READ_ROBIN_BC_INFLOW    #READ_FLUX_INF    #SET_CONC_BOUNDARY    #SET_CST_FLUX_BOUNDARY                                                                                        p          p            p                                                                                                                                                                    
                                                                      
                                                                       
                                                            (          
   1         Ā    $                                                       #SET_BCS_LABEL 8   1         Ā    $                                                       #SET_EVAP ;   1         Ā    $                                                  	     #READ_BCS >   1         Ā    $                                                  
     #READ_DIRICHLET_BCS A   1         Ā    $                                                       #READ_ROBIN_BC_INFLOW D   1         Ā    $                                                       #READ_FLUX_INF G   1         Ā    $                                                       #SET_CONC_BOUNDARY J   1         Ā    $                                                       #SET_CST_FLUX_BOUNDARY N          j      fn#fn '   
  H   J   STABILITY_PARAMETERS_M %   R  H   J   DIFF_PROPS_HETEROG_M $     H   J   SPATIAL_DISCR_RAD_M #   â  H   J   SPATIAL_DISCR_1D_M    *  H   J   VECTORS_M %   r        PROPS_C+PROPERTIES_M 7     P   a   PROPS_C%SOURCE_TERM_ORDER+PROPERTIES_M 1   Ũ     a   PROPS_C%SOURCE_TERM+PROPERTIES_M 6   s     a   PROPS_C%SOURCE_TERM_FLAG+PROPERTIES_M 0     P   a   PROPS_C%HOMOG_FLAG+PROPERTIES_M 5   _  e   a   PROPS_C%SET_SOURCE_TERM+PROPERTIES_M -   Ä  k       SET_SOURCE_TERM+PROPERTIES_M 2   /  ]   a   SET_SOURCE_TERM%THIS+PROPERTIES_M 9        a   SET_SOURCE_TERM%SOURCE_TERM+PROPERTIES_M ;      k   a   PROPS_C%SET_SOURCE_TERM_ORDER+PROPERTIES_M 3     q       SET_SOURCE_TERM_ORDER+PROPERTIES_M 8   ü  ]   a   SET_SOURCE_TERM_ORDER%THIS+PROPERTIES_M E   Y  H   a   SET_SOURCE_TERM_ORDER%SOURCE_TERM_ORDER+PROPERTIES_M :   Ą  j   a   PROPS_C%SET_SOURCE_TERM_FLAG+PROPERTIES_M 2   	  c       SET_SOURCE_TERM_FLAG+PROPERTIES_M 7   n	  ]   a   SET_SOURCE_TERM_FLAG%THIS+PROPERTIES_M 6   Ë	  [   a   SET_SOURCE_TERM_FLAG%BCS+PROPERTIES_M 0   &
  `   a   PROPS_C%READ_PROPS+PROPERTIES_M (   
  {       READ_PROPS+PROPERTIES_M -     ]   a   READ_PROPS%THIS+PROPERTIES_M 1   ^  T   a   READ_PROPS%FILENAME+PROPERTIES_M 6   ē  e   a   READ_PROPS%SPATIAL_DISCR+PROPERTIES_M 5     e   a   PROPS_C%ARE_PROPS_HOMOG+PROPERTIES_M -   |  Z       ARE_PROPS_HOMOG+PROPERTIES_M 2   Ö  ]   a   ARE_PROPS_HOMOG%THIS+PROPERTIES_M :   3  Ä       DIFF_PROPS_HETEROG_C+DIFF_PROPS_HETEROG_M B   ũ  e   a   DIFF_PROPS_HETEROG_C%PROPS_C+DIFF_PROPS_HETEROG_M C   \     a   DIFF_PROPS_HETEROG_C%POROSITY+DIFF_PROPS_HETEROG_M E   ø     a   DIFF_PROPS_HETEROG_C%DISPERSION+DIFF_PROPS_HETEROG_M Q     l   a   DIFF_PROPS_HETEROG_C%SET_PROPS_DIFF_HETEROG+DIFF_PROPS_HETEROG_M <      x       SET_PROPS_DIFF_HETEROG+DIFF_PROPS_HETEROG_M A   x  j   a   SET_PROPS_DIFF_HETEROG%THIS+DIFF_PROPS_HETEROG_M E   â     a   SET_PROPS_DIFF_HETEROG%POROSITY+DIFF_PROPS_HETEROG_M G   v     a   SET_PROPS_DIFF_HETEROG%DISPERSION+DIFF_PROPS_HETEROG_M E   
  m   a   DIFF_PROPS_HETEROG_C%READ_PROPS+DIFF_PROPS_HETEROG_M =   w  {       READ_PROPS_DIFF_HETEROG+DIFF_PROPS_HETEROG_M B   ō  j   a   READ_PROPS_DIFF_HETEROG%THIS+DIFF_PROPS_HETEROG_M F   \  T   a   READ_PROPS_DIFF_HETEROG%FILENAME+DIFF_PROPS_HETEROG_M K   °  e   a   READ_PROPS_DIFF_HETEROG%SPATIAL_DISCR+DIFF_PROPS_HETEROG_M J     j   a   DIFF_PROPS_HETEROG_C%ARE_PROPS_HOMOG+DIFF_PROPS_HETEROG_M :     Z       ARE_DIFF_PROPS_HOMOG+DIFF_PROPS_HETEROG_M ?   Ų  j   a   ARE_DIFF_PROPS_HOMOG%THIS+DIFF_PROPS_HETEROG_M ;   C         COMPUTE_STAB_PARAMS+STABILITY_PARAMETERS_M @   Ę  c   a   COMPUTE_STAB_PARAMS%THIS+STABILITY_PARAMETERS_M E   -  ]   a   COMPUTE_STAB_PARAMS%PROPS_OBJ+STABILITY_PARAMETERS_M E     H   a   COMPUTE_STAB_PARAMS%MESH_SIZE+STABILITY_PARAMETERS_M E   Ō  H   a   COMPUTE_STAB_PARAMS%TIME_STEP+STABILITY_PARAMETERS_M $     c       SET_BCS_LABEL+BCS_M )   }  [   a   SET_BCS_LABEL%THIS+BCS_M (   Ø     a   SET_BCS_LABEL%BCS+BCS_M    t  d       SET_EVAP+BCS_M $   Ø  [   a   SET_EVAP%THIS+BCS_M $   3  H   a   SET_EVAP%EVAP+BCS_M    {  h       READ_BCS+BCS_M $   ã  [   a   READ_BCS%THIS+BCS_M (   >  T   a   READ_BCS%FILENAME+BCS_M )     h       READ_DIRICHLET_BCS+BCS_M .   ú  [   a   READ_DIRICHLET_BCS%THIS+BCS_M 2   U  T   a   READ_DIRICHLET_BCS%FILENAME+BCS_M +   Đ  h       READ_ROBIN_BC_INFLOW+BCS_M 0     [   a   READ_ROBIN_BC_INFLOW%THIS+BCS_M 4   l  T   a   READ_ROBIN_BC_INFLOW%FILENAME+BCS_M $   Ā  h       READ_FLUX_INF+BCS_M )   (  [   a   READ_FLUX_INF%THIS+BCS_M -     T   a   READ_FLUX_INF%FILENAME+BCS_M (   Ũ  v       SET_CONC_BOUNDARY+BCS_M -   M  [   a   SET_CONC_BOUNDARY%THIS+BCS_M 1   Ļ  H   a   SET_CONC_BOUNDARY%CONC_INF+BCS_M 1   ð  H   a   SET_CONC_BOUNDARY%CONC_OUT+BCS_M ,   8  d       SET_CST_FLUX_BOUNDARY+BCS_M 1     [   a   SET_CST_FLUX_BOUNDARY%THIS+BCS_M 1   ũ  H   a   SET_CST_FLUX_BOUNDARY%FLUX+BCS_M ,   ?   u       SET_TARGETS+SPATIAL_DISCR_M 1   ī   e   a   SET_TARGETS%THIS+SPATIAL_DISCR_M 8   !  H   a   SET_TARGETS%NUM_TARGETS+SPATIAL_DISCR_M 1   a!  H   a   SET_TARGETS%FLAG+SPATIAL_DISCR_M ,   Đ!  g       SET_MEASURE+SPATIAL_DISCR_M 1   "  e   a   SET_MEASURE%THIS+SPATIAL_DISCR_M 4   u"  H   a   SET_MEASURE%MEASURE+SPATIAL_DISCR_M +   ―"  f       SET_SCHEME+SPATIAL_DISCR_M 0   ##  e   a   SET_SCHEME%THIS+SPATIAL_DISCR_M 2   #  H   a   SET_SCHEME%SCHEME+SPATIAL_DISCR_M *   Ð#  h       READ_MESH+SPATIAL_DISCR_M /   8$  e   a   READ_MESH%THIS+SPATIAL_DISCR_M 3   $  T   a   READ_MESH%FILENAME+SPATIAL_DISCR_M .   ņ$  i       GET_MESH_SIZE+SPATIAL_DISCR_M 3   Z%  e   a   GET_MESH_SIZE%THIS+SPATIAL_DISCR_M 0   ŋ%  H   a   GET_MESH_SIZE%I+SPATIAL_DISCR_M (   &  b       GET_DIM+SPATIAL_DISCR_M -   i&  e   a   GET_DIM%THIS+SPATIAL_DISCR_M 0   Î&  Z       COMPUTE_MEASURE+SPATIAL_DISCR_M 5   ('  e   a   COMPUTE_MEASURE%THIS+SPATIAL_DISCR_M ,   '         REFINE_MESH+SPATIAL_DISCR_M 1   (  e   a   REFINE_MESH%THIS+SPATIAL_DISCR_M 1   q(  Ž   a   REFINE_MESH%CONC+SPATIAL_DISCR_M 5   )  Ž   a   REFINE_MESH%CONC_EXT+SPATIAL_DISCR_M 4   É)  H   a   REFINE_MESH%REL_TOL+SPATIAL_DISCR_M #   *         STAB_PARAMS_DIFF_C 1   *  k   a   STAB_PARAMS_DIFF_C%STAB_PARAMS_C 5   
+         STAB_PARAMS_C+STABILITY_PARAMETERS_M B   +  P   a   STAB_PARAMS_C%DELTA_T_CRIT+STABILITY_PARAMETERS_M I   Ý+  i   a   STAB_PARAMS_C%COMPUTE_STAB_PARAMS+STABILITY_PARAMETERS_M (   F,  P   a   STAB_PARAMS_DIFF_C%BETA 7   ,  n   a   STAB_PARAMS_DIFF_C%COMPUTE_STAB_PARAMS )   -         COMPUTE_STAB_PARAMS_DIFF .   -  h   a   COMPUTE_STAB_PARAMS_DIFF%THIS 3   ó-  ]   a   COMPUTE_STAB_PARAMS_DIFF%PROPS_OBJ 3   P.  H   a   COMPUTE_STAB_PARAMS_DIFF%MESH_SIZE 3   .  H   a   COMPUTE_STAB_PARAMS_DIFF%TIME_STEP 0   ā.  C      SPATIAL_DISCR_C+SPATIAL_DISCR_M <   #0  P   a   SPATIAL_DISCR_C%NUM_TARGETS+SPATIAL_DISCR_M D   s0  P   a   SPATIAL_DISCR_C%NUM_TARGETS_DEFINED+SPATIAL_DISCR_M =   Ã0  P   a   SPATIAL_DISCR_C%TARGETS_FLAG+SPATIAL_DISCR_M 8   1  P   a   SPATIAL_DISCR_C%MEASURE+SPATIAL_DISCR_M 7   c1  P   a   SPATIAL_DISCR_C%SCHEME+SPATIAL_DISCR_M :   ģ1  P   a   SPATIAL_DISCR_C%ADAPT_REF+SPATIAL_DISCR_M <   2  a   a   SPATIAL_DISCR_C%SET_TARGETS+SPATIAL_DISCR_M <   d2  a   a   SPATIAL_DISCR_C%SET_MEASURE+SPATIAL_DISCR_M ;   Å2  `   a   SPATIAL_DISCR_C%SET_SCHEME+SPATIAL_DISCR_M :   %3  _   a   SPATIAL_DISCR_C%READ_MESH+SPATIAL_DISCR_M >   3  c   a   SPATIAL_DISCR_C%GET_MESH_SIZE+SPATIAL_DISCR_M 8   į3  ]   a   SPATIAL_DISCR_C%GET_DIM+SPATIAL_DISCR_M @   D4  e   a   SPATIAL_DISCR_C%COMPUTE_MEASURE+SPATIAL_DISCR_M <   Đ4  a   a   SPATIAL_DISCR_C%REFINE_MESH+SPATIAL_DISCR_M    
5  O      BCS_T+BCS_M &   Y6  Ī   a   BCS_T%BCS_LABEL+BCS_M !   ý6  P   a   BCS_T%EVAP+BCS_M %   M7  P   a   BCS_T%CONC_INF+BCS_M %   7  P   a   BCS_T%CONC_OUT+BCS_M %   í7  P   a   BCS_T%FLUX_INF+BCS_M %   =8  P   a   BCS_T%FLUX_OUT+BCS_M *   8  c   a   BCS_T%SET_BCS_LABEL+BCS_M %   ð8  ^   a   BCS_T%SET_EVAP+BCS_M %   N9  ^   a   BCS_T%READ_BCS+BCS_M /   Ž9  h   a   BCS_T%READ_DIRICHLET_BCS+BCS_M 1   :  j   a   BCS_T%READ_ROBIN_BC_INFLOW+BCS_M *   ~:  c   a   BCS_T%READ_FLUX_INF+BCS_M .   á:  g   a   BCS_T%SET_CONC_BOUNDARY+BCS_M 2   H;  k   a   BCS_T%SET_CST_FLUX_BOUNDARY+BCS_M 
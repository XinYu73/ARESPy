  ôT  ·   k820309    l          18.0        Ø3Lb                                                                                                          
       XC_functional.f90 XC_MODULE                   @                               
                @                                         
                                                           
                                                          
                @          @                              
                                                           
      IXCDF IXC NSPIN LCORE_VAL SNLCC          @                                         
      SPH GRID DVOL                  @                                '                    #BUFFER 	                D                              	                                     @               @           
     '@                   #COMM    #MYID    #NUMPROCS    #ROOTID    #ISROOT    #NSTATE_PROC    #SUB2SUM    #MYGRID_RANGE    #RECVCOUNTS    #DISPLS    #GLOBAL_GRIDRANGE    #COMM2D    #COMMX    #COMMY    #RANKX    #RANKY    #PERIODS    #REORDER    #REMAINX    #REMAINY    #NDIMS    #DIMS     #COMMFFT !   #LOCAL_Z "   #LOCAL_Z_START #   #FFT_GRID_RANGE $   #FFT_RCOUNT %   #FFT_RDISPLS &   #FFT_SCOUNT '   #FFT_SDISPLS (                                                                                                                                                                                                                                                                                                                                                                                                                                                                     &                   &                                                                                                  x                   p          p            p                                                                                            	               &                                                                                               Ð              
               &                                                                                                                           &                   &                                                                                           x                                                            |                                                                                                                                                                                                                                                                         p          p            p                                                                                                                                                               p          p            p                                                                                                    p          p            p                                                                          ¨                                                                    ¬                  p          p            p                                                                     !     ´                                                       "     ¸                                                       #     ¼                                                     $            À                            &                   &                                                                                    %                                         &                                                                                    &            h                            &                                                                                    '            °                            &                                                                                    (            ø                            &                                                             @               @           )     '@                   #NZ_MAP *   #MYCOMM_CORES +   #MYCOMM_SIZE ,   #MYSEND_SIZE -   #LOCAL_MAP .   #LOCAL_MAP1D /   #BOUNDARY 0   #BOUNDARY1D 1   #RCOUNT 2   #RDISPLS 3   #SCOUNT 4   #SDISPLS 5                                           *                                          &                   &                                                                                      +            `                   p          p            p                                                                  ,            h                             &                   &                                                                                   -            È                             &                   &                                                                                   .            (                            &                   &                                                                                   /                                        &                   &                   &                                                                                      0                               p          p          p            p          p                                                                     1                              p          p            p                                                                   2                          	               &                                                                                    3            h             
               &                                                                                    4            °                            &                                                                                    5            ø                            &                                                             @               D           6     '@                   #RHOS 7   #RHO 8   #VXCS 9   #VHXCD :   #VLPP ;   #VH <   #EVAL =   #RHOC >   #GVEC ?   #GMASK @   #RVEC A   #LSP B   #ONEDLENGTH C                                            7                              
            &                   &                                                                                    8            `                 
            &                                                                                    9            ¨                 
            &                   &                                                                                    :                            
            &                   &                                                                                    ;            h                
            &                                                                                    <            °                
            &                                                                                    =            ø                
            &                   &                   &                                                                                    >            p                
            &                                                                                    ?            ¸             	   
            &                   &                                                                                     @                         
               &                                                                                    A            `                
            &                   &                                                                                     B            À                            &                   &                   &                                                                                      C     8                       @ @                               D     @      #PARALLEL_TYPE 
              @                              E                                                       F                                                        G                                                       H     
                                                    I     @      #GRID_TYPE 6                                               J     
                                                    K                                                                                                      L                                                         #         @                                  M     	              #XC_F90_FUNC_INIT%XC_F90_POINTER_T N   #P P   #INFO Q   #FUNCTIONAL R   #NSPIN S                    @                           N     '                    #BUFFER O                D                             O                                                             P                    #XC_F90_FUNC_INIT%XC_F90_POINTER_T N                                             Q                    #XC_F90_FUNC_INIT%XC_F90_POINTER_T N             
                                 R                     
                                 S                                                        T                                                      1                                             U                                       	               9#         @                                  V     	              #XC_F90_LDA_VXC%XC_F90_POINTER_T W   #P Y   #NP Z   #RHO [   #V \                    @                           W     '                    #BUFFER X                D                             X                             
                                 Y                   #XC_F90_LDA_VXC%XC_F90_POINTER_T W             
                                 Z                     
                                [     
                                               \     
       #         @                                  ]     	              #XC_F90_LDA_EXC%XC_F90_POINTER_T ^   #P `   #NP a   #RHO b   #ZK c                    @                           ^     '                    #BUFFER _                D                             _                             
                                 `                   #XC_F90_LDA_EXC%XC_F90_POINTER_T ^             
                                 a                     
                                b     
                                               c     
                                                    d                                
       ) L            1275070505                                             e                                
         X            1476395011#         @                                  f     	              #XC_F90_FUNC_END%XC_F90_POINTER_T g   #P i                    @                           g     '                    #BUFFER h                D                             h                             
                                i                    #XC_F90_FUNC_END%XC_F90_POINTER_T g                                                j                                                      8                                             k                                       m               109                                             l                                                      134                                             m                                       e               101                                             n                                                      130                                             o                                       j               106                                             p                                                      131#         @                                  q     	              #XC_F90_GGA_VXC%XC_F90_POINTER_T r   #P t   #NP u   #RHO v   #SIGMA w   #VRHO x   #VSIGMA y                    @                           r     '                    #BUFFER s                D                             s                             
                                 t                   #XC_F90_GGA_VXC%XC_F90_POINTER_T r             
                                 u                     
                                v     
                
                                w     
                                               x     
                                                y     
       #         @                                  z     	              #XC_F90_GGA_EXC%XC_F90_POINTER_T {   #P }   #NP ~   #RHO    #SIGMA    #ZK                     @                           {     '                    #BUFFER |                D                             |                             
                                 }                   #XC_F90_GGA_EXC%XC_F90_POINTER_T {             
                                 ~                     
                                     
                
                                     
                                                    
                                                                                                                                                                                                                                                                                                                                                                                                                                                       #         @                                                       #NPS    #RHOS    #VXC    #EXC              
  @                                                  
  @                                                  
      p        5  p        r    p          5  p        r      5 r E       5  p        r      5 r E                              D @                                                  
       p        5  p        r    p          5  p        r      5 r E       5  p        r      5 r E                               F @                                   
       #         @                                                      #NPS    #RHOS    #VXCS    #EXC              
                                                     
                                                     
      p        5  p        r    p          5  p        r      5 r E       5  p        r      5 r E                              D                                                    
       p        5  p        r    p          5  p        r      5 r E       5  p        r      5 r E                               F @                                   
       #         @                                                     #LIBXC_GGA_SET%GLOBAL_N2    #LIBXC_GGA_SET%GLOBAL_N1    #LIBXC_GGA_SET%NG2    #LIBXC_GGA_SET%NG1    #LIBXC_GGA_SET%N    #LIBXC_GGA_SET%N3    #LIBXC_GGA_SET%N2    #LIBXC_GGA_SET%N1    #NPS    #RHOS    #VXCS    #EXC                                                                                                                                                                                                                                                                                                                                                                                                                                                       
                                                     
                                                     
      p        5  p        r    p          5  p        r      5 r E       5  p        r      5 r E                              D                                                    
       p        5  p        r    p          5  p        r      5 r E       5  p        r      5 r E                               F @                                    
       #         @                                   ¡                    #NPS ¢   #RHOS £   #EX ¤   #VXS ¥             
                                 ¢                    
                                 £                    
      p        5  p        r ¢   p          5  p        r ¢     5 r E       5  p        r ¢     5 r E                               D @                              ¤     
                F @                              ¥                    
 	      p        5  p        r ¢   p          5  p        r ¢     5 r E       5  p        r ¢     5 r E                     #         @                                   ¦                    #NPS §   #RHOS ¨   #EC ©   #VCS ª             
                                 §                    
                                 ¨                    
      p        5  p        r §   p          5  p        r §     5 r E       5  p        r §     5 r E                               D @                              ©     
                F @                              ª                    
       p        5  p        r §   p          5  p        r §     5 r E       5  p        r §     5 r E                     #         @                                   «                    #NPS ¬   #RHOS ­   #SIGMA ®   #EX ¯   #VXS °             
                                 ¬                    
                                 ­                    
      p        5  p        r ¬   p          5  p        r ¬     5 r E       5  p        r ¬     5 r E                            0  
 @                              ®                   
              &                   &                                                     D @                              ¯     
                F @                              °                    
       p        5  p        r ¬   p          5  p        r ¬     5 r E       5  p        r ¬     5 r E                     #         @                                   ±                    #NPS ²   #RHOS ³   #SIGMA ´   #EC µ   #VCS ¶             
                                 ²                    
                                 ³                    
 "     p        5  p        r ²   p          5  p        r ²     5 r E       5  p        r ²     5 r E                            0  
 @                              ´                   
 #             &                   &                                                     D @                              µ     
                F @                              ¶                    
 $      p        5  p        r ²   p          5  p        r ²     5 r E       5  p        r ²     5 r E                            $      fn#fn    Ä   @   J   CONSTANTS      @   J   XC_F90_TYPES_M    D  @   J   XC_F90_LIB_M      @   J   LIBXC_FUNCS_M !   Ä  @   J   SMPI_MATH_MODULE      `   j  PARAMETERS    d  N   J  GRID_MODULE 0   ²  \       XC_F90_POINTER_T+XC_F90_TYPES_M >     H   %   XC_F90_POINTER_T%BUFFER+XC_F90_TYPES_M=BUFFER /   V  ð      PARALLEL_TYPE+SMPI_MATH_MODULE 4   F  H   a   PARALLEL_TYPE%COMM+SMPI_MATH_MODULE 4     H   a   PARALLEL_TYPE%MYID+SMPI_MATH_MODULE 8   Ö  H   a   PARALLEL_TYPE%NUMPROCS+SMPI_MATH_MODULE 6     H   a   PARALLEL_TYPE%ROOTID+SMPI_MATH_MODULE 6   f  H   a   PARALLEL_TYPE%ISROOT+SMPI_MATH_MODULE ;   ®  H   a   PARALLEL_TYPE%NSTATE_PROC+SMPI_MATH_MODULE 7   ö  ¬   a   PARALLEL_TYPE%SUB2SUM+SMPI_MATH_MODULE <   ¢     a   PARALLEL_TYPE%MYGRID_RANGE+SMPI_MATH_MODULE :   >     a   PARALLEL_TYPE%RECVCOUNTS+SMPI_MATH_MODULE 6   Ò     a   PARALLEL_TYPE%DISPLS+SMPI_MATH_MODULE @   f	  ¬   a   PARALLEL_TYPE%GLOBAL_GRIDRANGE+SMPI_MATH_MODULE 6   
  H   a   PARALLEL_TYPE%COMM2D+SMPI_MATH_MODULE 5   Z
  H   a   PARALLEL_TYPE%COMMX+SMPI_MATH_MODULE 5   ¢
  H   a   PARALLEL_TYPE%COMMY+SMPI_MATH_MODULE 5   ê
  H   a   PARALLEL_TYPE%RANKX+SMPI_MATH_MODULE 5   2  H   a   PARALLEL_TYPE%RANKY+SMPI_MATH_MODULE 7   z     a   PARALLEL_TYPE%PERIODS+SMPI_MATH_MODULE 7     H   a   PARALLEL_TYPE%REORDER+SMPI_MATH_MODULE 7   ^     a   PARALLEL_TYPE%REMAINX+SMPI_MATH_MODULE 7   ú     a   PARALLEL_TYPE%REMAINY+SMPI_MATH_MODULE 5     H   a   PARALLEL_TYPE%NDIMS+SMPI_MATH_MODULE 4   Þ     a   PARALLEL_TYPE%DIMS+SMPI_MATH_MODULE 7   z  H   a   PARALLEL_TYPE%COMMFFT+SMPI_MATH_MODULE 7   Â  H   a   PARALLEL_TYPE%LOCAL_Z+SMPI_MATH_MODULE =   
  H   a   PARALLEL_TYPE%LOCAL_Z_START+SMPI_MATH_MODULE >   R  ¬   a   PARALLEL_TYPE%FFT_GRID_RANGE+SMPI_MATH_MODULE :   þ     a   PARALLEL_TYPE%FFT_RCOUNT+SMPI_MATH_MODULE ;        a   PARALLEL_TYPE%FFT_RDISPLS+SMPI_MATH_MODULE :   &     a   PARALLEL_TYPE%FFT_SCOUNT+SMPI_MATH_MODULE ;   º     a   PARALLEL_TYPE%FFT_SDISPLS+SMPI_MATH_MODULE 4   N         GRID_DIFF_MAP_TYPE+SMPI_MATH_MODULE ;   N  ¬   a   GRID_DIFF_MAP_TYPE%NZ_MAP+SMPI_MATH_MODULE A   ú     a   GRID_DIFF_MAP_TYPE%MYCOMM_CORES+SMPI_MATH_MODULE @     ¬   a   GRID_DIFF_MAP_TYPE%MYCOMM_SIZE+SMPI_MATH_MODULE @   B  ¬   a   GRID_DIFF_MAP_TYPE%MYSEND_SIZE+SMPI_MATH_MODULE >   î  ¬   a   GRID_DIFF_MAP_TYPE%LOCAL_MAP+SMPI_MATH_MODULE @     Ä   a   GRID_DIFF_MAP_TYPE%LOCAL_MAP1D+SMPI_MATH_MODULE =   ^  ¼   a   GRID_DIFF_MAP_TYPE%BOUNDARY+SMPI_MATH_MODULE ?        a   GRID_DIFF_MAP_TYPE%BOUNDARY1D+SMPI_MATH_MODULE ;   ¶     a   GRID_DIFF_MAP_TYPE%RCOUNT+SMPI_MATH_MODULE <   J     a   GRID_DIFF_MAP_TYPE%RDISPLS+SMPI_MATH_MODULE ;   Þ     a   GRID_DIFF_MAP_TYPE%SCOUNT+SMPI_MATH_MODULE <   r     a   GRID_DIFF_MAP_TYPE%SDISPLS+SMPI_MATH_MODULE &     Ö       GRID_TYPE+GRID_MODULE +   Ü  ¬   a   GRID_TYPE%RHOS+GRID_MODULE *        a   GRID_TYPE%RHO+GRID_MODULE +     ¬   a   GRID_TYPE%VXCS+GRID_MODULE ,   È  ¬   a   GRID_TYPE%VHXCD+GRID_MODULE +   t     a   GRID_TYPE%VLPP+GRID_MODULE )        a   GRID_TYPE%VH+GRID_MODULE +     Ä   a   GRID_TYPE%EVAL+GRID_MODULE +   `      a   GRID_TYPE%RHOC+GRID_MODULE +   ô   ¬   a   GRID_TYPE%GVEC+GRID_MODULE ,    !     a   GRID_TYPE%GMASK+GRID_MODULE +   4"  ¬   a   GRID_TYPE%RVEC+GRID_MODULE *   à"  Ä   a   GRID_TYPE%LSP+GRID_MODULE 1   ¤#  H   a   GRID_TYPE%ONEDLENGTH+GRID_MODULE *   ì#  S       PARALLEL+SMPI_MATH_MODULE !   ?$  @       NSPIN+PARAMETERS %   $  @       IXCDF+PARAMETERS=IXC %   ¿$  @       LCORE_VAL+PARAMETERS !   ÿ$  @       SNLCC+PARAMETERS %   ?%  O       SPH+GRID_MODULE=GRID !   %  @       DVOL+GRID_MODULE    Î%  p       I4B+CONSTANTS    >&  p       DP+CONSTANTS .   ®&         XC_F90_FUNC_INIT+XC_F90_LIB_M A   I'  \      XC_F90_FUNC_INIT%XC_F90_POINTER_T+XC_F90_TYPES_M O   ¥'  H   %   XC_F90_FUNC_INIT%XC_F90_POINTER_T%BUFFER+XC_F90_TYPES_M=BUFFER 0   í'  o   a   XC_F90_FUNC_INIT%P+XC_F90_LIB_M 3   \(  o   a   XC_F90_FUNC_INIT%INFO+XC_F90_LIB_M 9   Ë(  @   a   XC_F90_FUNC_INIT%FUNCTIONAL+XC_F90_LIB_M 4   )  @   a   XC_F90_FUNC_INIT%NSPIN+XC_F90_LIB_M '   K)  q       XC_LDA_X+LIBXC_FUNCS_M *   ¼)  q       XC_LDA_C_PZ+LIBXC_FUNCS_M ,   -*         XC_F90_LDA_VXC+XC_F90_LIB_M ?   ¹*  \      XC_F90_LDA_VXC%XC_F90_POINTER_T+XC_F90_TYPES_M M   +  H   %   XC_F90_LDA_VXC%XC_F90_POINTER_T%BUFFER+XC_F90_TYPES_M=BUFFER .   ]+  m   a   XC_F90_LDA_VXC%P+XC_F90_LIB_M /   Ê+  @   a   XC_F90_LDA_VXC%NP+XC_F90_LIB_M 0   
,  @   a   XC_F90_LDA_VXC%RHO+XC_F90_LIB_M .   J,  @   a   XC_F90_LDA_VXC%V+XC_F90_LIB_M ,   ,         XC_F90_LDA_EXC+XC_F90_LIB_M ?   -  \      XC_F90_LDA_EXC%XC_F90_POINTER_T+XC_F90_TYPES_M M   s-  H   %   XC_F90_LDA_EXC%XC_F90_POINTER_T%BUFFER+XC_F90_TYPES_M=BUFFER .   »-  m   a   XC_F90_LDA_EXC%P+XC_F90_LIB_M /   (.  @   a   XC_F90_LDA_EXC%NP+XC_F90_LIB_M 0   h.  @   a   XC_F90_LDA_EXC%RHO+XC_F90_LIB_M /   ¨.  @   a   XC_F90_LDA_EXC%ZK+XC_F90_LIB_M +   è.  z       MPI_REAL8+SMPI_MATH_MODULE )   b/  z       MPI_SUM+SMPI_MATH_MODULE -   Ü/  u       XC_F90_FUNC_END+XC_F90_LIB_M @   Q0  \      XC_F90_FUNC_END%XC_F90_POINTER_T+XC_F90_TYPES_M N   ­0  H   %   XC_F90_FUNC_END%XC_F90_POINTER_T%BUFFER+XC_F90_TYPES_M=BUFFER /   õ0  n   a   XC_F90_FUNC_END%P+XC_F90_LIB_M /   c1  q       XC_LDA_C_VWN_RPA+LIBXC_FUNCS_M ,   Ô1  s       XC_GGA_X_PW91+LIBXC_FUNCS_M ,   G2  s       XC_GGA_C_PW91+LIBXC_FUNCS_M +   º2  s       XC_GGA_X_PBE+LIBXC_FUNCS_M +   -3  s       XC_GGA_C_PBE+LIBXC_FUNCS_M +    3  s       XC_GGA_X_B88+LIBXC_FUNCS_M +   4  s       XC_GGA_C_LYP+LIBXC_FUNCS_M ,   4  ¦       XC_F90_GGA_VXC+XC_F90_LIB_M ?   ,5  \      XC_F90_GGA_VXC%XC_F90_POINTER_T+XC_F90_TYPES_M M   5  H   %   XC_F90_GGA_VXC%XC_F90_POINTER_T%BUFFER+XC_F90_TYPES_M=BUFFER .   Ð5  m   a   XC_F90_GGA_VXC%P+XC_F90_LIB_M /   =6  @   a   XC_F90_GGA_VXC%NP+XC_F90_LIB_M 0   }6  @   a   XC_F90_GGA_VXC%RHO+XC_F90_LIB_M 2   ½6  @   a   XC_F90_GGA_VXC%SIGMA+XC_F90_LIB_M 1   ý6  @   a   XC_F90_GGA_VXC%VRHO+XC_F90_LIB_M 3   =7  @   a   XC_F90_GGA_VXC%VSIGMA+XC_F90_LIB_M ,   }7         XC_F90_GGA_EXC+XC_F90_LIB_M ?   8  \      XC_F90_GGA_EXC%XC_F90_POINTER_T+XC_F90_TYPES_M M   q8  H   %   XC_F90_GGA_EXC%XC_F90_POINTER_T%BUFFER+XC_F90_TYPES_M=BUFFER .   ¹8  m   a   XC_F90_GGA_EXC%P+XC_F90_LIB_M /   &9  @   a   XC_F90_GGA_EXC%NP+XC_F90_LIB_M 0   f9  @   a   XC_F90_GGA_EXC%RHO+XC_F90_LIB_M 2   ¦9  @   a   XC_F90_GGA_EXC%SIGMA+XC_F90_LIB_M /   æ9  @   a   XC_F90_GGA_EXC%ZK+XC_F90_LIB_M    &:  @       N1+GRID_MODULE    f:  @       N2+GRID_MODULE    ¦:  @       N3+GRID_MODULE    æ:  @       N+GRID_MODULE     &;  @       NG1+GRID_MODULE     f;  @       NG2+GRID_MODULE &   ¦;  @       GLOBAL_N1+GRID_MODULE &   æ;  @       GLOBAL_N2+GRID_MODULE    &<  m       XC_FUNCTIONAL "   <  @   a   XC_FUNCTIONAL%NPS #   Ó<    a   XC_FUNCTIONAL%RHOS "   ×=    a   XC_FUNCTIONAL%VXC "   Û>  @   a   XC_FUNCTIONAL%EXC    ?  n       LIBXC_LDA_SET "   ?  @   a   LIBXC_LDA_SET%NPS #   É?    a   LIBXC_LDA_SET%RHOS #   Í@    a   LIBXC_LDA_SET%VXCS "   ÑA  @   a   LIBXC_LDA_SET%EXC    B  -      LIBXC_GGA_SET >   >C  @     LIBXC_GGA_SET%GLOBAL_N2+GRID_MODULE=GLOBAL_N2 >   ~C  @     LIBXC_GGA_SET%GLOBAL_N1+GRID_MODULE=GLOBAL_N1 2   ¾C  @     LIBXC_GGA_SET%NG2+GRID_MODULE=NG2 2   þC  @     LIBXC_GGA_SET%NG1+GRID_MODULE=NG1 .   >D  @     LIBXC_GGA_SET%N+GRID_MODULE=N 0   ~D  @     LIBXC_GGA_SET%N3+GRID_MODULE=N3 0   ¾D  @     LIBXC_GGA_SET%N2+GRID_MODULE=N2 0   þD  @     LIBXC_GGA_SET%N1+GRID_MODULE=N1 "   >E  @   a   LIBXC_GGA_SET%NPS #   ~E    a   LIBXC_GGA_SET%RHOS #   F    a   LIBXC_GGA_SET%VXCS "   G  @   a   LIBXC_GGA_SET%EXC    ÆG  l       LIBXC_LDA_X     2H  @   a   LIBXC_LDA_X%NPS !   rH    a   LIBXC_LDA_X%RHOS    vI  @   a   LIBXC_LDA_X%EX     ¶I    a   LIBXC_LDA_X%VXS     ºJ  l       LIBXC_VWN1RPA_C $   &K  @   a   LIBXC_VWN1RPA_C%NPS %   fK    a   LIBXC_VWN1RPA_C%RHOS #   jL  @   a   LIBXC_VWN1RPA_C%EC $   ªL    a   LIBXC_VWN1RPA_C%VCS    ®M  w       LIBXC_B88_X     %N  @   a   LIBXC_B88_X%NPS !   eN    a   LIBXC_B88_X%RHOS "   iO  ¤   a   LIBXC_B88_X%SIGMA    P  @   a   LIBXC_B88_X%EX     MP    a   LIBXC_B88_X%VXS    QQ  w       LIBXC_LYP_C     ÈQ  @   a   LIBXC_LYP_C%NPS !   R    a   LIBXC_LYP_C%RHOS "   S  ¤   a   LIBXC_LYP_C%SIGMA    °S  @   a   LIBXC_LYP_C%EC     ðS    a   LIBXC_LYP_C%VCS 
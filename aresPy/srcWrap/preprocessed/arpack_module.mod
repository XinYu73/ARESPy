  ÿ(  Y   k820309    l          18.0        ÛOb                                                                                                          
       Arpack_module.f90 ARPACK_MODULE                                                     
                                                                                                                                                                                                                                                                   `                     #LOGFIL    #NDIGIT    #MGETV0    #MSAUPD    #MSAUP2 	   #MSAITR 
   #MSEIGT    #MSAPPS    #MSGETS    #MSEUPD    #MNAUPD    #MNAUP2    #MNAITR    #MNEIGH    #MNAPPS    #MNGETS    #MNEUPD    #MCAUPD    #MCAUP2    #MCAITR    #MCEIGH    #MCAPPS    #MCGETS    #MCEUPD                                                                                                                                                                                                                                                                                              	                                                            
                                                                                                                                                                                                                                                      $                                                            (                                                            ,                                                            0                                                            4                                                            8                                                            <                                                            @                                                            D                                                            H                                                            L                                                            P                                                            T                                                            X                                                            \                                                                                                 @B             1000000                                                                                   ô              500                                                                                   '              10000#         @                                                    
   #REAL_DIAGH_ARPACK%DIFF_MAP !   #REAL_DIAGH_ARPACK%PARALLEL "   #N %   #VEFF &   #NEV '   #EVEC (   #EVAL )   #RESID_RESTART *   #NEC +   #INFO ,   #MAXMVS -   #TOL .                                             !     @      #REAL_DIAGH_ARPACK%GRID_DIFF_MAP_TYPE #                                             "     @      #REAL_DIAGH_ARPACK%PARALLEL_TYPE $             
                                 %                    
  @                              &                    
    p          5  p        r %       5  p        r %                               
@ @                              '                     D                                (                   
               &                   &                                                    D                                )                    
     p          5  p        r '       5  p        r '                               
D                                *                   
               &                                                     
D                                +                      
D @                              ,                      
D                                -                      
@ @                              .     
                       @               D           $     '@                   #COMM /   #MYID 0   #NUMPROCS 1   #ROOTID 2   #ISROOT 3   #NSTATE_PROC 4   #SUB2SUM 5   #MYGRID_RANGE 6   #RECVCOUNTS 7   #DISPLS 8   #GLOBAL_GRIDRANGE 9   #COMM2D :   #COMMX ;   #COMMY <   #RANKX =   #RANKY >   #PERIODS ?   #REORDER @   #REMAINX A   #REMAINY B   #NDIMS C   #DIMS D   #COMMFFT E   #LOCAL_Z F   #LOCAL_Z_START G   #FFT_GRID_RANGE H   #FFT_RCOUNT I   #FFT_RDISPLS J   #FFT_SCOUNT K   #FFT_SDISPLS L                                              /                                                              0                                                             1                                                             2                                                              3                                                             4                                                           5                                         &                   &                                                                                      6            x                   p          p            p                                                                  7                          	               &                                                                                   8            Ð              
               &                                                                                   9                                        &                   &                                                                                      :     x                                                       ;     |                                                       <                                                            =                                                            >                                                            ?                              p          p            p                                                                     @                                                            A                              p          p            p                                                                     B                               p          p            p                                                                     C     ¨                                                       D            ¬                  p          p            p                                                                     E     ´                                                       F     ¸                                                       G     ¼                                                     H            À                            &                   &                                                                                    I                                         &                                                                                    J            h                            &                                                                                    K            °                            &                                                                                    L            ø                            &                                                            @               D           #     '@                   #NZ_MAP M   #MYCOMM_CORES N   #MYCOMM_SIZE O   #MYSEND_SIZE P   #LOCAL_MAP Q   #LOCAL_MAP1D R   #BOUNDARY S   #BOUNDARY1D T   #RCOUNT U   #RDISPLS V   #SCOUNT W   #SDISPLS X                                           M                                          &                   &                                                                                      N            `                   p          p            p                                                                  O            h                             &                   &                                                                                   P            È                             &                   &                                                                                   Q            (                            &                   &                                                                                   R                                        &                   &                   &                                                                                      S                               p          p          p            p          p                                                                     T                              p          p            p                                                                   U                          	               &                                                                                    V            h             
               &                                                                                    W            °                            &                                                                                    X            ø                            &                                                  (      fn#fn    È   @   J   CONSTANTS      p       I4B+CONSTANTS    x  p       DP+CONSTANTS $   è  p     ARPACK_MODULE!DEBUG    X  H      LOGFIL       H      NDIGIT    è  H      MGETV0    0  H      MSAUPD    x  H      MSAUP2    À  H      MSAITR      H      MSEIGT    P  H      MSAPPS      H      MSGETS    à  H      MSEUPD    (  H      MNAUPD    p  H      MNAUP2    ¸  H      MNAITR       H      MNEIGH    H  H      MNAPPS      H      MNGETS    Ø  H      MNEUPD       H      MCAUPD    h  H      MCAUP2    °  H      MCAITR    ø  H      MCEIGH    @	  H      MCAPPS    	  H      MCGETS    Ð	  H      MCEUPD    
  w       MAXN    
  s       MAXNEV      u       MAXNCV "   w  ñ       REAL_DIAGH_ARPACK E   h  j     REAL_DIAGH_ARPACK%DIFF_MAP+SMPI_MATH_MODULE=DIFF_MAP E   Ò  e     REAL_DIAGH_ARPACK%PARALLEL+SMPI_MATH_MODULE=PARALLEL $   7  @   a   REAL_DIAGH_ARPACK%N '   w  ´   a   REAL_DIAGH_ARPACK%VEFF &   +  @   a   REAL_DIAGH_ARPACK%NEV '   k  ¤   a   REAL_DIAGH_ARPACK%EVEC '     ´   a   REAL_DIAGH_ARPACK%EVAL 0   Ã     a   REAL_DIAGH_ARPACK%RESID_RESTART &   O  @   a   REAL_DIAGH_ARPACK%NEC '     @   a   REAL_DIAGH_ARPACK%INFO )   Ï  @   a   REAL_DIAGH_ARPACK%MAXMVS &     @   a   REAL_DIAGH_ARPACK%TOL O   O  ð     REAL_DIAGH_ARPACK%PARALLEL_TYPE+SMPI_MATH_MODULE=PARALLEL_TYPE F   ?  H   a   REAL_DIAGH_ARPACK%PARALLEL_TYPE%COMM+SMPI_MATH_MODULE F     H   a   REAL_DIAGH_ARPACK%PARALLEL_TYPE%MYID+SMPI_MATH_MODULE J   Ï  H   a   REAL_DIAGH_ARPACK%PARALLEL_TYPE%NUMPROCS+SMPI_MATH_MODULE H     H   a   REAL_DIAGH_ARPACK%PARALLEL_TYPE%ROOTID+SMPI_MATH_MODULE H   _  H   a   REAL_DIAGH_ARPACK%PARALLEL_TYPE%ISROOT+SMPI_MATH_MODULE M   §  H   a   REAL_DIAGH_ARPACK%PARALLEL_TYPE%NSTATE_PROC+SMPI_MATH_MODULE I   ï  ¬   a   REAL_DIAGH_ARPACK%PARALLEL_TYPE%SUB2SUM+SMPI_MATH_MODULE N        a   REAL_DIAGH_ARPACK%PARALLEL_TYPE%MYGRID_RANGE+SMPI_MATH_MODULE L   7     a   REAL_DIAGH_ARPACK%PARALLEL_TYPE%RECVCOUNTS+SMPI_MATH_MODULE H   Ë     a   REAL_DIAGH_ARPACK%PARALLEL_TYPE%DISPLS+SMPI_MATH_MODULE R   _  ¬   a   REAL_DIAGH_ARPACK%PARALLEL_TYPE%GLOBAL_GRIDRANGE+SMPI_MATH_MODULE H     H   a   REAL_DIAGH_ARPACK%PARALLEL_TYPE%COMM2D+SMPI_MATH_MODULE G   S  H   a   REAL_DIAGH_ARPACK%PARALLEL_TYPE%COMMX+SMPI_MATH_MODULE G     H   a   REAL_DIAGH_ARPACK%PARALLEL_TYPE%COMMY+SMPI_MATH_MODULE G   ã  H   a   REAL_DIAGH_ARPACK%PARALLEL_TYPE%RANKX+SMPI_MATH_MODULE G   +  H   a   REAL_DIAGH_ARPACK%PARALLEL_TYPE%RANKY+SMPI_MATH_MODULE I   s     a   REAL_DIAGH_ARPACK%PARALLEL_TYPE%PERIODS+SMPI_MATH_MODULE I     H   a   REAL_DIAGH_ARPACK%PARALLEL_TYPE%REORDER+SMPI_MATH_MODULE I   W     a   REAL_DIAGH_ARPACK%PARALLEL_TYPE%REMAINX+SMPI_MATH_MODULE I   ó     a   REAL_DIAGH_ARPACK%PARALLEL_TYPE%REMAINY+SMPI_MATH_MODULE G     H   a   REAL_DIAGH_ARPACK%PARALLEL_TYPE%NDIMS+SMPI_MATH_MODULE F   ×     a   REAL_DIAGH_ARPACK%PARALLEL_TYPE%DIMS+SMPI_MATH_MODULE I   s  H   a   REAL_DIAGH_ARPACK%PARALLEL_TYPE%COMMFFT+SMPI_MATH_MODULE I   »  H   a   REAL_DIAGH_ARPACK%PARALLEL_TYPE%LOCAL_Z+SMPI_MATH_MODULE O     H   a   REAL_DIAGH_ARPACK%PARALLEL_TYPE%LOCAL_Z_START+SMPI_MATH_MODULE P   K  ¬   a   REAL_DIAGH_ARPACK%PARALLEL_TYPE%FFT_GRID_RANGE+SMPI_MATH_MODULE L   ÷     a   REAL_DIAGH_ARPACK%PARALLEL_TYPE%FFT_RCOUNT+SMPI_MATH_MODULE M        a   REAL_DIAGH_ARPACK%PARALLEL_TYPE%FFT_RDISPLS+SMPI_MATH_MODULE L        a   REAL_DIAGH_ARPACK%PARALLEL_TYPE%FFT_SCOUNT+SMPI_MATH_MODULE M   ³     a   REAL_DIAGH_ARPACK%PARALLEL_TYPE%FFT_SDISPLS+SMPI_MATH_MODULE Y   G         REAL_DIAGH_ARPACK%GRID_DIFF_MAP_TYPE+SMPI_MATH_MODULE=GRID_DIFF_MAP_TYPE M   G!  ¬   a   REAL_DIAGH_ARPACK%GRID_DIFF_MAP_TYPE%NZ_MAP+SMPI_MATH_MODULE S   ó!     a   REAL_DIAGH_ARPACK%GRID_DIFF_MAP_TYPE%MYCOMM_CORES+SMPI_MATH_MODULE R   "  ¬   a   REAL_DIAGH_ARPACK%GRID_DIFF_MAP_TYPE%MYCOMM_SIZE+SMPI_MATH_MODULE R   ;#  ¬   a   REAL_DIAGH_ARPACK%GRID_DIFF_MAP_TYPE%MYSEND_SIZE+SMPI_MATH_MODULE P   ç#  ¬   a   REAL_DIAGH_ARPACK%GRID_DIFF_MAP_TYPE%LOCAL_MAP+SMPI_MATH_MODULE R   $  Ä   a   REAL_DIAGH_ARPACK%GRID_DIFF_MAP_TYPE%LOCAL_MAP1D+SMPI_MATH_MODULE O   W%  ¼   a   REAL_DIAGH_ARPACK%GRID_DIFF_MAP_TYPE%BOUNDARY+SMPI_MATH_MODULE Q   &     a   REAL_DIAGH_ARPACK%GRID_DIFF_MAP_TYPE%BOUNDARY1D+SMPI_MATH_MODULE M   ¯&     a   REAL_DIAGH_ARPACK%GRID_DIFF_MAP_TYPE%RCOUNT+SMPI_MATH_MODULE N   C'     a   REAL_DIAGH_ARPACK%GRID_DIFF_MAP_TYPE%RDISPLS+SMPI_MATH_MODULE M   ×'     a   REAL_DIAGH_ARPACK%GRID_DIFF_MAP_TYPE%SCOUNT+SMPI_MATH_MODULE N   k(     a   REAL_DIAGH_ARPACK%GRID_DIFF_MAP_TYPE%SDISPLS+SMPI_MATH_MODULE 
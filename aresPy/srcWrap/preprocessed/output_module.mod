  Ō$  V   k820309    l          18.0        bQb                                                                                                          
       Energy_module.f90 OUTPUT_MODULE                                                     
                                                           
                @          @                              
                                                           
       ESHIFT_PS ESHIFT_TOT                   @               @                '@                   #COMM    #MYID    #NUMPROCS    #ROOTID 	   #ISROOT 
   #NSTATE_PROC    #SUB2SUM    #MYGRID_RANGE    #RECVCOUNTS    #DISPLS    #GLOBAL_GRIDRANGE    #COMM2D    #COMMX    #COMMY    #RANKX    #RANKY    #PERIODS    #REORDER    #REMAINX    #REMAINY    #NDIMS    #DIMS    #COMMFFT    #LOCAL_Z    #LOCAL_Z_START    #FFT_GRID_RANGE    #FFT_RCOUNT     #FFT_RDISPLS !   #FFT_SCOUNT "   #FFT_SDISPLS #                                                                                                                                                                                                                                      	                                                              
                                                                                                                                                                 &                   &                                                                                                  x                   p          p            p                                                                                            	               &                                                                                               Š              
               &                                                                                                                           &                   &                                                                                           x                                                            |                                                                                                                                                                                                                                                                         p          p            p                                                                                                                                                               p          p            p                                                                                                    p          p            p                                                                          Ø                                                                   ¬                  p          p            p                                                                          “                                                            ø                                                            ¼                                                                 Ą                            &                   &                                                                                                                              &                                                                                    !            h                            &                                                                                    "            °                            &                                                                                    #            ų                            &                                                             @               @           $     '@                   #NZ_MAP %   #MYCOMM_CORES &   #MYCOMM_SIZE '   #MYSEND_SIZE (   #LOCAL_MAP )   #LOCAL_MAP1D *   #BOUNDARY +   #BOUNDARY1D ,   #RCOUNT -   #RDISPLS .   #SCOUNT /   #SDISPLS 0                                           %                                          &                   &                                                                                      &            `                   p          p            p                                                                  '            h                             &                   &                                                                                   (            Č                             &                   &                                                                                   )            (                            &                   &                                                                                   *                                        &                   &                   &                                                                                      +                               p          p          p            p          p                                                                     ,                              p          p            p                                                                   -                          	               &                                                                                    .            h             
               &                                                                                    /            °                            &                                                                                    0            ų                            &                                                                                        1     @      #PARALLEL_TYPE                 @                              2     
                                                   3     
                                                    4                                                         #         @                                  5                    #INLABEL 6   #FLAG 7                                            6                     1                                            7                                                        8     
                                                   9     
                 
                 hĖm96;@        27.2113834279111D0                                            :     
                                                   ;     
                                                   <     
                                                   =     
                                                   >     
                                                   ?     
                                                   @     
                                                    A                                                                                                    B                                                       C                                                       D                                                       E            #         @                                   F                    #OUTPUT%GLOBAL_N2 G   #OUTPUT%GLOBAL_N1 H   #OUTPUT%NG2 I   #OUTPUT%NG1 J   #OUTPUT%N K                                            G                                                     H                                                     I                                                     J                                                     K            #         @                                  L                    #WRITE_DENSITY%N3 M   #WRITE_DENSITY%N2 N   #WRITE_DENSITY%N1 O   #WRITE_DENSITY%NG2 P   #WRITE_DENSITY%NG1 Q   #WRITE_DENSITY%NI3 R   #WRITE_DENSITY%NI2 S   #WRITE_DENSITY%NI1 T                                            M                                                     N                                                     O                                                     P                                                     Q                                                       R                                                       S                                                       T            #         @                                  U                            (      fn#fn    Č   @   J   CONSTANTS      @   J   ENERGY_MODULE !   H  @   J   SMPI_MATH_MODULE      U   J  STRUCT_MODULE /   Ż  š      PARALLEL_TYPE+SMPI_MATH_MODULE 4   Ķ  H   a   PARALLEL_TYPE%COMM+SMPI_MATH_MODULE 4     H   a   PARALLEL_TYPE%MYID+SMPI_MATH_MODULE 8   ]  H   a   PARALLEL_TYPE%NUMPROCS+SMPI_MATH_MODULE 6   „  H   a   PARALLEL_TYPE%ROOTID+SMPI_MATH_MODULE 6   ķ  H   a   PARALLEL_TYPE%ISROOT+SMPI_MATH_MODULE ;   5  H   a   PARALLEL_TYPE%NSTATE_PROC+SMPI_MATH_MODULE 7   }  ¬   a   PARALLEL_TYPE%SUB2SUM+SMPI_MATH_MODULE <   )     a   PARALLEL_TYPE%MYGRID_RANGE+SMPI_MATH_MODULE :   Å     a   PARALLEL_TYPE%RECVCOUNTS+SMPI_MATH_MODULE 6   Y     a   PARALLEL_TYPE%DISPLS+SMPI_MATH_MODULE @   ķ  ¬   a   PARALLEL_TYPE%GLOBAL_GRIDRANGE+SMPI_MATH_MODULE 6     H   a   PARALLEL_TYPE%COMM2D+SMPI_MATH_MODULE 5   į  H   a   PARALLEL_TYPE%COMMX+SMPI_MATH_MODULE 5   )	  H   a   PARALLEL_TYPE%COMMY+SMPI_MATH_MODULE 5   q	  H   a   PARALLEL_TYPE%RANKX+SMPI_MATH_MODULE 5   ¹	  H   a   PARALLEL_TYPE%RANKY+SMPI_MATH_MODULE 7   
     a   PARALLEL_TYPE%PERIODS+SMPI_MATH_MODULE 7   
  H   a   PARALLEL_TYPE%REORDER+SMPI_MATH_MODULE 7   å
     a   PARALLEL_TYPE%REMAINX+SMPI_MATH_MODULE 7        a   PARALLEL_TYPE%REMAINY+SMPI_MATH_MODULE 5     H   a   PARALLEL_TYPE%NDIMS+SMPI_MATH_MODULE 4   e     a   PARALLEL_TYPE%DIMS+SMPI_MATH_MODULE 7     H   a   PARALLEL_TYPE%COMMFFT+SMPI_MATH_MODULE 7   I  H   a   PARALLEL_TYPE%LOCAL_Z+SMPI_MATH_MODULE =     H   a   PARALLEL_TYPE%LOCAL_Z_START+SMPI_MATH_MODULE >   Ł  ¬   a   PARALLEL_TYPE%FFT_GRID_RANGE+SMPI_MATH_MODULE :        a   PARALLEL_TYPE%FFT_RCOUNT+SMPI_MATH_MODULE ;        a   PARALLEL_TYPE%FFT_RDISPLS+SMPI_MATH_MODULE :   ­     a   PARALLEL_TYPE%FFT_SCOUNT+SMPI_MATH_MODULE ;   A     a   PARALLEL_TYPE%FFT_SDISPLS+SMPI_MATH_MODULE 4   Õ         GRID_DIFF_MAP_TYPE+SMPI_MATH_MODULE ;   Õ  ¬   a   GRID_DIFF_MAP_TYPE%NZ_MAP+SMPI_MATH_MODULE A        a   GRID_DIFF_MAP_TYPE%MYCOMM_CORES+SMPI_MATH_MODULE @     ¬   a   GRID_DIFF_MAP_TYPE%MYCOMM_SIZE+SMPI_MATH_MODULE @   É  ¬   a   GRID_DIFF_MAP_TYPE%MYSEND_SIZE+SMPI_MATH_MODULE >   u  ¬   a   GRID_DIFF_MAP_TYPE%LOCAL_MAP+SMPI_MATH_MODULE @   !  Ä   a   GRID_DIFF_MAP_TYPE%LOCAL_MAP1D+SMPI_MATH_MODULE =   å  ¼   a   GRID_DIFF_MAP_TYPE%BOUNDARY+SMPI_MATH_MODULE ?   ”     a   GRID_DIFF_MAP_TYPE%BOUNDARY1D+SMPI_MATH_MODULE ;   =     a   GRID_DIFF_MAP_TYPE%RCOUNT+SMPI_MATH_MODULE <   Ń     a   GRID_DIFF_MAP_TYPE%RDISPLS+SMPI_MATH_MODULE ;   e     a   GRID_DIFF_MAP_TYPE%SCOUNT+SMPI_MATH_MODULE <   ł     a   GRID_DIFF_MAP_TYPE%SDISPLS+SMPI_MATH_MODULE *     S       PARALLEL+SMPI_MATH_MODULE (   ą  @       ESHIFT_PS+STRUCT_MODULE )      @       ESHIFT_TOT+STRUCT_MODULE    `  p       I4B+CONSTANTS ,   Š  _       WRITE_TIME+SMPI_MATH_MODULE 4   /  L   a   WRITE_TIME%INLABEL+SMPI_MATH_MODULE 1   {  @   a   WRITE_TIME%FLAG+SMPI_MATH_MODULE $   »  @       EBAND+ENERGY_MODULE "   ū         HART2EV+CONSTANTS "   }  @       EXC+ENERGY_MODULE !   ½  @       EH+ENERGY_MODULE #   ż  @       EEXT+ENERGY_MODULE *   =  @       EII+STRUCT_MODULE=EIONION #   }  @       ETOT+ENERGY_MODULE !   ½  @       FE+ENERGY_MODULE "   ż  @       FE0+ENERGY_MODULE    =  p       DP+CONSTANTS    ­  @       TIME_TOTAL0    ķ  @       TIME_TOTAL1    -  @       TIME_SCF0    m  @       TIME_SCF1    ­  ¢       OUTPUT 7   O   @     OUTPUT%GLOBAL_N2+GRID_MODULE=GLOBAL_N2 7      @     OUTPUT%GLOBAL_N1+GRID_MODULE=GLOBAL_N1 +   Ļ   @     OUTPUT%NG2+GRID_MODULE=NG2 +   !  @     OUTPUT%NG1+GRID_MODULE=NG1 '   O!  @     OUTPUT%N+GRID_MODULE=N    !  ż       WRITE_DENSITY 0   "  @     WRITE_DENSITY%N3+GRID_MODULE=N3 0   Ģ"  @     WRITE_DENSITY%N2+GRID_MODULE=N2 0   #  @     WRITE_DENSITY%N1+GRID_MODULE=N1 2   L#  @     WRITE_DENSITY%NG2+GRID_MODULE=NG2 2   #  @     WRITE_DENSITY%NG1+GRID_MODULE=NG1 8   Ģ#  @     WRITE_DENSITY%NI3+GRID_MODULE=GLOBAL_N3 8   $  @     WRITE_DENSITY%NI2+GRID_MODULE=GLOBAL_N2 8   L$  @     WRITE_DENSITY%NI1+GRID_MODULE=GLOBAL_N1    $  H       WRITE_BAND 
# 1 "Read_module.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "Read_module.f90"
MODULE read_module
   !#########################READ MODULE#########################!
   !* FOR    : read the input data                               !
   !* Author : Qiang Xu                                          !
   !* Date   : 2017/07/04                                        !
   !#############################################################!
   USE constants

   USE smpi_math_module

   IMPLICIT NONE
   !> character array
   !> reserve the attribute for UPF
   INTERFACE get_value
     MODULE PROCEDURE get_value_int,get_value_real
   END INTERFACE
   TYPE attribute
     CHARACTER(len=120)  :: value
   ENDTYPE attribute
   TYPE(attribute),ALLOCATABLE :: attribute_data(:)
CONTAINS
   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !------------------------read input.dat------------------------
   SUBROUTINE read_file(infile)
      USE parameters
      USE math, ONLY:change_case,find_keywords
      USE m_time_evaluate, only: memory_sum
      IMPLICIT NONE
      !-----------------------------------------------------------
      CHARACTER(*),INTENT(IN) :: infile
      !-----------------------------------------------------------
      INTEGER(I4B)    :: l_str
      INTEGER(I4B)    :: ios,id_ex,id_pound,id_key,id_value !id of key words
      INTEGER(I4B)    :: i,j,k
      INTEGER(I4B)    :: nele !the number of elements
      INTEGER(I4B)    :: vmajor,vminor,vmicro ! XC version
      !
      CHARACTER       :: ch_mark
      CHARACTER(100)  :: instr
      CHARACTER(100)  :: str
      !
      LOGICAL         :: lexist
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      IF(parallel%isroot)THEN

      !-------start input.dat-----------
      INQUIRE(FILE=infile,EXIST=lexist)
      !test the input.dat
      IF(.NOT.lexist)THEN
         WRITE(6,*) '>>>WARNING<<<:input file is not exist.stop!!!'
         CALL generate_infile(infile)
         STOP
      ENDIF
      !
      OPEN(100,FILE=infile,STATUS='old')
      ch_mark='='
      DO WHILE(.TRUE.)
         READ(100,'(A100)',IOSTAT=ios) instr
         IF( ios /= 0 ) EXIT
         instr=ADJUSTL(instr)
         !change the upper case to lower case
         CALL change_case(instr,str,2)
         !the length of strings
         l_str=LEN_TRIM(str)
         !mark '!' '#'
         id_ex=INDEX(str(1:l_str),'!')
         id_pound=INDEX(str(1:l_str),'#')
         !cancel the strings after '!' or '#'
         IF (id_ex > 0 .AND. id_pound > 0 ) THEN
            l_str=MIN(l_str,id_ex-1,id_pound-1)
         ELSEIF(id_ex > 0 .AND. id_pound == 0 ) THEN
            l_str=MIN(l_str,id_ex-1)
         ELSEIF(id_ex == 0 .AND. id_pound > 0 ) THEN
            l_str=MIN(l_str,id_pound-1)
         ENDIF
         l_str=LEN_TRIM(str(1:l_str))
         IF(l_str<2) CYCLE
         !
         CALL find_keywords(str,ch_mark,id_key,id_value)
         !System name,no use
         IF(str(1:id_key) =='system')THEN
            str=instr
            READ(str(id_value:l_str),*) system_name
            WRITE(6,*) 'Task name >>> ',system_name
         ENDIF
         !cell name
         IF(str(1:id_key)=='cellfile')THEN
            str=instr
            cellfile_name=TRIM(ADJUSTL(str(id_value:l_str)))
         ENDIF
         !output file name !!!> this does not mpi_bcast
         IF(str(1:id_key)=='outfile')THEN
            str=instr
            outfile=TRIM(ADJUSTL(str(id_value:l_str)))
         ENDIF
         !input elements
         IF( str(1:id_key) == 'elements')THEN
            nele=0
            k=id_value
            DO i=id_value,l_str
               IF(str(i:i) == " ")THEN
                  nele=nele+1
                  elements(nele)=str(k:i-1)
               ENDIF
               IF (str(i-1:i-1) == " " .AND. str(i:i) /= " ")THEN
                  k=i
               ENDIF
            ENDDO
            nele=nele+1
            elements(nele)=str(k:l_str)
         ENDIF
         !pseudopotential file name
         IF(str(1:id_key)=='ppfile')THEN
            str=instr
            ntype=0
            DO i=id_value,l_str
               IF(str(i:i) == " ")THEN
                  ntype=ntype+1
                  ppfile_name(ntype)=str(k:i-1)
               ENDIF
               IF (str(i-1:i-1) == " " .and. str(i:i) /= " ") then
                  k=i
               ENDIF
            ENDDO
            ntype=ntype+1
            ppfile_name(ntype)=str(k:l_str)
         ENDIF
         !Ecut
         IF(str(1:id_key) == 'ecut')THEN
            READ(str(id_value:l_str),*) Ecut
            IF(Ecut>0.d0)THEN
               Ecut=Ecut*ev2hart
               init_gap=SQRT((pi**2/Ecut)/2)
            ENDIF
         ENDIF
         !gridsize
         IF(str(1:id_key) == 'gridsize')THEN
            IF(Ecut<=0.d0)THEN
               READ(str(id_value:l_str),*) init_gap
               init_gap=init_gap*ang2bohr
            ENDIF
         ENDIF
         IF(str(1:id_key) == 'gridsizebohr')THEN
            IF(Ecut<=0.d0)THEN
               READ(str(id_value:l_str),*) init_gap
               !init_gap=init_gap*ang2bohr
            ENDIF
         ENDIF
         !For isolate system
         IF(str(1:id_key)=='lbvk')THEN
            IF(str(id_value:id_value)=='t')THEN
                WRITE(6,*) '[KS Task] BvK Cell'
                LBvK=.TRUE.
            ELSE
                WRITE(6,*) 'KS Task [k-Space representation Cell]'
                LBvK=.FALSE.
            ENDIF
         ENDIF
         !For subsystem DFT
         IF(str(1:id_key)=='lsub')THEN
            LSUB=.TRUE.
         ELSE
            LSUB=.FALSE.
         ENDIF
         !For FAT cycle
         IF(str(1:id_key)=='lfat')THEN
            IF(str(id_value:id_value)=='t')THEN
               LFAT=.TRUE.
            ELSE
               LFAT=.FALSE.
            ENDIF
         ENDIF
         !For only one subsystem
         IF(str(1:id_key)=='lone')THEN
            IF(str(id_value:id_value)=='t')THEN
               LONE=.TRUE.
            ELSE
               LONE=.FALSE.
            ENDIF
         ENDIF
         !The number of subsystem
         IF(str(1:id_key)=='nsub')THEN
            READ(str(id_value:l_str),*) Nsub
         ENDIF
         !The number of fat time
         IF(str(1:id_key)=='nfat')THEN
            READ(str(id_value:l_str),*) NFAT
            IF(NFAT<0) NFAT=0
         ENDIF
         !For OFDFT
         !For isolate system
         IF(str(1:id_key)=='lofdft')THEN
            IF(str(id_value:id_value)=='t')THEN
                LOFDFT=.TRUE.
            ELSE
                LOFDFT=.FALSE.
            ENDIF
         ENDIF
         !coe of aTFbvW
         IF(str(1:id_key)=='tfvw')THEN
            READ(str(id_value:l_str),*) TFVW
         ENDIF
         !k-grids
         IF(str(1:id_key)=='kspacing')THEN
            IF(LBvK)THEN
               !For BvK cell, Just need gamma-point
               kspacing=-1.d0
            ELSE
               READ(str(id_value:l_str),*) kspacing
            ENDIF
            if(kspacing.le.0)then
               Lgamma=.true.
            endif
         ENDIF
         !for k-point file
         IF(str(1:id_key)=='lkpt')THEN
            IF(str(id_value:id_value)=='t')THEN
                LKP=.TRUE.
                WRITE(6,*) "K-points by [KPOINTS]"
            ELSE
                LKP=.FALSE.
                WRITE(6,*) "K-points by [KSPACING]"
            ENDIF
         ENDIF
         !KEDF type
         IF( str(1:id_key) == 'kedf' )THEN
            READ(str(id_value:l_str),*) iKEDF
         ENDIF
         !exchange-correlation type
         IF( str(1:id_key) == 'xcdf' )THEN
            READ(str(id_value:l_str),*) iXCDF
            CALL xc_f90_version(vmajor,vminor,vmicro)
            WRITE(6,'(1X,"Libxc version:",I1,".",I1,".",I1)') &
             & vmajor,vminor,vmicro
         ENDIF
         !finite difference order
         IF( str(1:id_key) == 'norder' )THEN
            READ(str(id_value:l_str),"(i2)") finite_order
            finite_order=finite_order/2
            IF (finite_order>20 .or. finite_order<0) THEN
               WRITE(6,*) "STOP!!:the order of finite difference should in [1-10]"
               STOP
            ENDIF
            WRITE(6,*) "The order of finite difference is:",finite_order*2
         ENDIF
         !diag(H)
         IF(str(1:id_key)=='nadst')THEN
            READ(str(id_value:l_str),*) NADDSTATES
         ENDIF
         !>>>spin
         IF(str(1:id_key)=='lspin')THEN
            IF(str(id_value:id_value)=='t')THEN
                NSPIN=2
                WRITE(6,*) "Spin polarized"
            ELSE
                NSPIN=1
                WRITE(6,*) "Spin unpolarized"
            ENDIF
         ENDIF
         !---
         !>>>smearing for metal
         IF( str(1:id_key) == 'nsmear' )THEN
            READ(str(id_value:l_str),*) Nsmear
         ENDIF
         !---
         IF( str(1:id_key) == 'wsmear' )THEN
            READ(str(id_value:l_str),*) Wsmear
            !Wsmear to hartree
            Wsmear=Wsmear/hart2ev
         ENDIF
         !<<<smearing for metal
         IF( str(1:id_key) == 'istart' )THEN
            READ(str(id_value:l_str),*) ISTART
         ENDIF
         !Idiag
         IF(str(1:id_key) == 'idiag')THEN
            READ(str(id_value:l_str),*) Idiag
         ENDIF
         !CheM
         IF(str(1:id_key) == 'chem')THEN
            READ(str(id_value:l_str),*) CheM
            IF(Idiag/=0)THEN
                WRITE(6,*) "Chebyshev filter is used,order:",CheM
            ENDIF
         ENDIF
         !CheM0
         IF(str(1:id_key) == 'chem0')THEN
            READ(str(id_value:l_str),*) CheM0
            IF(Idiag/=0)THEN
                ! WRITE(6,*) "first Chebyshev filter order:",CheM0
            ENDIF
         ENDIF
         !first by arpack?
         IF(str(1:id_key)=='lfirst')THEN
            IF(str(id_value:id_value)=='t')THEN
                LFIRST=.TRUE.
                IF(Idiag==1)THEN
                   WRITE(6,*) 'Chebyshev:first step by ARPACK'
                ENDIF
            ELSE
                LFIRST=.FALSE.
                IF(Idiag==1)THEN
                   WRITE(6,*) 'Chebyshev:first step diag free'
                ENDIF
            ENDIF
         ENDIF
         !For partial RR
         IF(str(1:id_key) == 'nprr')THEN
            READ(str(id_value:l_str),*) NPRR
            IF(NPRR>=0.AND.LBvK)THEN
                WRITE(6,*) "[Partial Rayleigh-Ritz] Fraction Occupation:",NPRR
            ENDIF
         ENDIF
         !Lorthnorm
         IF(str(1:id_key)=='lorthnorm')THEN
            IF(str(id_value:id_value)=='t')THEN
                LRROrthNorm=.TRUE.
                WRITE(6,*) 'RayleighRitz need OrthNorm step'
            ELSE
                LRROrthNorm=.FALSE.
                WRITE(6,*) "RayleighRitz needn't OrthNorm step"
            ENDIF
         ENDIF
         !Lrandom
         IF(str(1:id_key)=='lrandom')THEN
            IF(str(id_value:id_value)=='t')THEN
                Lrandom=.TRUE.
                WRITE(6,*) 'Initialize the subspace by full random'
            ELSE
                Lrandom=.FALSE.
                WRITE(6,*) "Initialize the subspace by Slater's orbital"
            ENDIF
         ENDIF
         !linrho
         IF(str(1:id_key)=='linrho')THEN
            IF(str(id_value:id_value)=='t')THEN
                LINRHO=.TRUE.
                WRITE(6,*) 'initial density by input'
            ELSE
                LINRHO=.FALSE.
                WRITE(6,*) 'initial density by program'
            ENDIF
         ENDIF
         !>>>mixer data
         ! id of mix method select
         IF( str(1:id_key) == 'imixer' )THEN
            READ(str(id_value:l_str),*) IMIXER
         ENDIF
         !max step
         IF( str(1:id_key) == 'nmiter' )THEN
            READ(str(id_value:l_str),*) NMITER
         ENDIF
         !simple mixing step
         IF(str(1:id_key)=='nsmix')THEN
            READ(str(id_value:l_str),*) NSMIX
            IF(NSMIX>NMITER)THEN
               WRITE(6,*) 'STOP!!:the simple step shoule less than the max step'
            ENDIF
         ENDIF
         !history step number
         IF(str(1:id_key)=='nhmix')THEN
            READ(str(id_value:l_str),*) NHMIX
            IF(NHMIX>NMITER)THEN
               WRITE(6,*) 'STOP!!:the history step shoule less than the max step'
            ENDIF
         ENDIF
         !simple mixing factor
         IF(str(1:id_key)=='malpha')THEN
            READ(str(id_value:l_str),*) MALPHA
            IF(MALPHA>1.d0.OR.MALPHA<0.d0)THEN
               WRITE(6,*) 'STOP!!:the MALPHA should in [0,1] '
            ENDIF
         ENDIF
         !Anderson mixing factor
         IF(str(1:id_key)=='mbeta')THEN
            READ(str(id_value:l_str),*) mbeta
            IF(MALPHA>1.d0.OR.MALPHA<0.d0)THEN
               WRITE(6,*) 'STOP!!:the MBETA should in [0,1] '
            ENDIF
         ENDIF
         !kerker
         IF(str(1:id_key)=='amix')THEN
            READ(str(id_value:l_str),*) AMIX
            IF(AMIX<0.d0)THEN
              AMIX=0.4_DP
            ENDIF
         ENDIF
         !resta
         IF(str(1:id_key)=='resta')THEN
            !> resta(3)-> q0,epsilon0,Rs
            READ(str(id_value:l_str),*) resta(3)
         ENDIF

         IF(str(1:id_key)=='bmix')THEN
            READ(str(id_value:l_str),*) BMIX
            BMIX=BMIX*bohr2ang
         ENDIF
         !for ill matrix
         IF(str(1:id_key)=='w0am')THEN
            READ(str(id_value:l_str),*) W0AM
         ENDIF
         !tolerance of density convergence
         IF(str(1:id_key)=='rtol')THEN
            READ(str(id_value:l_str),*) RTOL
         ENDIF
         !tolerance of total energy
         IF(str(1:id_key)=='etol')THEN
            READ(str(id_value:l_str),*) ETOL
         ENDIF
         !<<<mixer data
         !For isolate system
         IF(str(1:id_key)=='lband')THEN
            IF(str(id_value:id_value)=='t')THEN
                Lband=.TRUE.
            ELSE
                Lband=.FALSE.
            ENDIF
         ENDIF
         !===================OPTIMIZE====================
         !relaxion step
         !max step
         IF( str(1:id_key) == 'nssp' )THEN
            READ(str(id_value:l_str),*) Nssp
            Nssp=MAX(Nssp,0)
            WRITE(6,*) 'Max simulate steps is:',Nssp
         ENDIF
         !> optimize method
         IF( str(1:id_key) == 'iopm' )THEN
            READ(str(id_value:l_str),*) iopm
            if(iopm==3)WRITE(6,*) 'Optimize Method is: ',"fire"
            if(iopm==4)WRITE(6,*) 'Max simulate steps is: ',"CG"
         ENDIF
         !> steps of SCF with fixed rho
         IF( str(1:id_key) == 'step_fixrho' )THEN
            READ(str(id_value:l_str),*) step_fixrho
            WRITE(6,*) 'Steps with fixed charge rho : ',step_fixrho
         ENDIF
         !tolerance of optimize force (ipom=4)
         IF(str(1:id_key)=='ediffg')THEN
            READ(str(id_value:l_str),*) ediffg
         ENDIF
         !> Pesudopotential formation option
         IF( str(1:id_key) == 'pp_opt' )THEN
           READ(str(id_value:l_str),*) PP_identifer
           IF( PP_identifer .gt. 1 )stop "please check input.dat: PP_OPT error"
         ENDIF
         !===========================================================
         !##isolate boundary condition
         !periodic? and lattice
         IF( str(1:id_key) == 'lradius_auto')THEN
           IF(str(id_value:id_value)=='t')THEN
               Lradius_auto=.TRUE.
           ELSE
               Lradius_auto=.FALSE.
           ENDIF
         ENDIF
         IF( str(1:id_key) == 'lpbc')THEN
           IF(str(id_value:id_value)=='t')THEN
               Lpbc=.TRUE.
           ELSE
               Lpbc=.FALSE.
           ENDIF
         ENDIF
         IF( str(1:id_key) == 'calforce')THEN
           IF(str(id_value:id_value)=='t')THEN
               CalForce=.TRUE.
           ELSE
               CalForce=.FALSE.
           ENDIF
         ENDIF
         IF( str(1:id_key) == 'cell_shape')THEN
           READ(str(id_value:l_str),*) cell_shape
         ENDIF
         IF( str(1:id_key) == 'cell_thick')THEN
           READ(str(id_value:l_str),*) cell_thick
         ENDIF
         ! IF( str(1:id_key) == 'hartree_direct')THEN
         !   IF(str(id_value:id_value)=='t')THEN
         !       Hartree_direct=.TRUE.
         !   ELSE
         !       Hartree_direct=.FALSE.
         !   ENDIF
         ! ENDIF
         ! IF( str(1:id_key) == 'mcm_flag')THEN
         !   IF(str(id_value:id_value)=='t')THEN
         !       MCM_flag=.TRUE.
         !   ELSE
         !       MCM_flag=.FALSE.
         !   ENDIF
         ! ENDIF
         IF(str(1:id_key)=='hartree_method')THEN
            READ(str(id_value:l_str),*) hartree_method
         ENDIF
         IF( str(1:id_key) == 'addcharge')THEN
           READ(str(id_value:l_str),*) ADDcharge
         ENDIF
         IF( str(1:id_key) == 'isormax')THEN
           READ(str(id_value:l_str),*) IsoRmax
           RadiusMax=IsoRmax*ang2bohr
           Lcellbohr=RadiusMax*2.d0
           Lcell=IsoRmax*2.d0
           IF(2.d0*IsoRmax>Lcell)THEN
             WRITE(6,*)"Isolation: rmax exceed Lcell"
             stop
           ENDIF
         ENDIF
         IF( str(1:id_key) == 'isormaxbohr')THEN
           READ(str(id_value:l_str),*) IsoRmaxbohr
           IsoRmax=Isormaxbohr*bohr2ang
           RadiusMax=Isormaxbohr
           Lcellbohr=IsoRmaxbohr*2.d0
           Lcell=IsoRmax*2.d0
           IF(2.d0*IsoRmax>Lcell)THEN
             WRITE(6,*)"Isolation: rmax exceed Lcell"
             stop
           ENDIF
         ENDIF
         !> Lcell : not necessary
         IF( str(1:id_key) == 'lcell')THEN
            READ(str(id_value:l_str),*) Lcell
         ENDIF
         IF( str(1:id_key) == 'lcellbohr')THEN
            READ(str(id_value:l_str),*) Lcellbohr
            Lcell=Lcellbohr*bohr2ang
         ENDIF
         !> Lcell END
         IF( str(1:id_key) == 'nvc')THEN
           READ(str(id_value:l_str),*) NVC
         ENDIF
         IF( str(1:id_key) == 'isonorder')THEN
           READ(str(id_value:l_str),*) ISOnorder
         ENDIF
         IF( str(1:id_key) == 'tolcg')THEN
           READ(str(id_value:l_str),*) TOLCG
         ENDIF
         IF( str(1:id_key) == 'nfcd')THEN
           READ(str(id_value:l_str),*) NFCD
         ENDIF
         IF( str(1:id_key) == 'isolmax')THEN
           READ(str(id_value:l_str),*) ISOLmax
         ENDIF
         IF( str(1:id_key) == 'fmmprec')THEN
           READ(str(id_value:l_str),*) iprec_fmm
         ENDIF
         ! IF( str(1:id_key) == 'fmm')THEN
         !   IF(str(id_value:id_value)=='t')THEN
         !       IFMM=.TRUE.
         !   ELSE
         !       IFMM=.FALSE.
         !   ENDIF
         ! ENDIF
         !======================Double-grid=====================
         !LDG
         IF(str(1:id_key)=='ldg')THEN
            IF(str(id_value:id_value)=='t')THEN
                LDG=.TRUE.
                WRITE(6,*) 'Use time-saving double-grid Method'
            ELSE
                LDG=.FALSE.
            ENDIF
         ENDIF
         !NDG
         IF( str(1:id_key) == 'ndg' )THEN
            READ(str(id_value:l_str),*) NDG
            IF(LDG)THEN
               WRITE(6,*) '[Order of double-grid]',NDG
            ENDIF
         ENDIF
         !N_near
         IF( str(1:id_key) == 'n_near' )THEN
            READ(str(id_value:l_str),*) n_near
            IF(LDG)THEN
               WRITE(6,*) '[Num near corse-grid(each side)]',n_near
            ENDIF
         ENDIF
         IF( str(1:id_key) == 'inpol' )THEN
            READ(str(id_value:l_str),*) inpol
            IF(LDG)THEN
               WRITE(6,*) '[interpolation order for dense-grid]',inpol
            ENDIF
         ENDIF
         !=======================INIT_MO==================
         !parameter=0,sto;parameter=1,random;2,MO
         IF( str(1:id_key) == 'idinit' )THEN
            READ(str(id_value:l_str),*) IDinit
         ENDIF
         IF( str(1:id_key) == 'mo_file' )THEN
            READ(str(id_value:l_str),*) MO_file
            MO_file=trim(adjustl(MO_file))
            !element first character is capitalization
            CALL change_case(MO_file(1:1),MO_file(1:1),1)
         ENDIF
         !> fix atom in x,y,z,xy,xz,yz
         IF(str(1:id_key) == 'fix_xyz')THEN
            READ(str(id_value:l_str),*)fix_xyz
            write(6,*)'## optimize with fixed direction ##'
            write(6,*)'## x    y    z    xy    xz    yz ##'
            write(6,'(6I5)')fix_xyz
         ENDIF

         !> LGamma for gamma point calculation
         IF(str(1:id_key)=='lgamma')THEN
            IF(str(id_value:id_value)=='t')THEN
                LGamma=.TRUE.
                KSPACING=-1
            ELSE
                LGamma=.FALSE.
            ENDIF
         ENDIF

      ENDDO

      CLOSE(100)
      !--------end input.dat------------
      IF(nele/=ntype)THEN
        WRITE(6,*) 'STOP!!!!!!'
        WRITE(6,*) 'Nelements',nele
        WRITE(6,*) 'Ntype',ntype
        WRITE(6,*) 'the pseudopotential files number or elements is wrong'
        STOP
      ENDIF

      ENDIF !> read on process 0 finished
      !> broadcast varables necessary
      !>> 
      CALL MPI_BCAST(finite_order ,1,MPI_INTEGER4 ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(KSPACING     ,1,MPI_REAL8    ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(LKP          ,1,MPI_LOGICAL  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(Lpbc         ,1,MPI_LOGICAL  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(Nspin        ,1,MPI_INTEGER4 ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(init_gap     ,1,MPI_REAL8    ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(LBvK         ,1,MPI_LOGICAL  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(Idiag        ,1,MPI_INTEGER4 ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(NPRR         ,1,MPI_INTEGER4 ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(LRadRho      ,1,MPI_LOGICAL  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(NMITER       ,1,MPI_INTEGER4 ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(IMIXER       ,1,MPI_INTEGER4 ,parallel%rootid,parallel%comm,mpinfo)
      !>> converge precision
      CALL MPI_BCAST(RTOL         ,1,MPI_REAL8    ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(ETOL         ,1,MPI_REAL8    ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(Nstates      ,1,MPI_INTEGER4 ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(LFIRST       ,1,MPI_LOGICAL  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(LINRHO       ,1,MPI_LOGICAL  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(Lrandom      ,1,MPI_LOGICAL  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(Nsmear       ,1,MPI_INTEGER4 ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(Wsmear       ,1,MPI_REAL8    ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(Lsub         ,1,MPI_LOGICAL  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(Nsub         ,1,MPI_INTEGER4 ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(LFAT         ,1,MPI_LOGICAL  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(Lone         ,1,MPI_LOGICAL  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(NFAT         ,1,MPI_INTEGER4 ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(CheM         ,1,MPI_INTEGER4 ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(LRROrthNorm  ,1,MPI_LOGICAL  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(RadiusMax    ,1,MPI_REAL8    ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(NVC          ,1,MPI_INTEGER4  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(ISOnorder    ,1,MPI_INTEGER4  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(NFCD         ,1,MPI_INTEGER4  ,parallel%rootid,parallel%comm,mpinfo)
      !>> optimize
      CALL MPI_BCAST(ISTART       ,1,MPI_INTEGER4  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(Nssp         ,1,MPI_INTEGER4  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(step_fixrho         ,1,MPI_INTEGER4  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(Wexict        ,1,MPI_REAL8  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(LDG           ,1,MPI_LOGICAL  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(NHMIX         ,1,MPI_INTEGER4  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(MALPHA        ,1,MPI_REAL8  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(MBETA         ,1,MPI_REAL8  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(NSMIX         ,1,MPI_INTEGER4  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(AMIX       ,1,MPI_REAL8 ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(BMIX       ,1,MPI_REAL8 ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(resta       ,3,MPI_REAL8 ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(W0AM          ,1,MPI_REAL8  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(NDG           ,1,MPI_INTEGER4  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(n_near           ,1,MPI_INTEGER4  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(inpol           ,1,MPI_INTEGER4  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(IOPM          ,1,MPI_INTEGER4  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(IGOAL         ,1,MPI_INTEGER4  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(TIMES         ,1,MPI_REAL8  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(PRESS         ,1,MPI_REAL8  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(TOLF          ,1,MPI_REAL8  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(TOLP          ,1,MPI_REAL8  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(iXCDF         ,1,MPI_INTEGER4  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(ADDcharge     ,1,MPI_INTEGER4  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(Naddstates    ,1,MPI_INTEGER4  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(TOLCG         ,1,MPI_REAL8  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(BLOCK_MBNB    ,1,MPI_INTEGER4  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(ntype ,1,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(hartree_method,1,MPI_INTEGER4  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(CalForce      ,1,MPI_LOGICAL  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(cell_shape    ,1,MPI_INTEGER4 ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(cell_thick    ,1,MPI_REAL8  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(Lcell    ,1,MPI_REAL8  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(Lcellbohr   ,1,MPI_REAL8  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(Lradius_auto  ,1,MPI_LOGICAL  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(IDinit ,1,MPI_INTEGER4 ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(PP_identifer ,1,MPI_INTEGER4 ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(Lgamma           ,1,MPI_LOGICAL  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(fix_xyz ,6,MPI_INTEGER4 ,parallel%rootid,parallel%comm,mpinfo)
      CALL memory_sum("bcast_parameters",real(38,DP)*I4B+&
           &18*DP+14*DP)


   CALL smpi_init_comm(Lpbc)

      !-------set nstates---------
      !> read POSCAR
      CALL read_pos(ntype,cellfile_name)
      !> read local pseudopotential
      CALL read_pspot(ntype,ppfile_name)
      !
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    END SUBROUTINE read_file
    !---------------------------read poscar------------------------
    SUBROUTINE read_pos(nty,filename)
      USE parameters  , ONLY :  elements, Lpbc
      USE math,  ONLY : dir2car,car2dir,det,change_case &
           &, atom_mass,atom_effcharge
      USE struct_module
      USE m_time_evaluate, only: memory_sum
      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B),INTENT(IN)     :: nty
      CHARACTER(LEN=*),INTENT(IN) :: filename
      !



      INTEGER(I4B)  :: i,j,natom_3,nty_4,nty_3

      INTEGER(I4B)  :: ele_n(nty),ele_id(nty+1)
      INTEGER(I4B)  :: filestatu
      REAL(DP)      :: lat_ratio
      REAL(DP)      :: lat_read(3,3)
      CHARACTER(1)  :: chr
      CHARACTER(30) :: ctmp
      LOGICAL       :: ldir
      !
      call memory_sum("read_pos_local",real(nty,DP)*I4B+(nty+1)*I4B+9*DP)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      IF(parallel%isroot)THEN

         OPEN(101,FILE=filename,IOSTAT=filestatu)
         IF(filestatu/=0)THEN
            WRITE(6,*) 'STOP:Could not open',filename
            STOP
         ENDIF
         READ(101,*)          !title,no use
         READ(101,*) lat_ratio
         !read lattice matrix
         DO i=1,3
            READ(101,*) lat_read(:,i)
         ENDDO
         !turn to atomic unit
         lat_mat(:,:)=lat_read(:,:)*lat_ratio*ang2bohr
         !read over poscar
         DO
            READ(101,'(A)') ctmp
            j=INDEX(ctmp,'Dir')
            IF(j/=0)THEN
               ldir=.TRUE.
               EXIT
            ENDIF
            j=INDEX(ctmp,'Car')
            IF(j/=0)THEN
               ldir=.false.
               EXIT
            ENDIF
            READ(ctmp,*) chr
            IF(LGE(chr,'0').AND.LLE(chr,'9'))THEN
               BACKSPACE(101)
               READ(101,*) ele_n(1:nty)
            ENDIF
         ENDDO
         !store id of every elements
         ele_id(1)=1
         DO i=2,nty+1
            ele_id(i)=SUM(ele_n(1:i-1))+1
         ENDDO
         !total atom
         natom=ele_id(nty+1)-1
         !creat structural data
         CALL creat_struct(nty,natom)
         struct%eleid(:)=ele_id(:)
         struct%nati(:)=ele_n(:)
         ! store name
         DO i=1,nty
            CALL change_case(elements(i)(1:1),elements(i)(1:1),1)
         ENDDO
         struct%elements(1:nty)=elements(1:nty)(1:3)
         ! calculate mass and zeta for atom
         DO i=1,nty
            CALL atom_mass( struct%elements(i), struct%mass(i) )
            CALL atom_effcharge( struct%elements(i), &
                 & struct%Lmax(i),struct%prinq(:,i),struct%zeta(:,i))
         ENDDO
         !
         IF (ldir) THEN
            DO i=1,natom
               READ(101,*) struct%pos(:,i)
            ENDDO
            CALL dir2car(struct%pos,struct%poscar,lat_mat)
         ELSE
            DO i=1,natom
               READ(101,*) struct%poscar(:,i)
            ENDDO
            struct%poscar=struct%poscar*ang2bohr
            CALL car2dir(struct%poscar,struct%pos,lat_mat)
         ENDIF
         IF (.NOT.Lpbc)THEN
            CALL ResetLattice()
         ENDIF
         !
         !print*,'lat'
         !print*,lat_mat
         !print*,'pos'
         !print*,recip_lat
         ! print*,'poscar'
         ! print*,struct%poscar
         ! STOP
         CLOSE(101)

      ENDIF !> if parallel%isroot
      CALL MPI_BCAST(natom,1,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
      IF(.NOT.parallel%isroot)THEN
         CALL creat_struct(nty,natom)
      ENDIF
      ! CALL MPI_BCAST(struct%Zion ,nty  ,MPI_INTEGER4,parallel.rootid,parallel.comm,mpinfo)
      CALL MPI_BCAST(struct%nati ,nty  ,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(struct%eleid,nty+1,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
      natom_3=3*natom
      CALL MPI_BCAST(struct%pos  ,natom_3,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(struct%poscar,natom_3,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(struct%mass ,nty  ,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
      nty_4=4*nty
      nty_3=3*nty
      CALL MPI_BCAST(struct%zeta ,nty_4 ,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(struct%prinq,nty_4 ,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(struct%Lmax ,nty   ,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(struct%elements,nty_3,MPI_CHARACTER,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(lat_mat,9,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      call memory_sum("read_pos_local",real(nty,DP)*I4B+(nty+1)*I4B+9*DP)
    ENDSUBROUTINE read_pos
    !-----------------------PARTING-LINE----------------
    SUBROUTINE ResetLattice()
      USE parameters ,    ONLY: Lcell,IsoRmax,RadiusMax
      USE struct_module , ONLY: lat_mat, struct, natom
      IMPLICIT NONE
      INTEGER(I4B)           :: i
      REAL(DP)               :: Ave_Poscar(3) !=0.d0
      Ave_Poscar(:)=0.d0
      lat_mat=0.d0
      DO i=1,3,1
         lat_mat(i,i)=Lcell*ang2bohr
      ENDDO
      DO i=1,3,1
         Ave_Poscar(i)=SUM(struct%poscar(i,:))
      ENDDO
      Ave_Poscar(:)=Ave_Poscar(:)/natom
      DO i=1,3,1
         struct%poscar(i,:)=struct%poscar(i,:)-Ave_Poscar(i)
      ENDDO
      struct%poscar=struct%poscar+Lcell*ang2bohr*0.5d0
      struct%pos=struct%poscar/(Lcell*ang2bohr)
      !CALL out_CONCAR('new-pos')
    ENDSUBROUTINE ResetLattice
    !-----------------------read recip popt------------------------
    SUBROUTINE read_pspot_atom(Ity,filename,ps)
      USE parameters , ONLY : LRadRho
      USE MathSplines
      USE pspot_module , ONLY : pspot
      USE math , ONLY : dfdr
      IMPLICIT NONE
      !IN
      INTEGER(I4B),INTENT(IN)  :: Ity
      CHARACTER(30),INTENT(IN) :: filename
      !OUT
      TYPE(pspot),INTENT(OUT) :: ps
      !LOCAL
      REAL(DP) :: COARSE, MEDIUM, FINE, PRECISE, EXTREME
      CHARACTER(80) :: line
      CHARACTER(20) :: word
      INTEGER(I4B)  :: unit_recpot  , &
           &    filestatu    ,   &
           &    version(3)  ,  &
           &    Id
      INTEGER(I4B)  :: i,j,k,l,m,q,   &
           &     count_points , &
           &     num_points   , &
           &     ncp_num_projectors ,&
           &     nproj
      !---temp array
      REAL(DP),ALLOCATABLE :: dtemp(:,:,:) &
           &,  knots(:)
      INTEGER(I4B),DIMENSION(:),ALLOCATABLE :: nbl,nnn
      REAL(DP),DIMENSION(:,:),ALLOCATABLE :: ncp_D0
      INTEGER(I4B),DIMENSION(:),ALLOCATABLE :: ncp_projector_l
      REAL(DP),DIMENSION(:,:),ALLOCATABLE :: ncp_beta
      !---
      REAL(DP) :: ecut,ionic_charge &
           &,   fact1 ,factor
      LOGICAL :: lopen,lcount
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      PRINT*,'[Reading PseudoPP FILE]',filename
      ! Figer out what kind of file
      Id=INDEX(filename,".",back=.TRUE.)
      IF(TRIM(filename(Id:))/=".rrpot")THEN
         WRITE(6,*) "read_pspot_atom:file must be .rrpot"
         STOP
      ENDIF
      !Inquire file
      unit_recpot=Ity+1500
      !
      DO
         INQUIRE( UNIT=unit_recpot, OPENED=lopen )
         IF ( lopen ) THEN
            unit_recpot = unit_recpot + 1
         ELSE
            EXIT
         ENDIF
      ENDDO
      !open the file
      OPEN(UNIT=unit_recpot,FILE=filename,IOSTAT=filestatu)

      IF (filestatu/=0) THEN
         WRITE(6,*) "Could not open ", TRIM(filename),"."
         WRITE(6,*) "This file does not exit "
         STOP
      ENDIF

      !======================read comment=========================
      line=''
      DO WHILE(INDEX(line,'END COMMENT')==0)
         READ(unit_recpot,'(a)',err=100) line
         IF(  INDEX(line,'COARSE')/=0 .OR. &
              &   INDEX(line,'MEDIUM')/=0 .OR. &
              &   INDEX(line,'FINE')/=0          ) THEN
            !----------------------------------------
            READ(line,*,err=1,end=1) ecut,word
            !PRINT*,'Ecut:',word,ecut,'eV'
            SELECT CASE(trim(word))
            CASE('COARSE') ; COARSE = ecut
            CASE('MEDIUM') ; MEDIUM = ecut
            CASE('FINE')   ; FINE   = ecut
            END SELECT
            !----------------------------------------
            ! * Set PRECISE and EXTREME (like castep)
            PRECISE = 1.2d0 * FINE
            EXTREME = 1.6d0 * FINE
         ENDIF
1        CONTINUE
      ENDDO
      !we should add the accuracy,waiting...
      !
      !
      BACKSPACE(unit_recpot)
      BACKSPACE(unit_recpot)
      READ(unit_recpot,*) ps%Zion,ps%rcut
      READ(unit_recpot,*)


      ! * Read the version
      READ(unit_recpot,*,err=100,end=100) version(1),version(2),version(3)

      !IF((version(1)/=3).OR.(version(2)/=5 .AND.version(2)/=6)) THEN
      IF((version(1)/=3).OR.(version(2)/=5).OR.(version(3)/=1)) THEN
         WRITE(6,*) 'read_pspot_atom: expecting a version 3.5.1 .rrpot file'
         STOP
      ENDIF

      ! * Read the qmax

      READ(unit_recpot,*,err=100,end=100) ps%qmax,ps%rmax

      ! * Find out the number of projectors, and g-grid points

      lcount = .true.
      num_points     = 0
      ncp_num_projectors = 0
      line=''
      DO WHILE((len_trim(adjustl(line))/=4).and.(ncp_num_projectors<2*(3+1)))

         READ(unit_recpot,'(a)',err=100,end=100) line

         IF(len_trim(adjustl(line))==1)THEN
            ncp_num_projectors = ncp_num_projectors+1
            lcount = .false.
         ENDIF

         IF(lcount)THEN
            DO i=1,len(line)
               IF(line(i:i)=='.') num_points = num_points + 1
            ENDDO
         ENDIF

      ENDDO
      !store
      ps%numps=num_points
      !==========================local part===========================
      !array to store lpp
      ALLOCATE(ps%Vlocq(num_points))
      ! ** Read norm conserving pseudopotential into local arrays

      ! * Back to the top of the file

      REWIND(unit_recpot)

      ! * Get to the top of the pseudopotential proper

      line = ''
      DO WHILE(index(line,'END COMMENT')==0)
         READ(unit_recpot,'(a)',err=100) line  ! To the end of the comments
      ENDDO

      DO i=1,2
         READ(unit_recpot,'(a)',err=100) line  ! To the end of the header
      ENDDO

      ! * Read the local component of the pseudopotential

      READ(unit_recpot,*,err=100,end=100)  (ps%Vlocq(i),i=1,ps%numps)

      ! * Calculate the ionic charge by looking at the g -> 0 part
      !   of the local potential

      !azero = io_atomic_to_unit(1.0_dp,'ang')
      !hart2ev = io_unit_to_atomic(io_atomic_to_unit(1.0_dp,'ev'),'hartree')

      fact1 = 4.d0*pi
      !fact1 = hart2ev*bohr2ang
      ionic_charge = REAL(NINT(-ps%Vlocq(2)*(ps%qmax/REAL(ps%numps-1,DP))**2/fact1),DP)
      !============================
      ! print*,"ps%Vlocq",ps%Vlocq(2)
      ! print*,"ps%qmax",ps%qmax
      ! print*,"ps%numps",ps%numps
      ! print*,"fact1",fact1
      !============================
      ! * Double check ionic charge for this species
      IF(NINT(ps%Zion) /= NINT(ionic_charge)) THEN
         PRINT*,'read_pspot_atom: ionic charge mismatch,type:',Ity
         PRINT*,'Zion and ionic_charge',ps%Zion,ionic_charge
         STOP
      ENDIF
      !for check nonlocal part
      ps%nproj=ncp_num_projectors
      !=========================nonlocal part==============================
      IF(ps%nproj>0)THEN
         !allocate
         ALLOCATE(dtemp(ncp_num_projectors,ncp_num_projectors,0:3))
         ALLOCATE(nbl(0:3))
         ALLOCATE(nnn(ncp_num_projectors))

         ! * Read the non-local component of the pseudopotential
         ALLOCATE(ncp_projector_l(ncp_num_projectors))
         ALLOCATE(ncp_beta(ps%numps,ncp_num_projectors))

         nbl = 0
         dtemp = 0.d0
         DO k=1,ncp_num_projectors
            READ(unit_recpot,*,end=100,err=100) ncp_projector_l(k)
            nbl(ncp_projector_l(k)) = nbl(ncp_projector_l(k)) + 1
            nnn(k) = nbl(ncp_projector_l(k))
            READ(unit_recpot,*,end=100,err=100) &
                 &  (dtemp(nnn(k),i,ncp_projector_l(k)),i=1,nnn(k))
            READ(unit_recpot,*,end=100,err=100) (ncp_beta(i,k),i=1,ps%numps)
         ENDDO
         ! * Fill in other half of dtemp
         DO l=0,maxval(ncp_projector_l(1:ncp_num_projectors))
            DO i=1,nbl(l)
               DO j=1,i-1
                  dtemp(j,i,l) = dtemp(i,j,l)
               ENDDO
            ENDDO
         ENDDO

         ! * Convert dtemp to ncp_D0

         ALLOCATE(ncp_D0(ncp_num_projectors,ncp_num_projectors))

         DO i=1,ncp_num_projectors

            l = ncp_projector_l(i)

            DO j=1,ncp_num_projectors
               IF( ncp_projector_l(j) == l ) THEN
                  ncp_D0(i,j) = dtemp(nnn(i),nnn(j),l)
               ELSE
                  ncp_D0(i,j) = 0.0_dp
               ENDIF
            ENDDO

         ENDDO

      ENDIF
      !===============Radial Charge Density==================
      READ(unit_recpot,'(a)',err=100) line
      IF(trim(adjustl(line))/='1000')THEN
         LRadRho=.TRUE.
         ALLOCATE(ps%denr(ps%numps))
         READ(unit_recpot,*,end=100,err=100) (ps%denr(i),i=1,ps%numps)
         !
         READ(unit_recpot,'(a)',err=100) line
         IF(trim(adjustl(line))/='1000')THEN
            WRITE(6,*) 'read_pspot_atom:Problem reading psp file at ending'
            STOP
         ENDIF
      ENDIF
      !======================================================
      IF(ncp_num_projectors > 2*(3+1)) THEN
         WRITE(6,*) 'read_pspot_atom:Check the number of porjectors'
         STOP
      ENDIF

      CLOSE(unit_recpot)
      !============================Convert Pseudopotential=============================
      ! ************************************************
      ! ** Convert pseudopotential into required form **
      ! ************************************************
      !factor = 4.d0*pi*io_atomic_to_unit(1.0_dp,'eV')*io_atomic_to_unit(1.0_dp,'ang')**3
      !ps%Vlocq(:) = ps%Vlocq(:)/factor
      IF(ps%nproj>0)THEN
         ! * Set the number of projectors, and angular momentum indexing
         nproj=0
         DO k=1,ncp_num_projectors
            nproj=nproj+2*ncp_projector_l(k)+1
         ENDDO

         ALLOCATE(ps%proj_l(nproj))
         ALLOCATE(ps%proj_m(nproj))

         ps%nproj = 0
         DO k=1,ncp_num_projectors
            ps%nproj = ps%nproj+2*ncp_projector_l(k)+1 ! 2l+1
            ps%proj_l(ps%nproj-2*ncp_projector_l(k):ps%nproj) = ncp_projector_l(k)
            ps%proj_m(ps%nproj-2*ncp_projector_l(k):ps%nproj) = &
                 (/(m,m=-ncp_projector_l(k),ncp_projector_l(k))/)
         ENDDO

         ! * Set up some indexing
         DEALLOCATE(nnn)
         ALLOCATE(nnn(ps%nproj))

         m = 0
         DO k=1,ncp_num_projectors
            m = m+2*ncp_projector_l(k)+1 ! 2l+1
            nnn(m-2*ncp_projector_l(k):m) = k
         ENDDO

         ! * Set ps%q and ps%D0
         ALLOCATE(ps%D0(nproj,nproj))

         DO k=1,ps%nproj
            DO m=1,ps%nproj

               IF((ps%proj_l(k)==ps%proj_l(m)).AND.&
                    &   (ps%proj_m(k)==ps%proj_m(m))) THEN

                  IF(ABS(ncp_D0(nnn(k),nnn(m))) > 0.0_dp) THEN
                     !ps%D0(n,m) =io_unit_to_atomic(1.0_dp/ncp_D0(nnn(n),nnn(m)),'eV')
                     ps%D0(k,m) =1.0_dp/  ncp_D0(nnn(k),nnn(m))
                  ELSE
                     ps%D0(k,m) = 0.0_dp
                  ENDIF

               ELSE

                  ps%D0(k,m) = 0.0_dp

               ENDIF

            ENDDO
         ENDDO
         ! * Set the beta functions
         ALLOCATE(ps%beta_r(ps%numps,ps%nproj))
         DO k=1,ps%nproj
            ps%beta_r(1:ps%numps,k) = ncp_beta(1:ps%numps,nnn(k))
         ENDDO
      ENDIF
      !=====================parpare to interplote===============
      ALLOCATE(knots(ps%numps))
      DO i=1,ps%numps
         knots(i)=REAL(i,DP)
      ENDDO
      ps%qspacing=ps%qmax/REAL((ps%numps-1),DP)
      !local part
      ALLOCATE(ps%VlocqS(ps%numps))
      ALLOCATE(ps%ddVl_dq2(ps%numps))
      !subtract 1/r in q-space
      ps%VlocqS(1)=ps%Vlocq(1)
      ps%VlocqS(2:)=ps%Vlocq(2:)+ 4.d0*pi*ps%Zion/(ps%qspacing*(knots(2:)-1))**2
      !derivetive 2 order
      CALL spline_cubic_set ( ps%numps, knots, ps%VlocqS,         &
           &    1, 0._DP, 1, 0._DP, ps%ddVl_dq2)
      !nonlocal part
      IF(ps%nproj>0)THEN
         ps%rspacing=ps%rmax/REAL((ps%numps-1),DP)
         !
         !ALLOCATE(ps%ddbeta_dr2(ps%numps,ps%nproj))
         ALLOCATE(ps%dbeta_dr(ps%numps,ps%nproj))
         DO i=1,ps%nproj
            !CALL spline_cubic_set ( ps%numps, knots, ps%beta_r(:,i),      &
            !              &    1, 0._DP, 1, 0._DP, ps%ddbeta_dr2(:,i))
            CALL dfdr(ps%numps,ps%rspacing,ps%beta_r(:,i),ps%dbeta_dr(:,i))
            !xqtest
            !OPEN(1111,FILE='dfdr.dat')
            !   DO j=1,ps%numps
            !      WRITE(1111,*) (j-1)*ps%rspacing,ps%dbeta_dr(j,1)*2
            !   ENDDO
            !CLOSE(1111)
            !STOP
         ENDDO
      ENDIF
      !charge density
      IF(LRadRho)THEN
         ALLOCATE(ps%ddden_dr2(ps%numps))
         CALL spline_cubic_set ( ps%numps, knots, ps%denr,      &
              &    1, 0._DP, 1, 0._DP, ps%ddden_dr2)
      ENDIF
      !========================End========================
      !destroy>>>
      IF(ALLOCATED(dtemp))           DEALLOCATE(dtemp)
      IF(ALLOCATED(knots))           DEALLOCATE(knots)
      IF(ALLOCATED(nbl))             DEALLOCATE(nbl)
      IF(ALLOCATED(nnn))             DEALLOCATE(nnn)
      IF(ALLOCATED(ncp_D0))          DEALLOCATE(ncp_D0)
      IF(ALLOCATED(ncp_projector_l)) DEALLOCATE(ncp_projector_l)
      IF(ALLOCATED(ncp_beta))        DEALLOCATE(ncp_beta)
      !<<<destroy
      RETURN
100   WRITE(6,*) 'read_pspot_atom:error'
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ENDSUBROUTINE read_pspot_atom
    !--------------------------read_pspot--------------------------
    SUBROUTINE read_pspot(nty,filenames)
      USE parameters , ONLY : Nstates,Nstates_global, Nspin,Naddstates,Lpbc,ADDcharge,PP_identifer, &
           & LRadRho
      USE pspot_module
      USE struct_module, ONLY:struct,ncharge
      USE m_time_evaluate, ONLY: memory_sum
      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B),INTENT(IN)  :: nty
      CHARACTER(30),INTENT(IN) :: filenames(nty)
      !LOCAL
      INTEGER(I4B) :: Ity
      INTEGER(I4B) :: Ipj,l,m , i  ,mnp
      INTEGER(I4B) :: N_BetaR,N_D0
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !for every speice atom
      ALLOCATE(psp(nty))
      !

      IF(parallel%isroot)THEN

         DO Ity=1,nty
            IF( PP_identifer==0 .AND. Lpbc )THEN
               CALL read_pspot_atom(Ity,filenames(Ity),psp(Ity))
            ELSEIF( PP_identifer==0 .AND. (.NOT. Lpbc) )THEN
               CALL read_realpot_atom(Ity,filenames(Ity),psp(Ity))
            ELSEIF( PP_identifer==1 )THEN
               CALL read_upf(Ity,filenames(Ity),psp(Ity))
            ENDIF
            struct%Zion(Ity)=psp(Ity)%Zion
         ENDDO
!#ifdef DEBUG
! #ifdef 1
!          if(parallel%isroot)then
! #endif
!          print*,'psp%D0',psp(1)%D0
! #ifdef 1
!          endif
! #endif
!#endif


      ENDIF !> parallel%isroot
      !> parallel broadcast
      if(lpbc)then
         if(PP_identifer==0)then
            DO Ity=1,nty
               CALL MPI_BCAST(struct%Zion(Ity) ,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
               CALL MPI_BCAST(psp(Ity)%Zion ,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
               CALL MPI_BCAST(psp(Ity)%rcut ,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
               CALL MPI_BCAST(psp(Ity)%qmax ,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
               CALL MPI_BCAST(psp(Ity)%rmax ,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
               CALL MPI_BCAST(psp(Ity)%numps ,1,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
               IF(.NOT.ALLOCATED(psp(Ity)%Vlocq))THEN
                  ALLOCATE(psp(Ity)%Vlocq(psp(Ity)%numps))
                  call memory_sum('ps_Vlocq',real(psp(Ity)%numps,DP)*DP)
               ENDIF
               IF(.NOT.ALLOCATED(psp(Ity)%denr))THEN
                  ALLOCATE(psp(Ity)%denr(psp(Ity)%numps))
                  call memory_sum('ps_denr',real(psp(Ity)%numps,DP)*DP)
               ENDIF
               CALL MPI_BCAST(psp(Ity)%Vlocq ,psp(Ity)%numps,&
                    & MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
               CALL MPI_BCAST(psp(Ity)%denr,psp(Ity)%numps,&
                    & MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
               CALL MPI_BCAST(psp(Ity)%nproj ,1,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
               CALL MPI_BCAST(LRadRho ,1,MPI_LOGICAL,parallel%rootid,parallel%comm,mpinfo)
               IF(LRadRho)THEN
                  ! print *,psp(Ity)%numps_den,parallel.myid
                  IF(.NOT.ALLOCATED(psp(Ity)%denr))THEN
                     allocate(psp(Ity)%denr(psp(Ity)%numps))
                     call memory_sum('ps_denr',real(psp(Ity)%numps_den,DP)*DP)
                  ENDIF
                  CALL MPI_BCAST(psp(Ity)%denr ,psp(Ity)%numps,&
                       MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
               ENDIF
               !         print *,"ccc",parallel.myid
               ! if (allocated(psp(ity)%proj_l)) deallocate (psp(ity)%proj_l)
               IF(.NOT.ALLOCATED(psp(Ity)%proj_l))THEN
                  ALLOCATE(psp(Ity)%proj_l(psp(Ity)%nproj))
                  call memory_sum('ps_proj_l',real(psp(Ity)%nproj,DP)*I4B)
               ENDIF
               !psp(ity)%proj_l = 2 !longlongago,debug for communction used it,which induced by bcasted variable not allocated all
               IF(.NOT.ALLOCATED(psp(Ity)%proj_m))THEN
                  ALLOCATE(psp(Ity)%proj_m(psp(Ity)%nproj))
                  call memory_sum('ps_proj_m',real(psp(Ity)%nproj,DP)*I4B)
               ENDIF
               IF(.NOT.ALLOCATED(psp(Ity)%D0))THEN
                  ALLOCATE(psp(Ity)%D0(psp(Ity)%nproj,psp(Ity)%nproj))
                  call memory_sum('ps_r_real',real(size(psp(Ity)%D0),DP)*DP)
               ENDIF
               !print *,'debug 20181208',size(psp),size(psp(ity)%proj_l),ity,parallel%myid,psp(Ity)%nproj,'ll',psp(ity)%proj_l
               ! CALL MPI_Barrier(parallel%comm,mpinfo)
               CALL MPI_BCAST(psp(Ity)%proj_l,psp(Ity)%nproj,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
               ! CALL MPI_BCAST(psp(Ity)%proj_l(1),1,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
               ! CALL MPI_Barrier(parallel%comm,mpinfo)
               CALL MPI_BCAST(psp(Ity)%proj_m ,psp(Ity)%nproj,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
               n_D0=psp(Ity)%nproj**2
               CALL MPI_BCAST(psp(Ity)%D0     ,n_D0,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
               IF(.NOT.ALLOCATED(psp(Ity)%beta_r)) ALLOCATE(psp(Ity)%beta_r(psp(Ity)%numps,psp(Ity)%nproj))
               N_BetaR=psp(Ity)%numps*psp(Ity)%nproj
               CALL MPI_BCAST(psp(Ity)%beta_r ,N_BetaR,&
                    MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
               call MPI_BCAST(psp(Ity)%qspacing,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
               call MPI_BCAST(psp(Ity)%rspacing,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
               if(.not.allocated(psp(Ity)%VlocqS))allocate(psp(Ity)%VlocqS(psp(Ity)%numps))
               if(.not.allocated(psp(Ity)%ddVl_dq2))allocate(psp(Ity)%ddVl_dq2(psp(Ity)%numps))
               if(.not.allocated(psp(Ity)%dbeta_dr))allocate(psp(Ity)%dbeta_dr(psp(Ity)%numps,psp(Ity)%nproj))
               if(.not.allocated(psp(Ity)%ddden_dr2))allocate(psp(Ity)%ddden_dr2(psp(Ity)%numps))
               call MPI_BCAST(psp(Ity)%VlocqS,psp(Ity)%numps,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
               call MPI_BCAST(psp(Ity)%ddVl_dq2,psp(Ity)%numps,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
               call MPI_BCAST(psp(Ity)%dbeta_dr,N_BetaR,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
               call MPI_BCAST(psp(Ity)%ddden_dr2,psp(Ity)%numps,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
               if(.not.allocated(psp(Ity)%r_real))allocate(psp(Ity)%r_real(psp(Ity)%numps))
               call MPI_BCAST(psp(Ity)%r_real,psp(Ity)%numps,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
            ENDDO
         else
            DO Ity=1,nty
               CALL MPI_BCAST(struct%Zion(Ity) ,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
               CALL MPI_BCAST(psp(Ity)%Zion ,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
               CALL MPI_BCAST(psp(Ity)%rcut ,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
               CALL MPI_BCAST(psp(Ity)%qmax ,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
               CALL MPI_BCAST(psp(Ity)%rmax ,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
               CALL MPI_BCAST(psp(Ity)%numps ,1,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
               CALL MPI_BCAST(psp(Ity)%qnumps ,1,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
               IF(.NOT.ALLOCATED(psp(Ity)%Vlocq))THEN
                  ALLOCATE(psp(Ity)%Vlocq(psp(Ity)%qnumps))
                  call memory_sum('ps_Vlocq',real(psp(Ity)%numps,DP)*DP)
               ENDIF
               IF(.NOT.ALLOCATED(psp(Ity)%denr))THEN
                  ALLOCATE(psp(Ity)%denr(psp(Ity)%numps))
                  call memory_sum('ps_denr',real(psp(Ity)%numps,DP)*DP)
               ENDIF
               CALL MPI_BCAST(psp(Ity)%Vlocq ,psp(Ity)%qnumps,&
                    & MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
               CALL MPI_BCAST(psp(Ity)%denr,psp(Ity)%numps,&
                    & MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
               CALL MPI_BCAST(psp(Ity)%nproj ,1,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
               CALL MPI_BCAST(LRadRho ,1,MPI_LOGICAL,parallel%rootid,parallel%comm,mpinfo)
               IF(LRadRho)THEN
                  ! print *,psp(Ity)%numps_den,parallel.myid
                  IF(.NOT.ALLOCATED(psp(Ity)%denr))THEN
                     allocate(psp(Ity)%denr(psp(Ity)%numps))
                     call memory_sum('ps_denr',real(psp(Ity)%numps_den,DP)*DP)
                  ENDIF
                  CALL MPI_BCAST(psp(Ity)%denr ,psp(Ity)%numps,&
                       MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
               ENDIF
               !         print *,"ccc",parallel.myid
               ! if (allocated(psp(ity)%proj_l)) deallocate (psp(ity)%proj_l)
               IF(.NOT.ALLOCATED(psp(Ity)%proj_l))THEN
                  ALLOCATE(psp(Ity)%proj_l(psp(Ity)%nproj))
                  call memory_sum('ps_proj_l',real(psp(Ity)%nproj,DP)*I4B)
               ENDIF
               !psp(ity)%proj_l = 2 !longlongago,debug for communction used it,which induced by bcasted variable not allocated all
               IF(.NOT.ALLOCATED(psp(Ity)%proj_m))THEN
                  ALLOCATE(psp(Ity)%proj_m(psp(Ity)%nproj))
                  call memory_sum('ps_proj_m',real(psp(Ity)%nproj,DP)*I4B)
               ENDIF
               IF(.NOT.ALLOCATED(psp(Ity)%D0))THEN
                  ALLOCATE(psp(Ity)%D0(psp(Ity)%nproj,psp(Ity)%nproj))
                  call memory_sum('ps_r_real',real(size(psp(Ity)%D0),DP)*DP)
               ENDIF
               !print *,'debug 20181208',size(psp),size(psp(ity)%proj_l),ity,parallel%myid,psp(Ity)%nproj,'ll',psp(ity)%proj_l
               ! CALL MPI_Barrier(parallel%comm,mpinfo)
               CALL MPI_BCAST(psp(Ity)%proj_l,psp(Ity)%nproj,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
               ! CALL MPI_BCAST(psp(Ity)%proj_l(1),1,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
               ! CALL MPI_Barrier(parallel%comm,mpinfo)
               CALL MPI_BCAST(psp(Ity)%proj_m ,psp(Ity)%nproj,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
               n_D0=psp(Ity)%nproj**2
               CALL MPI_BCAST(psp(Ity)%D0     ,n_D0,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
               IF(.NOT.ALLOCATED(psp(Ity)%beta_r)) ALLOCATE(psp(Ity)%beta_r(psp(Ity)%numps,psp(Ity)%nproj))
               N_BetaR=psp(Ity)%numps*psp(Ity)%nproj
               CALL MPI_BCAST(psp(Ity)%beta_r ,N_BetaR,&
                    MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
               IF(.NOT.ALLOCATED(psp(Ity)%dbeta_dr))THEN
                  ALLOCATE(psp(Ity)%dbeta_dr(psp(Ity)%qnumps,psp(Ity)%nproj))
                  call memory_sum('ps_dbeta_dr',real(psp(Ity)%qnumps,DP)*DP)
               ENDIF
               CALL MPI_BCAST(psp(Ity)%dbeta_dr ,psp(Ity)%qnumps*psp(Ity)%nproj,&
                    MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
               call MPI_BCAST(psp(Ity)%qspacing,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
               call MPI_BCAST(psp(Ity)%rspacing,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
               if(.not.allocated(psp(Ity)%VlocqS))allocate(psp(Ity)%VlocqS(psp(Ity)%qnumps))
               if(.not.allocated(psp(Ity)%ddVl_dq2))allocate(psp(Ity)%ddVl_dq2(psp(Ity)%qnumps))
               call MPI_BCAST(psp(Ity)%VlocqS,psp(Ity)%qnumps,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
               call MPI_BCAST(psp(Ity)%ddVl_dq2,psp(Ity)%qnumps,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
               if(.not.allocated(psp(Ity)%r_real))allocate(psp(Ity)%r_real(psp(Ity)%numps))
               call MPI_BCAST(psp(Ity)%r_real,psp(Ity)%numps,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
            ENDDO
         endif
      else
         DO Ity=1,nty
            CALL MPI_BCAST(struct%Zion(Ity) ,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
            CALL MPI_BCAST(psp(Ity)%Zion ,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
            CALL MPI_BCAST(psp(Ity)%rcut ,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
            CALL MPI_BCAST(psp(Ity)%qmax ,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
            CALL MPI_BCAST(psp(Ity)%rmax ,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
            CALL MPI_BCAST(psp(Ity)%numps ,1,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
            IF(.NOT.ALLOCATED(psp(Ity)%V_loc))THEN
               ALLOCATE(psp(Ity)%V_loc(psp(Ity)%numps))
               call memory_sum('ps_V_loc',real(psp(Ity)%numps,DP)*DP)
            ENDIF
            IF(.NOT.ALLOCATED(psp(Ity)%r_real))THEN
               ALLOCATE(psp(Ity)%r_real(psp(Ity)%numps))
               call memory_sum('ps_r_real',real(psp(Ity)%numps,DP)*DP)
            ENDIF
            CALL MPI_BCAST(psp(Ity)%V_loc ,psp(Ity)%numps,&
                 & MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
            CALL MPI_BCAST(psp(Ity)%r_real,psp(Ity)%numps,&
                 & MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
            CALL MPI_BCAST(psp(Ity)%nproj ,1,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
            CALL MPI_BCAST(LRadRho ,1,MPI_LOGICAL,parallel%rootid,parallel%comm,mpinfo)
            IF(LRadRho)THEN
               ! print *,psp(Ity)%numps_den,parallel.myid
               CALL MPI_BCAST(psp(Ity)%numps_den ,1,&
                    MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
               IF(.NOT.ALLOCATED(psp(Ity)%denr))THEN
                  allocate(psp(Ity)%denr(psp(Ity)%numps_den))
                  call memory_sum('ps_denr',real(psp(Ity)%numps_den,DP)*DP)
               ENDIF
               CALL MPI_BCAST(psp(Ity)%denr ,psp(Ity)%numps_den,&
                    MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
            ENDIF
            !         print *,"ccc",parallel.myid
            ! if (allocated(psp(ity)%proj_l)) deallocate (psp(ity)%proj_l)
            IF(.NOT.ALLOCATED(psp(Ity)%proj_l))THEN
               ALLOCATE(psp(Ity)%proj_l(psp(Ity)%nproj))
               call memory_sum('ps_proj_l',real(psp(Ity)%nproj,DP)*I4B)
            ENDIF
            !psp(ity)%proj_l = 2 !longlongago,debug for communction used it,which induced by bcasted variable not allocated all
            IF(.NOT.ALLOCATED(psp(Ity)%proj_m))THEN
               ALLOCATE(psp(Ity)%proj_m(psp(Ity)%nproj))
               call memory_sum('ps_proj_m',real(psp(Ity)%nproj,DP)*I4B)
            ENDIF
            IF(.NOT.ALLOCATED(psp(Ity)%D0))THEN
               ALLOCATE(psp(Ity)%D0(psp(Ity)%nproj,psp(Ity)%nproj))
               call memory_sum('ps_r_real',real(size(psp(Ity)%D0),DP)*DP)
            ENDIF
            !print *,'debug 20181208',size(psp),size(psp(ity)%proj_l),ity,parallel%myid,psp(Ity)%nproj,'ll',psp(ity)%proj_l
            ! CALL MPI_Barrier(parallel%comm,mpinfo)
            CALL MPI_BCAST(psp(Ity)%proj_l,psp(Ity)%nproj,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
            ! CALL MPI_BCAST(psp(Ity)%proj_l(1),1,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
            ! CALL MPI_Barrier(parallel%comm,mpinfo)
            CALL MPI_BCAST(psp(Ity)%proj_m ,psp(Ity)%nproj,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
            n_D0=psp(Ity)%nproj**2
            CALL MPI_BCAST(psp(Ity)%D0     ,n_D0,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
            IF(.NOT.ALLOCATED(psp(Ity)%beta_r)) ALLOCATE(psp(Ity)%beta_r(psp(Ity)%numps,psp(Ity)%nproj))
            N_BetaR=psp(Ity)%numps*psp(Ity)%nproj
            CALL MPI_BCAST(psp(Ity)%beta_r ,N_BetaR,&
                 MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
         ENDDO
      endif

      ncharge=SUM(struct%nati(:)*struct%Zion(:))+ADDcharge
      !states we need
      Nstates=NINT(real(ncharge)*NSPIN/2)+Naddstates

      if(parallel%isroot)print*,'[ Total states used ]',Nstates



      Nstates_global=Nstates

      call array_split(Nstates)
      if(Lpbc)Nstates=parallel%nstate_proc

      !
      max_nproj=MAXVAL(psp(:)%nproj)
      max_rcut=MAXVAL(psp(:)%rcut)
      !
      mnp=MAXVAL(psp(:)%numps)
      ALLOCATE(tknots(mnp))
      DO i=1,mnp
         tknots(i)=REAL(i,DP)
      ENDDO
      call memory_sum('tknots',real(size(tknots),DP)*DP)
      !OPEN(10010,FILE='test.dat')
      !  WRITE(10010,*) psp(1)%VlocqS
      !CLOSE(10010)
      !STOP
      !OPEN(10010,FILE='nolocal.dat')
      !   ! DO Ipj=1,psp(1)%nproj
      !   !    l=psp(1)%proj_l(Ipj)
      !   !    m=psp(1)%proj_m(Ipj)
      !   !    WRITE(10010,*) l,m
      !   !    WRITE(10010,*) psp(1)%D0(Ipj,Ipj)
      !   !    WRITE(10010,*) psp(1)%beta_r(:,Ipj)
      !   ! ENDDO
      !   WRITE(10010,*) psp(1)%denr(:)
      !CLOSE(10010)
      !STOP
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ENDSUBROUTINE read_pspot
    !-----------------------read recip popt------------------------
    SUBROUTINE read_realpot_atom(Ity,filename,ps)
      USE parameters , ONLY : LRadRho
      USE MathSplines
      USE pspot_module , ONLY : pspot
      USE math , ONLY : dfdr
      USE m_time_evaluate, ONLY: memory_sum,memory_free
      IMPLICIT NONE
      !IN
      INTEGER(I4B),INTENT(IN)  :: Ity
      CHARACTER(30),INTENT(IN) :: filename
      !OUT
      TYPE(pspot),INTENT(OUT) :: ps
      !LOCAL
      REAL(DP) :: COARSE, MEDIUM, FINE, PRECISE, EXTREME
      CHARACTER(80) :: line
      CHARACTER(20) :: word
      INTEGER(I4B)  :: unit_recpot  , &
           &    filestatu    ,   &
           &    version(3)  ,  &
           &    Id
      INTEGER(I4B)  :: i,j,k,l,m,q,   &
           &     count_points , &
           &     num_points   , &
           &     ncp_num_projectors ,&
           &     nproj
      !---temp array
      REAL(DP),ALLOCATABLE :: dtemp(:,:,:) &
           &,  knots(:)
      INTEGER(I4B),DIMENSION(:),ALLOCATABLE :: nbl,nnn
      REAL(DP),DIMENSION(:,:),ALLOCATABLE :: ncp_D0
      INTEGER(I4B),DIMENSION(:),ALLOCATABLE :: ncp_projector_l
      REAL(DP),DIMENSION(:,:),ALLOCATABLE :: ncp_beta
      !---
      REAL(DP) :: ecut,ionic_charge &
           &,   fact1 ,factor
      LOGICAL :: lopen,lcount
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      PRINT*,'[Reading PseudoPP FILE]',filename
      ! Figer out what kind of file
      Id=INDEX(filename,".",back=.TRUE.)
      IF(TRIM(filename(Id:))/=".realpot")THEN
         WRITE(6,*) "read_pspot_atom:file must be .realpot"
         STOP
      ENDIF
      !Inquire file
      unit_recpot=Ity+1500
      !
      DO
         INQUIRE( UNIT=unit_recpot, OPENED=lopen )
         IF ( lopen ) THEN
            unit_recpot = unit_recpot + 1
         ELSE
            EXIT
         ENDIF
      ENDDO
      !open the file
      OPEN(UNIT=unit_recpot,FILE=filename,IOSTAT=filestatu)

      IF (filestatu/=0) THEN
         WRITE(6,*) "Could not open ", TRIM(filename),"."
         WRITE(6,*) "This file does not exit "
         STOP
      ENDIF

      !======================read comment=========================
      line=''
      DO WHILE(INDEX(line,'END COMMENT')==0)
         READ(unit_recpot,'(a)',err=101) line
         IF(  INDEX(line,'COARSE')/=0 .OR. &
              &   INDEX(line,'MEDIUM')/=0 .OR. &
              &   INDEX(line,'FINE')/=0          ) THEN
            !----------------------------------------
            READ(line,*,err=1,end=1) ecut,word
            !PRINT*,'Ecut:',word,ecut,'eV'
            SELECT CASE(trim(word))
            CASE('COARSE') ; COARSE = ecut
            CASE('MEDIUM') ; MEDIUM = ecut
            CASE('FINE')   ; FINE   = ecut
            END SELECT
            !----------------------------------------
            ! * Set PRECISE and EXTREME (like castep)
            PRECISE = 1.2d0 * FINE
            EXTREME = 1.6d0 * FINE
         ENDIF
1        CONTINUE
      ENDDO
      !we should add the accuracy,waiting...
      !
      !
      BACKSPACE(unit_recpot)
      BACKSPACE(unit_recpot)
      READ(unit_recpot,*) ps%Zion,ps%rcut
      READ(unit_recpot,*)


      ! * Read the version
      READ(unit_recpot,*,err=101,end=101) version(1),version(2),version(3)

      !IF((version(1)/=3).OR.(version(2)/=5 .AND.version(2)/=6)) THEN
      IF((version(1)/=3).OR.(version(2)/=5).OR.(version(3)/=1)) THEN
         WRITE(6,*) 'read_pspot_atom: expecting a version 3.5.1 .rrpot file'
         STOP
      ENDIF

      ! * Read the qmax

      READ(unit_recpot,*,err=101,end=101)ps%qmax, ps%rmax
      !print*,"ll",ps%rmax
      ! * Find out the number of projectors, and g-grid points

      !lcount = .true.
      num_points     = 0
      ncp_num_projectors = 0
      line=''
      !=======================================================
      !##READ ARRAY SIZE
      READ(unit_recpot,*)num_points
      !=======================================================
      !##READ NUMBER OF PROJECTORS
      DO WHILE((len_trim(adjustl(line))/=5).and.(ncp_num_projectors<2*(3+1)))
         READ(unit_recpot,'(a)',err=101,end=101) line
         IF(len_trim(adjustl(line))==1)THEN
            ncp_num_projectors = ncp_num_projectors+1
            !lcount = .false.
         ENDIF
         !IF(lcount)THEN
         !   DO i=1,len(line)
         !      IF(line(i:i)=='.') num_points = num_points + 1
         !   ENDDO
         !ENDIF
      ENDDO
      !store
      ps%numps=num_points
      !==============================
      !##TEST
      ! print*,"ps%numps",ps%numps,ncp_num_projectors
      !==============================
      !==========================local part===========================
      !array to store lpp
      ALLOCATE(ps%V_loc(num_points))
      ALLOCATE(ps%r_real(num_points))
      call memory_sum('ps_local_V_r',real(2,DP)*num_points*DP)
      ! ** Read norm conserving pseudopotential into local arrays

      ! * Back to the top of the file

      REWIND(unit_recpot)

      ! * Get to the top of the pseudopotential proper

      line = ''
      DO WHILE(index(line,'END COMMENT')==0)
         READ(unit_recpot,'(a)',err=101) line  ! To the end of the comments
      ENDDO

      DO i=1,3
         READ(unit_recpot,'(a)',err=101) line  ! To the end of the header
      ENDDO

      ! * Read the local component of the pseudopotential
      READ(unit_recpot,*,err=101,end=101)  (ps%r_real(i),i=1,ps%numps)
      READ(unit_recpot,*,err=101,end=101)  (ps%V_loc(i),i=1,ps%numps)

      ! * Calculate the ionic charge by looking at the g -> 0 part
      !   of the local potential

      !azero = io_atomic_to_unit(1.0_dp,'ang')
      !hart2ev = io_unit_to_atomic(io_atomic_to_unit(1.0_dp,'ev'),'hartree')

      !fact1 = 4.d0*pi
      !!fact1 = hart2ev*bohr2ang
      !ionic_charge = REAL(NINT(-ps%Vlocq(2)*(ps%qmax/REAL(ps%numps-1,DP))**2/fact1),DP)
      !! * Double check ionic charge for this species
      !IF(ps%Zion /= NINT(ionic_charge)) THEN
      !    PRINT*,'read_pspot_atom: ionic charge mismatch,type:',Ity
      !    PRINT*,'Zion and ionic_charge',ps%Zion,ionic_charge
      !    STOP
      !ENDIF
      !for check nonlocal part
      ps%nproj=ncp_num_projectors
      !==============================
      !##TEST
      ! print*,"ps%nproj",ps%nproj
      !==============================
      !##PLOT Vloc {{{
      !open(1119,file="frxyz_0.4")
      !write(1119,*)ps%V_loc
      !close(1119)
      !open(1119,file="rxyz_0.4")
      !write(1119,*)ps%r_real
      !close(1119)
      !##PLOT Vloc }}}
      !==============================
      !##READ NONLOCAL PART
      !=========================nonlocal part==============================
      IF(ps%nproj>0)THEN
         !allocate
         ALLOCATE(dtemp(ncp_num_projectors,ncp_num_projectors,0:3))
         ALLOCATE(nbl(0:3))
         ALLOCATE(nnn(ncp_num_projectors))

         ! * Read the non-local component of the pseudopotential
         ALLOCATE(ncp_projector_l(ncp_num_projectors))
         ALLOCATE(ncp_beta(ps%numps,ncp_num_projectors))
         call memory_sum("ps_temp_array",&
              &(real(size(dtemp),DP)+size(ncp_beta))*DP+&
              &(size(nbl)+size(nnn)+size(ncp_projector_l))*I4B)

         nbl = 0
         dtemp = 0.d0
         ! print*,'norm conserved pseudopotential projectors',ncp_num_projectors
         DO k=1,ncp_num_projectors
            !> the k_th value: l
            READ(unit_recpot,*,end=101,err=101) ncp_projector_l(k)
            !> NumBer of projectors with l, the l map for order of total l to order of different l
            nbl(ncp_projector_l(k)) = nbl(ncp_projector_l(k)) + 1   !> the fucked number equal to one forever, use this to confused everyone
            nnn(k) = nbl(ncp_projector_l(k))                        !> in fact, it was the ith number in D0 in l dimension
            !> a normalized const
            READ(unit_recpot,*,end=101,err=101) &
                 &  (dtemp(nnn(k),i,ncp_projector_l(k)),i=1,nnn(k)) 




            READ(unit_recpot,*,end=101,err=101) (ncp_beta(i,k),i=1,ps%numps)
         ENDDO
         !> in norm conserved pseudopotential dtemp just a monodrome
         ! * Fill in other half of dtemp
         DO l=0,maxval(ncp_projector_l(1:ncp_num_projectors))
            DO i=1,nbl(l)
               DO j=1,i-1
                  dtemp(j,i,l) = dtemp(i,j,l)
               ENDDO
            ENDDO
         ENDDO

         ! * Convert dtemp to ncp_D0

         ALLOCATE(ncp_D0(ncp_num_projectors,ncp_num_projectors))
         call memory_sum("ps_temp_array_ncp_D0",real(size(ncp_D0),DP)*DP)

         DO i=1,ncp_num_projectors

            l = ncp_projector_l(i)

            DO j=1,ncp_num_projectors
               IF( ncp_projector_l(j) == l ) THEN
                  ncp_D0(i,j) = dtemp(nnn(i),nnn(j),l)
                  ! print *,'ncp_D0',ncp_D0(i,j)
               ELSE
                  ncp_D0(i,j) = 0.0_dp
               ENDIF
            ENDDO

         ENDDO

      ENDIF
      !===============Radial Charge Density==================
      READ(unit_recpot,'(a)',err=101) line
      ! print*,line
      IF(trim(adjustl(line))/='10000')THEN
         LRadRho=.TRUE.
         !=======================================================
         !##READ ARRAY SIZE
         READ(unit_recpot,*)num_points
         !ps%numps=num_points
         !=======================================================
         ALLOCATE(ps%denr(num_points))
         call memory_sum("ps_denr",real(size(ps%denr),DP)*DP)
         ps%numps_den=num_points
         READ(unit_recpot,*,end=101,err=101) (ps%denr(i),i=1,num_points)
         !
         READ(unit_recpot,'(a20)',err=123) line
         IF(trim(adjustl(line))/='10000')THEN
            WRITE(6,*) 'read_pspot_atom:Problem reading psp file at ending'
            STOP
         ENDIF
      ELSE
         LRadRho=.FALSE.
      ENDIF
      !======================================================
      IF(ncp_num_projectors > 2*(3+1)) THEN
         WRITE(6,*) 'read_pspot_atom:Check the number of porjectors'
         STOP
      ENDIF
      CLOSE(unit_recpot)
      !============================Convert Pseudopotential=============================
      ! ************************************************
      ! ** Convert pseudopotential into required form **
      ! ************************************************
      !factor = 4.d0*pi*io_atomic_to_unit(1.0_dp,'eV')*io_atomic_to_unit(1.0_dp,'ang')**3
      !ps%Vlocq(:) = ps%Vlocq(:)/factor
      IF(ps%nproj>0)THEN
         ! * Set the number of projectors, and angular momentum indexing
         nproj=0
         DO k=1,ncp_num_projectors
            nproj=nproj+2*ncp_projector_l(k)+1
         ENDDO
         ALLOCATE(ps%proj_l(nproj))
         ALLOCATE(ps%proj_m(nproj))
         call memory_sum("ps_proj",real(2,DP)*nproj*I4B)
         ps%nproj = 0
         DO k=1,ncp_num_projectors
            ps%nproj = ps%nproj+2*ncp_projector_l(k)+1 ! 2l+1
            ps%proj_l(ps%nproj-2*ncp_projector_l(k):ps%nproj) = ncp_projector_l(k)
            ps%proj_m(ps%nproj-2*ncp_projector_l(k):ps%nproj) = &
                 (/(m,m=-ncp_projector_l(k),ncp_projector_l(k))/)
         ENDDO
         ! * Set up some indexing
         call memory_free('ps_reallocate',real(size(nnn),DP)*I4B)
         DEALLOCATE(nnn)
         ALLOCATE(nnn(ps%nproj))
         call memory_sum('ps_reallocate',real(size(nnn),DP)*I4B)
         m = 0
         DO k=1,ncp_num_projectors
            m = m+2*ncp_projector_l(k)+1 ! 2l+1
            nnn(m-2*ncp_projector_l(k):m) = k
         ENDDO
         ! * Set ps%q and ps%D0
         ! print*,'shenjingbing##',ncp_D0
         ALLOCATE(ps%D0(nproj,nproj))
         call memory_sum('ps_D0',real(size(ps%D0),DP)*DP)
         DO k=1,ps%nproj
            DO m=1,ps%nproj
               IF((ps%proj_l(k)==ps%proj_l(m)).AND.&
                    &   (ps%proj_m(k)==ps%proj_m(m))) THEN
                  IF(ABS(ncp_D0(nnn(k),nnn(m))) > 0.0_dp) THEN
                     !ps%D0(n,m) =io_unit_to_atomic(1.0_dp/ncp_D0(nnn(n),nnn(m)),'eV')
                     ps%D0(k,m) =1.0_dp/  ncp_D0(nnn(k),nnn(m))
                  ELSE
                     ps%D0(k,m) = 0.0_dp
                  ENDIF
               ELSE
                  ps%D0(k,m) = 0.0_dp
               ENDIF
            ENDDO
         ENDDO
         ! * Set the beta functions
         ALLOCATE(ps%beta_r(ps%numps,ps%nproj))
         call memory_sum('ps_beta_r',real(size(ps%beta_r),DP)*I4B)
         DO k=1,ps%nproj
            ps%beta_r(1:ps%numps,k) = ncp_beta(1:ps%numps,nnn(k))
         ENDDO
      ENDIF
      !=====================parpare to interplote===============
      ALLOCATE(knots(num_points))
      call memory_sum('ps_local_knots',real(size(knots),DP)*DP)
      DO i=1,num_points
         knots(i)=REAL(i,DP)
      ENDDO
      ps%rmax=ps%r_real(num_points)
      ps%rspacing=ps%rmax/REAL((num_points-1),DP)
      IF(LRadRho)THEN
         ALLOCATE(ps%ddden_dr2(num_points))
         call memory_sum('ps_ddden_dr2',real(size(ps%ddden_dr2),DP)*DP)
         CALL spline_cubic_set ( num_points, knots, ps%denr,      &
              &    1, 0._DP, 1, 0._DP, ps%ddden_dr2)
      ENDIF
      !========================End========================
      !destroy>>>
      call memory_free('ps_local',&
              &(real(size(dtemp),DP)+size(ncp_beta)+size(knots)&
              &+size(ncp_D0))*DP+&
              &(size(nbl)+size(nnn)+size(ncp_projector_l))*I4B)
      IF(ALLOCATED(dtemp))           DEALLOCATE(dtemp)
      IF(ALLOCATED(knots))           DEALLOCATE(knots)
      IF(ALLOCATED(nbl))             DEALLOCATE(nbl)
      IF(ALLOCATED(nnn))             DEALLOCATE(nnn)
      IF(ALLOCATED(ncp_D0))          DEALLOCATE(ncp_D0)
      IF(ALLOCATED(ncp_projector_l)) DEALLOCATE(ncp_projector_l)
      IF(ALLOCATED(ncp_beta))        DEALLOCATE(ncp_beta)
      RETURN
222   WRITE(6,*) 'read_pspot_atom:error'
123   print*,"err read_mode line",line
101   WRITE(6,*) 'read_pspot_atom:error'
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ENDSUBROUTINE read_realpot_atom
    !>added by Lantian Xue < 9 aug 2018 >-------------------------
    !-----------------------DIVIDER-LINE--------------------------
    FUNCTION exist_in(string1,string2)
      IMPLICIT NONE
      LOGICAL :: exist_in
      CHARACTER(len=*) :: string1,string2
      !> if string1 existed in string2, return TRUE
      !> else, return FALSE
      INTEGER(I4B) :: len_str1, len_str2
      INTEGER(I4B) :: i_ei
      !> simplified from "i in exist_in"
      !====== ====== ======
      len_str1=len(string1)
      len_str2=len(string2)
      IF( (len_str1 - len_str2) .GT. 0)THEN
         WRITE(6,*) "read_upf: ERROR; read_module; exist_in"
         stop
      ENDIF
      DO i_ei=1, len_str2 - len_str1 +1, 1
         IF(string1 == string2(i_ei:i_ei+len_str1-1))THEN
            exist_in=.TRUE.
            return
         ENDIF
      ENDDO
      exist_in=.FALSE.
      !====== ====== ======
    END FUNCTION exist_in
    FUNCTION exist_ibegin(string1,string2)
      IMPLICIT NONE
      INTEGER(I4B) :: exist_ibegin
      CHARACTER(len=*) :: string1,string2
      !> if string1 existed in string2, return TRUE
      !> else, return FALSE
      INTEGER(I4B) :: len_str1, len_str2
      INTEGER(I4B) :: i_ei
      !> simplified from "i in exist_in"
      !====== ====== ======
      len_str1=len(string1)
      len_str2=len(string2)
      IF( (len_str1 - len_str2) .GT. 0)THEN
         WRITE(6,*) "read_upf: ERROR; read_module; exist_ibegin"
         stop
      ENDIF
      DO i_ei=1, len_str2 - len_str1 +1, 1
         IF(string1 == string2(i_ei:i_ei+len_str1-1))THEN
            exist_ibegin=i_ei
            return
         ENDIF
      ENDDO
      exist_ibegin=0
      !====== ====== ======
    END FUNCTION exist_ibegin
    !-----------------------DIVIDER-LINE--------------------------
    SUBROUTINE scan_head(file_unit,title,start_old)
      IMPLICIT NONE
      INTEGER(I4B) :: file_unit
      !> Unit of the input file
      CHARACTER(len=*) :: title
      !> Label to be matched
      LOGICAL :: start_old,on_off=.FALSE.
      !> Flag: if .false. rewind the file
      CHARACTER(len=100) :: read_string
      !> The whole line  from file
      INTEGER(I4B) :: ios
      !> I/O state
      INTEGER(I4B) :: error_id
      !> the id in error Information, not necessary
      INTEGER(I4B) :: i,counter=0
      !====== ====== ======
      ! print*,"read"//title
      ios=0
      IF (.NOT.start_old) REWIND(file_unit)
      !> obtain the linage of attribute
      DO WHILE (ios==0)
         READ (file_unit, '(A)', iostat = ios, err = 313) read_string
         !  print *,"read in",read_string,ios
         IF (exist_in("<PP_"//title, read_string) )THEN
            counter=0
            on_off=.TRUE.
         ENDIF
         IF (on_off) THEN
            counter=counter+1
            ! print *,"cycle 1",read_string
         ENDIF
         IF (exist_in(">", read_string).AND.on_off)THEN
            on_off=.FALSE.
            ! print*,"GET DATA END"
            exit
         ENDIF
      ENDDO
      ! print*,"attribute counter ->",on_off,counter
      !> go back "<PP_..."
      DO i=1,counter,1
         backspace(file_unit)
      ENDDO
      !> read
      ALLOCATE(attribute_data(1:counter))
      DO WHILE (ios==0)
         READ (file_unit, '(A)', iostat = ios, err = 313) read_string
         IF (exist_in("<PP_"//title, read_string) ) THEN
            backspace(file_unit)
            ! print *,"cycle 2",read_string
            DO i=1,counter,1
               READ (file_unit, '(A)', iostat = ios, err = 313) attribute_data(i)%value
               ! print *,attribute_data(i)%value
            ENDDO
            return
         ENDIF
      ENDDO
313   print*, 'ERROR:read upf, iostate is ', abs(ios)
      !====== ====== ======
    END SUBROUTINE scan_head
    !-----------------------DIVIDER-LINE--------------------------
    SUBROUTINE scan_tail(file_unit,title)
      !>For test! so not just read to cross
      IMPLICIT NONE
      INTEGER(I4B) :: file_unit
      !> Unit of the input file
      CHARACTER(len=*) :: title
      !> Label to be matched
      CHARACTER(len=80) :: read_string
      !> The whole line  from file
      INTEGER(I4B) :: ios
      !> I/O state
      INTEGER(I4B) :: error_id
      !> the id in error Information, not necessary
      !====== ====== ======
      DEALLOCATE(attribute_data)
      ios=0
      DO WHILE (ios==0)
         READ (file_unit, '(A)', iostat = ios, err = 313) read_string
         IF (exist_in("</PP_"//title//">", read_string) ) RETURN
      ENDDO
      IF(ios/=0)print*,"scan end err"
      return
313   print*, 'ERROR:read upf, iostate is ', title
      !====== ====== ======
    END SUBROUTINE scan_tail
    !-----------------------PARTING-LINE--------------------------
    SUBROUTINE read_upf(Ity,filename,ps)
      USE pspot_module , ONLY : pspot
      USE parameters , ONLY : LRadRho,Lpbc
      USE Math,        ONLY : direct_productlm,fourier_1d,interp,dfdr
      USE MathSplines, ONLY : spline_cubic_set
      !> UPF Format
      !>   PP_INFO
      !>   PP_HEADER
      !>   PP_MESH  :  {PP_R|PP_RAB}
      !>   PP_NLCC (optional)
      !>   PP_LOCAL
      !>   PP_NONLOCAL :  {PP_BETA|PP_DIJ}
      !>   PP_SEMILOCAL (optional, only for norm-conserving)
      !>   PP_PSWFC (optional)
      !>   PP_FULL_WFC (only for PAW)
      !>   PP_RHOATOM
      !>   PP_PAW (only for PAW)
      IMPLICIT NONE

      INTEGER(I4B),INTENT(IN)   :: Ity
      CHARACTER(30),INTENT(IN) :: filename
      !>
      TYPE(pspot)         :: ps
      INTEGER :: is, ios, unit_UPF
      CHARACTER(len=256) :: line
      INTEGER(I4B),ALLOCATABLE :: temp_l(:)
      REAL(DP),ALLOCATABLE     :: temp_beta_r(:,:),tempD0(:,:)
      INTEGER(I4B)             :: temp_nproj,i,j,temp_a=0,m,Id
      INTEGER(I4B)             :: k
      REAL(DP)           :: temp_rcut,dq
      REAL(DP),ALLOCATABLE     :: temp_denr(:)
      !> for array periodic boundary condition
      REAL(DP),ALLOCATABLE     :: temp_rab(:),temp_g(:)&
           &,temp_beta_r_interp(:),temp_r(:)&
           &,temp_beta_r_uni(:,:)
      REAL(DP),ALLOCATABLE :: knots(:)
      !>===========================
      PRINT*,'[Reading PseudoPP FILE]',filename
      unit_UPF=Ity+1600
      OPEN(unit=unit_UPF,file=filename,status='old',form='formatted',iostat=ios)
      WRITE ( *, * ) " Reading pseudopotential file in UPF format..."

      !> Search for Header
      if(allocated(attribute_data))DEALLOCATE(attribute_data)
      CALL scan_head (unit_UPF, "HEADER", .true.)
      CALL read_pseudo_header (ps%Zion,ps%numps,temp_nproj)
      ! print *, "pa%numps",ps%numps
      !CALL scan_tail (unit_UPF, "HEADER")
      DEALLOCATE(attribute_data)

      !> Search for mesh information
      !CALL scan_head (iunps, "MESH", .true.)
      CALL scan_head (unit_UPF, "R ", .true.)
      ALLOCATE(ps%r_real(ps%numps))
      READ(unit_UPF,*) ps%r_real
      !> the max r
      ps%rmax=ps%r_real(ps%numps)
      !CALL read_pseudo_mesh (ps%r_real)
      CALL scan_tail (unit_UPF, "MESH")

      !> Search for Local potential
      CALL scan_head (unit_UPF, "LOCAL", .true.)
      ALLOCATE(ps%V_loc(ps%numps))
      READ(unit_UPF,*) ps%V_loc

      !> Ry to hartree
      ps%V_loc=ps%V_loc*0.5d0
      ! ps%V_loc=ps%V_loc

      !CALL read_pseudo_local( ps%V_loc )
      CALL scan_tail (unit_UPF, "LOCAL")

      !-------->Search for Nonlocal potential
      CALL scan_head (unit_UPF, "NONLOCAL", .true.)
      DEALLOCATE(attribute_data)
      ALLOCATE(temp_beta_r(ps%numps,temp_nproj))
      ALLOCATE(tempD0(temp_nproj,temp_nproj))
      ALLOCATE(temp_l(temp_nproj))
      !READ(unit_UPF,*) ps%beta_r
      CALL read_pseudo_nonlocal(unit_UPF,temp_nproj &
           & ,temp_beta_r,tempD0,temp_rcut,temp_l)
      !CALL scan_tail (unit_UPF, "NONLOCAL")
      !DEALLOCATE(attribute_data)

      !> nolocal cutoff radius
      if(ps%rcut<temp_rcut)ps%rcut=temp_rcut

      !>reset nonlocal pseudopotential
      temp_a=0
      DO j=1,temp_nproj,1
         temp_a=temp_a+2*temp_l(j)+1
      ENDDO

      ALLOCATE(ps%proj_l(temp_a))
      ALLOCATE(ps%proj_m(temp_a))
      ALLOCATE(ps%beta_r(ps%numps,temp_a))
      ALLOCATE(ps%D0(temp_a,temp_a))
      ! print *,'temp_a',temp_a
      ! print *,'shape(ps%beta_r)',shape(ps%beta_r)
      ps%beta_r=0.d0
      ps%D0=0.d0
      ! print*,'beta_r_test',ps%beta_r(1,5)

      ps%nproj = 0
      DO j=1,temp_nproj
         ! print*,'n',j
         ! print*,'temp_l',temp_l(j)
         ps%nproj = ps%nproj+2*temp_l(j)+1 ! 2l+1
         ps%proj_l(ps%nproj-2*temp_l(j):ps%nproj) = temp_l(j)
         DO k=ps%nproj-2*temp_l(j),ps%nproj,1
            ! print*,'l',k,j
            !> Ry to hartree
            if(ps%r_real(1)==0)then
               ps%beta_r(2:,k) = temp_beta_r(2:,j)/ps%r_real(2:)*0.5d0
               ! ps%beta_r(1,k) = temp_beta_r(1,j)/ps%r_real(1)*0.5d0
               ps%beta_r(1,k) = ps%beta_r(2,k)
            else
               ps%beta_r(:,k) = temp_beta_r(:,j)/ps%r_real(:)*0.5d0
            endif
         ENDDO
         ps%proj_m(ps%nproj-2*temp_l(j):ps%nproj) = &
              (/(m,m=-temp_l(j),temp_l(j))/)
      ENDDO

      !> Ry to hartree
      tempD0=tempD0*2.d0
      ! tempD0=1/(tempD0)
      !> assign D0
      CALL direct_productlm(temp_nproj,ps%nproj,temp_l &
           & ,ps%proj_l,tempD0,ps%D0)
      ! write(6,*)'shape ps%D0',shape(ps%D0)
      ! do i=1,size(ps%D0,2),1
      !    write(6,'(1000(F6.2,2X))')ps%D0(:,i)
      ! enddo

     !> Search for atomic charge
     CALL scan_head (unit_UPF, "RHOATOM", .true.)
     ALLOCATE(temp_denr(ps%numps))
     READ(unit_UPF,*) temp_denr
     CALL scan_tail (unit_UPF, "RHOATOM")
     ALLOCATE(ps%denr(ps%numps))
     LRadRho=.true.

     IF(LRadRho)THEN
        IF(ps%r_real(1)==0)then
           ps%denr(2:)= temp_denr(2:)/ps%r_real(2:)**2
           !ps%beta_r(1,k) = temp_beta_r(1,j)/ps%rad(1)*0.5d0
           ps%denr(1) = ps%denr(2)
        ELSE
           ps%denr(:) = temp_denr(:)/ps%r_real(:)**2
        ENDIF
        ps%denr(:)=ps%denr(:)/(4*pi)
     ENDIF

     if(sum(ps%denr)<=0.1d0)then
        LRadRho=.false.
        write(6,*)'read_upf: error, total RHOATOM small than 0.1,'
        write(6,*)'read_upf: start with uniform density'
     endif
     ps%numps_den=ps%numps

     !> reset for periodic boundary condition
     IF(Lpbc)THEN
        !> periodic boundary condition need varables:
        !> ps% Zion rcut
        !> numps nproj denr(numps)
        !> qmax rmax ,parameters for uniform distance
        !> Vlocq(numps)

        !> read rab for fourier transform
        CALL scan_head(unit_UPF,"RAB", .false.)
        ALLOCATE(temp_rab(ps%numps))
        read(unit_UPF,*)temp_rab
        ! do i =1,1258,1
        !    temp_rab(i)=(ps%r_real(i)+0.413125362778D-3)*0.01
        ! enddo

        !> set qmax and Vlocq
        ps%qmax=100.d0*bohr2ang
        ps%qnumps=6001
        dq=ps%qmax/(ps%qnumps-1)

        allocate(ps%Vlocq(ps%qnumps),temp_g(ps%qnumps))
        temp_g=(/(i*dq,i=0,ps%qnumps-1,1)/)
        CALL fourier_1d(ps%numps,ps%r_real,temp_rab,ps%V_loc &
             &,0,ps%qnumps,temp_g,ps%Vlocq,ps%Zion)
        ! open(1212,file='vloc_upf')
        ! write(1212,*)ps%Vlocq
        ! close(1212)

        !> reset the rmax
        ! ps%rmax=20

        !>> reset beta_r in a uniform grid
        ! deallocate(ps%beta_r)
        ! allocate(temp_beta_r_interp(ps%numps))
        ! allocate(ps%beta_r(4001,ps%nproj))
        ! allocate(temp_r(4001))
        ! temp_r=(/(i*0.005d0,i=0,4000)/)
        ! temp_a=0
        ! DO j=1,temp_nproj,1
        !    temp_a=temp_a+2*temp_l(j)+1
        !    !> vr to v and Ry to Hartree
        !    if(ps%r_real(1)==0)then
        !        temp_beta_r_interp(2:)=&
        !             & temp_beta_r(2:,j)!/ps%r_real(2:)!*0.5d0
        !        temp_beta_r_interp(1)=&
        !             & temp_beta_r(1,j)!/ps%r_real(1)!*0.5d0
        !    else
        !       temp_beta_r_interp(:) =&
        !            & temp_beta_r(:,j)!/ps%r_real(:)*0.5d0
        !    endif

        !    !> assign new beta_r
        !    do k=temp_a-2*temp_l(j),temp_a,1
        !    do i=1,4001,1
        !       ps%beta_r(i,k)= interp(ps%numps &
        !            &,temp_beta_r_interp &
        !            &,ps%r_real,(i-1)*0.005d0)
        !    enddo
        !    enddo
        ! ENDDO
        ! open(1111,file='nonlocal')
        ! write(1111,*)ps%beta_r
        ! close(1111)

        !> reset denr
        ! deallocate(ps%denr)
        ! allocate(ps%denr(4001))
        ! do i=1,4001,1
        !    ps%denr(i)= interp(ps%numps &
        !         &,temp_denr,ps%r_real,(i-1)*0.005d0)
        ! enddo

        !> release the array for confined system
        deallocate(ps%V_loc)
        ! deallocate(ps%r_real)
        allocate(ps%qmesh(ps%qnumps))
        ps%qmesh=temp_g

        !================parpare to interplote===============
        ALLOCATE(knots(ps%qnumps))
        DO i=1,ps%qnumps
           knots(i)=REAL(i,DP)
        ENDDO
        ps%qspacing=ps%qmax/real(ps%qnumps-1,DP)
        !> local part
        ALLOCATE(ps%VlocqS(ps%qnumps))
        ALLOCATE(ps%ddVl_dq2(ps%qnumps))
        !> subtract 1/r in q-space
        ps%VlocqS(1)=ps%Vlocq(1)
        ps%VlocqS(2:)=ps%Vlocq(2:)+ 4.d0*pi*ps%Zion/(ps%qspacing*(knots(2:)-1))**2
        !> derivetive 2 order
        CALL spline_cubic_set ( ps%qnumps, knots, ps%VlocqS, &
             &    1, 0._DP, 1, 0._DP, ps%ddVl_dq2)

        !> nonlocal part
        ps%rmax=20
        allocate(temp_beta_r_interp(ps%numps))
        allocate(temp_beta_r_uni(ps%qnumps,ps%nproj))
        allocate(temp_r(ps%qnumps))
        temp_r=(/(i*ps%qspacing,i=0,ps%qnumps)/)
        temp_a=0
        DO j=1,temp_nproj,1
           temp_a=temp_a+2*temp_l(j)+1
           !> vr to v and Ry to Hartree
           if(ps%r_real(1)==0)then
               temp_beta_r_interp(2:)=&
                    & temp_beta_r(2:,j)!/ps%r_real(2:)!*0.5d0
               temp_beta_r_interp(1)=&
                    & temp_beta_r(1,j)!/ps%r_real(1)!*0.5d0
           else
              temp_beta_r_interp(:) =&
                   & temp_beta_r(:,j)!/ps%r_real(:)*0.5d0
           endif

           !> assign new beta_r
           do k=temp_a-2*temp_l(j),temp_a,1
           do i=1,ps%qnumps,1
              temp_beta_r_uni(i,k)= interp(ps%numps &
                   &,temp_beta_r_interp &
                   &,ps%r_real,(i-1)*ps%qspacing)
           enddo
           enddo
        ENDDO
        IF(ps%nproj>0)THEN
           ps%rspacing=ps%rmax/REAL((ps%qnumps-1),DP)
           !
           ! ALLOCATE(ps%ddbeta_dr2(ps%qnumps,ps%nproj))
           ALLOCATE(ps%dbeta_dr(ps%qnumps,ps%nproj))
           DO i=1,ps%nproj
              CALL dfdr(ps%qnumps,ps%rspacing,temp_beta_r_uni(:,i),ps%dbeta_dr(:,i))
           ENDDO
        ENDIF

        !> charge density
        ! IF(LRadRho)THEN
        !    ALLOCATE(ps%ddden_dr2(ps%numps))
        !    CALL spline_cubic_set ( ps%numps, knots, ps%denr, &
        !         &    1, 0._DP, 1, 0._DP, ps%ddden_dr2)
        ! ENDIF
      !========================End========================
     ENDIF

     !>the END
     ! WRITE ( *, * ) "... read UPF file done"
     RETURN
     CLOSE (unit=unit_UPF)
     20 WRITE(6,*)"ERR in read UPF"
   END SUBROUTINE read_upf
   SUBROUTINE read_pseudo_header(Zion,mesh_size,nproj)
     IMPLICIT NONE
     LOGICAL :: find_flag
     REAL(DP) :: Zion
     INTEGER(I4B) :: mesh_size  !>numps
     !INTEGER(I4B) :: lmax
     !INTEGER(I4B) :: l_local
     INTEGER(I4B) :: nproj
     INTEGER(I4B) :: i

     !> STRUCTURE
     !> <PP\_HEADER attr1="value1" ... attrN="valueN"> ... </PP\_HEADER>
     !>    attr             value
     !>    generated       "Generation code"
     !>    author          "Author"
     !>    date            "Generation date"
     !>    comment         "Brief description"
     !>    element         "Chemical Symbol"
     !>    pseudo_type     "NC | SL | 1/r | US | PAW"
     !>    relativistic    "scalar | full | nonrelativistic"
     !>    is_ultrasoft    .F. | .T.
     !>    is_paw          .F. | .T.
     !>    is_coulomb      .F. | .T.
     !>    has_so          .F. | .T.
     !>    has_wfc         .F. | .T.
     !>    has_gipaw       .F. | .T.
     !>    paw\_as\_gipaw    .F. | .T.
     !>    core_correction .F. | .T.
     !>    functional      "dft"
     !>    z_valence        Zval
     !>    total_psenergy   etotps
     !>    wfc_cutoff       ecutwfc
     !>    rho_cutoff       ecutrho
     !>    l_max            lmax
     !>    l\_max\_rho        lmax_rho
     !>    l_local          lloc
     !>    mesh_size        mesh
     !>    number\_of\_wfc    nwfc
     !>    number\_of\_proj   nbeta
     !> </PP_HEADER>
     !>===========================
     DO i=1,size(attribute_data),1
       call get_value(attribute_data(i)%value,'z_valence=',Zion,find_flag)
       IF(find_flag)exit
     ENDDO

     DO i=1,size(attribute_data),1
       call get_value(attribute_data(i)%value,'mesh_size',mesh_size,find_flag)
       IF(find_flag)exit
     ENDDO
     !>==================================
     !>it's so dozzy that reading information make me
      !DO i=1,size(attribute_data),1
      !   call get_value(attribute_data(i)%value,'l_max',l_max,find_flag)
      !   IF(find_flag)exit
      !ENDDO
      !
      !DO i=1,size(attribute_data),1
      !   call get_value(attribute_data(i)%value,'l_local',l_local,find_flag)
      !   IF(find_flag)exit
      !ENDDO

     DO i=1,size(attribute_data),1
       call get_value(attribute_data(i)%value,'number_of_proj',nproj,find_flag)
       IF(find_flag)exit
     ENDDO
     !>==================================
   END SUBROUTINE read_pseudo_header
   SUBROUTINE read_pseudo_nonlocal(unit_UPF,nl,beta_r,D0,rcut,proj_l)
     IMPLICIT NONE
     LOGICAL :: find_flag
     INTEGER(I4B) :: unit_UPF,nl
     REAL(DP)     :: beta_r(:,:)
     INTEGER(I4B) :: proj_l(:)
     INTEGER(I4B) :: i,j
     REAL(DP)     :: D0(:,:)
     REAL(DP)     :: rcut
     CHARACTER(len=2) :: id_temp
     CHARACTER(len=20) :: title_in
     !>===========================
     DO i=1,nl,1
       write(title_in,'("BETA."I1)')i
       CALL scan_head (unit_UPF,trim(adjustl(title_in)),.FALSE.)
       READ(unit_UPF,*)beta_r(:,i)
       !> get angular_momentum
       DO j=1,size(attribute_data),1
         call get_value(attribute_data(j)%value,'angular_momentum',proj_l(i),find_flag)
         IF(find_flag)exit
       ENDDO
       !> get rcut
       DO j=1,size(attribute_data),1
         call get_value(attribute_data(j)%value,'cutoff_radius=',rcut,find_flag)
         IF(find_flag)exit
       ENDDO
       DEALLOCATE(attribute_data)
     ENDDO

     !DEALLOCATE(attribute_data)
     CALL scan_head(unit_UPF, "DIJ",.true.)
     READ(unit_UPF,*) D0
     CALL scan_tail(unit_UPF, "DIJ")
     !>==================================
   END SUBROUTINE read_pseudo_nonlocal
   !>=================================
   SUBROUTINE get_value_int(char_in,char_find,variable,find_flag)
     IMPLICIT NONE
     INTEGER(I4B)  :: i,j,ia,ib,foo
     CHARACTER(len=120),INTENT(IN) :: char_in
     CHARACTER(len=*),INTENT(IN) :: char_find
     INTEGER(I4B)      :: variable
     LOGICAL           :: find_flag
     !>read "char_find" from "char_in" and assign to variable
     IF(exist_in(char_find,char_in))THEN
       foo=exist_ibegin(char_find,char_in)
       find_value1: DO j=foo,120,1
         IF(char_in(j:j)=='"')THEN
           ia=j
           exit find_value1
         ENDIF
       ENDDO find_value1
       find_value2: DO j=ia+1,120,1
         IF(char_in(j:j)=='"')THEN
           ib=j
           exit find_value2
         ENDIF
       ENDDO find_value2
       READ(char_in(ia+1:ib-1),*)variable
       ! print *, char_find, variable, char_in(ia+1:ib-1)
       find_flag=.TRUE.
     ELSE
       find_flag=.FALSE.
     ENDIF
   END SUBROUTINE get_value_int
   SUBROUTINE get_value_real(char_in,char_find,variable,find_flag)
     IMPLICIT NONE
     INTEGER(I4B)  :: i,j,ia,ib,foo
     CHARACTER(len=120),INTENT(IN) :: char_in
     CHARACTER(len=*),INTENT(IN) :: char_find
     REAL(DP)      :: variable
     LOGICAL        :: find_flag
     !>read "char_find" from "char_in" and assign to variable
     IF(exist_in(char_find,char_in))THEN
       foo=exist_ibegin(char_find,char_in)
       find_value1: DO j=foo,120,1
         IF(char_in(j:j)=='"')THEN
           ia=j
           exit find_value1
         ENDIF
       ENDDO find_value1
       find_value2: DO j=ia+1,120,1
         IF(char_in(j:j)=='"')THEN
           ib=j
           exit find_value2
         ENDIF
       ENDDO find_value2
       READ(char_in(ia+1:ib-1),*)variable
       ! print *, char_find, variable, char_in(ia+1:ib-1)
       find_flag=.TRUE.
     ELSE
       find_flag=.FALSE.
     ENDIF
   END SUBROUTINE get_value_real
   !>==========================={UPF_END}============================
   !> read yaehmop output
   SUBROUTINE get_MO_coefficient()
     USE struct_module
     USE parameters,     ONLY:nstates,MO_file
     IMPLICIT NONE
     LOGICAL  :: logical_spec_line, read_flag
     CHARACTER(len=100) :: str
     CHARACTER(len=3)   :: atom_sign(3)
     CHARACTER(len=6)   :: atom_orbital(3)
     INTEGER(I4B)       :: atom_id(3), l(3),m(3)
     INTEGER(I4B)       :: MO_id,MO_id_test
     INTEGER(I4B)       :: parentheses_id(2),i_o,i,j  !> read occupation, counters

     read_flag=.TRUE.
     !read coefficient
     open(555,file=MO_file)
     read(555,'(A)')str
     !>find line which satisfy re(^;   > REAL:)
     logical_spec_line = .not.(exist_in(";	***> REAL:",str))
     DO WHILE(logical_spec_line)
       read(555,'(A)',err=256)str
       logical_spec_line = .not.(exist_in(";	***> REAL:",str))
     ENDDO

     !>obtain the MO coefficient in AOs of all atoms
     !>COEFF (n_atom, n_l, n_m, n_orbital)
     IF(.NOT.ALLOCATED(struct%coeff))&
      & ALLOCATE(struct%coeff(natom,0:2,-2:2,Nstates))
      struct%coeff=99999.d0
     READ(555,'(A)',err=257)str
     DO WHILE(read_flag)
       CALL parse_headline(str,atom_sign,atom_id,atom_orbital)
       CALL sign2lm(atom_orbital,l,m)
       !>read the sub segment
       IF(atom_id(2)==0)THEN
         do MO_id = 1,Nstates
           read(555,*,err=257) MO_id_test, &
                            &  struct%coeff(atom_id(1),l(1),m(1),MO_id)
         enddo
       ELSEIF(atom_id(3)==0)THEN
         do MO_id = 1,Nstates
           read(555,*,err=257) MO_id_test, &
                            &  struct%coeff(atom_id(1),l(1),m(1),MO_id),&
                            &  struct%coeff(atom_id(2),l(2),m(2),MO_id)
         enddo
       ELSE
         do MO_id = 1,Nstates
           read(555,*,err=257) MO_id_test, &
                            &  struct%coeff(atom_id(1),l(1),m(1),MO_id),&
                            &  struct%coeff(atom_id(2),l(2),m(2),MO_id),&
                            &  struct%coeff(atom_id(3),l(3),m(3),MO_id)
         enddo
       ENDIF
       !>resume the last id
       MO_id=MO_id-1
       !>judge out
       IF(atom_id(3)==0) read_flag=.FALSE.
       IF(MO_id/=MO_id_test) THEN
         stop "MO_id of MO_file read conflict"
       ENDIF

       !>locate the HEADLINE
       read(555,'(A)') str
       str=adjustl(str)
       DO WHILE(iachar('A')>=iachar(str(1:1)) &
            &  .or.iachar(str(1:1))>=iachar('Z'))
         read(555,'(A)') str
         str=adjustl(str)
         IF(str(1:1)==';'.or.str(1:1)=='#')THEN
           read_flag=.FALSE.
           exit
         ENDIF
       ENDDO
     ENDDO

     !> read the occupation info
     IF(.NOT.ALLOCATED(struct%occupy))&
      & ALLOCATE(struct%occupy(Nstates))
     logical_spec_line = .not.(exist_in("**** Energies (in eV)  and Occupa",str))
     DO WHILE(logical_spec_line)
        read(555,'(A)',err=256)str
       logical_spec_line = .not.(exist_in("**** Energies (in eV)  and Occupa",str))
     ENDDO
    struct%Noccupy=0
    struct%occupy=0.d0
     DO i_o=1,Nstates,1
        read(555,'(A)',err=256)str
       !> split string by []
       !> character id
       i=1
       !>counter of parentheses
       j=1
       do while(j.lt.3)
         do while(.true.)
           if(index(str(i:),'[')==1.or.index(str(i:),']')==1)exit
           if(index(str(i:),'[')==0.and.index(str(i:),']')==0)exit
           i=i+1
         enddo
         parentheses_id(j)=i
         i=i+1
         j=j+1
       enddo
       !> read occupy(i_o)
       READ(str(parentheses_id(1)+1:parentheses_id(2)-1),*,err=258)struct%occupy(i_o)
       struct%Noccupy=struct%Noccupy+1
       if(struct%occupy(i_o).lt.0.00001d0) exit
     ENDDO
     close(555)
     RETURN
     256 WRITE(6,*) "MO.info: can't find the line start with ';  &
                     & ***> REAL'"
     257 WRITE(6,*) "READ MO COEFFICIENT ERROR =>> Read_module ",MO_id,MO_id_test
     258 WRITE(6,*) "READ MO OCCUPATION ERROR =>> Read_module ","states_n",i_o
   END SUBROUTINE get_MO_coefficient
   !===================================================================
   SUBROUTINE parse_headline(str, atom_sign, atom_id, atom_orbital)
     IMPLICIT NONE
     CHARACTER(len=100) :: str
     CHARACTER(len=3)   :: atom_sign(3)
     INTEGER(I4B)       :: atom_id(3)
     CHARACTER(len=6)   :: atom_orbital(3)
     INTEGER(I4B)       :: parentheses_id(6),i,j
     !>split string by parentheses
     !>character id
     i=1
     !>counter of parentheses
     j=1
     do while(j.lt.7)
       do while(.true.)
         if(index(str(i:),'(')==1.or.index(str(i:),')')==1)exit
         if(index(str(i:),'(')==0.and.index(str(i:),')')==0)exit
         i=i+1
       enddo
       parentheses_id(j)=i
       i=i+1
       j=j+1
     enddo
     atom_id=0
     !>READ HEADLINE
     IF((parentheses_id(4)-parentheses_id(3))<2)THEN
       read(str(1:parentheses_id(1)-1),*)atom_sign(1)
       read(str(parentheses_id(1)+1:parentheses_id(2)-1),*)atom_id(1)
       read(str(parentheses_id(2)+1:),*) atom_orbital(1)
       atom_id(2:3)=0
       atom_sign(2:3)="NULL"
       atom_orbital(2:3)="NULL"
     ELSEIF((parentheses_id(6)-parentheses_id(5))<2)THEN
       read(str(1:parentheses_id(1)-1),*)atom_sign(1)
       read(str(parentheses_id(1)+1:parentheses_id(2)-1),*)atom_id(1)
       read(str(parentheses_id(2)+1:parentheses_id(3)-1),*) &
          & atom_orbital(1),atom_sign(2)
       read(str(parentheses_id(3)+1:parentheses_id(4)-1),*)atom_id(2)
       read(str(parentheses_id(4)+1:),*) atom_orbital(2)
       atom_id(3)=0
       atom_sign(3)="NULL"
       atom_orbital(3)="NULL"
     ELSE
       read(str(1:parentheses_id(1)-1),*)atom_sign(1)
       read(str(parentheses_id(1)+1:parentheses_id(2)-1),*)atom_id(1)
       read(str(parentheses_id(2)+1:parentheses_id(3)-1),*) &
           & atom_orbital(1),atom_sign(2)
       read(str(parentheses_id(3)+1:parentheses_id(4)-1),*)atom_id(2)
       read(str(parentheses_id(4)+1:parentheses_id(5)-1),*) &
           & atom_orbital(2),atom_sign(3)
       read(str(parentheses_id(5)+1:parentheses_id(6)-1),*)atom_id(3)
       read(str(parentheses_id(6)+1:),*)atom_orbital(3)
     ENDIF
   END SUBROUTINE parse_headline
   !===================================================================
   SUBROUTINE sign2lm(atom_orbital,l,m)
     IMPLICIT NONE
     CHARACTER(len=6),INTENT(IN) :: atom_orbital(3)
     INTEGER(I4B)                :: l(3),m(3),i
     !> identify the s p d f ...
     DO i=1,3,1
       IF(index(atom_orbital(i),'s')/=0)THEN
         l(i)=0
         m(i)=0
       ELSEIF(index(atom_orbital(i),'p')/=0)THEN
         l(i)=1
         IF(index(atom_orbital(i),'x')/=0)m(i)=-1
         IF(index(atom_orbital(i),'y')/=0)m(i)=1
         IF(index(atom_orbital(i),'z')/=0)m(i)=0
       ELSEIF(index(atom_orbital(i),'d')/=0)THEN
         l(i)=2
         IF(exist_in('dxy',atom_orbital(i)))   m(i)=-2
         IF(exist_in('dxz',atom_orbital(i)))   m(i)=-1
         IF(exist_in('dz2',atom_orbital(i)))   m(i)=0
         IF(exist_in('dyz',atom_orbital(i)))   m(i)=1
         IF(exist_in('dx2y2',atom_orbital(i))) m(i)=2
       ENDIF
     ENDDO
   END SUBROUTINE sign2lm
   !===================================================================
   SUBROUTINE destroy_MOinit()
     USE struct_module, ONLY:struct
     IMPLICIT NONE
     IF(ALLOCATED(struct%coeff))THEN
       DEALLOCATE(struct%coeff)
     ENDIF
   END SUBROUTINE destroy_MOinit
   !===================================================================
   SUBROUTINE out_CONCAR(filename)
     USE struct_module , ONLY: lat_mat, struct, natom
     IMPLICIT NONE
     CHARACTER(len=*) :: filename
     INTEGER(I4B) :: i


     open(1234,file=filename)
     write(1234,*)"structure"
     write(1234,*)lat_mat
     write(1234,*)"element"
     write(1234,*)natom
     write(1234,*)"Car"
     do i=1,natom,1
        write(1234,*)struct%poscar(:,i)
     enddo
     write(1234,*)"Dir"
     do i=1,natom,1
        write(1234,*)struct%pos(:,i)
     enddo
     close(1234)
   ENDSUBROUTINE out_CONCAR
   SUBROUTINE read_CHGCAR()
     !> read chgcar for set boundary
     IMPLICIT NONE
     INTEGER(I4B)  :: chg_n1,chg_n2,chg_n3
     REAL(DP),allocatable ::CHGCAR_in(:,:,:)

     open(1111,file='CHGCAR')
     read(1111,*)chg_n1,chg_n2,chg_n3
     allocate(CHGCAR_in(chg_n1,chg_n2,chg_n3))
     read(1111,*)CHGCAR_in
     close(1111)
   ENDSUBROUTINE read_CHGCAR
END MODULE read_module
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

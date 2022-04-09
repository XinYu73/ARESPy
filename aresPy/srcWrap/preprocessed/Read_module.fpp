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
CONTAINS
   !------------------------read input.dat------------------------
   SUBROUTINE read_file(infile)
      USE parameters
      USE math, ONLY:change_case,find_keywords
      USE pspot_module, ONLY : read_pspot
      IMPLICIT NONE
      !-----------------------------------------------------------
      CHARACTER(7),INTENT(IN) :: infile
      !-----------------------------------------------------------
      INTEGER(I4B)    :: l_str
      INTEGER(I4B)    :: ios,id_ex,id_pound,id_key,id_value !id of key words
      INTEGER(I4B)    :: i,j,k
      INTEGER(I4B)    :: vmajor,vminor,vmicro ! XC version
      CHARACTER       :: ch_mark
      CHARACTER(100)  :: instr
      CHARACTER(100)  :: str
      LOGICAL         :: lexist
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      IF(parallel%isroot)THEN

      !-------start input.dat-----------
      INQUIRE(FILE=infile,EXIST=lexist)
      !test the input.dat
      IF(.NOT.lexist)THEN
         WRITE(*,*) '>>>WARNING<<<:input file is not exist.stop!!!'
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
            WRITE(*,*) 'Task name >>> ',system_name
         ENDIF
         !cell name
         IF(str(1:id_key)=='cellfile')THEN
            str=instr
            cellfile_name=TRIM(ADJUSTL(str(id_value:l_str)))
         ENDIF
         !pseudopotential file name
         IF(str(1:id_key)=='ppfile')THEN
            str=instr
            ntype=0
            DO i=id_value,l_str
               IF(str(i:i) == " " .and.  str(i-1:i-1) /= " " )THEN
                  ntype=ntype+1
                  ppfile_name(ntype)=trim(adjustl(str(k:i-1)))
               ENDIF
               IF (str(i-1:i-1) == " " .and. str(i:i) /= " ") then
                  k=i
               ENDIF
            ENDDO
            ntype=ntype+1
            ppfile_name(ntype)=trim(adjustl(str(k:l_str)))
         ENDIF
         !Ecut
         IF(str(1:id_key) == 'ecut')THEN
            READ(str(id_value:l_str),*) Ecut
            IF(Ecut>0.d0)THEN
               Ecut=Ecut*ev2hart
               init_gap=SQRT((pi**2/Ecut)/2)
               PRINT*,'Ecut=',Ecut*hart2eV,'eV'
            ENDIF
         ENDIF 
         !gridsize
         IF(str(1:id_key) == 'gridspacing')THEN
            IF(Ecut<=0.d0)THEN
               READ(str(id_value:l_str),*) init_gap
               init_gap=init_gap*ang2bohr
               Ecut=(pi**2/init_gap**2)/2
               PRINT*,'Ecut=',Ecut*hart2eV,'eV'
            ENDIF
         ENDIF 
!gridmesh
         IF(str(1:id_key) == 'gridmesh')THEN
!IF(Ecut<=0.d0)THEN
               READ(str(id_value:l_str),*) gridn
!ENDIF
         ENDIF 
!k-grids
         IF(str(1:id_key)=='kspacing')THEN
!ang^-1 to bohr^-1
            READ(str(id_value:l_str),*) kspacing
            kspacing=kspacing*bohr2ang !/(2._DP*pi)
         ENDIF
!kgrid mesh
         IF(str(1:id_key)=='kgrids')THEN
            READ(str(id_value:l_str),*) kgrid
         ENDIF
!finite difference order
         IF( str(1:id_key) == 'norder' )THEN
            READ(str(id_value:l_str),*) finite_order
            IF(finite_order<0)THEN
               finite_order = -finite_order
            ENDIF
            finite_order=finite_order/2
            IF (finite_order>1000 .or. finite_order<0) THEN
               WRITE(*,*) "STOP!!:the order of finite difference should in [0-40]"
               STOP
            ENDIF
            WRITE(*,*) "The order of finite difference is:",finite_order*2
         ENDIF
!diag(H)
         IF(str(1:id_key)=='nadst')THEN
            READ(str(id_value:l_str),*) NADDSTATES
         ENDIF
!>>>spin
         IF(str(1:id_key)=='lspin')THEN
            IF(str(id_value:id_value)=='t')THEN
                NSPIN=2
                WRITE(*,*) "Spin polarized"
            ELSE
                NSPIN=1
                WRITE(*,*) "Spin unpolarized"
            ENDIF
         ENDIF
!---
         IF( str(1:id_key) == 'snlcc' )THEN
            READ(str(id_value:l_str),*) Snlcc
!Checking
            IF(Snlcc<0._DP.OR.Snlcc>1._DP)THEN
               PRINT*,'STOP: Snlcc must be in [0,1]'
               STOP
            ENDIF
         ENDIF
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
                WRITE(*,*) "Chebyshev filter is used,order:",CheM
            ENDIF
         ENDIF
!CheM0
         IF(str(1:id_key) == 'chem0')THEN
            READ(str(id_value:l_str),*) CheM0
            IF(Idiag/=0)THEN
                WRITE(*,*) "first Chebyshev filter order:",CheM0
            ENDIF
         ENDIF
!first by arpack?
         IF(str(1:id_key)=='lfirst')THEN
            IF(str(id_value:id_value)=='t')THEN
                LFIRST=.TRUE.
                IF(Idiag==1)THEN
                   WRITE(*,*) 'Chebyshev:first step by ARPACK'
                ENDIF
            ELSE
                LFIRST=.FALSE.
                IF(Idiag==1)THEN
                   WRITE(*,*) 'Chebyshev:first step diag free'
                ENDIF
            ENDIF
         ENDIF
!Lorthnorm
         IF(str(1:id_key)=='lorthnorm')THEN
            IF(str(id_value:id_value)=='t')THEN
                LRROrthNorm=.TRUE.
                WRITE(*,*) 'RayleighRitz need OrthNorm step'
            ELSE
                LRROrthNorm=.FALSE.
                WRITE(*,*) "RayleighRitz needn't OrthNorm step"
            ENDIF
         ENDIF
!Lrandom
         IF(str(1:id_key)=='lrandom')THEN
            IF(str(id_value:id_value)=='t')THEN
                Lrandom=.TRUE.
                WRITE(*,*) 'Initialize the subspace by full random'
            ELSE
                Lrandom=.FALSE.
                WRITE(*,*) "Initialize the subspace by Pseudo orbitals"
            ENDIF
         ENDIF
!linrho
         IF(str(1:id_key)=='linrho')THEN
            IF(str(id_value:id_value)=='t')THEN
                LINRHO=.TRUE.
                WRITE(*,*) 'initial density by input'
            ELSE
                LINRHO=.FALSE.
                WRITE(*,*) 'initial density by program'
            ENDIF
         ENDIF
!>>>mixer data
!Type for mixing
         IF( str(1:id_key) == 'imixer' )THEN
            READ(str(id_value:l_str),*) IMIXER
            SELECT CASE(IMIXER)
! CASE(1)
!     WRITE(*,*) 'Use Simple Mixing + Anderson Mixing'
            CASE(0)
                WRITE(*,*) 'Use Simple Mixing + (r)Pulay Mixing'
            CASE(2)
                WRITE(*,*) 'Use Simple Mixing + (r)Pulay Mixing + resta'
            CASE default
                WRITE(*,*) 'Use Simple Mixing + (r)Pulay Mixing + kerker'
            ENDSELECT
         ENDIF
!max step
         IF( str(1:id_key) == 'nmiter' )THEN
            READ(str(id_value:l_str),*) NMITER
         ENDIF
!simple mixing step
         IF(str(1:id_key)=='nsmix')THEN
            READ(str(id_value:l_str),*) NSMIX
            IF(NSMIX>NMITER)THEN
               WRITE(*,*) 'STOP!!:the simple step shoule less than the max step'
            ENDIF
         ENDIF
!history step number
         IF(str(1:id_key)=='nhmax')THEN
            READ(str(id_value:l_str),*) NHMAX
            IF(NHMAX>NMITER)THEN
               WRITE(*,*) 'STOP!!:the history step shoule less than the max step'
            ENDIF
         ENDIF
!min
         IF(str(1:id_key)=='nhmin')THEN
            READ(str(id_value:l_str),*) NHMIN
            IF(NHMIN<=0.AND.NHMIN>=ABS(NHMAX))THEN
               WRITE(*,*) 'STOP!!:the mininal history step > 0 and <ABS(NHMAX)'
            ENDIF
         ENDIF
!simple mixing factor
         IF(str(1:id_key)=='malpha')THEN
            READ(str(id_value:l_str),*) MALPHA
            IF(MALPHA>1.d0.OR.MALPHA<0.d0)THEN
               WRITE(*,*) 'STOP!!:the MALPHA should in [0,1] '
            ENDIF
         ENDIF
!Anderson mixing factor
         IF(str(1:id_key)=='mbeta')THEN
            READ(str(id_value:l_str),*) mbeta
            IF(MBETA>1.d0.OR.MBETA<0.d0)THEN
               WRITE(*,*) 'STOP!!:the MBETA should in [0,1] '
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
         IF(str(1:id_key)=='w0ma')THEN
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
         IF(str(1:id_key)=='lband')THEN
            IF(str(id_value:id_value)=='t')THEN
                Lband=.TRUE.
            ELSE
                Lband=.FALSE.
            ENDIF
         ENDIF
!relaxion step
!max step
         IF( str(1:id_key) == 'nssp' )THEN
            READ(str(id_value:l_str),*) Nssp
            Nssp=MAX(Nssp,0)
            WRITE(*,*) 'Max simulate steps is:',Nssp
         ENDIF
!Atomic rcut
         IF(str(1:id_key)=='atomrc')THEN
            READ(str(id_value:l_str),*) AtomRc
            AtomRc=AtomRc*ang2bohr
            WRITE(*,*) '[The Atomic Spherical Cutting Radius is (Angs)]',AtomRc*bohr2ang
         ENDIF

         IF(str(1:id_key)=='dcharge')THEN
            READ(str(id_value:l_str),*) dcharge
            WRITE(*,*) '[The Adding Charge]',dcharge
         ENDIF
!> Parallel dims
         IF(str(1:id_key)=='npara')THEN
            READ(str(id_value:l_str),*) NPARA
            WRITE(*,*) '[Parallel split by user]',NPARA
         ENDIF

!====================== Outputing ====================
         IF(str(1:id_key)=='lcharge')THEN
            IF(str(id_value:id_value)=='t')THEN
                LCHARGE=.TRUE.
            ELSE
                LCHARGE=.FALSE.
            ENDIF
         ENDIF

         IF(str(1:id_key)=='lwave')THEN
            IF(str(id_value:id_value)=='t')THEN
                LWAVE=.TRUE.
            ELSE
                LWAVE=.FALSE.
            ENDIF
         ENDIF

         IF(str(1:id_key)=='lmom')THEN
            IF(str(id_value:id_value)=='t')THEN
                LMOM=.TRUE.
                PRINT*,'[Maximum Overlap Method is Used]'
            ELSE
                LMOM=.FALSE.
            ENDIF
         ENDIF


         IF( str(1:id_key) == 'ixc' )THEN
            READ(str(id_value:l_str),*) ixc
         ENDIF

      ENDDO

      CLOSE(100)
!--------end input.dat------------
      IF(CheM<1)THEN
         CheM=16
      ENDIF
!-------set nstates---------
      IF(nspin==1)THEN
         IF(LMOM)THEN
            nspin=2
            PRINT*,'MOM:Note that reseted the nspin=2'
         ENDIF
      ENDIF


      ENDIF !> read on process 0 finished
!bcast the parameters
!INTEGER
      CALL MPI_BCAST(nspin        ,1,MPI_INTEGER4 ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(finite_order ,1,MPI_INTEGER4 ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(ntype        ,1,MPI_INTEGER4 ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(naddstates   ,1,MPI_INTEGER4 ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(ISTART       ,1,MPI_INTEGER4 ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(Idiag        ,1,MPI_INTEGER4 ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(CheM0        ,1,MPI_INTEGER4 ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(CheM         ,1,MPI_INTEGER4 ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(gridn        ,3,MPI_INTEGER4 ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(kgrid        ,3,MPI_INTEGER4 ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(nssp         ,1,MPI_INTEGER4 ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(Nsmear       ,1,MPI_INTEGER4 ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(IMIXER       ,1,MPI_INTEGER4 ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(NMITER       ,1,MPI_INTEGER4 ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(NSMIX        ,1,MPI_INTEGER4 ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(NHMAX        ,1,MPI_INTEGER4 ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(NHMIN        ,1,MPI_INTEGER4 ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(NMITER       ,1,MPI_INTEGER4 ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(BLOCK_MBNB   ,1,MPI_INTEGER4 ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(nwf0         ,1,MPI_INTEGER4 ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(NPARA        ,2,MPI_INTEGER4 ,parallel%rootid,parallel%comm,mpinfo)
!REAL
      CALL MPI_BCAST(dcharge      ,1,MPI_REAL8    ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(Snlcc        ,1,MPI_REAL8    ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(init_gap     ,1,MPI_REAL8    ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(Ecut         ,1,MPI_REAL8    ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(kspacing     ,1,MPI_REAL8    ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(Wsmear       ,1,MPI_REAL8    ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(MALPHA       ,1,MPI_REAL8    ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(MBETA        ,1,MPI_REAL8    ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(AMIX        ,1,MPI_REAL8    ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(BMIX        ,1,MPI_REAL8    ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(RESTA        ,3,MPI_REAL8    ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(W0AM         ,1,MPI_REAL8    ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(RTOL         ,1,MPI_REAL8    ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(ETOL         ,1,MPI_REAL8    ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(AtomRc       ,1,MPI_REAL8    ,parallel%rootid,parallel%comm,mpinfo)
!LOGICAL
      CALL MPI_BCAST(LFIRST       ,1,MPI_LOGICAL  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(LINRHO       ,1,MPI_LOGICAL  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(LAtomRho     ,1,MPI_LOGICAL  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(LRROrthNorm  ,1,MPI_LOGICAL  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(Lrandom      ,1,MPI_LOGICAL  ,parallel%rootid,parallel%comm,mpinfo)
!
!
      CALL MPI_BCAST(LWAVE        ,1,MPI_LOGICAL  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(LCHARGE      ,1,MPI_LOGICAL  ,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(Lband        ,1,MPI_LOGICAL  ,parallel%rootid,parallel%comm,mpinfo)
!for PBC init
      CALL smpi_init_pbc()

!read POSCAR
      CALL read_poscar(ntype,cellfile_name)
!read pseudopotentials
      CALL read_pspot(ntype,ppfile_name)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   END SUBROUTINE read_file
!---------------------------read poscar------------------------
   SUBROUTINE read_poscar(nty,filename)
      USE parameters  , ONLY :  elements,init_gap
      USE math,  ONLY : dir2car,car2dir,det,change_case &
                  &, atom_mass
      USE struct_module
      IMPLICIT NONE
!IN/OUT
      INTEGER(I4B),INTENT(IN)     :: nty
      CHARACTER(LEN=*),INTENT(IN) :: filename
!

      INTEGER(I4B)  :: i,j,natom_3,nty_3

      INTEGER(I4B)  :: ele_n(nty),ele_id(nty+1)
      INTEGER(I4B)  :: filestatu
      REAL(DP)      :: lat_ratio
      REAL(DP)      :: lat_read(3,3)
      CHARACTER(1)  :: chr
      CHARACTER(30) :: ctmp
      LOGICAL       :: ldir
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      IF(parallel%isroot)THEN


      OPEN(101,FILE=filename,IOSTAT=filestatu)

          IF(filestatu/=0)THEN
             WRITE(*,*) 'STOP:Could not open',filename
             STOP
          ENDIF

          READ(101,*)          !title,no use
          READ(101,*) lat_ratio
!read lattice matrix
          DO i=1,3
             READ(101,*) lat_read(:,i)
          ENDDO
!filter the 1e-7
!turn to atomic unit
          lat_mat(:,:)=lat_read(:,:)*lat_ratio*ang2bohr
!read over poscar
          DO
             READ(101,'(A)') ctmp
             j=INDEX(ctmp,'Dir')
             i=INDEX(ctmp,'dir')
             IF(j/=0.OR.i/=0)THEN
                ldir=.TRUE.
                EXIT
             ENDIF
             j=INDEX(ctmp,'Car')
             i=INDEX(ctmp,'car')
             IF(j/=0.OR.i/=0)THEN
                ldir=.FALSE.
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
          struct%eleid(:) = ele_id(:)
          struct%nati(:) = ele_n(:)
!   ! store name
!   DO i=1,nty
!      CALL change_case(elements(i)(1:1),elements(i)(1:1),1)
!   ENDDO
!   struct%elements(1:nty)=elements(1:nty)(1:3)
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
      CLOSE(101)
      CALL ResetPOS(natom,lat_mat,struct%pos,struct%poscar)
      OPEN(102,FILE='NEWCAR')
         WRITE(102,*) 'Cell'
         WRITE(102,*) '1.0'
         WRITE(102,*) lat_mat(:,:)*bohr2ang
         WRITE(102,*) struct%nati(:)
         WRITE(102,*) 'Direct'
         WRITE(102,*)  struct%pos(:,:)
      CLOSE(102)

      ENDIF !> if parallel%isroot
!parallel bcast
      CALL MPI_BCAST(natom,1,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
!nati
      IF(.NOT.parallel%isroot)THEN
         CALL creat_struct(nty,natom)
      ENDIF
      CALL MPI_BCAST(struct%nati     ,nty  ,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(struct%mass     ,nty  ,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(struct%eleid    ,nty+1,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
      natom_3=3*natom
      CALL MPI_BCAST(struct%pos      ,natom_3,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(struct%poscar   ,natom_3,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
      nty_3=3*nty
      CALL MPI_BCAST(struct%elements ,nty_3,MPI_CHARACTER,parallel%rootid,parallel%comm,mpinfo)
      CALL MPI_BCAST(lat_mat         ,9,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
!at one time
      CALL MPI_Barrier(parallel%comm,mpinfo)

!print*,struct%eleid
!STOP
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE read_poscar
!---------------------------Reset POS--------------------------
   SUBROUTINE ResetPOS(natom,lat,pos,poscar)
      !set z-axis as tallest
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN) :: natom
      REAL(DP),INTENT(INOUT) :: lat(3,3)
      REAL(DP),INTENT(INOUT) :: pos(3,natom),poscar(3,natom)
      !LOCAL
      REAL(DP) :: alat(3),apos(natom),aposcar(natom)
      REAL(DP) :: x2,y2,z2
      LOGICAL :: lflag
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      x2=SUM(lat(:,1))**2
      y2=SUM(lat(:,2))**2
      z2=SUM(lat(:,3))**2
!z-
      lflag= (z2>=y2).AND.(z2>=x2)
      IF(lflag) RETURN
!y->z, z->x, x->y
      lflag= (y2>=x2).AND.(y2>=z2)
      IF(lflag)THEN
!z
         alat(:)=lat(:,3)
         apos(:)=pos(3,:)
         aposcar(:)=poscar(3,:)
!y->z
         lat(:,3)=lat(:,2)
         pos(3,:)=pos(2,:)
         poscar(3,:)=poscar(2,:)
!x->y
         lat(:,2)=lat(:,1)
         pos(2,:)=pos(1,:)
         poscar(2,:)=poscar(1,:)
!z->x
         lat(:,1)=alat(:)
         pos(1,:)=apos(:)
         poscar(1,:)=aposcar(:)

         RETURN
      ENDIF 
!x->z, z->y, y->x
      lflag= (x2>=y2).AND.(x2>=z2)
      IF(lflag)THEN
!z
         alat(:)=lat(:,3)
         apos(:)=pos(3,:)
         aposcar(:)=poscar(3,:)
!x->z
         lat(:,3)=lat(:,1)
         pos(3,:)=pos(1,:)
         poscar(3,:)=poscar(1,:)
!y->x
         lat(:,1)=lat(:,2)
         pos(1,:)=pos(2,:)
         poscar(1,:)=poscar(2,:)
!z->y
         lat(:,2)=alat(:)
         pos(2,:)=apos(:)
         poscar(2,:)=aposcar(:)
         RETURN
      ENDIF

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE ResetPOS
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
END MODULE read_module   

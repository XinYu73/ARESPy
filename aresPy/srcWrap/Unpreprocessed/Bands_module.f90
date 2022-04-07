!############################################!
!*For    : draw the bands                    !
!*Author : Qiang Xu                          !
!*Date   : 2017/11/13                        !
!############################################!
MODULE band_structure
   USE grid_module , ONLY : n,n1,n2,n3
   USE constants
   IMPLICIT NONE
   !TYPE kpPATH_type
   !     CHARACTER(70)        :: SNAMK !name
   !     INTEGER(I4B) :: NKPTS ,  &!Number of k-points tot path
   !              &      NLine ,  &!number of lines
   !              &      NKPT  ,  &!number of k-points per line
   !              &      NEV       !number of states
   !     REAL(DP),ALLOCATABLE :: VKPT(:,:) !vector of vec(k)
   !     REAL(DP),ALLOCATABLE :: EIGVAL(:,:)
   !     REAL(DP),ALLOCATABLE :: WTK(:)
   !ENDTYPE kpPATH_type

   !
   !INTEGER(I4B),PARAMETER :: NKDIMD=1000

   !TYPE(kpPATH_type)  :: KPTB !k-points for band structure
CONTAINS
   !-------------------band_Structural-----------------------
   SUBROUTINE Init_BandStruct(numk,kvec)
      USE grid_module , ONLY : KPT,nk,destroy_KPT
      USE struct_module , ONLY : recip_lat
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN) :: numk
      REAL(DP),INTENT(IN) :: kvec(:,:)
      !LOCAL
      INTEGER(I4B) :: Ik
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      nk=numk
      !destroy arrays
      CALL destroy_KPT
      ALLOCATE(KPT%vec(3,nk))
      ALLOCATE(KPT%vcar(3,nk))
      ALLOCATE(KPT%wk(nk))
      !kvec
      KPT%vec(:,:)=kvec(:,:)
      DO Ik=1,nk
         KPT%vcar(:,Ik)=MATMUL(recip_lat,KPT%vec(:,Ik))
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE Init_BandStruct
   !--------------------set_bands_data-----------------------
   SUBROUTINE Band_Begin(numk,kvec)
      USE parameters, ONLY : finite_order
      USE Grid_module , ONLY : Build_rgrid, gap , nk&
#ifdef MPI
                          &,global_n3,global_n2,global_n1&
#endif
                          &,FillQTable,FillRTable
      USE struct_module , ONLY : lat_mat,recip_lat,volume &
                           &,charge_ave , ncharge ,struct, &
                           & lat_para,reclat_para
      USE math , ONLY : inv_33,det,dir2car,lat2matrix
      USE potential_module , ONLY : vlpp
#ifdef MPI
      USE nlpot_module , ONLY : initialize_nlpot_band
#else
      USE nlpot_module , ONLY : initialize_nlpot
#endif
      USE FOURIER, ONLY : FFT,PlanFFT,CleanFFT
      USE finite_module , ONLY : init_finite
#ifdef MPI
      USE smpi_math_module, ONLY: parallel
#endif
      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B),INTENT(IN) :: numk
      REAL(DP),INTENT(IN) :: kvec(:,:) !dir vectors
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     !pos to poscar
     CALL dir2car(struct%pos,struct%poscar,lat_mat)
     !set cell
     recip_lat(:,:)=2.d0*pi*(inv_33(TRANSPOSE(lat_mat)))
     volume = ABS(det(lat_mat))
     !cell_para
     !PRINT*,'Building real space grid ...'
     CALL lat2matrix(lat_para,lat_mat,2)
     CALL lat2matrix(reclat_para,recip_lat,2)
!print*,'lat_mat'
!print*,lat_mat
!STOP
     !average charge density
     charge_ave=ncharge/volume
      !creat grids
      CALL Build_rgrid()
      !creat k-points
      CALL Init_BandStruct(numk,kvec)
      !>>>Finite difference
      CALL init_finite(finite_order,gap)
      !<<<Finite difference
      !Plan for FFT
#ifdef MPI
      CALL PlanFFT(global_n1,global_n2,global_n3)
#else
      CALL PlanFFT(n1,n2,n3)
#endif
      !recip grid mesh data
      CALL FillQTable()
      !real grid mesh data
      CALL FillRTable()
      !perpare vlpp
      CALL vlpp()
      !perpare nolocal part
#ifdef MPI
      CALL initialize_nlpot_band()
#else
      CALL initialize_nlpot()
#endif
      ! CALL initialize_nlpot()
#ifdef MPI
      if(parallel%isroot)THEN
#endif
      WRITE(*,*)'R-space-GRIDS:',n1,n2,n3
      WRITE(*,*)'K-points-used:',nk
#ifdef MPI
      endif
#endif
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE Band_Begin
   !--------------------band structure-----------------------
   SUBROUTINE build_bands()
      USE Struct_module , ONLY : ncharge,struct,volume,lat_para
#ifdef MPI
      USE parameters , ONLY : Nspin,ETOL,system_name &
                      &, nev=>Nstates_global, ISTART
      USE smpi_math_module, ONLY:parallel
#else
      USE parameters , ONLY : Nspin,ETOL,system_name &
           &, nev=>Nstates, ISTART
#endif
      USE grid_module , ONLY : grid,nk,KPT
      IMPLICIT NONE
      INTEGER(I4B) :: Ik,Ien
      REAL(DP),ALLOCATABLE :: eigval(:,:)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !KPTB%NEV=ncharge/2+NADDSTATES
      !read band path
      IF(ISTART==1)THEN
         !vasp type
         CALL read_bandpath('KPOINTS.dat')
      ELSEIF(ISTART==2)THEN
         !castep type
         CALL Band_pathread('KPOINTS.dat')
      ELSE
         STOP 'BandStruct: ISTART should be 1 or 2'
      ENDIF
      !read charge density
      CALL read_density('ares.chg',grid%rhoS(:,:,:,1))
      !caleigenval
      ALLOCATE(eigval(nev,nk))
#ifdef MPI
      if(parallel%isroot)PRINT*,'BUILDING:EIGENVAL FILE ...'
#else
      PRINT*,'BUILDING:EIGENVAL FILE ...'
#endif
      CALL cal_band(grid%rhoS,nev,eigval)
#ifdef MPI
      if(parallel%isroot)then
         OPEN(105,FILE='ares.band.dat')
         WRITE(105,11) struct%nati,Nspin,1
         WRITE(105,12) volume*bohr2ang**3,lat_para(1:3)*bohr2ang*0.1d0,0.5e-15
         WRITE(105,13) ETOL
         WRITE(105,14) 'CAR'
         WRITE(105,"(2X,A20)")  system_name
         WRITE(105,15) ncharge,nk,nev
         WRITE(105,*)
         DO Ik=1,nk

            WRITE(105,16) KPT%vec(:,Ik) !,KPTB%WTK(IPT)
            !output eigenvalue
            DO Ien=1,nev
               WRITE(105,17) Ien,eigval(Ien,Ik)*hart2ev
            ENDDO
            WRITE(105,*)

         ENDDO
         CLOSE(105)
11       FORMAT(1X,I4,1X,I4,1X,I4)
12       FORMAT(E15.7,1X,E15.7,1X,E15.7,1X,E15.7,1X,E15.7)
13       FORMAT(E15.7)
14       FORMAT(2X,A3)
15       FORMAT(2X,I4,1X,I4,1X,I4)
         !16 FORMAT(1X,F10.7,1X,F10.7,1X,F10.7,1X,F10.7)
16       FORMAT(1X,F10.7,1X,F10.7,1X,F10.7)
17       FORMAT(1X,I4,1X,F15.7)
         PRINT*,'DONE!'
         ! OPEN(106,FILE='BNDKPOINTS')
         ! DO Ik=1,nk
         !    WRITE(106,*) KPT%vec(:,Ik)
         ! ENDDO
         ! CLOSE(106)
      endif
#else
      OPEN(105,FILE='EIGENVAL')
          WRITE(105,11) struct%nati,Nspin,1
          WRITE(105,12) volume*bohr2ang**3,lat_para(1:3)*bohr2ang*0.1d0,0.5e-15
          WRITE(105,13) ETOL
          WRITE(105,14) 'CAR'
          WRITE(105,"(2X,A20)")  system_name
          WRITE(105,15) ncharge,nk,nev
          WRITE(105,*)
          DO Ik=1,nk

             WRITE(105,16) KPT%vec(:,Ik) !,KPTB%WTK(IPT)
             !output eigenvalue
             DO Ien=1,nev
                WRITE(105,17) Ien,eigval(Ien,Ik)*hart2ev
             ENDDO
             WRITE(105,*)

          ENDDO
      CLOSE(105)
      11 FORMAT(1X,I4,1X,I4,1X,I4)
      12 FORMAT(E15.7,1X,E15.7,1X,E15.7,1X,E15.7,1X,E15.7)
      13 FORMAT(E15.7)
      14 FORMAT(2X,A3)
      15 FORMAT(2X,I4,1X,I4,1X,I4)
      !16 FORMAT(1X,F10.7,1X,F10.7,1X,F10.7,1X,F10.7)
      16 FORMAT(1X,F10.7,1X,F10.7,1X,F10.7)
      17 FORMAT(1X,I4,1X,F15.7)
      PRINT*,'DONE!'
      OPEN(106,FILE='BNDKPOINTS')
         DO Ik=1,nk
            WRITE(106,*) KPT%vec(:,Ik)
         ENDDO
      CLOSE(106)
#endif
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE build_bands
   !---------------------read BANDPATH-----------------------
   SUBROUTINE read_bandpath(infile)
      USE parameters , ONLY : NEV=>Nstates
#ifdef MPI
      USE smpi_math_module, ONLY:parallel
#endif
      IMPLICIT NONE
      !IN/OUT
      CHARACTER(11),INTENT(IN) :: infile !band path
      !LOCAL
      CHARACTER(70)        :: SNAMK !name
      CHARACTER(1) :: CSEL,CLINE
      INTEGER(I4B) :: NINTER,NKP,IERR,IINDEX,N
      REAL(DP),ALLOCATABLE :: VKPT2(:,:) , VKPT(:,:)
      REAL(DP) :: SHIFT(3)
      !
      INTEGER(I4B) :: NKPTS ,  &!Number of k-points tot path
               &      NLine ,  &!number of lines
               &      NKPT !number of k-points per line
      LOGICAL  :: lexist
      INTEGER(I4B),PARAMETER :: NKDIMD=1000 !maxk
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !-------start input.dat-----------
      INQUIRE(FILE=infile,EXIST=lexist)
      !test the input.dat
      IF(.NOT.lexist)THEN
         WRITE(*,*) '>>>WARNING<<<: KPOINTS.dat file is not exist.stop!!!'
         STOP
      ENDIF
      !
      OPEN(103,FILE=infile,STATUS='OLD',IOSTAT=IERR)

         !READ(103,'(A40)',ERR=70111,END=70111) KPTB%SNAMK
         READ(103,'(A70)') SNAMK
#ifdef MPI
         if(parallel%isroot)WRITE(*,*) 'KPOINTS :',SNAMK
#else
         WRITE(*,*) 'KPOINTS :',SNAMK
#endif


         READ(103,*) NKPT
         IF(NKPT>100)THEN
             WRITE(*,*) 'KPONITS : k-points number can not exceed 100'
             STOP
         ENDIF

         READ(103,'(A1)') CSEL
         IF (CSEL=='L'.OR.CSEL=='l') THEN

            CLINE='L'
            READ(103,'(A1)') CSEL
            NKPT=MAX( NKPT,2)
#ifdef MPI
         if(parallel%isroot)WRITE(*,*)'KPONITS : interpolating k-points between supplied coordinates'
#else
         WRITE(*,*)'KPONITS : interpolating k-points between supplied coordinates'
#endif

         ELSE
            CLINE=" "
#ifdef MPI
         if(parallel%isroot)WRITE(*,*) 'KPONITS.dat: Line must,STOP'
#else
            WRITE(*,*) 'KPONITS.dat: Line mode,STOP'
#endif
            WRITE(*,*) 'KPONITS.dat: Line mode,STOP'
            STOP
         ENDIF
         !-------------------------------------------------------------
         ALLOCATE(VKPT2(3,NKDIMD))

         NINTER = NKPT
         NKP=0  !count the number of line-end-points
         DO 
            NKP=NKP+1
            IF (NKP>NKDIMD) THEN
#ifdef MPI
               if(parallel%isroot)WRITE(*,*)'ERROR in RD_KPOINTS (mkpoints.F): increase NKDIMD to',NKP
#else
               WRITE(*,*)'ERROR in RD_KPOINTS (mkpoints.F): increase NKDIMD to',NKP
#endif
               STOP
            ENDIF

            READ(103,*,IOSTAT=IERR) &
     &           VKPT2(1,NKP),VKPT2(2,NKP),VKPT2(3,NKP)
            IF (IERR/=0) EXIT
         ENDDO
         !NKPTS is used to store the lines-end number now
         NKPTS=NKP-1

         IINDEX=0
         ! make NKPTS even
         NKPTS=(NKPTS/2)*2
         Nline=NKPTS/2
         !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
         ALLOCATE(VKPT(3,Nline*NKPT))
         !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
         DO NKP=1,NKPTS-1,2
            SHIFT=(VKPT2(:,NKP+1)-VKPT2(:,NKP))/(NINTER-1)
            DO N=0,NINTER-1
               IINDEX=IINDEX+1
               IF (IINDEX>NKDIMD) THEN
#ifdef MPI
                  if(parallel%isroot)WRITE(*,*)'ERROR in RD_KPOINTS (mkpoints.F): increase NKDIMD to',IINDEX
#else
                  WRITE(*,*)'ERROR in RD_KPOINTS (mkpoints.F): increase NKDIMD to',IINDEX
#endif
                  STOP
               ENDIF
               VKPT(:,IINDEX)=VKPT2(:,NKP)+SHIFT*N
            ENDDO
         ENDDO

         NKPTS=Nline*NINTER
      CLOSE(103)
      CALL Band_Begin(NKPTS,VKPT)
      !destory the array
      DEALLOCATE(VKPT,VKPT2)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !ALLOCATE(KPTB%WTK(KPTB%NKPTS))
      !!nouse,just suit for vasp
      !KPTB%WTK(:)=1.d0/KPTB%NKPTS
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !print*,KPTB%WTK
#ifdef MPI
      if(parallel%isroot)PRINT*,'KPOINTS : PATH READ DONE!!!'
#else
      PRINT*,'KPOINTS : PATH READ DONE!!!'
#endif
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE read_bandpath
   !---------------------Band_pathread-----------------------
   SUBROUTINE Band_pathread(infile)
      USE parameters , ONLY : NEV=>Nstates
#ifdef MPI
      USE smpi_math_module, ONLY: parallel
#endif
      IMPLICIT NONE
      CHARACTER(11),INTENT(IN) :: infile !band path
      !LOCAL
      CHARACTER(70)        :: SNAMK !name
      INTEGER(I4B) :: numk,Ik,Itmp,IERR
      REAL(DP),ALLOCATABLE :: VKPT(:,:)
      LOGICAL  :: lexist
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      INQUIRE(FILE=infile,EXIST=lexist)
      !test the input.dat
      IF(.NOT.lexist)THEN
         WRITE(*,*) '>>>WARNING<<<: KPOINTS.dat file is not exist.stop!!!'
         STOP
      ENDIF

      OPEN(103,FILE=infile,STATUS='OLD',IOSTAT=IERR)
         !nouse : just title
         READ(103,'(A70)') SNAMK
#ifdef MPI
         if(parallel%isroot)WRITE(*,*) 'KPOINTS :',SNAMK
#else
         WRITE(*,*) 'KPOINTS :',SNAMK
#endif
         !read k number
         READ(103,*) numk
         !ALLOCATE
         ALLOCATE(VKPT(3,numk))
         DO Ik=1,numk
            READ(103,*) Itmp,VKPT(:,Ik)
         ENDDO
      CLOSE(103)
      CALL Band_Begin(numk,VKPT)
      !deallocate
      DEALLOCATE(VKPT)
#ifdef MPI
      if(parallel%isroot)PRINT*,'KPOINTS : PATH READ DONE!!!'
#else
      PRINT*,'KPOINTS : PATH READ DONE!!!'
#endif
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE Band_pathread
   !----------------------read CHGCAR------------------------
   SUBROUTINE read_density(infile,rho)
      USE constants
#ifdef MPI
      USE smpi_math_module, ONLY:parallel,mpi_real8,mpinfo
      USE grid_module, ONLY:global_n3,global_n2,global_n1,global_n
#endif
      IMPLICIT NONE
      !IN/OUT
      CHARACTER(8),INTENT(IN) :: infile !CHGCAR
      REAL(DP),INTENT(OUT)    :: rho(:,:,:) !density
#ifdef MPI
      REAL(DP)    :: rho_global(global_n1,global_n2,global_n3) !density
      INTEGER(I4B) :: ix,iy
#endif
      !LOCAL
      INTEGER(I4B) :: IERR ,&
                 &    n(3) 
      LOGICAL :: lexist 
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      INQUIRE(FILE=infile,EXIST=lexist)
      print*,'infile',infile,'exist',Lexist
      !test the input.dat
      IF(.NOT.lexist)THEN
#ifdef MPI
         if(parallel%isroot)WRITE(*,*) '>>>WARNING<<<: CHGCAR-file is not exist.stop!!!'
#else
         WRITE(*,*) '>>>WARNING<<<: CHARGE-file(ares.chg) is not exist.stop!!!'
#endif
         STOP
      ENDIF

#ifdef MPI
      if(parallel%isroot)PRINT*, 'READ CHGCAR ... '
#else
      PRINT*, 'READ CHGCAR ... '
#endif
#ifdef MPI
      IF(parallel%isroot)THEN
         OPEN(104,FILE=infile,STATUS='OLD',IOSTAT=IERR)

         READ(104,*) n(:)
         !check
         IF( ( global_n1/=n(1) ) .OR. ( global_n2/=n(2) ) .OR. ( global_n3/=n(3) ) )THEN
             WRITE(*,*) 'CHGCAR: grids number is not suit for GRIDGAP in input.dat'
             STOP
         ENDIF

         READ(104,*) rho_global(:,:,:)

         CLOSE(104)
      ENDIF
      CALL MPI_BCAST(rho_global,global_n,MPI_REAL8,parallel%rootid&
           &,parallel%comm,mpinfo)
      rho=0.d0
      ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
      iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
      call MPI_SCATTERV(rho_global,parallel%recvcounts,parallel%displs,MPI_REAL8&
           & ,rho(ix,iy,1),parallel%mygrid_range(3),MPI_REAL8,parallel%rootid&
           & ,parallel%commx,mpinfo)
#else
      OPEN(104,FILE=infile,STATUS='OLD',IOSTAT=IERR)

         READ(104,*) n(:)
         !check
         IF( ( n1/=n(1) ) .OR. ( n2/=n(2) ) .OR. ( n3/=n(3) ) )THEN
             WRITE(*,*) 'CHGCAR: grids number is not suit for GRIDGAP in input.dat'
             STOP
         ENDIF

         READ(104,*) rho(:,:,:)

         CLOSE(104)
#endif
#ifdef MPI
         if(parallel%isroot)PRINT*, 'DONE!'
#else
         PRINT*, 'DONE!'
#endif
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE read_density
   !-----------------calculate band struct-------------------
   SUBROUTINE cal_band(rhoS,nev,eigval)
      !
      USE constants
      !USE parameters, ONLY : Lsub,LONE,Nsub
      USE parameters, ONLY : Nspin,ETOL,hart2ev
      !USE subsystem_module , ONLY : sub
      !USE subscf_module , ONLY : v_ksced
      USE math , ONLY : sort_eigval
      USE potential_module , ONLY : cal_veff
      USE grid_module , ONLY : nk
#ifdef MPI
      USE smpi_math_module, ONLY:parallel,mpi_real8,mpinfo,mpi_sum
      USE grid_module, ONLY: global_n1,global_n2,global_n3,global_n
      USE Arpack_module , ONLY : diagH_arpack_band,diagH_arpack
      USE pspot_module, ONLY: max_nproj
      USE m_time_evaluate, ONLY: filename
#else
      USE Arpack_module , ONLY : diagH_arpack
#endif
      !IN/OUT
      REAL(DP),INTENT(IN) :: rhoS(:,:,:,:)
      INTEGER(I4B),INTENT(IN) :: nev
      REAL(DP),INTENT(OUT) :: EIGVAL(:,:)
      !LOCAL
      INTEGER(I4B) :: IKP,info,isc,maxmvs,nec
      REAL(DP) :: veff(n1,n2,n3,Nspin) ! , &
             !&    vtemp(n1,n2,n3)
#ifdef MPI
      COMPLEX(DCP) :: resid_int(global_n)
      COMPLEX(DCP) :: psi(global_n,nev)
      REAL(DP)    :: veff_global(global_n1,global_n2,global_n3,Nspin)
      INTEGER(I4B) :: ix,iy
      REAL(DP)     :: EIGVAL_local(nev,nk)
#else
      COMPLEX(DCP) :: resid_int(n)
      COMPLEX(DCP) :: psi(n,nev)
#endif
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL cal_veff(rhoS,veff)

      info=0
      resid_int(:)=(0.d0,0.d0)
#ifdef MPI
     ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
     iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
     call MPI_ALLGATHERV(veff(ix,iy,1,1),parallel%mygrid_range(3),&
          & MPI_REAL8,veff_global,parallel%recvcounts,parallel%displs&
          &, MPI_REAL8, parallel%commx,mpinfo)
      EIGVAL_local=0.d0
#endif
      DO IKP=1,nk
#ifdef MPI
         if(mod(IKP,parallel%numprocs)==parallel%myid)then
         diagH:DO isc=1,5
                  maxmvs=30000*isc
                  CALL diagH_arpack_band(  veff_global(:,:,:,1),IKP,nev,psi,EIGVAL_local(:,IKP), &
                     &     resid_int,nec,info,maxmvs,ETOL)
               ENDDO diagH
         CALL sort_eigval(nev,EIGVAL_local(:,IKP))
         endif
#else
         diagH:DO isc=1,5
            maxmvs=30000*isc
            CALL diagH_arpack(  veff(:,:,:,1),IKP,nev,psi,EIGVAL(:,IKP), &
                 &     resid_int,nec,info,maxmvs,ETOL)
         ENDDO diagH
         CALL sort_eigval(nev,EIGVAL(:,IKP))
#endif
      ENDDO
#ifdef MPI
      call mpi_allreduce(EIGVAL_local, EIGVAL, size(EIGVAL), mpi_real8, mpi_sum, parallel%comm, mpinfo)
#endif
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
   ENDSUBROUTINE cal_band
   !---------------------------------------------------------
ENDMODULE band_structure

# 1 "Struct_module.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "Struct_module.f90"
!!##########################################################!!
!!*************  author: Qiang Xu & Lantian Xue  ***********!!
!!*************  date  :      2021-07-31         ***********!!
!!*************      All Data  Information       ***********!!
!!##########################################################!!
!---------------------------PARTING LINE--------------------------
MODULE struct_module
!!####################################################!!
!!*************  author: Qiang Xu          ***********!!
!!*************  date  : 2017-07-11        ***********!!
!!*************  Save Structure Information        ***!!
!!####################################################!!
   USE constants
   IMPLICIT NONE
   TYPE struct_type
      REAL(DP),ALLOCATABLE      :: Zion (:)
      INTEGER(I4B),ALLOCATABLE  :: nati(:)  !number of atoms for each atom type
      INTEGER(I4B),ALLOCATABLE  :: eleid(:)   !atoms for each atom type id
      REAL(DP),ALLOCATABLE      :: pos(:,:)  !direct
      REAL(DP),ALLOCATABLE      :: poscar(:,:) ! car
      REAL(DP)                  :: stress(3,3)
      REAL(DP),ALLOCATABLE      :: forces(:,:)       ! total forces
      REAL(DP),ALLOCATABLE      :: mass(:)       ! atomic mass
!for STO
      REAL(DP),ALLOCATABLE      :: zeta(:,:)    ! zeta for STO
      INTEGER(I4B),ALLOCATABLE  :: prinq(:,:)   !atom's principle quantum
      INTEGER(I4B),ALLOCATABLE  :: Lmax(:)   !max L need to consider
      CHARACTER(len=3),ALLOCATABLE :: elements(:)
   END TYPE struct_type
!-------------------------------------------------------------
   INTEGER(I4B)              :: natom,nzion! total number of atoms
   INTEGER(I4B)              :: naty ! number of atom types
   REAL(DP)                  :: ncharge    ! total change of the system
   REAL(DP)                  :: charge_ave ! ave change of the system
   REAL(DP)                  :: volume,volsp
   REAL(DP)                  :: lat_mat(3,3)
   REAL(DP)                  :: lat_para(6)
   REAL(DP)                  :: recip_lat(3,3) !!! the lattice matrix in reciprocal space
   REAL(DP)                  :: reclat_para(6)
   REAL(DP)                  :: Eionion ! energy of ion-ion
   REAL(DP)                  :: Eshift_ps,Eshift_tot !total energy shifts
!cell symmetry
   REAL(DP)                  :: Opsym(3,3,48) !operator of symmetry
   REAL(DP)                  :: Otrans(3,48)
   INTEGER(I4B)              :: nsym  !number of symmetry we used
   INTEGER(I4B)              :: num_t  !number of translations we used
   INTEGER(I4B)              :: c_i(8)
   INTEGER(I4B)              :: Odet(48)
   TYPE(struct_type) :: struct
CONTAINS
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!--------------------------------------------------------------
   SUBROUTINE creat_struct(numtyp,numatom)
      USE m_time_evaluate, only: memory_sum
      IMPLICIT NONE
!
      INTEGER(I4B),INTENT(IN) :: numtyp
      INTEGER(I4B),INTENT(IN) :: numatom
!TYPE(struct_type),INTENT(OUT) :: struct
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      naty=numtyp
      natom=numatom
      CALL destroy_struct()
      ALLOCATE(struct%Zion(numtyp))
      ALLOCATE(struct%nati(numtyp))
      ALLOCATE(struct%eleid(numtyp+1))
      ALLOCATE(struct%pos(3,natom))
      ALLOCATE(struct%poscar(3,natom))
      ALLOCATE(struct%forces(3,natom))
      ALLOCATE(struct%mass(numtyp))
      ALLOCATE(struct%elements(numtyp))
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      call memory_sum("creat_struct",real(size(struct%Zion),DP)*DP&
           &+size(struct%nati)*I4B+size(struct%eleid)*I4B&
           &+size(struct%pos)*DP+size(struct%poscar)*DP&
           &+size(struct%forces)*DP+size(struct%mass)*DP&
           &+size(struct%Lmax)*I4B+size(struct%elements)*3)
   ENDSUBROUTINE creat_struct
!--------------------------------------------------------------
   SUBROUTINE destroy_struct()
     USE m_time_evaluate, ONLY: memory_free
      IMPLICIT NONE
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF(ALLOCATED(struct%Zion))      DEALLOCATE(struct%Zion)
      IF(ALLOCATED(struct%nati))      DEALLOCATE(struct%nati)
      IF(ALLOCATED(struct%eleid))     DEALLOCATE(struct%eleid)
      IF(ALLOCATED(struct%pos))       DEALLOCATE(struct%pos)
      IF(ALLOCATED(struct%poscar))    DEALLOCATE(struct%poscar)
      IF(ALLOCATED(struct%forces))    DEALLOCATE(struct%forces)
      IF(ALLOCATED(struct%mass))    DEALLOCATE(struct%mass)
      IF(ALLOCATED(struct%zeta))    DEALLOCATE(struct%zeta)
      IF(ALLOCATED(struct%prinq))    DEALLOCATE(struct%prinq)
      IF(ALLOCATED(struct%Lmax))    DEALLOCATE(struct%Lmax)
      IF(ALLOCATED(struct%elements))  DEALLOCATE(struct%elements)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE destroy_struct
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDMODULE struct_module
!---------------------------PARTING LINE--------------------------
MODULE pspot_module
   USE constants
   IMPLICIT NONE
!defined
   TYPE pspot
      REAL(DP) :: Zion
      INTEGER(I4B) :: numps
      INTEGER(I4B) :: qnumps
      CHARACTER(len=3) :: elename
!> local part
      REAL(DP)     :: qmax
      REAL(DP)     :: qspacing
      REAL(DP),ALLOCATABLE :: qmesh(:)
      REAL(DP),ALLOCATABLE :: Vlocq(:)
      REAL(DP),ALLOCATABLE :: VlocqS(:)
      REAL(DP),ALLOCATABLE :: ddVl_dq2(:)
!> for real local pseudopotential used in confined system
      REAL(DP),ALLOCATABLE :: r(:),rab(:)
      REAL(DP),ALLOCATABLE :: Vlocr(:)
!> nonlocal real space
      INTEGER(I4B) :: nproj
      INTEGER(I4B),ALLOCATABLE :: proj_l(:)
      INTEGER(I4B),ALLOCATABLE :: proj_m(:),indx(:)
      REAL(DP) :: rcut
      REAL(DP) :: rmax
      REAL(DP) :: rspacing
      REAL(DP),ALLOCATABLE :: Dij(:,:)
      REAL(DP),ALLOCATABLE :: beta_r(:,:)
      REAL(DP),ALLOCATABLE :: dbeta_dr(:,:)
!radial charge density
      LOGICAL :: lden
      REAL(DP),ALLOCATABLE :: denr(:)
      REAL(DP),ALLOCATABLE :: dden_dr(:)
!partial-core charge density
      LOGICAL :: lcore
      REAL(DP),ALLOCATABLE :: denc(:) &
                           &, ddenc_dr(:)
!For Self Energy correction
      REAL(DP):: rnoverlap
!Total PS/AE energy
      REAL(DP) :: eps,eae
!Total atomic PS-wave functions
      INTEGER(I4B) :: nwfa
      INTEGER(I4B),ALLOCATABLE :: wfal(:)
      REAL(DP),ALLOCATABLE :: wfar(:,:)
   ENDTYPE pspot
!
   TYPE(pspot),ALLOCATABLE :: psp(:)
!
   INTEGER(I4B) :: max_nproj,max_nwfa
   REAL(DP)     :: max_rcut
!> character array
!> reserve the attribute for UPF
   INTERFACE get_value
     MODULE PROCEDURE get_value_int,get_value_real &
                & ,get_value_char,get_value_logic
   END INTERFACE
   TYPE attribute
     CHARACTER(len=150)  :: value 
   ENDTYPE attribute
   TYPE(attribute),ALLOCATABLE :: attribute_data(:)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
CONTAINS
!---------------------------PARTING LINE--------------------------
   SUBROUTINE read_pspot(nty,filenames)
     USE parameters , ONLY : Nstates,Nstates_global,nspin,Naddstates & 
          & ,dcharge,LAtomRho,Lcore_val,ixc
     USE struct_module, ONLY:struct,nzion,ncharge,Eshift_ps,Eshift_tot

     USE smpi_math_module
     USE m_time_evaluate, ONLY: memory_sum

     IMPLICIT NONE
!IN/OUT
     INTEGER(I4B),INTENT(IN)  :: nty
     CHARACTER(30),INTENT(IN) :: filenames(nty)
!LOCAL
     INTEGER(I4B) :: Ity
     INTEGER(I4B) :: Ipj,l,m , i  ,mnp

     INTEGER(I4B) :: N_comm

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!for every speice atom
     ALLOCATE(psp(nty))

      IF(parallel%isroot)THEN

!read atomic pseudopotentials
     DO Ity=1,nty
        CALL read_psupf_atom(Ity,filenames(Ity),psp(Ity))
        struct%Zion(Ity)=NINT(psp(Ity)%Zion)
        IF(psp(Ity)%lcore) Lcore_val=.TRUE.
        IF((.NOT.psp(Ity)%lden)) LAtomRho=.FALSE.
     ENDDO
!max rcut and nproj
     max_rcut=MAXVAL(psp(:)%rcut)
     max_nproj=MAXVAL(psp(:)%nproj)
     max_nwfa=MAXVAL(psp(:)%nwfa)
!calculate the total charge in systems
     nzion=SUM(struct%nati(:)*struct%Zion(:))
!(dcharge: remove/add charge)
     ncharge=REAL(nzion,DP) + dcharge
!states we need
     IF(Naddstates>0)THEN
        Nstates=NINT(ncharge/2)+Naddstates
     ELSE
        Nstates=MAX(NINT(ncharge/2*1.1),NINT(ncharge/2+15))
     ENDIF

     Nstates=MAX(Nstates,parallel%numprocs+1)

!Energy Shifts
     Eshift_ps=-SUM(struct%nati(:)*psp(:)%eps)
     Eshift_tot=Eshift_ps+SUM(struct%nati(:)*psp(:)%eae)
     WRITE(*,*) '[Ion and Charge Num.]', nzion,ncharge
     WRITE(*,*) '[Eigen States Num.]', Nstates
!Check XC functional
     IF(ixc==1)THEN
        WRITE(*,*) 'XC: PZ-LDA (Perdew-Zunger 1981)'
     ELSEIF(ixc==-1)THEN
        WRITE(*,*) 'XC: PW91-GGA (Perdew-Wang 1991)'
     ELSEIF(ixc==-2)THEN
        WRITE(*,*) 'XC: PBE-GGA (Pewdew-Burke-Ernzerhof 1996)'
     ELSEIF(ixc==-3)THEN
        WRITE(*,*) 'XC: BLYP-GGA (B88 + LYP)'
     ELSE
        WRITE(*,*) 'STOP: Please Check The XC Used in Pseudopotential Files'
        STOP
     ENDIF

     ENDIF
!parallel boardcast
     CALL MPI_BCAST(Lcore_val        ,1,MPI_LOGICAL,parallel%rootid,parallel%comm,mpinfo)
     CALL MPI_BCAST(LAtomRho         ,1,MPI_LOGICAL,parallel%rootid,parallel%comm,mpinfo)
     CALL MPI_BCAST(max_rcut         ,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
     CALL MPI_BCAST(ncharge          ,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
     CALL MPI_BCAST(Eshift_ps        ,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
     CALL MPI_BCAST(Eshift_tot       ,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
     CALL MPI_BCAST(max_nproj        ,1,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
     CALL MPI_BCAST(max_nwfa         ,1,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
     CALL MPI_BCAST(nzion            ,1,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
     CALL MPI_BCAST(Nstates          ,1,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
     CALL MPI_BCAST(ixc              ,1,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)

!split the states for MatMat
     Nstates_global=Nstates
     CALL array_split(Nstates_global)
     Nstates=parallel%nstate_proc

     DO Ity=1,nty
        CALL MPI_BCAST(struct%Zion(Ity)   ,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
        CALL MPI_BCAST(psp(Ity)%Zion      ,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
        CALL MPI_BCAST(psp(Ity)%rcut      ,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
        CALL MPI_BCAST(psp(Ity)%rmax      ,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
        CALL MPI_BCAST(psp(Ity)%qmax      ,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
        CALL MPI_BCAST(psp(Ity)%qspacing  ,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
        CALL MPI_BCAST(psp(Ity)%rnoverlap ,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
        CALL MPI_BCAST(psp(Ity)%numps     ,1,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
        CALL MPI_BCAST(psp(Ity)%qnumps    ,1,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
        CALL MPI_BCAST(psp(Ity)%nproj     ,1,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
        CALL MPI_BCAST(psp(Ity)%nwfa      ,1,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
        CALL MPI_BCAST(psp(Ity)%lden      ,1,MPI_LOGICAL,parallel%rootid,parallel%comm,mpinfo)
        CALL MPI_BCAST(psp(Ity)%lcore     ,1,MPI_LOGICAL,parallel%rootid,parallel%comm,mpinfo)
!real-space mesh (r)
        IF(.NOT.ALLOCATED(psp(Ity)%r))THEN
           ALLOCATE(psp(Ity)%r(psp(Ity)%numps))
           call memory_sum('ps_r_real',real(psp(Ity)%numps,DP)*DP)
        ENDIF
        CALL MPI_BCAST(psp(Ity)%r,psp(Ity)%numps,&
                 & MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
!real-space mesh (rab)
        IF(.NOT.ALLOCATED(psp(Ity)%rab))THEN
           ALLOCATE(psp(Ity)%rab(psp(Ity)%numps))
           call memory_sum('ps_rab_real',real(psp(Ity)%numps,DP)*DP)
        ENDIF
        CALL MPI_BCAST(psp(Ity)%rab,psp(Ity)%numps,&
                 & MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
!recip.-space mesh (qmesh) for PBC
        IF(.NOT.ALLOCATED(psp(Ity)%qmesh))THEN
           ALLOCATE(psp(Ity)%qmesh(psp(Ity)%qnumps))
           call memory_sum('ps_qmesh_real',real(psp(Ity)%qnumps,DP)*DP)
        ENDIF
        CALL MPI_BCAST(psp(Ity)%qmesh,psp(Ity)%qnumps,&
                 & MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
!local pseudopotentials, recip. for PBC
        IF(.NOT.ALLOCATED(psp(Ity)%Vlocq))THEN
           ALLOCATE(psp(Ity)%Vlocq(psp(Ity)%qnumps))
           call memory_sum('ps_V_loc',real(psp(Ity)%qnumps,DP)*DP)
        ENDIF
        CALL MPI_BCAST(psp(Ity)%Vlocq,psp(Ity)%qnumps,&
                 & MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
        IF(.not.allocated(psp(Ity)%VlocqS))   ALLOCATE(psp(Ity)%VlocqS(psp(Ity)%qnumps))
        IF(.not.allocated(psp(Ity)%ddVl_dq2)) ALLOCATE(psp(Ity)%ddVl_dq2(psp(Ity)%qnumps))
        CALL MPI_BCAST(psp(Ity)%VlocqS,psp(Ity)%qnumps,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
        CALL MPI_BCAST(psp(Ity)%ddVl_dq2,psp(Ity)%qnumps,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
!nonlocal pseudopotentials
        IF(psp(Ity)%nproj>0)THEN
!l-index
           IF(.NOT.ALLOCATED(psp(Ity)%proj_l))THEN
              ALLOCATE(psp(Ity)%proj_l(psp(Ity)%nproj))
              call memory_sum('ps_proj_l',real(psp(Ity)%nproj,DP)*I4B)
           ENDIF
           CALL MPI_BCAST(psp(Ity)%proj_l,psp(Ity)%nproj,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
!m-index
           IF(.NOT.ALLOCATED(psp(Ity)%proj_m))THEN
              ALLOCATE(psp(Ity)%proj_m(psp(Ity)%nproj))
              call memory_sum('ps_proj_m',real(psp(Ity)%nproj,DP)*I4B)
           ENDIF
           CALL MPI_BCAST(psp(Ity)%proj_m,psp(Ity)%nproj,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
!Dij
           IF(.NOT.ALLOCATED(psp(Ity)%Dij))THEN
              ALLOCATE(psp(Ity)%Dij(psp(Ity)%nproj,psp(Ity)%nproj))
              call memory_sum('ps_Dij',real(size(psp(Ity)%Dij),DP)*DP)
           ENDIF
           N_comm=psp(Ity)%nproj**2
           CALL MPI_BCAST(psp(Ity)%Dij,N_comm,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
!Vnl
           IF(.NOT.ALLOCATED(psp(Ity)%beta_r))THEN
              ALLOCATE(psp(Ity)%beta_r(psp(Ity)%numps,psp(Ity)%nproj))
              call memory_sum('ps_betar',real(size(psp(Ity)%beta_r),DP)*DP)
           ENDIF
           N_comm=psp(Ity)%numps*psp(Ity)%nproj
           CALL MPI_BCAST(psp(Ity)%beta_r ,N_comm,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)

        ENDIF
        CALL MPI_Barrier(parallel%comm,mpinfo)
!core charge
        IF(psp(Ity)%lcore)THEN
!core density
           IF(.NOT.ALLOCATED(psp(Ity)%denc))THEN
              ALLOCATE(psp(Ity)%denc(psp(Ity)%numps))
              call memory_sum('ps_denc',real(psp(Ity)%numps,DP)*DP)
           ENDIF
           CALL MPI_BCAST(psp(Ity)%denc ,psp(Ity)%numps,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
   
        ENDIF
        CALL MPI_Barrier(parallel%comm,mpinfo)
!density
        IF(psp(Ity)%lden)THEN
!ps-valence density
           IF(.NOT.ALLOCATED(psp(Ity)%denr))THEN
              ALLOCATE(psp(Ity)%denr(psp(Ity)%numps))
              call memory_sum('ps_denr',real(psp(Ity)%numps,DP)*DP)
           ENDIF
           CALL MPI_BCAST(psp(Ity)%denr ,psp(Ity)%numps,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
        ENDIF
        CALL MPI_Barrier(parallel%comm,mpinfo)
!wvfs
        IF(psp(Ity)%nwfa>0)THEN !wfa
!l-index
           IF(.NOT.ALLOCATED(psp(Ity)%wfal))THEN
              ALLOCATE(psp(Ity)%wfal(psp(Ity)%nwfa))
              call memory_sum('ps_wfal',real(psp(Ity)%nwfa,DP)*I4B)
           ENDIF
           CALL MPI_BCAST(psp(Ity)%wfal ,psp(Ity)%nwfa,MPI_INTEGER4,parallel%rootid,parallel%comm,mpinfo)
!wvf(r)
           IF(.NOT.ALLOCATED(psp(Ity)%wfar))THEN
              ALLOCATE(psp(Ity)%wfar(psp(Ity)%numps,psp(Ity)%nwfa))
              call memory_sum('ps_wfar',real( size(psp(Ity)%wfar),DP)*DP )
           ENDIF
           N_comm=psp(Ity)%numps*psp(Ity)%nwfa
           CALL MPI_BCAST(psp(Ity)%wfar ,N_comm,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)

        ENDIF !wfa
        CALL MPI_Barrier(parallel%comm,mpinfo)

     ENDDO

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE read_pspot
!---------------------------PARTING LINE--------------------------
   SUBROUTINE read_psupf_atom(Ity,filename,ps)
      USE Math,        ONLY : dfdr,fourier_1d
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
!IN/OUT
      INTEGER(I4B),INTENT(IN)   :: Ity
      CHARACTER(30),INTENT(IN) :: filename
!LOCAL
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
           &,temp_beta_r_interp(:),temp_r(:)
      REAL(DP),ALLOCATABLE :: knots(:)
      LOGICAL :: lpsp !quire the psp file
      INTEGER(I4B) :: ztmp
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      PRINT*,'[Reading PseudoPP FILE]',filename
      INQUIRE(FILE=TRIM(ADJUSTL(filename)),EXIST=lpsp)
      IF(.NOT.lpsp)THEN
         STOP 'Check the Pseudopotential Name!!!'
         READ(filename,*) ztmp
         WRITE(*,*) '>>>WARNING<<<: Pseudopotential Not Exist, changed by All-Electron Potential'
         WRITE(*,*) '>>>Element<<<:',ztmp
         CALL AEP_generator(ztmp,ps)
         RETURN
      ENDIF
      unit_UPF=Ity+1600
      OPEN(unit=unit_UPF,file=filename,status='old',form='formatted',iostat=ios)
!WRITE ( *, * ) " Reading pseudopotential file in UPF format..."

!Initialize data
      ps%rcut=0._DP 
!> Search for Header
      CALL scan_head(unit_UPF, "HEADER", .true.)
      CALL read_pseudo_header(ps,temp_nproj)
!> Search for mesh information
!radius
      CALL scan_head (unit_UPF, "R ", .true.)
      ALLOCATE(ps%r(ps%numps))
      READ(unit_UPF,*) ps%r
      ps%rmax=ps%r(ps%numps) !the max radius
!CALL scan_tail (unit_UPF, "MESH")
!rab calculation
      CALL scan_head (unit_UPF, "RAB ", .true.)
      ALLOCATE(ps%rab(ps%numps))
      READ(unit_UPF,*) ps%rab
      CALL scan_tail (unit_UPF, "MESH")
!-------->Search for Local potential
      IF(ps%lcore)THEN
         CALL scan_head (unit_UPF, "NLCC", .true.)
         ALLOCATE(ps%denc(ps%numps))
         READ(unit_UPF,*) ps%denc
!ps%denc=ps%denc/2._DP
         PRINT*,'total core:',SUM(4*pi*ps%denc*ps%r**2*ps%rab)
      ENDIF
!-------->Search for Local potential
      CALL scan_head (unit_UPF, "LOCAL", .true.)
      ALLOCATE(ps%Vlocr(ps%numps))
      READ(unit_UPF,*) ps%Vlocr
!> Ry to hartree
      ps%Vlocr=ps%Vlocr*0.5_DP
      CALL scan_tail (unit_UPF, "LOCAL")
!for recipical space
      ps%qmax=100._DP*bohr2ang
      ps%qnumps=6001
      dq=ps%qmax/(ps%qnumps-1)
      ps%qspacing=dq

      ALLOCATE(ps%Vlocq(ps%qnumps),ps%qmesh(ps%qnumps))
      ps%qmesh=(/(i*dq ,i=0,ps%qnumps-1,1)/)
      CALL fourier_1d(ps%numps,ps%r,ps%rab,ps%Vlocr &
           &,0,ps%qnumps,ps%qmesh,ps%Vlocq,ps%Zion)
!nouse for PBC
      DEALLOCATE(ps%vlocr)
!================parpare to interplote===============
      ALLOCATE(knots(ps%qnumps))
      DO i=1,ps%qnumps
         knots(i)=REAL(i,DP)
      ENDDO
!> local part
      ALLOCATE(ps%VlocqS(ps%qnumps))
      ALLOCATE(ps%ddVl_dq2(ps%qnumps))
!> subtract 1/r in q-space
      ps%VlocqS(1)=ps%Vlocq(1)
      ps%VlocqS(2:)=ps%Vlocq(2:)+ 4.d0*pi*ps%Zion/(ps%qmesh(2:))**2
!> derivetive 2 order
      CALL spline_cubic_set ( ps%qnumps, knots, ps%VlocqS, &
           &    1, 0._DP, 1, 0._DP, ps%ddVl_dq2)
!-------->Search for Nonlocal potential
      CALL scan_head (unit_UPF, "NONLOCAL", .true.)
      DEALLOCATE(attribute_data)
      ALLOCATE(temp_beta_r(ps%numps,temp_nproj))
      ALLOCATE(tempD0(temp_nproj,temp_nproj))
      ALLOCATE(temp_l(temp_nproj))
!READ(unit_UPF,*) ps%beta_r
      CALL read_pseudo_nonlocal(unit_UPF,temp_nproj &
           & ,temp_beta_r,tempD0,temp_rcut,temp_l)
      CALL scan_tail (unit_UPF, "NONLOCAL")
!> nolocal cutoff radius
      if(ps%rcut<temp_rcut) ps%rcut=temp_rcut
!For self-energy correction
      ps%rnoverlap=0.5_DP*ps%rcut
!>reset nonlocal pseudopotential
      temp_a=0
      DO j=1,temp_nproj,1
         temp_a=temp_a+2*temp_l(j)+1
      ENDDO

      ALLOCATE(ps%proj_l(temp_a),ps%proj_m(temp_a),ps%indx(temp_a))
      ALLOCATE(ps%beta_r(ps%numps,temp_a))
      ALLOCATE(ps%Dij(temp_a,temp_a))
! print *,'temp_a',temp_a
! print *,'shape(ps%beta_r)',shape(ps%beta_r)
      ps%beta_r=0.d0
      ps%Dij=0.d0
      ps%nproj = 0
      DO j=1,temp_nproj
!adding the 2l+1 for every l-index
         ps%nproj = ps%nproj+2*temp_l(j)+1 ! 2l+1
!set l-index for projectors
         ps%proj_l(ps%nproj-2*temp_l(j):ps%nproj) = temp_l(j)
!set m-index,m=(-l:l)
         ps%proj_m(ps%nproj-2*temp_l(j):ps%nproj) = &
              (/(m,m=-temp_l(j),temp_l(j))/)
!set origin-index
         ps%indx(ps%nproj-2*temp_l(j):ps%nproj)=j
!set the beta(r) for every projectors
         DO k=ps%nproj-2*temp_l(j),ps%nproj,1
! print*,'l',k,j
!> Ry to hartree
            IF(ps%r(1)==0)then
               ps%beta_r(2:,k) = temp_beta_r(2:,j)/ps%r(2:)*0.5d0
!ps%beta_r(1,k) = temp_beta_r(1,j)/ps%rad(1)*0.5d0
               ps%beta_r(1,k) = ps%beta_r(2,k)
            ELSE
               ps%beta_r(:,k) = temp_beta_r(:,j)/ps%r(:)*0.5d0
            ENDIF

         ENDDO
      ENDDO

!> Ry to hartree
      tempD0=tempD0*2.d0
!> assign DIJ
!CALL direct_productlm(temp_nproj,ps%nproj,temp_l &
!     & ,ps%proj_l,tempD0,ps%Dij)
!
      DO j=1,ps%nproj
      DO k=1,ps%nproj
         IF(ps%proj_l(j)==ps%proj_l(k).and.&
           &ps%proj_m(j)==ps%proj_m(k)) THEN
!set the Dij from D0
            ps%Dij(j,k) = tempD0(ps%indx(j),ps%indx(k))
         ENDIF
      ENDDO
      ENDDO
!> Search for atomic wavefunctions
      IF(ps%nwfa>0)THEN
         CALL scan_head (unit_UPF, "PSWFC", .true.)
         DEALLOCATE(attribute_data)
         ALLOCATE(ps%wfal(ps%nwfa),ps%wfar(ps%numps,ps%nwfa))
         CALL read_pseudo_pswfc(unit_UPF,ps%nwfa,ps%numps,ps%wfal,ps%wfar) 
         CALL scan_tail (unit_UPF, "PSWFC")
!over r
         DO k=1,ps%nwfa
            IF(ps%r(1)==0)THEN
               ps%wfar(2:,k)=ps%wfar(2:,k)/ps%r(2:)
               ps%wfar(1,k)=ps%wfar(2,k)
            ELSE
               ps%wfar(:,k)=ps%wfar(:,k)/ps%r(:)
            ENDIF
         ENDDO
      ENDIF
      
!> Search for atomic charge
      CALL scan_head (unit_UPF, "RHOATOM", .true.)
      ALLOCATE(temp_denr(ps%numps))
!print *,shape(temp_denr)
      READ(unit_UPF,*) temp_denr
      CALL scan_tail (unit_UPF, "RHOATOM")
      ALLOCATE(ps%denr(ps%numps))
      ps%lden=.TRUE.
      if(sum(temp_denr)<=0.1d0)then
         ps%lden=.FALSE.
         write(6,*)'read_upf: error, total RHOATOM small than 0.1,'
         write(6,*)'read_upf: start with uniform density'
      endif
      IF(ps%lden)THEN
         IF(ps%r(1)==0)then
            ps%denr(2:)= temp_denr(2:)/ps%r(2:)**2
!ps%beta_r(1,k) = temp_beta_r(1,j)/ps%rad(1)*0.5d0
            ps%denr(1) = ps%denr(2)
         ELSE
            ps%denr(:) = temp_denr(:)/ps%r(:)**2
         ENDIF
         ps%denr(:)=ps%denr(:)/(4*pi)
      ENDIF
!>the END
!WRITE ( *, * ) "... read UPF file done"
      RETURN
      CLOSE (unit=unit_UPF)
      IF( ALLOCATED(attribute_data) ) DEALLOCATE(attribute_data)
      20 WRITE(6,*)"ERR in read UPF"
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE read_psupf_atom
!---------------------------PARTING LINE--------------------------
    SUBROUTINE scan_head(file_unit,title,start_old)
      IMPLICIT NONE
      INTEGER(I4B) :: file_unit
!> Unit of the input file
      CHARACTER(len=*) :: title
!> Label to be matched
      LOGICAL :: start_old,on_off=.FALSE.
!> Flag: if .false. rewind the file
      CHARACTER(len=150) :: read_string
!> The whole line  from file
      INTEGER(I4B) :: ios
!> I/O state
      INTEGER(I4B) :: error_id
!> the id in error Information, not necessary
      INTEGER(I4B) :: i,counter=0
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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
      IF( ALLOCATED(attribute_data) ) DEALLOCATE(attribute_data)
      ALLOCATE(attribute_data(1:counter))
!start
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
313   print*, 'ERROR:Reading upf, iostate is ', abs(ios)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    END SUBROUTINE scan_head
!---------------------------PARTING LINE--------------------------
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
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF ( ALLOCATED(attribute_data) ) DEALLOCATE(attribute_data)
      ios=0
      DO WHILE (ios==0)
         READ (file_unit, '(A)', iostat = ios, err = 313) read_string
         IF (exist_in("</PP_"//title//">", read_string) ) RETURN
      ENDDO
      IF(ios/=0)print*,"scan end err"
      RETURN
313   print*, 'ERROR:read upf, iostate is ', title
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    END SUBROUTINE scan_tail
!---------------------------PARTING LINE--------------------------
   SUBROUTINE read_pseudo_header(ps,nproj)
     USE parameters, ONLY : ixc
     IMPLICIT NONE
     TYPE(pspot),INTENT(OUT) :: ps
!INTEGER(I4B) :: lmax
!INTEGER(I4B) :: l_local
     INTEGER(I4B),INTENT(OUT) :: nproj
     LOGICAL :: find_flag
     INTEGER(I4B) :: i
     CHARACTER(len=80) :: tempc
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
!>    total_aeenergy   etotae
!>    wfc_cutoff       ecutwfc
!>    rho_cutoff       ecutrho
!>    l_max            lmax
!>    l\_max\_rho        lmax_rho
!>    l_local          lloc
!>    mesh_size        mesh
!>    number\_of\_wfc    nwfc
!>    number\_of\_proj   nbeta
!> </PP_HEADER>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     ps%eps=0._DP
     ps%eae=0._DP
!Store the element name.
     DO i=1,size(attribute_data),1
       call get_value(attribute_data(i)%value,'element=',ps%elename,find_flag)
       IF(find_flag)exit
     ENDDO
!pseudopotential type
     DO i=1,size(attribute_data),1
       call get_value(attribute_data(i)%value,'pseudo_type=',tempc,find_flag)
       IF(find_flag)exit
     ENDDO
!check
     IF( .NOT.(tempc(1:2)=='NC'.OR.tempc(1:2)=='SL') )THEN
        WRITE(*,*) 'STOP: Only NCPP can be used now!!!'
        STOP
     ENDIF
!Store the valance num.
     DO i=1,size(attribute_data),1
       call get_value(attribute_data(i)%value,'z_valence=',ps%Zion,find_flag)
       IF(find_flag)exit
     ENDDO
!Store the num. of points.
     DO i=1,size(attribute_data),1
       call get_value(attribute_data(i)%value,'mesh_size',ps%numps,find_flag)
       IF(find_flag)exit
     ENDDO
!Store the num. of projectors.
     DO i=1,size(attribute_data),1
       call get_value(attribute_data(i)%value,'number_of_proj',nproj,find_flag)
       IF(find_flag)exit
     ENDDO
!XC functional?
     DO i=1,size(attribute_data),1
       call get_value(attribute_data(i)%value,'functional=',tempc,find_flag)
       IF(find_flag)exit
     ENDDO
!XC select
     IF(tempc(1:3)=='LDA'.OR.tempc(1:2)=='PZ')THEN
        ixc=1 !PZ81
     ELSEIF(tempc(1:4)=='PW91')THEN
        ixc=-1 !PW91
     ELSEIF(tempc(1:3)=='PBE'.OR.tempc(1:3)=='GGA')THEN
        ixc=-2 !PBE96
     ELSEIF(tempc(1:4)=='BLYP')THEN
        ixc=-3 !BLYP
     ELSE
        WRITE(*,*) 'STOP: Only PZ and PBE XCs are supported in ARES now'
        STOP
     ENDIF
!Nonlinear core-valence?
     DO i=1,size(attribute_data),1
       call get_value(attribute_data(i)%value,'core_correction=',ps%lcore,find_flag)
       IF(find_flag)exit
     ENDDO
!# of wf
     DO i=1,size(attribute_data),1
       call get_value(attribute_data(i)%value,'number_of_wfc=',ps%nwfa,find_flag)
       IF(find_flag)exit
     ENDDO  
!Energy of pseudoatom
     DO i=1,size(attribute_data),1
       call get_value(attribute_data(i)%value,'total_psenergy=',ps%eps,find_flag)
       IF(find_flag)exit
     ENDDO
!Energy of ae-atom
     DO i=1,size(attribute_data),1
       call get_value(attribute_data(i)%value,'total_aeenergy=',ps%eae,find_flag)
       IF(find_flag)exit
     ENDDO
!Ry. to Ha.
     ps%eps=0.5_DP*ps%eps
     ps%eae=0.5_DP*ps%eae
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   END SUBROUTINE read_pseudo_header
!---------------------------PARTING LINE--------------------------
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
       READ(char_in(ia+1:ib-1),*) variable
! print *, char_find, variable, char_in(ia+1:ib-1)
       find_flag=.TRUE.
     ELSE
       find_flag=.FALSE.
     ENDIF
   END SUBROUTINE get_value_int
!---------------------------PARTING LINE--------------------------
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
!---------------------------PARTING LINE--------------------------
   SUBROUTINE get_value_char(char_in,char_find,variable,find_flag)
     IMPLICIT NONE
     INTEGER(I4B)  :: i,j,ia,ib,foo
     CHARACTER(len=120),INTENT(IN) :: char_in
     CHARACTER(len=*),INTENT(IN) :: char_find
     CHARACTER(len=*)      :: variable
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
       READ(char_in(ia+1:ib-1),*) variable
! print *, char_find, variable, char_in(ia+1:ib-1)
       find_flag=.TRUE.
     ELSE
       find_flag=.FALSE.
     ENDIF
   END SUBROUTINE get_value_char
!---------------------------PARTING LINE--------------------------
   SUBROUTINE get_value_logic(char_in,char_find,variable,find_flag)
     IMPLICIT NONE
     INTEGER(I4B)  :: j,ia,foo
     CHARACTER(len=120),INTENT(IN) :: char_in
     CHARACTER(len=*),INTENT(IN) :: char_find
     LOGICAL      :: variable
     CHARACTER(len=1) :: tmpc
     LOGICAL      :: find_flag
!>read "char_find" from "char_in" and assign to variable
     IF(exist_in(char_find,char_in))THEN
       foo=exist_ibegin(char_find,char_in)
       find_value1: DO j=foo,120,1
         IF(char_in(j:j)=='"')THEN
           ia=j
           exit find_value1
         ENDIF
       ENDDO find_value1

       READ(char_in(ia+1:ia+1),*) tmpc
! print *, char_find, variable, char_in(ia+1:ib-1)
       IF(tmpc=='T')THEN
         variable=.TRUE.
       ELSE
         variable=.FALSE.
       ENDIF
       find_flag=.TRUE.
     ELSE
       find_flag=.FALSE.
     ENDIF
   END SUBROUTINE get_value_logic
!---------------------------PARTING LINE--------------------------
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
!---------------------------PARTING LINE--------------------------
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
!---------------------------PARTING LINE--------------------------
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
         call get_value(attribute_data(j)%value,'angular_momentum=',proj_l(i),find_flag)
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
!---------------------------PARTING LINE--------------------------
   SUBROUTINE read_pseudo_pswfc (unit_UPF,nwfc,nps,wfcl,wfcr)
!---------------------------------------------------------------------
!
     IMPLICIT NONE
!
     INTEGER(I4B),INTENT(IN) :: unit_UPF,nwfc ,nps
     INTEGER(I4B),INTENT(OUT) :: wfcl(:)
     REAL(DP),INTENT(OUT) :: wfcr(:,:)
!LOCAL
     INTEGER :: i,j
     LOGICAL :: find_flag
     CHARACTER(len=20) :: title_in
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
     DO i = 1, nwfc
        write(title_in,'("CHI."I1)')i
        CALL scan_head (unit_UPF,trim(adjustl(title_in)),.FALSE.)
!> get angular_momentum
        DO j=1,size(attribute_data),1
           call get_value(attribute_data(j)%value,' l=',wfcl(i),find_flag)
           IF(find_flag) exit
        ENDDO
        
!read wvfs
        READ(unit_UPF,*) wfcr(:,i)
        DEALLOCATE(attribute_data)
     ENDDO
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   END SUBROUTINE read_pseudo_pswfc
!---------------------------PARTING LINE--------------------------
SUBROUTINE AEP_generator(nz,ps)
   IMPLICIT NONE
   INTEGER(I4B),INTENT(IN) ::  nz ! # of nuclear charge
   TYPE(pspot),INTENT(OUT) :: ps
!LOCAL
   INTEGER(I4B) :: ip
   INTEGER(I4B),PARAMETER :: nps=10 !gridmesh
   REAL(DP),PARAMETER :: par_a=0.000413125362778 &
                    &,par_b=0.01_DP !mesh parameters
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   ps%zion=REAL(nz,DP)
   ps%numps=nps
   ps%nproj=0
   ps%rcut=0
   ps%lden=.FALSE.
   ps%lcore=.FALSE.
   ps%eps=0._DP
   ps%eae=0._DP
   ALLOCATE(ps%r(nps),ps%Vlocr(nps))
   DO ip=nps,1,-1
      ps%r(ip)=par_a*( EXP(par_b*(ip-1)) - 1._DP)
      IF(ip/=1)THEN
         ps%Vlocr(ip)=-nz/ps%r(ip)
      ELSE
         ps%Vlocr(1)=ps%Vlocr(2)
      ENDIF
   ENDDO
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDSUBROUTINE AEP_generator
!---------------------------PARTING LINE--------------------------
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDMODULE pspot_module
!---------------------------PARTING LINE--------------------------
MODULE grid_module
!##############################################!
!*For    : real space grid mesh and relate mesh!
!*Author : Qiang Xu                            !
!*Date   : 2017-7-12                           !
!##############################################!
  USE constants
  IMPLICIT NONE
!defined
  TYPE grid_type
!3D real space grid mesh
     REAL(DP),ALLOCATABLE :: rhoS(:,:)
     REAL(DP),ALLOCATABLE :: rho(:)
     REAL(DP),ALLOCATABLE :: vxcS(:,:)
     REAL(DP),ALLOCATABLE :: vhxcd(:,:)
     REAL(DP),ALLOCATABLE :: vlpp(:)
     REAL(DP),ALLOCATABLE :: vh(:)
     REAL(DP),ALLOCATABLE :: eval(:,:,:) !eigen-values
!partial core density
     REAL(DP),ALLOCATABLE :: rhoc(:)
!recip q table
     REAL(DP),ALLOCATABLE :: gVec(:,:)
     LOGICAL ,ALLOCATABLE :: gMask(:)
!r table ktable
     REAL(DP),ALLOCATABLE :: rVec(:,:)
!> logical for judge if in sphere around atoms
     LOGICAL,ALLOCATABLE :: lsp(:,:,:)
!> sphere grid
     INTEGER(I4B) :: oneDlength
  ENDTYPE grid_type
!k-points data
  TYPE kgrid_type
!
     REAL(DP),ALLOCATABLE :: vec(:,:)
     REAL(DP),ALLOCATABLE :: vcar(:,:)
     REAL(DP),ALLOCATABLE :: wk(:)
  ENDTYPE kgrid_type
!eigen data
  TYPE eigen_type
!eigen-states
     COMPLEX(DCP),ALLOCATABLE :: wvf(:,:,:,:)
!eigen-values
     REAL(DP),ALLOCATABLE :: val(:,:,:)
     REAL(DP),ALLOCATABLE :: wvfG(:,:,:)
  ENDTYPE eigen_type
!------------------Basic data--------------------
  INTEGER(I4B) :: n1,n2,n3,n
  INTEGER(I4B) :: ng1,ng2,ng3,ng
  REAL(DP)     :: gap(3)
  REAL(DP)     :: dvol
!-------------------kgrids-----------------------
  INTEGER(I4B) :: nk1,nk2,nk3,nk
  REAL(DP)     :: kdispl(3)
!basic grid mesh
  TYPE(grid_type)  ::  grid
!band structrue kpoints
  TYPE(kgrid_type) ::  KPT
!eigen states and value(k-represent)
  TYPE(eigen_type) :: eigen
!

  INTEGER(I4B) :: global_n1,global_n2,global_n3,global_n
  INTERFACE fft_sph
     MODULE PROCEDURE fft_sph_r2c
     MODULE PROCEDURE fft_sph_c2r
  END INTERFACE fft_sph

CONTAINS
!-------------------PARTING LINE--------------------
  SUBROUTINE Build_rgrid()
!USE math , ONLY : lat2matrix
     USE parameters , ONLY : NSPIN,init_gap,gridn,Lcore_val
     USE struct_module , ONLY : lat_mat,lat_para  &
                       &,  recip_lat ,reclat_para &
                       &, ncharge, volume,volsp,charge_ave
     IMPLICIT NONE
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!PRINT*,'Building real space grid ...'
!CALL lat2matrix(lat_para,lat_mat,2)
!CALL lat2matrix(reclat_para,recip_lat,2)
!calculate the grid number
     IF(gridn(1)>0.AND.gridn(2)>0.AND.gridn(3)>0)THEN
        n1=gridn(1)
        n2=gridn(2)
        n3=gridn(3)
     ELSE
        n1=NINT(lat_para(1)/init_gap)
        n2=NINT(lat_para(2)/init_gap)
        n3=NINT(lat_para(3)/init_gap)
     ENDIF
!avoiding odd grid
     n1=(n1+1)/2*2
     n2=(n2+1)/2*2
     n3=(n3+1)/2*2
!grid size
     gap(1)=lat_para(1)/n1
     gap(2)=lat_para(2)/n2
     gap(3)=lat_para(3)/n3
!total mesh points in real space
     n=n1*n2*n3
!
     ng1=n1/2+1
     ng2=n2
     ng3=n3
     ng=ng1*ng2*ng3
!
!volume of every grid point
     dvol=volume/n
     CALL destroy_rgrid()

!CALL build_parallel_3d_grid()
     CALL build_parallel_sph_grid()

     volsp=global_n*dvol
!charge_ave=ne/(global_n*dvol)
     charge_ave=ncharge/volsp
!3D
     ALLOCATE(grid%rhoS(n,nspin))
     ALLOCATE(grid%rho(n))
     ALLOCATE(grid%vxcS(n,nspin))
     ALLOCATE(grid%vhxcd(n,nspin))
     ALLOCATE(grid%vlpp(n))
     ALLOCATE(grid%vh(n))
!partial core charge
     IF(Lcore_val)THEN
        ALLOCATE(grid%rhoc(n))
     ENDIF
!recip grids
     ALLOCATE(grid%gVec(4,ng))
     ALLOCATE(grid%gMask(ng))
!r-space
     ALLOCATE(grid%rVec(4,n))

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE Build_rgrid
!-------------------PARTING LINE--------------------
  SUBROUTINE destroy_rgrid()
     IMPLICIT NONE
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!3D
     IF(ALLOCATED(grid%rhoS))   DEALLOCATE(grid%rhoS)
     IF(ALLOCATED(grid%rho))   DEALLOCATE(grid%rho)
     IF(ALLOCATED(grid%vxcS))   DEALLOCATE(grid%vxcS)
     IF(ALLOCATED(grid%vhxcd))  DEALLOCATE(grid%vhxcd)
     IF(ALLOCATED(grid%vlpp))   DEALLOCATE(grid%vlpp)
     IF(ALLOCATED(grid%vh))  DEALLOCATE(grid%vh)
     IF(ALLOCATED(grid%rhoc))   DEALLOCATE(grid%rhoc)
!recip
     IF(ALLOCATED(grid%gVec))   DEALLOCATE(grid%gVec)
     IF(ALLOCATED(grid%gMask))  DEALLOCATE(grid%gMask)
!r-space , k-mesh
     IF(ALLOCATED(grid%rVec))   DEALLOCATE(grid%rVec)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE destroy_rgrid
!-------------------PARTING LINE--------------------
  SUBROUTINE Build_kgrid()
     USE parameters , ONLY : NSPIN,KSPACING,kgrid,IGamma,Isym
     USE struct_module , ONLY : recip_lat,lat_para
     IMPLICIT NONE
!
     INTEGER(I4B) :: k1,k2,k3,I,IERR,Itmp
     REAL(DP) :: wk0,kvec(3)
     LOGICAL  :: Lexist,lkmesh
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!PRINT*,'Building k-points grid ...'
     CALL destroy_KPT()
     lkmesh=((kgrid(1)>0).AND.(kgrid(2)>0).AND.(kgrid(3)>0))
!case for kspacing work
     IF(KSPACING>0.d0.OR.lkmesh)THEN
!symmetry operator
        CALL Symm_kgrid()
     ELSE
!only gamma point
        IGamma=1
        nk1=1
        nk2=1
        nk3=1
        nk=1
!>>>>>>>>>>>>>>>>>>>>>>>>>>
        ALLOCATE(KPT%vec(3,nk))
        ALLOCATE(KPT%vcar(3,nk))
        ALLOCATE(KPT%wk(nk))
!<<<<<<<<<<<<<<<<<<<<<<<<<<
        KPT%wk(:)=1
        KPT%vec(:,:)=0.d0
        KPT%vcar(:,:)=0.d0
     ENDIF
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE Build_kgrid
!-------------------PARTING LINE--------------------
  SUBROUTINE destroy_KPT()
     IMPLICIT NONE
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     IF(ALLOCATED(KPT%wk))    DEALLOCATE(KPT%wk)
     IF(ALLOCATED(KPT%vec))   DEALLOCATE(KPT%vec)
     IF(ALLOCATED(KPT%vcar))  DEALLOCATE(KPT%vcar)
!IF(ALLOCATED(KPT%nksym))  DEALLOCATE(KPT%nksym)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE destroy_KPT
!-------------------PARTING LINE--------------------
  SUBROUTINE Build_eigen()
     USE parameters , ONLY : Nspin,Nstates &
                        &,   ISTART,IGamma,Nstates_global
     IMPLICIT NONE
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!PRINT*,'Building eigen data structural ...'
     CALL destroy_eigen()
     IF(IGamma<=0)THEN
        ALLOCATE(eigen%wvf(n,Nstates,nk,Nspin))
     ELSE !have gamma point
        IF(nk>=2)THEN
           ALLOCATE(eigen%wvf(n,Nstates,nk-1,Nspin))
        ENDIF
!
        ALLOCATE(eigen%wvfG(n,Nstates,Nspin))
     ENDIF
!eigen value

     ALLOCATE(eigen%val(Nstates_global,nk,Nspin))

!PRINT*,'Building eigen data structural done'
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE Build_eigen
!-------------------PARTING LINE--------------------
  SUBROUTINE destroy_eigen()
     IMPLICIT NONE
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!k
     IF(ALLOCATED(eigen%wvf)) DEALLOCATE(eigen%wvf)
     IF(ALLOCATED(eigen%val)) DEALLOCATE(eigen%val)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE destroy_eigen
!---------------------------------------------------
  SUBROUTINE FillQTable()
!##########################################################
!* CREATED_TIME  : 2015-05-08
!* AUTHOR        : Yuan Zhou
!* CHANGE        : Xuecheng Shao
!* ADD           :
!* DESCRIPTION   :
!     This routine updates the table of norms that tells energy calculators in
!     reciprocal space whether to include a given q-point in the sum or leave it
!     out (based on whether it is within the energyCutoff sphere in recip. sp.)
!     It should be called every time the cell dimensions are altered.
!* REFERENCES    :
!     ------
!* LOG           :
!     2015-05-08 :
!* LAST MODIFIED : 2015-05-08 08:04:24 PM
!##########################################################
     USE struct_module , ONLY : recip_lat
     USE math , ONLY : Norm
     IMPLICIT NONE
!
     INTEGER(I4B) :: ix,iy,iz,Ig
     REAL(DP)     :: mVector(3),qPoint(3)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      Ig=0
      DO iz = 1, ng3
         mVector(3) = iz - 1
! Wrap around if index is greater than 1/2 the table size
         IF (mVector(3)>ng3/2) mVector(3)=mVector(3)-ng3
      DO iy = 1, ng2
         mVector(2) = iy - 1
! Wrap around if index is greater than 1/2 the table size
         IF (mVector(2) > ng2/2) mVector(2)=mVector(2)-ng2
      DO ix = 1, ng1
         mVector(1) = ix - 1
!
         Ig=Ig+1
! Calculate the qPoint cartesian coordinates given by this mVector.
!qPoint = MATMUL(mVector,recip_lat)
         qPoint = MATMUL(recip_lat,mVector)
! Assign this point in the qTable and qVectors.
         grid%gVec(1:3,Ig) = qPoint(:)
         grid%gVec(4,Ig) = Norm(qPoint)
         IF(((ix >= 2) &
            .OR. (ix == 1 .AND. iy >= 2 .AND. iy <= ng2/2+1) &
            .OR. (ix == 1 .AND. iy == 1 &
            .AND. iz >= 2 .AND. iz <= ng3/2+1))) THEN
            grid%gMask(Ig) = .TRUE.
         ELSE
            grid%gMask(Ig) = .FALSE.
         END IF
      ENDDO
      ENDDO
      ENDDO
!open(unit=121, file="qVectors_omp.dat",status="replace", access="sequential",action="readwrite")
!open(unit=121, file="qVectors.dat",status="replace", access="sequential",action="readwrite")
!write(121,*) qVectors
!close(121)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE FillQTable
!------------------Fill r_Table---------------------
  SUBROUTINE FillRTable()
     USE struct_module , ONLY : lat_mat,recip_lat
     USE math , ONLY : Norm

     USE smpi_math_module

     IMPLICIT NONE
!
     INTEGER(I4B) :: ix,iy,iz,I
     REAL(DP)     :: mVector(3),rPoint(3)
     REAL(DP)     :: slat(3,3)
     INTEGER(I4B) :: shiftn
!Initial


     shiftn=parallel%mygrid_range(1)-1
!Initial
     slat(:,1)=lat_mat(:,1)/global_n1
     slat(:,2)=lat_mat(:,2)/global_n2
     slat(:,3)=lat_mat(:,3)/global_n3
!r mesh
     DO I = 1,n
!get id
        ix=diff_map%local_map(1,shiftn+I)
        iy=diff_map%local_map(2,shiftn+I)
        iz=diff_map%local_map(3,shiftn+I)
        mVector(3) = iz - 1
        mVector(2) = iy - 1
        mVector(1) = ix - 1
!
        rPoint = MATMUL(slat,mVector)
! Assign this point in the qTable and qVectors.
        grid%rVec(1:3,I) = rPoint(:)
        grid%rVec(4,I) = Norm(rPoint)
     ENDDO

  ENDSUBROUTINE FillRTable
!-----------------------symm-kgrid------------------
  SUBROUTINE Symm_kgrid()
     USE math , ONLY : thr2mat,norm
     USE parameters , ONLY : KSPACING,Isym,kgrid,IGamma
     USE struct_module , ONLY : Opsym,nsym,recip_lat,lat_para,reclat_para
     IMPLICIT NONE
     INTEGER(I4B) :: nkr  !total k-points
     REAL(DP),ALLOCATABLE :: xk(:,:),wk(:),wk0(:),kvec(:,:)
     INTEGER(I4B),ALLOCATABLE :: equiv(:)
     REAL(DP) :: xkr(3),xx,yy,zz,fact
!
     INTEGER(I4B) :: k1,k2,k3,Ik,Is,I,J,K,Il,offsetk(3)
     INTEGER(I4B) :: nk0
     LOGICAL :: time_reversal,linlist
     REAL(DP) :: tmpwk
     REAL(DP),PARAMETER :: eps=1D-5
     LOGICAL :: LGamma=.FALSE.
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     IF(Isym>=0)THEN
!print*,'11111'
        time_reversal=.TRUE.
     ELSE
!print*,'11112'
        time_reversal=.FALSE.
     ENDIF
!kgrid
     IF(kgrid(1)>0.AND.kgrid(2)>0.AND.kgrid(3)>0)THEN
       nk1=kgrid(1)
       nk2=kgrid(2)
       nk3=kgrid(3)
     ELSE
!nk1=MAX(1,CEILING(2*pi/lat_para(1)/KSPACING))
!nk2=MAX(1,CEILING(2*pi/lat_para(2)/KSPACING))
!nk3=MAX(1,CEILING(2*pi/lat_para(3)/KSPACING))
       nk1=MAX(1,CEILING(reclat_para(1)/KSPACING))
       nk2=MAX(1,CEILING(reclat_para(2)/KSPACING))
       nk3=MAX(1,CEILING(reclat_para(3)/KSPACING))
     ENDIF
     nkr=nk1*nk2*nk3
     ALLOCATE(xk(3,nkr),kvec(3,nkr),wk(nkr),wk0(nkr))
     ALLOCATE(equiv(nkr))
!generate the k-points by Monkhorst and Pack sampling
     Ik=0
     DO k3=1,nk3
     DO k2=1,nk2
     DO k1=1,nk1
        Ik=Ik+1 
        xk(1,Ik)=(2.0_DP*k1-nk1-1)/2.0_DP/nk1
        xk(2,Ik)=(2.0_DP*k2-nk2-1)/2.0_DP/nk2
        xk(3,Ik)=(2.0_DP*k3-nk3-1)/2.0_DP/nk3
     ENDDO
     ENDDO
     ENDDO
     offsetk(1)=1-nk1
     offsetk(2)=1-nk2
     offsetk(3)=1-nk3
!find the Irreducible Brillouin Zon(IBZ) k-points
     DO Ik=1,nkr
        equiv(Ik)=Ik
     ENDDO
!
     DO Ik=1,nkr
!  check if this k-point has already been found equivalent to another
        IF(equiv(Ik)==Ik)THEN
           wk(Ik)=1.0_DP
!  check if there are equivalent k-point to this in the list
!  (excepted those previously found to be equivalent to another)
!  check both k and -k
           DO Is=1,nsym
              xkr(:)=MATMUL(TRANSPOSE(Opsym(:,:,Is)),xk(:,Ik))
              xkr(:)=xkr(:)-NINT(xkr(:))
!IF(Odet(Ik)==-1)
              xx=xkr(1)*nk1-0.5_DP*offsetk(1)
              yy=xkr(2)*nk2-0.5_DP*offsetk(2)
              zz=xkr(3)*nk3-0.5_DP*offsetk(3)
              linlist= ABS(xx-NINT(xx))<=eps .AND. &
                    &  ABS(yy-NINT(yy))<=eps .AND. &
                    &  ABS(zz-NINT(zz))<=eps
              IF(linlist)THEN
                 I=MOD(NINT(xkr(1)*nk1-0.5_DP*offsetk(1)+2*nk1),nk1)+1
                 J=MOD(NINT(xkr(2)*nk2-0.5_DP*offsetk(2)+2*nk2),nk2)+1
                 K=MOD(NINT(xkr(3)*nk3-0.5_DP*offsetk(3)+2*nk3),nk3)+1
!find the index
                 CALL thr2mat(nk1,nk2,nk3,I,J,K,Il)
!havn't been used?
                 IF(Il>Ik.AND.equiv(Il)==Il.AND.(Isym>0))THEN
                    equiv(Il)=Ik
                    wk(Ik)=wk(Ik)+1._DP
                 ELSE

                    IF(equiv(Il)/=Ik.OR.Il<Ik)THEN
                       WRITE(*,*) 'Symmetry k-points: Some thing wrong 1'
                       STOP
                    ENDIF

                 ENDIF

              ENDIF
!time reversal
              IF(time_reversal)THEN
                 xx=-xkr(1)*nk1-0.5_DP*offsetk(1)
                 yy=-xkr(2)*nk2-0.5_DP*offsetk(2)
                 zz=-xkr(3)*nk3-0.5_DP*offsetk(3)
                 linlist=ABS(xx-NINT(xx))<=eps .AND. &
                     &   ABS(yy-NINT(yy))<=eps .AND. &
                     &   ABS(zz-NINT(zz))<=eps
                 IF(linlist)THEN
                    I=MOD(NINT(-xkr(1)*nk1-0.5_DP*offsetk(1)+2*nk1),nk1)+1
                    J=MOD(NINT(-xkr(2)*nk2-0.5_DP*offsetk(2)+2*nk2),nk2)+1
                    K=MOD(NINT(-xkr(3)*nk3-0.5_DP*offsetk(3)+2*nk3),nk3)+1
                    CALL thr2mat(nk1,nk2,nk3,I,J,K,Il)
                    IF(Il>Ik.AND.equiv(Il)==Il)THEN
                       equiv(Il)=Ik
                       wk(Ik)=wk(Ik)+1._DP
                    ELSE

                       IF(equiv(Il)/=Ik.OR.Il<Ik)THEN
                          WRITE(*,*) 'Symmetry k-points: Some thing wrong 2'
                          STOP
                       ENDIF

                    ENDIF

                 ENDIF

              ENDIF
!all symmetry operator
           ENDDO

        ENDIF
!cycle all k
     ENDDO
!
     nk0=0
     fact=0._DP
     wk0(:)=0._DP
     kvec(:,:)=0._DP
     DO Ik=1,nkr
        IF(equiv(Ik)==Ik)THEN
          nk0=nk0+1
          IF(nk0>nkr) STOP 'Symmetry k-points: Some thing wrong 3'
          wk0(nk0)=wk(Ik)
          fact=fact+wk0(nk0)
!store kvec
          kvec(:,nk0)=xk(:,Ik)
        ENDIF
     ENDDO
!normalize
     wk0(1:nk0)=wk0(1:nk0)/fact
!build data
     nk=nk0
!-----------------------
     ALLOCATE(KPT%vec(3,nk))
     ALLOCATE(KPT%vcar(3,nk))
     ALLOCATE(KPT%wk(nk))
!ALLOCATE(KPT%nksym(nk))
!-----------------------
     KPT%vec(:,1:nk)=kvec(:,1:nk)
     KPT%wk(1:nk)=wk0(1:nk)
!print*,'kmesh used',sum(KPT%wk(:)),fact
     DO Ik=1,nk

        KPT%vcar(:,Ik)=MATMUL(recip_lat,KPT%vec(:,Ik))
!print*,KPT%vec(:,Ik),KPT%wk(Ik)
        IF(norm(KPT%vcar(:,Ik))<xtiny)THEN
           LGamma=.TRUE.
           IGamma=Ik
        ENDIF

     ENDDO
!set Gamma Point in end
     IF(LGamma.AND.nk>1)THEN
!set IGamma to 1
        KPT%vec(:,IGamma)=KPT%vec(:,nk)
        KPT%vcar(:,IGamma)=KPT%vcar(:,nk)
        tmpwk=KPT%wk(IGamma)
        KPT%wk(IGamma)=KPT%wk(nk)
!replace
        KPT%vec(:,nk)=0._DP
        KPT%vcar(:,nk)=0._DP
        KPT%wk(nk)=tmpwk
        IGamma=nk
     ENDIF
!deallocate something
     DEALLOCATE(xk,kvec,wk,wk0,equiv)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE Symm_kgrid
!--------------------symm-density-------------------
  SUBROUTINE Symm_density(rho)
!USE parameters , ONLY : nspin
     USE struct_module , ONLY : nsym !,Opsym,Otrans,lat_mat
     IMPLICIT NONE
!IN/OUT
     REAL(DP),INTENT(INOUT) :: rho(n1,n2,n3)
!LOCAL
     INTEGER(I4B) :: Ix,Iy,Iz,Ip,Isy,Idir
!REAL(DP) :: rt(3),offset(3) !temp veriables
     INTEGER(I4B) :: Ixyz(3),srxyz(3)  &
                  & , sizen(3) , indexs(3,nsym) ! , Nop
     REAL(DP) :: weight,rhot !weight
!REAL(DP) :: rhoin(n1,n2,n3) !,trans(3)
     LOGICAL  :: lflag(n1,n2,n3),ltmp ! ,lsym(nsym)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     IF(nsym<=1) RETURN
!start symmetry the density
!store rhoin
!rhoin(:,:,:)=rho(:,:,:)
!store
     sizen(1)=n1
     sizen(2)=n2
     sizen(3)=n3
!initialize
     lflag(:,:,:)=.FALSE.
!Nop=0
!DO Isy=1,nsym
!   CALL Isymm_apply(Isy,sizen,sizen,srxyz,lsym(Isy))
!   IF(lsym(Isy))THEN
!      Nop=Nop+1
!   ENDIF
!ENDDO
!print*,'finally find symm:',Nop
!set weight
     weight=1._DP/nsym
!all real space points
     DO Iz=1,n3
     DO Iy=1,n2
     DO Ix=1,n1
        IF(lflag(Ix,Iy,Iz)) CYCLE
!set integer mesh coordinates and
!move to center of cell in real coordinates
        Ixyz(1)=Ix
        Ixyz(2)=Iy
        Ixyz(3)=Iz
!set rhot=0
        rhot=0._DP
!cycle all symmetry operators
        DO Isy=1,nsym
!operate on vector
           CALL Isymm_apply(Isy,sizen,Ixyz,srxyz,ltmp)
!store index
           indexs(:,Isy)=srxyz(:)
           rhot=rhot+rho(srxyz(1),srxyz(2),srxyz(3))
           IF(.NOT.(rhot<1000.d0.AND.rhot>0.d0))THEN
               print*,rhot
               print*,srxyz
               STOP
           ENDIF
!print*,'xyz',Ixyz
!print*,'newxyz',srxyz
!pause
        ENDDO
!apply weight
        rhot=weight*rhot
        DO Isy=1,nsym
           rho(indexs(1,Isy),indexs(2,Isy),indexs(3,Isy)) &
          &  = rhot
           lflag(indexs(1,Isy),indexs(2,Isy),indexs(3,Isy))=.TRUE.
        ENDDO
!        print*,'origin density',rhoin(Ix,Iy,Iz)
!        print*,'density',rho(Ix,Iy,Iz)
!pause
     ENDDO
     ENDDO
     ENDDO
!test
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE Symm_density
!---------------------------------------------------
  SUBROUTINE Isymm_apply(Isy,nsize,vin,vout,lsymm)
     USE math , ONLY : inv_33
     USE struct_module , ONLY : nsym,Opsym,Otrans,lat_mat
     IMPLICIT NONE
!IN/OUT
     INTEGER(I4B),INTENT(IN)  :: Isy,nsize(3),vin(3)
     INTEGER(I4B),INTENT(OUT) :: vout(3) !out point
     LOGICAL,INTENT(OUT) :: lsymm
!LOCAL
     REAL(DP) :: vec(3),svec(3) !temp
     REAL(DP) :: ntrans(3),Op(3,3),offset(3) !operator
     LOGICAL :: linlist !list of point
     REAL(DP),PARAMETER :: eps=1D-5 !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!set vec(in)
     offset(:)=1._DP
     vec(:)=(vin(:)-offset(:))/nsize(:)
!point operater
!vout=NINT(MATMUL(TRANSPOSE(Opsym(:,:,Isym)),vin(:)))
!svec(:)=MATMUL(TRANSPOSE(Opsym(:,:,Isym)),vec(:))+Otrans(:,Isym)
     svec(:)=MATMUL(Opsym(:,:,Isy),vec(:))+Otrans(:,Isy)
!translation
!ntrans(1)=NINT(Otrans(1,Isym)*n1)
!ntrans(2)=NINT(Otrans(2,Isym)*n2)
!ntrans(3)=NINT(Otrans(3,Isym)*n3)
!move back
!print*,Opsym(:,:,Isym)
!print*,'vecin',vec
!print*,'svec',svec
!pause
!move back
     svec(:)=svec(:)*nsize(:)+offset(:)
!print*,'vecin',vin
!print*,'svec',svec
!Does it in list?
     linlist=ABS(svec(1)-NINT(svec(1)))<=eps .AND. &
           & ABS(svec(2)-NINT(svec(2)))<=eps .AND. &
           & ABS(svec(3)-NINT(svec(3)))<=eps
     IF(linlist)THEN
!find the grid
       lsymm=.TRUE.
!output new mesh
       vout(:) = MOD( NINT(svec(:)+10*nsize(:)) , nsize(:) )
       IF(vout(1)==0) vout(1)=n1 
       IF(vout(2)==0) vout(2)=n2
       IF(vout(3)==0) vout(3)=n3
!test
       IF(vout(1)<0.OR.vout(2)<0.OR.vout(3)<0)THEN
          WRITE(*,*) 'Bug in symmetry density'
          STOP
       ENDIF
     ELSE
!didn't find the mesh
       lsymm=.FALSE.
!output origin mesh
       vout(:)=vin(:)
     ENDIF
!print*,'vout',vout
!pause
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE Isymm_apply
!---------------------------------------------------
  SUBROUTINE build_parallel_3d_grid()
    USE parameters, ONLY: NSPIN,finite_order
    USE smpi_math_module, ONLY: grid_split,parallel,diff_map,destroy_diff_map
    implicit none
    INTEGER(I4B) :: i,j,k

    call destroy_diff_map()
    allocate(diff_map%nz_map(2,n3))
!> border in z
    do k=1,n3,1
       diff_map%nz_map(1,k)=(k-1)*n1*n2+1
       diff_map%nz_map(2,k)=k*n1*n2
    enddo

!> divide n grid in dims(1) cores
    global_n1=n1
    global_n2=n2
    global_n3=n3
    global_n=n
    CALL grid_split(global_n, parallel%dims(2),parallel%commx,parallel%rankx&
         & , parallel%mygrid_range, parallel%recvcounts&
         & , parallel%displs, parallel%global_gridrange&
         & , n1, n2, n3, n)
  ENDSUBROUTINE build_parallel_3d_grid
!---------------------------------------------------

  SUBROUTINE build_parallel_sph_grid()
    USE parameters, ONLY: NSPIN,finite_order
    USE smpi_math_module, ONLY: grid_split,parallel,destroy_diff_map
    implicit none
    INTEGER(I4B) :: i,j,k

!> divide n grid in dims(1) cores
    global_n1=n1
    global_n2=n2
    global_n3=n3
    global_n=n
    CALL grid_split_sph(global_n, parallel%dims(2),parallel%commx,parallel%rankx&
         & , parallel%mygrid_range, parallel%recvcounts&
         & , parallel%displs, parallel%global_gridrange&
         & , n1, n2, n3, n)
  ENDSUBROUTINE build_parallel_sph_grid
!---------------------------------------------------
  SUBROUTINE grid_split_sph(ngrid,ncore,comm,id,grid_range,recvcounts,displs,gridrange_sum,n1,n2,n3,n)
    USE parameters, ONLY: atomRc
    USE struct_module, ONLY:lat_mat,struct,naty
    USE math,  ONLY:norm
    USE smpi_math_module, ONLY: MPI_INTEGER4,mpinfo,diff_map,destroy_diff_map
    implicit none
!> in/out
    INTEGER(I4B),INTENT(IN)  :: ngrid,ncore,comm,id
    INTEGER(I4B),INTENT(OUT) :: grid_range(3)
    INTEGER(I4B),INTENT(OUT) :: recvcounts(ncore)
    INTEGER(I4B),INTENT(OUT) :: displs(ncore)
    INTEGER(I4B),INTENT(OUT) :: gridrange_sum(3,ncore)
    INTEGER(I4B),INTENT(INOUT) :: n1,n2,n3,n
!> local
    REAL(DP) :: slat(3,3),mVector(3),rpoint(3),temp_r
    INTEGER(I4B) :: i,j,k,In,Ip,Ity,Ia
    INTEGER(I4B) :: average,redund
    INTEGER(I4B) :: recvcounts_temp(ncore)
    INTEGER(I4B) :: displs_temp(ncore)
    INTEGER(I4B) :: nrepl=1,icx,icy,icz

!> border in z
    call destroy_diff_map()
    allocate(diff_map%nz_map(2,n3))

    slat(:,1)=lat_mat(:,1)/n1
    slat(:,2)=lat_mat(:,2)/n2
    slat(:,3)=lat_mat(:,3)/n3

!r mesh
    In=0
    if(.not.allocated(grid%lsp))allocate(grid%lsp(n1,n2,n3))
    grid%lsp=.false.

    if(atomRc>0.d0)then
       CALL sphere_region(n1,n2,n3,grid%lsp,diff_map%nz_map,in)
    else
       DO k = 1,n3
          mVector(3) = k - 1
          diff_map%nz_map(1,k)=In+1 !> for parallel
       DO j = 1,n2
          mVector(2) = j - 1
       DO i = 1,n1
          mVector(1) = i - 1
          rPoint = MATMUL(slat,mVector)

          in=in+1
          grid%lsp(i,j,k)=.true.
       ENDDO
       ENDDO
       diff_map%nz_map(2,k)=In
       ENDDO
    endif

    grid%Onedlength=In
    Ip=0; In=0
    allocate(diff_map%local_map(3,grid%Onedlength))
    do k=1,n3,1
    do j=1,n2,1
    do i=1,n1,1
       Ip=Ip+1
       if(grid%lsp(i,j,k))then
          in=in+1
          diff_map%local_map(:,in)=(/i,j,k/)
       endif
    enddo
    enddo
    enddo


!> split the grid into "ncore" cores
    average=grid%oneDlength/ncore
    redund=mod(in,ncore)
    if(id.lt.redund)then
!> id*(l+1) + 1, ..., (id*(l+1)+1) + (l+1) - 1
       grid_range(1)=id*(average+1)+1
       grid_range(2)=(id+1)*(average+1)
       grid_range(3)=average+1
    else
!> redund*(l+1)+(id-redund)*l+1, ..., (redund*(l+1)+(id-redund)*l+1) + l - 1
       grid_range(1)=redund+id*average+1
       grid_range(2)=redund+(id+1)*average
       grid_range(3)=average
    endif

!> set the recvcounts and displs
    do i=0,ncore-1,1
       j=i+1
       if(i.lt.redund)then
!> id*(l+1) + 1, ..., (id*(l+1)+1) + (l+1) - 1
          recvcounts(j)=average+1
          displs(j)=i*(average+1)+1-1
       else
!> redund*(l+1)+(id-redund)*l+1, ..., (redund*(l+1)+(id-redund)*l+1) + l - 1
          recvcounts(j)=average
          displs(j)=redund+i*average+1-1
       endif
    enddo

!> set the global_gridrange
    recvcounts_temp=3
    displs_temp=(/(i*3,i=0,ncore-1,1)/)
    CALL MPI_ALLGATHERV(grid_range,3,MPI_INTEGER4,gridrange_sum,recvcounts_temp,&
         displs_temp,MPI_INTEGER4,comm,mpinfo)

!> set the local cubic grid
! n3=(grid_range(2)-1)/n1/n2+1-(grid_range(1)-1)/n1/n2-1+1
    n3=diff_map%local_map(3,grid_range(2))&
         & -diff_map%local_map(3,grid_range(1))
    n=grid_range(3)
    global_n=grid%oneDlength
  ENDSUBROUTINE grid_split_sph

!---------------------------------------------------

  SUBROUTINE sphere_region(n1,n2,n3,Lsphere,nz_map,nsphere)
    USE parameters, ONLY: atomRc
    USE struct_module, ONLY:lat_mat,struct,naty,natom
    USE math,  ONLY:norm,cross,inv_33
    implicit none
!> in/out
    INTEGER(I4B) :: n1,n2,n3
    LOGICAL :: Lsphere(n1,n2,n3)
    INTEGER(I4B) :: nz_map(2,n3)
    INTEGER(I4B) :: nsphere
!> local
    INTEGER(I4B) :: rb2,nb(3),nbh(3),Ia,Id(3,2)
    INTEGER(I4B) :: Ix,Iy,Iz,In,ngrid_z(n3)
    INTEGER(I4B) :: mVector(3),i,j,k
    REAL(DP) :: slat(3,3),r(3),rpoint(3),dist

!> init
    Lsphere=.false.
    ngrid_z=0

    slat(:,1)=lat_mat(:,1)/n1
    slat(:,2)=lat_mat(:,2)/n2
    slat(:,3)=lat_mat(:,3)/n3

!> inner cur ball radius
    dist=dot_product(cross(slat(:,1),slat(:,2)),slat(:,3))/&
         &(norm(cross(slat(:,1),slat(:,2)))+&
         & norm(cross(slat(:,2),slat(:,3)))+&
         & norm(cross(slat(:,3),slat(:,1))))*3

    rb2 = 2.0*atomRc !box length
    nb(1) = min(INT(rb2*norm(slat(:,1))/dist/gap(1))+2,n1)
    nb(2) = min(INT(rb2*norm(slat(:,2))/dist/gap(2))+2,n2)
    nb(3) = min(INT(rb2*norm(slat(:,3))/dist/gap(3))+2,n3)
!xq avoid odd grids
    do i =1,3
       IF(MOD(nb(i),2)/=0)THEN
          nb(i)=nb(i)+1
       ENDIF
       nbh(i)=nb(i)/2
    enddo

!mesh index
    DO Ia=1,natom

!1) find atomic begin and end point
       r(:)=struct%poscar(:,Ia)
       rpoint=MATMUL(inv_33(slat),r)
       ix=INT(rpoint(1))+1
       iy=INT(rpoint(2))+1
       iz=INT(rpoint(3))+1
!x
       id(1,1)=ix-nbh(1)+1
       id(1,2)=ix+nbh(1)
!
       id(2,1)=iy-nbh(2)+1
       id(2,2)=iy+nbh(2)
!z
       id(3,1)=iz-nbh(3)+1
       id(3,2)=iz+nbh(3)

!2) judge in a small cubic grid
       do Iz=id(3,1),id(3,2),1
       do Iy=id(2,1),id(2,2),1
       do Ix=id(1,1),id(1,2),1

          mVector(3) = Iz - 1
          mVector(2) = Iy - 1
          mVector(1) = Ix - 1
          rPoint = MATMUL(slat,mVector)
          dist=norm(rpoint-r)

          if( dist < atomRc)then
!> set periodic boundary
             i=mod(mod(Ix-1,n1)+n1,n1)+1
             j=mod(mod(Iy-1,n2)+n2,n2)+1
             k=mod(mod(Iz-1,n3)+n3,n3)+1

             if(Lsphere(i,j,k))then
                cycle
             else
                Lsphere(i,j,k)=.true.
                nsphere=nsphere+1
                ngrid_z(k)=ngrid_z(k)+1
             endif
          endif

       enddo !> Ix
       enddo !> Iy
       enddo !> Iz

     ENDDO

!> assign nz_map
     in=0
     do i = 1,n3,1
        nz_map(1,i)=in+1
        in=in+ngrid_z(i)
        nz_map(2,i)=in
     enddo
     if(in/=nsphere)then
        stop "error: n in sphere inconsistence nz_map "
     endif

  ENDSUBROUTINE sphere_region

!---------------------------------------------------
  SUBROUTINE Sphere2Cubic(nps,f1d,f3d,rfill)

     USE smpi_math_module

     IMPLICIT NONE
     INTEGER(I4B),INTENT(IN) :: nps
     REAL(DP),INTENT(IN) :: f1d(nps)
     REAL(DP),OPTIONAL   :: rfill
     REAL(DP),INTENT(OUT) :: f3d(:,:,:)
     INTEGER(I4B) :: ix,iy,iz,id,nx,ny,nz,nn

     REAL(DP), ALLOCATABLE :: f3d_local(:,:,:)
     INTEGER(I4B) :: shiftn

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     nx=SIZE(f3d,1)
     ny=SIZE(f3d,2)
     nz=SIZE(f3d,3)
     nn=nx*ny*nz


     ALLOCATE(f3d_local(nx,ny,nz))
     f3d_local=0._DP
     shiftn=parallel%mygrid_range(1)-1


     id=0
     IF(PRESENT(rfill))THEN
        f3d(:,:,:)=rfill
     ELSE
        f3d(:,:,:)=0._DP
     ENDIF
!set values
     DO id=1,nps
!get id
        ix=diff_map%local_map(1,shiftn+id)
        iy=diff_map%local_map(2,shiftn+id)
        iz=diff_map%local_map(3,shiftn+id)
!get values

        f3d_local(ix,iy,iz)=f1d(id)

     ENDDO

!communcations

     CALL MPI_ALLREDUCE(f3d_local,f3d,nn,MPI_REAL8,MPI_SUM,parallel%commx,mpinfo)
     DEALLOCATE(f3d_local)

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE Sphere2Cubic
!---------------------------------------------------
  SUBROUTINE Cubic2Sphere(nps,f3d,f1d)

     USE smpi_math_module

     IMPLICIT NONE
     INTEGER(I4B),INTENT(IN) :: nps
     REAL(DP),INTENT(IN) :: f3d(:,:,:)
     REAL(DP),INTENT(OUT) :: f1d(nps)
     INTEGER(I4B) :: ix,iy,iz,id!,nx,ny,nz

     INTEGER(I4B) :: shiftn

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

     shiftn=parallel%mygrid_range(1)-1

!set values
     DO id=1,nps
!get id
        ix=diff_map%local_map(1,shiftn+id)
        iy=diff_map%local_map(2,shiftn+id)
        iz=diff_map%local_map(3,shiftn+id)
!get values
        f1d(id)=f3d(ix,iy,iz) 
     ENDDO
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE Cubic2Sphere
!---------------------------------------------------

  SUBROUTINE Cubic2Sphere_fft(nps,f3d,f1d,shiftn,shiftz)
     USE smpi_math_module
     IMPLICIT NONE
     INTEGER(I4B),INTENT(IN) :: nps
     REAL(DP),INTENT(IN) :: f3d(:,:,:)
     REAL(DP),INTENT(OUT) :: f1d(nps)
     INTEGER(I4B) :: shiftn,shiftz
     INTEGER(I4B) :: ix,iy,iz,id!,nx,ny,nz
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!set values
     DO id=1,nps
!get id
        ix=diff_map%local_map(1,shiftn+id)
        iy=diff_map%local_map(2,shiftn+id)
        iz=diff_map%local_map(3,shiftn+id)-shiftz
!get values
        f1d(id)=f3d(ix,iy,iz)
     ENDDO
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE Cubic2Sphere_fft
!---------------------------------------------------
  SUBROUTINE Sphere2Cubic_fft(nps,f1d,f3d,shiftn,shiftz)
     USE smpi_math_module
     IMPLICIT NONE
     INTEGER(I4B),INTENT(IN) :: nps
     REAL(DP),INTENT(IN) :: f1d(nps)
     REAL(DP),INTENT(OUT) :: f3d(:,:,:)
     INTEGER(I4B) :: shiftn,shiftz
     INTEGER(I4B) :: ix,iy,iz,id
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!set values
     DO id=1,nps
!get id
        ix=diff_map%local_map(1,shiftn+id)
        iy=diff_map%local_map(2,shiftn+id)
        iz=diff_map%local_map(3,shiftn+id)-shiftz
!get values
        f3d(ix,iy,iz)=f1d(id)
     ENDDO
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE Sphere2Cubic_fft
!-------------------divide line ---------------------
  FUNCTION fft_sph_r2c(array_r) result(array_c)
    USE smpi_math_module, ONLY: parallel,MPI_REAL8,mpinfo,smpi_exit
    USE FOURIER, ONLY: FFT
    implicit none
!> in/out
    REAL(DP),intent(in)  :: array_r(:)
!DIMENSION(ng1,ng2,parallel%local_z)
    COMPLEX(DCP) :: array_c(ng1,ng2,parallel%local_z)
!> local
    REAL(DP)  :: array_FFT(global_n1,global_n2,parallel%local_z)
    INTEGER(I4B) :: fft1d_size
    REAL(DP),allocatable  :: array_fft1d(:)

    fft1d_size=parallel%fft_grid_range(2,parallel%rankx+1)&
         & -parallel%fft_grid_range(1,parallel%rankx+1)+1
    allocate(array_fft1d(fft1d_size))
    CALL MPI_ALLTOALLV(array_r,parallel%fft_scount&
         & ,parallel%fft_sdispls,MPI_REAL8&
         & ,array_fft1d,parallel%fft_rcount&
         & ,parallel%fft_rdispls,MPI_REAL8&
         & ,parallel%commx,mpinfo)
!> sphere to cubic
    array_FFT=0.d0
    CALL sphere2cubic_fft(fft1d_size,array_fft1d,array_FFT&
         & ,parallel%fft_grid_range(1,parallel%rankx+1)-1&
         & ,parallel%local_z_start)

    array_c(:,:,:)=FFT(array_FFT)
    deallocate(array_fft1d)
  END FUNCTION fft_sph_r2c
!-------------------divide line ---------------------
  FUNCTION fft_sph_c2r(array_c)result(array_r)
    USE smpi_math_module, ONLY: parallel,MPI_REAL8,mpinfo,smpi_exit
    USE FOURIER, ONLY: FFT
    implicit none
!> in/out
    REAL(DP)  :: array_r(n)
!DIMENSION(ng1,ng2,parallel%local_z)
    COMPLEX(DCP),intent(in) :: array_c(:,:,:)
!> local
    REAL(DP)  :: array_FFT(global_n1,global_n2,parallel%local_z)
    INTEGER(I4B) :: fft1d_size
    REAL(DP),allocatable  :: array_fft1d(:)

    fft1d_size=parallel%fft_grid_range(2,parallel%rankx+1)&
         & -parallel%fft_grid_range(1,parallel%rankx+1)+1
    allocate(array_fft1d(fft1d_size))
    array_FFT=FFT(array_c)
    CALL cubic2sphere_fft(fft1d_size,array_FFT,array_fft1d&
         & ,parallel%fft_grid_range(1,parallel%rankx+1)-1 &
         & ,parallel%local_z_start)
    CALL MPI_ALLTOALLV(array_fft1d,parallel%fft_rcount&
         & ,parallel%fft_rdispls,MPI_REAL8&
         & ,array_r,parallel%fft_scount&
         & ,parallel%fft_sdispls,MPI_REAL8&
         & ,parallel%commx,mpinfo)
    deallocate(array_fft1d)
  END FUNCTION fft_sph_c2r
!-------------------divide line ---------------------

!---------------------------------------------------
   SUBROUTINE sumrhoS(nps,rhoS,rho)
      USE parameters , ONLY : nspin
      IMPLICIT NONE
!IN/OUT
      INTEGER(I4B)            :: nps
      REAL(DP),INTENT(IN)     :: rhoS(nps,nspin)
      REAL(DP),INTENT(OUT)    :: rho(nps)
!LOCAL
      INTEGER(I4B) :: Is
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF(nspin==1)THEN
         rho(:)=rhoS(:,1)
      ELSEIF(nspin==2)THEN
         rho(:)=rhoS(:,1)+rhoS(:,2)
      ELSE
         STOP 'sumrhoS errors'
      ENDIF
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE sumrhoS
!---------------------------------------------------
ENDMODULE grid_module

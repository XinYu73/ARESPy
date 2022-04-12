# 1 "Grid_module.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "Grid_module.f90"
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
     REAL(DP),ALLOCATABLE :: rhoS(:,:,:,:)
     REAL(DP),ALLOCATABLE :: vlpp(:,:,:)
     REAL(DP),ALLOCATABLE :: eval(:,:,:) !eigen-values
     !recip q table
     REAL(DP),ALLOCATABLE :: gVec(:,:)
     LOGICAL ,ALLOCATABLE :: gMask(:)
     !r table ktable
     REAL(DP),ALLOCATABLE :: rVec(:,:)
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
     COMPLEX(DP),ALLOCATABLE :: wvf(:,:,:,:)
     !eigen-values
     REAL(DP),ALLOCATABLE :: val(:,:,:)
  ENDTYPE eigen_type
  !real eigen_type
  TYPE eigen_type_r
     !eigen-states
     REAL(DP),ALLOCATABLE :: wvf(:,:,:,:)
     !eigen-values
     REAL(DP),ALLOCATABLE :: val(:,:,:)
  ENDTYPE eigen_type_r
  !ISO_density_sphere----xlt
# 56 "Grid_module.f90"
  TYPE charge_sphere
     INTEGER(I4B)             :: OneDLength
     REAL(DP)                 :: volume
     INTEGER(I4B),ALLOCATABLE :: x(:),y(:),z(:),n(:)
     REAL(DP),ALLOCATABLE     :: OneDSphere(:,:)
     REAL(DP),ALLOCATABLE     :: OneDVeff(:)
     REAL(DP),ALLOCATABLE     :: OneDwvf(:,:,:)
     REAL(DP),ALLOCATABLE     :: OneDval(:,:)
     REAL(DP),ALLOCATABLE     :: initMO(:,:)
     REAL(DP),ALLOCATABLE     :: vlpp(:)
  ENDTYPE charge_sphere
!   TYPE grid_diff_map_type
!      INTEGER(I4B),allocatable :: nz_map(:,:) !> the up id and down id for per nz
!      INTEGER(I4B) :: mycomm_cores(2)         !> number of cores for communcation
!      INTEGER(I4B),allocatable :: mycomm_size(:,:)
!      INTEGER(I4B),allocatable :: mysend_size(:,:)
!      INTEGER(I4B),allocatable :: local_map(:,:)
!      INTEGER(I4B)             :: boundary(2,3)
!   ENDTYPE grid_diff_map_type
  ! TYPE grid_diff_map_type
  !    INTEGER(I4B),allocatable :: local_nz(:,:)
  ! ENDTYPE grid_diff_map_type

  !------------------Basic data--------------------
  INTEGER(I4B) :: n1,n2,n3,n,nsn

  INTEGER(I4B) :: global_n1,global_n2,global_n3,global_n

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
  TYPE(eigen_type_r) :: eigen_r
  !isolation: sphere coordinates
  REAL(DP),ALLOCATABLE     :: dr(:,:,:)
  REAL(DP),ALLOCATABLE     :: cost(:,:,:)
  REAL(DP),ALLOCATABLE     :: sint(:,:,:)
  COMPLEX(DCP),ALLOCATABLE :: cns(:,:,:)
  INTEGER(I4B),ALLOCATABLE :: indx(:,:,:)
  !isolation: density sphere
  TYPE(charge_sphere)  :: rho_calc
  !> fft_sph
  INTERFACE fft_sph
     MODULE PROCEDURE fft_sph_r2c
     MODULE PROCEDURE fft_sph_c2r
  END INTERFACE fft_sph

CONTAINS
  !-------------------PARTING LINE--------------------
  SUBROUTINE Build_rgrid()
     USE math , ONLY : lat2matrix
     USE parameters , ONLY : NSPIN,init_gap,Lpbc,Lpbc2iso,isp
     USE struct_module , ONLY : lat_mat,lat_para  &
                       &,  recip_lat ,reclat_para &
                       &,  volume
     USE m_time_evaluate, only: memory_sum
     USE succeed, ONLY:Llastrho
     IMPLICIT NONE
     INTEGER(I4B),save :: n1_old,n2_old,n3_old
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     !PRINT*,'Building real space grid ...'
     CALL lat2matrix(lat_para,lat_mat,2)
     CALL lat2matrix(reclat_para,recip_lat,2)
     !calculate the grid number
     if(isp>1)then

        n1_old=global_n1
        n2_old=global_n2
        n3_old=global_n3





     endif
     n1=NINT(lat_para(1)/init_gap)
     n2=NINT(lat_para(2)/init_gap)
     n3=NINT(lat_para(3)/init_gap)
     !print*,"lat_para",lat_para(1),lat_para(2),lat_para(3),init_gap
     !=============================================================
     !isolate system and odd grid
     ! IF(.NOT.Lpbc.AND.MOD(n1,2)==0 )n1=n1+1
     ! IF(.NOT.Lpbc.AND.MOD(n2,2)==0 )n2=n2+1
     ! IF(.NOT.Lpbc.AND.MOD(n3,2)==0 )n3=n3+1
     !=============================================================
     !grid size
     gap(1)=lat_para(1)/n1
     gap(2)=lat_para(2)/n2
     gap(3)=lat_para(3)/n3

     !> restart succeed
     if((n1_old/=n1).or.(n2_old/=n2).or.(n3_old/=n3))then
        Llastrho=.false.
     endif

     !total mesh points in real space
     n=n1*n2*n3
     nsn=n*NSPIN
     !
     ng1=n1/2+1
     ng2=n2
     ng3=n3
     ng=ng1*ng2*ng3
     !
     !volume of every grid point
     dvol=volume/n
     CALL destroy_rgrid()

     CALL build_parallel_cubic_grid()

     !3D
     ALLOCATE(grid%rhoS(n1,n2,n3,NSPIN))
     ALLOCATE(grid%vlpp(n1,n2,n3))
     !recip grids
     ALLOCATE(grid%gVec(4,ng))
     ALLOCATE(grid%gMask(ng))
     !r-space
     ALLOCATE(grid%rVec(4,n))
     call memory_sum('build_grid',&
          &(real(size(grid%rhoS),DP)+size(grid%vlpp)&
          &+size(grid%gVec)+size(grid%rVec))*DP&
          &+size(grid%gMask)*kind(grid%gMask))
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE Build_rgrid
  !-----------------------PARTING-LINE----------------
  SUBROUTINE Build_rgrid_iso()
     USE math , ONLY : lat2matrix
     USE parameters , ONLY : NSPIN,init_gap,Lpbc,Lpbc2iso
     USE struct_module , ONLY : lat_mat,lat_para  &
                       &,  recip_lat ,reclat_para &
                       &,  volume
     USE m_time_evaluate, only: memory_sum
     IMPLICIT NONE
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     !PRINT*,'Building real space grid ...'
     CALL lat2matrix(lat_para,lat_mat,2)
     CALL lat2matrix(reclat_para,recip_lat,2)
     !calculate the grid number
     n1=NINT(lat_para(1)/init_gap)
     n2=NINT(lat_para(2)/init_gap)
     n3=NINT(lat_para(3)/init_gap)
     !print*,"lat_para",lat_para(1),lat_para(2),lat_para(3),init_gap
     !=============================================================
     !isolate system and odd grid
     IF(.NOT.Lpbc.AND.MOD(n1,2)==0 )n1=n1+1
     IF(.NOT.Lpbc.AND.MOD(n2,2)==0 )n2=n2+1
     IF(.NOT.Lpbc.AND.MOD(n3,2)==0 )n3=n3+1
     !=============================================================
     !grid size
     gap(1)=lat_para(1)/n1
     gap(2)=lat_para(2)/n2
     gap(3)=lat_para(3)/n3
     !total mesh points in real space
     n=n1*n2*n3
     nsn=n*NSPIN
     !
     ng1=n1/2+1
     ng2=n2
     ng3=n3
     ng=ng1*ng2*ng3
     !
     !volume of every grid point
     dvol=volume/n
     CALL destroy_rgrid()
     !3D
     ALLOCATE(grid%rhoS(n1,n2,n3,NSPIN))
     ALLOCATE(grid%vlpp(n1,n2,n3))
     !recip grids
     ALLOCATE(grid%gVec(4,ng))
     ALLOCATE(grid%gMask(ng))
     !r-space
     ALLOCATE(grid%rVec(4,n))
     call memory_sum('build_grid',&
          &(real(size(grid%rhoS),DP)+size(grid%vlpp)&
          &+size(grid%gVec)+size(grid%rVec))*DP&
          &+size(grid%gMask)*kind(grid%gMask))
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE Build_rgrid_iso
  !-----------------------PARTING-LINE----------------
  SUBROUTINE destroy_rgrid()
     USE m_time_evaluate, ONLY: memory_free
     IMPLICIT NONE
     !
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     !3D
     IF(ALLOCATED(grid%rhoS))THEN
        call memory_free('build_grid_rhos',real(size(grid%rhoS),DP)*DP)
        DEALLOCATE(grid%rhoS)
     ENDIF
     IF(ALLOCATED(grid%vlpp))THEN
        call memory_free('build_grid_vlpp',real(size(grid%vlpp),DP)*DP)
        DEALLOCATE(grid%vlpp)
     ENDIF
     !recip
     IF(ALLOCATED(grid%gVec))THEN
        call memory_free('build_grid_gVec',real(size(grid%gVec),DP)*DP)
        DEALLOCATE(grid%gVec)
     ENDIF
     IF(ALLOCATED(grid%gMask))THEN
        call memory_free('build_grid_gMask',real(size(grid%gMask),DP)&
             &*kind(grid%gMask))
        DEALLOCATE(grid%gMask)
     ENDIF
     !r-space , k-mesh
     IF(ALLOCATED(grid%rVec))THEN
        call memory_free('build_grid_rVec',real(size(grid%rVec),DP)*DP)
        DEALLOCATE(grid%rVec)
     ENDIF
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE destroy_rgrid
  !-------------------PARTING LINE--------------------
  SUBROUTINE Build_kgrid()
     USE parameters , ONLY : NSPIN,KSPACING,LKP,LBvK,Lpbc,LGamma
     USE struct_module , ONLY : recip_lat,lat_para
     USE m_time_evaluate, ONLY: memory_sum
     IMPLICIT NONE
     !
     INTEGER(I4B) :: k1,k2,k3,I,IERR
     REAL(DP) :: wk0,kvec(3)
     LOGICAL  :: Lexist
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     !PRINT*,'Building k-points grid ...'
     CALL destroy_KPT()
     IF(LKP.AND.(.NOT.LBvK).AND.Lpbc)THEN
        INQUIRE(FILE='KPOINTS',EXIST=lexist)
        !test
        IF(.NOT.lexist)THEN
           WRITE(6,*) '>>>WARNING<<<: KPOINTS.dat file is not exist.stop!!!'
           STOP
        ENDIF
        !grneralize k-points
        OPEN(106,FILE='KPOINTS',STATUS='OLD',IOSTAT=IERR)
           READ(106,*)
           READ(106,*) nk !number of kpoints
           READ(106,*)
           READ(106,*)
           !>>>>>>>>>>>>>>>>>>>>>>>>>>
           ALLOCATE(KPT%vec(3,nk))
           ALLOCATE(KPT%vcar(3,nk))
           ALLOCATE(KPT%wk(nk))
           call memory_sum('KPT',(real(size(KPT%vec),DP)+size(KPT%vcar)&
                &+size(KPT%wk))*DP)
           !<<<<<<<<<<<<<<<<<<<<<<<<<<
           DO I=1,nk
              !Read the kpoints data
              READ(106,*) kvec(1),kvec(2),kvec(3),wk0
              !store the data
              KPT%vec(:,I)=kvec(:)
              !cartesion coo
              KPT%vcar(:,I)=MATMUL(recip_lat,KPT%vec(:,I))
              !weight
              KPT%wk(I)=wk0
           ENDDO
        CLOSE(106)
     ELSE
        !case for kspacing work
        IF(KSPACING>0.d0.AND.Lpbc.AND.(.NOT.(LBvK)))THEN
           !kgrids
           nk1=MAX(2,CEILING(2*pi/lat_para(1)/KSPACING))
           nk2=MAX(2,CEILING(2*pi/lat_para(2)/KSPACING))
           nk3=MAX(2,CEILING(2*pi/lat_para(3)/KSPACING))
           !KSHIFT
           kdispl(:)=0.5d0
           !make k-grids even
           nk1=(nk1/2)*2
           nk2=(nk2/2)*2
           nk3=(nk3/2)*2
           !time inverse symmetry
           nk=nk1*nk2*nk3/2
           !>>>>>>>>>>>>>>>>>>>>>>>>>>
           ALLOCATE(KPT%vec(3,nk))
           ALLOCATE(KPT%vcar(3,nk))
           ALLOCATE(KPT%wk(nk))
           call memory_sum('KPT',(real(size(KPT%vec),DP)+size(KPT%vcar)&
                &+size(KPT%wk))*DP)
           !<<<<<<<<<<<<<<<<<<<<<<<<<<
           wk0=1.d0/nk
           KPT%wk(:)=wk0
           I=0
           DO k3=1,nk3/2
           DO k2=1,nk2
           DO k1=1,nk1
              I=I+1
              KPT%vec(1,I)=(k1-nk1/2-0.5d0)/nk1
              KPT%vec(2,I)=(k2-nk2/2-0.5d0)/nk2
              KPT%vec(3,I)=(k3-0.5d0)/nk3
              !cartentisien coo
              KPT%vcar(:,I)=MATMUL(recip_lat,KPT%vec(:,I))
           ENDDO
           ENDDO
           ENDDO
        ELSE
           !only gamma point
           nk1=1
           nk2=1
           nk3=1
           nk=1
           !>>>>>>>>>>>>>>>>>>>>>>>>>>
           ALLOCATE(KPT%vec(3,nk))
           ALLOCATE(KPT%vcar(3,nk))
           ALLOCATE(KPT%wk(nk))
           call memory_sum('KPT',(real(size(KPT%vec),DP)+size(KPT%vcar)&
                &+size(KPT%wk))*DP)
           !<<<<<<<<<<<<<<<<<<<<<<<<<<
           KPT%wk(:)=1
           KPT%vec(:,:)=0.d0
           KPT%vcar(:,:)=0.d0
        ENDIF
        !
     ENDIF
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE Build_kgrid
  !-------------------PARTING LINE--------------------
  SUBROUTINE destroy_KPT()
    USE m_time_evaluate, ONLY: memory_free
     IMPLICIT NONE
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     IF(ALLOCATED(KPT%wk))THEN
        call memory_free('KPT_wk',real(size(KPT%wk),DP)*DP)
        DEALLOCATE(KPT%wk)
     ENDIF
     IF(ALLOCATED(KPT%vec))THEN
        call memory_free('KPT_vec',real(size(KPT%vec),DP)*DP)
        DEALLOCATE(KPT%vec)
     ENDIF
     IF(ALLOCATED(KPT%vcar))THEN
        call memory_free('KPT_vcar',real(size(KPT%vcar),DP)*DP)
        DEALLOCATE(KPT%vcar)
     ENDIF
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE destroy_KPT
  !-------------------PARTING LINE--------------------
  SUBROUTINE Build_eigen()
     USE parameters , ONLY : Nspin,Nstates,Nstates_global,LBvK,NPRR,Lpbc,Lbvk,LGamma
     USE m_time_evaluate, ONLY: memory_sum
     USE struct_module, ONLY:natom
     USE succeed , ONLY : init_succeed_rho_real,init_succeed_rho_cmplx,Llastrho,destroy_succeed
     USE smpi_math_module, ONLY: parallel
     IMPLICIT NONE
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     !PRINT*,'Building eigen data structural ...'
     CALL destroy_eigen()
     IF(LBvK)THEN
        ALLOCATE(eigen_r%wvf(n,Nstates,nk,Nspin))
        call memory_sum('eigen_r_wvf',real(size(eigen_r%wvf),DP)*DP)
        IF(NPRR<=0)THEN
           ALLOCATE(eigen_r%val(Nstates,nk,Nspin))
           call memory_sum('eigen_r_val',real(size(eigen_r%val),DP)*DP)
        ENDIF
     ELSE
        IF(LGamma)THEN
           ALLOCATE(eigen_r%wvf(n,Nstates,nk,Nspin))
           ALLOCATE(eigen_r%val(Nstates_global,nk,Nspin))
           call memory_sum('eigen_r',real(size(eigen_r%wvf),DP)*DP&
                &+size(eigen_r%val)*DP)
        ELSE
           ALLOCATE(eigen%wvf(n,Nstates,nk,Nspin))
           ALLOCATE(eigen%val(Nstates_global,nk,Nspin))
           call memory_sum('eigen',real(size(eigen_r%wvf),DP)*DP&
                &+size(eigen%val)*DP)
      ENDIF
     ENDIF
     if(LGamma)then
        if(.not.Llastrho)then
           if(parallel%isroot)write(6,*)'[init succeed] new grid size is store'
           call destroy_succeed()
           CALL init_succeed_rho_real(n,natom,Nstates,Nspin,dvol)
        endif
     else
        if(.not.Llastrho)then
           if(parallel%isroot)write(6,*)'[init succeed] new grid size is store'
           call destroy_succeed()
           CALL init_succeed_rho_cmplx(n,natom,Nstates,nk,Nspin,dvol)
        endif
     endif
     !PRINT*,'Building eigen data structural done'
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE Build_eigen
  !-------------------PARTING LINE--------------------
  SUBROUTINE destroy_eigen()
     USE m_time_evaluate, ONLY: memory_free
     IMPLICIT NONE
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     !k
     IF(ALLOCATED(eigen%wvf))THEN
        call memory_free('eigen_wvf',real(size(eigen%wvf),DP)*DP)
        DEALLOCATE(eigen%wvf)
     ENDIF
     IF(ALLOCATED(eigen%val))THEN
        call memory_free('eigen_val',real(size(eigen%val),DP)*DP)
        DEALLOCATE(eigen%val)
     ENDIF
     !BvK
     IF(ALLOCATED(eigen_r%wvf))THEN
        call memory_free('eigen_r_wvf',real(size(eigen_r%wvf),DP)*DP)
        DEALLOCATE(eigen_r%wvf)
     ENDIF
     IF(ALLOCATED(eigen_r%val))THEN
        call memory_free('eigen_r_val',real(size(eigen_r%val),DP)*DP)
        DEALLOCATE(eigen_r%val)
     ENDIF
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

     USE smpi_math_module, ONLY: parallel

     IMPLICIT NONE
     !
     INTEGER(I4B) :: ix,iy,iz,I
     REAL(DP)     :: mVector(3),rPoint(3)
     REAL(DP)     :: slat(3,3)

     INTEGER(I4B) :: j,k,in,z_up
     INTEGER(I4B) :: d_z,y_down

     !Initial





     slat(:,1)=lat_mat(:,1)/global_n1
     slat(:,2)=lat_mat(:,2)/global_n2
     slat(:,3)=lat_mat(:,3)/global_n3

     !r mesh

     ix=mod(mod(parallel%mygrid_range(1)-1,global_n1*global_n2),n1)+1
     iy=mod((parallel%mygrid_range(1)-1)/global_n1,global_n2)+1
     iz=(parallel%mygrid_range(1)-1)/global_n1/global_n2+1
     ! print*,'xlt-test',ix,iy,iz,'myid',parallel%myid
     in=0 !parallel%mygrid_range(1)-1
     do i=ix,global_n1
        mVector(3) = iz - 1
        mVector(2) = iy - 1
        mVector(1) = i - 1
        in=in+1
        rPoint = MATMUL(slat,mVector)
        grid%rVec(1:3,in) = rPoint(:)
        grid%rVec(4,in) = Norm(rPoint)
     enddo
     z_up=(parallel%mygrid_range(2)-1)/global_n1/global_n2+1
     d_z=ceiling((z_up-iz)/(z_up-iz+1.d0))
     do j=iy+1,global_n2*d_z
        do i=1,global_n1
           mVector(3) = iz - 1
           mVector(2) = j - 1
           mVector(1) = i - 1
           in=in+1
           rPoint = MATMUL(slat,mVector)
           grid%rVec(1:3,in) = rPoint(:)
           grid%rVec(4,in) = Norm(rPoint)
        enddo
     enddo
     do k=iz+1,z_up-1
        do j=1,global_n2
           do i=1,global_n1
              mVector(3) = k - 1
              mVector(2) = j - 1
              mVector(1) = i - 1
              in=in+1
              rPoint = MATMUL(slat,mVector)
              grid%rVec(1:3,in) = rPoint(:)
              grid%rVec(4,in) = Norm(rPoint)
           enddo
        enddo
     enddo
     ix=mod(mod(parallel%mygrid_range(2)-1,global_n1*global_n2),global_n1)+1
     !> in{0|1},out{iy|1}
     y_down=(1-iy-1)*d_z+iy+1
     iy=mod((parallel%mygrid_range(2)-1)/global_n1,global_n2)+1
     iz=(parallel%mygrid_range(2)-1)/global_n1/global_n2+1
     do j=y_down,iy-1
        do i=1,global_n1
           mVector(3) = iz - 1
           mVector(2) = j - 1
           mVector(1) = i - 1
           in=in+1
           rPoint = MATMUL(slat,mVector)
           grid%rVec(1:3,in) = rPoint(:)
           grid%rVec(4,in) = Norm(rPoint)
        enddo
     enddo
     do i=1,ix
        mVector(3) = iz - 1
        mVector(2) = iy - 1
        mVector(1) = i - 1
        in=in+1
        rPoint = MATMUL(slat,mVector)
        grid%rVec(1:3,in) = rPoint(:)
        grid%rVec(4,in) = Norm(rPoint)
     enddo
# 640 "Grid_module.f90"
  ENDSUBROUTINE FillRTable
  !------------------Fill iso r_Table---------------------
  SUBROUTINE FillRTable_iso()
     USE struct_module , ONLY : lat_mat,recip_lat
     USE math , ONLY : Norm
     IMPLICIT NONE
     !
     INTEGER(I4B) :: ix,iy,iz,I
     REAL(DP)     :: mVector(3),rPoint(3)
     REAL(DP)     :: slat(3,3)
     !Initial
     slat(:,1)=lat_mat(:,1)/n1
     slat(:,2)=lat_mat(:,2)/n2
     slat(:,3)=lat_mat(:,3)/n3
     !r mesh
     I=0
     DO iz = 1,n3
        mVector(3) = iz - 1
     DO iy = 1,n2
        mVector(2) = iy - 1
     DO ix = 1,n1
        mVector(1) = ix - 1
        !
        I=I+1
        rPoint = MATMUL(slat,mVector)
        ! Assign this point in the qTable and qVectors.
        grid%rVec(1:3,I) = rPoint(:)
        grid%rVec(4,I) = Norm(rPoint)
     ENDDO
     ENDDO
     ENDDO
   ENDSUBROUTINE FillRTable_Iso
  !---------------------------------------------------
  SUBROUTINE build_ISO_sphere_grid()
      !USE struct_info, ONLY: Rermax
      USE struct_module
      USE parameters , ONLY : RadiusMax,NSPIN,Nstates,IDinit,finite_order,cell_shape,cell_thick
      USE succeed , ONLY : init_succeed_rho_real,Llastrho
      USE m_time_evaluate, ONLY: memory_sum,memory_free

      USE smpi_math_module , ONLY : grid_split,parallel,sphere,diff_map,&
           & grid_sphere_init,destroy_diff_map
      !> local
      INTEGER(I4B),ALLOCATABLE :: n_local(:)

      INTEGER(I4B) :: ix,iy,iz
      REAL(DP) :: orig(3),cosp,sinp
      INTEGER(I4B),ALLOCATABLE :: xyz(:,:)
      INTEGER(I4B) :: len

      !> set the z range
      INTEGER(I4B) :: z_range(2)

      !>>>>>>>>>>>>>>>>>>>>>>>>>START >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL matchposcar_ISO()
! #ifndef 1
      CALL destroy_ISO_sphere_grid()
! #endif
      ALLOCATE(dr(n1,n2,n3))    !r for grid
      ALLOCATE(cost(n1,n2,n3))  !cos(theta)
      ALLOCATE(sint(n1,n2,n3))  !sin(theta)
      ALLOCATE(cns(n1,n2,n3))   !exp(i*phi)
      ALLOCATE(indx(n1,n2,n3))  !logical for judge grid points which is in sphere
      call memory_sum('dr,cost,sint,cns,indx',5*real(size(dr),DP)*DP)
     !locale variables for assignment
      ALLOCATE(xyz(3,n1*n2*n3))
      call memory_sum('build_sphere_local_xyz',real(size(xyz),DP)*DP)

      CALL destroy_diff_map()
      ALLOCATE(n_local(n1*n2*n3))
      call memory_sum('build_sphere_local_n_local',real(size(n_local),DP)*I4B)
      ALLOCATE(diff_map%nz_map(2,n3))
      call memory_sum('diff_map_nz_map',real(size(diff_map%nz_map),DP)*I4B)

      !cubic cell?
      IF(n1/=n2.OR.n2/=n3.OR.n1/=n3) THEN
         WRITE(6,*) 'init_ISO: STOP!!! need cubic cell'
         STOP
      ENDIF
      orig(:)=0.5d0*gap(1)*(n1+2)
      !======================================
      !##XLT OBTAIN GRID INFORMATION
      !print*,"orig",orig
      !print*,"struct%poscar",struct%poscar
      !======================================
      rho_calc%OneDLength=0

      CALL get_z_range(z_range,gap(3),n3,cell_thick,cell_shape)

      DO iz=1,n3

      diff_map%nz_map(1,iz)=rho_calc%OneDLength+1

      DO iy=1,n2
      DO ix=1,n1
         CALL set_car2spe(ORIG,ix,iy,iz,dr(ix,iy,iz),cost(ix,iy,iz),sint(ix,iy,iz),cosp,sinp)
         cns(ix,iy,iz)=CMPLX(cosp,sinp)

         IF(dr(ix,iy,iz) <= RadiusMax)THEN

            if (ix >= z_range(1).and.ix <= z_range(2))then

            indx(ix,iy,iz) = 1
            rho_calc%OneDLength=rho_calc%OneDLength+1
            xyz(:,rho_calc%OneDLength)=(/ix,iy,iz/)

            n_local(rho_calc%OneDLength)=ix+(iy-1)*n1+(iz-1)*n1*n2
            endif

         ELSE
            indx(ix,iy,iz) = 0
         ENDIF

      ENDDO
      ENDDO

      diff_map%nz_map(2,iz)=rho_calc%OneDLength
      ! if(diff_map%nz_map(2,iz) == 0)diff_map%nz_map(2,iz)=1
      if(diff_map%nz_map(1,iz) > diff_map%nz_map(2,iz))diff_map%nz_map(1,iz)=diff_map%nz_map(2,iz)

      ENDDO
! #ifdef 1
!     ! if(parallel%isroot)print*,diff_map%nz_map
!       CALL destroy_ISO_sphere_grid()
! #endif
      !allocate density sphere
# 784 "Grid_module.f90"

      ! if(parallel%mygrid_range(1) == -2)CALL grid_split(rho_calc%OneDLength,parallel%numprocs&
      !      & , parallel%myid, parallel%mygrid_range,parallel%recvcounts,parallel%displs,&
      !      & parallel%global_gridrange)
      CALL grid_split(rho_calc%OneDLength,parallel%numprocs,parallel%comm&
           & , parallel%myid, parallel%mygrid_range,parallel%recvcounts,parallel%displs,&
           & parallel%global_gridrange)
      len=parallel%mygrid_range(3)
      !> init the rho to succeed the rho calculated by last ion step
      if(.not.Llastrho)CALL init_succeed_rho_real(len,natom,Nstates,Nspin,dvol)
      !> init sphere grid
      ALLOCATE(rho_calc%x(len))
      ALLOCATE(rho_calc%y(len))
      ALLOCATE(rho_calc%z(len))
      ALLOCATE(rho_calc%OneDVeff(len))
      ALLOCATE(rho_calc%OneDSphere(len,Nspin))
      ALLOCATE(rho_calc%OneDwvf(len,Nstates,NSPIN))
      ALLOCATE(rho_calc%n(len))
      ALLOCATE(rho_calc%OneDval(Nstates,NSPIN))
      ALLOCATE(sphere%map3d(4,rho_calc%OneDLength))
      call memory_sum('rho_calc',real(size(rho_calc%x),DP)*I4B*4&
           &+size(rho_calc%OneDVeff)*DP&
           &+size(rho_calc%OneDSphere)*DP&
           &+size(rho_calc%OneDwvf)*DP&
           &+size(rho_calc%OneDval)*DP&
           &+size(sphere%map3d)*I4B)

      IF(IDinit==1)ALLOCATE(rho_calc%initMO(rho_calc%OneDLength,Nstates))
      !> assign the grid map





      rho_calc%x=xyz(1,parallel%mygrid_range(1):parallel%mygrid_range(2))
      rho_calc%y=xyz(2,parallel%mygrid_range(1):parallel%mygrid_range(2))
      rho_calc%z=xyz(3,parallel%mygrid_range(1):parallel%mygrid_range(2))
      rho_calc%n=n_local(parallel%mygrid_range(1):parallel%mygrid_range(2))
      sphere%map3d(1:3,:)=xyz(:,1:rho_calc%OneDLength)
      sphere%map3d(4,:)=n_local(1:rho_calc%OneDLength)
      sphere%length=rho_calc%OneDLength
      !> grid communction parameters set
      !>> get the cores need to communction with i-th core and
      !>> the size of their data to be transport
      CALL grid_sphere_init(n1,n2,n3,finite_order)
      ! allocate(grid_diff_map%local_nz(2,numprocs))
      ! grid_diff_map%local_nz=0
      ! grid_diff_map(1,parallel%myid+1)=xyz(parallel%mygrid_range(1),3)
      ! grid_diff_map(2,parallel%myid+1)=xyz(parallel%mygrid_range(2),3)
      ! parallel
      ! CALL grid_comm_init(grid_diff_map)
      ALLOCATE(rho_calc%vlpp(len))
      call memory_sum('rho_calc%vlpp',real(len,DP)*DP)
      call memory_free('build_sphere_local_n_local',real(size(n_local),DP)*I4B)
      DEALLOCATE(n_local)

      call memory_free('build_sphere_local_xyz',real(size(xyz),DP)*DP)
      DEALLOCATE(xyz)
  !<<<<<<<<<<<<<<<<<<<<<<<<< END  F_DR<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDSUBROUTINE build_ISO_sphere_grid
!-----------------------DIVIDER-LINE--------------------------
  SUBROUTINE matchposcar_ISO()
    USE struct_module
    IMPLICIT NONE
    INTEGER(I4B) :: i
    !=========================================================
    !##move atoms in line with the grid coordinate system
    DO i=1,natom,1
      !## MOVE ATOMS POSITION
      struct%poscar(:,i)=struct%poscar(:,i)+gap(:)
    ENDDO
    DO i=1,n,1
      !## MOVE GRID r POSITION
      grid%rVec(1:3,i)=grid%rVec(1:3,i)+gap(:)
    ENDDO
    !=========================================================
  END SUBROUTINE matchposcar_ISO
!-----------------------PARTING-LINE--------------------------
  !-----------------------PARTING-LINE--------------------------
  !-----------------------DIVIDER-LINE--------------------------
  SUBROUTINE set_car2spe(ORIG,x,y,z,r,cost,sint,cosp,sinp)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: x,y,z
    REAL(8),INTENT(IN) :: ORIG(3)
    REAL(8),INTENT(OUT) :: r,cost,sint,cosp,sinp
    REAL(8) ::rr,rx,ry,rz
    INTEGER :: l
    !>>>>>>>>>>>>>>>>>>>>>>>start CAR2SPE>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    rx=gap(1)*REAL(x,8)-ORIG(1)
    ry=gap(2)*REAL(y,8)-ORIG(2)
    rz=gap(3)*REAL(z,8)-ORIG(3)
    r=SQRT(rx**2+ry**2+rz**2)
    rr=SQRT(rx**2+ry**2)
    IF(r>0.d0)THEN
       cost=rz/r
       sint=rr/r
    ELSE
       !cost=1.d0
       !sint=0.d0
       WRITE(6,*) 'init_ISO_sphere:the origin on grid'
       STOP
    ENDIF

    IF(rr>0.d0)THEN
      cosp=rx/rr
      sinp=ry/rr
    ELSE
      !cosp=1.d0
      !sinp=0.d0
       WRITE(6,*) 'init_ISO_sphere:the origin on grid'
       STOP
    ENDIF
    !<<<<<<<<<<<<<<<<<<<<<<<<<end CAR2SPE<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE set_car2spe
  !-----------------------PARTING-LINE--------------------------
  !-----------------------DIVIDER-LINE--------------------------
  SUBROUTINE destroy_ISO_sphere_grid()
    USE m_time_evaluate, ONLY: memory_free
    IMPLICIT NONE
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !coordinates in sphere
    if(allocated(dr))then
       call memory_free('dr cost sint cns indx',real(size(dr),DP)*5*DP)
       DEALLOCATE(dr, cost, sint, cns, indx)
    endif
    ! DEALLOCATE(dr, cost, sint, cns, indx,xyz)
    !coordinates in grid and OneDimension
    IF(ALLOCATED(rho_calc%x))THEN
       call memory_free('rho_calc_x',real(size(rho_calc%x),DP)*I4B)
       DEALLOCATE(rho_calc%x)
    ENDIF
    IF(ALLOCATED(rho_calc%y))THEN
       call memory_free('rho_calc_y',real(size(rho_calc%y),DP)*I4B)
       DEALLOCATE(rho_calc%y)
    ENDIF
    IF(ALLOCATED(rho_calc%z))THEN
       call memory_free('rho_calc_z',real(size(rho_calc%z),DP)*I4B)
       DEALLOCATE(rho_calc%z)
    ENDIF
    IF(ALLOCATED(rho_calc%OneDVeff))THEN
       call memory_free('rho_calc_OneDVeff',real(size(rho_calc%OneDVeff),DP)*DP)
       DEALLOCATE(rho_calc%OneDVeff)
    ENDIF
    IF(ALLOCATED(rho_calc%OneDSphere))THEN
       call memory_free('rho_calc_OneDSphere',real(size(rho_calc%OneDSphere),DP)*DP)
       DEALLOCATE(rho_calc%OneDSphere)
    ENDIF
    IF(ALLOCATED(rho_calc%OneDwvf))THEN
       call memory_free('rho_calc_OneDwvf',real(size(rho_calc%OneDwvf),DP)*DP)
       DEALLOCATE(rho_calc%OneDwvf)
    ENDIF
    IF(ALLOCATED(rho_calc%OneDval))THEN
       call memory_free('rho_calc_OneDval',real(size(rho_calc%OneDval),DP)*DP)
       DEALLOCATE(rho_calc%OneDval)
    ENDIF
    IF(ALLOCATED(rho_calc%initMO))THEN
       call memory_free('rho_calc_initMO',real(size(rho_calc%initMO),DP)*DP)
       DEALLOCATE(rho_calc%initMO)
    ENDIF

    IF(ALLOCATED(rho_calc%n))THEN
       call memory_free('rho_calc_n',real(size(rho_calc%n),DP)*I4B)
       DEALLOCATE(rho_calc%n)
    ENDIF
    IF(ALLOCATED(rho_calc%vlpp))THEN
       call memory_free('rho_calc_vlpp',real(size(rho_calc%vlpp),DP)*DP)
       DEALLOCATE(rho_calc%vlpp)
    ENDIF

    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


  END SUBROUTINE destroy_ISO_sphere_grid
  !-----------------------PARTING-LINE--------------------------
  !-----------------------DIVIDER-LINE--------------------------
  SUBROUTINE ISO_Vsphere(Vgrid,Vsphere)
    USE parameters ,  ONLY : NSPIN
    IMPLICIT NONE
    INTEGER(I4B)         :: i, j
    REAL(DP),INTENT(IN)  :: Vgrid(:,:,:,:)
    REAL(DP),INTENT(OUT) :: Vsphere(:,:)
    Vsphere=0.d0
    DO j=1,NSPIN,1
    DO i=1,rho_calc%OneDLength
      Vsphere(i,j)=Vgrid(rho_calc%x(i),rho_calc%y(i),rho_calc%z(i),j)
    ENDDO
    ENDDO
  END SUBROUTINE ISO_Vsphere
  !-----------------------PARTING-LINE--------------------------
  !-----------------------DIVIDER-LINE--------------------------
  SUBROUTINE ISO_Rho2grid(Rhosphere,Rhogrid)
    USE parameters ,  ONLY : NSPIN

    USE smpi_math_module , ONLY : parallel

    IMPLICIT NONE
    INTEGER(I4B)         :: i, j
    REAL(DP),INTENT(OUT) :: Rhogrid(:,:,:,:)
    REAL(DP),INTENT(IN)  :: Rhosphere(:,:)
    Rhogrid=0.d0
    DO j=1,NSPIN,1



    DO i=1,parallel%mygrid_range(3)

      Rhogrid(rho_calc%x(i),rho_calc%y(i),rho_calc%z(i),j)=Rhosphere(i,j)
    ENDDO
    ENDDO
  END SUBROUTINE ISO_Rho2grid
  !-----------------------DIVIDER-LINE--------------------------

  SUBROUTINE parallel_s2g(p,thrq)
    USE smpi_math_module, ONLY : parallel,mpi_real8,mpi_sum,mpinfo
    IMPLICIT NONE
    REAL(DP) :: p(parallel%mygrid_range(3))
    REAL(DP) :: thrq(n1,n2,n3)
    !> local
    REAL(DP)     :: sphere(rho_calc%OneDLength)
    INTEGER(I4B) :: i,ix,iy,iz
    REAL(DP) :: thrq_local(n1,n2,n3)

    !> obtain the contact sphere
    !call MPI_ALLGATHERV(p,parallel%mygrid_range(3),MPI_REAL8,sphere,&
    !     parallel%recvcounts,parallel%displs,MPI_REAL8,parallel%comm,mpinfo)

    !> convert sphere to grid
    thrq_local=0.d0
    thrq=0.d0
    do i=1,parallel%mygrid_range(3),1
       ix=rho_calc%x(i)
       iy=rho_calc%y(i)
       iz=rho_calc%z(i)
       thrq_local(ix,iy,iz)=p(i)
    enddo
    !> obtain the contact grid
    CALL MPI_ALLREDUCE(thrq_local,thrq,n,MPI_REAL8,MPI_SUM,parallel%comm,mpinfo)
    ! CALL smpi_reduce_sum_real(thrq,n)
  ENDSUBROUTINE parallel_s2g
  !-------------------------------------------------------------------------
  SUBROUTINE parallel_g2s(p,thrq)
    !> assign the sphere snippet array by a global grid array
    USE smpi_math_module, ONLY :parallel
    IMPLICIT NONE
    REAL(DP) :: p(parallel%mygrid_range(3))
    REAL(DP) :: thrq(n1,n2,n3)
    !> local
    INTEGER(I4B) :: i,ix,iy,iz

    !> convert grid to sphere
    p=0.d0
    do i=1,parallel%mygrid_range(3),1
       ix=rho_calc%x(i)
       iy=rho_calc%y(i)
       iz=rho_calc%z(i)
       p(i)=thrq(ix,iy,iz)
    enddo

  ENDSUBROUTINE parallel_g2s
  !-------------------------------------------------------------------------
  SUBROUTINE get_z_range(my_z_range,delta_z,ngrid_z,cell_thick,cell_shape)
    IMPLICIT NONE
    INTEGER(I4B) :: my_z_range(2)
    REAL(DP)     :: delta_z
    INTEGER(I4B) :: ngrid_z
    REAL(DP)     :: cell_thick
    INTEGER(I4B) :: cell_shape

    if(cell_shape==0)then
       my_z_range(1)=1
       my_z_range(2)=ngrid_z
    elseif(cell_shape==1)then
       my_z_range(1)=int((ngrid_z*delta_z/2.d0+delta_z-cell_thick*ang2bohr)/delta_z,I4B)
       my_z_range(2)=int((ngrid_z*delta_z/2.d0+delta_z+cell_thick*ang2bohr)/delta_z,I4B)
       if(my_z_range(1) < 1)my_z_range(1)=1
       if(my_z_range(2) > ngrid_z)my_z_range(2)=ngrid_z
    else
       print *,'grid_moduel::grt_z_range:701, parameter error'
       stop
    endif
  ENDSUBROUTINE get_z_range

  !===================================================================
  SUBROUTINE confirm_iso_radius()

    USE smpi_math_module
    USE read_module, ONLY: out_CONCAR

    USE struct_module ,ONLY :struct,natom,lat_mat
    USE parameters    ,ONLY :Lcell,Lcellbohr,ang2bohr,IsoRmaxbohr,IsoRmax&
         & ,RadiusMax,Lpbc,bohr2ang,Lradius_auto
    USE pspot_module  , ONLY :max_rcut
    IMPLICIT NONE
    INTEGER(I4B)  :: i
    REAL(DP)      :: max_r,r_temp,r
    REAL(DP)      :: len_trans
    !========================================
    !> find MAX(r)



       IF(.not.Lpbc.and.parallel%isroot)THEN

          max_r=0
          do i=1,natom,1
             r=sqrt(sum((struct%poscar(:,i)-Lcell*ang2bohr*0.5d0)**2))
             if(r>max_r)max_r=r
          enddo
          !> set the radius
          r_temp=max_r+max_rcut+2.5d0*ang2bohr
       ENDIF

       CALL MPI_BCAST(r_temp,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
       CALL MPI_BCAST(max_r,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
       CALL MPI_BCAST(max_rcut,1,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)

       !print *,r_temp,max_r,max_rcut
       ! IF(RadiusMax<r_temp.and.Lradius_auto)THEN
       IF(Lradius_auto)THEN
          IsoRmaxbohr=r_temp
          IsoRmax=IsoRmaxbohr*bohr2ang
          RadiusMax=IsoRmaxbohr
          Lradius_auto=.False.
          !> prompt the setting



          if(parallel%isroot)then
             print*,'PROMPT : the radius is too small,',IsoRmax,'(ang) is properly'
          endif

          len_trans=(IsoRmaxbohr*2.d0-Lcellbohr)/2.d0
          if(len_trans/=0)then
             Lcell=IsoRmax*2.d0
             Lcellbohr=RadiusMax*2.d0
             struct%poscar=struct%poscar+len_trans
             struct%pos=struct%poscar/Lcellbohr
             lat_mat=0.d0
             do i=1,3,1
                lat_mat(i,i)=Lcellbohr
             enddo
          endif
          ! print*,'poscar',struct%poscar
          ! print*,'pos',struct%pos



       if(parallel%isroot)call out_CONCAR('pos-2')
    ENDIF

  ENDSUBROUTINE confirm_iso_radius
  !>-----------------------------------------
  SUBROUTINE reshape_center() !{{{
    use parameters  , only : lcell,Lcellbohr
    USE struct_module, ONLY : struct,natom,lat_mat
    ! USE math ,ONLY : norm
    implicit none
    real(dp)                ::  lat_mat1(3,3), xyzbar(3)
    integer(i4b)  ::  i
    !> for judge atom position and sphere radius
    REAL(DP) :: posdis(natom)
    !-----------------------------------------------
    !> reset the direct pos
    do i = 1, 3
       xyzbar(i) = sum(struct%pos(i,:))/natom
       struct%pos(i,:) = struct%pos(i,:) - xyzbar(i)
    enddo

    ! do i = 1, natom
    !    posdis(i) = norm(matmul(lat_mat,struct%pos(:,i) ) )
    ! enddo
    struct%pos = 0.5d0 + struct%pos

    !> reset the poscar
    lat_mat1 = 0.d0
    do i = 1, 3
       lat_mat1(i,i) = lcell*ang2bohr
    enddo

    do i = 1 , natom
       struct%poscar(:,i) = matmul(lat_mat1(:,:) , struct%pos(:,i) )
    enddo
  ENDSUBROUTINE reshape_center!}}}
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine reset_poscar(new_gap,new_ni) !{{{
    !> reset the lcell by multiply new_gap and new_ni
    USE struct_module, ONLY : struct,natom
    implicit none
    !> IN/OUT
    REAL(DP)      :: new_gap
    INTEGER(I4B)  :: new_ni
    !> LOCAL
    REAL(DP)      ::  lat_mat1(3,3), xyzbar(3)
    REAL(DP)      ::  new_lcell
    integer(i4b)  ::  i
    !> for judge atom position and sphere radius
    REAL(DP) :: posdis(natom)
    !-----------------------------------------------

    !> reset the direct pos
    new_lcell = new_gap * new_ni * ang2bohr
    do i = 1, natom
       struct%pos(:,i) = struct%poscar(:,i)/new_lcell
    enddo
    do i = 1, 3
       xyzbar(i) = sum(struct%pos(i,:))/natom
       struct%pos(i,:) = struct%pos(i,:) - xyzbar(i)
    enddo

    ! do i = 1, natom
    !    posdis(i) = norm(matmul(lat_mat,struct%pos(:,i) ) )
    ! enddo
    struct%pos = 0.5d0 + struct%pos

    !> reset the poscar
    lat_mat1 = 0.d0
    do i = 1, 3
       lat_mat1(i,i) = new_lcell
    enddo

    do i = 1 , natom
       struct%poscar(:,i) = matmul(lat_mat1(:,:) , struct%pos(:,i) )
    enddo
  end subroutine reset_poscar!}}}


  SUBROUTINE build_parallel_cubic_grid()
    USE smpi_math_module, ONLY: grid_split,parallel,diff_map,grid_cubic_init,destroy_diff_map
    USE parameters, ONLY: NSPIN,finite_order,Nstates
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
    nsn=n*NSPIN

    !> reallocate rhoS for parallel
    ! if(allocated(grid%rhoS))deallocate(grid%rhoS)
    ! if(allocated(grid%vlpp))deallocate(grid%vlpp)
    ! allocate(grid%rhoS(n1,n2,n3,Nspin))
    ! allocate(grid%vlpp(local_n1,local_n2,local_n3))
    ! print *,'rank',parallel%rankx,'local_n[123]',n1,n2,n3
    ! print *,'rank',parallel%rankx,'my_grid_range',parallel%mygrid_range
    ! i=mod(mod(parallel%mygrid_range(1)-1,global_n1*global_n2),global_n1)+1
    ! j=mod((parallel%mygrid_range(1)-1)/global_n1,global_n2)+1
    ! k=(parallel%mygrid_range(1)-1)/global_n1/global_n2+1
    ! ! print *,'rank',parallel%rankx,'grid(ijk) start',i,j,k
    ! i=mod(mod(parallel%mygrid_range(2)-1,global_n1*global_n2),global_n1)+1
    ! j=mod((parallel%mygrid_range(2)-1)/global_n1,global_n2)+1
    ! k=(parallel%mygrid_range(2)-1)/global_n1/global_n2+1
    ! print *,'rank',parallel%rankx,'grid(ijk) end',i,j,k
    call grid_cubic_init(global_n1,global_n2,global_n3,global_n,finite_order)
  ENDSUBROUTINE build_parallel_cubic_grid


  SUBROUTINE rho_trans1D(rho3d,rho1d)
    USE smpi_math_module, ONLY: parallel
    implicit none
    REAL(DP),intent(in) :: rho3d(n1,n2,n3)
    REAL(DP),intent(out):: rho1d(n)
    !>local
    INTEGER(I4B) :: ix,iy,iz,in,i,j,k
    INTEGER(I4B) :: d_z,y_down

    ix=mod(mod(parallel%mygrid_range(1)-1,global_n1*global_n2),n1)+1
    iy=mod((parallel%mygrid_range(1)-1)/global_n1,global_n2)+1
    iz=1
    in=0 !parallel%mygrid_range(1)-1
    do i=ix,global_n1
       in=in+1
       rho1d(in) = rho3d(i,iy,iz)
    enddo
    d_z=ceiling((n3-iz)/(n3-iz+1.d0))
    do j=iy+1,global_n2*d_z
       do i=1,global_n1
          in=in+1
          rho1d(in) = rho3d(i,j,iz)
       enddo
    enddo
    do k=2,n3-1
       do j=1,global_n2
          do i=1,global_n1
             in=in+1
             rho1d(in) = rho3d(i,j,k)
          enddo
       enddo
    enddo
    ix=mod(mod(parallel%mygrid_range(2)-1,global_n1*global_n2),global_n1)+1
    y_down=(1-iy-1)*d_z+iy+1
    iy=mod((parallel%mygrid_range(2)-1)/global_n1,global_n2)+1
    iz=n3
    do j=y_down,iy-1
       do i=1,global_n1
          in=in+1
          rho1d(in) = rho3d(i,j,iz)
       enddo
    enddo
    do i=1,ix
       in=in+1
       rho1d(in) = rho3d(i,iy,iz)
    enddo
  END SUBROUTINE rho_trans1D
  SUBROUTINE rho_trans3D(rho1d,rho3d)
    USE smpi_math_module, ONLY: parallel
    implicit none
    REAL(DP),intent(in):: rho1d(n)
    REAL(DP),intent(out) :: rho3d(n1,n2,n3)
    !>local
    INTEGER(I4B) :: ix,iy,iz,in,i,j,k
    INTEGER(I4B) :: d_z,y_down

    ix=mod(mod(parallel%mygrid_range(1)-1,global_n1*global_n2),n1)+1
    iy=mod((parallel%mygrid_range(1)-1)/global_n1,global_n2)+1
    iz=1
    in=0 !parallel%mygrid_range(1)-1
    do i=ix,global_n1
       in=in+1
       rho3d(i,iy,iz) = rho1d(in)
    enddo
    d_z=ceiling((n3-iz)/(n3-iz+1.d0))
    do j=iy+1,global_n2*d_z
       do i=1,global_n1
          in=in+1
          rho3d(i,j,iz) = rho1d(in)
       enddo
    enddo
    do k=2,n3-1
       do j=1,global_n2
          do i=1,global_n1
             in=in+1
             rho3d(i,j,k) = rho1d(in)
          enddo
       enddo
    enddo
    ix=mod(mod(parallel%mygrid_range(2)-1,global_n1*global_n2),global_n1)+1
    y_down=(1-iy-1)*d_z+iy+1
    iy=mod((parallel%mygrid_range(2)-1)/global_n1,global_n2)+1
    iz=n3
    do j=y_down,iy-1
       do i=1,global_n1
          in=in+1
          rho3d(i,j,iz) = rho1d(in)
       enddo
    enddo
    do i=1,ix
       in=in+1
       rho3d(i,iy,iz) = rho1d(in)
    enddo
  END SUBROUTINE rho_trans3D

  !-------------------divide line ---------------------
  FUNCTION fft_sph_r2c(array_r) result(array_c)
    USE smpi_math_module, ONLY: parallel,MPI_REAL8,mpinfo,smpi_exit
    USE FOURIER, ONLY: FFT
    implicit none
    !> in/out
    REAL(DP),intent(in)  :: array_r(:,:,:)
    !DIMENSION(ng1,ng2,parallel%local_z)
    COMPLEX(DP) :: array_c(ng1,ng2,parallel%local_z)
    !> local
    INTEGER(I4B) :: ix,iy
    REAL(DP)  :: array_FFT(global_n1,global_n2,parallel%local_z)

    ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
    iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1

    CALL MPI_ALLTOALLV(array_r(ix,iy,1),parallel%fft_scount&
         & ,parallel%fft_sdispls,MPI_REAL8&
         & ,array_FFT,parallel%fft_rcount&
         & ,parallel%fft_rdispls,MPI_REAL8&
         & ,parallel%commx,mpinfo)
    array_c(:,:,:)=FFT(array_FFT)

  END FUNCTION fft_sph_r2c
  !-------------------divide line ---------------------
  FUNCTION fft_sph_c2r(array_c)result(array_r)
    USE smpi_math_module, ONLY: parallel,MPI_REAL8,mpinfo,smpi_exit
    USE FOURIER, ONLY: FFT
    implicit none
    !> in/out
    REAL(DP)  :: array_r(n1,n2,n3)
    !DIMENSION(ng1,ng2,parallel%local_z)
    COMPLEX(DP),intent(in) :: array_c(:,:,:)
    !> local
    INTEGER(I4B) :: ix,iy
    REAL(DP)  :: array_FFT(global_n1,global_n2,parallel%local_z)

    array_FFT=FFT(array_c)
    array_r=0.d0
    ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
    iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
    CALL MPI_ALLTOALLV(array_FFT,parallel%fft_rcount&
         & ,parallel%fft_rdispls,MPI_REAL8&
         & ,array_r(ix,iy,1),parallel%fft_scount&
         & ,parallel%fft_sdispls,MPI_REAL8&
         & ,parallel%commx,mpinfo)

  END FUNCTION fft_sph_c2r
  !-------------------divide line ---------------------
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<END<MODULE<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDMODULE grid_module

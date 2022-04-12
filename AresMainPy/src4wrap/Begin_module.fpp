# 1 "Begin_module.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "Begin_module.f90"
MODULE Begin_module
   USE constants
   IMPLICIT NONE
CONTAINS
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !-------------------PARTING LINE--------------------
  SUBROUTINE Initial_grid()
    use parameters, only:Lpbc
    implicit none

    if(lpbc)then
       call Initial_grid_per()
    else
       call Initial_grid_iso()
    endif

  ENDSUBROUTINE Initial_grid
  !-------------------PARTING LINE--------------------
  SUBROUTINE Initial_grid_per()
     USE math ,ONLY : lat2matrix,inv_33,det,dir2car
     USE struct_module , ONLY : lat_mat,recip_lat &
                            & , volume,charge_ave,ncharge &
                            & , struct
     USE parameters, ONLY : finite_order,KSPACING,LKP,Lpbc
     USE Grid_module , ONLY : Build_rgrid,Build_kgrid &
                          &  ,Build_eigen,n1,n2,n3 &
                          & ,nk1,nk2,nk3,nk,gap &
                          &,FillQTable,FillRTable &

                          &, global_n1,global_n2,global_n3 &
                          &,build_ISO_sphere_grid &
                          &,build_parallel_cubic_grid
     USE smpi_math_module, ONLY: parallel



     USE FOURIER, ONLY : FFT,PlanFFT,CleanFFT
     USE finite_module , ONLY : init_finite
     IMPLICIT NONE
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     CALL dir2car(struct%pos,struct%poscar,lat_mat)
     !set cell
     recip_lat(:,:)=2.d0*pi*TRANSPOSE(inv_33(lat_mat))
     volume = ABS(det(lat_mat))
     !average charge density
     charge_ave=ncharge/volume
     !creat grids
     CALL Build_rgrid()
     !creat k-points
     CALL Build_kgrid()
     !creat eigen-data
     CALL Build_eigen()
     !>>>Finite difference
     CALL init_finite(finite_order,gap)
     !<<<Finite difference
     !Plan for FFT

     CALL PlanFFT(global_n1,global_n2,global_n3)



     !recip grid mesh data
     CALL FillQTable()
     !real grid mesh data
     CALL FillRTable()
     !initial density
     ! CALL Initial_density()
     !

     if(parallel%isroot)then
        WRITE(6,*)'R-space-GRIDS:',global_n1,global_n2,global_n3
     !WRITE(6,'(A21,3I5,A6,I5)')'local R-space-GRIDS:',n1,n2,n3,',rank:',parallel%myid
     IF(LKP)THEN
        WRITE(6,*)'K-points-used:',nk
        !WRITE(6,*)'Waiting for this'
        !STOP
     ELSE
        WRITE(6,*)'K-space-GRIDS:',nk1,nk2,nk3
        IF(KSPACING>0.d0)THEN
           WRITE(6,*)'Num of K-used:',nk
        ELSE
           WRITE(6,*)'Only gamma point is used'
        ENDIF
     ENDIF
     endif      
# 101 "Begin_module.f90"
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE Initial_grid_per
  !-------------------PARTING LINE--------------------
  SUBROUTINE Initial_grid_iso()

    USE smpi_math_module

    USE math ,ONLY : lat2matrix,inv_33,det,dir2car
    USE struct_module , ONLY : lat_mat,recip_lat &
         & , volume,charge_ave,ncharge &
         & , struct
    USE parameters, ONLY : finite_order,KSPACING,LKP,Lpbc
    USE Grid_module , ONLY : Build_rgrid_iso,Build_kgrid &
         &  ,Build_eigen,n1,n2,n3 &
         & ,nk1,nk2,nk3,nk,gap &
         &,FillQTable,FillRTable_iso &
         &,build_ISO_sphere_grid &
         &,reshape_center,confirm_iso_radius,rho_calc
    USE FOURIER, ONLY : FFT,PlanFFT,CleanFFT
    USE finite_module , ONLY : init_finite
    IMPLICIT NONE
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call reshape_center()
    CALL confirm_iso_radius()
    CALL dir2car(struct%pos,struct%poscar,lat_mat)
    !set cell
    recip_lat(:,:)=2.d0*pi*TRANSPOSE(inv_33(lat_mat))
    volume = ABS(det(lat_mat))
    !average charge density
    charge_ave=ncharge/volume
    !creat grids
    CALL Build_rgrid_iso()
    !creat k-points
    CALL Build_kgrid()
    !creat eigen-data
    CALL Build_eigen()
    !>>>Finite difference
    CALL init_finite(finite_order,gap)
    !<<<Finite difference
    !Plan for FFT
    ! CALL PlanFFT_iso(n1,n2,n3)
    !recip grid mesh data
    CALL FillQTable()
    !real grid mesh data
    CALL FillRTable_iso()
    !initial density
    !CALL Initial_density()
    !
    !ISO shpere grid
    IF(.NOT.Lpbc)CALL build_ISO_sphere_grid()

    if(parallel%isroot)then

       WRITE(6,*)'R-space-GRIDS:',rho_calc%OneDLength
       ! IF(LKP)THEN
       !    WRITE(6,*)'K-points-used:',nk
       !    !WRITE(6,*)'Waiting for this'
       !    !STOP
       ! ELSE
       !    WRITE(6,*)'K-space-GRIDS:',nk1,nk2,nk3
       !    IF(KSPACING>0.d0)THEN
       !       WRITE(6,*)'Num of K-used:',nk
       !    ELSE
       !       WRITE(6,*)'Only gamma point is used'
       !    ENDIF
       ! ENDIF

    endif

    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE Initial_grid_iso
  !----------------Initial Charge density---------------------
  SUBROUTINE Initial_density()
     USE grid_module , ONLY : grid
     USE parameters , ONLY : LINRHO,NSPIN,LRadRho
     USE struct_module , ONLY : charge_ave

     USE smpi_math_module, ONLY:parallel

     IMPLICIT NONE
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     IF(LINRHO)THEN
        WRITE(6,*) 'Initialize density : by INRHO file'
        !read the density
        IF(NSPIN==1)THEN
           OPEN(107,FILE='INRHO')
               READ(107,*)
               READ(107,*) grid%rhoS
           CLOSE(107)
        ELSE
           OPEN(107,FILE='INRHO')
               READ(107,*)
               READ(107,*) grid%rhoS(:,:,:,1)
               READ(107,*)
               READ(107,*) grid%rhoS(:,:,:,2)
           CLOSE(107)
        ENDIF

     !=================================================
     !SHIELD
     ELSEIF(LRadRho)THEN

       if(parallel%isroot)WRITE(6,*) 'Initialize density : by pseudopotential files'



       CALL inichrg()
     !=================================================
     ELSE

        if(parallel%isroot)WRITE(6,*) 'Initialize density : by uniform density'



        !uniform density
        IF(NSPIN==1)THEN
           grid%rhoS(:,:,:,:)=charge_ave
        ELSE
           grid%rhoS(:,:,:,1)=charge_ave*0.48d0
           grid%rhoS(:,:,:,2)=charge_ave*0.52d0
        ENDIF

     ENDIF
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE Initial_density
  !------------------inichrg---------------------------
  SUBROUTINE inichrg()
     USE parameters , ONLY : NSPIN,PP_identifer
     USE math , ONLY : CubicSplineInterp,kahan_sum
     USE struct_module , ONLY : struct,naty,lat_mat,ncharge
     USE pspot_module , ONLY : psp
     USE Math, ONLY : interp

     USE smpi_math_module, ONLY: parallel,MPI_REAL8,mpinfo,MPI_SUM,smpi_exit
     USE m_time_evaluate, ONLY: filename
     USE grid_module , ONLY : n1,n2,n3,n,grid,dvol,global_n1,global_n2,global_n3



     IMPLICIT NONE
     !LOCAL
     REAL(DP) :: rho_t(n)
     INTEGER(I4B) :: nrep=1
     INTEGER(I4B) :: Ity,Ia,icx,icy,icz,I,ix,iy,iz
     REAL(DP) :: rat(3),ra(3),dist  &
          &  , acha , tele
     REAL(DP) :: r(3),f(3),c(3)
     INTEGER(I4B) :: m

     INTEGER(I4B) :: j,k,in
     REAL(DP),allocatable :: rho_temp(:,:,:)
     REAL(DP)     :: tele_local
     INTEGER(I4B) :: d_z,y_down
     REAL(DP) :: drcut=4*ang2bohr

     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     rho_t(:)=0.d0
     !Cycle all species of atoms
     DO Ity=1,naty
        !Cycle all atoms in this species
        DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
           ra(:)=struct%poscar(:,Ia)
           !all points
! #ifdef 1
!            DO I=parallel%mygrid_range(1),parallel%mygrid_range(2)
! #else
           DO I=1,n
! #endif
              !replicas
              DO icz=-nrep,nrep
              DO icy=-nrep,nrep
              DO icx=-nrep,nrep
                 rat(:)= icx*lat_mat(:,1)  &
                    & +  icy*lat_mat(:,2)  &
                    & +  icz*lat_mat(:,3)  &
                    & +  grid%rVec(1:3,I)  &
                    & -  ra(:)
                 dist=SQRT( rat(1)**2 + rat(2)**2 +rat(3)**2 )
                 IF(dist>drcut) CYCLE
                 !interpolate density
              ! m=0
              ! DO WHILE (psp(Ity)%r_real(m+1).lt.dist .and. (m.le.psp(Ity)%numps))!( m*psp(Ity)%rspacing.lt.dist )
              !    m=m+1
              ! ENDDO
              ! if(m<=1)m=2
              ! IF(dist.lt.psp(Ity)%rmax)THEN
              !    r=psp(Ity)%r_real(m-1:m+1)
              !    f=psp(Ity)%denr(m-1:m+1)
              !    acha=polynom(0,3,r,f,c,dist)
              ! ELSE
              !    r=psp(Ity)%r_real(psp(Ity)%numps-2:psp(Ity)%numps)
              !    f=psp(Ity)%denr(psp(Ity)%numps-2:psp(Ity)%numps)
              !    acha=polynom(0,3,r,f,c,dist)
              ! ENDIF
                 if(PP_identifer==0)then
                 acha=CubicSplineInterp(psp(Ity)%denr,psp(Ity)%ddden_dr2  &
                      &   ,psp(Ity)%rmax,psp(Ity)%rspacing,dist)
                 else
                 acha=interp(psp(Ity)%numps,psp(Ity)%denr &
                        &,psp(Ity)%r_real,dist)
                 ! acha=CubicSplineInterp(psp(Ity)%denr,psp(Ity)%ddden_dr2  &
                 !      &   ,psp(Ity)%rmax,psp(Ity)%rspacing,dist)
                 endif
                 rho_t(I) = rho_t(I) + acha
              ENDDO
              ENDDO
              ENDDO
           ENDDO !> I
           !
        ENDDO !> Ia
     ENDDO !> Ity
     !

     ! allocate(rho_temp(global_n1,global_n2,global_n3))
     ! call MPI_ALLGATHERV(rho_t(:),parallel%mygrid_range(3),&
     !      & MPI_REAL8,rho_temp,parallel%recvcounts,parallel%displs&
     !      &, MPI_REAL8, parallel%commx,mpinfo)
     ! tele=SUM(rho_temp)
     tele_local=kahan_SUM(size(rho_t),rho_t(:))*dvol
     ! tele_local=SUM(rho_t)*dvol
     CALL MPI_ALLREDUCE(tele_local,tele,1,MPI_REAL8,MPI_SUM,parallel%commx,mpinfo)



     !
     ! print*,'ncharge',ncharge,'tele',tele
     rho_t(:)=rho_t(:)*ncharge/tele
     !decide the initial density
     IF(NSPIN==1)THEN

        grid%rhoS=0.d0
        ! print*,'z2-z1',(parallel%mygrid_range(2)-1)/n1/n2+1-(parallel%mygrid_range(1)-1)/n1/n2-1+1,'local_n3',local_n3
        ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
        iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
        iz=1 !(parallel%mygrid_range(1)-1)/n1/n2+1
        in=0 !parallel%mygrid_range(1)-1
        do i=ix,n1
           in=in+1
           grid%rhoS(i,iy,iz,1)= rho_t(in)
           ! if(in==1)print*,'first',i,iy,iz
        enddo
        d_z=ceiling((n3-iz)/(n3-iz+1.d0))
        do j=iy+1,n2*d_z
           do i=1,n1
              in=in+1
              grid%rhoS(i,j,iz,1)= rho_t(in)
           enddo
        enddo
        do k=2,n3-1
           do j=1,n2
              do i=1,n1
                 in=in+1
                 grid%rhoS(i,j,k,1)= rho_t(in)
              enddo
           enddo
        enddo
        ix=mod(mod(parallel%mygrid_range(2)-1,n1*n2),n1)+1
        y_down=(1-iy-1)*d_z+iy+1
        iy=mod((parallel%mygrid_range(2)-1)/n1,n2)+1
        iz=n3 !(parallel%mygrid_range(2)-1)/n1/n2+1
        do j=y_down,iy-1
           do i=1,n1
              in=in+1
              grid%rhoS(i,j,iz,1)= rho_t(in)
           enddo
        enddo
        do i=1,ix
           in=in+1
           grid%rhoS(i,iy,iz,1)= rho_t(in)
        enddo
        !> end of local and end of global
        !>print*,'in',in,'mygrid_range(2)',parallel%mygrid_range(2)
# 384 "Begin_module.f90"
     ELSE

        ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
        iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
        iz=1 !(parallel%mygrid_range(1)-1)/n1/n2+1
        in=0 !parallel%mygrid_range(1)-1
        do i=ix,n1
           in=in+1
           grid%rhoS(i,iy,iz,1)=rho_t(in)*0.48d0
           grid%rhoS(i,iy,iz,2)=rho_t(in)*0.52d0
        enddo
        d_z=ceiling((n3-iz)/(n3-iz+1.d0))
        do j=iy+1,n2*d_z
           do i=1,n1
              in=in+1
              grid%rhoS(i,j,iz,1)=rho_t(in)*0.48d0
              grid%rhoS(i,j,iz,2)=rho_t(in)*0.52d0
           enddo
        enddo
        do k=2,n3-1
           do j=1,n2
              do i=1,n1
                 in=in+1
                 grid%rhoS(i,j,k,1)=rho_t(in)*0.48d0
                 grid%rhoS(i,j,k,2)=rho_t(in)*0.52d0
              enddo
           enddo
        enddo
        ix=mod(mod(parallel%mygrid_range(2)-1,n1*n2),n1)+1
        y_down=(1-iy-1)*d_z+iy+1
        iy=mod((parallel%mygrid_range(2)-1)/n1,n2)+1
        iz=n3 !(parallel%mygrid_range(1)-1)/n1/n2+1
        do j=y_down,iy-1
           do i=1,n1
              in=in+1
              grid%rhoS(i,j,iz,1)=rho_t(in)*0.48d0
              grid%rhoS(i,j,iz,2)=rho_t(in)*0.52d0
           enddo
        enddo
        do i=1,ix
           in=in+1
           grid%rhoS(i,iy,iz,1)=rho_t(in)*0.48d0
           grid%rhoS(i,iy,iz,2)=rho_t(in)*0.52d0
        enddo
        ! print*,'in',in,'mygrid_range(2)',parallel%mygrid_range(2)
# 441 "Begin_module.f90"
     ENDIF
# 457 "Begin_module.f90"
!OPEN(10086,FILE='intrho.dat')
!    WRITE(10086,*) grid%rhoS(:,:,:,:)
!CLOSE(10086)
!STOP
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE inichrg
  !-----------------------DIVIDER-LINE--------------------------
  SUBROUTINE init_ISO_density()

    USE smpi_math_module

    USE parameters    , ONLY : NSPIN,LRadRho,IDinit
    USE struct_module , ONLY : volume, ncharge
    USE grid_module   , ONLY : n1, n2, n3, rho_calc, grid
    IMPLICIT NONE
    INTEGER(I4B)      :: i, ix, iy, iz ,j
    !=====================================================================
    !##JUDGE THE INITIALIZE METHOD
# 484 "Begin_module.f90"
    IF(LRadRho.and.IDinit==0)THEN
      IF(parallel%isroot)WRITE(6,*) 'Initialize ISO density : by pseudopotential files'

      CALL inichrg_iso()
      !===============================================================
    ELSE
      !##uniform density



       IF(parallel%isroot)WRITE(6,*) 'Initialize ISO density : by uniform files'

      rho_calc%volume=volume/REAL(n1*n2*n3,DP)*rho_calc%OneDLength
      IF(NSPIN==1)THEN
        rho_calc%OneDSphere(:,:)=ncharge/rho_calc%volume
      ELSE
        rho_calc%OneDSphere(:,1)=ncharge/rho_calc%volume*0.48d0
        rho_calc%OneDSphere(:,2)=ncharge/rho_calc%volume*0.52d0
      ENDIF
      !> assignment to grid
# 513 "Begin_module.f90"
    ENDIF
    !=====================================================================
  END SUBROUTINE init_ISO_density
  !-----------------------DIVIDER-LINE--------------------------
  SUBROUTINE inichrg_iso()
     USE parameters , ONLY : NSPIN
     USE math , ONLY : CubicSplineInterp
     USE grid_module , ONLY : n1,n2,n3,n,grid,dvol,rho_calc,gap
     USE struct_module , ONLY : struct,naty,lat_mat,ncharge,natom
     USE pspot_module , ONLY : psp
     USE MathSplines   ,     ONLY: polynom

     USE smpi_math_module



     USE m_time_evaluate, ONLY: memory_sum,memory_free
     IMPLICIT NONE
     !LOCAL



     REAL(DP) :: rho_t(parallel%mygrid_range(3))

     INTEGER(I4B) :: nrep=0
     INTEGER(I4B) :: Ity,Ia,icx,icy,icz,I,ix,iy,iz,m
     REAL(DP) :: rat(3),ra(3),dist  &
              &  , acha , tele, coefficient(2), &
              &   threedot(3), c(3),r(3),f(3)

     REAL(DP)     :: tele_local
     integer(I4B) :: global_id
     INTEGER(I4B),allocatable :: atom_index(:)
     !> temp
     REAL(DP)    :: grid_rho(n)





     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

     call memory_sum('init_density_local',real(size(rho_t)+size(grid_rho),DP)*DP)



     rho_t=0.d0




     !> Cycle all species of atoms
     DO Ity=1,naty
        !Cycle all atoms in this species
        DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
           ra(:)=struct%poscar(:,Ia)
           !> all points



           DO I=1,parallel%mygrid_range(3)

              !> replicas
              rat(1)= rho_calc%x(I)*gap(1) - ra(1)
              rat(2)= rho_calc%y(I)*gap(1) - ra(2)
              rat(3)= rho_calc%z(I)*gap(1) - ra(3)
              dist=SQRT( rat(1)**2 + rat(2)**2 +rat(3)**2 )
              !> interpolate density
              m=0
              DO WHILE (psp(Ity)%r_real(m+1).lt.dist .and. (m.le.psp(Ity)%numps))!( m*psp(Ity)%rspacing.lt.dist )
                 m=m+1
              ENDDO
              IF(dist.lt.psp(Ity)%rmax)THEN
                 r=psp(Ity)%r_real(m-1:m+1)
                 f=psp(Ity)%denr(m-1:m+1)
                 acha=polynom(0,3,r,f,c,dist)
              ELSE
                 r=psp(Ity)%r_real(psp(Ity)%numps-2:psp(Ity)%numps)
                 f=psp(Ity)%denr(psp(Ity)%numps-2:psp(Ity)%numps)
                 acha=polynom(0,3,r,f,c,dist)
              ENDIF
              !> gauss distribution
              ! acha = 1.d0/(4*sqrt(2.d0*pi*a))**3*exp(-dist**2/(2*pi*a))
              rho_t(I)=rho_t(I) + acha



           ENDDO
        ENDDO !> Ia
     ENDDO !> Ity
# 648 "Begin_module.f90"
     !> assignment to grid
     ! grid%rhoS=0.d0
     ! IF(NSPIN==1)THEN
     !   DO i=1,parallel%mygrid_range(3)
     !      grid%rhoS(rho_calc%x(i),rho_calc%y(i),rho_calc%z(i),1)= &
     !         & rho_t(i)
     !   ENDDO
     !   grid_rho=reshape(grid%rhoS,(/n/))
     !   CALL smpi_reduce_sum_real_1d(grid_rho)
     !   grid%rhoS=reshape(grid_rho,(/n1,n2,n3,1/))
     ! ELSE
     !   DO i=1,parallel%mygrid_range(3)
     !     grid%rhoS(rho_calc%x(i),rho_calc%y(i),rho_calc%z(i),1)= &
     !         & rho_t(i)*0.48d0
     !     grid%rhoS(rho_calc%x(i),rho_calc%y(i),rho_calc%z(i),2)= &
     !         & rho_t(i)*0.52d0
     !   ENDDO
     !   grid_rho=0.d0
     !   grid_rho=reshape(grid%rhoS,(/2*n/))
     !   CALL smpi_reduce_sum_real_1d(grid_rho)
     !   grid%rhoS=reshape(grid_rho,(/n1,n2,n3,2/))
     ! ENDIF

     !> Normalize the charge density
     tele_local=SUM(rho_t(:))*dvol
     ! print*,"myid",parallel%myid,'tele_local',tele_local
     CALL MPI_ALLREDUCE(tele_local,tele,1,MPI_REAL8,MPI_SUM,parallel%comm,mpinfo)
     ! print*,"myid",parallel%myid,'tele',tele
     rho_t(:)=rho_t(:)*ncharge/tele
     rho_calc%OneDSphere=0.d0

     !> assignment the charge density
     IF(NSPIN==1)THEN
        CALL dcopy(parallel%mygrid_range(3),rho_t,1,rho_calc%OneDSphere,1)
     ELSE
        CALL daxpy(parallel%mygrid_range(3),0.48d0,rho_t,1,rho_calc%OneDSphere(:,1),1)
        CALL daxpy(parallel%mygrid_range(3),0.52d0,rho_t,1,rho_calc%OneDSphere(:,2),1)
     ENDIF



     call memory_free('init_density_local',(real(size(rho_t),DP)+size(grid_rho))*DP)



     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE inichrg_iso
  !-----------------------PARTING-LINE--------------------------
  SUBROUTINE init_MO_Density_iso()
    USE parameters
    USE struct_module
    USE read_module, ONLY:get_MO_coefficient, parse_headline, sign2lm, destroy_MOinit
    USE grid_module, ONLY:rho_calc,grid,dvol,n1,n2,n3
    IMPLICIT NONE
    REAL(DP) :: rho_t(rho_calc%OneDLength)
    INTEGER(I4B)   :: i,j
    REAL(DP)       :: tele !,Noccupy=9
    !REAL(DP)       :: occupy(9)=(/2,2,2,2,2,2,2,1,1/)

    CALL get_MO_coefficient()
    CALL ISO_init_MO(rho_calc%initMO)

    rho_t=0.d0
    DO i=1,struct%Noccupy
      DO j=1,rho_calc%OneDLength
        rho_t(j)=rho_t(j)+rho_calc%initMO(j,i)**2.d0*struct%occupy(i)
      ENDDO
    ENDDO
    !DO i=1,rho_calc%OneDLength,1
    !  rho_t(i)=SUM(rho_calc%initMO(i,1:Noccupy)**2*occupy(:))
      !rho_t(i)=1
    !ENDDO
    tele=SUM(rho_t(:))*dvol
    print*,"tele",tele
    print*,"ncharge",ncharge
    !
    rho_t(:)=rho_t(:)*ncharge/tele
    print*,"aaaaaaaaaaaaa",sum(rho_t)*dvol
    !decide the initial density
    !===============================================================
    !##assignment to grid
    grid%rhoS=0.d0
    IF(NSPIN==1)THEN
      DO i=1,rho_calc%OneDLength
        grid%rhoS(rho_calc%x(i),rho_calc%y(i),rho_calc%z(i),1)= &
            & rho_t(i)
      ENDDO
    ELSE
      DO i=1,rho_calc%OneDLength
        grid%rhoS(rho_calc%x(i),rho_calc%y(i),rho_calc%z(i),1)= &
            & rho_t(i)*0.48d0
        grid%rhoS(rho_calc%x(i),rho_calc%y(i),rho_calc%z(i),1)= &
            & rho_t(i)*0.52d0
      ENDDO
    ENDIF
    CALL destroy_MOinit()
!    rho_calc%volume=volume/REAL(n1*n2*n3,DP)*rho_calc%OneDLength
!    IF(NSPIN==1)THEN
!      rho_calc%OneDSphere(:,:)=ncharge/rho_calc%volume
!    ELSE
!      rho_calc%OneDSphere(:,1)=ncharge/rho_calc%volume*0.48d0
!      rho_calc%OneDSphere(:,2)=ncharge/rho_calc%volume*0.52d0
!    ENDIF
    !===============================================================
    !##assignment to grid
    !DO j=1,NSPIN,1
    !DO i=1,rho_calc%OneDLength
    !  grid%rhoS(rho_calc%x(i),rho_calc%y(i),rho_calc%z(i),j)= &
    !      & rho_calc%OneDSphere(i,j)
    !ENDDO
    !ENDDO
    !===============================================================
    !STOP
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE init_MO_Density_iso
  !-----------------------PARTING-LINE--------------------------
  SUBROUTINE ISO_init_MO(initX_sto)
    !initialize the subsystem by MO
    USE parameters , ONLY : Nstates
    USE struct_module , ONLY : struct,lat_mat,naty
    USE grid_module , ONLY : grid,n,n1,n2,n3,rho_calc,gap
    USE math , ONLY : atom_sto
    IMPLICIT NONE
    REAL(DP),INTENT(INOUT) :: initX_sto(:,:)
    !LOCAL
    INTEGER(I4B) :: Ity,Ia,Ip,icx,icy,icz,Ii &
         & ,nrep=1 &
         &,Il,Im,Nstart
    REAL(DP) :: ra(3),rat(4),f,randt
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !open(22222,file='STO_f')
    !open(22223,file='STO_r')
    initX_sto(:,:)=0.d0
    !all MO states
    DO Ii=1,Nstates,1
    !all type
    DO Ity=1,naty
    !all atom
    DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
      !store the position
      ra(:)=struct%poscar(:,Ia)
      !call the angle momentum
      DO Il=0,struct%Lmax(Ity)
      DO Im=-Il,Il
        !all points
        DO Ip=1,rho_calc%OneDLength
          !replicas
          !DO icz=-nrep,nrep
          !DO icy=-nrep,nrep
          !DO icx=-nrep,nrep
          !  rat(1:3)= icx*lat_mat(:,1)  &
          !          & +  icy*lat_mat(:,2)  &
          !          & +  icz*lat_mat(:,3)  &
          !          & +  grid%rVec(1:3,Ip)  &
          !          & -  ra(:)
            rat(1)=rho_calc%x(Ip)*gap(1)-ra(1)
            rat(2)=rho_calc%y(Ip)*gap(2)-ra(2)
            rat(3)=rho_calc%z(Ip)*gap(3)-ra(3)
            rat(4)=SQRT( rat(1)**2 + rat(2)**2 +rat(3)**2 )
            !>STO
            !print*,'Il',Il,'np~n',struct%prinq(Il+1,Ity)
            !print*,struct.prinq(:,Ity)
            CALL atom_STO(struct%prinq(Il+1,Ity),Il,Im, &
                        & struct%zeta(Il+1,Ity),rat,f)
            !if(Ia==1.and.Il==0.and.Im==0)THEN
            !  if(rat(1)-rat(2).lt.0.0000001d0.and.rat(1)-rat(2).lt.0.0000001d0)THEN
            !    write(22222,*) f
            !    write(22223,*) rat(4)
            !  endif
            !endif
            !>init
            !print*,'MO',Ii,struct.coeff(Ia,Il,Im,Ii),'atom',Ia
            initX_sto(Ip,Ii)=initX_sto(Ip,Ii)+struct%coeff(Ia,Il,Im,Ii)*f
          !ENDDO
          !ENDDO
          !ENDDO
        ENDDO !>grid
      ENDDO !>m
      ENDDO !>l
    ENDDO !>atoms of element "foo"
    ENDDO !>elements
  ENDDO !>MO states
    !open(22222,file='STO_f')
    !open(22223,file='STO_r')
    !   do Ii=1,10000,1
    !        rat(4)=real(Ii,DP)*0.0005
    !        CALL atom_STO(struct%prinq(0+1,1),0,0, &
    !                    & struct%zeta(0+1,1),rat(4),f)
    !         write(22222,*) f
    !         write(22223,*) rat(4)
    !   enddo
    !close(22222)
    !close(22223)
  ENDSUBROUTINE ISO_init_MO
  !-----------------------DIVIDER-LINE--------------------------
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDMODULE Begin_module

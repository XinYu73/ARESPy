MODULE Begin_module
   USE constants
   IMPLICIT NONE
CONTAINS
  !-------------------PARTING LINE--------------------
  SUBROUTINE Initial_grid_pbc()
     USE math ,ONLY : lat2matrix,inv_33,det,dir2car
     USE struct_module , ONLY : lat_mat,recip_lat,lat_para &
                            & , volume,charge_ave,ncharge &
                            & , struct
     USE parameters, ONLY : finite_order,KSPACING,atomRc
     USE Grid_module , ONLY : Build_rgrid,Build_kgrid &
                          &  ,Build_eigen,n1,n2,n3 &
                          & ,nk1,nk2,nk3,nk,gap &
                          &,FillQTable,FillRTable &
#ifdef MPI
                          &, global_n1,global_n2,global_n3 &
                          &, global_n,build_parallel_3d_grid,grid
     USE smpi_math_module, ONLY: parallel,smpi_diff_init_sph,mpinfo
     USE finite_module , ONLY : init_finite,cell_mu
#else
                          &,build_ISO_sphere_grid
     USE finite_module , ONLY : init_finite
#endif
     USE FOURIER, ONLY : FFT,PlanFFT,CleanFFT
     USE nlpot_module, ONLY : initialize_nlpot
     IMPLICIT NONE
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     CALL dir2car(struct%pos,struct%poscar,lat_mat)
     !cell_para
     CALL lat2matrix(lat_para,lat_mat,2)
     !set cell
     recip_lat(:,:)=2.d0*pi*TRANSPOSE(inv_33(lat_mat))
     volume = ABS(det(lat_mat))
     !creat grids
     CALL Build_rgrid()
     !creat k-points
     CALL Build_kgrid()
     !creat eigen-data
     CALL Build_eigen()
     !>>>Finite difference
     CALL init_finite(gap)
#ifdef MPI
      ! call smpi_diff_init(global_n1,global_n2,global_n3,global_n,finite_order,cell_mu)
      call smpi_diff_init_sph(global_n1,global_n2,global_n3,global_n,finite_order,cell_mu,grid%lsp)
#endif
     !<<<Finite difference
     !Plan for FFT
#ifdef MPI
     CALL PlanFFT(global_n1,global_n2,global_n3,grid%lsp)
#else
     CALL PlanFFT(n1,n2,n3)
#endif
     !recip grid mesh data
     CALL FillQTable()
     !real grid mesh data
     CALL FillRTable()
     !initial density
     CALL Initial_density()
!print*,'density',parallel%myid
!CALL MPI_BARRIER(parallel%comm,mpinfo)
     !initial nonlocal pseudopotentials
     CALL initialize_nlpot()
!print*,'initialize',parallel%myid
!CALL MPI_BARRIER(parallel%comm,mpinfo)
!STOP
#ifdef MPI
     if(parallel%isroot)then
        WRITE(6,*)'R-space-GRIDS:',global_n1,global_n2,global_n3
        IF(atomRc > 0.d0) WRITE(6,*)'R-space-sphere:',global_n1*global_n2*global_n3, &
                   & '->',grid%oneDlength
     !WRITE(6,'(A21,3I5,A6,I5)')'local R-space-GRIDS:',n1,n2,n3,',rank:',parallel%myid
        WRITE(6,*)'K-space-GRIDS:',nk1,nk2,nk3
        IF(KSPACING>0.d0)THEN
           WRITE(6,*)'Num of K-used:',nk
        ELSE
           WRITE(6,*)'Only gamma point is used'
        ENDIF
     endif
#else
     WRITE(6,*)'R-space-GRIDS:',n1,n2,n3
     WRITE(6,*)'K-space-GRIDS:',nk1,nk2,nk3
     IF(KSPACING>0.d0)THEN
        WRITE(6,*)'Num of K-used:',nk
     ELSE
        WRITE(6,*)'Only gamma point is used'
     ENDIF
#endif
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE Initial_grid_pbc
  !----------------Initial Charge density---------------------
  SUBROUTINE Initial_density()
     USE parameters , ONLY : LINRHO,LAtomRho,nspin
     USE grid_module , ONLY : sph=>grid,dvol &
           &,ni1=>global_n1,ni2=>global_n2,ni3=>global_n3 &
           &,ni=>global_n,nps=>n,Cubic2Sphere,sumrhoS
     USE struct_module , ONLY : charge_ave,ncharge
#ifdef MPI
     USE smpi_math_module
#endif
     IMPLICIT NONE
     INTEGER(I4B) :: in1,in2,in3,i,isp
     LOGICAL :: lin
     REAL(DP) :: rhoin(ni1,ni2,ni3,nspin)!,rhos(ni1,ni2,ni3)
     REAL(DP) :: tele,tele1,tele2
#ifdef MPI
     REAL(DP) :: tele_local,tele_local1,tele_local2
#endif
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     IF(LINRHO)THEN
#ifdef MPI
        IF(parallel%isroot)THEN
#endif
        rhoin=0._DP
        OPEN(111,FILE='INRHO')
           !DO i=1,16
           !    READ(111,*)
           !ENDDO
           DO isp=1,nspin
              READ(111,*) in1,in2,in3
              lin=(in1==ni1.AND.in2==ni2.AND.in3==ni3)
              IF(.NOT.lin)THEN
                 PRINT*,'Check the INRHO file meshs,should be',ni1,ni2,ni3
                 STOP
              ENDIF
              READ(111,*) rhoin(:,:,:,isp)
           ENDDO
        CLOSE(111)
#ifdef MPI
        ENDIF
        CALL MPI_BCAST(rhoin ,nspin*ni,MPI_REAL8,parallel%rootid,parallel%comm,mpinfo)
#endif
        IF(nspin==1)THEN
          CALL Cubic2Sphere(nps,rhoin(:,:,:,1),sph%rhoS(:,1))
        ELSE !nspin==2
          CALL Cubic2Sphere(nps,rhoin(:,:,:,1)+rhoin(:,:,:,2),sph%rhoS(:,1))
          CALL Cubic2Sphere(nps,rhoin(:,:,:,1)-rhoin(:,:,:,2),sph%rhoS(:,2))
          sph%rhoS(:,:)=sph%rhoS(:,:)/2._DP
        ENDIF

     ELSEIF(LAtomRho)THEN
        CALL inichrg_sp()
     ELSE
#ifdef MPI
        IF(parallel%isroot)THEN
#endif
        PRINT*,'Charge_ave in sphere:',charge_ave
#ifdef MPI
        ENDIF
#endif
        IF(nspin==1)THEN
           sph%rhoS(:,1)=charge_ave
        ELSE
           sph%rhoS(:,1)=0.52*charge_ave
           sph%rhoS(:,2)=0.48*charge_ave
        ENDIF
     ENDIF

#ifdef MPI
     tele_local=SUM(sph%rhoS)*dvol
     CALL MPI_ALLREDUCE(tele_local,tele,1,MPI_REAL8,MPI_SUM,parallel%commx,mpinfo)
     IF(nspin==2)THEN
        tele_local1=SUM(sph%rhoS(:,1))*dvol
        tele_local2=SUM(sph%rhoS(:,2))*dvol
        CALL MPI_ALLREDUCE(tele_local1,tele1,1,MPI_REAL8,MPI_SUM,parallel%commx,mpinfo)
        CALL MPI_ALLREDUCE(tele_local2,tele2,1,MPI_REAL8,MPI_SUM,parallel%commx,mpinfo)
     ENDIF
#else
     tele=SUM(sph%rhoS)*dvol
     IF(nspin==2)THEN
        tele1=SUM(sph%rhoS(:,1))*dvol
        tele2=SUM(sph%rhoS(:,2))*dvol
     ENDIF
#endif
     
#ifdef MPI
     IF(parallel%isroot)THEN
#endif
     PRINT*,'Total Charge# in Sphere is', tele
     IF(nspin==2) PRINT*,'Charge# alpha,beta',tele1,tele2
#ifdef MPI
     ENDIF

#endif
     CALL sumrhoS(nps,sph%rhoS,sph%rho)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE Initial_density
  !------------------inichrg---------------------------
  SUBROUTINE inichrg_sp()
     USE parameters , ONLY : nspin
     USE math , ONLY : interp
     USE grid_module , ONLY : n1,n2,n3,nsp=>n,sph=>grid,dvol
     USE struct_module , ONLY : struct,naty,lat_mat,ncharge
     USE pspot_module , ONLY : psp
#ifdef MPI
     USE smpi_math_module
#endif
     IMPLICIT NONE
     !LOCAL
     REAL(DP) :: rho_t(nsp)
     INTEGER(I4B) :: Ity,Ia,Ip,icx,icy,icz
     REAL(DP) :: rat(3),ra(3),rRnorm  &
              &  , acha , tele
     INTEGER(I4B) :: nrep=1
#ifdef MPI
     REAL(DP) :: tele_local
#endif
     REAL(DP) :: drcut=4*ang2bohr
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     rho_t(:)=0.d0
     !Cycle all species of atoms
     DO Ity=1,naty

        !Cycle all atoms in this species
        DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
           ra(:)=struct%poscar(:,Ia)
           !all points
           IF(nsp<1) CYCLE
           DO ip=1,nsp
              DO icz=-nrep,nrep
              DO icy=-nrep,nrep
              DO icx=-nrep,nrep
                 rat(:)= icx*lat_mat(:,1)  &
                    & +  icy*lat_mat(:,2)  &
                    & +  icz*lat_mat(:,3)  &
                    & +  sph%rVec(1:3,ip)  -  ra(:)
                 rRnorm=SQRT(DOT_PRODUCT(rat,rat))
                 IF(rRnorm>drcut) CYCLE
                 !interpolate density
                 acha=interp(psp(Ity)%numps,psp(Ity)%denr,psp(Ity)%r,rRnorm)

                 rho_t(ip) = rho_t(ip) + acha
              ENDDO
              ENDDO
              ENDDO
           ENDDO
           !
        ENDDO
     ENDDO
     !
#ifndef MPI
     tele=SUM(rho_t(:))*dvol
#else
     tele_local=SUM(rho_t(:))*dvol
     CALL MPI_ALLREDUCE(tele_local,tele,1,MPI_REAL8,MPI_SUM,parallel%commx,mpinfo)
#endif
     !re-scaled
     rho_t(:)=rho_t(:)*(ncharge/tele)
     !decide the initial density
     sph%rho(:)=rho_t(:)
     IF(nspin==1)THEN
        sph%rhoS(:,1)=rho_t(:)
     ELSE
        sph%rhoS(:,1)=rho_t(:)*0.48d0
        sph%rhoS(:,2)=rho_t(:)*0.52d0
     ENDIF
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE inichrg_sp
  !-----------------------DIVIDER-LINE--------------------------
ENDMODULE Begin_module

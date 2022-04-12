# 1 "Energy_module.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "Energy_module.f90"
MODULE energy_module
   USE constants
   USE parameters ,  ONLY : Lpbc
   USE grid_module , ONLY : dvol
   USE ewald , ONLY : ewald_energy,ISO_ewald_energy
   USE struct_module , ONLY : lat_mat,struct
   IMPLICIT NONE
   REAL(DP) :: Etot=0.d0  &
           &,  Ekine=0.d0 &
           &,  Ehart=0.d0 &
           &,  Exc=0.d0   &
           &,  Eband=0.d0 &
           &,  Eie=0.d0   &
           &,  Eienl=0.d0   &
           &,  Eewald=0.d0 &
           &,  Efm=0.d0    &
           &,  Tnad=0.d0   &
           &,  FE=0.d0    &
           &,  FE0=0.d0
   REAL(DP) :: WmaxL= 0.d0 &
           &,  WminL= 0.d0

CONTAINS
   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !###########################################################!
   !k-space calculate                                          !
   !###########################################################!
   !------------------total energy calculation------------------
   SUBROUTINE totalenergy(psi,rhoS,eval,fmenergy,ets,Llast)
      USE parameters , ONLY : NSPIN,Lpbc,Lke_pot
      USE smearing_module , ONLY : wke
      USE grid_module , ONLY : grid    ,n1,n2,n3
      USE potential_module
      USE pspot_module , ONLY : max_nproj
      USE libxc_module
      USE IsolateSet ,   ONLY : IsolateVCoulomb_a

      USE smpi_math_module, ONLY: parallel,mpi_real8,mpi_sum,mpinfo

      use m_time_evaluate, only:filename
      IMPLICIT NONE
      !IN/OUT
      COMPLEX(DP),INTENT(IN) :: psi(:,:,:,:)
      REAL(DP),INTENT(IN) :: rhoS(:,:,:,:)
      REAL(DP),INTENT(IN) :: eval(:,:,:),fmenergy,ets
      LOGICAL,INTENT(IN)  :: Llast
      !LOCAL
      REAL(DP) :: Evxc

      REAL(DP) :: Evxc_local
      INTEGER(I4B) :: ix,iy,iz,i

      REAL(DP),DIMENSION(n1,n2,n3) :: &
               vh  &!hartree
        &  ,  rho !total density
      REAL(DP) :: vxc(n1,n2,n3,NSPIN)
      INTEGER(I4B) :: Is
      REAL(DP),external :: ddot
      !Etot=Eband-Ehart-\int(rho*Vxc)d^3r+E_xc
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !total rho
      rho(:,:,:)=0.d0
      DO Is=1,NSPIN
         rho(:,:,:)=rho(:,:,:)+rhoS(:,:,:,Is)
      ENDDO
      !now
! #ifdef 1
!          ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
!          iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
!          call MPI_ALLGATHERV(rhoS(ix,iy,1,1),parallel%mygrid_range(3)&
!               &,MPI_REAL8,rho_global,parallel%recvcounts&
!               &,parallel%displs,MPI_REAL8&
!               & ,parallel%commx,mpinfo)
         ! open(1111,file="rho"//filename(10:11))
         ! ! do i=1,size(psi,2),1
         ! !    write(1111+parallel%myid,'(10000(F8.2,1X))')real(psi(:,i,1,1))
         ! ! enddo
         ! write(1111,*)rho!_global
         ! close(1111)
         ! print*,'sum rho_global',sum(rho)!_global)
         ! stop
! #endif
      CALL Ebands(eval,wke,Eband)
      !hartree
      CALL vhartree(rho,vh)
      CALL Ehartree(rho,vh,Ehart)
      !exchange correlation
      CALL LDAlib_energy(rhoS,Exc)
      CALL LDAlib_potential(rhoS,vxc)
      ! Evxc=SUM(vxc*rhoS)*dvol



      ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
      iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
      iz=1
      Evxc_local=0.d0
      do i=1,NSPIN
         Evxc_local=Evxc_local+ddot(parallel%mygrid_range(3),vxc(ix,iy,iz,i),1&
              & ,rhoS(ix,iy,iz,i),1)*dvol
      enddo
      call mpi_allreduce(Evxc_local, Evxc, 1, mpi_real8, mpi_sum, parallel%commx, mpinfo)

      !Ion-Ion
      IF(Llast)THEN
         Eewald=ewald_energy(lat_mat,struct%pos,struct%eleid,struct%Zion)
         !ion-e local part
         CALL E_ie(rho,grid%vlpp,Eie)
         !ion-e nonlocal part
         !IF(max_nproj>0)THEN
         !   CALL E_ienl(psi,wke,Eienl)
         !ENDIF
      ELSE
         Eewald=0.d0
         Eie=0.d0
         Eienl=0.d0
      ENDIF
      !total energy
      ! print*,'Eband-Ehart-Evxc+Exc+Eewald',Eband,Ehart,Evxc,Exc,Eewald
      Etot=Eband-Ehart-Evxc+Exc+Eewald
      !fermi energy
      Efm=fmenergy
      !Free energy
      FE=Etot-ets
      FE0=Etot-0.5d0*ets
      ! print*,'Efm,Ets',Efm,ets
      ! stop

      if(Llast.and.Lke_pot)then
         if(parallel%isroot)print*,"out Kenl Potential"
         CALL out_Ke_potential(Efm,grid%vlpp,vh,Vxc(:,:,:,1))
      endif

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE totalenergy

   SUBROUTINE out_Ke_potential(Efm,Vie,Vh,Vxc)
     USE grid_module , ONLY : grid,n1,n2,n3,global_n1&
          & ,global_n2,global_n3,global_n
      USE smpi_math_module, ONLY: parallel,mpi_real8,mpi_sum,mpinfo
     implicit none
     !> in/out
     REAL(DP) :: Efm, Vie(n1,n2,n3),Vh(n1,n2,n3)&
          & ,Vxc(n1,n2,n3)
     !> local
     INTEGER(I4B) :: ix,iy
     REAL(DP) :: Vke(n1,n2,n3)
     REAL(DP) :: Vke_global(global_n1,global_n2,global_n3)

     !> Vke=Efm-Vxc-Vh-Vie
     Vke=Efm-Vie-Vh-Vxc

     !> gather
     ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
     iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
     call MPI_ALLGATHERV(Vke(ix,iy,1),parallel%mygrid_range(3)&
          &,MPI_REAL8,Vke_global,parallel%recvcounts&
          &,parallel%displs,MPI_REAL8&
          & ,parallel%commx,mpinfo)

     !> out
     if(parallel%isroot)then
        open(1698,file='KE_POT')
        write(1698,*)global_n1,global_n2,global_n3
        write(1698,*)Vke_global
        close(1698)
     endif

   ENDSUBROUTINE out_Ke_potential

   !---------------total energy gamma calculation----------------
   SUBROUTINE totalenergy_gamma(psi,rhoS,eval,fmenergy,ets,Llast)
      USE parameters , ONLY : NSPIN,Lpbc
      USE smearing_module , ONLY : wke
      USE grid_module , ONLY : grid    ,n1,n2,n3
      USE potential_module
      USE pspot_module , ONLY : max_nproj
      USE libxc_module
      USE IsolateSet ,   ONLY : IsolateVCoulomb_a

      USE smpi_math_module, ONLY: parallel,mpi_real8,mpi_sum,mpinfo

      use m_time_evaluate, only:filename
      IMPLICIT NONE
      !IN/OUT
      REAL(DP),INTENT(IN) :: psi(:,:,:,:)
      REAL(DP),INTENT(IN) :: rhoS(:,:,:,:)
      REAL(DP),INTENT(IN) :: eval(:,:,:),fmenergy,ets
      LOGICAL,INTENT(IN)  :: Llast
      !LOCAL
      REAL(DP) :: Evxc

      REAL(DP) :: Evxc_local
      INTEGER(I4B) :: ix,iy,iz,i

      REAL(DP),DIMENSION(n1,n2,n3) :: &
               vh  &!hartree
        &  ,  rho !total density
      REAL(DP) :: vxc(n1,n2,n3,NSPIN)
      INTEGER(I4B) :: Is
      REAL(DP),external :: ddot
      !Etot=Eband-Ehart-\int(rho*Vxc)d^3r+E_xc
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !total rho
      rho(:,:,:)=0.d0
      DO Is=1,NSPIN
         rho(:,:,:)=rho(:,:,:)+rhoS(:,:,:,Is)
      ENDDO
      !now
! #ifdef 1
!          ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
!          iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
!          call MPI_ALLGATHERV(rhoS(ix,iy,1,1),parallel%mygrid_range(3)&
!               &,MPI_REAL8,rho_global,parallel%recvcounts&
!               &,parallel%displs,MPI_REAL8&
!               & ,parallel%commx,mpinfo)
         ! open(1111,file="rho"//filename(10:11))
         ! ! do i=1,size(psi,2),1
         ! !    write(1111+parallel%myid,'(10000(F8.2,1X))')real(psi(:,i,1,1))
         ! ! enddo
         ! write(1111,*)rho!_global
         ! close(1111)
         ! print*,'sum rho_global',sum(rho)!_global)
         ! stop
! #endif
      CALL Ebands(eval,wke,Eband)
      !hartree
      CALL vhartree(rho,vh)
      CALL Ehartree(rho,vh,Ehart)
      !exchange correlation
      CALL LDAlib_energy(rhoS,Exc)
      CALL LDAlib_potential(rhoS,vxc)
      ! Evxc=SUM(vxc*rhoS)*dvol



      ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
      iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
      iz=1
      Evxc_local=0.d0
      do i=1,NSPIN
         Evxc_local=Evxc_local+ddot(parallel%mygrid_range(3),vxc(ix,iy,iz,i),1&
              & ,rhoS(ix,iy,iz,i),1)*dvol
      enddo
      call mpi_allreduce(Evxc_local, Evxc, 1, mpi_real8, mpi_sum, parallel%commx, mpinfo)

      !Ion-Ion
      IF(Llast)THEN
         Eewald=ewald_energy(lat_mat,struct%pos,struct%eleid,struct%Zion)
         !ion-e local part
         CALL E_ie(rho,grid%vlpp,Eie)
         !ion-e nonlocal part
         !IF(max_nproj>0)THEN
         !   CALL E_ienl(psi,wke,Eienl)
         !ENDIF
      ELSE
         Eewald=0.d0
         Eie=0.d0
         Eienl=0.d0
      ENDIF
      !total energy
      ! print*,'Eband-Ehart-Evxc+Exc+Eewald',Eband,Ehart,Evxc,Exc,Eewald
      Etot=Eband-Ehart-Evxc+Exc+Eewald
      !fermi energy
      Efm=fmenergy
      !Free energy
      FE=Etot-ets
      FE0=Etot-0.5d0*ets
      ! print*,'Efm,Ets',Efm,ets
      ! stop
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ENDSUBROUTINE totalenergy_gamma
   !---------------------hartree energy-------------------------
   SUBROUTINE Ehartree(rho,vhart,Ehart)

      USE smpi_math_module, ONLY: parallel,mpi_real8,mpi_sum,mpinfo
      USE grid_module , ONLY : n1,n2,n3

      !
      IMPLICIT NONE
      REAL(DP),INTENT(OUT)      :: Ehart
      real(DP),dimension(:,:,:),INTENT(IN) :: rho,vhart
      real(DP),external :: ddot

      INTEGER(I4B) :: ix,iy,iz,i
      REAL(DP)     :: Ehart_local

      !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>



      ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
      iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
      iz=1
      Ehart_local=0.d0
      Ehart_local=0.5d0 * ddot(parallel%mygrid_range(3),vhart(ix,iy,iz),1&
              & ,rho(ix,iy,iz),1)*dvol
      call mpi_allreduce(Ehart_local, Ehart, 1, mpi_real8, mpi_sum, parallel%commx, mpinfo)

   ENDSUBROUTINE Ehartree
   !---------------------hartree energy-------------------------
   SUBROUTINE Ehartree_iso(rho,vhart,Ehart)

     USE smpi_math_module, ONLY: parallel,mpi_real8,mpi_sum,mpinfo

      !
      IMPLICIT NONE
      REAL(DP),INTENT(OUT)      :: Ehart



      real(DP),dimension(:),INTENT(IN) :: rho,vhart
      real(DP) :: Ehart_local

      real(DP),external :: ddot
      !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>



      Ehart_local=0.5d0 * ddot(size(vhart),vhart,1,rho,1)*dvol
      call mpi_allreduce(Ehart_local, Ehart, 1, mpi_real8, mpi_sum, parallel%comm, mpinfo)

      !<<<<<<<<<<<<<<<<<<<<< End Body  <<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE Ehartree_iso
   !----------------------ion-e energy--------------------------
   SUBROUTINE E_ie(rho,vion,Eie)

     USE smpi_math_module, ONLY : parallel,mpi_real8,mpi_sum,mpinfo
     USE grid_module, ONLY: n1,n2

     !
     IMPLICIT NONE
     REAL(DP),INTENT(IN) :: rho(:,:,:),vion(:,:,:)
     REAL(DP),INTENT(OUT) :: Eie
     REAL(DP),external :: ddot

     REAL(DP) :: Eie_local
     INTEGER(I4B) :: ix,iy,iz

     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



     ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
     iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
     iz=1
     Eie_local=0.d0
     Eie_local=ddot(parallel%mygrid_range(3),rho(ix,iy,iz),1,vion(ix,iy,iz),1)*dvol
     call mpi_allreduce(Eie_local, Eie, 1, mpi_real8, mpi_sum, parallel%commx, mpinfo)

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE
   !---------------------For  band energy-----------------------
   SUBROUTINE Ebands(eval,wke,Eband)
      !
      USE parameters , ONLY : NSPIN
      USE grid_module , ONLY : nk
      IMPLICIT NONE
      REAL(DP),INTENT(IN) :: eval(:,:,:)
      REAL(DP),INTENT(IN) :: wke(:,:,:)
      REAL(DP),INTENT(OUT) :: Eband
      !
      INTEGER(I4B) :: nev,Ioc,Ik,Ispin
      REAL(DP) :: Sn
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      Eband=0.d0
      nev=SIZE(eval,1)

      DO Ispin=1,NSPIN
         !integral in FBZ
         DO Ik=1,nk
            !every band
            DO Ioc=1,nev
                  !CALL smearSN(eval(Ioc,k1,k2,k3),fme,Wsmear,Sn)
                  Eband=Eband + eval(Ioc,Ik,Ispin)*wke(Ioc,Ik,Ispin)

            ENDDO
         ENDDO
      ENDDO
      !test last band
      WmaxL=MAXVAL(wke(nev,:,:))
      WminL=MINVAL(wke(nev,:,:))
      !>>>>>>>>>
      !print*,wke
      !STOP
      !<<<<<<<<<
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE Ebands
   !###########################################################!
   !For BvK cell calculate                                     !
   !###########################################################!
   !------------------total energy calculation------------------
   SUBROUTINE BvK_totalenergy(rhoS,eval,fmenergy,ets,Llast)
      USE parameters , ONLY : Nspin,Nstates
      USE grid_module , ONLY : grid , n1,n2,n3
      USE potential_module
      USE smearing_module , ONLY : wke
      USE libxc_module
      USE IsolateSet ,      ONLY : IsolateVCoulomb_a
      IMPLICIT NONE
      REAL(DP),INTENT(IN) :: rhoS(:,:,:,:)
      REAL(DP),INTENT(IN) :: eval(:,:),fmenergy,ets
      LOGICAL,INTENT(IN)  :: Llast
      !
      REAL(DP) :: Evxc
      REAL(DP),DIMENSION(n1,n2,n3) :: rho & !total density
                              &, vh  !hartree potential
      REAL(DP) :: vxc(n1,n2,n3,Nspin)
      INTEGER(I4B) :: Is
      !Etot=Eband-Ehart-\int(rho*Vxc)d^3r+E_xc
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !total rho
      rho(:,:,:)=0.d0
      DO Is=1,NSPIN
         rho(:,:,:)=rho(:,:,:)+rhoS(:,:,:,Is)
      ENDDO
      !calculate the band energy
      CALL BvK_Ebands(eval,wke(:,1,:),Eband)
      !hartree
      CALL vhartree(rho,vh)
      CALL Ehartree(rho,vh,Ehart)
      !exchange correlation
      CALL LDAlib_energy(rhoS,Exc)
      CALL LDAlib_potential(rhoS,vxc)
      Evxc=SUM(vxc*rhoS)*dvol
      !Ion-Ion
      IF(Llast)THEN
         Eewald=ewald_energy(lat_mat,struct%pos,struct%eleid,struct%Zion)
         !ion-e local part
         CALL E_ie(rho,grid%vlpp,Eie)
      ELSE
         Eewald=0.d0
         Eie=0.d0
         Eienl=0.d0
      ENDIF
      !total energy
      Etot=Eband-Ehart-Evxc+Exc+Eewald
      !fermi energy
      Efm=fmenergy
      !Free energy
      FE=Etot-ets
      FE0=Etot-0.5d0*ets
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE BvK_totalenergy
   !---------------BvK band energy calculation--------------
   SUBROUTINE BvK_Ebands(eval,focc,Eband)
      USE parameters , ONLY : Nspin,nev=>Nstates
      !
      IMPLICIT NONE
      !
      REAL(DP),INTENT(IN) :: eval(:,:)
      REAL(DP),INTENT(IN) :: focc(:,:)
      REAL(DP),INTENT(OUT) :: Eband
      !
      INTEGER(I4B) :: Ioc,Is
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      Eband=0.d0
      !
      DO Is=1,Nspin
         DO Ioc=1,nev
            !2.d0  double occ
            Eband=Eband + focc(Ioc,Is)*eval(Ioc,Is)
         ENDDO
      ENDDO
      !
      !test last band
      WmaxL=MAXVAL(focc(nev,:))
      WminL=MINVAL(focc(nev,:))
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE BvK_Ebands
   !###########################################################!
   !pRR For BvK cell calculate                                 !
   !###########################################################!
   !---------------pRR band energy calculation--------------
   SUBROUTINE pRR_Ebands(veff,Phi,Pbar,Eband)
      USE parameters    , ONLY :NSPIN,Nstates
      USE Grid_module , ONLY : n1,n2,n3,n
      USE chebyshev_module , ONLY : cal_HX_real
      USE Lapack_module , ONLY : matmat_real
      IMPLICIT NONE
      !IN/OUT
      REAL(DP),INTENT(IN)  :: veff(:,:,:,:) & !effective potential
                         &,   Phi(:,:,:)    & !subspace
                         &,   Pbar(:,:,:)     !project density matrix
      REAL(DP),INTENT(OUT) :: Eband
      !LOCAL
      REAL(DP) :: HPhi(n,Nstates,Nspin),PbarPhi(Nstates,n,Nspin)
      REAL(DP) :: nei
      INTEGER(I4B) :: Is,I
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      DO Is=1,Nspin
         !HPhi
         CALL cal_HX_real(veff(:,:,:,Is),Nstates,Phi(:,:,Is),HPhi(:,:,Is))
         !Pbar Phi
         CALL matmat_real(Pbar(:,:,Is),Phi(:,:,Is),'N','T',PbarPhi(:,:,Is))
      ENDDO
      !
      Eband=0.d0
      DO Is=1,Nspin
         DO I=1,n
            nei=DOT_PRODUCT(HPhi(I,:,Is),PbarPhi(:,I,Is))
            Eband=Eband+nei
         ENDDO
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE pRR_Ebands
   !------------------total energy calculation------------------
   SUBROUTINE pRR_totalenergy(rhoS,Phi,Pbar,fmenergy,ets,Llast)
      USE parameters , ONLY : Nspin,Nstates
      USE grid_module , ONLY : grid , n1,n2,n3
      USE potential_module
      USE smearing_module , ONLY : wke
      USE libxc_module
      USE IsolateSet ,      ONLY : IsolateVCoulomb_a
      IMPLICIT NONE
      REAL(DP),INTENT(IN) :: rhoS(:,:,:,:)
      REAL(DP),INTENT(IN) :: Phi(:,:,:)  &
                       &  ,  Pbar(:,:,:)
      REAL(DP),INTENT(IN) :: fmenergy,ets
      LOGICAL,INTENT(IN)  :: Llast
      !
      REAL(DP) :: Evxc
      REAL(DP),DIMENSION(n1,n2,n3) :: rho & !total density
                              &, vh  !hartree potential
      REAL(DP) :: vtmp(n1,n2,n3,Nspin)
      INTEGER(I4B) :: Is
      !Etot=Eband-Ehart-\int(rho*Vxc)d^3r+E_xc
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !total rho
      rho(:,:,:)=0.d0
      DO Is=1,NSPIN
         rho(:,:,:)=rho(:,:,:)+rhoS(:,:,:,Is)
      ENDDO
      !hartree
      CALL vhartree(rho,vh)
      CALL Ehartree(rho,vh,Ehart)
      !exchange correlation
      CALL LDAlib_energy(rhoS,Exc)
      CALL LDAlib_potential(rhoS,vtmp)
      Evxc=SUM(vtmp*rhoS)*dvol
      !Ion-Ion
      IF(Llast)THEN
         Eewald=ewald_energy(lat_mat,struct%pos,struct%eleid,struct%Zion)
         !ion-e local part
         CALL E_ie(rho,grid%vlpp,Eie)
      ELSE
         Eewald=0.d0
         Eie=0.d0
         Eienl=0.d0
      ENDIF
      !veff=vh+vlpp+vxc
      !CALL cal_veff(rhoS,vtmp)
      DO Is=1,Nspin
         vtmp(:,:,:,Is)=vtmp(:,:,:,Is)+vh(:,:,:)+grid%vlpp(:,:,:)
      ENDDO
      !calculate the band energy
      CALL pRR_Ebands(vtmp,Phi,Pbar,Eband)

      !total energy
      Etot=Eband-Ehart-Evxc+Exc+Eewald


      !fermi energy
      Efm=fmenergy
      !Free energy
      FE=Etot-ets
      FE0=Etot-0.5d0*ets
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ENDSUBROUTINE pRR_totalenergy
   !###########################################################!
   !ISO For cubic cell calculate                                 !
   !###########################################################!
   !------------------total energy calculation(ISO)-------------
   SUBROUTINE ISO_totalenergy(rhoS,eval,fmenergy,ets,Llast)
     !###########################################################!
     !For Isolate cell calculate                                 !
     !###########################################################!
     !------------------total energy calculation------------------
      USE parameters , ONLY : NSPIN,Lpbc
      USE smearing_module , ONLY : wke
      USE grid_module , ONLY : grid    ,n1,n2,n3,rho_calc
      USE potential_module
      USE libxc_module
      USE IsolateSet ,   ONLY : IsolateVCoulomb_a

      USE smpi_math_module, ONLY: parallel,mpi_real8,mpi_sum,mpinfo

      USE m_time_evaluate, ONLY: memory_sum,memory_free
      IMPLICIT NONE
      !IN/OUT



      REAL(DP),INTENT(IN) :: rhoS(:,:)

      REAL(DP),INTENT(IN) :: eval(:,:),fmenergy,ets
      LOGICAL,INTENT(IN)  :: Llast
      !LOCAL
      REAL(DP) :: Evxc







      REAL(DP) :: Evxc_local
      REAL(DP),DIMENSION(parallel%mygrid_range(3)) :: &
               vh  &!hartree
        &  ,  rho  !total density
      REAL(DP) :: vxc(parallel%mygrid_range(3),NSPIN)
      REAL(DP) :: vxc_ii(parallel%mygrid_range(3))

      INTEGER(I4B) :: Is,t1,t2,i,kk
      REAL(DP),external :: ddot
      !Etot=Eband-Ehart-\t(rho*Vxc)d^3r+E_xc
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      call memory_sum('iso_totalenergy_local',real(size(vh),DP)*2*DP+(size(vxc)+size(vxc_ii))*DP)
      !total rho
      rho=0.d0
      DO Is=1,NSPIN



         rho = rho + rhoS(:,Is)

      ENDDO
      !now
      CALL ISO_Ebands(eval,wke(:,1,:),Eband)
      !exchange correlation
      CALL LDAlib_energy_iso(rhoS,Exc)
      CALL LDAlib_potential_iso(rhoS,vxc)
      ! Evxc=SUM(vxc*rhoS)*dvol



      Evxc_local=ddot(size(rhoS),vxc,1,rhoS,1)*dvol
      call mpi_allreduce(Evxc_local, Evxc, 1, mpi_real8, mpi_sum, parallel%comm, mpinfo)

      !==================================================
      !##ANOTHER Evxc CALCULATE FUNCTION FROM ATALAS
      !Evxc=lda_energy(rhoS(:,:,:,1))
      !====================================================
      !##OUT Exc
      !open(unit=1118,file="Exc_ISO")
      !write(1118,*),size(rhoS,1),size(rhoS,2),size(rhoS,3)
      !write(1118,*)vxc*rhoS*dvol
      !close(1118)
      !====================================================
      !##Ion-Ion and LOCAL Vion ENERGY
      IF(Llast)THEN
        !================================
        !print*,"struct%Zion",struct%Zion
        Eewald=ISO_ewald_energy(lat_mat,struct%poscar,struct%eleid,struct%Zion)
        !print*,'Eewald',Eewald
        !================================
        !ion-e local part



        CALL ISO_E_ie(rho,rho_calc%vlpp,Eie)

        !print*,"lat_mat",lat_mat
        !ion-e nonlocal part
        !IF(max_nproj>0)THEN
        !   CALL E_ienl(psi,wke,Eienl)
        !ENDIF
        !print*,"SUM_pho",sum(rho)*dvol
        !================================
      ELSE
        !==================
        !IN SCF, NOT OUTPUT
        Eewald=0.d0
        Eie=0.d0
        Eienl=0.d0
        !==================
      ENDIF
      !========================================
      !##CAL Hartree ENERGY
      ! CALL system_clock(t1)
      CALL IsolateVCoulomb_a(rho,vh)
      ! CALL system_clock(t2)
      !print*,'Vhartree time',(t2-t1)/10000.d0
      CALL Ehartree_iso(rho,vh,Ehart)
      !========================================
      !>CORRECT FOR SCF (TERM3)
      ! IF(Llast)THEN
      !   IF(ALLOCATED(V_hxc_new)) DEALLOCATE(V_hxc_new)
      !   ALLOCATE(V_hxc_new(n1,n2,n3))
      !   vxc_ii=0.d0
      !   DO i=1,NSPIN,1
      !     vxc_ii=vxc_ii+vxc(:,:,:,i)/2.d0
      !   ENDDO
      !   V_hxc_new=vxc_ii/REAL(NSPIN,DP)+vh
      ! ELSE
      !   IF(ALLOCATED(V_hxc_new))THEN
      !     IF(ALLOCATED(V_hxc_old)) DEALLOCATE(V_hxc_old)
      !     ALLOCATE(V_hxc_old(n1,n2,n3))
      !     V_hxc_old=V_hxc_new
      !     vxc_ii=0.d0
      !     DO i=1,NSPIN,1
      !       vxc_ii=vxc_ii+vxc(:,:,:,i)/2.d0
      !     ENDDO
      !     V_hxc_new=vxc_ii/REAL(NSPIN,DP)+vh
      !   ELSE
      !     ALLOCATE(V_hxc_new(n1,n2,n3))
      !   ENDIF
      ! ENDIF
      !========================================
      !##ACCELERATE BY Vhartree
      V_accelerate=vh
      !print*,'set accelerate ok'
      ACCELERATE=.TRUE.
      !========================================
      !total energy
      Etot=Eband-Ehart-Evxc+Exc+Eewald
      !fermi energy
      Efm=fmenergy
      !Free energy
      FE=Etot-ets
      FE0=Etot-0.5d0*ets

      ! print*,'riguai','Eband',Eband,'Ehart',Ehart,'Evxc',Evxc,'Exc',Exc,'Eewald',Eewald
      call memory_free('iso_totalenergy_local',real(size(vh),DP)*2*DP+(size(vxc)+size(vxc_ii))*DP)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE ISO_totalenergy
   !-----------------------DIVIDER-LINE--------------------------
   !-----------------------PARTING-LINE--------------------------
   SUBROUTINE ISO_Ebands(eval,focc,Eband)
      USE parameters , ONLY : Nspin,nev=>Nstates
      !
      IMPLICIT NONE
      !
      REAL(DP),INTENT(IN) :: eval(:,:)
      REAL(DP),INTENT(IN) :: focc(:,1:)
      REAL(DP),INTENT(OUT) :: Eband
      !
      INTEGER(I4B) :: Ioc,Is
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      Eband=0.d0
      !
      DO Is=1,Nspin
         DO Ioc=1,nev
            !2.d0  double occ
            Eband=Eband + focc(Ioc,Is)*eval(Ioc,Is)
         ENDDO
      ENDDO
      !
      !test last band
      WmaxL=MAXVAL(focc(nev,:))
      WminL=MINVAL(focc(nev,:))
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE ISO_Ebands
   !-----------------------DIVIDER-LINE--------------------------
   SUBROUTINE ISO_E_ie(rho,Vionlpp,Eie)
     USE matvec_module , ONLY : nlocmatvec_r  !##cal nonlocal potential
     USE grid_module   , ONLY : n1,n2,n3,n
     USE parameters    , ONLY : NSPIN

     USE smpi_math_module, ONLY : parallel,mpi_real8,mpi_sum,mpinfo

     !
     IMPLICIT NONE



     REAL(DP),INTENT(IN) :: rho(:),vionlpp(:)  !##sum rho,Vion_local

     REAL(DP),INTENT(OUT) :: Eie
     !##Local
     INTEGER(I4B)         :: i,j,k
     !REAL(DP)             :: rho1D(n) &
     !     & ,vionnlp1D(n) &
     !     & ,vionnlp(n1,n2,n3)
     REAL(DP),external    :: ddot

     REAL(DP) :: Eie_local

     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     ! rho1D=(/rho(:,:,:)/)
     ! vionnlp1D=0.d0
     ! !##calculate the nonlocal E_ie
     ! ! CALL nlocmatvec_r(rho1D,vionnlp1D)
     ! print*,"sum(rho_circle)",SUM(rho1D)*dvol
     ! Eie=SUM(rho*Vionlpp)*dvol



     Eie_local=ddot(size(rho),rho,1,Vionlpp,1)*dvol
     call mpi_allreduce(Eie_local, Eie, 1, mpi_real8, mpi_sum, parallel%comm, mpinfo)

     !print*,"Eie",Eie*rydberg
     !print*,"dvol",dvol
     !print*,"nonlocal",SUM(vionnlp1D)*dvol
     !open(1118,file="Vionnlp")
     !write(1118,*)size(Vionnlp,1),size(Vionnlp,2),size(Vionnlp,3)
     !write(1118,*)vionnlp
     !close(1118)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE ISO_E_ie
   !-----------------------PARTING-LINE--------------------------
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   FUNCTION lda_energy(rhoReal)!{{{
      USE constants
      IMPLICIT NONE

      REAL(kind=DP) :: &
         LDA_Energy             ! The XC energy.

      REAL(kind=DP), DIMENSION(:,:,:), INTENT(IN) :: &
         rhoReal               ! Electron density in real space, spin DEPENDENT

      !>> INTERNAL VARIABLES <<!
      REAL(kind=DP), PARAMETER :: &
         ft = 4.0_DP / 3.0_DP, &              ! four thirds
         ot = 1.0_DP/3.0_DP, &                ! one third
         mot = -1.0_DP/3.0_DP, &              ! negative (minus) one third
         cX = -0.73855876638202234_DP, &      ! -0.75_DP * (3._DP/pi)**ot
         cC = 0.62035049089940009_DP          ! (7.5E-1_DP / pi)**ot

      REAL(kind=DP) :: &
         exch, &                              ! Running exchange energy
         eC, &                                ! Correlation energy coefficient
         rs, &                                ! (3/(4.pi.rho))^1/3
         corr                                 ! Running correlation energy

      REAL(kind=DP), DIMENSION(2), PARAMETER :: &
         a = (3.11E-2_DP,1.555E-2_DP), &      ! Parameter A for rs < 1 in [1]
         b = (-4.8E-2_DP,-2.69E-2_DP), &      ! Parameter B for rs < 1 in [1]
         c = (2E-3_DP,7E-4_DP), &             ! Parameter C (Ceperley-Alder) [1]
         d = (-1.16E-2_DP,-4.8E-3_DP), &      ! Parameter D (C-A) from [1]
         g = (-1.423E-1_DP,-8.43E-2_DP), &    ! Parameter gamma for rs > 1 [1]
         b1 = (1.0529_DP,1.3981_DP), &        ! Parameter beta1 for rs > 1 [1]
         b2 = (3.334E-1_DP,2.611E-1_DP)       ! Parameter beta2 for rs > 1 [1]

      INTEGER(I4B) :: &
         ix, iy, iz                           ! Dummy counters
      REAL(DP)                      :: rhoft(size(rhoReal,1),size(rhoReal,2),size(rhoReal,3))
      !>>>>>>>>>>>>>>>>>>>>> Main Body >>>>>>>>>>>>>>>>>>>>>>>>
      corr = 0.d0
      ! Exchange energy.
      exch = cX * SUM(rhoReal**ft)
      ! Correlation energy.
      !$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:corr) PRIVATE(iz,iy,ix,rs,eC)
      DO iz=1, SIZE(rhoReal,3)
         DO iy=1, SIZE(rhoReal,2)
            DO ix=1, SIZE(rhoReal,1)

               ! Calculate r subscript s [1] in atomic units.
               rs = cC * rhoReal(ix,iy,iz)**mot
               IF (rs<1) THEN
                  eC = LOG(rs)*(a(1) + c(1) * rs) + b(1) + d(1) * rs
               ELSE
                  eC = g(1) / (1 + b1(1) * SQRT(rs) + b2(1) * rs)
               END IF
               corr = corr + rhoReal(ix,iy,iz) * eC
            END DO
         END DO
      END DO
      !$OMP END PARALLEL DO
      LDA_Energy= exch + corr
   END function!}}}
ENDMODULE energy_module

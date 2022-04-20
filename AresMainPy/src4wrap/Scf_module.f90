MODULE scf_module
!############################################################!
!*For :self consistent                                       !
!*Author : Qiang Xu                                          !
!*Date   : 2017-7-20                                         !
!############################################################!
   USE constants
#ifdef MPI
   USE smpi_math_module, ONLY:parallel
#endif
   IMPLICIT NONE
CONTAINS
   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!############################################################!
!*For :k-space self consistent                               !
!*Author : Qiang Xu                                          !
!*Date   : 2017-07-20                                        !
!############################################################!
   SUBROUTINE ArpackSCF(rhoS,psi,eval)
      USE parameters , ONLY : NMITER,RTOL,ETOL,Nstates,NSPIN,MALPHA
      USE grid_module , ONLY : n1,n2,n3,n,nsn,nk,dvol
      USE mixer_module
      USE energy_module
      USE Struct_module , ONLY : ncharge,natom
      USE smearing_module , ONLY : wke,smear_init,destroy_smear,Fermilevel

      IMPLICIT NONE
      !IN/OUT
      REAL(DP),INTENT(INOUT)  :: rhoS(:,:,:,:)
      COMPLEX(DP),INTENT(OUT) :: psi(:,:,:,:)
      REAL(DP),INTENT(OUT)    :: eval(:,:,:)
      !LOCAL
      REAL(DP) :: rhoSd(n1,n2,n3,NSPIN)
      REAL(DP) :: fme,ne,ets
      REAL(DP) :: drho,dtoten,toten,totend,res
      INTEGER(I4B) :: Iter,nev
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ifndef MPI
      nev=Nstates
      ne=REAL(ncharge,8)
      !inirial smear
      CALL smear_init(nev)
      !initial mixer
      CALL init_mixer_data()
      !++++++++++++++++++++++++++++++First diag+++++++++++++++++++++++++++++++
      rhoSd(:,:,:,:)=rhoS(:,:,:,:)
      !First diag
      CALL solver_spin(rhoS,psi,eval,nev,100.d0)
      !update rhoS
      CALL smear_updaterho(nev,ne,psi,eval,wke,fme,ets,rhoS)
      CALL totalenergy(psi,rhoS,eval,fme,ets,.FALSE.)
      toten=Etot
      !first mixing by simple mixing
      rhoS=rhoSd+MALPHA*rhoS
      !cycle for SCF
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO Iter=1,NMITER

         !store old data
         rhoSd(:,:,:,:)=rhoS(:,:,:,:)
         totend=toten

         !diag H
         CALL solver_spin(rhoS,psi,eval,nev,5.0D-5)

         !update rhoS
         CALL smear_updaterho(nev,ne,psi,eval,wke,fme,ets,rhoS)
         CALL totalenergy(psi,rhoS,eval,fme,ets,.FALSE.)
         toten=Etot

         !check for exit
         drho=SUM(ABS(rhoS-rhoSd))/nsn
         dtoten=(toten-totend)*hart2ev/natom
         !===============================================================
         WRITE(6,"(1X,A9,I4,1X,A6,E15.7,1X,A4,E15.7)") &
              &      '>ArpkSCF:',iter,'dTOTEN',dtoten,      &
              &  'dRHO',drho
         !===============================================================
         IF(drho<RTOL.AND.ABS(dtoten)<ETOL) EXIT
         !mixing for update phi and rho
         CALL mixing(Iter,rhoS,rhoSd,res)
      ENDDO
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL destroy_mixer()
      IF(Iter==NMITER) THEN
         WRITE(6,*) 'SCF:failed STOP'
         STOP
      ENDIF

      WRITE(6,*)'SCF DONE!STEP is:',Iter
      !call total energy
      CALL totalenergy(psi,rhoS,eval,fme,ets,.TRUE.)
#else
      print *,"only support CheFSI method in confined system"
      psi=0.d0
      eval=0.d0
#endif
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE ArpackSCF
   !-----------------------Ispin solver-----------------------
   SUBROUTINE solver_spin(rhoS,psi,eval,nev,diagTOL)
      !
      USE grid_module , ONLY : n1,n2,n3
      USE parameters , ONLY : NSPIN
      USE potential_module , ONLY : cal_veff
      IMPLICIT NONE
      REAL(DP),INTENT(IN) :: rhoS(:,:,:,:)
      REAL(DP),INTENT(IN) :: diagTOL !for diag
      COMPLEX(DP),INTENT(OUT) :: psi(:,:,:,:)
      REAL(DP),INTENT(OUT) :: eval(:,:,:)
      INTEGER,INTENT(IN) :: nev
      !LOCAL
      INTEGER(I4B) :: Is
      REAL(DP) :: veff(n1,n2,n3,NSPIN)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      DO Is=1,NSPIN
         CALL cal_veff(rhoS,veff)
         ! print *,"only support isolated system"
         CALL ksolver(veff(:,:,:,Is),psi(:,:,:,Is),eval(:,:,Is),nev,diagTOL)
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE solver_spin
   !----------------------k-space solver----------------------
   SUBROUTINE ksolver(veff,kpsi,keval,nev,TOL)
      USE Grid_module , ONLY : n,n1,n2,n3,nk,KPT
      USE Arpack_module , ONLY : diagH_arpack
      USE math , ONLY : csort_eigen
      IMPLICIT NONE
      REAL(DP),INTENT(IN) :: veff(:,:,:)
      REAL(DP),INTENT(IN) :: TOL  !for diag
      COMPLEX(DP),INTENT(OUT) :: kpsi(:,:,:)
      REAL(DP),INTENT(OUT) :: keval(:,:)
      INTEGER(I4B),INTENT(IN) :: nev  !,Ispin
      !
      INTEGER(I4B) :: Ik,nec,isc,info
      INTEGER(I4B) :: maxmvs
      COMPLEX(DP)  :: resid_int(n)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !center uniform sampling
      !!$OMP PARALLEL DO PRIVATE(Ik,Fk,info,resid_int,isc,maxmvs)
      info=0
      resid_int=(0.d0,0.d0)
      DO Ik=1,nk
         !resid_int=(0.d0,0.d0)
         diagH:DO isc=1,5
                  maxmvs=30000*isc
                  CALL diagH_arpack(  veff, Ik ,nev,kpsi(:,:,Ik),keval(:,Ik), &
                     &     resid_int,nec,info,maxmvs,TOL)
                  IF(nec>=nev)THEN
                     EXIT diagH
                  ELSE
                     WRITE(6,*)'diag(H) failed,restarting'
                     !STOP          !I don't want to do this
                     !kinfo(k1,k2,k3)=0
                  ENDIF
              ENDDO diagH
         IF(isc>=5)THEN
            WRITE(6,*) 'diag(H) failed,STOP,increase NADST and try again'
            STOP
         ENDIF
         !sort eigenvalue and eigenstates
         CALL csort_eigen(nev,keval(:,Ik),kpsi(:,:,Ik))
      ENDDO
!print*,'ksolver'
!print*,keval(:,:)
!STOP
      !!$OMP END PARALLEL DO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE ksolver
   !--------------------smear update  rho---------------------
!#ifndef MPI
   SUBROUTINE smear_updaterho(nev,ne,psi,eval,wke_l,fme,ets,rhoS)
      USE smearing_module
      USE parameters , ONLY : NSPIN,Nstates_global
      USE grid_module , ONLY : n,nsn,dvol,nk,KPT
#ifdef MPI
      USE smpi_math_module, ONLY: MPI_REAL8, MPI_SUM, mpinfo
      use m_time_evaluate, only:filename
#endif
      !
      IMPLICIT NONE
      !INPUT
      INTEGER(I4B),INTENT(IN) :: nev   !number of states
      REAL(DP),INTENT(IN)     :: ne !number of charge
      COMPLEX(DP),INTENT(IN) :: psi(:,:,:,:) !eigen-states
      REAL(DP),INTENT(IN)     :: eval(:,:,:) !eigen-values
      REAL(DP),INTENT(OUT)     :: wke_l(:,:,:)  !weight of k
      REAL(DP),INTENT(OUT) :: ets !energy of electronic entropy
      !OUT PUT
      REAL(DP),INTENT(OUT) :: fme !fermi leval
      REAL(DP),INTENT(OUT) :: rhoS(nsn)
#ifdef MPI
      REAL(DP) ::rhoS_local(nsn)
#endif
      !
      INTEGER(I4B) :: Ik,Ispin,Id,lft,rit !for index
      INTEGER(I4B) :: Iocc !for occupation state index
      REAL(DP) :: prho(n) ! 1D density in each k
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !calculate the weight of k points
#ifndef MPI
      CALL Fermilevel(ne,nev,nk,KPT%wk,eval,wke_l,fme,ets)
#else
      CALL Fermilevel(ne,Nstates_global,nk,KPT%wk,eval,wke_l,fme,ets)
      if(parallel%isroot)then
         ! print*,'ne',ne,'nev',nev
         ! print*,'eval',eval
         ! print*,'fme',fme
         ! print*,'ets',ets
      endif
#endif
      !
      DO Ispin=1,NSPIN
         prho(:)=0.d0
         DO Ik=1,nk
            DO Iocc=1,nev
               !occ number
               IF(ABS(wke_l(parallel%sub2sum(Iocc,parallel%ranky+1),Ik,Ispin))<1e-14) CYCLE
               !
               DO Id=1,n
#ifndef MPI
                  prho(Id)=prho(Id)+ABS(psi(Id,Iocc,Ik,Ispin))**2*wke_l(Iocc,Ik,Ispin)
#else
                  prho(Id)=prho(Id)+ABS(psi(Id,Iocc,Ik,Ispin))**2&
                       &*wke_l(parallel%sub2sum(Iocc,parallel%ranky+1),Ik,Ispin)
#endif
               ENDDO
            ENDDO
         ENDDO
         !
         lft=(Ispin-1)*n+1
         rit=Ispin*n
         !
         rhoS(lft:rit)=prho(:)
      ENDDO
      !dvol for psi normal
#ifdef MPI
       ! rhoS=rhoS/dvol
       rhoS_local=rhoS/dvol
       CALL MPI_ALLREDUCE(rhoS_local, rhoS, nsn, MPI_REAL8,&
            & MPI_SUM, parallel%commy, mpinfo)
#else
      rhoS=rhoS/dvol
#endif
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE smear_updaterho
!############################################################!
!*For :k-space Chebyshev Filtering                           !
!*Author : Qiang Xu                                          !
!*Date   : 2017-11-28                                        !
!############################################################!
   !-------------------------CheFSI---------------------------
   SUBROUTINE CheFSI(rhoS,psi,eval)
      USE parameters , ONLY : NMITER,RTOL,ETOL,Nstates,Nspin &
                            &,LFIRST,LINRHO,Lrandom,MALPHA &
                            &,Nstates_global
      USE succeed , ONLY: Llastrho,get_rho,get_psi
#ifndef MPI
      USE grid_module , ONLY : n1,n2,n3,n,nk,nsn,dvol,KPT
#else
      USE grid_module , ONLY : n1,n2,n3,n,nk,nsn,dvol,KPT&
           &,global_n1,global_n2,global_n3,global_n,rho_trans1D&
           &, rho_trans3D
      use smpi_math_module, only: mpinfo,mpi_real8,mpi_complex16&
           & ,mpi_sum,start_time,end_time,write_time
#endif
      USE potential_module, ONLY:cal_veff
      USE mixer_module
      USE energy_module
      USE struct_module , ONLY : ncharge,natom
      USE smearing_module , ONLY : wke,smear_init,destroy_smear,fermilevel
      USE chebyshev_module , ONLY : first_CheSCF_random,first_CheSCF_sto
      use m_time_evaluate, only: filename
      !
      IMPLICIT NONE
      REAL(DP),INTENT(INOUT)  :: rhoS(:,:,:,:)
      COMPLEX(DP),INTENT(OUT) :: psi(:,:,:,:)
      REAL(DP),INTENT(OUT)    :: eval(:,:,:)
      !LOCAL
      ! REAL(DP) :: rhoSd(n1,n2,n3,NSPIN)
      REAL(DP) :: veff(n1,n2,n3,nspin),veffd(n1,n2,n3,Nspin)
      REAL(DP) :: fme,ne,ets
      REAL(DP) :: drho,dtoten,toten,totend,res
      INTEGER(I4B) :: Iter , nev, Iwd
      integer(i4b) :: i
#ifdef MPI
      REAL(DP),allocatable :: rho_1D(:,:),rhod_1D(:,:)
      REAL(DP) :: rho_global(global_n1,global_n2,global_n3)
      COMPLEX(DP) :: psi_global(global_n1,global_n2,global_n3)
      integer(i4b) :: ix,iy,iz
#endif
!!!!xinyu
write(*,*)"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
write(*,*)struct%pos
write(*,*)lat_mat
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !number of state we need
      nev=Nstates
      ne=REAL(ncharge,8)
      !inirial smear
#ifndef MPI
      CALL smear_init(nev)
#else
      CALL smear_init(Nstates_global)
#endif
      !initial mixer
      CALL init_mixer_data()
      !Calculate the effective potential
      if(Llastrho)then
         CALL get_rho(natom,struct%poscar,parallel%mygrid_range(3),Nspin,rhoS,dvol)
         CALL get_psi(parallel%mygrid_range(3),nev,nk,Nspin,psi)
      CALL cal_veff(rhoS,veff)
      veffd=veff
      else
      CALL cal_veff(rhoS,veff)
      !first step,don't consider spin
      IF(LFIRST)THEN
         !First step by Arpack
         print*,'WARNING: parallel arpack some times err'
         CALL solver_spin(rhoS,psi,eval,nev,100.d0)
      ELSE
         IF(Lrandom)THEN
            !first step by random subspace
            CALL first_CheSCF_random(rhoS,nev,psi,eval)
            ! print *,"only support isolated system"
         ELSE
            !first step by Rayleigh-Ritz step(STO+PW)
            CALL first_CheSCF_sto(veff,nev,psi,eval)
            ! print *,"only support isolated system"
         ENDIF
      ENDIF
! #ifdef MPI
!       ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
!       iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
!       open(1111+parallel%myid,file="rho"//filename(10:11))
!       do i=1,nev,1
!       call MPI_ALLGATHERV(psi(:,1,1,1),parallel%mygrid_range(3)&
!            &,MPI_COMPLEX16,psi_global,parallel%recvcounts&
!            &,parallel%displs,MPI_COMPLEX16&
!            & ,parallel%commx,mpinfo)
!       ! do i=1,size(psi,2),1
!       !    write(1111+parallel%myid,'(10000(F8.2,1X))')real(psi(:,i,1,1))
!       ! enddo
!       write(1111+parallel%myid,*)psi_global
!       enddo
!       close(1111+parallel%myid)
!       stop
! #endif
#ifdef MPI
      allocate(rho_1D(n,NSPIN))
      ! do i=1,NSPIN
      !    call rho_trans1D(rhoS(:,:,:,i),rho_1D(:,i))
      ! enddo
      CALL smear_updaterho(nev,ne,psi,eval,wke,fme,ets,rho_1D)
      do i=1,NSPIN
         call rho_trans3D(rho_1D(:,i),rhoS(:,:,:,i))
      enddo
      deallocate(rho_1D)
#else
      CALL smear_updaterho(nev,ne,psi,eval,wke,fme,ets,rhoS)
#endif
      ! print *,'fme',fme
      CALL totalenergy(psi,rhoS,eval,fme,ets,.FALSE.)
      toten=Etot
      ! print *,'toten',toten
      ! veffd=veff
      veffd=veff
      CALL cal_veff(rhoS,veff)
      !first step use simple mixing
      ! rhoS=rhoSd+MALPHA*rhoS
      CALL mixing(0,veff,veffd,res)
      endif
! #ifdef MPI
!       ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
!       iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
!       call MPI_ALLGATHERV(rhoS(ix,iy,1,1),parallel%mygrid_range(3)&
!            &,MPI_REAL8,rho_global,parallel%recvcounts&
!            &,parallel%displs,MPI_REAL8&
!            & ,parallel%commx,mpinfo)
!       open(1111+parallel%myid,file="rho"//filename(10:11))
!       ! do i=1,size(psi,2),1
!       !    write(1111+parallel%myid,'(10000(F8.2,1X))')real(psi(:,i,1,1))
!       ! enddo
!       write(1111+parallel%myid,*)rho_global
!       close(1111+parallel%myid)
!       print*,'sum rho_global1',sum(rho_global)
!       stop
! #endif
! #ifdef MPI
!       print *,'shape(rho)',shape(rhoS)
!       open(1111+parallel%myid,file="rho"//filename(10:11))
!       ! do i=1,size(psi,2),1
!       !    write(1111+parallel%myid,'(10000(F8.2,1X))')real(psi(:,i,1,1))
!       ! enddo
!       write(1111+parallel%myid,*)rhoS
!       close(1111+parallel%myid)
!       stop
! #endif
      !toten=509612.d0
      !iter
      DO Iter=1,NMITER
#ifdef MPI
      CALL start_time('scf',.true.)
#endif
         ! rhoSd=rhoS
         totend=toten
         !call chebyshev filter
         CALL filter_spin(veff,psi,eval)

         !update rho
#ifdef MPI
         allocate(rho_1D(n,NSPIN))
         ! do i=1,NSPIN
         !    call rho_trans1D(rhoS(:,:,:,i),rho_1D(:,i))
         ! enddo
         CALL smear_updaterho(nev,ne,psi,eval,wke,fme,ets,rho_1D)
         ! print*,'update rho'
         do i=1,NSPIN
            call rho_trans3D(rho_1D(:,i),rhoS(:,:,:,i))
         enddo
         deallocate(rho_1D)
#else
         CALL smear_updaterho(nev,ne,psi,eval,wke,fme,ets,rhoS)
#endif

         !> calculate energy
         CALL totalenergy(psi,rhoS,eval,fme,ets,.TRUE.)
         toten=Etot

         !check
         dtoten=ABS(toten-totend)/natom*hart2ev
#ifndef MPI
            
!!!!
         WRITE(6,"(1X,A8,I4,1X,A6,ES18.9,1X,A6,ES15.7,1X,A8,ES15.7)")     &
              &      '>CheFSI:',iter,'Energy',toten*hart2ev,'Res(E)',dtoten ,  &
              &  'Res(rho)',res
#else
         IF(parallel%isroot) WRITE(6,"(1X,A8,I4,1X,A6,ES18.9,1X,A6,ES15.7,1X,A8,ES15.7)")     &
              &      '>CheFSI:',iter,'Energy',toten*hart2ev,'Res(E)',dtoten ,  &
              &  'Res(rho)',res
#endif
         IF(res<RTOL.AND.dtoten<ETOL)THEN
            IWD=IWD+1
            IF(IWD==3) EXIT
         ELSE
            IWD=0
         ENDIF
         ! drho=SUM(ABS(rhoS-rhoSd))/nsn
         ! dtoten=(toten-totend)*hart2ev/natom
         ! !===============================================================
         ! WRITE(6,"(1X,A8,I4,1X,A6,E15.7,1X,A4,E15.7)") &
         !      &      '>CheFSI:',iter,'dTOTEN',dtoten,      &
         !      &  'dRHO',drho
         ! IF(drho<RTOL.AND.ABS(dtoten)<ETOL) EXIT
         !===============================================================
         CALL Cal_Veff(rhoS,veff)
         !mixing potential
         CALL mixing(Iter,veff,veffd,res)
#ifdef MPI
      CALL end_time('scf',.true.)
      if(parallel%isroot)CALL write_time('scf',.true.)
#endif
      ENDDO
      CALL destroy_mixer()
      IF(Iter==NMITER) THEN
         WRITE(6,*) 'SCF:failed STOP'
         !STOP
      ENDIF

#ifdef MPI
      if(parallel%isroot)WRITE(6,*)'SCF DONE!STEP is:',Iter
#else
      WRITE(6,*)'SCF DONE!STEP is:',Iter
#endif
      !call total energy
      CALL totalenergy(psi,rhoS,eval,fme,ets,.TRUE.)
      !print*,'Etot',Etot
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !>> debug, out the total psi
! #ifdef DEBUG
!       ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
!       iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
!       call MPI_ALLGATHERV(psi(:,1,1,1),parallel%mygrid_range(3)&
!            &,MPI_COMPLEX16,psi_global,parallel%recvcounts&
!            &,parallel%displs,MPI_COMPLEX16&
!            & ,parallel%commx,mpinfo)
!       open(1111+parallel%myid,file="rho"//filename(10:11))
!       ! do i=1,size(psi,2),1
!       !    write(1111+parallel%myid,'(10000(F8.2,1X))')real(psi(:,i,1,1))
!       ! enddo
!       write(1111+parallel%myid,*)psi_global
!       close(1111+parallel%myid)
!       stop
! #endif
   ENDSUBROUTINE CHeFSI
   !-----------------------filter_spin------------------------
   SUBROUTINE filter_spin(veff,psi,eval)
      !need to be improve for spin
      USE parameters , ONLY : NSPIN
      USE Grid_module , ONLY : n1,n2,n3, nk
      ! USE potential_module , ONLY : cal_veff
      USE chebyshev_module , ONLY : cheby_filter_RR
      IMPLICIT NONE
      REAL(DP),INTENT(IN) :: veff(:,:,:,:)
      COMPLEX(DP),INTENT(INOUT) :: psi(:,:,:,:)
      REAL(DP),INTENT(INOUT) :: eval(:,:,:)
      !LOCAL
      ! REAL(DP) :: veff(n1,n2,n3,NSPIN)
      INTEGER(I4B) :: Is,Ik
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! CALL cal_veff(rhoS,veff)
      ! print *,"only support isolated system"
      DO Is=1,NSPIN
      DO Ik=1,nk !> kfilter
         CALL cheby_filter_RR(Ik,veff(:,:,:,Is),psi(:,:,Ik,Is)&
              &,eval(:,Ik,Is))
      ENDDO
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE filter_spin
   !-------------------------kfilter--------------------------
   ! SUBROUTINE kfilter(veff,X,RRv)
   !    USE parameters , ONLY : LRROrthNorm
   !    USE Grid_module , ONLY : nk
   !    USE chebyshev_module , ONLY : cheby_filter_RR
   !    IMPLICIT NONE
   !    REAL(DP),INTENT(IN) :: veff(:,:,:)
   !    COMPLEX(DP),INTENT(INOUT) :: X(:,:,:)
   !    REAL(DP),INTENT(INOUT) :: RRv(:,:)
   !    !
   !    INTEGER(I4B) :: Ik
   !    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !    DO Ik=1,nk
   !       CALL cheby_filter_RR(Ik,veff,X(:,:,Ik),RRv(:,Ik))
   !    ENDDO
   !    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ! ENDSUBROUTINE kfilter
!############################################################!
!*For :BvK cell self consistent                              !
!*Author : Qiang Xu                                          !
!*Date   : 2018-03-01                                        !
!############################################################!
   SUBROUTINE ArpackSCF_r(rhoS,psi,eval)
      USE parameters , ONLY : NMITER,RTOL,ETOL,Nstates,NSPIN,MALPHA
      USE grid_module , ONLY : n1,n2,n3,n,nsn,nk,dvol
      USE mixer_module
      USE energy_module
      USE Struct_module , ONLY : ncharge,natom
      USE smearing_module , ONLY : wke,smear_init,destroy_smear,Fermilevel

      IMPLICIT NONE
      !IN/OUT
      REAL(DP),INTENT(INOUT)  :: rhoS(:,:,:,:)
      REAL(DP),INTENT(OUT) :: psi(:,:,:,:)
      REAL(DP),INTENT(OUT)    :: eval(:,:,:)
      !LOCAL
      REAL(DP) :: rhoSd(n1,n2,n3,NSPIN)
      REAL(DP) :: fme,ne,ets
      REAL(DP) :: drho,dtoten,toten,totend,res
      INTEGER(I4B) :: Iter,nev
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      nev=Nstates
      ne=REAL(ncharge,8)
      !inirial smear
      CALL smear_init(nev)
      !initial mixer
      CALL init_mixer_data()
      !++++++++++++++++++++++++++++++First diag+++++++++++++++++++++++++++++++
      rhoSd(:,:,:,:)=rhoS(:,:,:,:)
      !First diag
      CALL solver_spin_r(rhoS,psi(:,:,1,:),eval(:,1,:),nev,100.D0)
      !update rhoS
      CALL real_smear_updaterho(nev,ne,psi,eval,wke,fme,ets,rhoS)
      CALL BvK_totalenergy(rhoS,eval(:,1,:),fme,ets,.FALSE.)
      toten=Etot
      !first mixing by simple mixing
      rhoS=rhoSd+MALPHA*rhoS
      !cycle for SCF
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO Iter=1,NMITER

         !store old data
         rhoSd(:,:,:,:)=rhoS(:,:,:,:)
         totend=toten

         !diag H
         CALL solver_spin_r(rhoS,psi(:,:,1,:),eval(:,1,:),nev,5D-5)

         !update rhoS
         CALL real_smear_updaterho(nev,ne,psi,eval,wke,fme,ets,rhoS)
         CALL BvK_totalenergy(rhoS,eval(:,1,:),fme,ets,.FALSE.)
         toten=Etot

         !check for exit
         drho=SUM(ABS(rhoS-rhoSd))/nsn
         dtoten=(toten-totend)*hart2ev/natom
         !===============================================================
         WRITE(6,"(1X,A9,I4,1X,A6,E15.7,1X,A4,E15.7)") &
              &      '>ArpkSCF:',iter,'dTOTEN',dtoten,      &
              &  'dRHO',drho
         !===============================================================
         IF(drho<RTOL.AND.ABS(dtoten)<ETOL) EXIT
         !mixing for update phi and rho
         CALL mixing(Iter,rhoS,rhoSd,res)
      ENDDO
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL destroy_mixer()
      IF(Iter==NMITER) THEN
         WRITE(6,*) 'SCF:failed STOP'
         STOP
      ENDIF

      WRITE(6,*)'SCF DONE!STEP is:',Iter
      !call total energy
      CALL BvK_totalenergy(rhoS,eval(:,1,:),fme,ets,.TRUE.)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE ArpackSCF_r
   !-----------------------Ispin solver-----------------------
   SUBROUTINE solver_spin_r(rhoS,psi,eval,nev,diagTOL)
      !
      USE grid_module , ONLY : n1,n2,n3,n
      USE parameters , ONLY : NSPIN
      USE potential_module , ONLY : cal_veff
      !arpk
      USE Arpack_module , ONLY : real_diagH_arpack
      USE math, ONLY : rsort_eigen
      IMPLICIT NONE
      REAL(DP),INTENT(IN) :: rhoS(:,:,:,:)
      REAL(DP),INTENT(IN) :: diagTOL !for diag
      REAL(DP),INTENT(OUT) :: psi(:,:,:)
      REAL(DP),INTENT(OUT) :: eval(:,:)
      INTEGER,INTENT(IN) :: nev
      !LOCAL
      INTEGER(I4B) :: Is
      REAL(DP) :: veff(n1,n2,n3,NSPIN)
      !Arpk
      INTEGER(I4B) :: nec,isc,info
      REAL(DP)     :: resid_int(n)
      INTEGER(I4B) :: maxmvs
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !
      DO Is=1,NSPIN
         CALL cal_veff(rhoS,veff)
         ! print *,"only support isolated system"
         !initialize arpk
         info=0
         resid_int(:)=0.d0
         diagH:DO isc=1,5
                   maxmvs=isc*30000
                   CALL real_diagH_arpack(veff(:,:,:,Is),nev,psi(:,:,Is),eval(:,Is), &
                      &    resid_int,nec,info,maxmvs,diagTOL)
                   IF(nec>=nev)THEN
                      EXIT diagH
                   ELSE
                      WRITE(6,*)'diag(H) failed,restarting'
                      !STOP          !I don't want to do this
                      !info=0
                   ENDIF

              ENDDO diagH

         IF(isc>=5)THEN
            WRITE(6,*) 'diag(H) failed,STOP,increase NADST and try again'
            STOP
         ENDIF

         IF(nev>1)THEN
             CALL rsort_eigen(nev,eval(:,Is),psi(:,:,Is))
         ENDIF

      ENDDO
!print*,'solver'
!print*,eval(:,:)
!pause
!STOP
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE solver_spin_r
   !--------------------smear update  rho---------------------
   SUBROUTINE real_smear_updaterho(nev,ne,psi,eval,wke_l,fme,ets,rhoS)
      USE smearing_module
      USE parameters , ONLY : NSPIN
      USE grid_module , ONLY : n,nsn,dvol,nk,KPT
      !
      IMPLICIT NONE
      !INPUT
      INTEGER(I4B),INTENT(IN) :: nev   !number of states
      REAL(DP),INTENT(IN)     :: ne !number of charge
      REAL(DP),INTENT(IN) :: psi(:,:,:,:) !eigen-states
      REAL(DP),INTENT(IN)     :: eval(:,:,:) !eigen-values
      REAL(DP),INTENT(OUT)     :: wke_l(:,:,:)  !weight of k
      REAL(DP),INTENT(OUT) :: ets !energy of electronic entropy
      !OUT PUT
      REAL(DP),INTENT(OUT) :: fme !fermi leval
      REAL(DP),INTENT(OUT) :: rhoS(nsn)
      !
      INTEGER(I4B) :: Ik,Ispin,Id,lft,rit !for index
      INTEGER(I4B) :: Iocc !for occupation state index
      REAL(DP) :: prho(n) ! 1D density in each k
      REAL(DP) :: small=1e-14
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !calculate the weight of k points
      CALL Fermilevel(ne,nev,nk,KPT%wk,eval,wke_l,fme,ets)
      !
      DO Ispin=1,NSPIN
         prho(:)=0.d0
         DO Ik=1,nk
            DO Iocc=1,nev
               !occ number
               IF(ABS(wke_l(Iocc,Ik,Ispin))<small) CYCLE
               !
               DO Id=1,n
                  prho(Id)=prho(Id)+psi(Id,Iocc,Ik,Ispin)**2*wke_l(Iocc,Ik,Ispin)
               ENDDO
            ENDDO
         ENDDO
         !
         lft=(Ispin-1)*n+1
         rit=Ispin*n
         !
         rhoS(lft:rit)=prho(:)
      ENDDO
      !dvol for psi normal
      rhoS(:)=rhoS(:)/dvol
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE real_smear_updaterho
!############################################################!
!*For    : BvK cell ChebyShev filtering + RR                 !
!*Author : Qiang Xu                                          !
!*Date   : 2018-03-06                                        !
!############################################################!
   SUBROUTINE BvK_CheFSI(rhoS,psi,eval)
      USE parameters , ONLY : NMITER,RTOL,ETOL,Nstates,Nspin &
                            &,LFIRST,LINRHO,Lrandom,MALPHA
      USE grid_module , ONLY : n1,n2,n3,n,nk,nsn,dvol
      USE mixer_module
      USE energy_module
      USE struct_module , ONLY : ncharge,natom
      USE smearing_module , ONLY : wke,smear_init,destroy_smear,fermilevel
#ifndef MPI
      USE chebyshev_module , ONLY : BvK_first_CheSCF_sto_rand
#endif
      !USE Potential_module , ONLY : cal_veff !Xltdelect
      !
      IMPLICIT NONE
      REAL(DP),INTENT(INOUT)  :: rhoS(:,:,:,:)
      REAL(DP),INTENT(OUT) :: psi(:,:,:,:)
      REAL(DP),INTENT(OUT)    :: eval(:,:,:)
      !LOCAL
      REAL(DP) :: rhoSd(n1,n2,n3,NSPIN)
      REAL(DP) :: fme,ne,ets
      REAL(DP) :: drho,dtoten,toten,totend,res
      INTEGER(I4B) :: Iter , nev
      !REAL(DP) :: veff(n1,n2,n3,NSPIN) ! Xltdelect
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !number of state we need
      nev=Nstates
      ne=REAL(ncharge,8)
      !initial smear
      CALL smear_init(nev)
      !================================================
      !!##XLT test I/O FOR DETECT ISO
      !open(unit=1118,file="RHO_ISO")
      !READ(1118,*)
      !READ(1118,*)rhoS
      !close(1118)
      !goto 1188
      !================================================
      !initial mixer
      CALL init_mixer_data()
      !first step,by STO+random
#ifndef MPI
      CALL BvK_first_CheSCF_sto_rand(rhoS,nev,psi(:,:,1,:),eval(:,1,:))
      ! print *,"only support isolated system"
#endif
      rhoSd(:,:,:,:)=rhoS(:,:,:,:)
      CALL real_smear_updaterho(nev,ne,psi,eval,wke,fme,ets,rhoS)
      CALL BvK_totalenergy(rhoS,eval(:,1,:),fme,ets,.FALSE.)
      toten=Etot
      !first step use simple mixing
      rhoS=rhoSd+MALPHA*rhoS
      !iter
      DO Iter=1,NMITER
         rhoSd(:,:,:,:)=rhoS(:,:,:,:)
         totend=toten
         !call chebyshev filter
         CALL BvK_filter_spin(rhoS,psi(:,:,1,:),eval(:,1,:))

         !update rho
         CALL real_smear_updaterho(nev,ne,psi,eval,wke,fme,ets,rhoS)
         CALL BvK_totalenergy(rhoS,eval(:,1,:),fme,ets,.FALSE.)
         toten=Etot

         !check
         drho=SUM(ABS(rhoS-rhoSd))/nsn
         dtoten=(toten-totend)*hart2ev/natom
         !===============================================================
         WRITE(6,"(1X,A8,I4,1X,A6,E15.7,1X,A4,E15.7)") &
              &      '>CheFSI:',iter,'dTOTEN',dtoten,      &
              &  'dRHO',drho
         IF(drho<RTOL.AND.ABS(dtoten)<ETOL) EXIT
         !===============================================================
         !mixing potential
         CALL mixing(Iter,rhoS,rhoSd,res)
      ENDDO
      CALL destroy_mixer()
      IF(Iter==NMITER) THEN
         WRITE(6,*) 'SCF:failed STOP'
         !STOP
      ENDIF

      WRITE(6,*)'SCF DONE!STEP is:',Iter
      !call total energy
      !1188 &
      CALL BvK_totalenergy(rhoS,eval(:,1,:),fme,ets,.TRUE.)
      print*,"sum rho",sum(rhoS)
      print*,"dvol",dvol
      !open(unit=111, file="RHO")
      !write(unit=111, fmt=*)rhoS
      !close(111)
      !open(unit=1118,file="EVAL")
      !write(1118,*)eval(:,1,:)
      !close(1118)
      !open(unit=1118,file="FME")
      !write(1118,*)fme
      !close(1118)
      !open(unit=1118,file="ETS")
      !write(1118,*)ets
      !close(1118)
      !open(unit=1118,file="WKE")
      !write(1118,*)wke
      !close(1118)
      !CALL cal_veff(rhoS,veff)
      !open(unit=1118,file="VEFF_PERIOD")
      !write(1118,*),size(veff,1),size(veff,2),size(veff,3)
      !write(1118,*)veff
      !close(1118)
      !print*,'Etot',Etot
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE BvK_CHeFSI
   !-----------------------filter_spin------------------------
   SUBROUTINE BvK_filter_spin(rhoS,X,D)
      !need to be improve for spin
      USE parameters , ONLY : Nspin
      USE Grid_module , ONLY : n1,n2,n3
      USE potential_module , ONLY : cal_veff
      USE chebyshev_module , ONLY : cheby_filtering_GRRr
      IMPLICIT NONE
      REAL(DP),INTENT(IN) :: rhoS(:,:,:,:)
      REAL(DP),INTENT(INOUT) :: X(:,:,:)
      REAL(DP),INTENT(INOUT) :: D(:,:)
      !LOCAL
      REAL(DP) :: veff(n1,n2,n3,NSPIN)
      INTEGER(I4B) :: Is
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL cal_veff(rhoS,veff)
      ! print *,"only support isolated system"
      DO Is=1,Nspin
         CALL cheby_filtering_GRRr(veff(:,:,:,Is),X(:,:,Is),D(:,Is))
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE BvK_filter_spin
!############################################################!
!*For    : BvK cell Chebyshev filtering + pRR                !
!*Author : Qiang Xu                                          !
!*Date   : 2018-03-07                                        !
!############################################################!
   SUBROUTINE pRR_filter_spin(Nfs,Nfe,rhoS,X,Efr,fme,ets,Pbar)
      !need to be improve for spin
      USE parameters , ONLY : Nspin,nev=>Nstates
      USE Grid_module , ONLY : n1,n2,n3,n,nk,KPT,dvol
      USE potential_module , ONLY : cal_veff
      USE chebyshev_module , ONLY : Cheby_filtering_pRRr
      USE smearing_module , ONLY : wke,fermilevel
      USE Lapack_module , ONLY : matmat_real
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN) :: Nfs  !number of f-occ states
      REAL(DP),INTENT(INOUT)  :: rhoS(:,:,:,:) &
                            &, Nfe
      REAL(DP),INTENT(OUT) ::  Pbar(:,:,:) &
                            &, X(:,:,:) &
                            &, fme &
                            &, ets
      REAL(DP),INTENT(INOUT) :: Efr(:,:,:)
      !LOCAL
      REAL(DP) :: veff(n1,n2,n3,NSPIN)
      REAL(DP) :: Cfr(nev,Nfs,Nspin),CfrN(nev,Nfs),PbarX(nev,n)
      INTEGER(I4B) :: Is,Ifr,I,Ix,Iy,Iz
      !
      REAL(DP) :: a,b,al,Cnsp
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL cal_veff(rhoS,veff)
      ! print *,"only support isolated system"
      Cnsp=2.d0/Nspin
      !
      DO Is=1,Nspin
         !filtering
         CALL Cheby_filtering_pRRr(veff(:,:,:,Is)   &
             & ,X(:,:,Is),Nfs,Cfr(:,:,Is),Efr(:,1,Is))
      ENDDO
      !calculate the project density matrix
      !1)find fermi-level and occupation number
      CALL Fermilevel(Nfe,Nfs,nk,KPT%wk,Efr,wke,fme,ets)
      !2)calculate Pbar
      DO Is=1,Nspin
         !
         DO Ifr=1,Nfs
            CfrN(:,Ifr)=Cfr(:,Ifr,Is)*( wke(Ifr,1,Is)-Cnsp )
         ENDDO
         !
         CALL matmat_real(CfrN(:,:),Cfr(:,:,Is),'N','T',Pbar(:,:,Is))
         !Pbar(:,:,Is)=MATMUL(CfrN(:,:),TRANSPOSE(Cfr(:,:,Is)))
         !print*,'Pbar',SIZE(Pbar,1),SIZE(Pbar,2)
         !print*,'CfrN',SIZE(CFRN,1),SIZE(CfrN,2)
         !print*,'Cfr',SIZE(Cfr,1),SIZE(Cfr,2)
         !
         !OPEN(1111,FILE='Pbar.dat')
         !   WRITE(1111,*) Pbar(:,:,1)
         !CLOSE(1111)
         !OPEN(1112,FILE='Cfrn.dat')
         !   WRITE(1112,*) CfrN(:,:)
         !CLOSE(1112)
         !OPEN(1113,FILE='Cfr.dat')
         !   WRITE(1113,*) Cfr(:,:,1)
         !CLOSE(1113)
         !STOP
         DO I=1,nev
            Pbar(I,I,Is)=Cnsp+Pbar(I,I,Is)
         ENDDO
      ENDDO
      !update the density
      DO Is=1,Nspin
         CALL matmat_real(Pbar(:,:,Is),X(:,:,Is),'N','T',PbarX(:,:))
         !PbarX(:,:)=MATMUL(Pbar(:,:,Is),TRANSPOSE(X(:,:,Is)))
         !OPEN(1111,FILE='Pbar.dat')
         !   WRITE(1111,*) Pbar
         !CLOSE(1111)
         !OPEN(1112,FILE='X.dat')
         !   WRITE(1112,*) X
         !CLOSE(1112)
         !OPEN(1113,FILE='PX.dat')
         !   WRITE(1113,*) PbarX
         !CLOSE(1113)
         !calculate spin density
         I=0
         DO Iz=1,n3
         DO Iy=1,n2
         DO Ix=1,n1
            I=I+1
            rhoS(Ix,Iy,Iz,Is)=DOT_PRODUCT(X(I,:,Is),PbarX(:,I))
         ENDDO
       ENDDO
         ENDDO
      ENDDO
      rhoS(:,:,:,:)=rhoS(:,:,:,:)/dvol
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE pRR_filter_spin
   !-----------------------filter_spin------------------------
   SUBROUTINE pRR_CheFSI(rhoS,Phi)
      USE parameters , ONLY : NMITER,RTOL,ETOL,nev=>Nstates,Nspin &
                            &,LFIRST,LINRHO,Lrandom,MALPHA,Nfs=>NpRR
      USE grid_module , ONLY : n1,n2,n3,n,nk,nsn,dvol
      USE mixer_module
      USE energy_module
      USE struct_module , ONLY : ncharge,natom
      USE smearing_module , ONLY : wke,smear_init,destroy_smear,fermilevel
#ifndef MPI
      USE chebyshev_module , ONLY : BvK_first_CheSCF_sto_rand
#endif
      !
      IMPLICIT NONE
      REAL(DP),INTENT(INOUT)  :: rhoS(:,:,:,:)
      REAL(DP),INTENT(OUT) :: Phi(:,:,:,:)
      !LOCAL
      REAL(DP) :: rhoSd(n1,n2,n3,NSPIN)
      REAL(DP) :: fme,ets,Ne,Nfe,eval(nev,nk,Nspin)
      REAL(DP) :: drho,dtoten,toten,totend,res
      REAL(DP),ALLOCATABLE,DIMENSION(:,:,:) :: Efr , Pbar
      INTEGER(I4B) :: Iter
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !number of state we need
      !-------------------------------------------------------
      Ne=REAL(ncharge,8)             !Number of total electron
      Nfe=ne-((nev-Nfs)*2.d0)        !Number of frac-electron
      IF(Nfe<=0.d0.OR.Nfe>Ne)THEN
         WRITE(6,*) 'pRR: Electron number wroung,Check NPRR'
         WRITE(6,*) 'fraction electron',Nfe,'total electron',ne
         STOP
      ENDIF
      WRITE(6,*) '[pRR-Number offracation electron]',Nfe
      !-------------------------------------------------------
      !allocate some array
      !-------------------------------------------------------
      ALLOCATE(Efr(Nfs,nk,Nspin),Pbar(nev,nev,Nspin))
      !-------------------------------------------------------
      !inirial smear
      CALL smear_init(nev)
      !initial mixer
      CALL init_mixer_data()
      !first step,by STO+random
#ifndef MPI
      CALL BvK_first_CheSCF_sto_rand(rhoS,nev,Phi(:,:,1,:),eval(:,1,:))
      ! stop "only support isolated system"
#endif
      rhoSd(:,:,:,:)=rhoS(:,:,:,:)
      CALL real_smear_updaterho(nev,ne,Phi,eval,wke,fme,ets,rhoS)
      CALL BvK_totalenergy(rhoS,eval(:,1,:),fme,ets,.FALSE.)
      !store up eigen-value
      Efr(1:Nfs,:,:)=eval(nev-Nfs+1:nev,:,:)
      !store old data
      toten=Etot
      !first step use simple mixing
      rhoS=rhoSd+MALPHA*rhoS
      !iter
      DO Iter=1,NMITER
         rhoSd(:,:,:,:)=rhoS(:,:,:,:)
         totend=toten
         !filter
         CALL pRR_filter_spin(Nfs,Nfe,rhoS,Phi(:,:,1,:),Efr,fme,ets,Pbar)
         CALL pRR_totalenergy(rhoS,Phi(:,:,1,:),Pbar,fme,ets,.FALSE.)
         toten=Etot
         !check
         drho=SUM(ABS(rhoS-rhoSd))/nsn
         dtoten=(toten-totend)*hart2ev/natom
         !===============================================================
         WRITE(6,"(1X,A8,I4,1X,A6,E15.7,1X,A4,E15.7)") &
              &      '>CheFSI:',iter,'dTOTEN',dtoten,      &
              &  'dRHO',drho
         IF(drho<RTOL.AND.ABS(dtoten)<ETOL) EXIT
         !===============================================================
         !mixing potential
         CALL mixing(Iter,rhoS,rhoSd,res)
      ENDDO
      CALL destroy_mixer()
      IF(Iter==NMITER) THEN
         WRITE(6,*) 'SCF:failed STOP'
         !STOP
      ENDIF

      WRITE(6,*)'SCF DONE!STEP is:',Iter
      !call total energy
      CALL pRR_totalenergy(rhoS,Phi(:,:,1,:),Pbar,fme,ets,.TRUE.)
      !print*,'Etot',Etot
      DEALLOCATE(Efr,Pbar)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE pRR_CHeFSI
!#endif
!############################################################!
!*For :Isolate cell self consistent                          !
!*Author : qiangxu, modified by xlt                          !
!*Date   : 2018-04-10                                        !
!############################################################!
#ifndef MPI
   SUBROUTINE ISO_ArpackSCF_r(rhoS,psi,eval)
      USE parameters , ONLY : NMITER,RTOL,ETOL,Nstates,NSPIN,MALPHA
      USE grid_module , ONLY : n1,n2,n3,n,nsn,nk,dvol
      USE mixer_module
      USE energy_module
      USE Struct_module , ONLY : ncharge,natom
      USE smearing_module , ONLY : wke,smear_init,destroy_smear,Fermilevel_iso
      USE potential_module , ONLY : V_accelerate

      IMPLICIT NONE
      !IN/OUT
      REAL(DP),INTENT(INOUT)  :: rhoS(:,:,:,:)
      REAL(DP),INTENT(OUT) :: psi(:,:,:,:)
      REAL(DP),INTENT(OUT)    :: eval(:,:,:)
      !LOCAL
      REAL(DP) :: rhoSd(n1,n2,n3,NSPIN)
      REAL(DP) :: fme,ne,ets
      REAL(DP) :: drho,dtoten,toten,totend,res
      INTEGER(I4B) :: Iter,nev,t1,t2
      CHARACTER(20) :: name
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ALLOCATE(V_accelerate(n1,n2,n3))
      nev=Nstates
      ne=REAL(ncharge,8)
      !inirial smear
      CALL smear_init(nev)
      !initial mixer
      CALL init_mixer_data()
      !++++++++++++++++++++++++++++++First diag+++++++++++++++++++++++++++++++
      !==================================================
      !#PLOT
      !open(1156,file="CHGCAR_0")
      !WRITE(1156,*)size(rhoS,1),size(rhoS,2),size(rhoS,3)
      !WRITE(1156,*)rhoS
      !CLOSE(1156)
      !==================================================
      rhoSd(:,:,:,:)=rhoS(:,:,:,:)
      !First diag
      CALL ISO_solver_spin_r(rhoS,psi(:,:,1,:),eval(:,1,:),nev,100.D0)
      !update rhoS
      CALL ISO_smear_updaterho(nev,ne,psi,eval,wke,fme,ets,rhoS)
      CALL ISO_totalenergy(rhoS,eval(:,1,:),fme,ets,.FALSE.)
      toten=Etot
      !first mixing by simple mixing
      rhoS=rhoSd+MALPHA*rhoS
      !==================================================
      !#PLOT
      !open(1156,file="CHGCAR_1")
      !WRITE(1156,*)size(rhoS,1),size(rhoS,2),size(rhoS,3)
      !WRITE(1156,*)rhoS
      !CLOSE(1156)
      !==================================================
      !cycle for SCF
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO Iter=1,NMITER

         !store old data
         rhoSd(:,:,:,:)=rhoS(:,:,:,:)
         totend=toten

         !diag H
         CALL system_clock(t1)
         CALL ISO_solver_spin_r(rhoS,psi(:,:,1,:),eval(:,1,:),nev,5D-5)
         !update rhoS
         CALL ISO_smear_updaterho(nev,ne,psi,eval,wke,fme,ets,rhoS)
         CALL mixing(Iter,rhoS,rhoSd,res)
         CALL ISO_totalenergy(rhoS,eval(:,1,:),fme,ets,.FALSE.)
         toten=Etot
         CALL system_clock(t2)
         print*,'#######sum time',(t2-t1)/10000.d0
         !check for exit
         drho=SUM(ABS(rhoS-rhoSd))/nsn*5
         dtoten=(toten-totend)*hart2ev/natom
         !==========================================
         !#PLOT
         !write(name,*)Iter+1
         !open(1156,file="CHGCAR_"//TRIM(ADJUSTL(name)))
         !WRITE(1156,*)size(rhoS,1),size(rhoS,2),size(rhoS,3)
         !WRITE(1156,*)rhoS
         !CLOSE(1156)
         !==========================================
         !===============================================================
         WRITE(6,"(1X,A9,I4,1X,A6,E15.7,1X,A4,E15.7)") &
             &      '>ArpkSCF:',iter,'dTOTEN',dtoten,      &
             &  'dRHO',drho
         !===============================================================
         IF(drho<RTOL.AND.ABS(dtoten)<ETOL) EXIT
         !mixing for update phi and rho
         !!!!!!#######drho=SUM(ABS(rhoS-rhoSd))/nsn
         !!!!!!#######IF(drho<RTOL)THEN
         !!!!!!#######   CALL ISO_totalenergy(rhoS,eval(:,1,:),fme,ets,.FALSE.)
         !!!!!!#######   toten=Etot
         !!!!!!#######ELSE
         !!!!!!#######   toten=0.d0
         !!!!!!#######ENDIF
         !!!!!!#######dtoten=(toten-totend)*hart2ev/natom
         !!!!!!#######!===============================================================
         !!!!!!#######WRITE(6,"(1X,A9,I4,1X,A6,E15.7,1X,A4,E15.7)") &
         !!!!!!#######     &      '>ArpkSCF:',iter,'dTOTEN',dtoten,      &
         !!!!!!#######     &  'dRHO',drho
         !!!!!!#######!===============================================================
         !!!!!!#######IF(drho<RTOL.AND.ABS(dtoten)<ETOL) EXIT
         !!!!!!#######!mixing for update phi and rho
         !!!!!!#######CALL mixing(Iter,rhoS,rhoSd)


      ENDDO
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL destroy_mixer()
      IF(Iter==NMITER) THEN
         WRITE(6,*) 'SCF:failed STOP'
         STOP
      ENDIF

      WRITE(6,*)'SCF DONE!STEP is:',Iter
      !call total energy
      CALL ISO_totalenergy(rhoS,eval(:,1,:),fme,ets,.TRUE.)
      !WRITE(6,"(F15.9/)"),((eval(t1,1,t2)*hart2ev,t1=1,size(eval,1)),t2=1,size(eval,3))
      DEALLOCATE(V_accelerate)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE ISO_ArpackSCF_r
#endif
   !-----------------------Ispin solver-----------------------
   SUBROUTINE ISO_solver_spin_r(rhoS,psi,eval,nev,diagTOL)
      !
      USE grid_module , ONLY : n1,n2,n3,n,rho_calc,ISO_Vsphere
      USE parameters , ONLY : NSPIN
      USE potential_module , ONLY : cal_veff
      !arpk
      USE Arpack_module , ONLY : ISO_diagH_arpack
      USE math, ONLY : rsort_eigen
      IMPLICIT NONE
      REAL(DP),INTENT(IN) :: rhoS(:,:,:,:)
      REAL(DP),INTENT(IN) :: diagTOL !for diag
      REAL(DP),INTENT(OUT) :: psi(:,:,:)
      REAL(DP),INTENT(OUT) :: eval(:,:)
      INTEGER,INTENT(IN) :: nev
      !LOCAL
      INTEGER(I4B) :: Is
      REAL(DP) :: veff(n1,n2,n3,NSPIN)
      !Arpk
      INTEGER(I4B) :: nec,isc,info
      REAL(DP)     :: resid_int(rho_calc%OneDLength)
      INTEGER(I4B) :: maxmvs,t1,t2
      REAL(DP)     :: Veff_sphere(rho_calc%OneDLength,NSPIN)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !
      DO Is=1,NSPIN
#ifndef MPI
         CALL cal_veff(rhoS,veff)
         ! print *,"only support isolated system with CheFSI"
#endif
         CALL ISO_Vsphere(veff,Veff_sphere)
         !initialize arpk
         info=0
         resid_int(:)=0.d0
         diagH:DO isc=1,5
                   maxmvs=isc*30000
                   CALL ISO_diagH_arpack(Veff_sphere(:,Is),nev,psi(:,:,Is),eval(:,Is), &
                      &    resid_int,nec,info,maxmvs,diagTOL)
                   IF(nec>=nev)THEN
                      EXIT diagH
                   ELSE
                      WRITE(6,*)'diag(H) failed,restarting'
                      !STOP          !I don't want to do this
                      !info=0
                   ENDIF

              ENDDO diagH

         IF(isc>=5)THEN
            WRITE(6,*) 'diag(H) failed,STOP,increase NADST and try again'
            STOP
         ENDIF

         IF(nev>1)THEN
             CALL rsort_eigen(nev,eval(:,Is),psi(:,:,Is))
         ENDIF

      ENDDO
      !print*,'solver'
      !print*,eval(:,:)
      !pause
      !STOP
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE ISO_solver_spin_r
   !--------------------smear update  rho---------------------
   SUBROUTINE ISO_smear_updaterho(nev,ne,psi,eval,wke_l,fme,ets,rhoS_out)
      USE smearing_module
      USE parameters , ONLY : NSPIN
      USE grid_module , ONLY : dvol,nk,KPT,rho_calc,n,nsn,n1,n2,n3,ISO_Rho2grid
#ifdef MPI
      USE smpi_math_module
#endif
      USE m_time_evaluate, ONLY: memory_sum,memory_free
      !
      IMPLICIT NONE
      !INPUT
      INTEGER(I4B),INTENT(IN) :: nev   !number of states
      REAL(DP),INTENT(IN)     :: ne !number of charge
#ifndef MPI
      REAL(DP),INTENT(IN) :: psi(:,:,:,:) !eigen-states
      REAL(DP),INTENT(IN)     :: eval(:,:,:) !eigen-values
#else
      REAL(DP),INTENT(IN) :: psi(:,:,:) !eigen-states
      REAL(DP),INTENT(IN)     :: eval(:,:) !eigen-values
#endif
      REAL(DP),INTENT(OUT)     :: wke_l(:,:,:)  !weight of k
      REAL(DP),INTENT(OUT) :: ets !energy of electronic entropy
      !OUT PUT
      REAL(DP),INTENT(OUT) :: fme !fermi leval
#ifndef MPI
      REAL(DP),INTENT(OUT) :: rhoS_out(n1,n2,n3,NSPIN)
#else
      REAL(DP),INTENT(OUT) :: rhoS_out(parallel%mygrid_range(3),NSPIN)
#endif
      !
#ifndef MPI
      REAL(DP)     :: rhoS(rho_calc%OneDLength,NSPIN)
#else
      REAL(DP) :: rhoS(parallel%mygrid_range(3),NSPIN) ! 1D density in each k
#endif
      INTEGER(I4B) :: Ik,Ispin,Id,lft,rit !for index
      INTEGER(I4B) :: Iocc !for occupation state index
#ifndef MPI
      REAL(DP) :: prho(rho_calc%OneDLength) ! 1D density in each k
#else
      REAL(DP) :: prho(parallel%mygrid_range(3)) ! 1D density in each k
#endif
      REAL(DP) :: small=1e-14
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      call memory_sum('updaterho_local',real(size(rhoS),DP)*DP+size(prho)*DP)
      !calculate the weight of k points
      CALL Fermilevel_iso(ne,nev,nk,KPT%wk,eval,wke_l,fme,ets)
      !
      DO Ispin=1,NSPIN
         prho(:)=0.d0
#ifndef MPI
         DO Ik=1,nk
#endif
            DO Iocc=1,nev
               !occ number
               IF(ABS(wke_l(Iocc,1,Ispin))<small) CYCLE
               !
#ifndef MPI
               DO Id=1,rho_calc%OneDLength
                  prho(Id)=prho(Id)+psi(Id,Iocc,Ik,Ispin)**2*wke_l(Iocc,Ik,Ispin)
#else
               DO Id=1,parallel%mygrid_range(3)
                  prho(Id)=prho(Id)+psi(Id,Iocc,Ispin)**2*wke_l(Iocc,1,Ispin)
#endif
               ENDDO
#ifndef MPI
            ENDDO
#endif
         ENDDO
         !
         !ft=(Ispin-1)*rho_calc%OneDLength+1
         !rit=Ispin*rho_calc%OneDLength
         !
         rhoS(:,Ispin)=prho(:)
      ENDDO
      !dvol for psi normal
      rhoS(:,:)=rhoS(:,:)/dvol
#ifndef MPI
      CALL ISO_Rho2grid(rhoS,rhoS_out)
#else
      CALL dcopy(parallel%mygrid_range(3),rhoS,1,rhoS_out,1)
#endif
      call memory_free('updaterho_local',real(size(rhoS),DP)*DP+size(prho)*DP)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE ISO_smear_updaterho
!############################################################!
!*For :Isolate cell self consistent                          !
!*Author : qiangxu, modified by xlt                          !
!*Date   : 2018-04-10                                        !
!############################################################!
   SUBROUTINE ISO_CheFSI(rhoS,psi,eval)
#include "my_macro.h"
#ifdef MPI
      USE smpi_math_module
      USE scalapack_module
      USE succeed , ONLY: Llastrho,get_rho,get_psi
#else
      USE succeed , ONLY: Llastrho,get_rho,get_psi !get_new_rho_psi
#endif
      USE parameters , ONLY : NMITER,RTOL,ETOL,Nstates,Nspin &
           &,LFIRST,LINRHO,Lrandom,MALPHA &
           &, step_fixrho
#ifdef MPI
      USE grid_module , ONLY : n1,n2,n3,n,nk,nsn,dvol,parallel_s2g,rho_calc
#else
      USE grid_module , ONLY : n1,n2,n3,n,nk,nsn,dvol,eigen_r,rho_calc,grid
#endif
      USE mixer_module
      USE energy_module
      USE struct_module , ONLY : ncharge,natom,struct
      USE smearing_module , ONLY : wke,smear_init,destroy_smear,fermilevel_iso
      USE chebyshev_module , ONLY : ISO_first_CheSCF_sto_rand
      USE potential_module , ONLY : V_accelerate,ACCELERATE
      !================================================
      !xlt 2018-5-18 USED FOR TEST
      ! USE grid_module , ONLY : ISO_Rho2grid,rho_calc
      !================================================
      USE array_io ,ONLY: output
      USE m_time_evaluate, ONLY: memory_sum, memory_free
      IMPLICIT NONE
#ifndef MPI
      REAL(DP),INTENT(INOUT)  :: rhoS(:,:,:,:)
      REAL(DP),INTENT(OUT) :: psi(:,:,:,:)
      REAL(DP),INTENT(OUT)    :: eval(:,:,:)
#else
      REAL(DP),INTENT(INOUT)  :: rhoS(:,:)
      REAL(DP),INTENT(OUT) :: psi(:,:,:)
      REAL(DP),INTENT(OUT)    :: eval(:,:)
#endif
      !LOCAL
#ifndef MPI
      REAL(DP) :: rhoSd(n1,n2,n3,NSPIN)
      REAL(DP) :: rhoSr(n1,n2,n3,NSPIN)
#else
      REAL(DP) :: rhoSd(parallel%mygrid_range(3),NSPIN)
      REAL(DP) :: rhoSr(parallel%mygrid_range(3),NSPIN)
#endif
      REAL(DP) :: fme,ne,ets,res
      REAL(DP) :: drho,dtoten,toten,totend
      INTEGER(I4B) :: Iter , nev, IWD
      CHARACTER(20) :: name
      REAL(DP) :: rho_out(n1,n2,n3)
      INTEGER(I4B) :: i
#ifdef MPI
      REAL(DP) :: t
#endif
      ! LOGICAL :: l1=.false.
      LOGICAL :: l1=.false.
#ifndef MPI
      INTEGER(I4B) :: In
#endif
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      call memory_sum('chefsi_local',real(size(rhoSd),DP)*DP*2+&
           & size(rho_out)*DP)
      !> ACCELERATE BY Vhartree
#ifndef MPI
      ALLOCATE(V_accelerate(n1,n2,n3))
#else
      ALLOCATE(V_accelerate(parallel%mygrid_range(3)))
#endif
      call memory_sum('V_accelerate',real(size(V_accelerate),DP)*DP)
      !> number of state we need
#ifdef MPI
      CALL start_time('init smear mix',l1)
#endif


#ifdef DEBUG
#ifndef MPI
#else
      if(parallel%isroot)then
      open(111,file="rhoS_Scf")
      write(111,*)rhoS
      close(111)
      endif
#endif
#endif


      nev=Nstates
      ne=REAL(ncharge,8)
      !> initial smear
      CALL smear_init(nev)
#ifdef NOSCF
      !> FOR TEST
      goto 1118
#endif
      !> initial mixer
      CALL init_mixer_data()
      !> first step,by STO+random
   IF(.not.Llastrho)THEN
#ifndef MPI
      CALL ISO_first_CheSCF_sto_rand(rhoS,nev,psi(:,:,1,:),eval(:,1,:))
#else
      ! CALL get_rho(natom,struct%poscar,parallel%mygrid_range(3),Nspin,rhoS)
      ! print*,'parallel%mygrid_range',parallel%mygrid_range
      CALL end_time('init smear mix',l1)
      CALL write_time('init smear mix',l1)
      ! CALL end_time('Init_data',l1)
      ! CALL write_time('Init_data',l1)
      CALL start_time('scf',l1)
      CALL start_time('init_space',l1)
      CALL ISO_first_CheSCF_sto_rand(rhoS,nev,psi,eval)
      call end_time('init_space',l1)
      ! call print_time('init_space',t)
      ! IF(parallel%isroot) WRITE(6,'(A20,1X,SP,F15.7,A3)')'*Init space&
           ! & time ->',t,'s'
#endif
      rhoSd=rhoS
      CALL ISO_smear_updaterho(nev,ne,psi,eval,wke,fme,ets,rhoS)
      ! print *,"fme",fme
      ! print *,"ets",ets
      ! print*,'eval',eval
!#ifndef MPI
!      open(111,file="rhoS_Scf")
!      write(111,*)rhoSd
!      close(111)
!#else
!      rho_out=0.d0
!      do i=1,rho_calc%OneDLength
!         rho_out(rho_calc%x(i),rho_calc%y(i),rho_calc%z(i))=rhoSd(i,1)
!      enddo
!      if(parallel%isroot)then
!      open(111,file="rhoS_Scf")
!      write(111,*)rho_out
!      close(111)
!      endif
!#endif
!      CALL get_new_rho_psi(natom,struct%poscar,n,n1,n2,n3,Nspin,rhoS,Nstates,eigen_r%wvf,grid%gVec)
!      psi=0.d0
!      DO i=1,rho_calc%OneDLength
!         In=rho_calc%x(i)+n1*(rho_calc%y(i)-1)+n1*n2*(rho_calc%z(i)-1)
!         psi(i,:,1,:)=eigen_r%wvf(In,:,1,:)
!      ENDDO
#ifndef MPI
      CALL ISO_totalenergy(rhoS,eval(:,1,:),fme,ets,.FALSE.)
#else
      CALL ISO_totalenergy(rhoS,eval,fme,ets,.FALSE.)
#endif
      toten=Etot
       ! print *,"toten",toten
       ! print*,'sumrhos',sum(rhoS)
       ! stop
      !first step use simple mixing
      !rhoS=rhoSd+MALPHA*rhoS
      !==========================================
      CALL mixing_iso(0,rhoS,rhoSd,res)
#ifdef MPI
      CALL end_time('scf',l1)
      CALL write_time('scf',l1)
#endif
   ELSE
#ifdef MPI
      CALL get_psi(parallel%mygrid_range(3),nev,Nspin,psi)
      CALL get_rho(natom,struct%poscar,parallel%mygrid_range(3),Nspin,rhoS,dvol)
#else
      CALL get_psi(rho_calc%OneDLength,nev,Nspin,psi)
      CALL get_rho(natom,struct%poscar,n,Nspin,rhoS,dvol)
!      CALL get_new_rho_psi(natom,struct%poscar,n,n1,n2,n3,Nspin,rhoS,Nstates,eigen_r%wvf,grid%gVec)
!      psi=0.d0
!      DO i=1,rho_calc%OneDLength
!         In=rho_calc%x(i)+n1*(rho_calc%y(i)-1)+n1*n2*(rho_calc%z(i)-1)
!         psi(i,:,1,:)=eigen_r%wvf(In,:,1,:)
!      ENDDO
!      ! open(1234,file='rho_succeed')
!      ! write(1234,*)rhoS
!      ! close(1234)
!      ! print*,'scf',rhoS(59328:59338,1,1,1)
!      ! print*,'scf',psi(59328:59338,1,1,1)
!      ! open(1234,file='rho_scf')
!      ! write(1234,*)grid%rhoS
!      ! close(1234)
!      ! open(1234,file='psi_scf')
!      ! write(1234,*)eigen_r%wvf(:,1,1,1)
!      ! close(1234)
#endif
      toten=Etot
      rhoSd=rhoS
      CALL mixing_iso(0,rhoS,rhoSd,res)
      ACCELERATE=.FALSE.
      ! print *,"scf_module",rhoS(1:10,1)
   ENDIF
      !iter
      DO Iter=1,NMITER
         rhoSd=rhoS
         totend=toten
#ifdef MPI
         MPI_SUMTIME("iso_filter",l1,ISO_filter_spin(rhoS,psi(:,:,:),eval(:,:)))
         if(step_fixrho<=Iter.or.(.not.Llastrho))then
            ! if(parallel%isroot)print *,'update rho'
            MPI_SUMTIME("iso_updaterho",l1,ISO_smear_updaterho(nev,ne,psi,eval,wke,fme,ets,rhoS))
         endif
         MPI_SUMTIME("iso_mixing",l1,mixing_iso(Iter,rhoS,rhoSd,res))
         MPI_SUMTIME("iso_totalenergy",l1,ISO_totalenergy(rhoS,eval(:,:),fme,ets,.FALSE.))
#else
         call ISO_filter_spin(rhoS,psi(:,:,1,:),eval(:,1,:))
         call ISO_smear_updaterho(nev,ne,psi,eval,wke,fme,ets,rhoS)
         call mixing_iso(Iter,rhoS,rhoSd,res)
         call ISO_totalenergy(rhoS,eval(:,1,:),fme,ets,.FALSE.)
#endif
         toten=Etot
         !> check
         dtoten=ABS(toten-totend) !/ABS(toten+totend) !*hart2ev/natom*50.d0
#ifndef MPI
         WRITE(6,"(1X,A8,I4,1X,A6,ES18.9,1X,A6,ES15.7,1X,A8,ES15.7)")     &
              &      '>CheFSI:',iter,'Energy',toten*hart2ev,'Res(E)',dtoten ,  &
              &  'Res(rho)',res
#else
         IF(parallel%isroot) WRITE(6,"(1X,A8,I4,1X,A6,ES18.9,1X,A6,ES15.7,1X,A8,ES15.7)")     &
              &      '>CheFSI:',iter,'Energy',toten*hart2ev,'Res(E)',dtoten ,  &
              &  'Res(rho)',res
#endif
         IF(res<RTOL.AND.dtoten<ETOL)THEN
            IWD=IWD+1
            IF(IWD==3) EXIT
         ELSE
            IWD=0
         ENDIF
         !===============================================================
         !mixing potential
         !=================================================
      ENDDO
      CALL destroy_mixer()
#ifndef MPI
      IF(Iter==NMITER) THEN
         WRITE(6,*) 'SCF:failed STOP'
         !STOP
      ENDIF

      WRITE(6,*)'SCF DONE!STEP is:',Iter
#else
      IF(Iter==NMITER+1) THEN
          IF(parallel%isroot)WRITE(6,*) 'SCF:failed STOP'
         ! STOP
      ENDIF

      IF(parallel%isroot)WRITE(6,*)'SCF DONE!STEP is:',Iter
      CALL write_sum_time("iso_filter",l1)
      CALL write_sum_time("iso_updaterho",l1)
      CALL write_sum_time("iso_mixing",l1)
      CALL write_sum_time("iso_totalenergy",l1)
      CALL write_sum_time("chebyshev_filter_scaled_ISO",l1)
      CALL write_sum_time("GRayleigh_Ritz_ISO",l1)
      CALL write_sum_time("parallel_communcation",.true.)
#endif
      !call total energy
      !=====================================================
#ifdef NOSCF
      1118 &!,FOR TEST POTENTIAL CALCULATE
      if(Llastrho)then
      CALL get_new_rho_psi(natom,struct%poscar,n,n1,n2,n3,Nspin,rhoS,Nstates,eigen_r%wvf,grid%gVec)
      open(1234,file='rho_scf')
      write(1234,*)grid%rhoS
      close(1234)
      stop
      endif
#endif
      !==================================================
#ifndef MPI
      CALL ISO_totalenergy(rhoS,eval(:,1,:),fme,ets,.TRUE.)
#else
      CALL ISO_totalenergy(rhoS,eval(:,:),fme,ets,.TRUE.)
      CALL parallel_s2g(psi(:,nev,1),rho_out)
      !open(1173,file='psi_nev')
      !write(1173,*)rho_out
      !close(1173)
#endif
      call memory_free('V_accelerate',real(size(V_accelerate),DP)*DP)
      DEALLOCATE(V_accelerate)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      call memory_free('chefsi_local',real(size(rhoSd),DP)*DP*2+&
           & size(rho_out)*DP)
   ENDSUBROUTINE ISO_CHeFSI
   !-----------------------filter_spin------------------------
   SUBROUTINE ISO_filter_spin(rhoS,X,D)
      !need to be improve for spin
      USE parameters , ONLY : Nspin
      USE Grid_module , ONLY : n1,n2,n3,rho_calc,ISO_Vsphere
      USE potential_module , ONLY : cal_veff_iso
      USE chebyshev_module , ONLY : cheby_filtering_GRRiso
#ifdef MPI
      USE smpi_math_module , ONLY : parallel
#endif
      USE m_time_evaluate, ONLY: memory_sum, memory_free
      IMPLICIT NONE
#ifndef MPI
      REAL(DP),INTENT(IN) :: rhoS(:,:,:,:)
#else
      REAL(DP),INTENT(IN) :: rhoS(:,:)
#endif
      REAL(DP),INTENT(INOUT) :: X(:,:,:)
      REAL(DP),INTENT(INOUT) :: D(:,:)
      !LOCAL
#ifndef MPI
      REAL(DP) :: veff(n1,n2,n3,NSPIN),veff_ISO(rho_calc%OneDLength,NSPIN)
#else
      REAL(DP) :: veff(parallel%mygrid_range(3),NSPIN)
#endif
      INTEGER(I4B) :: Is,te, tf, tg
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ifndef MPI
      call memory_sum('iso_filter_spin',real(size(veff),DP)*DP+size(veff)*DP)
#else
      call memory_sum('iso_filter_spin',real(size(veff),DP)*DP)
#endif
      CALL system_clock(te)
      CALL cal_veff_iso(rhoS,veff)
      CALL system_clock(tf)
      !CALL ISO_Vsphere(veff,veff_ISO)
      DO Is=1,Nspin
#ifndef MPI
         CALL cheby_filtering_GRRiso(veff(:,:,:,Is),X(:,:,Is),D(:,Is))
#else
         CALL cheby_filtering_GRRiso(veff(:,Is),X(:,:,Is),D(:,Is))
#endif
      END DO
      CALL system_clock(tg)
      !=============================================
      !print *,'calculate effectial time -->',(tf-te)/10000.d0
      !print *,'filter GeneralRealyRitz  time -->',(tg-tf)/10000.d0
      !print *,'total chebyshev filter time -->',(tg-te)/10000.d0
      !=============================================
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#ifndef MPI
      call memory_free('iso_filter_spin',real(size(veff),DP)*DP+size(veff)*DP)
#else
      call memory_free('iso_filter_spin',real(size(veff),DP)*DP)
#endif
   ENDSUBROUTINE ISO_filter_spin
   !----------------------CheFSI_Gamma-------------------------
   SUBROUTINE CheFSI_gamma(rhoS,psi,eval)
      USE parameters , ONLY : NMITER,RTOL,ETOL,Nstates,Nspin &
                            &,LFIRST,LINRHO,Lrandom,MALPHA &
                            &,Nstates_global
      USE succeed , ONLY: Llastrho,get_rho,get_psi
#ifndef MPI
      USE grid_module , ONLY : n1,n2,n3,n,nk,nsn,dvol,KPT
#else
      USE grid_module , ONLY : n1,n2,n3,n,nk,nsn,dvol,KPT&
           &,global_n1,global_n2,global_n3,global_n,rho_trans1D&
           &, rho_trans3D
      use smpi_math_module, only: mpinfo,mpi_real8,mpi_complex16&
           & ,mpi_sum,start_time,end_time,write_time
#endif
      USE potential_module, ONLY:cal_veff
      USE mixer_module
      USE energy_module
      USE struct_module , ONLY : ncharge,natom
      USE smearing_module , ONLY : wke,smear_init,destroy_smear,fermilevel
      USE chebyshev_module , ONLY : first_CheSCF_sto_gamma
      use m_time_evaluate, only: filename
      !
      IMPLICIT NONE
      REAL(DP),INTENT(INOUT)  :: rhoS(:,:,:,:)
      REAL(DP),INTENT(OUT) :: psi(:,:,:,:)
      REAL(DP),INTENT(OUT)    :: eval(:,:,:)
      !LOCAL
      ! REAL(DP) :: rhoSd(n1,n2,n3,NSPIN)
      REAL(DP) :: veff(n1,n2,n3,NSPIN),veffd(n1,n2,n3,NSPIN)
      REAL(DP) :: fme,ne,ets
      REAL(DP) :: drho,dtoten,toten,totend,res
      INTEGER(I4B) :: Iter , nev, Iwd
      integer(i4b) :: i
#ifdef MPI
      REAL(DP),allocatable :: rho_1D(:,:),rhod_1D(:,:)
      REAL(DP) :: rho_global(global_n1,global_n2,global_n3)
      REAL(DP) :: psi_global(global_n1,global_n2,global_n3)
      integer(i4b) :: ix,iy,iz
#endif
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !number of state we need
      nev=Nstates
      ne=REAL(ncharge,8)
      !inirial smear
#ifndef MPI
      CALL smear_init(nev)
#else
      CALL smear_init(Nstates_global)
#endif
      ! print*,'sum(rho^2)',sum(rhoS**2),'size(rhoS)',size(rhoS)
      ! open(1111,file='rho_init')
      ! write(1111,*)rhoS
      ! close(1111)
      !initial mixer
      CALL init_mixer_data()
      !Calculate the effective potential
      if(Llastrho)then
         CALL get_rho(natom,struct%poscar,parallel%mygrid_range(3),Nspin,rhoS,dvol)
         CALL get_psi(parallel%mygrid_range(3),nev,Nspin,psi(:,:,1,:))
      CALL cal_veff(rhoS,veff)
      veffd=veff
      else
      CALL cal_veff(rhoS,veff)
      !first step,don't consider spin
      IF(LFIRST)THEN
         !First step by Arpack
         print*,'WARNING: Gamma no arpack now'
         stop
         ! CALL solver_spin(rhoS,psi,eval,nev,100.d0)
      ELSE
         IF(Lrandom)THEN
            !first step by random subspace
            ! CALL first_CheSCF_random(rhoS,nev,psi,eval)
            print *,"only support isolated system"
            stop
         ELSE
            !first step by Rayleigh-Ritz step(STO+PW)
            CALL first_CheSCF_sto_gamma(veff,nev,psi,eval)
            ! print *,"only support isolated system"
         ENDIF
      ENDIF
      ! rhoSd(:,:,:,:)=rhoS(:,:,:,:)
#ifdef MPI
      allocate(rho_1D(n,NSPIN))
      ! do i=1,NSPIN
      !    call rho_trans1D(rhoS(:,:,:,i),rho_1D(:,i))
      ! enddo
      CALL smear_updaterho_gamma(nev,ne,psi,eval,wke,fme,ets,rho_1D)
      do i=1,NSPIN
         call rho_trans3D(rho_1D(:,i),rhoS(:,:,:,i))
      enddo
      deallocate(rho_1D)
#else
      CALL smear_updaterho_gamma(nev,ne,psi,eval,wke,fme,ets,rhoS)
#endif
      ! rhoSd=rhoS
      ! print *,'fme',fme
      CALL totalenergy_gamma(psi,rhoS,eval,fme,ets,.FALSE.)
      toten=Etot
      veffd=veff
      CALL cal_veff(rhoS,veff)
      !first step use simple mixing
      CALL mixing(0,veff,veffd,res)
      endif ! Llastrho}}
      !iter
      DO Iter=1,NMITER
#ifdef MPI
      CALL start_time('scf',.true.)
#endif
         totend=toten
         !call chebyshev filter
         CALL filter_spin_gamma(veff,psi,eval)

         !update rho
#ifdef MPI
         allocate(rho_1D(n,NSPIN))
         ! do i=1,NSPIN
         !    call rho_trans1D(rhoS(:,:,:,i),rho_1D(:,i))
         ! enddo
         CALL smear_updaterho_gamma(nev,ne,psi,eval,wke,fme,ets,rho_1D)
         ! print*,'update rho'
         do i=1,NSPIN
            call rho_trans3D(rho_1D(:,i),rhoS(:,:,:,i))
         enddo
         deallocate(rho_1D)
#else
         CALL smear_updaterho_gamma(nev,ne,psi,eval,wke,fme,ets,rhoS)
#endif

         !> calculate energy
         CALL totalenergy_gamma(psi,rhoS,eval,fme,ets,.TRUE.)
         toten=Etot

         !check
         dtoten=ABS(toten-totend)*hart2ev/natom
#ifndef MPI
         WRITE(6,"(1X,A8,I4,1X,A6,ES18.9,1X,A6,ES15.7,1X,A8,ES15.7)")     &
              &      '>CheFSI:',iter,'Energy',toten*hart2ev,'Res(E)',dtoten ,  &
              &  'Res(rho)',res
#else
         IF(parallel%isroot) WRITE(6,"(1X,A8,I4,1X,A6,ES18.9,1X,A6,ES15.7,1X,A8,ES15.7)")     &
              &      '>CheFSI:',iter,'Energy',toten*hart2ev,'Res(E)',dtoten ,  &
              &  'Res(rho)',res
#endif
         IF(res<RTOL.AND.dtoten<ETOL)THEN
            IWD=IWD+1
            IF(IWD==3) EXIT
         ELSE
            IWD=0
         ENDIF

         !> veff next step
         CALL Cal_Veff(rhoS,veff)
         !mixing potential
         CALL mixing(Iter,veff,veffd,res)
#ifdef MPI
      CALL end_time('scf',.true.)
      if(parallel%isroot)CALL write_time('scf',.true.)
#endif
      ENDDO
      CALL destroy_mixer()
      IF(Iter==NMITER) THEN
         WRITE(6,*) 'SCF:failed STOP'
         !STOP
      ENDIF

#ifdef MPI
      if(parallel%isroot)WRITE(6,*)'SCF DONE!STEP is:',Iter
#else
      WRITE(6,*)'SCF DONE!STEP is:',Iter
#endif
      !call total energy
      CALL totalenergy_gamma(psi,rhoS,eval,fme,ets,.TRUE.)
      !print*,'Etot',Etot
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ENDSUBROUTINE CheFSI_Gamma
   !--------------------smear update  rho---------------------
   SUBROUTINE smear_updaterho_gamma(nev,ne,psi,eval,wke_l,fme,ets,rhoS)
      USE smearing_module
      USE parameters , ONLY : NSPIN,Nstates_global
      USE grid_module , ONLY : n,nsn,dvol,nk,KPT
#ifdef MPI
      USE smpi_math_module, ONLY: MPI_REAL8, MPI_SUM, mpinfo
      use m_time_evaluate, only:filename
#endif
      !
      IMPLICIT NONE
      !INPUT
      INTEGER(I4B),INTENT(IN) :: nev   !number of states
      REAL(DP),INTENT(IN)     :: ne !number of charge
      REAL(DP),INTENT(IN) :: psi(:,:,:,:) !eigen-states
      REAL(DP),INTENT(IN)     :: eval(:,:,:) !eigen-values
      REAL(DP),INTENT(OUT)     :: wke_l(:,:,:)  !weight of k
      REAL(DP),INTENT(OUT) :: ets !energy of electronic entropy
      !OUT PUT
      REAL(DP),INTENT(OUT) :: fme !fermi leval
      REAL(DP),INTENT(OUT) :: rhoS(nsn)
#ifdef MPI
      REAL(DP) ::rhoS_local(nsn)
#endif
      !
      INTEGER(I4B) :: Ik,Ispin,Id,lft,rit !for index
      INTEGER(I4B) :: Iocc !for occupation state index
      REAL(DP) :: prho(n) ! 1D density in each k
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !calculate the weight of k points
#ifndef MPI
      CALL Fermilevel(ne,nev,nk,KPT%wk,eval,wke_l,fme,ets)
#else
      CALL Fermilevel(ne,Nstates_global,nk,KPT%wk,eval,wke_l,fme,ets)
      if(parallel%isroot)then
         ! print*,'ne',ne,'nev',nev
         ! print*,'eval',eval
         ! print*,'fme',fme
         ! print*,'ets',ets
      endif
#endif
      !
      DO Ispin=1,NSPIN
         prho(:)=0.d0
         DO Ik=1,nk
            DO Iocc=1,nev
               !occ number
               IF(ABS(wke_l(Iocc,Ik,Ispin))<1e-14) CYCLE
               !
               DO Id=1,n
#ifndef MPI
                  prho(Id)=prho(Id)+ABS(psi(Id,Iocc,Ik,Ispin))**2*wke_l(Iocc,Ik,Ispin)
#else
                  prho(Id)=prho(Id)+ABS(psi(Id,Iocc,Ik,Ispin))**2&
                       &*wke_l(parallel%sub2sum(Iocc,parallel%ranky+1),Ik,Ispin)
#endif
               ENDDO
            ENDDO
         ENDDO
         !
         lft=(Ispin-1)*n+1
         rit=Ispin*n
         !
         rhoS(lft:rit)=prho(:)
      ENDDO
      !dvol for psi normal
#ifdef MPI
       ! rhoS=rhoS/dvol
       rhoS_local=rhoS/dvol
       CALL MPI_ALLREDUCE(rhoS_local, rhoS, nsn, MPI_REAL8,&
            & MPI_SUM, parallel%commy, mpinfo)
#else
      rhoS=rhoS/dvol
#endif
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ENDSUBROUTINE smear_updaterho_gamma
   !-----------------------filter_spin------------------------
   SUBROUTINE filter_spin_gamma(veff,psi,eval)
      !need to be improve for spin
      USE parameters , ONLY : NSPIN
      USE Grid_module , ONLY : n1,n2,n3, nk
      USE potential_module , ONLY : cal_veff
      USE chebyshev_module , ONLY : cheby_filter_RR_gamma
      IMPLICIT NONE
      REAL(DP),INTENT(IN) :: veff(:,:,:,:)
      REAL(DP),INTENT(INOUT) :: psi(:,:,:,:)
      REAL(DP),INTENT(INOUT) :: eval(:,:,:)
      !LOCAL
      INTEGER(I4B) :: Is,Ik
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! print *,"only support isolated system"
      DO Is=1,NSPIN
      DO Ik=1,nk !> kfilter
         CALL cheby_filter_RR_gamma(Ik,veff(:,:,:,Is),psi(:,:,Ik,Is)&
              &,eval(:,Ik,Is))
      ENDDO
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ENDSUBROUTINE filter_spin_gamma
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDMODULE scf_module

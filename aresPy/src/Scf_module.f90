MODULE scf_module
!############################################################!
!*For :self-consistent-field module                          !
!*Author : Qiang Xu                                          !
!*Date   : 2017-7-20                                         !
!############################################################!
   USE constants
#ifdef MPI
   USE smpi_math_module
#endif
   IMPLICIT NONE
   INTEGER(I4B) :: IWD !windows for test scf
   LOGICAL :: LSCF !test for scf
CONTAINS
   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   SUBROUTINE electronicSCF()
      USE parameters , ONLY : Idiag
      USE grid_module, ONLY : grid,eigen,n
      USE struct_module , ONLY : struct
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ifdef MPI
      IF(parallel%isroot) THEN
#endif
      WRITE(*,*) '>----SCF Iterations for Solving KS Equations----<'
#ifdef MPI
      ENDIF
      CALL MPI_Barrier(parallel%comm,mpinfo)
#endif
      IF(Idiag==0)THEN
         CALL ArpackSCF(n,grid%rhoS,grid%rho,eigen)
      ELSE
         CALL CheFSI(n,grid%rhoS,grid%rho,eigen)
      ENDIF
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE electronicSCF
!############################################################!
!*For :k-space self consistent                               !
!*Author : Qiang Xu                                          !
!*Date   : 2017-07-20                                        !
!############################################################!
   SUBROUTINE ArpackSCF(nps,rhoS,rho,eig)
      USE parameters , ONLY : NMITER,RTOL,ETOL,nev=>Nstates,nev_tot=>Nstates_global
      USE grid_module , ONLY : dvol,eigen_type
      USE mixer_module
      USE energy_module, ONLY : Etot,TotalEnergy
      USE Struct_module , ONLY : ne=>ncharge,natom
      USE smearing_module , ONLY : smear_init,Fermilevel,smear_updaterho
      USE potential_module , ONLY : CalVeff

      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B),INTENT(IN) :: nps
      REAL(DP),INTENT(INOUT)  :: rhoS(nps,nspin),rho(nps)
      TYPE(eigen_type),INTENT(INOUT) :: eig
      !LOCAL
      REAL(DP),DIMENSION(nps,nspin) :: Veff,Veffd
      REAL(DP) :: dEtot,Etotd,Res
      INTEGER(I4B) :: Iter, Is
      REAL(DP),PARAMETER :: diagTOL=5e-5
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IWD=0
      !inirial smear
#ifdef MPI
      CALL smear_init(nev_tot)
#else
      CALL smear_init(nev)
#endif
      !initial mixer
      CALL init_mixer(nps)
      !++++++++++++++++++++++++++++++First diag+++++++++++++++++++++++++++++++
      !First diag
      CALL CalVeff(nps,rhoS,rho,Veff)
      CALL EigenSolver_real(nps,nev,Veff,eig%wvfG(:,:,:),eig%val(:,1,:),100._DP)
      !update rhoS
      CALL smear_updaterho(nps,nev,ne,eig,rhoS,rho)
      CALL TotalEnergy(nps,eig,rhoS,rho)
      !Updatting veff
      Veffd=Veff
      CALL CalVeff(nps,rhoS,rho,Veff)
      !first mixing by simple mixing
      CALL mixing(0,Veff,Veffd,res)
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO Iter=1,NMITER

         Etotd=Etot

         !diag H
         CALL EigenSolver_real(nps,nev,Veff,eig%wvfG(:,:,:),eig%val(:,1,:),diagTOL)

         !update rhoS/Energy
         CALL smear_updaterho(nps,nev,ne,eig,rhoS,rho)
         CALL TotalEnergy(nps,eig,rhoS,rho)
     
         !check for exit
         dEtot=ABS(Etot-Etotd)*hart2ev/natom
         !===============================================================
         WRITE(*,"(1X,A8,I4,1X,A6,ES15.7,1X,A6,ES15.7,1X,A8,ES15.7)")     &
              &      '>Arpack:',iter,'Energy',Etot*hart2ev,'Res(E)',dEtot ,  &
              &  'Res(POT)',res
         !Test for exit
         IF(Res<RTOL.AND.ABS(dEtot)<ETOL)THEN
            IWD=IWD+1
            IF(IWD==3) EXIT
         ELSE
            IWD= 0
         ENDIF
         !===============================================================
         !Updating effective potential
         CALL CalVeff(nps,rhoS,rho,Veff)
         !mixing for update rho and store the old rho
         CALL mixing(Iter,Veff,Veffd,res)

      ENDDO
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF(Iter==NMITER) THEN
         WRITE(*,*) 'SCF:failed STOP'
         STOP
      ENDIF

      WRITE(*,*) 'SCF DONE!STEP is:',Iter
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE ArpackSCF
   !-----------------------solver of Gamma-----------------------
   SUBROUTINE EigenSolver_real(nps,nev,veff,psi,eval,diagTOL)
 
      USE parameters , ONLY : nspin,nev_tot=>Nstates_global
      !arpk
      USE Arpack_module , ONLY : real_diagH_arpack
      USE math, ONLY : rsort_eigen
      IMPLICIT NONE
      REAL(DP),INTENT(IN) :: veff(nps,nspin)
      REAL(DP),INTENT(IN) :: diagTOL !for diag 
      REAL(DP),INTENT(OUT) :: psi(nps,nev,nspin)
      REAL(DP),INTENT(OUT) :: eval(nev,nspin)
      INTEGER,INTENT(IN) :: nps,nev !num. of points and states
      !LOCAL
      INTEGER(I4B) :: Is,Ik
      !Arpk
      INTEGER(I4B) :: nec,isc,info
      REAL(DP)     :: resid_int(nps)
      INTEGER(I4B) :: maxmvs
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !
      DO Is=1,nspin
         !initialize arpk
         info=0
         resid_int(:)=0.d0
         diagH:DO isc=1,5
                   maxmvs=isc*30000
                   CALL real_diagH_arpack(nps,veff(:,Is),nev,psi(:,:,Is),eval(:,Is), &
                      &    resid_int,nec,info,maxmvs,diagTOL)
                   IF(nec>=nev)THEN
                      EXIT diagH
                   ELSE
                      WRITE(*,*)'diag(H) failed,restarting'
                      !STOP          !I don't want to do this
                      !info=0
                   ENDIF

              ENDDO diagH

         IF(isc>=5)THEN
            WRITE(*,*) 'diag(H) failed,STOP,increase NADST and try again'
            STOP
         ENDIF

         IF(nev>1)THEN
             CALL rsort_eigen(nev,eval(:,Is),psi(:,:,Is))
         ENDIF
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE EigenSolver_real
!############################################################!
!*For :k-space Chebyshev Filtering                           !
!*Author : Qiang Xu                                          !
!*Date   : 2017-11-28                                        !
!############################################################!
   !-------------------------CheFSI---------------------------
   SUBROUTINE CheFSI(nps,rhoS,rho,eig)
      USE parameters , ONLY : NMITER,RTOL,ETOL,nspin &
             & ,nev=>Nstates,nev_tot=>Nstates_global &
             & ,LWAVE
      USE grid_module , ONLY : dvol,eigen_type
      USE mixer_module
      USE energy_module, ONLY : Etot,TotalEnergy
      USE Struct_module , ONLY : ne=>ncharge,natom !,struct
      USE smearing_module , ONLY : smear_init,Fermilevel,smear_updaterho
      USE potential_module , ONLY : CalVeff
      USE chebyshev_module , ONLY : BuildSubspace
#ifdef MPI
      USE smpi_math_module
#endif
      !
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN) :: nps
      REAL(DP),INTENT(INOUT)  :: rhoS(nps,nspin),rho(nps)
      TYPE(eigen_type),INTENT(INOUT) :: eig
      !LOCAL
      REAL(DP) :: dEtot,Etotd,res
      REAL(DP) :: veff(nps,nspin),veffd(nps,nspin)
      INTEGER(I4B) :: Iter
!REAL(DP) :: tmp,tmp_loc
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !CALL MPI_barrier
      IWD=0
      !inirial smear
#ifdef MPI
      CALL smear_init(nev_tot)
#else
      CALL smear_init(nev)
#endif
      !initial mixer
      CALL init_mixer(nps)
      !first step by Rayleigh-Ritz step(STO+random)
      CALL CalVeff(nps,rhoS,rho,Veff)
      !initialize the subspace
      CALL BuildSubspace(nps,nev,veff,eig)
      !updating density and energy
      CALL smear_updaterho(nps,nev,ne,eig,rhoS,rho)
!tmp_loc=SUM(rhoS*Veff)*dvol
!CALL MPI_ALLREDUCE(tmp_loc,tmp,1,MPI_REAL8,MPI_SUM,parallel%comm,mpinfo)
!IF(parallel%isroot) print*,tmp
!STOP
      !evaluate the total energy
      CALL TotalEnergy(nps,eig,rhoS,rho)
      !store Veff
      Veffd=Veff
      CALL CalVeff(nps,rhoS,rho,Veff)
!tmp_loc=SUM(Veff)
!CALL MPI_ALLREDUCE(tmp_loc,tmp,1,MPI_REAL8,MPI_SUM,parallel%comm,mpinfo)
!IF(parallel%isroot) print*,tmp
      !first step use simple mixing
      CALL mixing(0,veff,veffd,res)
!tmp_loc=SUM(rhoS*Veff)*dvol
!CALL MPI_ALLREDUCE(tmp_loc,tmp,1,MPI_REAL8,MPI_SUM,parallel%comm,mpinfo)
!IF(parallel%isroot) print*,tmp
      !iteration
      DO Iter=1,NMITER
#ifdef MPI
      CALL start_time('scf',.true.)
#endif
         Etotd=Etot
         !call chebyshev filter
         CALL filter_spin_Gamma(nps,nev,veff,eig%wvfG,eig%val(:,1,:))
!IF(parallel%isroot) print*,eig%val
!STOP
         !update rho
         CALL smear_updaterho(nps,nev,ne,eig,rhoS,rho)
         CALL TotalEnergy(nps,eig,rhoS,rho)
         !check
         dEtot=ABS(Etot-Etotd)*hart2ev/natom
!print*,'ne',ne
!pause
         !===============================================================
#ifdef MPI
         IF(parallel%isroot) THEN
#endif
         WRITE(*,"(1X,A8,I4,1X,A6,ES15.7,1X,A6,ES15.7,1X,A8,ES15.7)")     &
              &      '>CheFSI:',iter,'Energy',Etot*hart2ev,'Res(E)',dEtot ,  &
              &  'Res(POT)',res
#ifdef MPI
         ENDIF
#endif
!print*,'TEST',dEtot,res,RTOL,ETOL
         !Test for exit
         IF(Res<RTOL.AND.ABS(dEtot)<ETOL)THEN
            IWD=IWD+1
            IF(IWD==3) EXIT
         ELSE
            IWD= 0
         ENDIF
         !===============================================================
         CALL CalVeff(nps,rhoS,rho,veff)
         !mixing potential
         CALL mixing(Iter,veff,veffd,res)
#ifdef MPI
      CALL end_time('scf',.true.)
      ! if(parallel%isroot)CALL write_time('scf',.true.)
#endif
      ENDDO
!PRINT'111111111'
#ifdef MPI
      IF(parallel%isroot)THEN
#endif
      IF(Iter==NMITER) THEN
         WRITE(*,*) 'SCF:failed STOP'
         !STOP
      ENDIF

      WRITE(*,*) 'SCF DONE!STEP is:',Iter
#ifdef MPI
      ENDIF
#endif
!PRINT'2222222222'
      !Re-Diag(H) for outputing wave functionals
      !IF(LWAVE) &
      ! &  CALL EigenSolver_real(nps,nev,Veff,eig%vecG,eig%val(:,1,:),5e-5)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE CheFSI
   !-------------------------CheFSI---------------------------
   SUBROUTINE filter_spin_Gamma(nps,nev,veff,X,D)
      !need to be improve for spin
      USE parameters , ONLY : nspin
      USE chebyshev_module , ONLY : cheby_filtering_GRRr
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN) :: nps,nev
      REAL(DP),INTENT(IN) :: veff(nps,nspin)
      REAL(DP),INTENT(INOUT) :: X(nps,nev,nspin)
      REAL(DP),INTENT(INOUT) :: D(:,:)
      !LOCAL
      INTEGER(I4B) :: Is
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      DO Is=1,nspin
         CALL cheby_filtering_GRRr(nps,nev,veff(:,Is),X(:,:,Is),D(:,Is))
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE filter_spin_Gamma
   !-------------------------kfilter--------------------------
ENDMODULE scf_module

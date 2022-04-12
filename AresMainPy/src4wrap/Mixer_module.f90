MODULE mixer_module
!###################################################!
!*For   : Mixing the charge density to convergence  !
!*Author: Qiang Xu                                  !
!*Data  : 2017/07/30                                !
!*IMIXER:   0   For simple mixing + Anderson mixing !
!         other For simple mixing + (r)Pulay mixing !
!###################################################!
   USE constants
   USE parameters , ONLY : IMIXER
   IMPLICIT NONE
   TYPE mixer_data
        !Basic data (Simple Mixing and Anderson,real space)
        INTEGER :: &
          NHMIX,NHMIN, & ! Number of history vector for Anderson mixing
          NAM, &      ! Dimensions of the vector for Anderson mixing
          NITRA       ! Number of simple mixing method
        REAL(DP) :: &
          alpha, &    ! simple mixing factor
          beta, &     ! non-simple mixing factor
          w0          ! for ill matrix
        REAL(DP), dimension(:,:), allocatable :: &
          DXL, &      !delta vector
          DFL, &      !delta rediuse vector
          VOMA
        REAL(DP) :: SP !factor for scalar product

        !For Kerker mixing (recip-space)
        REAL(DP),ALLOCATABLE :: kerker(:)
        COMPLEX(DP),ALLOCATABLE ::  & 
                              &    DXGL(:,:) & !input vector
                              &,   DFGL(:,:)   !residue
        REAL(DP) :: BMIX &  !For  G^2 Cut
                & , AMIX   !For minimial kerker
   ENDTYPE mixer_data
   !data for mixer
   TYPE(mixer_data) :: mixer
CONTAINS
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!-------------------------------------------------------------------------
  SUBROUTINE init_mixer_data
     USE parameters, ONLY: Lpbc
     implicit none

     if(Lpbc)then
        call init_mixer_data_per()
     else
        call init_mixer_data_iso()
     endif

   END SUBROUTINE init_mixer_data
  !-------------------init_mixer_data_per---------------------
   SUBROUTINE init_mixer_data_per()!{{{
      USE grid_module , ONLY : nsn,dvol,ng,ng1,ng2,ng3
      USE parameters , ONLY : NHMIX,NHMIN,MALPHA,MBETA,NSMIX,W0AM,Nspin,AMIX,BMIX
#ifdef MPI
      USE smpi_math_module, ONLY : parallel
      USE grid_module , ONLY :n1,n2,n3
#endif
      USE m_time_evaluate, ONLY: memory_sum
      IMPLICIT NONE
      !
      INTEGER(I4B) :: NUH,NAM
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !Basisc data
      !first times for simpling mixing
      mixer%NITRA = NSMIX
      !coe for simple mixing
      mixer%alpha = MALPHA
      !coe for non-simple mixing
      mixer%beta = MBETA
      !number of history for non-simpling mixing
      mixer%NHMIX = ABS(NHMIX)
      mixer%NHMIN = NHMIN
      NUH = mixer%NHMIX
      !for ill-matrix
      mixer%w0 = W0AM
      !
      !dimension for mixing vector
      IF(IMIXER==0)THEN
#ifdef MPI
         mixer%NAM = n1*n2*n3*Nspin
#else
         mixer%NAM = nsn
#endif
      ELSE
#ifdef MPI
         mixer%NAM = nspin*ng1*ng2*parallel%local_z
#else
         mixer%NAM=ng*Nspin
#endif
      ENDIF
      NAM = mixer%NAM

      CALL destroy_mixer()
      !
      SELECT CASE(IMIXER)
      CASE(0)
#ifdef MPI
         if(parallel%isroot)print*,'[ Mixing ] Simple + rPulay mixing'
#else
         print*,'[ Mixing ] Simple + rPulay mixing'
#endif
         !real space mixing
         ALLOCATE(mixer%DXL(NAM,NUH))
         ALLOCATE(mixer%DFL(NAM,NUH))
         ALLOCATE(mixer%VOMA(NUH,NUH))
         call memory_sum('mixing1',real(size(mixer%DXL),DP)*DP&
              &+size(mixer%DFL)*DP+size(mixer%VOMA)*DP)
         mixer%DXL = 0.d0
         mixer%DFL = 0.d0
         mixer%VOMA = 0.d0
         !WEIGHTS FOR SCALAR PRODUCT
         mixer%SP=dvol
      CASE(2)
#ifdef MPI
         if(parallel%isroot)print*,'[ Mixing ] Simple + resta mixing'
#else
         print*,'[ Mixing ] Simple + resta mixing'
#endif
         !recip-space mixing
#ifdef MPI
         NAM=nspin*ng1*ng2*parallel%local_z
         ALLOCATE(mixer%kerker(ng1*ng2*parallel%local_z))
#else
         NAM=nspin*ng
         ALLOCATE(mixer%kerker(ng))
#endif
         ALLOCATE(mixer%DXGL(NAM,NUH))
         ALLOCATE(mixer%DFGL(NAM,NUH))
         call memory_sum('mixing2',real(size(mixer%kerker),DP)*DP&
              &+size(mixer%DXGL)*DP+size(mixer%VOMA)*DP)
         mixer%AMIX=AMIX
         mixer%BMIX=BMIX
         mixer%DXGL(:,:)=0.d0
         mixer%DFGL(:,:)=0.d0
         CALL init_resta()
      CASE default
#ifdef MPI
         if(parallel%isroot)print*,'[ Mixing ] Simple + rPulayK mixing'
#else
         print*,'[ Mixing ] Simple + rPulayK mixing'
#endif
         !recip-space mixing
#ifdef MPI
         NAM=nspin*ng1*ng2*parallel%local_z
         ALLOCATE(mixer%kerker(ng1*ng2*parallel%local_z))
#else
         NAM=nspin*ng
         ALLOCATE(mixer%kerker(ng))
#endif
         ALLOCATE(mixer%DXGL(NAM,NUH))
         ALLOCATE(mixer%DFGL(NAM,NUH))
         call memory_sum('mixing0',real(size(mixer%DXL),DP)*DP&
              &+size(mixer%DFL)*DP)
         mixer%DXGL = 0.d0
         mixer%DFGL = 0.d0
         mixer%AMIX=AMIX
         mixer%BMIX=BMIX
         CALL init_kerker()
      ENDSELECT
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE init_mixer_data_per!}}}
  !-------------------init_mixer_data_iso---------------------
   SUBROUTINE init_mixer_data_iso()!{{{
      USE grid_module , ONLY : nsn,dvol,ng,ng1,ng2,ng3
      USE parameters , ONLY : NHMIX,NHMIN,MALPHA,MBETA,NSMIX,W0AM,Nspin,AMIX,BMIX
#ifdef MPI
      USE smpi_math_module, ONLY : parallel
#endif
      USE m_time_evaluate, ONLY: memory_sum
      IMPLICIT NONE
      !
      INTEGER(I4B) :: NUH,NAM
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !Basisc data
      !first times for simpling mixing
      mixer%NITRA = NSMIX
      !coe for simple mixing
      mixer%alpha = MALPHA
      !coe for non-simple mixing
      mixer%beta = MBETA
      !number of history for non-simpling mixing
      mixer%NHMIX = ABS(NHMIX)
      mixer%NHMIN = NHMIN
      NUH = mixer%NHMIX
      !for ill-matrix
      mixer%w0 = W0AM
      !
      !dimension for mixing vector
      IF(IMIXER==2)THEN
         mixer%NAM=ng*Nspin
      ELSE
#ifndef MPI
         mixer%NAM = nsn
#else
         mixer%NAM = parallel%mygrid_range(3)*Nspin
#endif
      ENDIF
      NAM = mixer%NAM

      CALL destroy_mixer()
      !
      SELECT CASE(IMIXER)
      CASE(1)
         !real space mixing
         ALLOCATE(mixer%DXL(NAM,NUH))
         ALLOCATE(mixer%DFL(NAM,NUH))
         ALLOCATE(mixer%VOMA(NUH,NUH))
         call memory_sum('mixing1',real(size(mixer%DXL),DP)*DP&
              &+size(mixer%DFL)*DP+size(mixer%VOMA)*DP)
         mixer%DXL = 0.d0
         mixer%DFL = 0.d0
         mixer%VOMA = 0.d0
         !WEIGHTS FOR SCALAR PRODUCT
         mixer%SP=dvol
      CASE(2)
         !recip-space mixing
         ALLOCATE(mixer%kerker(ng))
         ALLOCATE(mixer%DXGL(NAM,NUH))
         ALLOCATE(mixer%DFGL(NAM,NUH))
         call memory_sum('mixing2',real(size(mixer%kerker),DP)*DP&
              &+size(mixer%DXGL)*DP+size(mixer%VOMA)*DP)
         mixer%AMIX=AMIX
         mixer%BMIX=BMIX
         mixer%DXGL(:,:)=0.d0
         mixer%DFGL(:,:)=0.d0
         CALL init_kerker()
      CASE default
         !real-space mixing
         ALLOCATE(mixer%DXL(NAM,NUH))
         ALLOCATE(mixer%DFL(NAM,NUH))
         call memory_sum('mixing0',real(size(mixer%DXL),DP)*DP&
              &+size(mixer%DFL)*DP)
         mixer%DXL = 0.d0
         mixer%DFL = 0.d0
         mixer%AMIX=0.4d0
         mixer%BMIX=1.d0
      ENDSELECT
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE init_mixer_data_iso!}}}
!-----------------------------PARTING LINE---------------------------------
   SUBROUTINE destroy_mixer()!{{{
     USE m_time_evaluate, ONLY: memory_free
      !
      IMPLICIT NONE
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !real space mixing
      IF(ALLOCATED(mixer%DXL))THEN
         DEALLOCATE(mixer%DXL)
         CALL memory_free('mixing_DXL',real(size(mixer%DXL),DP)*DP)
      ENDIF
      IF(ALLOCATED(mixer%DFL))THEN
         DEALLOCATE(mixer%DFL)
         CALL memory_free('mixing_DFL',real(size(mixer%DFL),DP)*DP)
      ENDIF
      IF(ALLOCATED(mixer%VOMA))THEN
         DEALLOCATE(mixer%VOMA)
         CALL memory_free('mixing_VOMA',real(size(mixer%VOMA),DP)*DP)
      ENDIF
      !recip-space mixing
      IF(ALLOCATED(mixer%kerker))THEN
         DEALLOCATE(mixer%kerker)
         CALL memory_free('mixing_kerker',real(size(mixer%kerker),DP)*DP)
      ENDIF
      IF(ALLOCATED(mixer%DXGL))THEN
         DEALLOCATE(mixer%DXGL)
         CALL memory_free('mixing_DXGL',real(size(mixer%DXGL),DP)*DP)
      ENDIF
      IF(ALLOCATED(mixer%DFGL))THEN
         DEALLOCATE(mixer%DFGL)
         CALL memory_free('mixing_DFGL',real(size(mixer%DFGL),DP)*DP)
      ENDIF
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE destroy_mixer!}}}
!--------------------------------mixing------------------------------------
   SUBROUTINE mixing(Iter,XOUT,XIN,Res)!{{{
      USE parameters , ONLY : Nspin
#ifdef MPI
     USE smpi_math_module , ONLY : parallel,mpi_real8,mpi_sum,mpinfo
     USE grid_module , ONLY : dvol,n,ng1,ng2,fft_sph
#else
      USE grid_module , ONLY : ng1,ng2,ng3,ng,dvol,nsn
#endif
      USE struct_module , ONLY : ncharge
      USE FOURIER
      USE m_time_evaluate, ONLY: memory_sum,memory_free
      IMPLICIT NONE
      !
      INTEGER,INTENT(IN) :: Iter !present iteration step
      REAL(DP),INTENT(INOUT) :: XIN(:,:,:,:),XOUT(:,:,:,:)
      REAL(DP),INTENT(OUT) :: Res !Norm of residual of density
      !LOCAL
      COMPLEX(DP),DIMENSION(:,:,:,:),ALLOCATABLE :: &
                           & XINg,RLg
#ifdef MPI
     REAL(DP) :: res_local
#endif
      INTEGER(I4B) :: Is
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF(Iter<1)THEN
         res=-1.00
         XOUT=XIN
         RETURN
      ENDIF
      SELECT CASE(IMIXER)
      CASE(0)
         !real-space mixing
         CALL Anderson_mixing(Iter,XOUT,XIN,Res)
         ! res=res/ncharge
      CASE(2)
         !recip-space mixing
         !allocate
#ifdef MPI
         ALLOCATE(XINg(ng1,ng2,parallel%local_z,Nspin)&
              &,RLg(ng1,ng2,parallel%local_z,Nspin))
#else
         ALLOCATE(XINg(ng1,ng2,ng3,Nspin),RLg(ng1,ng2,ng3,Nspin))
#endif
         call memory_sum('mixing_recip',real(size(XINg),DP)*2*DP+size(RLg)*2*DP)
         !store the res
         XOUT(:,:,:,:)=XOUT(:,:,:,:)-XIN(:,:,:,:)
         !residual
#ifdef MPI
        res_local=SUM(XOUT*XOUT)
        CALL MPI_ALLREDUCE(res_local,res,1,MPI_REAL8,MPI_SUM,parallel%commx,mpinfo)
#else
         res=SUM(XOUT*XOUT)
#endif
         res=SQRT(res)
         !to recip-space
         DO Is=1,Nspin
#ifdef MPI
            XINg(:,:,:,Is)  = FFT_sph(XIN(:,:,:,Is))
            RLg(:,:,:,Is)   = FFT_sph(XOUT(:,:,:,Is))
#else
            XINg(:,:,:,Is)  = FFT(XIN(:,:,:,Is))
            RLg(:,:,:,Is)   = FFT(XOUT(:,:,:,Is))
#endif
         ENDDO
!        WRITE(6,'(9x,I2,1X,A6,F18.10,3x,F18.10)') iter,'RLG   ', real(rlg(2,2,1,1)),aimag(rlg(2,2,1,1))
!        WRITE(6,'(9x,I2,1X,A6,F18.10,3x,F18.10)') iter,'XING  ', real(xing(2,2,1,1)),aimag(xing(2,2,1,1))
         CALL resta_mixing(Iter,RLg,XINg)
         !to real-space
         DO Is=1,Nspin
#ifdef MPI
            XIN(:,:,:,Is)  = FFT_sph(XINg(:,:,:,Is))
#else
            XIN(:,:,:,Is)  = FFT(XINg(:,:,:,Is))
#endif
         ENDDO
         XOUT(:,:,:,:)=XIN(:,:,:,:)
         !scale
         call memory_free('mixing_recip',real(size(XINg),DP)*2*DP+size(RLg)*2*DP)
         DEALLOCATE(XINg,RLg)
      CASE default
         !recip-space mixing
         !allocate
#ifdef MPI
         ALLOCATE(XINg(ng1,ng2,parallel%local_z,Nspin)&
              &,RLg(ng1,ng2,parallel%local_z,Nspin))
#else
         ALLOCATE(XINg(ng1,ng2,ng3,Nspin),RLg(ng1,ng2,ng3,Nspin))
#endif
         call memory_sum('mixing_recip',real(size(XINg),DP)*2*DP+size(RLg)*2*DP)
         !store the res
         XOUT(:,:,:,:)=XOUT(:,:,:,:)-XIN(:,:,:,:)
         !residual
#ifdef MPI
        res_local=SUM(XOUT*XOUT)
        CALL MPI_ALLREDUCE(res_local,res,1,MPI_REAL8,MPI_SUM,parallel%commx,mpinfo)
#else
         res=SUM(XOUT*XOUT)
#endif
         res=SQRT(res)
         !to recip-space
         DO Is=1,Nspin
#ifdef MPI
            XINg(:,:,:,Is)  = FFT_sph(XIN(:,:,:,Is))
            RLg(:,:,:,Is)   = FFT_sph(XOUT(:,:,:,Is))
#else
            XINg(:,:,:,Is)  = FFT(XIN(:,:,:,Is))
            RLg(:,:,:,Is)   = FFT(XOUT(:,:,:,Is))
#endif
         ENDDO
         CALL rPulayK_mixing(Iter,RLg,XINg)
         !to real-space
         DO Is=1,Nspin
#ifdef MPI
            XIN(:,:,:,Is)  = FFT_sph(XINg(:,:,:,Is))
#else
            XIN(:,:,:,Is)  = FFT(XINg(:,:,:,Is))
#endif
         ENDDO
         XOUT(:,:,:,:)=XIN(:,:,:,:)
         !scale
         call memory_free('mixing_recip',real(size(XINg),DP)*2*DP+size(RLg)*2*DP)
         DEALLOCATE(XINg,RLg)
      ENDSELECT
      res=res/ncharge
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE mixing!}}}
!--------------------------------mixing------------------------------------
   SUBROUTINE mixing_iso(Iter,XOUT,XIN,Res)!{{{
     USE parameters , ONLY : Nspin
     USE grid_module , ONLY : ng1,ng2,ng3,ng,dvol,nsn
     USE struct_module , ONLY : ncharge
     USE FOURIER
     USE m_time_evaluate, ONLY: memory_sum,memory_free
     IMPLICIT NONE
     !
     INTEGER,INTENT(IN) :: Iter !present iteration step
#ifndef MPI
     REAL(DP),INTENT(INOUT) :: XIN(:,:,:,:),XOUT(:,:,:,:)
#else
     REAL(DP),INTENT(INOUT) :: XIN(:,:),XOUT(:,:)
#endif
     REAL(DP),INTENT(OUT) :: Res !Norm of residual of density
     !LOCAL
     COMPLEX(DP),DIMENSION(:,:,:,:),ALLOCATABLE :: &
          & XINg,RLg
     INTEGER(I4B) :: Is
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     IF(Iter<1)THEN
        res=-1.00
        XOUT=XIN
        RETURN
     ENDIF
     SELECT CASE(IMIXER)
     CASE(1)
        !real-space mixing
        CALL Anderson_mixing(Iter,XOUT,XIN,Res)
        res=res/ncharge
#ifndef MPI
     CASE(2)
        !recip-space mixing
        !allocate
        ALLOCATE(XINg(ng1,ng2,ng3,Nspin),RLg(ng1,ng2,ng3,Nspin))
        call memory_sum('mixing_recip',real(size(XINg),DP)*2*DP+size(RLg)*2*DP)
        !store the res
        XOUT(:,:,:,:)=XOUT(:,:,:,:)-XIN(:,:,:,:)
        !residual
        res=SUM(XOUT*XOUT)/ncharge
        !to recip-space
        DO Is=1,Nspin
           XINg(:,:,:,Is)  = FFT(XIN(:,:,:,Is))
           RLg(:,:,:,Is)   = FFT(XOUT(:,:,:,Is))
        ENDDO
        !        WRITE(6,'(9x,I2,1X,A6,F18.10,3x,F18.10)') iter,'RLG   ', real(rlg(2,2,1,1)),aimag(rlg(2,2,1,1))
        !        WRITE(6,'(9x,I2,1X,A6,F18.10,3x,F18.10)') iter,'XING  ', real(xing(2,2,1,1)),aimag(xing(2,2,1,1))
        CALL rPulayK_mixing(Iter,RLg,XINg)
        !to real-space
        DO Is=1,Nspin
           XIN(:,:,:,Is)  = FFT(XINg(:,:,:,Is))
        ENDDO
        XOUT(:,:,:,:)=XIN(:,:,:,:)
        !scale
        call memory_free('mixing_recip',real(size(XINg),DP)*2*DP+size(RLg)*2*DP)
        DEALLOCATE(XINg,RLg)
#endif
        ! CASE(3)
        !    !recip-space mixing
        !    !allocate
        !    ALLOCATE(XINg(ng1,ng2,ng3,Nspin),RLg(ng1,ng2,ng3,Nspin))
        !    ALLOCATE(XIN_grid(n1,n2,n3,Nspin),XOUT_grid(n1,n2,n3,Nspin))
        !    !> resize the rho
        !    do Is=1,Nspin
        !    do in=1,rho_calc%OneDLength,1
        !       XIN_grid(rho_calc%x(in),rho_calc%y(in),rho_calc%z(in),Is)=XIN()
        !    !store the res
        !    XOUT(:,:,:,:)=XOUT(:,:,:,:)-XIN(:,:,:,:)
        !    !residual
        !    res=SUM(XOUT*XOUT)/ncharge
        !    !to recip-space
        !       DO Is=1,Nspin
        !          XINg(:,:,:,Is) = FFT(XIN_grid(:,:,:,Is))
        !          RLg(:,:,:,Is)   = FFT(XOUT_grid(:,:,:,Is))
        !       ENDDO
        !    CALL rPulayK_mixing(Iter,RLg,XINg)
        !    !to real-space
        !       DO Is=1,Nspin
        !          XINg(:,:,:,Is) = FFT(XIN_grid(:,:,:,Is))
        !          RLg(:,:,:,Is)   = FFT(XOUTi_grid(:,:,:,Is))
        !       ENDDO
        !    !scale
        !    DEALLOCATE(XINg,RLg)
     CASE default
        !real-space r-pulay mixing
        CALL rPulay_mixing(Iter,XOUT,XIN,Res)
        res=res/ncharge
     ENDSELECT
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE mixing_iso!}}}
!########################################################!
!*For    : Simple Mixing + Anderson mixing               !
!########################################################!
!-----------------------------PARTING LINE----------------------------------
   SUBROUTINE Anderson_mixing(IITER,XOUT,XIN,err)!{{{
      USE grid_module , ONLY : nsn,n1,n2
      USE parameters , ONLY : NHMIX,Lpbc
#ifdef MPI
      USE smpi_math_module , ONLY : parallel,mpi_real8,mpi_sum,mpinfo
      use m_time_evaluate, ONLY:filename
#endif
      !
      IMPLICIT NONE
      !
      INTEGER,INTENT(IN) :: IITER !present iteration step
#ifndef MPI
      REAL(DP),INTENT(INOUT) :: XIN(nsn),XOUT(nsn) !
#else
      REAL(DP),INTENT(INOUT) :: XIN(mixer%NAM),XOUT(mixer%NAM) !mixer%NAM
#endif
      REAL(DP),INTENT(OUT) :: err
      !
      REAL(DP) :: alpha, beta !factor for mixing
      REAL(DP) :: XL(mixer%NAM) , FL(mixer%NAM), XN(mixer%NAM)
      integer :: &
        NAM, NUHMAX,NUHMIN,NITRA, &
        J, I1, I2, I3, &
        IH, JH, &
        NUH1,RITER,DH
      REAL(DP) :: dnrm2
#ifdef MPI
      REAL(DP) :: err_local
      ! !> kahan sum
      ! REAL(DP) ::kahan_y,kahan_t,kahan_eps
      ! INTEGER(I4B) :: i,ix
      ! REAL(DP) :: FL_global(global_n)
#endif
      REAL(DP),external :: ddot
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! kahan_eps=0.d0
      NAM   = mixer%NAM
      NUHMAX = mixer%NHMIX
      NUHMIN = mixer%NHMIN
      NITRA = mixer%NITRA
      alpha = mixer%alpha
      beta  = mixer%beta
      !>>>input rho
      XL=XIN
      FL=XOUT-XIN
      !<<<input rho
! #ifdef MPI
!       ix=mod(parallel%mygrid_range(1)-1,n1*n2)+1
!       print*,'size(XIN)',size(XIN),parallel%mygrid_range(3),ix
!       call MPI_ALLGATHERV(FL(ix),parallel%mygrid_range(3)&
!            &,MPI_REAL8,FL_global,parallel%recvcounts&
!            &,parallel%displs,MPI_REAL8&
!            & ,parallel%commx,mpinfo)
!       if(parallel%myid==parallel%numprocs-1)then
!          open(1212,file='Fl'//filename(10:11))
!          write(1212,*)FL_global
!          close(1212)
!       endif
!       print*,'sum FL_global',sum(FL_global**2)
!       call MPI_ALLGATHERV(XIN(ix),parallel%mygrid_range(3)&
!            &,MPI_REAL8,FL_global,parallel%recvcounts&
!            &,parallel%displs,MPI_REAL8&
!            & ,parallel%commx,mpinfo)
!       if(parallel%myid==parallel%numprocs-1)then
!          open(1212,file='Xin'//filename(10:11))
!          write(1212,*)FL_global
!          close(1212)
!       endif
!       print*,'sum XIN_global',sum(FL_global**2)
!       call MPI_ALLGATHERV(XOUT(ix),parallel%mygrid_range(3)&
!            &,MPI_REAL8,FL_global,parallel%recvcounts&
!            &,parallel%displs,MPI_REAL8&
!            & ,parallel%commx,mpinfo)
!       if(parallel%myid==parallel%numprocs-1)then
!          open(1212,file='Xout'//filename(10:11))
!          write(1212,*)FL_global
!          close(1212)
!       endif
!       print*,'sum XOUT_global',sum(FL_global**2)
! #endif
#ifndef MPI
      err=ddot(size(FL),FL,1,FL,1)
#else
      err_local=ddot(size(FL),FL,1,FL,1)
      ! do i=1,size(FL),1
      !    kahan_y=FL(i)*FL(i)-kahan_eps
      !    kahan_t=err_local+kahan_y
      !    kahan_eps=(kahan_t-err_local)-kahan_y
      !    err_local=kahan_t
      ! enddo
      call mpi_allreduce(err_local, err, 1, mpi_real8, mpi_sum, parallel%commx, mpinfo)
#endif
      err=sqrt(err)

      IF(IITER>1) THEN
          mixer%DXL(:,1)=mixer%DXL(:,1)-XL(:)
          mixer%DFL(:,1)=mixer%DFL(:,1)-FL(:)
          CALL OM1C(NAM, NUHMAX, mixer%SP, mixer%DFL, mixer%VOMA)
      END IF

      IF(IITER<=NITRA.OR.IITER<=2) THEN
          XN=XL+alpha*FL
      ELSE
          DH=NUHMAX-NUHMIN+1
          IF(NHMIX>0.OR.(IITER-1)<=NUHMIN)THEN
             !Anderson mixing
             NUH1=MIN(NUHMAX,IITER-1)
          ELSE
             !restart Anderson mixing
             RITER=IITER-NUHMIN-1
             NUH1=MOD(RITER,DH)+NUHMIN
          ENDIF
          !print*,'Iter,NUH1',IITER,NUH1,NUHMAX

          !call Anderson
          !IF(NUH1==0)THEN
          !ELSE
          CALL AMST(beta,mixer%w0,NAM,NUH1,&
                    mixer%DXL,mixer%DFL,mixer%SP,&
                    XL,FL,mixer%VOMA,XN)
          !ENDIF
      ENDIF

      IF(IITER>0)THEN

         DO JH=NUHMAX,2,-1
         DO IH=JH,NUHMAX
            mixer%VOMA(IH,JH)=mixer%VOMA(IH-1,JH-1) 
         END DO
         END DO

         !ONE STEP BACKWARD
         DO IH=NUHMAX,2,-1
            mixer%DXL(:,IH)=mixer%DXL(:,IH-1)
            mixer%DFL(:,IH)=mixer%DFL(:,IH-1)
         END DO
         mixer%DXL(:,1)=XL
         mixer%DFL(:,1)=FL

      ENDIF
      !>>>output rho
      XIN(:)=XN(:)
      XOUT(:)=XN(:)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE Anderson_mixing!}}}
!-----------------------PARTING LINE---------------------------
   SUBROUTINE OM1C(NAM,NUH,SP,DFP,VOMA)!{{{

#ifdef MPI
      USE smpi_math_module , ONLY : parallel,mpi_real8,mpi_sum,mpinfo
#endif
       IMPLICIT NONE

       INTEGER, INTENT(IN) :: & 
          NAM, NUH

       !REAL*8, DIMENSION(:), INTENT(IN) :: & 
       !   SP
       REAL(8),INTENT(IN) :: SP
 
       REAL(8), DIMENSION(:,:), INTENT(IN) :: & 
          DFP
  
       REAL(8), DIMENSION(:,:), INTENT(OUT) :: &
          VOMA

!C*************************************************
!C     THE FIRST COLUMN OF OVERLAP MATRIX FOR
!C     THE ANDERSON MIXING PROCEDURE
!C*************************************************
!C     NAM - NUMBER OF VARIABLES
!C     NUH - NUMBER OF PREVIOUS VECTORS
!C-------------------------------------------------
!C   ON INPUT:
!C     SP(I), I=1,...,NAM,  -  WEIGHTS FOR THE
!C                             SCALAR PRODUCT
!C     DFP(I,IH), I=1,...,NAM, IH=1,...,NUH,  -
!C              - DIFFERENCES OF PREVIOUS VECTORS F
!C
!C   ON OUTPUT:
!C     VOMA(IH,1),  IH=1,...,NUH,  -
!C                 - 1ST COLUMN OF OVERLAP MATRIX
!C*************************************************

       REAL(8), DIMENSION(NAM) :: & 
          VECT

#ifdef MPI
       REAL(8), DIMENSION(NUH) :: VOMA_local
       ! !> kahan sum
       ! REAL(DP) ::kahan_y,kahan_t,kahan_eps
#endif

       INTEGER :: & 
          IH, I

! #ifdef MPI
!        kahan_eps=0.d0
! #endif
      DO 310 IH=1, NUH
        VOMA(IH,1)=0.d0
#ifdef MPI
        VOMA_local(IH)=0.d0
#endif
310   CONTINUE

      DO 312 I=1,NAM
        !VECT(I) = SP(I) * DFP(I,1)
        VECT(I) = SP * DFP(I,1)
312   CONTINUE

      DO 320 IH=1,NUH
      DO 325 I=1,NAM
#ifndef MPI
         VOMA(IH,1) = VOMA(IH,1) + DFP(I,IH) * VECT(I)
#else
         VOMA_local(IH) = VOMA_local(IH) + DFP(I,IH) * VECT(I)
         ! kahan_y=DFP(I,IH)*VECT(I)-kahan_eps
         ! kahan_t=VOMA_local(IH)+kahan_y
         ! kahan_eps=(kahan_t-VOMA_local(IH))-kahan_y
         ! VOMA_local(IH)=kahan_t
#endif
325   CONTINUE
320   CONTINUE
#ifdef MPI
      call mpi_allreduce(VOMA_local, VOMA(:,1), NUH, mpi_real8, mpi_sum, parallel%commx, mpinfo)
#endif

   END SUBROUTINE OM1C!}}}
!-------------------------PARTING LINE--------------------------
   SUBROUTINE AMST(BETA,W0,NAM,NUH,DXP,DFP,SP,XL,FL,VOMA,XN)!{{{
        USE math , ONLY : SOPO
#ifdef MPI
      USE smpi_math_module , ONLY : parallel,mpi_real8,mpi_sum,mpinfo
#endif
        IMPLICIT NONE

        REAL(8), INTENT(IN) :: & 
           BETA, W0

        INTEGER, INTENT(IN) ::& 
           NAM, NUH

        REAL(8), DIMENSION(:,:), INTENT(IN) :: &
           DXP, DFP 

        REAL(8), DIMENSION(:), INTENT(IN) :: &
            XL, FL !,SP
        REAL(8),INTENT(IN) :: SP

        REAL(8), DIMENSION(:,:), INTENT(IN) :: &
           VOMA

        REAL(8), DIMENSION(:), INTENT(OUT) :: &
           XN

!C*************************************************
!C     ONE STEP OF THE ANDERSON MIXING PROCEDURE
!C     TO SOLVE NON-LINEAR EQUATIONS  F(X)=0
!C*************************************************
!C    BETA - MIXING PARAMETER
!C      W0 - A SMALL QUANTITY TO IMPROVE STABILITY
!C     NAM - NUMBER OF VARIABLES
!C     NUH - NUMBER OF PREVIOUS VECTORS (NUH.GE.2)
!C-------------------------------------------------
!C   ON INPUT:
!C     DXP(I,IH), I=1,...,NAM, IH=1,...,NUH,  -
!C              - DIFFERENCES OF PREVIOUS VECTORS X
!C     DFP(I,IH), I=1,...,NAM, IH=1,...,NUH,  -
!C              - DIFFERENCES OF PREVIOUS VECTORS F
!C     SP(I), I=1,...,NAM,  -  WEIGHTS FOR THE
!C                             SCALAR PRODUCT
!C     XL(I), I=1,...,NAM,  -  THE LAST VECTOR X
!C     FL(I), I=1,...,NAM,  -  THE LAST VECTOR F
!C     VOMA(IH,JH),  IH=1,...,NUH, JH=1,...,IH, -
!C              - LOWER TRIANGLE OF OVERLAP MATRIX
!C
!C   ON OUTPUT:
!C     XN(I), I=1,...,NAM,   -  THE NEW VECTOR X
!C*************************************************

      REAL(8), DIMENSION(NAM) :: & 
         VECT

      REAL(8), DIMENSION(NUH) :: & 
         WORK

      REAL(8), DIMENSION(NUH, 1) :: &
           T
#ifdef MPI
      REAL(8), DIMENSION(NUH) :: T_local
      ! !> kahan sum
      ! REAL(DP) ::kahan_y,kahan_t,kahan_eps
#endif

      REAL(8), DIMENSION(NUH, NUH) :: &
         A

       INTEGER :: & 
          IH, I, JH, MNUH

       REAL(8) :: & 
          DUM, BT

      MNUH=NUH
! #ifdef MPI
!       kahan_eps=0.d0
! #endif

      DO 301 IH=1,NUH
        T(IH,1)=0.d0
#ifdef MPI
        T_local(IH)=0.d0
#endif
301   CONTINUE
      DO 302 I=1,NAM
        XN(I)=0.d0
302   CONTINUE
!C                                    OVERLAP MATRIX
      DO 310 JH=1,NUH
      DO 311 IH=JH,NUH
        A(IH,JH)=VOMA(IH,JH)
311   CONTINUE
310   CONTINUE
      DUM=1.d0+W0**2
      DO 314 IH=1,NUH
        A(IH,IH)=DUM*A(IH,IH)
314   CONTINUE
!C                                        R.H.S.
      DO 320 I=1,NAM
        !VECT(I) = SP(I) * FL(I)
        VECT(I) = SP * FL(I)
320   CONTINUE


      DO 322 IH=1,NUH
      DO 324 I=1,NAM
#ifndef MPI
         T(IH,1) = T(IH,1) + DFP(I,IH) * VECT(I)
#else
         T_local(IH) = T_local(IH) + DFP(I,IH) * VECT(I)
         ! kahan_y=DFP(I,IH)*VECT(I)-kahan_eps
         ! kahan_t=T_local(IH)+kahan_y
         ! kahan_eps=(kahan_t-T_local(IH))-kahan_y
         ! T_local(IH)=kahan_t
#endif
324   CONTINUE
322   CONTINUE
#ifdef MPI
         CALL mpi_allreduce(T_local, T(:,1), NUH, mpi_real8, mpi_sum, parallel%commx, mpinfo)
#endif
!C                            SOLUTION OF LINEAR SYSTEM
      CALL SOPO(A,MNUH,NUH,T,MNUH,1,WORK)
!C                                    NEW VECTOR X
      DO 340 IH=1,NUH
        BT=BETA*T(IH,1)

      DO 341 I=1,NAM
        XN(I) = XN(I) + DFP(I,IH)*BT
341   CONTINUE

340   CONTINUE
      DO 342 IH=1,NUH
      DO 343 I=1,NAM
        XN(I) = XN(I) + DXP(I,IH)*T(IH,1)
343   CONTINUE
342   CONTINUE
      DO 345 I=1,NAM
        XN(I) = - XN(I) + BETA*FL(I)
345   CONTINUE
      DO 346 I=1,NAM
        XN(I) = XN(I) + XL(I)
346   CONTINUE

   END SUBROUTINE AMST!}}}
!-------------------------PARTING LINE--------------------------
!########################################################!
!*For    :    Simple Mixing + (r)Pulay                   !
!########################################################!
!--------------------Pulay Mixing---------------------------
   SUBROUTINE rPulay_mixing(Iter,XOUT,XIN,err)!{{{
      USE parameters  , ONLY : Nspin,NHMIX
      USE grid_module , ONLY : nsn,dvol
      IMPLICIT NONE
      !INOUT
      INTEGER(I4B),INTENT(IN) :: Iter !number of persent iteration
      REAL(DP),INTENT(INOUT)  :: XOUT(nsn),XIN(nsn)
                           !IN  : persent XIN and XOUT
                           !OUT : XIN/XOUT = next XIN
      REAL(DP),INTENT(OUT) :: err
      !LOCAL
      REAL(DP),DIMENSION(nsn) :: FL,XL,XN
      INTEGER(I4B) :: dime,nuhmax,nuhmin,niter &
                &,    Is ,Ig,I , IH , NUH1 ,DH , RITER
      REAL(DP) :: alpha,beta
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      dime=mixer%NAM
      nuhmax=mixer%NHMIX
      nuhmin=mixer%NHMIN
      niter=mixer%NITRA
      alpha=mixer%alpha
      beta=mixer%beta
      !----------------------
      !store RL
      XL(:)=XIN(:)
      FL(:)=XOUT(:)-XIN(:)
      !----------------------
      err=SQRT(DOT_PRODUCT(FL,FL))
      !----------------------
      IF(Iter>1) THEN
          mixer%DXL(:,1)=mixer%DXL(:,1)-XL(:)
          mixer%DFL(:,1)=mixer%DFL(:,1)-FL(:)
      END IF
      !----------------------
      IF(Iter<=MAX(niter,2))THEN
         !simple mixing
         XN(:)=XL(:)+alpha*FL(:)
      ELSE
          DH=NUHMAX-NUHMIN+1
          IF(NHMIX>0.OR.(Iter-1)<=NUHMIN)THEN
             !Anderson mixing
             NUH1=MIN(NUHMAX,Iter-1)
          ELSE
             !restart Anderson mixing
             RITER=Iter-NUHMIN-1
             NUH1=MOD(RITER,DH)+NUHMIN
          ENDIF

         !(r)Pulay mixing
         CALL rPulay_mix(beta,mixer%w0,dime,NUH1,&
                   & mixer%DXL,mixer%DFL,XL,FL,XN)
      ENDIF
      !-----------------------
      IF(Iter>0)THEN


         !ONE STEP BACKWARD
         DO IH=NUHMAX,2,-1
            mixer%DXL(:,IH)=mixer%DXL(:,IH-1)
            mixer%DFL(:,IH)=mixer%DFL(:,IH-1)
         END DO
         mixer%DXL(:,1)=XL
         mixer%DFL(:,1)=FL

      ENDIF
      !------------------------
      !output the densiy
      XOUT(:)=XN(:)
      XIN(:)=XN(:)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE rpulay_mixing!}}}
!-------------------------PARTING LINE--------------------------
   SUBROUTINE rpulay_mix(beta,w0,dime,NH,DXL,DRL,XL,RL,XN)!{{{
      USE Lapack_module , ONLY : invmat_real
      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B),INTENT(IN) :: dime,NH
      REAL(DP),INTENT(IN)  :: beta & !coe 
               &, w0  &   !for ill matrix
               &, DXL(:,:),DRL(:,:) &
               &, XL(:),RL(:)
      REAL(DP),INTENT(OUT) :: XN(:)
      !LOCAL
      INTEGER(I4B) :: I
      REAL(DP) :: invAM(NH,NH),DRR(NH),cl(NH) !,w2
      REAL(DP) :: Xbar(dime),Rbar(dime)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF(NH<=0)THEN
         WRITE(6,*) 'rpulay_mix: NH must > 0',NH
      ENDIF
      !AM
      invAM(:,:)=MATMUL(TRANSPOSE(DRL(:,1:NH)),DRL(:,1:NH))
      !for ill
      !w2=w0**2
      !DO I=1,NH
      !   invAM(I,I)=invAM(I,I) + w2
      !ENDDO
      CALL invmat_real(invAM)
      !solve the value
      DO I=1,NH
         DRR(I)=DOT_PRODUCT(DRL(:,I),RL(:))
      ENDDO
      cl(:)= -1.d0*MATMUL(invAM,DRR) 
      !XN(:)=Xbar(:)+beta*Rbar(:)
      Xbar(:)=XL(:)
      Rbar(:)=RL(:)
      DO I=1,NH
         Xbar(:)=Xbar(:)+cl(I)*DXL(:,I)
         Rbar(:)=Rbar(:)+cl(I)*DRL(:,I)
      ENDDO
      XN(:)=Xbar(:)+beta*Rbar(:)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE rpulay_mix!}}}
!-------------------real-space-kerker---------------------------
!   SUBROUTINE real_kerker(XL,RL,alpha,XN)!{{{
!      USE parameters , ONLY : Nspin,finite_order
!      USE grid_module , ONLY : n1,n2,n3,nsn,dvol
!      USE struct_module , ONLY : ncharge
!      USE finite_module , ONLY : real_nabla2 
!      IMPLICIT NONE
!      REAL(DP),INTENT(IN)  :: XL(n1,n2,n3,Nspin) &
!                         &,   RL(n1,n2,n3,Nspin) &
!                         &,   alpha
!      REAL(DP),INTENT(OUT) :: XN(n1,n2,n3,Nspin)
!      !
!      REAL(DP),DIMENSION(n1,n2,n3,Nspin) :: &
!                             &LR,LLR,LLLR,GR
!      INTEGER(I4B) :: Is,ix,iy,iz
!      REAL(DP) :: a1,a2,a3,t1,t2,t0
!      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!      a1 = alpha/mixer%BMIX**2
!      a2 = alpha/mixer%BMIX**4
!      a3 = alpha/mixer%BMIX**6
!      GR=0.d0
!      !1st-order
!      DO Is=1,Nspin
!         CALL real_nabla2(RL(:,:,:,Is),finite_order,LR(:,:,:,Is))
!         !update
!         DO iz=1,n3
!         DO iy=1,n2
!         DO ix=1,n1
!            GR(ix,iy,iz,Is) =      GR(ix,iy,iz,Is)   &
!                  &         + a1 * LR(ix,iy,iz,Is)
!         ENDDO
!         ENDDO
!         ENDDO
!      ENDDO
!      !2nd-0rder
!      DO Is=1,Nspin
!         CALL realnabla2(LR(:,:,:,Is),finite_order,LLR(:,:,:,Is))
!         !update
!         DO iz=1,n3
!         DO iy=1,n2
!         DO ix=1,n1
!            GR(ix,iy,iz,Is) =      GR(ix,iy,iz,Is)   &
!                  &         + a2 * LLR(ix,iy,iz,Is)
!         ENDDO
!         ENDDO
!         ENDDO
!      ENDDO
!      !3rd-order   
!      DO Is=1,Nspin
!         CALL realnabla2(LLR(:,:,:,Is),finite_order,LLLR(:,:,:,Is))
!         !update
!         DO iz=1,n3
!         DO iy=1,n2
!         DO ix=1,n1
!            GR(ix,iy,iz,Is) =      GR(ix,iy,iz,Is)   &
!                  &         + a3 * LLLR(ix,iy,iz,Is)
!         ENDDO
!         ENDDO
!         ENDDO
!      ENDDO
!      !--------------------------------------------------
!      !Normalize
!      GR=GR-SUM(GR)/nsn
!      !Output XN
!      DO Is=1,Nspin
!         !update
!         DO iz=1,n3
!         DO iy=1,n2
!         DO ix=1,n1
!            t1=alpha*(RL(ix,iy,iz,Is))
!            t2=GR(ix,iy,iz,Is)
!            IF(ABS(t1)>ABS(t2))THEN
!                t0=t2
!            ELSE
!                t0=t1
!            ENDIF
!            XN(ix,iy,iz,Is) = XL(ix,iy,iz,Is) + t0 
!         ENDDO
!         ENDDO
!         ENDDO
!      ENDDO
!      !normalize
!      XN(:,:,:,:)=XN/(SUM(XN)*dvol)*ncharge
!      !--------------------------------------------------
!      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!   ENDSUBROUTINE real_kerker!}}}
!########################################################!
!*For    : Simple Mixing + (r)Pulay + kerker             !
!########################################################!
!--------------------------init kerker--------------------------
   SUBROUTINE init_kerker()!{{{
#ifdef MPI
      USE smpi_math_module, ONLY:parallel
      USE grid_module , ONLY : ng1,ng2,grid
#else
      USE grid_module , ONLY : ng,grid
#endif
      IMPLICIT NONE
      INTEGER(I4B) :: Ig
      REAL(DP) :: kerk,g2,BMIX2
#ifdef MPI
      INTEGER(I4B) :: g_shift
#endif
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !cycle all recip-space
      BMIX2=mixer%BMIX**2
#ifdef MPI
      g_shift=parallel%local_z_start*ng2*ng1
      DO Ig=1,parallel%local_z*ng2*ng1
         g2=grid%gVec(4,Ig+g_shift)**2
         kerk=g2/(g2+BMIX2)
         mixer%kerker(Ig)=MAX(kerk,mixer%AMIX)
      ENDDO
#else
      DO Ig=1,ng
         g2=grid%gVec(4,Ig)**2
         kerk=g2/(g2+BMIX2)
         !kerk=g2/(g2+BMIX2)
         mixer%kerker(Ig)=MAX(kerk,mixer%AMIX)
         !mixer%kerker(Ig)=kerk
      ENDDO
#endif
      !mixer%kerker(1)=0.d0
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE init_kerker!}}}
!--------------------------init kerker--------------------------
!########################################################!
!* mixer%NUH > 0 : Pulay mixing + kerker                 !
!* mixer%NUH < 0 : rPulay mixing                         ! 
!########################################################!
   SUBROUTINE rPulayK_mixing(Iter,RLg,XINg)!{{{
      USE parameters  , ONLY : Nspin,NHMIX
#ifdef MPI
      USE smpi_math_module, ONLY:parallel
      USE grid_module , ONLY : ng1,ng2
#else
      USE grid_module , ONLY : ng
#endif
      IMPLICIT NONE
      !INOUT
      INTEGER(I4B),INTENT(IN) :: Iter !number of persent iteration
      COMPLEX(DP),INTENT(INOUT)  :: RLg(mixer%NAM),XINg(mixer%NAM)
                           !IN  : persent XIN and RLg
                           !OUT : XIN/XOUT = next XIN
      !LOCAL
      COMPLEX(DP),DIMENSION(mixer%NAM) :: RL,XL,XN
      INTEGER(I4B) :: dime,nuhmax,nuhmin,niter &
                &,    Is ,Ig,I , IH , NUH1 , RITER , DH
      REAL(DP) :: alpha,beta
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      dime=mixer%NAM
      nuhmax=mixer%NHMIX
      nuhmin=mixer%NHMIN
      niter=mixer%NITRA
      alpha=mixer%alpha
      beta=mixer%beta
      !----------------------
      !store RL
      XL(:)=XINg(:)
      RL(:)=RLg(:)
      !----------------------
      IF(Iter>1) THEN
          mixer%DXGL(:,1)=mixer%DXGL(:,1)-XL(:)
          mixer%DFGL(:,1)=mixer%DFGL(:,1)-RL(:)
      END IF
      !----------------------
      IF(Iter<=MAX(niter,2))THEN
         !Kerker simple mixing
         I=0
         DO Is=1,Nspin
#ifdef MPI
            DO Ig=1,parallel%local_z*ng2*ng1
               I=I+1
               XN(I)=XL(I)+alpha*mixer%kerker(Ig)*RL(I)
               !XN(I)=XL(I)+alpha*RL(I)
            ENDDO
#else
            DO Ig=1,ng
               I=I+1
               XN(I)=XL(I)+alpha*mixer%kerker(Ig)*RL(I)
               !XN(I)=XL(I)+alpha*RL(I)
            ENDDO
#endif
         ENDDO
      ELSE
         NUH1=MIN(NUHMAX,Iter-1)
         !(r)Pulay mixing
         CALL rPulayK_mix(beta,mixer%w0,dime,NUH1,&
                   & mixer%DXGL,mixer%DFGL,XL,RL,XN)

      ENDIF
      !-----------------------
      IF(Iter>0)THEN


         !ONE STEP BACKWARD
         DO IH=NUHMAX,2,-1
            mixer%DXGL(:,IH)=mixer%DXGL(:,IH-1)
            mixer%DFGL(:,IH)=mixer%DFGL(:,IH-1)
         END DO
         mixer%DXGL(:,1)=XL
         mixer%DFGL(:,1)=RL
      
      ENDIF
      !------------------------
      !output the densiy
      RLg(:)=XN(:)
      XINg(:)=XN(:)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE rpulayK_mixing!}}}
!----------------------------rPulay_mix---------------------------
   SUBROUTINE rpulayK_mix(beta,w0,dime,NH,DXL,DRL,XL,RL,XN)!{{{
      USE parameters , ONLY : Nspin
#ifdef MPI
       USE smpi_math_module , ONLY : parallel,mpi_complex16&
            &,mpi_sum,mpinfo
      USE grid_module , ONLY : ng1,ng2
#else
      USE grid_module , ONLY : ng
#endif
      USE Lapack_module , ONLY : matmat,invmat
      IMPLICIT NONE
      !IN/OUT
      REAL(DP),INTENT(IN) :: beta,w0
      INTEGER(I4B),INTENT(IN) :: dime , NH
      COMPLEX(DP),INTENT(IN)  :: DXL(:,:),DRL(:,:)  &
                           &,     XL(:),RL(:)
      COMPLEX(DP),INTENT(OUT) :: XN(:)
      !LOCAL
      COMPLEX(DP) :: invAM(NH,NH),DRR(NH),cl(NH) &
                   &,Xbar(dime),Rbar(dime)
#ifdef MPI
      COMPLEX(DP) :: invAM_local(NH,NH),DRR_local
#endif
      INTEGER(I4B) :: I,Is,Ig,I_end
      !REAL(DP) :: w2
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF(NH<=0)THEN
         WRITE(6,*) 'rpulay_mix: NH must > 0',NH
         STOP
      ENDIF
      !AM
#ifdef MPI
      if(size(DRL)==0)then
         invAM_local=cmplx(0.d0,0.d0)
      else
      CALL matmat(DRL(:,1:NH),DRL(:,1:NH) &
              & ,'C','N',invAM_local)
      endif
      CALL MPI_ALLREDUCE(invAM_local,invAM,NH**2,MPI_COMPLEX16,MPI_SUM,parallel%commx,mpinfo)
#else
      CALL matmat(DRL(:,1:NH),DRL(:,1:NH) &
              & ,'C','N',invAM)
#endif
      !CALL invAM
      CALL invmat(invAM)
      !solve
      DO I=1,NH
#ifdef MPI
         DRR_local=DOT_PRODUCT(DRL(:,I),RL(:))
         CALL MPI_ALLREDUCE(DRR_local,DRR(I),1,MPI_COMPLEX16,MPI_SUM,parallel%commx,mpinfo)
#else
         DRR(I)=DOT_PRODUCT(DRL(:,I),RL(:))
#endif
      ENDDO
      cl(:)=-1.d0*MATMUL(invAM,DRR)
      !call vector
      Xbar(:)=XL(:)
      Rbar(:)=RL(:)
      DO I=1,NH
         Xbar(:)=Xbar(:)+cl(I)*DXL(:,I)
         Rbar(:)=Rbar(:)+cl(I)*DRL(:,I)
      ENDDO
      !output
      I=0
#ifdef MPI
      I_end=parallel%local_z*ng2*ng1
#else
      I_end=ng
#endif
      DO Is=1,Nspin
         DO Ig=1,I_end
            I=I+1
            XN(I)=Xbar(I)+beta*mixer%kerker(Ig)*Rbar(I)
         ENDDO
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE rpulayK_mix!}}}
!-------------------------PARTING LINE--------------------------
!* mixer%NUH > 0 : Pulay mixing + resta                 !
!* mixer%NUH < 0 : rPulay mixing                         ! 
   SUBROUTINE resta_mixing(Iter,RLg,XINg)
     USE parameters  , ONLY : Nspin,NHMIX
#ifdef MPI
     USE smpi_math_module, ONLY:parallel
     USE grid_module , ONLY : ng1,ng2
#else
     USE grid_module , ONLY : ng
#endif
     IMPLICIT NONE
     !INOUT
     INTEGER(I4B),INTENT(IN) :: Iter !number of persent iteration
     COMPLEX(DP),INTENT(INOUT)  :: RLg(mixer%NAM),XINg(mixer%NAM)
     !IN  : persent XIN and RLg
     !OUT : XIN/XOUT = next XIN
     !LOCAL
     COMPLEX(DP),DIMENSION(mixer%NAM) :: RL,XL,XN
     INTEGER(I4B) :: dime,nuhmax,nuhmin,niter &
          &,    Is ,Ig,I , IH , NUH1 , RITER , DH
     REAL(DP) :: alpha ,beta
     INTEGER(I4B) :: I_end
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      dime=mixer%NAM
      nuhmax=mixer%NHMIX
      nuhmin=mixer%NHMIN
      niter=mixer%NITRA
      alpha=mixer%alpha
      beta=mixer%beta
     !----------------------
     !store RL
     XL(:)=XINg(:)
     RL(:)=RLg(:)
     !----------------------
     IF(Iter>1) THEN
        mixer%DXGL(:,1)=mixer%DXGL(:,1)-XL(:)
        mixer%DFGL(:,1)=mixer%DFGL(:,1)-RL(:)
     END IF
     !----------------------
     IF(Iter<=MAX(niter,2))THEN
        !Kerker simple mixing
        I=0
        DO Is=1,Nspin
#ifdef MPI
           I_end=parallel%local_z*ng2*ng1
#else
           I_end=ng
#endif
           DO Ig=1,I_end
              I=I+1
              XN(I)=XL(I)+alpha*mixer%kerker(Ig)*RL(I)
              !XN(I)=XL(I)+alpha*RL(I)
           ENDDO
        ENDDO
     ELSE
        NUH1=MIN(NUHMAX,Iter-1)
        !(r)Pulay mixing
        CALL rPulayK_mix(beta,mixer%w0,dime,NUH1,&
             & mixer%DXGL,mixer%DFGL,XL,RL,XN)


     ENDIF
     !-----------------------
     IF(Iter>0)THEN


        !ONE STEP BACKWARD
        DO IH=NUHMAX,2,-1
           mixer%DXGL(:,IH)=mixer%DXGL(:,IH-1)
           mixer%DFGL(:,IH)=mixer%DFGL(:,IH-1)
        END DO
        mixer%DXGL(:,1)=XL
        mixer%DFGL(:,1)=RL

     ENDIF
     !------------------------
     !output the densiy
     RLg(:)=XN(:)
     XINg(:)=XN(:)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE resta_mixing
!-------------------------PARTING LINE--------------------------
   SUBROUTINE init_resta()!{{{
     USE parameters, ONLY:resta
#ifdef MPI
     USE smpi_math_module, ONLY:parallel
     USE grid_module , ONLY : ng1,ng2,grid
#else
     USE grid_module , ONLY : ng,grid
#endif
     IMPLICIT NONE
     INTEGER(I4B) :: Ig
     REAL(DP) :: kerk,g2,g_abs
#ifdef MPI
     INTEGER(I4B) :: g_shift
     INTEGER(I4B) :: Ig_start
#endif
     REAL(DP) :: q0,epsilon0,Rs
     REAL(DP) :: a1,q02
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     !cycle all recip-space
     ! q0=0.68
     ! epsilon0=10
     ! Rs=6.61
     q0=resta(1)
     epsilon0=resta(2)
     Rs=resta(3)

     q02=q0**2
#ifdef MPI
     g_shift=parallel%local_z_start*ng2*ng1
     if(g_shift==0)then
        mixer%kerker(1)=1/epsilon0
        Ig_start=2
     else
        Ig_start=1
     endif
     DO Ig=Ig_start,parallel%local_z*ng2*ng1
        g2=grid%gVec(4,Ig+g_shift)**2
        g_abs=sqrt(g2)
        a1=q02*sin(g_abs*Rs)/(epsilon0*g_abs*Rs)
        kerk=(a1+g2)/(g2+q02)
        mixer%kerker(Ig)=kerk
     ENDDO
#else
     mixer%kerker(1)=1/epsilon0
     DO Ig=2,ng
        g2=grid%gVec(4,Ig+g_shift)**2
        g_abs=sqrt(g2)
        a1=q02*sin(g_abs*Rs)/(epsilon0*g_abs*Rs)
        kerk=(a1+g2)/(g2+q02)
        mixer%kerker(Ig)=kerk
     ENDDO
#endif
     !mixer%kerker(1)=0.d0
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE init_resta!}}}
!-------------------------PARTING LINE--------------------------
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDMODULE mixer_module

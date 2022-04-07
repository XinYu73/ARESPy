MODULE mixer_module
!###################################################!
!*For   : Mixing the charge density to convergence  !
!*Author: Qiang Xu                                  !
!*Data  : 2017/07/30                                !
!*IMIXER:   0   For simple mixing + Anderson mixing !
!         other For simple mixing + (r)Pulay mixing !
!###################################################!
   USE constants
   USE parameters , ONLY : IMIXER,NMITER,NHMAX,NHMIN &
         &,MALPHA,MBETA,NSMIX,W0AM,nspin,AMIX,BMIX
   USE grid_module , ONLY : dvol,ng
   IMPLICIT NONE
   TYPE mixer_data
        REAL(DP),DIMENSION(:,:),ALLOCATABLE :: &
          DXL, &      !delta vector
          DFL, &      !delta rediuse vector
          VOMA
        !For Kerker mixing (recip-space)
        REAL(DP),ALLOCATABLE :: kerker(:)
        COMPLEX(DCP),ALLOCATABLE ::  &
          &    DXGL(:,:) & !input vector
          &,   DRGL(:,:)   !residue
   ENDTYPE mixer_data
   !data for mixer
   INTEGER(I4B) :: NAM,NUH 
   TYPE(mixer_data) :: mixer
CONTAINS
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!-------------------------------------------------------------------------
   SUBROUTINE init_mixer(nps)
#ifdef MPI
      USE grid_module, ONLY:ng1,ng2
      USE smpi_math_module, ONLY: parallel
#endif
      IMPLICIT NONE
      !
      INTEGER(I4B),INTENT(IN) :: nps
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      NAM=nps*nspin
      NUH=NHMAX
      CALL destroy_mixer()
      SELECT CASE(IMIXER)
      CASE(0)
         !real space mixing
         ALLOCATE(mixer%DXL(NAM,NUH))
         ALLOCATE(mixer%DFL(NAM,NUH))
         ALLOCATE(mixer%VOMA(NUH,NUH))
         mixer%DXL = 0.d0
         mixer%DFL = 0.d0
         mixer%VOMA = 0.d0
      CASE(2)
#ifdef MPI
         NAM=nspin*ng1*ng2*parallel%local_z
         !recip-space mixing
         ALLOCATE(mixer%kerker(ng1*ng2*parallel%local_z))
#else
         NAM=nspin*ng
         !recip-space mixing
         ALLOCATE(mixer%kerker(ng))
#endif
         ALLOCATE(mixer%DXGL(NAM,NUH))
         ALLOCATE(mixer%DRGL(NAM,NUH))
         mixer%DXGL(:,:)=0.d0
         mixer%DRGL(:,:)=0.d0
         CALL init_resta()
      CASE default
#ifdef MPI
         NAM=nspin*ng1*ng2*parallel%local_z
         !recip-space mixing
         ALLOCATE(mixer%kerker(ng1*ng2*parallel%local_z))
#else
         NAM=nspin*ng
         !recip-space mixing
         ALLOCATE(mixer%kerker(ng))
#endif
         ALLOCATE(mixer%DXGL(NAM,NUH))
         ALLOCATE(mixer%DRGL(NAM,NUH))
         mixer%DXGL(:,:)=0.d0
         mixer%DRGL(:,:)=0.d0
         CALL init_kerker()
      ENDSELECT
      IF(MBETA<0._DP) MBETA=0.8*MALPHA
      !WEIGHTS FOR SCALAR PRODUCT
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE init_mixer
!-----------------------------PARTING LINE---------------------------------
   SUBROUTINE destroy_mixer()
      !
      IMPLICIT NONE
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !real space mixing
      IF(ALLOCATED(mixer%DXL))  DEALLOCATE(mixer%DXL)
      IF(ALLOCATED(mixer%DFL))  DEALLOCATE(mixer%DFL)
      IF(ALLOCATED(mixer%VOMA)) DEALLOCATE(mixer%VOMA)
      IF(ALLOCATED(mixer%kerker)) DEALLOCATE(mixer%kerker)
      IF(ALLOCATED(mixer%DXGL)) DEALLOCATE(mixer%DXGL)
      IF(ALLOCATED(mixer%DRGL)) DEALLOCATE(mixer%DRGL)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE destroy_mixer
!--------------------------------mixing------------------------------------
   SUBROUTINE mixing(Iter,XOUT,XIN,Res)
     USE parameters , ONLY : nspin
#ifdef MPI
     USE smpi_math_module , ONLY : parallel,mpi_real8,mpi_sum,mpinfo
     USE grid_module , ONLY : dvol,n,ng1,ng2,fft_sph
#else
     USE grid_module , ONLY : dvol,n,ng
#endif
     USE struct_module , ONLY : ncharge
     IMPLICIT NONE
     !
     INTEGER,INTENT(IN) :: Iter !present iteration step
     REAL(DP),INTENT(INOUT) :: XIN(n,nspin),XOUT(n,nspin)
     REAL(DP),INTENT(OUT) :: Res !Norm of residual of density
     COMPLEX(DCP),ALLOCATABLE :: XINg(:,:,:,:),RLg(:,:,:,:)
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
        !real space
        CALL Anderson_mixing(Iter,XIN,XOUT,Res)


     CASE(2) !resta
        !recip. space
        !allocate
#ifdef MPI
        ALLOCATE(XINg(ng1,ng2,parallel%local_z,Nspin)&
             &,RLg(ng1,ng2,parallel%local_z,Nspin))
#else
        ALLOCATE(XINg(ng1,ng2,ng3,Nspin),RLg(ng1,ng2,ng3,Nspin))
#endif
        !store the res
        XOUT=XOUT-XIN
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
           XINg(:,:,:,Is)  = FFT_sph(XIN(:,Is))
           RLg(:,:,:,Is)   = FFT_sph(XOUT(:,Is))
#else
           !shape err     XINg(:,:,:,Is)  = FFT(XIN(:,:,:,Is))
           !              RLg(:,:,:,Is)   = FFT(XOUT(:,:,:,Is))
#endif
        ENDDO
        CALL resta_mixing(Iter,RLg,XINg)
        !to real-space
        DO Is=1,Nspin
#ifdef MPI
           XIN(:,Is)  = FFT_sph(XINg(:,:,:,Is))
#else
           !shape err     XIN(:,:,:,Is)  = FFT(XINg(:,:,:,Is))
#endif
        ENDDO
        XOUT=XIN
        DEALLOCATE(XINg,RLg)


     CASE default
        !recip. space
        !allocate
#ifdef MPI
        ALLOCATE(XINg(ng1,ng2,parallel%local_z,Nspin)&
             &,RLg(ng1,ng2,parallel%local_z,Nspin))
#else
        ALLOCATE(XINg(ng1,ng2,ng3,Nspin),RLg(ng1,ng2,ng3,Nspin))
#endif
        !store the res
        XOUT=XOUT-XIN
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
           XINg(:,:,:,Is)  = FFT_sph(XIN(:,Is))
           RLg(:,:,:,Is)   = FFT_sph(XOUT(:,Is))
#else
           !shape err     XINg(:,:,:,Is)  = FFT(XIN(:,:,:,Is))
           !              RLg(:,:,:,Is)   = FFT(XOUT(:,:,:,Is))
#endif
        ENDDO
        CALL rPulayK_mixing(Iter,RLg,XINg)
        !to real-space
        DO Is=1,Nspin
#ifdef MPI
           XIN(:,Is)  = FFT_sph(XINg(:,:,:,Is))
#else
           !shape err     XIN(:,:,:,Is)  = FFT(XINg(:,:,:,Is))
#endif
        ENDDO
        XOUT=XIN
        !scale
        DEALLOCATE(XINg,RLg)
     ENDSELECT
     !res=res
     res=res*dvol
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE mixing
!########################################################!
!*For    : Simple Mixing + Anderson mixing               !
!########################################################!
!-----------------------------PARTING LINE----------------------------------
   SUBROUTINE Anderson_mixing(IITER,XIN,XOUT,err)
      USE parameters , ONLY : NHMAX
#ifdef MPI
      USE smpi_math_module , ONLY : parallel,mpi_real8,mpi_sum,mpinfo
#endif
      !
      IMPLICIT NONE
      !
      INTEGER,INTENT(IN) :: IITER !present iteration step
      REAL(DP),INTENT(INOUT) :: XIN(NAM),XOUT(NAM) !
      REAL(DP),INTENT(OUT) :: err
      !
      REAL(DP) :: alpha, beta !factor for mixing
      REAL(DP), DIMENSION(NAM) :: XL, FL, XN
      integer :: &
        NUHMAX,NUHMIN,NSIMPLE, &
        J, I1, I2, I3, &
        IH, JH, &
        NUH1,RITER,DH
#ifdef MPI
      REAL(DP) :: err_loc
#endif
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      NUHMAX = NHMAX
      NUHMIN = NHMIN
      NSIMPLE = NSMIX
      alpha = MALPHA
      beta  = MBETA
      !>>>input rho
      XL(:)=XIN
      FL(:)=XOUT(:)-XIN(:)
      !<<<input rho 
#ifdef MPI
      err_loc=DOT_PRODUCT(FL,FL)
      CALL MPI_ALLREDUCE(err_loc,err,1,MPI_REAL8,MPI_SUM,parallel%commx,mpinfo)
#else
      err=DOT_PRODUCT(FL,FL)
#endif
      err=sqrt(err)

      IF(IITER>1) THEN
          mixer%DXL(:,1)=mixer%DXL(:,1)-XL(:)
          mixer%DFL(:,1)=mixer%DFL(:,1)-FL(:)
          CALL OM1C(NAM, NUHMAX, dvol, mixer%DFL, mixer%VOMA)
      END IF

      IF(IITER<=NSIMPLE.OR.IITER<=2) THEN
          !simple mixing
          XN=XL+alpha*FL
      ELSE
          !Anderson mixing
          NUH1=MIN(NUHMAX,IITER-1)
          !print*,'Iter,NUH1',IITER,NUH1,NUHMAX
          
          !call Anderson
          !IF(NUH1==0)THEN
          !ELSE
          CALL AMST(beta,W0AM,NAM,NUH1,&
                    mixer%DXL,mixer%DFL,dvol,&
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
   ENDSUBROUTINE Anderson_mixing
!-----------------------PARTING LINE---------------------------
   SUBROUTINE OM1C(NAM,NUH,SP,DFP,VOMA)
#ifdef MPI
       USE smpi_math_module , ONLY : parallel,mpi_real8,mpi_sum,mpinfo
#endif

       IMPLICIT NONE

       INTEGER, INTENT(IN) :: & 
          NAM, NUH

       !REAL*8, DIMENSION(:), INTENT(IN) :: & 
       !   SP
       REAL*8,INTENT(IN) :: SP
 
       REAL*8, DIMENSION(:,:), INTENT(IN) :: & 
          DFP
  
       REAL*8, DIMENSION(:,:), INTENT(OUT) :: &
          VOMA
#ifdef MPI
       REAL*8,DIMENSION(NUH) :: VOMA_loc
#endif
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

       REAL*8, DIMENSION(NAM) :: & 
          VECT

       INTEGER :: & 
          IH, I

      DO 310 IH=1, NUH
        VOMA(IH,1)=0.d0
#ifdef MPI
        VOMA_loc(IH)=0._DP
#endif
310   CONTINUE

      DO 312 I=1,NAM
        !VECT(I) = SP(I) * DFP(I,1)
        VECT(I) = SP * DFP(I,1)
312   CONTINUE

      DO 320 IH=1,NUH
      DO 325 I=1,NAM
#ifdef MPI
        VOMA_loc(IH) = VOMA_loc(IH) + DFP(I,IH) * VECT(I)
#else
        VOMA(IH,1) = VOMA(IH,1) + DFP(I,IH) * VECT(I)
#endif
325   CONTINUE
320   CONTINUE
#ifdef MPI
      CALL MPI_ALLREDUCE(VOMA_loc,VOMA(:,1),NUH,MPI_REAL8,MPI_SUM,parallel%commx,mpinfo)
#endif
   END SUBROUTINE OM1C
!-------------------------PARTING LINE--------------------------
   SUBROUTINE AMST(BETA,W0,NAM,NUH,DXP,DFP,SP,XL,FL,VOMA,XN)
        USE math , ONLY : SOPO
#ifdef MPI
        USE smpi_math_module , ONLY : parallel,mpi_real8,mpi_sum,mpinfo
#endif
        IMPLICIT NONE

        REAL*8, INTENT(IN) :: & 
           BETA, W0

        INTEGER, INTENT(IN) ::& 
           NAM, NUH

        REAL*8, DIMENSION(:,:), INTENT(IN) :: &
           DXP, DFP 

        REAL*8, DIMENSION(:), INTENT(IN) :: &
            XL, FL !,SP
        REAL*8,INTENT(IN) :: SP

        REAL*8, DIMENSION(:,:), INTENT(IN) :: &
           VOMA

        REAL*8, DIMENSION(:), INTENT(OUT) :: &
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

      REAL*8, DIMENSION(NAM) :: & 
         VECT

      REAL*8, DIMENSION(NUH) :: & 
         WORK

      REAL*8, DIMENSION(NUH, 1) :: &
         T
#ifdef MPI
      REAL*8,DIMENSION(NUH) :: T_loc
#endif
        
      REAL*8, DIMENSION(NUH, NUH) :: &
         A

       INTEGER :: & 
          IH, I, JH, MNUH

       REAL*8 :: & 
          DUM, BT

      MNUH=NUH

      DO 301 IH=1,NUH
        T(IH,1)=0.d0
#ifdef MPI
        T_loc(IH)=0._DP
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
#ifdef MPI
        T_loc(IH) = T_loc(IH) + DFP(I,IH) * VECT(I)
#else
        T(IH,1) = T(IH,1) + DFP(I,IH) * VECT(I)
#endif
324   CONTINUE
322   CONTINUE
!C
#ifdef MPI
      CALL MPI_ALLREDUCE(T_loc,T(:,1),NUH,MPI_REAL8,MPI_SUM,parallel%commx,mpinfo)
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

   END SUBROUTINE AMST
!-------------------------PARTING LINE--------------------------
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
      BMIX2=BMIX**2
#ifdef MPI
      g_shift=parallel%local_z_start*ng2*ng1
      DO Ig=1,parallel%local_z*ng2*ng1
         g2=grid%gVec(4,Ig+g_shift)**2
         kerk=g2/(g2+BMIX2)
         mixer%kerker(Ig)=MAX(kerk,AMIX)
      ENDDO
#else
      DO Ig=1,ng
         g2=grid%gVec(4,Ig)**2
         kerk=g2/(g2+BMIX2)
         !kerk=g2/(g2+BMIX2)
         mixer%kerker(Ig)=MAX(kerk,AMIX)
         !mixer%kerker(Ig)=kerk
      ENDDO
#endif
      !mixer%kerker(1)=0.d0
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE init_kerker!}}}
!-------------------------PARTING LINE--------------------------
!* mixer%NUH > 0 : Pulay mixing + kerker                 !
!* mixer%NUH < 0 : rPulay mixing                         ! 
   SUBROUTINE rPulayK_mixing(Iter,RLg,XINg)
     USE parameters  , ONLY : Nspin,NHMAX
#ifdef MPI
     USE smpi_math_module, ONLY:parallel
     USE grid_module , ONLY : ng1,ng2
#else
     USE grid_module , ONLY : ng
#endif
      IMPLICIT NONE
      !INOUT
      INTEGER(I4B),INTENT(IN) :: Iter !number of persent iteration
      COMPLEX(DCP),INTENT(INOUT)  :: RLg(NAM),XINg(NAM)
                           !IN  : persent XIN and RLg
                           !OUT : XIN/XOUT = next XIN
      !LOCAL
      COMPLEX(DCP),DIMENSION(NAM) :: RL,XL,XN
      INTEGER(I4B) :: dime,nuhmax,nuhmin,niter &
                &,    Is ,Ig,I , IH , NUH1 , RITER , DH
      REAL(DP) :: alpha ,beta
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      dime=NAM
      nuhmax=NHMAX
      nuhmin=NHMIN
      niter=NSMIX
      alpha=malpha
      beta=mbeta
      !----------------------
      !store RL
      XL(:)=XINg(:)
      RL(:)=RLg(:)
      !----------------------
      IF(Iter>1) THEN
          mixer%DXGL(:,1)=mixer%DXGL(:,1)-XL(:)
          mixer%DRGL(:,1)=mixer%DRGL(:,1)-RL(:)
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
         CALL rPulayK_mix(beta,W0AM,dime,NUH1,&
                   & mixer%DXGL,mixer%DRGL,XL,RL,XN)


      ENDIF
      !-----------------------
      IF(Iter>0)THEN


         !ONE STEP BACKWARD
         DO IH=NUHMAX,2,-1
            mixer%DXGL(:,IH)=mixer%DXGL(:,IH-1)
            mixer%DRGL(:,IH)=mixer%DRGL(:,IH-1)
         END DO
         mixer%DXGL(:,1)=XL
         mixer%DRGL(:,1)=RL
      
      ENDIF
      !------------------------
      !output the densiy
      RLg(:)=XN(:)
      XINg(:)=XN(:)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE rpulayK_mixing
!-------------------------PARTING LINE--------------------------
   SUBROUTINE rpulayK_mix(beta,w0,dime,NH,DXL,DRL,XL,RL,XN)
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
      REAL(DP),INTENT(IN) :: w0,beta
      INTEGER(I4B),INTENT(IN) :: dime , NH
      COMPLEX(DCP),INTENT(IN)  :: DXL(:,:),DRL(:,:)  &
                           &,     XL(:),RL(:)
      COMPLEX(DCP),INTENT(OUT) :: XN(:)
      !LOCAL
      COMPLEX(DCP) :: invAM(NH,NH),DRR(NH),cl(NH) &
                   &,Xbar(dime),Rbar(dime)
#ifdef MPI
      COMPLEX(DCP) :: invAM_local(NH,NH),DRR_local
#endif
      INTEGER(I4B) :: I,Is,Ig,I_end
      !REAL(DP) :: w2
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF(NH<=0)THEN
         WRITE(*,*) 'rpulay_mix: NH must > 0',NH
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
            !XN(I)=Xbar(I)+mixer%kerker(Ig)*Rbar(I)
            !XN(I)=Xbar(I)+beta*mixer%kerker(Ig)*RL(I)
         ENDDO
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE rpulayK_mix
!-------------------------PARTING LINE--------------------------
!* mixer%NUH > 0 : Pulay mixing + resta                 !
!* mixer%NUH < 0 : rPulay mixing                         ! 
   SUBROUTINE resta_mixing(Iter,RLg,XINg)
     USE parameters  , ONLY : Nspin,NHMAX
#ifdef MPI
     USE smpi_math_module, ONLY:parallel
     USE grid_module , ONLY : ng1,ng2
#else
     USE grid_module , ONLY : ng
#endif
     IMPLICIT NONE
     !INOUT
     INTEGER(I4B),INTENT(IN) :: Iter !number of persent iteration
     COMPLEX(DCP),INTENT(INOUT)  :: RLg(NAM),XINg(NAM)
     !IN  : persent XIN and RLg
     !OUT : XIN/XOUT = next XIN
     !LOCAL
     COMPLEX(DCP),DIMENSION(NAM) :: RL,XL,XN
     INTEGER(I4B) :: dime,nuhmax,nuhmin,niter &
          &,    Is ,Ig,I , IH , NUH1 , RITER , DH
     REAL(DP) :: alpha ,beta
     INTEGER(I4B) :: I_end
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     dime=NAM
     nuhmax=NHMAX
     nuhmin=NHMIN
     niter=NSMIX
     alpha=malpha
     beta=mbeta
     !----------------------
     !store RL
     XL(:)=XINg(:)
     RL(:)=RLg(:)
     !----------------------
     IF(Iter>1) THEN
        mixer%DXGL(:,1)=mixer%DXGL(:,1)-XL(:)
        mixer%DRGL(:,1)=mixer%DRGL(:,1)-RL(:)
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
        CALL rPulayK_mix(beta,W0AM,dime,NUH1,&
             & mixer%DXGL,mixer%DRGL,XL,RL,XN)


     ENDIF
     !-----------------------
     IF(Iter>0)THEN


        !ONE STEP BACKWARD
        DO IH=NUHMAX,2,-1
           mixer%DXGL(:,IH)=mixer%DXGL(:,IH-1)
           mixer%DRGL(:,IH)=mixer%DRGL(:,IH-1)
        END DO
        mixer%DXGL(:,1)=XL
        mixer%DRGL(:,1)=RL

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

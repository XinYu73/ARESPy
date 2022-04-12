# 1 "Chebyshev_fliter.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "Chebyshev_fliter.f90"
!#############################################################!
!        Chebyshev filter Method for diag(H)                  !
!*Author : Qiang Xu                                           !
!*Date   : 2017/11/28                                         !
!#############################################################!
MODULE chebyshev_module
   USE constants

   USE smpi_math_module

   USE parameters , ONLY : CheM,sn=>Nstates
   USE grid_module , ONLY : n,rho_calc,n1,n2,n3
   IMPLICIT NONE
   INTEGER(I4B),save :: iiii=1
CONTAINS
   !##############################################################!
   !*For   : k-space chebyshev filter                             !
   !*Date  : 2017/11/28                                           !
   !##############################################################!
   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !-----------------Chebyshev_filter subspace---------------------
   SUBROUTINE cheby_filter_RR(Ik,veff,X,D)
      USE parameters , ONLY : LRROrthNorm
      USE Lapack_module ,  ONLY : OrthNorm

      USE parameters, ONLY: Nstates_global

      IMPLICIT NONE
      !INOUT
      INTEGER(I4B),INTENT(IN) :: Ik
      REAL(DP),INTENT(IN) :: veff(:,:,:)
      COMPLEX(DP),INTENT(INOUT) :: X(:,:)
      REAL(DP),INTENT(OUT) :: D(:)  !rayleigh-ritz value
      !LOCAL
      REAL(DP) :: a,b,al

      INTEGER(I4B) :: bcastid,i,j

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !Boundary
      a=MAXVAL(D)
      al=MINVAL(D)



      !> find the core deal with the highest state
      aa:DO j=parallel%dims(1),1,-1 !> all core
      DO i=size(parallel%sub2sum,1),1,-1  !> all states per core
         IF(parallel%sub2sum(i,j)==Nstates_global)THEN
            bcastid=j-1
         CALL Estupb(7,Ik,veff,X(:,parallel%nstate_proc),b)
            exit aa
         ENDIF
      ENDDO
      ENDDO aa
      CALL MPI_BCAST(b,1,MPI_REAL8,bcastid,parallel%commy,mpinfo)

      !filtering
      CALL chebyshev_filter_scaled(Ik,veff,X,CheM,a,b,al)
      !CALL chebyshev_filter(kp,veff,X,CheM,a,b)
      !Rayleigh-Ritz step
      IF(LRROrthNorm)THEN
         !Standard RR

         if(parallel%isroot)print*,'Lorthnorm need to be false in mpi'
         stop

         CALL OrthNorm(X)
         CALL Rayleigh_Ritz(Ik,veff,X,D)
      ELSE
         !Generalize RR
         CALL GRayleigh_Ritz(Ik,veff,X,D)
      ENDIF
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE cheby_filter_RR
   !-----------------------GRayleigh_Ritz--------------------------
   SUBROUTINE GRayleigh_Ritz(Ik,veff,X,D)
      USE Lapack_module , ONLY : GeneralizeEigen,matmat

      USE parameters, ONLY: Nstates_global, BLOCK_MBNB
      USE ScaLapack_module, ONLY: SL_GeneralizeEigen,twoD_map,&
           & SL_matmat_gridcn,SL_matmat_gridnn
           ! & SL_matmat
      USE Grid_module, ONLY: global_n
      ! USE m_time_evaluate, ONLY: filename

      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B),INTENT(IN) :: Ik
      REAL(DP),INTENT(IN) :: veff(:,:,:)
      COMPLEX(DP),INTENT(INOUT) :: X(:,:)
      REAL(DP),INTENT(OUT) :: D(:)
      !



      COMPLEX(DP),DIMENSION(twoD_map(1,parallel%rankx,&
           &parallel%ranky),sn) :: S_hat,H_hat,Q

      COMPLEX(DP) :: Xnew(n,sn)

      INTEGER(I4B) :: mb,nb,cmb,cnb

      ! integer(i4b) :: i
!xq
!INTEGER(I4B) :: t1,t2,t3,t4,t5,t6
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !Calculate the overlap matrix
!call system_clock(t1)
      ! x=cmplx(0.d0,0.d0)






      ! do i=1,sn,1
      !    X(parallel%sub2sum(i,parallel%myid+1),i)=&
      !         & cmplx(parallel%sub2sum(i,parallel%myid+1),0)
      ! enddo
      mb=1
      nb=1!BLOCK_MBNB
      cmb=1!BLOCK_MBNB
      cnb=1!BLOCK_MBNB
      ! CALL SL_matmat('c','n',X,X,S_hat,&
      CALL SL_matmat_gridcn('c','n',X,X,S_hat,&
           &global_n,Nstates_global,global_n,Nstates_global,mb,nb,mb,nb,cmb,cnb)

!call system_clock(t2)
      !Calculate the project hamiltion
!call system_clock(t3)
      CALL Rayleigh_quotient(Ik,veff,sn,X,H_hat)
!call system_clock(t4)
      !solve the generalized eigenvalue problem
!call system_clock(t5)



      CALL SL_GeneralizeEigen(Nstates_global,H_hat,S_hat,Nstates_global&
           &,Nstates_global,Nstates_global,Nstates_global&
           &,cmb,cnb,cmb,cnb,Q,Nstates_global,Nstates_global&
           &,D,cmb,cnb)

!call system_clock(t6)
      !X=XQ



      ! CALL SL_matmat('n','n',X,Q,Xnew,global_n,Nstates_global&
      CALL SL_matmat_gridnn('n','n',X,Q,Xnew,global_n,Nstates_global&
           &,Nstates_global,Nstates_global,mb,nb,cmb,cnb,cmb,cnb)

!print*,'Overlap',(t2-t1)/10000.d0
!print*,'Hhat',(t4-t3)/10000.d0
!print*,'Genera',(t6-t5)/10000.d0
!STOP
      X=Xnew
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE GRayleigh_Ritz
   !----------------------Chebyshev_filter-------------------------
   SUBROUTINE chebyshev_filter(Ik,veff,X,m,a,b)
      !
      IMPLICIT NONE
      !INOUT
      INTEGER(I4B),INTENT(IN) :: Ik !k-point
      REAL(DP),INTENT(IN) :: veff(:,:,:) !veff
      COMPLEX(DP),INTENT(INOUT) :: X(:,:)
      INTEGER(I4B),INTENT(IN)   :: m ! the m degree Chebyshev polynomial we used
      REAL(DP),INTENT(IN) :: a,b  !interval [a,b] to be filter
      !LOCAL
      REAL(DP) :: e &   !(b-a)/2
            & ,   c    !(b+a)/2
      !
      COMPLEX(DP),DIMENSION(n,sn)  :: HV,Y,Ynew
      INTEGER(I4B) :: Ic
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF(m<2)THEN
         WRITE(6,*) 'Chebyshev_filter: the degree m must larger than 1'
      ENDIF

      e = ( b - a ) / 2
      c = ( b + a ) / 2

      CALL cal_HX(Ik,veff,sn,X,HV)
      Y(:,:)=( HV(:,:) - c*X(:,:) ) / e

      !Chebyshev filtering
      DO Ic=2,m

         !CALL HY
         CALL cal_HX(Ik,veff,sn,Y,HV)
         Ynew=2.d0*( HV-c*Y ) / e - X
         !store the Y

         X=Y

         Y=Ynew


      ENDDO
      !out put filted space
      X=Ynew
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE chebyshev_filter
   !-------------------chebyshev_filter_scaled---------------------
   SUBROUTINE chebyshev_filter_scaled(Ik,veff,X,m,a,b,al)
      !
      IMPLICIT NONE
      !INOUT
      INTEGER(I4B),INTENT(IN) :: Ik !k-point
      REAL(DP),INTENT(IN) :: veff(:,:,:) !veff
      COMPLEX(DP),INTENT(INOUT) :: X(:,:)
      INTEGER(I4B),INTENT(IN)   :: m ! the m degree Chebyshev polynomial we used
      REAL(DP),INTENT(IN) :: a,b,al  !interval [a,b] to be filter
      !LOCAL
      REAL(DP) :: e &   !(b-a)/2
            & ,   c &   !(b+a)/2
            & , sigma,sigmanew,tau
      !
      COMPLEX(DP),DIMENSION(n,sn)  :: HV,Y,Ynew
      INTEGER(I4B) :: Ic
      REAL(DP) :: temp1,temp2
!xq
!INTEGER(I4B) :: t1,t2
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF(m<2)THEN
         WRITE(6,*) 'Chebyshev_filter: the degree m must larger than 1'
      ENDIF

      e = ( b - a ) / 2
      c = ( b + a ) / 2
      sigma=e / (c-al)
      tau=2.d0 / sigma
!CALL system_clock(t1)
      CALL cal_HX(Ik,veff,sn,X,HV)
!CALL system_clock(t2)
!print*,'HX time:',(t2-t1)/10000.d0
      temp1=sigma / e
      Y= ( HV - c*X ) * temp1
      DO Ic=2,m
         sigmanew=1.d0 / (tau-sigma)
         CALL cal_HX(Ik,veff,sn,Y,HV)
         temp1=2.d0*sigmanew/e
         temp2=sigma*sigmanew
         Ynew=( HV - c*Y )*temp1 - temp2*X
         X=Y
         Y=Ynew
         sigma=sigmanew
      ENDDO

      X=Ynew
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE chebyshev_filter_scaled
   !-------------------------------HX------------------------------
   SUBROUTINE cal_HX(Ik,veff,nst,V,HV)
      USE matvec_module , ONLY : cmatvec
      IMPLICIT NONE
      !INOUT
      INTEGER(I4B),INTENT(IN) :: Ik,nst
      REAL(DP),INTENT(IN)  :: veff(:,:,:) !effective potential
      COMPLEX(DP),INTENT(IN)  :: V(:,:)
      COMPLEX(DP),INTENT(OUT) :: HV(:,:)
      !LOCAL
      INTEGER(I4B) :: Is
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! print*,'xlt-test     nst',nst
      DO Is=1,nst
      ! print*,'start',parallel%myid
         CALL cmatvec(veff,Ik,V(:,Is),HV(:,Is),n)
      ! print*,'end',parallel%myid
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE cal_HX
   !---------------------Rayleigh-quotient-------------------------
   SUBROUTINE Rayleigh_quotient(Ik,veff,nst,x,xhx)
     !xhx=(X,H,X)



      ! USE ScaLapack_module, ONLY: SL_matmat
      USE ScaLapack_module, ONLY: SL_matmat_gridcn
      USE parameters, ONLY: BLOCK_MBNB,Nstates_global
      USE Grid_module, ONLY: global_n

      use m_time_evaluate ,only: memory_sum,memory_free,filename
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN) :: Ik &
                              &,nst
      REAL(DP),INTENT(IN) :: veff(:,:,:)
      COMPLEX(DP) :: x(:,:),xhx(:,:)
      !LOCAL
      COMPLEX(DP) :: hx(n,nst)
      integer(I4B) :: i
      ! COMPLEX(DP) :: hx_global(global_n)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      call memory_sum("first_Rayleigh_quotient",real(size(hx),DP)*DP)
      CALL cal_HX(Ik,veff,nst,x,hx)
      !xhx(:,:)=MATMUL( TRANSPOSE(CONJG(x)) , hx )



      ! CALL SL_matmat('c','n',x,hx,xhx,global_n,Nstates_global,&
      CALL SL_matmat_gridcn('c','n',x,hx,xhx,global_n,Nstates_global,&
           &global_n,Nstates_global,1,BLOCK_MBNB,1,BLOCK_MBNB,&
           & BLOCK_MBNB,BLOCK_MBNB)

      call memory_free("first_Rayleigh_quotient",real(size(hx),DP)*DP)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! #ifndef 1
!        open(1111,file="Hhat")
!        do i=1,size(x,2),1
!           write(1111,'(10000(F25.16,1X))')real(hx(:,i))
!        enddo
!        close(1111)
!        stop 'serial'
! #else
!       open(1111+parallel%myid,file="Hx"//filename(10:11))
!       ! call MPI_ALLGATHERV(hx(:,1),parallel%mygrid_range(3)&
!       !      &,MPI_COMPLEX16,hx_global,parallel%recvcounts&
!       !      &,parallel%displs,MPI_COMPLEX16&
!       !      & ,parallel%commx,mpinfo)
!       ! write(1111+parallel%myid,*)real(hx_global)

!       do i=1,size(xhx,2),1
!          ! write(1111+parallel%myid,'(10000(F25.16,1X))')real(xhx(:,i))
!          write(1111+parallel%myid,*)xhx(:,i)
!       enddo
!       close(1111+parallel%myid)
!       call MPI_BARRIER(parallel%comm,mpinfo)
!       stop 'parallel'
! #endif
   ENDSUBROUTINE Rayleigh_quotient
   !---------------------Rayleigh-Ritz step------------------------
   SUBROUTINE Rayleigh_Ritz(Ik,veff,X,D)
      USE Lapack_module ,  ONLY : diagM,matmat
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN) :: Ik
      REAL(DP),INTENT(IN) :: veff(:,:,:)
      COMPLEX(DP),INTENT(INOUT) :: X(:,:)
      REAL(DP),INTENT(OUT) :: D(:)
      !LOCAL
      COMPLEX(DP),DIMENSION(sn,sn) :: Hhat,Q
      COMPLEX(DP) :: Xnew(n,sn)
      INTEGER(I4B) :: I
!xq
!integer(I4B) :: t1,t2,t3,t4,t5,t6
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !F1:Hhat=(X,H,X)
!call system_clock(t1)
      CALL Rayleigh_quotient(Ik,veff,sn,X,Hhat)
!call system_clock(t2)
      !F2:eigen-decomposition Q,D
!call system_clock(t3)
      CALL diagM(Hhat,Q,D)
!call system_clock(t4)
      !F3:X=XQ
!call system_clock(t5)
      !X=MATMUL( X , Q )
      CALL matmat(X,Q,'N','N',Xnew)
      X=Xnew
!call system_clock(t6)
!print*,'R-Q time',(t2-t1)/10000.d0
!print*,'E-D time',(t4-t3)/10000.d0
!print*,'MAT time',(t6-t5)/10000.d0
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE
   !-------------------------upper bound---------------------------
   SUBROUTINE Estupb(k,Ik,veff,vec,b)
      !
      USE matvec_module , ONLY : cmatvec
      USE Lapack_module , ONLY : Norm_2,diagM
      IMPLICIT NONE
      !INOUT
      INTEGER(I4B),INTENT(IN) :: k , Ik
      !hamiltonian
      REAL(DP),INTENT(IN) :: veff(:,:,:)
      COMPLEX(DP),INTENT(IN) :: vec(:)
      REAL(DP),INTENT(OUT) :: b
      !LOCAL
      COMPLEX(DP) :: v0(n),v(n),f(n),alpha
      REAL(DP) :: beta,eval(k),mz
      COMPLEX(DP) :: T(k,k),evec(k,k)
      INTEGER(I4B) :: J,Nj
      COMPLEX(DP),external  :: Zdotc

      COMPLEX(DP) :: alpha_local,b_local
      REAL(DP) :: beta_local

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      v=vec !should normalize
      T(:,:)=0.d0
      Nj=MIN(k,10)
      !
      CALL cmatvec(veff,Ik,v,f,n)



      alpha_local=zdotc(size(f),f,1,v,1)
      call mpi_allreduce(alpha_local,alpha,1,mpi_complex16&
           & ,mpi_sum, parallel%commx, mpinfo)

      f=f-alpha*v
      T(1,1)=alpha
      DO J=2,Nj
         ! beta=SQRT(REAL(DOT_PRODUCT(f,f),8))



         beta_local=zdotc(size(f),f,1,f,1)
         call mpi_allreduce(beta_local, beta,1,mpi_real8&
              &, mpi_sum,parallel%commx,mpinfo)

         beta=sqrt(beta)
         v0=v
         v=f/beta
         CALL cmatvec(veff,Ik,v,f,n)
         f=f-beta*v0
         ! alpha=DOT_PRODUCT(f,v)



      alpha_local=zdotc(size(f),f,1,v,1)
      call mpi_allreduce(alpha_local,alpha,1,mpi_complex16&
           & ,mpi_sum, parallel%commx, mpinfo)

         f=f-alpha*v
         T(J,J-1)=beta
         T(J-1,J)=beta
         T(J,J)=alpha
      ENDDO
      !NORM2(T)
      !b=Norm_2(T,k) + SQRT(REAL(DOT_PRODUCT(f,f),8))
      CALL diagM(T,evec,eval)
      !mz=MAX( ABS(evec(k,k)), ABS(evec(k,k-1)) , ABS(evec(k,k-2)) )
      mz=ABS(evec(k,k))



      b_local=zdotc(size(f),f,1,f,1)
      call mpi_allreduce(b_local,b,1,mpi_complex16,mpi_sum&
           &, parallel%commx, mpinfo)

      b=eval(k) + SQRT(REAL(b,DP))*mz
      !STOP
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE Estupb
   !##############################################################!
   !*First :                Random Method                         !
   !         Ref: Y. K. Zhou etc , JCP , 274(2014)770-782         !
   !##############################################################!
   !-----------------------init up low bound-----------------------
   SUBROUTINE inituplow(k,Ik,veff,v,eval,a,b,al)
      !
      USE matvec_module , ONLY : cmatvec
      USE Lapack_module , ONLY : diagM
      USE math , ONLY : Norm
      IMPLICIT NONE
      !INOUT
      INTEGER(I4B),INTENT(IN) :: k,Ik
      !hamiltonian
      REAL(DP),INTENT(IN) :: veff(:,:,:)
      COMPLEX(DP),INTENT(INOUT) :: v(:)
      REAL(DP),INTENT(OUT) :: a,b,al,eval(:)
      !LOCAL
      COMPLEX(DP) :: v0(n),f(n),alpha
      REAL(DP) :: beta, &
              &   fbeta=0.5d0
      COMPLEX(DP) :: T(k,k),evec(k,k)
      INTEGER(I4B) :: J,Nj
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      T(:,:)=0.d0
      Nj=MIN(k,10)
      !
      CALL cmatvec(veff,Ik,v,f,n)
      alpha=DOT_PRODUCT(f,v)
      f=f-alpha*v
      T(1,1)=alpha
      DO J=1,Nj
         beta=SQRT(REAL(DOT_PRODUCT(f,f),8))
         v0=v
         v=f/beta
         CALL cmatvec(veff,Ik,v,f,n)
         f=f-beta*v0
         alpha=DOT_PRODUCT(f,v)
         f=f-alpha*v
         T(J,J-1)=beta
         T(J-1,J)=beta
         T(J,J)=alpha
      ENDDO
      !Rayleigh-Ritz value
      CALL diagM(T,evec,eval)
      a=fbeta*eval(1)+(1.d0-fbeta)*eval(k)
      al=eval(1)
      b=eval(k) + SQRT(REAL(DOT_PRODUCT(f,f),8))*ABS(evec(k,k))!+1e-10
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE inituplow
   !------------------------first_step_filter----------------------
   SUBROUTINE first_SCFstep_filter(Ik,veff,X,eval)
      USE parameters , ONLY : CF0=>CheM0,LRROrthNorm
      USE Lapack_module , ONLY : diagM,OrthNorm,matmat &
                       &, GeneralizeEigen
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN) :: Ik
      REAL(DP),INTENT(IN) :: veff(:,:,:)
      COMPLEX(DP),INTENT(INOUT) :: X(:,:)
      REAL(DP),INTENT(OUT) :: eval(:)
      !LOCAL
      REAL(DP) :: a,b,al
      INTEGER(I4B) :: I,Niter=4
      REAL(DP) :: evald(sn),deval,TOL=1e-8
      COMPLEX(DP),DIMENSION(sn,sn)  :: Hhat,Q,Shat
      COMPLEX(DP)  :: Xnew(n,sn),vec(n)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !low up bound
      vec(:)=X(:,sn)
      CALL inituplow(7,Ik,veff,vec,eval,a,b,al)
      !
      evald=eval
      DO I=1,Niter
         !
         CALL chebyshev_filter_scaled(Ik,veff,X,CF0,a,b,al)
         IF(LRROrthNorm)THEN
            !OrthNorm
            CALL OrthNorm(X)
            !xHx
            CALL Rayleigh_quotient(Ik,veff,sn,X,Hhat)
            !eigen-decomposion
            CALL diagM(Hhat,Q,eval)
         ELSE
            !Overlap matrix
            CALL matmat(X,X,'C','N',Shat)
            !projected hamiltonian
            CALL Rayleigh_quotient(Ik,veff,sn,X,Hhat)
            !eigen-decomposion
            CALL GeneralizeEigen(sn,Hhat,Shat,Q,eval)
         ENDIF
         deval=SUM(ABS(eval-evald))/sn
         !-----------------
         IF(deval<TOL) EXIT
         !-----------------
         !store old eigenvalue
         evald(:)=eval(:)
         !update the new bound
         a=eval(sn)
         al=eval(1)
      !psi=MATMUL( psi , Q )
      ENDDO

      !rotation
      CALL matmat(X,Q,'N','N',Xnew)
      X=Xnew
      !PRINT*,'Free diag>>>iter:',I
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE first_SCFstep_filter
   !-------------------------random subspace-----------------------
   SUBROUTINE randomfast(psi,nev,radius)
      USE grid_module , ONLY : n
      USE Lapack_module , ONLY : OrthNorm
      !
      IMPLICIT NONE
      !INOUT
      COMPLEX(DP),INTENT(OUT) :: psi(:,:)
      REAL(DP),INTENT(IN) :: radius
      INTEGER(I4B),INTENT(IN) :: nev
      !LOCAL
      REAL(DP) :: rant,rantr
      INTEGER(I4B) :: Is,Id
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CALL random_seed()
      DO Is=1,nev
         DO Id=1,n
            CALL random_number(rant)
            CALL random_number(rantr)
            rant = -radius+rant*2.d0*radius
            rantr= -radius+rantr*2.d0*radius
            psi(Id,Is)=CMPLX(rant,rantr)
         ENDDO
      ENDDO
      CALL OrthNorm(psi)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE randomfast
   !------------------first_step spin-------------------------
   SUBROUTINE first_CheSCF_random(rhoS,nev,psi,eval)
      USE parameters , ONLY : Nspin
      USE potential_module , ONLY : cal_veff
      USE grid_module , ONLY : n,n1,n2,n3
      IMPLICIT NONE
      !IN/OUT
      REAL(DP),INTENT(IN) :: rhoS(:,:,:,:)
      INTEGER(I4B),INTENT(IN) :: nev
      COMPLEX(DP),INTENT(OUT) :: psi(:,:,:,:)
      REAL(DP),INTENT(OUT) :: eval(:,:,:)
      !LOCAL
      INTEGER(I4B) :: Is
      REAL(DP) :: veff(n1,n2,n3,NSPIN)
      COMPLEX(DP) :: psi_ran(n,nev)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !
      CALL cal_veff(rhoS,veff)
      !initialize the subspace
      CALL randomfast(psi_ran,nev,50.d0)
      !
      DO Is=1,Nspin
         CALL first_CheSCF(veff(:,:,:,Is),psi_ran,nev,psi(:,:,:,Is),eval(:,:,Is))
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE first_CheSCF_random
   !----------------------first step----------------------------
   SUBROUTINE first_CheSCF(veff,psi_ran,nev,psi,eval)
      USE grid_module , ONLY : nk,KPT
      USE struct_module , ONLY : ncharge
      IMPLICIT NONE
      !IN/OUT
      REAL(DP),INTENT(IN) :: veff(:,:,:)
      COMPLEX(DP),INTENT(IN) :: psi_ran(:,:)
      COMPLEX(DP),INTENT(OUT) :: psi(:,:,:)
      REAL(DP),INTENT(OUT) :: eval(:,:)
      INTEGER(I4B),INTENT(IN) :: nev
      !LOCAL
      REAL(DP) :: a,b,al
      INTEGER(I4B) :: Ik
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !
      DO Ik=1,nk
         psi(:,:,Ik)=psi_ran(:,:)
         CALL first_SCFstep_filter(Ik,veff,psi(:,:,Ik),eval(:,Ik))
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE first_CheSCF
   !##############################################################!
   !*First :                STO+PW Method                         !
   !##############################################################!
   !------------initialize_subspace by random STO-------------
   SUBROUTINE cheby_init_sto(Nmax,Npw,initX_sto)
      !initialize the subsystem by random STO
      USE parameters , ONLY : Nstates,Nstates_global
      USE struct_module , ONLY : struct,lat_mat,naty
      USE grid_module , ONLY : grid,n,n1,n2,n3
      USE math , ONLY : atom_sto
      USE Lapack_module , ONLY : OrthNorm
      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B),INTENT(IN) :: Nmax & !max dimen 'super'-subspace
                             &, Npw  !number of plane wave used
      COMPLEX(DP),INTENT(OUT) :: initX_sto(:,:)
      !LOCAL
      INTEGER(I4B) :: Ity,Ia,Ip,icx,icy,icz,Ii &
            & ,nrep=1 &



            &,Il,Im,Nstart &
            &,Is,nassign

      REAL(DP) :: ra(3),rat(4),f,gdotr
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      initX_sto(:,:)=cmplx(0.d0,0.d0)
      !all states
      Ii=0

      !> init the count independently in per process
      nassign=0

      !all type
      DO Ity=1,naty
         !all atom
         DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
            !store the position
            ra(:)=struct%poscar(:,Ia)
            !call the angle momentum
            DO Il=0,struct%Lmax(Ity)
               DO Im=-Il,Il
                  !call the states
                  Ii=Ii+1
                  if(Ii.gt.Nstates_global)Ii=1

                  do Is=1,parallel%nstate_proc,1
                  If(Ii==parallel%sub2sum(Is,parallel%ranky+1))THEN
                     nassign=nassign+1
                     IF(nassign.gt.parallel%nstate_proc)nassign=1

                  !all points
                  DO Ip=1,n
                     !replicas
                     DO icz=-nrep,nrep
                     DO icy=-nrep,nrep
                     DO icx=-nrep,nrep
                        rat(1:3)= icx*lat_mat(:,1)  &
                          & +  icy*lat_mat(:,2)  &
                          & +  icz*lat_mat(:,3)  &
                          & +  grid%rVec(1:3,Ip)  &
                          & -  ra(:)
                        rat(4)=SQRT( rat(1)**2 + rat(2)**2 +rat(3)**2 )
                        !STO
                    CALL atom_STO(struct%prinq(Il+1,Ity),Il,Im,struct%zeta(Il+1,Ity),rat,f)
                    !init



                    initX_sto(Ip,nassign)=initX_sto(Ip,nassign)+CMPLX(f,0.d0)

                     ENDDO
                     ENDDO
                     ENDDO
                     !
                  ENDDO

                  ENDIF
                  ENDDO !> Is

                  !
               ENDDO
            ENDDO
            !
         ENDDO
         !
      ENDDO
! OPEN(file_unit,FILE=filename(8:11)//'test.dat')
! WRITE(file_unit,*) '14 14 14'
! WRITE(file_unit,*) initX_sto
! CLOSE(file_unit)
! CALL MPI_finalize(mpinfo)
! stop
      !For plane-wave if needed
      IF(Npw>0)THEN
        Nstart=Nmax-Npw+1
        DO Ii=Nstart,Nmax
           DO Ip=1,n

              gdotr=DOT_PRODUCT(grid%gVec(1:3,parallel%sub2sum(Ii,parallel%ranky+1)),grid%rVec(1:3,Ip))



              initX_sto(Ip,Ii)=EXP(IMAG*gdotr)
           ENDDO
        ENDDO
        !orth-Norm if need
        CALL OrthNorm(initX_sto)
      ENDIF
      !orth-Norm if need
      !CALL OrthNorm(X)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE cheby_init_sto
   !---------------------cheby_subspace----------------------------
   SUBROUTINE first_CheSCF_sto(veff,nev,psi,eval)
      USE parameters , ONLY : Nspin
      ! USE potential_module , ONLY : cal_veff
      USE grid_module , ONLY : n,n1,n2,n3,nk
      USE struct_module , ONLY : naty,struct
      use m_time_evaluate ,only: memory_sum,memory_free

      USE smpi_math_module, ONLY: states_split,parallel,mpinfo
      ! USE ScaLapack_module, ONLY: twoD_map

      IMPLICIT NONE
      !INOUT
      REAL(DP),INTENT(IN) :: veff(:,:,:,:)
      INTEGER(I4B),INTENT(IN) :: nev
      COMPLEX(DP) :: psi(:,:,:,:)
      REAL(DP)  :: eval(:,:,:)
      !LOCAL
      INTEGER(I4B) :: Is,Ik,Ity &
               &,Nmax,Npw
      ! REAL(DP) :: veff(n1,n2,n3,Nspin)
      COMPLEX(DP),ALLOCATABLE :: initX(:,:)!,Shat(:,:),Hhat(:,:)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !Calculate the effective potential
      ! CALL cal_veff(rhoS,veff)
      call memory_sum("first_Chebyshev_local",real(size(veff),DP)*DP)
      !initialize the temp subspace
      !estimate number of the states
      Nmax=0
      Npw=0
      DO Ity=1,naty
         Nmax=Nmax+struct%nati(Ity)*(struct%Lmax(Ity)+1)**2
      ENDDO

      call states_split(Nmax)
      ! NMAX=twoD_map(2,parallel%rankx,parallel%ranky)

      !test
      IF(Nmax<=nev)THEN
         Npw=nev-Nmax
         Nmax=nev

         PRINT*,'[Add plane wave for subspace in rank',parallel%myid,']',Npw



      ELSE
         Nmax=nev
      ENDIF
      !allocate
      ALLOCATE(initX(n,nev))
      call memory_sum("first_Chebyshev",real(size(initX),DP)*DP)
      !init a 'super'-subspace
      CALL cheby_init_sto(nev,Npw,initX)
      !Raleigh-Ritz step
      !all spin
      DO Is=1,Nspin
         !all kpoints
         DO Ik=1,nk
            CALL first_subspace_stopw(Ik,veff(:,:,:,Is), &
               &   Nmax,nev,initX,psi(:,:,Ik,Is),eval(:,Ik,Is))
            !filter
        !    IF(Npw>0)THEN
        !       CALL first_SCFstep_filter(Ik,veff(:,:,:,Is), &
        !       &   psi(:,:,Ik,Is),eval(:,Ik,Is))
        !    ENDIF
         ENDDO
      ENDDO
      call memory_free("first_Chebyshev",real(size(initX),DP)*DP)
      DEALLOCATE(initX)
      call memory_free("first_Chebyshev_local",real(size(veff),DP)*DP)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ENDSUBROUTINE first_CheSCF_sto
   !----------------------------------------------------------
   SUBROUTINE first_subspace_stopw(Ik,veff,Nmax,nev,initX,X,D)
      USE Lapack_module , ONLY : diagM,OrthNorm,matmat &
                       &, GeneralizeEigen
      ! use m_time_evaluate, only:memory_sum,memory_free
      use m_time_evaluate, only:memory_sum,memory_free,filename

      use parameters, only: Nstates_global, BLOCK_MBNB
      use ScaLapack_module, only: SL_GeneralizeEigen &
           & ,SL_matmat_gridcn, SL_matmat_gridnn,twoD_map
           ! &, SL_matmat,twoD_map
      USE Grid_module, ONLY: global_n

      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B),INTENT(IN) :: Ik &
                             &, Nmax,nev
      REAL(DP),INTENT(IN) :: veff(:,:,:)
      COMPLEX(DP),INTENT(INOUT)  :: initX(:,:)
      COMPLEX(DP),INTENT(OUT) :: X(:,:) !subspace
      REAL(DP),INTENT(OUT) :: D(:)      !eigenvalue
      !LOCAL



      COMPLEX(DP),DIMENSION(twoD_map(1,parallel%rankx,&
           &parallel%ranky),Nmax) :: Shat,Hhat,Qs

      COMPLEX(DP) :: Xnew(n,Nmax)
      REAL(DP) :: Ds(Nmax)

      INTEGER(I4B) :: mb,nb,cmb,cnb
      INTEGER(I4B) :: ix,iy
      ! COMPLEX(DCP) :: initx_global(global_n)
      REAL(DP) :: veff_global(global_n)

      integer(I4B) :: i
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !Raleigh-Ritz step
      !Overlap matrix
      call memory_sum("first_subspace_sto_local",(real(size(Shat),DP)*3+size(Xnew)+size(Ds))*DP)



      mb=1 !n1*n2*n3
      nb=1!BLOCK_MBNB
      cmb=1!BLOCK_MBNB
      cnb=1!BLOCK_MBNB
      ! print *,'shape initX',shape(initX)
      ! CALL SL_matmat('c','n',initX,initX,Shat,&
      CALL SL_matmat_gridcn('c','n',initX,initX,Shat,&
           &global_n,Nstates_global,global_n,Nstates_global,mb,nb,mb,nb,cmb,cnb)
      ! ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
      ! iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
      ! call MPI_ALLGATHERV(veff(ix,iy,1),parallel%mygrid_range(3)&
      !      &,MPI_REAL8,veff_global,parallel%recvcounts&
      !      &,parallel%displs,MPI_REAL8&
      !      & ,parallel%commx,mpinfo)
      ! ! ! print *,'shape initX2',shape(initX)
      ! open(1111+parallel%myid,file="veff_cheb"//filename(10:11))
      ! write(1111+parallel%myid,*)veff_global
      ! do i=1,size(initX,2),1
      !    write(1111+parallel%myid,'(10000(F25.16,1X))')real(initX(:,i))
      ! enddo
      ! close(1111+parallel%myid)
      ! open(1111+parallel%myid,file="Shat_p"//filename(10:11))
      ! do i=1,size(Shat,2),1
      !    write(1111+parallel%myid,'(1000(F25.16,2X))')real(Shat(:,i))
      ! enddo
      ! close(1111+parallel%myid)
      ! call smpi_exit()
      ! CALL MPI_finalize(mpinfo)
      ! stop

      !projected hamiltonian
      CALL Rayleigh_quotient(Ik,veff,Nmax,initX,Hhat)
      !eigen-decomposion



      CALL SL_GeneralizeEigen(Nstates_global,Hhat,Shat,Nstates_global&
           &,Nstates_global,Nstates_global,Nstates_global&
           &,cmb,cnb,cmb,cnb,Qs,Nstates_global,Nstates_global&
           &,D,cmb,cnb)

! #ifndef 1
!        open(1111,file="Hhat")
!        do i=1,size(Hhat,2),1
!           write(1111,'(10000(F25.15,1X))')real(Hhat(:,i))
!        enddo
!           write(1111,*)D
!        close(1111)
!        stop
! #else
!      open(1111+parallel%myid,file="Hhat"//filename(10:11))
!      do i=1,size(Hhat,2),1
!         write(1111+parallel%myid,'(10000(F25.15,1X))')real(Hhat(:,i))
!      enddo
!          write(1111+parallel%myid,*)D
!      close(1111+parallel%myid)
!      call smpi_exit()
! #endif
      !rotation



      ! CALL SL_matmat('n','n',initX,Qs,X,global_n,Nstates_global&
      CALL SL_matmat_gridnn('n','n',initX,Qs,X,global_n,Nstates_global&
           &,Nstates_global,Nstates_global,mb,nb,cmb,cnb,cmb,cnb)
     ! open(1111+parallel%myid,file="Hhat"//filename(10:11))
     ! do i=1,size(X,2),1
     !    write(1111+parallel%myid,'(10000(F25.15,1X))')real(X(:,i))
     ! enddo
     !     write(1111+parallel%myid,*)D
     ! close(1111+parallel%myid)
     ! call smpi_exit()

      !eigen-value
      !D(1:nev)=Ds(1:nev)
      !X(:,1:nev)=Xnew(:,1:nev)
      call memory_free("first_subspace_sto_local",(real(size(Shat),DP)*3+size(Xnew)+size(Ds))*DP)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE first_subspace_stopw
   !##############################################################!
   !*For    :              BvK Cell Method First Step             !
   !*Author : Qiang Xu                                            !
   !*Date   : 2018-03-05                                          !
   !*First  :               STO+PW Method                         !
   !##############################################################!
   !---------------------cheby_subspace----------------------------
   SUBROUTINE BvK_first_CheSCF_sto_rand(rhoS,nev,psi,eval)
      USE parameters , ONLY : Nspin
      USE potential_module , ONLY : cal_veff
      USE grid_module , ONLY : n,n1,n2,n3,nk
      USE struct_module , ONLY : naty,struct
      IMPLICIT NONE
      !INOUT
      REAL(DP),INTENT(IN) :: rhoS(:,:,:,:)
      INTEGER(I4B),INTENT(IN) :: nev
      REAL(DP) :: psi(:,:,:)
      REAL(DP)  :: eval(:,:)
      !LOCAL
      INTEGER(I4B) :: Is,Ity &
               &,Nmax,Nrand
      REAL(DP) :: veff(n1,n2,n3,Nspin)
      REAL(DP),ALLOCATABLE :: initX(:,:)!,Shat(:,:),Hhat(:,:)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !Calculate the effective potential
      CALL cal_veff(rhoS,veff)
      !initialize the temp subspace
      !estimate number of the states
      Nmax=0
      Nrand=0
      DO Ity=1,naty
         Nmax=Nmax+struct%nati(Ity)*(struct%Lmax(Ity)+1)**2
      ENDDO
      !test
      IF(Nmax<=nev)THEN
         Nrand=nev-Nmax
         Nmax=nev
         PRINT*,'[Add random wave for subspace]',Nrand
      ENDIF
      !allocate
      ALLOCATE(initX(n,Nmax))
      !init a 'super'-subspace
      CALL BvK_cheby_init_sto_rand(Nmax,Nrand,initX)
      !Raleigh-Ritz step
      !all spin
      DO Is=1,Nspin
         CALL first_subspace_sto_rand(veff(:,:,:,Is), &
            &   Nmax,nev,initX,psi(:,:,Is),eval(:,Is))
      ENDDO
      DEALLOCATE(initX)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE BvK_first_CheSCF_sto_rand
   !---------------------cheby_subspace----------------------------
   SUBROUTINE BvK_cheby_init_sto_rand(Nmax,Nrand,initX_sto)
      !initialize the subsystem by random STO
      USE parameters , ONLY : Nstates
      USE struct_module , ONLY : struct,lat_mat,naty
      USE grid_module , ONLY : grid,n,n1,n2,n3
      USE math , ONLY : atom_sto
      USE Lapack_module , ONLY : OrthNorm_real
      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B),INTENT(IN) :: Nmax & !max dimen 'super'-subspace
                             &, Nrand  !number of random wave used
      REAL(DP),INTENT(OUT) :: initX_sto(:,:)
      !LOCAL
      INTEGER(I4B) :: Ity,Ia,Ip,icx,icy,icz,Ii &
            & ,nrep=1 &
            &,Il,Im,Nstart
      REAL(DP) :: ra(3),rat(4),f,randt
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      initX_sto(:,:)=0.d0
      !all states
      Ii=0
      !all type
      DO Ity=1,naty
         !all atom
         DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
            !store the position
            ra(:)=struct%poscar(:,Ia)
            !call the angle momentum
            DO Il=0,struct%Lmax(Ity)
               DO Im=-Il,Il
                  !call the states
                  Ii=Ii+1
                  !all points
                  DO Ip=1,n
                     !replicas
                     DO icz=-nrep,nrep
                     DO icy=-nrep,nrep
                     DO icx=-nrep,nrep
                        rat(1:3)= icx*lat_mat(:,1)  &
                          & +  icy*lat_mat(:,2)  &
                          & +  icz*lat_mat(:,3)  &
                          & +  grid%rVec(1:3,Ip)  &
                          & -  ra(:)
                        rat(4)=SQRT( rat(1)**2 + rat(2)**2 +rat(3)**2 )
                        !STO
                        CALL atom_STO(struct%prinq(Il+1,Ity),Il,Im,struct%zeta(Il+1,Ity),rat,f)
                        !init
                        initX_sto(Ip,Ii)=initX_sto(Ip,Ii)+f
                     ENDDO
                     ENDDO
                     ENDDO
                     !
                  ENDDO
                  !
               ENDDO
            ENDDO
            !
         ENDDO
         !
      ENDDO
!OPEN(10086,FILE='test.dat')
!WRITE(10086,*) '14 14 14'
!WRITE(10086,*) REAL(initX_sto(:,1),DP)
!CLOSE(10086)
!STOP
      !For plane-wave if needed
      IF(Nrand>0)THEN
        CALL random_seed()
        Nstart=Nmax-Nrand+1
        DO Ii=Nstart,Nmax
           DO Ip=1,n
              CALL random_number(randt)
              initX_sto(Ip,Ii)=randt
           ENDDO
        ENDDO
        !orth-Norm if need
        CALL OrthNorm_real(initX_sto)
      ENDIF
      !orth-Norm if need
      !CALL OrthNorm(X)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE BvK_cheby_init_sto_rand
   !----------------------------------------------------------
   SUBROUTINE first_subspace_sto_rand(veff,Nmax,nev,initX,X,D)
      USE Lapack_module , ONLY : matmat_real &
                       &, GeneralizeEigen_real
      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B),INTENT(IN) :: Nmax,nev
      REAL(DP),INTENT(IN) :: veff(:,:,:)
      REAL(DP),INTENT(IN)  :: initX(:,:)
      REAL(DP),INTENT(OUT) :: X(:,:) !subspace
      REAL(DP),INTENT(OUT) :: D(:)      !eigenvalue
      !LOCAL
      REAL(DP),DIMENSION(Nmax,Nmax) :: Shat,Hhat,Qs
      REAL(DP) :: Xnew(n,Nmax)
      REAL(DP) :: Ds(Nmax)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !Raleigh-Ritz step
      !Overlap matrix
      CALL matmat_real(initX,initX,'T','N',Shat)
      !projected hamiltonian
      CALL Rayleigh_quotient_real(veff,Nmax,initX,Hhat)
      !eigen-decomposion
      CALL GeneralizeEigen_real(Nmax,Hhat,Shat,Qs,Ds)
      !rotation
      CALL matmat_real(initX,Qs,'N','N',Xnew)
      !eigen-value
      D(1:nev)=Ds(1:nev)
      X(:,1:nev)=Xnew(:,1:nev)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE first_subspace_sto_rand
   !##############################################################!
   !*For    :     BvK Cell Method Filtering (Standard RR)         !
   !*Author : Qiang Xu                                            !
   !*Date   : 2018-03-05                                          !
   !##############################################################!
   !---------------------Rayleigh-quotient-------------------------
   SUBROUTINE Rayleigh_quotient_real(veff,nst,x,xhx)
      !xhx=(X,H,X)
      USE Lapack_module , ONLY : matmat_real
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN) :: nst
      REAL(DP),INTENT(IN) :: veff(:,:,:)
      REAL(DP) :: x(:,:),xhx(:,:)
      !LOCAL
      REAL(DP) :: hx(n,nst)
      integer(I4b) :: i
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! print *,'output H --> dimensions: ',n,n
      ! open(1122,file="Hamiltonian")
      ! write(1122,*)'dimensions',n, n
      ! do i = 1,n,1
      !    x=0.d0
      !    x(i,1)=1.d0
      !    CALL cal_HX_real(veff,1,x(:,1:1),hx(:,1:1))
      !    write(1122,'(20000(G16.9,X))')hx(:,1)
      ! enddo
      ! close(1122)
      ! stop "out Hamiltonian end"
      CALL cal_HX_real(veff,nst,x,hx)
      !xhx(:,:)=MATMUL( TRANSPOSE(CONJG(x)) , hx )
      CALL matmat_real(x,hx,'T','N',xhx)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE Rayleigh_quotient_real
   !-------------------------------HX------------------------------
   SUBROUTINE cal_HX_real(veff,nst,V,HV)
      USE matvec_module , ONLY : rmatvec
      IMPLICIT NONE
      !INOUT
      INTEGER(I4B),INTENT(IN) :: nst
      REAL(DP),INTENT(IN)  :: veff(:,:,:) !effective potential
      REAL(DP),INTENT(IN)  :: V(:,:)
      REAL(DP),INTENT(OUT) :: HV(:,:)
      !LOCAL
      INTEGER(I4B) :: Is,ts,tt
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      DO Is=1,nst
         CALL system_clock(ts)
         CALL rmatvec(veff,V(:,Is),HV(:,Is),n)
         CALL system_clock(tt)
      ENDDO
      !=========================================
      !print *,"nst",nst
      !print *,'cal HX time -->',(ts-tt)/10000.d0
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE cal_HX_real
   !-------------------------upper bound---------------------------
   SUBROUTINE Estupb_real(k,veff,vec,b)
      !
      USE matvec_module , ONLY : rmatvec
      USE Lapack_module , ONLY : diagM_real
      IMPLICIT NONE
      !INOUT
      INTEGER(I4B),INTENT(IN) :: k
      !hamiltonian
      REAL(DP),INTENT(IN) :: veff(:,:,:)
      REAL(DP),INTENT(IN) :: vec(:)
      REAL(DP),INTENT(OUT) :: b
      !LOCAL
      REAL(DP) :: v0(n),v(n),f(n),alpha
      REAL(DP) :: beta,eval(k),mz
      REAL(DP) :: T(k,k),evec(k,k)
      INTEGER(I4B) :: J,Nj
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      v=vec !should normalize
      T(:,:)=0.d0
      Nj=MIN(k,10)
      !
      CALL rmatvec(veff,v,f,n)
      alpha=DOT_PRODUCT(f,v)
      f=f-alpha*v
      T(1,1)=alpha
      DO J=2,Nj
         beta=SQRT(REAL(DOT_PRODUCT(f,f),8))
         v0=v
         v=f/beta
         CALL rmatvec(veff,v,f,n)
         f=f-beta*v0
         alpha=DOT_PRODUCT(f,v)
         f=f-alpha*v
         T(J,J-1)=beta
         T(J-1,J)=beta
         T(J,J)=alpha
      ENDDO
      !NORM2(T)
      !b=Norm_2(T,k) + SQRT(REAL(DOT_PRODUCT(f,f),8))
      CALL diagM_real(T,evec,eval)
      !mz=MAX( ABS(evec(k,k)), ABS(evec(k,k-1)) , ABS(evec(k,k-2)) )
      mz=ABS(evec(k,k))
      b=eval(k) + SQRT(REAL(DOT_PRODUCT(f,f),8))*mz
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE Estupb_real
   !-------------------chebyshev_filter_scaled---------------------
   SUBROUTINE chebyshev_filter_scaled_real(veff,X,m,a,b,al)
      !
      IMPLICIT NONE
      !INOUT
      REAL(DP),INTENT(IN) :: veff(:,:,:) !veff
      REAL(DP),INTENT(INOUT) :: X(:,:)
      INTEGER(I4B),INTENT(IN)   :: m ! the m degree Chebyshev polynomial we used
      REAL(DP),INTENT(IN) :: a,b,al  !interval [a,b] to be filter
      !LOCAL
      REAL(DP) :: e &   !(b-a)/2
            & ,   c &   !(b+a)/2
            & , sigma,sigmanew,tau
      !
      REAL(DP),DIMENSION(n,sn)  :: HV,Y,Ynew
      INTEGER(I4B) :: Ic
      REAL(DP) :: temp1,temp2
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF(m<2)THEN
         WRITE(6,*) 'Chebyshev_filter: the degree m must larger than 1'
      ENDIF

      !===========================
      !print*,"{-------------------"
      !===========================
      e = ( b - a ) / 2
      c = ( b + a ) / 2
      sigma=e / (c-al)
      tau=2.d0 / sigma
      CALL cal_HX_real(veff,sn,X,HV)
      temp1=sigma / e
      Y= ( HV - c*X ) * temp1
      DO Ic=2,m
         sigmanew=1.d0 / (tau-sigma)
         CALL cal_HX_real(veff,sn,Y,HV)
         temp1=2.d0*sigmanew/e
         temp2=sigma*sigmanew
         Ynew=( HV - c*Y )*temp1 - temp2*X
         X=Y
         Y=Ynew
         sigma=sigmanew
      ENDDO
      !===========================
      !print*,"-------------------}"
      !===========================
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE chebyshev_filter_scaled_real
   !-----------------------GRayleigh_Ritz--------------------------
   SUBROUTINE GRayleigh_Ritz_real(veff,X,D)
      USE Lapack_module , ONLY : GeneralizeEigen_real,matmat_real
      IMPLICIT NONE
      !IN/OUT
      REAL(DP),INTENT(IN) :: veff(:,:,:)
      REAL(DP),INTENT(INOUT) :: X(:,:)
      REAL(DP),INTENT(OUT) :: D(:)
      !
      REAL(DP),DIMENSION(sn,sn) :: S_hat,H_hat,Q
      REAL(DP) :: Xnew(n,sn)
      INTEGER(I4B) :: i
      CHARACTER(len=8):: str_id
!xq
!INTEGER(I4B) :: t1,t2,t3,t4,t5,t6
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !Calculate the overlap matrix
!call system_clock(t1)
      CALL matmat_real(X,X,'T','N',S_hat)
!call system_clock(t2)
      !Calculate the project hamiltion
!call system_clock(t3)
      CALL Rayleigh_quotient_real(veff,sn,X,H_hat)
      ! print *,'output H'
      ! write(str_id,'(I8)')iiii
      ! open(1122,file="Hamiltonian"//trim(adjustl(str_id))//".txt")
      ! write(1122,*)'dimensions',sn, sn
      ! do i = 1,sn,1
      !    write(1122,'(2000(G16.9,X))')H_hat(:,i)
      ! enddo
      ! close(1122)
      ! iiii=iiii+1
!call system_clock(t4)
      !solve the generalized eigenvalue problem
!call system_clock(t5)
      CALL GeneralizeEigen_real(sn,H_hat,S_hat,Q,D)
!call system_clock(t6)
      !X=XQ
      CALL matmat_real(X,Q,'N','N',Xnew)
!print*,'Overlap',(t2-t1)/10000.d0
!print*,'Hhat',(t4-t3)/10000.d0
!print*,'Genera',(t6-t5)/10000.d0
!STOP
      X=Xnew
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE GRayleigh_Ritz_real
   !-----------------non-OrthNorm Chebyshev_filter ----------------
   SUBROUTINE cheby_filtering_GRRr(veff,X,D)
      !To avoid OrthNorm
      IMPLICIT NONE
      REAL(DP),INTENT(IN) :: veff(:,:,:)
      REAL(DP),INTENT(INOUT) :: X(:,:)
      REAL(DP),INTENT(INOUT) :: D(:)  !rayleigh-ritz value
      !LOCAL
      REAL(DP) :: a,b,al
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !call (al,a,b)
      a=MAXVAL(D)
      al=MINVAL(D)
      !up boundary
      CALL Estupb_real(7,veff,X(:,sn),b)
      !filtering (a,b)
      CALL chebyshev_filter_scaled_real(veff,X,CheM,a,b,al)
      !GRR (Rayleigh-Ritz Step)
      CALL GRayleigh_Ritz_real(veff,X,D)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE cheby_filtering_GRRr
   !##############################################################!
   !*For    :     BvK Cell Method Filtering (pRR)                 !
   !*Author : Qiang Xu                                            !
   !*Date   : 2018-03-06                                          !
   !##############################################################!
   !-----------------non-OrthNorm Chebyshev_filter ----------------
   SUBROUTINE cheby_filtering_pRRr(veff,X,Nfs,Cfr,Efr)
      !To avoid OrthNorm
      IMPLICIT NONE
      REAL(DP),INTENT(IN) :: veff(:,:,:)
      INTEGER(I4B),INTENT(IN) :: Nfs
      REAL(DP),INTENT(INOUT) :: X(:,:) &
                        &,      Efr(:) !pRR value
      REAL(DP),INTENT(OUT) :: Cfr(:,:) !pRR vector
      !LOCAL
      REAL(DP) :: a,b,al
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !call (al,a,b)
      a=MAXVAL(Efr)
      al=MINVAL(Efr)
      !up boundary
      CALL Estupb_real(7,veff,X(:,sn),b)
      !filtering (a,b)
      CALL chebyshev_filter_scaled_real(veff,X,CheM,a,b,al)
      !rRR (Rayleigh-Ritz Step)
      CALL PartialRayleighRitz(veff,X,Nfs,Cfr,Efr)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE cheby_filtering_pRRr
   !--------------------Partial Rayleigh Ritz----------------------
   SUBROUTINE PartialRayleighRitz(veff,X,Nfs,Cfr,Efr)
      USE parameters, ONLY : LRROrthNorm
      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B) :: &
                      & Nfs    !fraction states
      REAL(DP),INTENT(IN) :: veff(:,:,:)
      REAL(DP),INTENT(INOUT) :: X(:,:)
      REAL(DP),INTENT(OUT) :: Efr(:),Cfr(:,:)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF(LRROrthNorm)THEN
          CALL pRR_OrthNorm(veff,X,Nfs,Cfr,Efr)
      ELSE
          !CALL pRR_General(veff,X,D)
          print*,'waiting for pRR general'
          STOP
      ENDIF
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE PartialRayleighRitz
   !--------------------Partial Rayleigh Ritz----------------------
   SUBROUTINE pRR_OrthNorm(veff,X,Nfs,Cfr,Efr)
      USE Lapack_module , ONLY : OrthNorm_real
      USE Arpack_module , ONLY :rdiagM_arpk
      IMPLICIT NONE
      INTEGER(I4B) :: &
                      & Nfs    !fraction states
      REAL(DP),INTENT(IN) :: veff(:,:,:)
      REAL(DP),INTENT(INOUT) :: X(:,:)
      REAL(DP),INTENT(OUT) :: Efr(:),Cfr(:,:)
      !LOCAL
      REAL(DP),DIMENSION(sn,sn) :: Hhat,Q
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !Orth-Norm the subspace
      CALL OrthNorm_real(X)
      !xhx
      CALL Rayleigh_quotient_real(veff,sn,X,Hhat)
      !Partial diagM
      CALL rdiagM_arpk(Hhat,sn,Nfs,Cfr,Efr)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE pRR_OrthNorm
   !############################################################!
   !*For :Isolate cell self consistent                          !
   !*Author : qiangxu, modified by xlt                          !
   !*Date   : 2018-04-10                                        !
   !############################################################!
   !-----------------------DIVIDER-LINE--------------------------
   SUBROUTINE ISO_first_CheSCF_sto_rand(rhoS,nev,psi,eval)
      USE parameters , ONLY : Nspin,IDinit
      USE potential_module , ONLY : cal_veff_iso,ACCELERATE
      USE grid_module , ONLY : n,n1,n2,n3,nk,rho_calc,ISO_Vsphere
      USE struct_module , ONLY : naty,struct,natom





      USE m_time_evaluate, ONLY: memory_sum, memory_free
      IMPLICIT NONE
      !INOUT



      REAL(DP),INTENT(IN) :: rhoS(:,:)

      INTEGER(I4B),INTENT(IN) :: nev
      REAL(DP) :: psi(:,:,:)
      REAL(DP)  :: eval(:,:)
      !LOCAL
      INTEGER(I4B) :: Is,Ity &
               &,Nmax,Nrand,Ia




      REAL(DP) :: veff(parallel%mygrid_range(3),Nspin)
      REAL(DP) :: veff_local(parallel%mygrid_range(3))
      REAL(DP) :: psi_local(parallel%mygrid_range(3),nev)
      REAL(DP) :: eval_local(nev)

      REAL(DP),ALLOCATABLE :: initX(:,:)!,Shat(:,:),Hhat(:,:)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      call memory_sum('Iso_first_chescf_local',(real(size(veff),DP)+&
           & size(veff_local)+size(psi_local)+size(eval_local))*DP)




      ! print *,'bbbb'
      !Calculate the effective potential

# 1 "my_macro.h" 1
!output


!output time



!deallocate and allocate


!output file and line


!VERSION information
# 1440 "Chebyshev_fliter.f90" 2
      ACCELERATE=.FALSE.
      CALL cal_veff_iso(rhoS,veff)

# 1457 "Chebyshev_fliter.f90"

      !initialize the temp subspace
      !estimate number of the states
      IF(IDinit==0)THEN  !Slater orbital
      Nmax=0
      Nrand=0
      DO Ity=1,naty
            Nmax=Nmax+struct%nati(Ity)*(struct%Lmax(Ity)+1)**2
      ENDDO
      !test
      IF(Nmax<=nev)THEN
         Nrand=nev-Nmax
         Nmax=nev

         if(parallel%isroot)PRINT*,'[Add random wave for subspace]',Nrand



      ELSE
        Nmax=nev
      ENDIF
      !allocate



      ALLOCATE(initX(parallel%mygrid_range(3),Nmax))
      call memory_sum('Iso_first_chescf_initx_local',real(size(initX),DP)*DP)

      !init a 'super'-subspace
      ! print *,'storand1'
      CALL ISO_cheby_init_sto_rand(Nmax,Nrand,initX)
# 1499 "Chebyshev_fliter.f90"
      ! print *,'storand2'
      !> Raleigh-Ritz step
      !> all spin
      DO Is=1,Nspin




         Veff_local=Veff(:,Is)
         psi_local=psi(:,:,Is)
         eval_local=eval(:,Is)
         CALL ISO_first_subspace_sto_rand(Veff_local, &
              &   Nmax,nev,initX,psi_local,eval_local)
         Veff(:,Is)=Veff_local
         psi(:,:,Is)=psi_local
         eval(:,Is)=eval_local

      ENDDO
   ELSEIF(IDinit==1)THEN  !random initial orbital
      Nmax=nev
      Nrand=nev
      !allocate



      ALLOCATE(initX(parallel%mygrid_range(3),Nmax))

      call memory_sum('Iso_first_chescf_initx_local',real(size(initX),DP)*DP)
      CALL ISO_cheby_init_rand(Nmax,Nrand,initX)
      DO Is=1,Nspin




           Veff_local=Veff(:,Is)
           psi_local=psi(:,:,Is)
           eval_local=eval(:,Is)
           CALL ISO_first_subspace_sto_rand(Veff_local, &
                &   Nmax,nev,initX,psi_local,eval_local)
           Veff(:,Is)=Veff_local
           psi(:,:,Is)=psi_local
           eval(:,Is)=eval_local

        ENDDO
# 1561 "Chebyshev_fliter.f90"
      ENDIF
      IF(ALLOCATED(initX))THEN
         DEALLOCATE(initX)
         call memory_free('Iso_first_chescf_initx_local',real(size(initX),DP)*DP)
      ENDIF

      call memory_free('Iso_first_chescf_local',(real(size(veff),DP)+&
           & size(veff_local)+size(psi_local)+size(eval_local))*DP)




      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE ISO_first_CheSCF_sto_rand
   !-----------------------PARTING-LINE--------------------------
   !---------------------cheby_subspace----------------------------
   SUBROUTINE ISO_cheby_init_sto_rand(Nmax,Nrand,initX_sto)
      !initialize the subsystem by random STO
      USE parameters , ONLY : Nstates,Lpbc
      USE struct_module , ONLY : struct,lat_mat,naty
      USE grid_module , ONLY : grid,n,n1,n2,n3,rho_calc,gap
      USE math , ONLY : atom_sto
      USE Lapack_module , ONLY : OrthNorm_real

      USE smpi_math_module, ONLY: parallel

      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B),INTENT(IN) :: Nmax & !max dimen 'super'-subspace
                             &, Nrand  !number of random wave used
      REAL(DP),INTENT(OUT) :: initX_sto(:,:)
      !LOCAL
      INTEGER(I4B) :: Ity,Ia,Ip,icx,icy,icz,Ii &
            & ,nrep=1 &
            &,Il,Im,Nstart
      REAL(DP) :: ra(3),rat(4),f,randt
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      initX_sto(:,:)=0.d0
      !all states
      Ii=0
      !all type
     DO Ity=1,naty
        !all atom
        DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
           !store the position
           ra(:)=struct%poscar(:,Ia)
           !call the angle momentum
           DO Il=0,struct%Lmax(Ity)
              DO Im=-Il,Il
                 !call the states
                 Ii=Ii+1
                 IF(Ii.gt.Nmax)Ii=1
                 !all points



                 DO Ip=1,parallel%mygrid_range(3),1

                    !replicas
                      rat(1)=rho_calc%x(Ip)*gap(1)-ra(1)
                      rat(2)=rho_calc%y(Ip)*gap(2)-ra(2)
                      rat(3)=rho_calc%z(Ip)*gap(3)-ra(3)
                      rat(4)=SQRT( rat(1)**2 + rat(2)**2 +rat(3)**2 )
                      !STO
                      CALL atom_STO(struct%prinq(Il+1,Ity),Il,Im, &
                                  & struct%zeta(Il+1,Ity),rat,f)
                      !init
                      initX_sto(Ip,Ii)=initX_sto(Ip,Ii)+f
                 ENDDO
              ENDDO !> Im
           ENDDO !> Il
        ENDDO !> Ia
     ENDDO !> Ity

      !For plane-wave if needed
      IF(Nrand>0)THEN
        CALL random_seed()
        Nstart=Nmax-Nrand+1
        DO Ii=Nstart,Nmax



           DO Ip=1,parallel%mygrid_range(3),1

              CALL random_number(randt)
              initX_sto(Ip,Ii)=randt
           ENDDO
        ENDDO
        !orth-Norm if need
        ! CALL OrthNorm_real(initX_sto)
      ENDIF
      !orth-Norm if need
      !CALL OrthNorm(X)
      !> check
! #ifndef 1
!       open(111,file="initX_s")
!       write(111,*)initX_sto(:,1)
!       close(111)
! #else
!       if ( parallel%isroot ) then
!          open(111,file="initXi_p")
!          write(111,*)initX_sto(:,1)
!          close(111)
!       end if
! #endif
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE ISO_cheby_init_sto_rand
   !-----------------------DIVIDER-LINE--------------------------
   SUBROUTINE ISO_cheby_init_rand(Nmax,Nrand,initX_sto)
      !initialize the subsystem by random STO
      USE parameters , ONLY : Nstates,Lpbc
      USE struct_module , ONLY : struct,lat_mat,naty
      USE grid_module , ONLY : grid,n,n1,n2,n3,rho_calc,gap
      USE math , ONLY : atom_sto
      USE Lapack_module , ONLY : OrthNorm_real

      USE smpi_math_module, ONLY: parallel

      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B),INTENT(IN) :: Nmax & !max dimen 'super'-subspace
                             &, Nrand  !number of random wave used
      REAL(DP),INTENT(OUT) :: initX_sto(:,:)
      !LOCAL
      INTEGER(I4B) :: Ity,Ia,Ip,icx,icy,icz,Ii &
            & ,nrep=1 &
            &,Il,Im,Nstart
      REAL(DP) :: ra(3),rat(4),f,randt
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      initX_sto(:,:)=0.d0
      !all states
      Ii=0

      !For plane-wave if needed
      IF(Nrand>0)THEN
        CALL random_seed()
        Nstart=Nmax-Nrand+1
        DO Ii=Nstart,Nmax



           DO Ip=1,parallel%mygrid_range(3),1

              CALL random_number(randt)
              initX_sto(Ip,Ii)=randt !cos(real(Ii*Ip))
           ENDDO
        ENDDO
        !orth-Norm if need
        ! CALL OrthNorm_real(initX_sto)
      ENDIF
      !orth-Norm if need
      !CALL OrthNorm(X)
      !> check
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE ISO_cheby_init_rand
   !-----------------------DIVIDER-LINE--------------------------
   SUBROUTINE ISO_first_subspace_sto_rand(veff_3D,Nmax,nev,initX,X,D)
      USE Lapack_module , ONLY : matmat_real &
                       &, GeneralizeEigen_real

      USE parameters    , ONLY : BLOCK_MBNB
      USE ScaLapack_module, ONLY:SL_GeneralizeEigen_real,SL_matmat_real

      USE m_time_evaluate, ONLY: memory_sum, memory_free
      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B),INTENT(IN) :: Nmax,nev



      REAL(DP),INTENT(IN) :: veff_3D(:)

      REAL(DP),INTENT(IN)  :: initX(:,:)
      REAL(DP),INTENT(OUT) :: X(:,:) !subspace
      REAL(DP),INTENT(OUT) :: D(:)      !eigenvalue
      !LOCAL
      REAL(DP),DIMENSION(Nmax,Nmax) :: Shat,Hhat,Qs
      REAL(DP) :: Xnew(rho_calc%OneDLength,Nmax)
      REAL(DP) :: Ds(Nmax)

      REAL(DP),ALLOCATABLE :: Shat_local(:,:), Hhat_local(:,:), qs_local(:,:)
      INTEGER(I4B)         :: m,n,mb,nb,cm,cn,cmb,cnb
      REAL(DP),ALLOCATABLE :: rshp(:,:),QT(:,:),XT(:,:)!,rshp3(:,:),rShat_local(:,:),rHhat_local(:,:)

      INTEGER(I4B)         :: i
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      !> for matmat in parallel
      m=rho_calc%OneDLength  !>> X dimension 1
      n=nev                    !>> X dimension 2
      mb=BLOCK_MBNB !>> X splited dimension 1
      nb=nev !parallel%nstate_proc  !>> X splited dimension 2
      cm=nev!parallel%nstate_proc !>> { X^T [* H]* X } dimension 1
      cn=nev !>> { X^T [* H]* X } dimension 2
      cmb=BLOCK_MBNB!nev !>> { X^T [* H]* X } splited dimension 1
      cnb=BLOCK_MBNB!parallel%nstate_proc !>> { X^T [* H]* X } splited dimension 2
      !> allocate the array of parallel calculation
      IF(allocated(Shat_local))THEN
         call memory_free('iso_first_subspace_shat_local',real(size(Shat_local),DP)*DP)
         deallocate(Shat_local)
      ENDIF
      IF(allocated(Hhat_local))THEN
         call memory_free('iso_first_subspace_hhat_local',real(size(Hhat_local),DP)*DP)
         deallocate(Hhat_local)
      ENDIF
      IF(allocated(qs_local))THEN
         call memory_free('iso_first_subspace_qs_local',real(size(qs_local),DP)*DP)
         deallocate(qs_local)
      ENDIF
      allocate(Shat_local(parallel%nstate_proc,nev))
      allocate(Hhat_local(parallel%nstate_proc,nev))
      allocate(qs_local(parallel%nstate_proc,nev))
      call memory_sum('iso_first_subspace_local',real(size(Shat_local),DP)*DP&
           &+size(Hhat_local)*DP+size(qs_local)*DP)

      !Raleigh-Ritz step
      !Overlap matrix
# 1787 "Chebyshev_fliter.f90"
      ! print *,"nstate",parallel%nstate_proc
      ! print*,"initX",shape(initX)
      ! print*,"shat_local",shape(Shat_local)
      ! print *,"m,n,mb,nb,cmb,cnb",m,n,mb,nb,cmb,cnb
      ! if(parallel%isroot)then
      ! open(111,file="initx_p")
      ! write(111,*)initX(:,1)
      ! close(111)
      ! endif
      CALL SL_matmat_real('t','n',initx,initx,Shat_local,m,n,m,n,mb,nb,mb,nb,cmb,cnb)
      ! if(parallel%isroot)then
      ! print *,"shape Shat",shape(Shat_local)
      ! open(111,file="Shat_p")
      !   ! do i=1,size(Shat_local,2),1
      !   !    write(111,'(1000(F6.2,2X))')Shat_local(:,i)
      !   ! enddo
      ! write(111,*)shat_local
      ! close(111)
      ! endif

      ! stop
      !projected hamiltonian



      CALL Rayleigh_quotient_ISO(veff_3D,nev,initX,Hhat_local)

      !eigen-decomposion
# 1833 "Chebyshev_fliter.f90"
      CALL SL_GeneralizeEigen_real(nmax,hhat_local,shat_local,nmax,nmax,nmax,nmax, &
           &cmb,cnb,cmb,cnb,qs_local,nmax,nmax,ds,cmb,cnb)
      ! if(parallel%isroot)then
      ! print *,"shape Qs",shape(qs_local)
      ! open(111,file="Hhat_p")
      !   do i=1,size(hhat_local,2),1
      !      write(111,'(1000(F8.4,2X))')hhat_local(:,i)
      !   enddo
      ! close(111)
      ! open(111,file="Shat_p")
      !   do i=1,size(qs_local,2),1
      !      write(111,'(1000(F8.4,2X))')shat_local(:,i)
      !   enddo
      ! close(111)
      ! open(111,file="Qs_p")
      !   do i=1,size(qs_local,2),1
      !      write(111,'(1000(F8.4,2X))')qs_local(:,i)
      !   enddo
      ! close(111)
      ! endif

      !rotation



      if(allocated(shat_local))then
         call memory_free('iso_first_subspace_shat_local',real(size(Shat_local),DP)*DP)
         deallocate(shat_local)
      endif
      if(allocated(hhat_local))then
         call memory_free('iso_first_subspace_hhat_local',real(size(Hhat_local),DP)*DP)
         deallocate(hhat_local)
      endif
      CALL SL_matmat_real('n','n',initx,qs_local,X,m,n,nmax,nmax,mb,nb,cmb,cmb,cmb,cnb)
      !> reshape
      ! allocate(rshp(nev,parallel%nstate_proc))
      ! allocate(QT(nev,parallel%nstate_proc))
      ! allocate(XT(m,parallel%nstate_proc))
      ! rshp=0
      ! do i=1,parallel%nstate_proc,1
      !    rshp(parallel%sub2sum(i,parallel%myid+1),i)=1
      ! enddo
      ! !> X * Q * P
      ! CALL SL_matmat_real('n','n',initx,qs_local,XT,m,n,nmax,nmax,mb,nb,cmb,cmb,mb,cnb)
      ! CALL SL_matmat_real('n','n',XT,rshp,x,m,n,nmax,nmax,mb,nb,cmb,cmb,mb,cnb)
      ! !> release memory
      if(allocated(qs_local))then
         call memory_free('iso_first_subspace_hhat_local',real(size(qs_local),DP)*DP)
         deallocate(qs_local)
      endif
      ! if(allocated(rshp)) deallocate(rshp)
      ! if(allocated(QT)) deallocate(QT)
      ! if(allocated(XT)) deallocate(XT)

      !eigen-value
      D(1:nev)=Ds(1:nev)



      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE ISO_first_subspace_sto_rand
   !-----------------------PARTING-LINE--------------------------
   !-----------------------DIVIDER-LINE--------------------------
   SUBROUTINE Rayleigh_quotient_ISO(veff_3D,nst,x,xhx)
      !xhx=(X,H,X)
      USE Lapack_module , ONLY : matmat_real

      USE ScaLapack_module , ONLY : SL_matmat_real
      USE parameters    , ONLY : BLOCK_MBNB

      USE m_time_evaluate, ONLY: memory_sum,memory_free
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN) :: nst



      REAL(DP),INTENT(IN) :: veff_3D(:)

      REAL(DP) :: x(:,:),xhx(:,:)
      !LOCAL



      REAL(DP) :: hx(parallel%mygrid_range(3),nst)

      INTEGER(I4B) :: i,j
      character(len=4) :: str_id,str_id2
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      call memory_sum('first_subspace_Rayleigh_local',real(size(hx),DP)*DP)
      CALL cal_HX_ISO(veff_3D,nst,x,hx)
# 1941 "Chebyshev_fliter.f90"
      !xhx(:,:)=MATMUL( TRANSPOSE(CONJG(x)) , hx )



     CALL SL_matmat_real('t','n',x,hx,xhx,rho_calc%OneDLength,nst,&
          &rho_calc%OneDLength,nst,&
          &BlOCK_MBNB,nst,&
          &BLOCK_MBNB,nst,&
          &BLOCK_MBNB,BLOCK_MBNB)


      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      call memory_free('first_subspace_Rayleigh_local',real(size(hx),DP)*DP)
   ENDSUBROUTINE Rayleigh_quotient_ISO
   !-----------------------PARTING-LINE--------------------------
   !-----------------------DIVIDER-LINE--------------------------
   SUBROUTINE cal_HX_ISO(veff,nst,V,HV)
      USE matvec_module , ONLY : ISO_rmatvec
      USE m_time_evaluate, ONLY: memory_sum,memory_free
      IMPLICIT NONE
      !INOUT
      INTEGER(I4B),INTENT(IN) :: nst



      REAL(DP),INTENT(IN)  :: veff(:) !effective potential

      REAL(DP),INTENT(IN)  :: V(:,:)
      REAL(DP),INTENT(OUT) :: HV(:,:)
      !LOCAL
      INTEGER(I4B) :: Is,tq,tr
      LOGICAL :: l1=.false.




      REAL(DP) :: V_local(parallel%mygrid_range(3))
      REAL(DP) :: HV_local(parallel%mygrid_range(3))
      character(len=4) :: str_id,str_id2

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      call memory_sum('cal_HX_local',(real(size(V_local),DP)+size(HV_local))*DP)
      DO Is=1,nst,1

         CALL start_time('sub_calHX_rmatvec',l1)

         V_local=V(:,Is)
         HV_local=HV(:,Is)



         CALL ISO_rmatvec(veff,V_local,HV_local,parallel%mygrid_range(3))

         HV(:,Is)=HV_local

         CALL end_time('sub_calHX_rmatvec',l1)
      ! CALL MPI_BARRIER(parallel%comm,mpinfo)

      ENDDO
      !=========================================
      call memory_free('cal_HX_local',(real(size(V_local),DP)+size(HV_local))*DP)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE cal_HX_ISO
   !-----------------------PARTING-LINE--------------------------
   !-----------------------DIVIDER-LINE--------------------------
   !-----------------non-OrthNorm Chebyshev_filter ----------------
   SUBROUTINE cheby_filtering_GRRiso(veff_3D,X_sphe,D)
      !To avoid OrthNorm
     ! USE array_io
      IMPLICIT NONE



      REAL(DP),INTENT(IN) :: veff_3D(:)

      REAL(DP),INTENT(INOUT) :: X_sphe(:,:)
      REAL(DP),INTENT(INOUT) :: D(:)  !rayleigh-ritz value
      !LOCAL
      REAL(DP) :: a,b,al

      INTEGER(I4B) :: bcastid,i,j
      LOGICAL :: l1=.true.

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !call (al,a,b)
      a=MAXVAL(D)
      al=MINVAL(D)
      !> up boundary
      CALL Estupb_ISO(7,veff_3D,X_sphe(:,sn),b)
      !> filtering (a,b)

      CALL start_time("chebyshev_filter_scaled_ISO",l1); CALL chebyshev_filter_scaled_ISO(veff_3D,X_sphe,CheM,a,b,al) ; CALL end_time("chebyshev_filter_scaled_ISO",l1)



      !> GRR (Rayleigh-Ritz Step)

      CALL start_time("GRayleigh_Ritz_ISO",l1); CALL GRayleigh_Ritz_ISO(veff_3D,X_sphe,D) ; CALL end_time("GRayleigh_Ritz_ISO",l1)



      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE cheby_filtering_GRRiso
   !-----------------------PARTING-LINE--------------------------
   !-----------------------DIVIDER-LINE--------------------------
   !-------------------------upper bound---------------------------
   SUBROUTINE Estupb_ISO(k,veff_3D,vec,b)
      !
      USE Lapack_module , ONLY : diagM_real
      USE matvec_module , ONLY : ISO_rmatvec
      USE array_io

      USE smpi_math_module , ONLY : parallel,mpi_real8,mpi_sum,mpinfo

      USE m_time_evaluate, ONLY: memory_sum,memory_free
      IMPLICIT NONE
      !INOUT
      INTEGER(I4B),INTENT(IN) :: k
      !hamiltonian



      REAL(DP),INTENT(IN) :: veff_3D(:)

      REAL(DP),INTENT(IN) :: vec(:)
      REAL(DP),INTENT(OUT) :: b
      !LOCAL




      REAL(DP) :: v0(parallel%mygrid_range(3)),v(parallel%mygrid_range(3)),f(parallel%mygrid_range(3))
      REAL(DP) :: alpha_local,alpha,beta_local,beta,eval(k),mz,b_local

      REAL(DP) :: T(k,k),evec(k,k)
      INTEGER(I4B) :: Nj
      INTEGER(I4B) :: J
      REAL(DP),external :: ddot
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      call memory_sum('Estub_iso',(real(size(v0),DP)+size(v)+size(f))*DP+2*k**2*DP)
      v=vec !should normalize
      T(:,:)=0.d0
      Nj=MIN(k,10)
      !



      CALL ISO_rmatvec(veff_3D,v,f,parallel%mygrid_range(3))

      ! CALL output(size(f),f,'f_estupb_serial')
      ! CALL output(size(veff_3D),reshape(f,(/size(veff_3D)/)),'v_estupb_serial')
      ! alpha=DOT_PRODUCT(f,v)



      alpha_local=ddot(size(f),f,1,v,1)
      call mpi_allreduce(alpha_local, alpha, 1, mpi_real8, mpi_sum, parallel%comm, mpinfo)

      f=f-alpha*v
      T(1,1)=alpha
      DO J=2,Nj
         ! beta=SQRT(REAL(DOT_PRODUCT(f,f),8))



         beta_local=ddot(size(f),f,1,f,1)
         call mpi_allreduce(beta_local, beta, 1, mpi_real8, mpi_sum, parallel%comm, mpinfo)

         beta=sqrt(beta)
         v0=v
         v=f/beta



         CALL ISO_rmatvec(veff_3D,v,f,parallel%mygrid_range(3))

         f=f-beta*v0
         !>>
         ! alpha=DOT_PRODUCT(f,v)



         alpha_local=ddot(size(f),f,1,v,1)
         call mpi_allreduce(alpha_local, alpha, 1, mpi_real8, mpi_sum, parallel%comm, mpinfo)

         f=f-alpha*v
         T(J,J-1)=beta
         T(J-1,J)=beta
         T(J,J)=alpha
      ENDDO
      !NORM2(T)
      !b=Norm_2(T,k) + SQRT(REAL(DOT_PRODUCT(f,f),8))
      CALL diagM_real(T,evec,eval)
      !mz=MAX( ABS(evec(k,k)), ABS(evec(k,k-1)) , ABS(evec(k,k-2)) )
      mz=ABS(evec(k,k))
      ! b=eval(k) + SQRT(REAL(DOT_PRODUCT(f,f),8))*mz



      b_local=ddot(size(f),f,1,f,1)
      call mpi_allreduce(b_local, b, 1, mpi_real8, mpi_sum, parallel%comm, mpinfo)

      b=eval(k) + SQRT(REAL(b,DP))*mz
      call memory_free('Estub_iso',(real(size(v0),DP)+size(v)+size(f))*DP+2*k**2*DP)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE Estupb_ISO
   !-----------------------PARTING-LINE--------------------------
   !-----------------------DIVIDER-LINE--------------------------
   !-------------------chebyshev_filter_scaled---------------------
   SUBROUTINE chebyshev_filter_scaled_ISO(veff_3D,X,m,a,b,al)
      !
     USE m_time_evaluate, ONLY: memory_sum,memory_free
      IMPLICIT NONE
      !INOUT



      REAL(DP),INTENT(IN) :: veff_3D(:) !veff

      REAL(DP),INTENT(INOUT) :: X(:,:)
      INTEGER(I4B),INTENT(IN)   :: m ! the m degree Chebyshev polynomial we used
      REAL(DP),INTENT(IN) :: a,b,al  !interval [a,b] to be filter
      !LOCAL
      REAL(DP) :: e &   !(b-a)/2
            & ,   c &   !(b+a)/2
            & , sigma,sigmanew,tau
      !



      REAL(DP),DIMENSION(parallel%mygrid_range(3),sn)  :: HV,Y,Ynew

      INTEGER(I4B) :: Ic,tl,tm,tn,to,tp
      REAL(DP) :: temp1,temp2
      INTEGER(I4B) :: array_size
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      call memory_sum('filter_iso_local',real(size(HV),DP)*3*DP)
      IF(m<2)THEN
         WRITE(6,*) 'Chebyshev_filter: the degree m must larger than 1'
      ENDIF

      e = ( b - a ) / 2
      c = ( b + a ) / 2
      sigma=e / (c-al)
      tau=2.d0 / sigma



      CALL cal_HX_ISO(veff_3D,sn,X,HV)

      temp1=sigma / e
      Y= ( HV - c*X ) * temp1
      DO Ic=2,m
         sigmanew=1.d0 / (tau-sigma)



         CALL cal_HX_ISO(veff_3D,sn,Y,HV)

         temp1=2.d0*sigmanew/e
         temp2=sigma*sigmanew
         Ynew=( HV - c*Y )*temp1 - temp2*X
         X=Y
         Y=Ynew
         ! > temp2=sigma*sigmanew
         ! temp2=-sigma*sigmanew
         !> Ynew=( HV - c*Y )*temp1 - temp2*X
         ! array_size=size(HV)
         ! CALL dscal(array_size,temp1,HV,1)
         ! temp1=-temp1*c
         ! CALL daxpy(array_size,temp1,Y,1,HV,1)
         ! CALL daxpy(array_size,temp2,X,1,HV,1)
         ! ! X=Y
         ! CALL dcopy(array_size,Y,1,X,1)
         ! ! Y=Ynew
         ! CALL dcopy(array_size,HV,1,Y,1)
         sigma=sigmanew
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      call memory_free('filter_iso_local',real(size(HV),DP)*3*DP)
   ENDSUBROUTINE chebyshev_filter_scaled_ISO
   !-----------------------PARTING-LINE--------------------------
   !-----------------------DIVIDER-LINE--------------------------
   SUBROUTINE GRayleigh_Ritz_ISO(veff_3D,X,D)
      USE Lapack_module , ONLY : GeneralizeEigen_real,matmat_real

      USE parameters    , ONLY : BLOCK_MBNB
      USE ScaLapack_module

      USE m_time_evaluate, ONLY: memory_sum,memory_free
      IMPLICIT NONE
      !IN/OUT



      REAL(DP),INTENT(IN) :: veff_3D(:)

      REAL(DP),INTENT(INOUT) :: X(:,:)
      REAL(DP),INTENT(OUT) :: D(:)
      !




      REAL(DP) :: Xnew(parallel%mygrid_range(3),sn)
      REAL(DP),ALLOCATABLE :: Shat_local(:,:), Hhat_local(:,:), Q(:,:)
      INTEGER(I4B)         :: m,n,mb,nb,cm,cn,cmb,cnb
      REAL(DP),ALLOCATABLE :: rshp(:,:),QPT(:,:)
      INTEGER(I4B)         :: i

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



      call memory_sum('Rayleigh_ritz_local',real(size(Xnew),DP)*DP)


      !> for matmat in parallel
      m=rho_calc%OneDLength  !>> X dimension 1
      n=sn                    !>> X dimension 2
      mb=BLOCK_MBNB !>> X splited dimension 1
      nb=sn !parallel%nstate_proc  !>> X splited dimension 2
      cm=sn!parallel%nstate_proc !>> { X^T [* H]* X } dimension 1
      cn=sn !>> { X^T [* H]* X } dimension 2
      cmb=BLOCK_MBNB!nev !>> { X^T [* H]* X } splited dimension 1
      cnb=BLOCK_MBNB!parallel%nstate_proc !>> { X^T [* H]* X } splited dimension 2
      !> allocate the array of parallel calculation
      IF(allocated(Shat_local))THEN
         call memory_free('Rayleigh_ritz_local',real(size(Shat_local),DP)*DP)
         deallocate(Shat_local)
      ENDIF
      IF(allocated(Hhat_local))THEN
         call memory_free('Rayleigh_ritz_local',real(size(Hhat_local),DP)*DP)
         deallocate(Hhat_local)
      ENDIF
      IF(allocated(Q))THEN
         call memory_free('Rayleigh_ritz_local',real(size(Q),DP)*DP)
         deallocate(Q)
      ENDIF
      allocate(Shat_local(parallel%nstate_proc,sn))
      allocate(Hhat_local(parallel%nstate_proc,sn))
      allocate(Q(parallel%nstate_proc,sn))
      call memory_sum('Rayleigh_ritz_local',(real(size(Shat_local),DP)+size(Hhat_local)+size(Q))*DP)

      !> Raleigh-Ritz step
      !Calculate the overlap matrix



      CALL SL_matmat_real('t','n',X,X,Shat_local,m,n,m,n,mb,nb,mb,nb,cmb,cnb)

      !Calculate the project hamiltion



      CALL Rayleigh_quotient_ISO(veff_3D,sn,X,Hhat_local)

      !solve the generalized eigenvalue problem



      CALL SL_GeneralizeEigen_real(sn,Hhat_local,Shat_local,sn,sn,sn,sn, &
           &cmb,cnb,cmb,cnb,Q,sn,sn,D,cmb,cnb)

      !X=XQ



      if(allocated(shat_local))THEN
         call memory_free('Rayleigh_ritz_local',real(size(Shat_local),DP)*DP)
         deallocate(shat_local)
      ENDIF
      if(allocated(hhat_local))THEN
         call memory_free('Rayleigh_ritz_local',real(size(hhat_local),DP)*DP)
         deallocate(hhat_local)
      ENDIF
      CALL SL_matmat_real('n','n',X,Q,Xnew,m,n,sn,sn,mb,nb,cmb,cmb,mb,cnb)
      ! allocate(rshp(sn,parallel%nstate_proc))
      ! allocate(QPT(sn,parallel%nstate_proc))
      ! rshp=0
      ! do i=1,parallel%nstate_proc,1
      !    rshp(parallel%sub2sum(i,parallel%myid+1),i)=1
      ! enddo
      ! !> X * Q * P
      ! CALL SL_matmat_real('n','n',Q,rshp,QPT,n,n,sn,sn,mb,nb,cmb,cmb,mb,cnb)
      ! CALL SL_matmat_real('n','n',X,QPT,Xnew,m,n,sn,sn,mb,nb,cmb,cmb,mb,cnb)
      ! if(allocated(Q)) deallocate(Q)
      ! if(allocated(rshp)) deallocate(rshp)
      ! if(allocated(QPT)) deallocate(QPT)

      X=Xnew

      IF(allocated(Q))THEN
         call memory_free('Rayleigh_ritz_local',real(size(Q),DP)*DP)
         deallocate(Q)
      ENDIF

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



      call memory_free('Rayleigh_ritz_local',real(size(Xnew),DP)*DP)

   ENDSUBROUTINE GRayleigh_Ritz_ISO
   !-----------------------PARTING-LINE--------------------------
   !---------------------cheby_subspace----------------------------
   SUBROUTINE first_CheSCF_sto_gamma(veff,nev,psi,eval)
      USE parameters , ONLY : Nspin
      USE potential_module , ONLY : cal_veff
      USE grid_module , ONLY : n,n1,n2,n3,nk
      USE struct_module , ONLY : naty,struct
      use m_time_evaluate ,only: memory_sum,memory_free

      USE smpi_math_module, ONLY: states_split,parallel,mpinfo
      ! USE ScaLapack_module, ONLY: twoD_map

      IMPLICIT NONE
      !INOUT
      REAL(DP),INTENT(IN) :: veff(:,:,:,:)
      INTEGER(I4B),INTENT(IN) :: nev
      REAL(DP) :: psi(:,:,:,:)
      REAL(DP)  :: eval(:,:,:)
      !LOCAL
      INTEGER(I4B) :: Is,Ik,Ity &
               &,Nmax,Npw
      REAL(DP),ALLOCATABLE :: initX(:,:)!,Shat(:,:),Hhat(:,:)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !initialize the temp subspace
      !estimate number of the states
      Nmax=0
      Npw=0
      DO Ity=1,naty
         Nmax=Nmax+struct%nati(Ity)*(struct%Lmax(Ity)+1)**2
      ENDDO

      call states_split(Nmax)
      ! NMAX=twoD_map(2,parallel%rankx,parallel%ranky)

      !test
      IF(Nmax<=nev)THEN
         Npw=nev-Nmax
         Nmax=nev
         !PRINT*,'[Add plane wave for subspace]',Npw
      ELSE
         Nmax=nev
      ENDIF
      !allocate
      ALLOCATE(initX(n,nev))
      call memory_sum("first_Chebyshev",real(size(initX),DP)*DP)
      !init a 'super'-subspace
      CALL cheby_init_sto_gamma(nev,Npw,initX)
      ! open(1111,file='initX',status='old')
      ! read(1111,*)initX
      ! close(1111)
      !Raleigh-Ritz step
      !all spin
      DO Is=1,Nspin
         !all kpoints
         DO Ik=1,nk
            CALL first_subspace_stopw_gamma(Ik,veff(:,:,:,Is), &
               &   Nmax,nev,initX,psi(:,:,Ik,Is),eval(:,Ik,Is))
            !filter
        !    IF(Npw>0)THEN
        !       CALL first_SCFstep_filter(Ik,veff(:,:,:,Is), &
        !       &   psi(:,:,Ik,Is),eval(:,Ik,Is))
        !    ENDIF
         ENDDO
      ENDDO
      call memory_free("first_Chebyshev",real(size(initX),DP)*DP)
      DEALLOCATE(initX)
      call memory_free("first_Chebyshev_local",real(size(veff),DP)*DP)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ENDSUBROUTINE first_CheSCF_sto_gamma
   !------------initialize_subspace by random STO-------------
   SUBROUTINE cheby_init_sto_gamma(Nmax,Npw,initX_sto)
      !initialize the subsystem by random STO
      USE parameters , ONLY : Nstates,Nstates_global
      USE struct_module , ONLY : struct,lat_mat,naty
      USE grid_module , ONLY : grid,n,n1,n2,n3
      USE math , ONLY : atom_sto
      USE Lapack_module , ONLY : OrthNorm_real
      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B),INTENT(IN) :: Nmax & !max dimen 'super'-subspace
                             &, Npw  !number of plane wave used
      REAL(DP),INTENT(OUT) :: initX_sto(:,:)
      !LOCAL
      INTEGER(I4B) :: Ity,Ia,Ip,icx,icy,icz,Ii &
            & ,nrep=1 &



            &,Il,Im,Nstart &
            &,Is,nassign

      REAL(DP) :: ra(3),rat(4),f,gdotr
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      initX_sto(:,:)=cmplx(0.d0,0.d0)
      !all states
      Ii=0

      !> init the count independently in per process
      nassign=0

      !all type
      DO Ity=1,naty
         !all atom
         DO Ia=struct%eleid(Ity),struct%eleid(Ity+1)-1
            !store the position
            ra(:)=struct%poscar(:,Ia)
            !call the angle momentum
            DO Il=0,struct%Lmax(Ity)
               DO Im=-Il,Il
                  !call the states
                  Ii=Ii+1
                  if(Ii.gt.Nstates_global)Ii=1

                  do Is=1,parallel%nstate_proc,1
                  If(Ii==parallel%sub2sum(Is,parallel%ranky+1))THEN
                     nassign=nassign+1
                     IF(nassign.gt.parallel%nstate_proc)nassign=1

                  !all points
                  DO Ip=1,n
                     !replicas
                     DO icz=-nrep,nrep
                     DO icy=-nrep,nrep
                     DO icx=-nrep,nrep
                        rat(1:3)= icx*lat_mat(:,1)  &
                          & +  icy*lat_mat(:,2)  &
                          & +  icz*lat_mat(:,3)  &
                          & +  grid%rVec(1:3,Ip)  &
                          & -  ra(:)
                        rat(4)=SQRT( rat(1)**2 + rat(2)**2 +rat(3)**2 )
                        !STO
                    CALL atom_STO(struct%prinq(Il+1,Ity),Il,Im,struct%zeta(Il+1,Ity),rat,f)
                    !init



                    initX_sto(Ip,nassign)=initX_sto(Ip,nassign)+CMPLX(f,0.d0)

                     ENDDO
                     ENDDO
                     ENDDO
                     !
                  ENDDO

                  ENDIF
                  ENDDO !> Is

                  !
               ENDDO
            ENDDO
            !
         ENDDO
         !
      ENDDO
! OPEN(file_unit,FILE=filename(8:11)//'test.dat')
! WRITE(file_unit,*) '14 14 14'
! WRITE(file_unit,*) initX_sto
! CLOSE(file_unit)
! CALL MPI_finalize(mpinfo)
! stop
      !For plane-wave if needed
      IF(Npw>0)THEN
        Nstart=Nmax-Npw+1
        DO Ii=Nstart,Nmax
           DO Ip=1,n
              gdotr=DOT_PRODUCT(grid%gVec(1:3,Ii),grid%rVec(1:3,Ip))
              initX_sto(Ip,Ii)=EXP(IMAG*gdotr)
           ENDDO
        ENDDO
        !orth-Norm if need
        CALL OrthNorm_real(initX_sto)
      ENDIF
      !orth-Norm if need
      !CALL OrthNorm(X)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ENDSUBROUTINE cheby_init_sto_gamma
   !----------------------------------------------------------
   SUBROUTINE first_subspace_stopw_gamma(Ik,veff,Nmax,nev,initX,X,D)
      USE Lapack_module , ONLY : diagM,OrthNorm,matmat &
                       &, GeneralizeEigen
      ! use m_time_evaluate, only:memory_sum,memory_free
      use m_time_evaluate, only:memory_sum,memory_free,filename

      use parameters, only: Nstates_global, BLOCK_MBNB
      use ScaLapack_module, only: SL_GeneralizeEigen_real2 &
           & ,SL_matmat_realtn, SL_matmat_realnn,twoD_map
           ! &, SL_matmat,twoD_map
      USE Grid_module, ONLY: global_n

      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B),INTENT(IN) :: Ik &
                             &, Nmax,nev
      REAL(DP),INTENT(IN) :: veff(:,:,:)
      REAL(DP),INTENT(INOUT)  :: initX(:,:)
      REAL(DP),INTENT(OUT) :: X(:,:) !subspace
      REAL(DP),INTENT(OUT) :: D(:)      !eigenvalue
      !LOCAL



      REAL(DP),DIMENSION(twoD_map(1,parallel%rankx,&
           &parallel%ranky),Nmax) :: Shat,Hhat,Qs

      REAL(DP) :: Xnew(n,Nmax)
      REAL(DP) :: Ds(Nmax)

      INTEGER(I4B) :: mb,nb,cmb,cnb
      INTEGER(I4B) :: ix,iy
      ! COMPLEX(DCP) :: initx_global(global_n)
      REAL(DP) :: veff_global(global_n)

      integer(I4B) :: i
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !Raleigh-Ritz step
      !Overlap matrix
      call memory_sum("first_subspace_sto_local",(real(size(Shat),DP)*3+size(Xnew)+size(Ds))*DP)



      mb=1 !n1*n2*n3
      nb=1!BLOCK_MBNB
      cmb=1!BLOCK_MBNB
      cnb=1!BLOCK_MBNB
      ! CALL SL_matmat_real('c','n',initX,initX,Shat,&
      CALL SL_matmat_realtn('t','n',initX,initX,Shat,&
           &global_n,Nstates_global,global_n,Nstates_global,mb,nb,mb,nb,cmb,cnb)
      ! ix=mod(mod(parallel%mygrid_range(1)-1,n1*n2),n1)+1
      ! iy=mod((parallel%mygrid_range(1)-1)/n1,n2)+1
      ! call MPI_ALLGATHERV(veff(ix,iy,1),parallel%mygrid_range(3)&
      !      &,MPI_REAL8,veff_global,parallel%recvcounts&
      !      &,parallel%displs,MPI_REAL8&
      !      & ,parallel%commx,mpinfo)
      ! ! ! print *,'shape initX2',shape(initX)
      ! open(1111+parallel%myid,file="veff_cheb"//filename(10:11))
      ! write(1111+parallel%myid,*)veff_global
      ! do i=1,size(initX,2),1
      !    write(1111+parallel%myid,'(10000(F25.16,1X))')real(initX(:,i))
      ! enddo
      ! close(1111+parallel%myid)
      ! open(1111+parallel%myid,file="Shat_p"//filename(10:11))
      ! do i=1,size(Shat,2),1
      !    write(1111+parallel%myid,'(1000(F25.16,2X))')real(Shat(:,i))
      ! enddo
      ! close(1111+parallel%myid)
      ! call smpi_exit()
      ! CALL MPI_finalize(mpinfo)
      ! stop

      !projected hamiltonian
      CALL Rayleigh_quotient_gamma(Ik,veff,Nmax,initX,Hhat)
      !eigen-decomposion



      CALL SL_GeneralizeEigen_real2(Nstates_global,Hhat,Shat,Nstates_global&
           &,Nstates_global,Nstates_global,Nstates_global&
           &,cmb,cnb,cmb,cnb,Qs,Nstates_global,Nstates_global&
           &,D,cmb,cnb)

! #ifndef 1
!        open(1111,file="Hhat")
!        do i=1,size(Hhat,2),1
!           write(1111,'(10000(F25.15,1X))')real(Hhat(:,i))
!        enddo
!           write(1111,*)D
!        close(1111)
!        stop
! #else
!      open(1111+parallel%myid,file="Hhat"//filename(10:11))
!      do i=1,size(Hhat,2),1
!         write(1111+parallel%myid,'(10000(F25.15,1X))')real(Hhat(:,i))
!      enddo
!          write(1111+parallel%myid,*)D
!      close(1111+parallel%myid)
!      call smpi_exit()
! #endif
      !rotation



      ! CALL SL_matmat('n','n',initX,Qs,X,global_n,Nstates_global&
      CALL SL_matmat_realnn('n','n',initX,Qs,X,global_n,Nstates_global&
           &,Nstates_global,Nstates_global,mb,nb,cmb,cnb,cmb,cnb)
     ! open(1111+parallel%myid,file="Hhat"//filename(10:11))
     ! do i=1,size(X,2),1
     !    write(1111+parallel%myid,'(10000(F25.15,1X))')real(X(:,i))
     ! enddo
     !     write(1111+parallel%myid,*)D
     ! close(1111+parallel%myid)
     ! call smpi_exit()

      !eigen-value
      !D(1:nev)=Ds(1:nev)
      !X(:,1:nev)=Xnew(:,1:nev)
      call memory_free("first_subspace_sto_local",(real(size(Shat),DP)*3+size(Xnew)+size(Ds))*DP)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ENDSUBROUTINE first_subspace_stopw_gamma
    !---------------------------------------------------------
   SUBROUTINE Rayleigh_quotient_gamma(Ik,veff,nst,x,xhx)
     !xhx=(X,H,X)



      ! USE ScaLapack_module, ONLY: SL_matmat
      USE ScaLapack_module, ONLY: SL_matmat_realtn
      USE parameters, ONLY: BLOCK_MBNB,Nstates_global
      USE Grid_module, ONLY: global_n

      use m_time_evaluate ,only: memory_sum,memory_free,filename
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN) :: Ik &
                              &,nst
      REAL(DP),INTENT(IN) :: veff(:,:,:)
      REAL(DP) :: x(:,:),xhx(:,:)
      !LOCAL
      REAL(DP) :: hx(n,nst)
      integer(I4B) :: i
      ! COMPLEX(DP) :: hx_global(global_n)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      call memory_sum("first_Rayleigh_quotient",real(size(hx),DP)*DP)
      CALL cal_HX_gamma(Ik,veff,nst,x,hx)
      !xhx(:,:)=MATMUL( TRANSPOSE(CONJG(x)) , hx )



      ! CALL SL_matmat('c','n',x,hx,xhx,global_n,Nstates_global,&
      CALL SL_matmat_realtn('t','n',x,hx,xhx,global_n,Nstates_global,&
           &global_n,Nstates_global,1,BLOCK_MBNB,1,BLOCK_MBNB,&
           & BLOCK_MBNB,BLOCK_MBNB)

      call memory_free("first_Rayleigh_quotient",real(size(hx),DP)*DP)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ENDSUBROUTINE Rayleigh_quotient_gamma
   !-------------------------------HX------------------------------
   SUBROUTINE cal_HX_gamma(Ik,veff,nst,V,HV)
      USE matvec_module , ONLY : rmatvec_gamma
      IMPLICIT NONE
      !INOUT
      INTEGER(I4B),INTENT(IN) :: Ik,nst
      REAL(DP),INTENT(IN)  :: veff(:,:,:) !effective potential
      REAL(DP),INTENT(IN)  :: V(:,:)
      REAL(DP),INTENT(OUT) :: HV(:,:)
      !LOCAL
      INTEGER(I4B) :: Is
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! print*,'xlt-test     nst',nst
      DO Is=1,nst
      ! print*,'start',parallel%myid
         CALL rmatvec_gamma(veff,Ik,V(:,Is),HV(:,Is),n)
      ! print*,'end',parallel%myid
      ENDDO
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ENDSUBROUTINE cal_HX_gamma
   SUBROUTINE cheby_filter_RR_gamma(Ik,veff,X,D)
      USE parameters , ONLY : LRROrthNorm
      USE Lapack_module ,  ONLY : OrthNorm

      USE parameters, ONLY: Nstates_global

      IMPLICIT NONE
      !INOUT
      INTEGER(I4B),INTENT(IN) :: Ik
      REAL(DP),INTENT(IN) :: veff(:,:,:)
      REAL(DP),INTENT(INOUT) :: X(:,:)
      REAL(DP),INTENT(OUT) :: D(:)  !rayleigh-ritz value
      !LOCAL
      REAL(DP) :: a,b,al

      INTEGER(I4B) :: bcastid,i,j

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !Boundary
      a=MAXVAL(D)
      al=MINVAL(D)



      !> find the core deal with the highest state
      aa:DO j=parallel%dims(1),1,-1 !> all core
      DO i=size(parallel%sub2sum,1),1,-1  !> all states per core
         IF(parallel%sub2sum(i,j)==Nstates_global)THEN
            bcastid=j-1
         CALL Estupb_gamma(7,Ik,veff,X(:,parallel%nstate_proc),b)
            exit aa
         ENDIF
      ENDDO
      ENDDO aa
      CALL MPI_BCAST(b,1,MPI_REAL8,bcastid,parallel%commy,mpinfo)

      !filtering
      CALL chebyshev_filter_scaled_gamma(Ik,veff,X,CheM,a,b,al)
      !CALL chebyshev_filter(kp,veff,X,CheM,a,b)
      !Rayleigh-Ritz step
      !Generalize RR
      CALL GRayleigh_Ritz_gamma(Ik,veff,X,D)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ENDSUBROUTINE cheby_filter_RR_gamma
   SUBROUTINE Estupb_gamma(k,Ik,veff,vec,b)
      !
      USE matvec_module , ONLY : rmatvec_gamma
      USE Lapack_module , ONLY : Norm_2,diagM_real
      IMPLICIT NONE
      !INOUT
      INTEGER(I4B),INTENT(IN) :: k , Ik
      !hamiltonian
      REAL(DP),INTENT(IN) :: veff(:,:,:)
      REAL(DP),INTENT(IN) :: vec(:)
      REAL(DP),INTENT(OUT) :: b
      !LOCAL
      REAL(DP) :: v0(n),v(n),f(n),alpha
      REAL(DP) :: beta,eval(k),mz
      REAL(DP) :: T(k,k),evec(k,k)
      INTEGER(I4B) :: J,Nj
      REAL(DP),external  :: ddot

      REAL(DP) :: alpha_local,b_local
      REAL(DP) :: beta_local

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      v=vec !should normalize
      T(:,:)=0.d0
      Nj=MIN(k,10)
      !
      CALL rmatvec_gamma(veff,Ik,v,f,n)



      alpha_local=ddot(size(f),f,1,v,1)
      call mpi_allreduce(alpha_local,alpha,1,mpi_real8&
           & ,mpi_sum, parallel%commx, mpinfo)

      f=f-alpha*v
      T(1,1)=alpha
      DO J=2,Nj
         ! beta=SQRT(REAL(DOT_PRODUCT(f,f),8))



         beta_local=ddot(size(f),f,1,f,1)
         call mpi_allreduce(beta_local, beta,1,mpi_real8&
              &, mpi_sum,parallel%commx,mpinfo)

         beta=sqrt(beta)
         v0=v
         v=f/beta
         CALL rmatvec_gamma(veff,Ik,v,f,n)
         f=f-beta*v0
         ! alpha=DOT_PRODUCT(f,v)



      alpha_local=ddot(size(f),f,1,v,1)
      call mpi_allreduce(alpha_local,alpha,1,mpi_real8&
           & ,mpi_sum, parallel%commx, mpinfo)

         f=f-alpha*v
         T(J,J-1)=beta
         T(J-1,J)=beta
         T(J,J)=alpha
      ENDDO
      !NORM2(T)
      !b=Norm_2(T,k) + SQRT(REAL(DOT_PRODUCT(f,f),8))
      CALL diagM_real(T,evec,eval)
      !mz=MAX( ABS(evec(k,k)), ABS(evec(k,k-1)) , ABS(evec(k,k-2)) )
      mz=ABS(evec(k,k))



      b_local=ddot(size(f),f,1,f,1)
      call mpi_allreduce(b_local,b,1,mpi_real8,mpi_sum&
           &, parallel%commx, mpinfo)

      b=eval(k) + SQRT(REAL(b,DP))*mz
      !STOP
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ENDSUBROUTINE Estupb_Gamma
   SUBROUTINE chebyshev_filter_scaled_gamma(Ik,veff,X,m,a,b,al)
      !
      IMPLICIT NONE
      !INOUT
      INTEGER(I4B),INTENT(IN) :: Ik !k-point
      REAL(DP),INTENT(IN) :: veff(:,:,:) !veff
      REAL(DP),INTENT(INOUT) :: X(:,:)
      INTEGER(I4B),INTENT(IN)   :: m ! the m degree Chebyshev polynomial we used
      REAL(DP),INTENT(IN) :: a,b,al  !interval [a,b] to be filter
      !LOCAL
      REAL(DP) :: e &   !(b-a)/2
            & ,   c &   !(b+a)/2
            & , sigma,sigmanew,tau
      !
      REAL(DP),DIMENSION(n,sn)  :: HV,Y,Ynew
      INTEGER(I4B) :: Ic
      REAL(DP) :: temp1,temp2
!xq
!INTEGER(I4B) :: t1,t2
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF(m<2)THEN
         WRITE(6,*) 'Chebyshev_filter: the degree m must larger than 1'
      ENDIF

      e = ( b - a ) / 2
      c = ( b + a ) / 2
      sigma=e / (c-al)
      tau=2.d0 / sigma
!CALL system_clock(t1)
      CALL cal_HX_gamma(Ik,veff,sn,X,HV)
!CALL system_clock(t2)
!print*,'HX time:',(t2-t1)/10000.d0
      temp1=sigma / e
      Y= ( HV - c*X ) * temp1
      DO Ic=2,m
         sigmanew=1.d0 / (tau-sigma)
         CALL cal_HX_gamma(Ik,veff,sn,Y,HV)
         temp1=2.d0*sigmanew/e
         temp2=sigma*sigmanew
         Ynew=( HV - c*Y )*temp1 - temp2*X
         X=Y
         Y=Ynew
         sigma=sigmanew
      ENDDO

      X=Ynew
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE chebyshev_filter_scaled_gamma
   SUBROUTINE GRayleigh_Ritz_gamma(Ik,veff,X,D)
      USE Lapack_module , ONLY : GeneralizeEigen,matmat

      USE parameters, ONLY: Nstates_global, BLOCK_MBNB
      USE ScaLapack_module, ONLY: SL_GeneralizeEigen_real2,twoD_map,&
           & SL_matmat_realtn,SL_matmat_realnn
           ! & SL_matmat
      USE Grid_module, ONLY: global_n
      ! USE m_time_evaluate, ONLY: filename

      IMPLICIT NONE
      !IN/OUT
      INTEGER(I4B),INTENT(IN) :: Ik
      REAL(DP),INTENT(IN) :: veff(:,:,:)
      REAL(DP),INTENT(INOUT) :: X(:,:)
      REAL(DP),INTENT(OUT) :: D(:)
      !



      REAL(DP),DIMENSION(twoD_map(1,parallel%rankx,&
           &parallel%ranky),sn) :: S_hat,H_hat,Q

      REAL(DP) :: Xnew(n,sn)

      INTEGER(I4B) :: mb,nb,cmb,cnb

      ! integer(i4b) :: i
!xq
!INTEGER(I4B) :: t1,t2,t3,t4,t5,t6
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !Calculate the overlap matrix
!call system_clock(t1)
      ! x=cmplx(0.d0,0.d0)






      ! do i=1,sn,1
      !    X(parallel%sub2sum(i,parallel%myid+1),i)=&
      !         & cmplx(parallel%sub2sum(i,parallel%myid+1),0)
      ! enddo
      mb=1
      nb=1!BLOCK_MBNB
      cmb=1!BLOCK_MBNB
      cnb=1!BLOCK_MBNB
      ! CALL SL_matmat('c','n',X,X,S_hat,&
      CALL SL_matmat_realtn('t','n',X,X,S_hat,&
           &global_n,Nstates_global,global_n,Nstates_global,mb,nb,mb,nb,cmb,cnb)

!call system_clock(t2)
      !Calculate the project hamiltion
!call system_clock(t3)
      CALL Rayleigh_quotient_gamma(Ik,veff,sn,X,H_hat)
!call system_clock(t4)
      !solve the generalized eigenvalue problem
!call system_clock(t5)



      CALL SL_GeneralizeEigen_real2(Nstates_global,H_hat,S_hat,Nstates_global&
           &,Nstates_global,Nstates_global,Nstates_global&
           &,cmb,cnb,cmb,cnb,Q,Nstates_global,Nstates_global&
           &,D,cmb,cnb)

!call system_clock(t6)
      !X=XQ



      ! CALL SL_matmat('n','n',X,Q,Xnew,global_n,Nstates_global&
      CALL SL_matmat_realnn('n','n',X,Q,Xnew,global_n,Nstates_global&
           &,Nstates_global,Nstates_global,mb,nb,cmb,cnb,cmb,cnb)

!print*,'Overlap',(t2-t1)/10000.d0
!print*,'Hhat',(t4-t3)/10000.d0
!print*,'Genera',(t6-t5)/10000.d0
!STOP
      X=Xnew
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ENDSUBROUTINE GRayleigh_Ritz_Gamma
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDMODULE chebyshev_module

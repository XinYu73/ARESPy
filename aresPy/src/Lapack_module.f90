MODULE Lapack_module
!###########################################################!
!*For: Lapack                                               !
!*Author:Qiang Xu                                           !
!*Date:2017-7-20                                            !
!###########################################################!
  USE constants
  !USE MKL95_PRECISION
  IMPLICIT NONE
  INTEGER(I4B) :: MMAX=1000000,NMAX=10000
CONTAINS
  !----------------------------------------------------------
  !#########################################################!
  !  CALL ZHEEV(jobs,uplo,n,a,lad,w,work,lwork,rwork,info)  !
  !#########################################################!
  SUBROUTINE diagM(mat,evec,eval)
     !
     IMPLICIT NONE
     COMPLEX(DCP),INTENT(IN)  :: mat(:,:) !IN:martix
     COMPLEX(DCP),INTENT(OUT) :: evec(:,:) !OUT:engivector
     REAL(DP),INTENT(OUT) :: eval(:)     !engivalue
     !
     INTEGER(I4B) :: lda            !see Lapack's ZHEEV
     INTEGER(I4B) :: lwork !,lwmax=1000
     COMPLEX(DCP) :: work(1000)
     REAL(DP),ALLOCATABLE     :: rwork(:)
     INTEGER(I4B) :: info
     INTEGER(I4B) :: dime
     !
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     dime=SIZE(mat,1)
     IF(dime>MMAX)THEN
        WRITE(6,*) 'diagM:Check the dimension of MMAX'
        STOP
     ENDIF
     lda=dime
     !---creat arrays
     ALLOCATE(rwork(3*dime-2))
     !use mat for lapack
     evec(:,:)=mat(:,:)
     !query the optimal workspace
     lwork=-1
     CALL ZHEEV('V','U',dime,evec,lda,eval,work,lwork,rwork,info)
     !lwork=MIN(lwmax,INT(work(1)))
     lwork=MIN(2*dime,INT(work(1)))

     !solve engivalue
     CALL ZHEEV('V','U',dime,evec,lda,eval,work,lwork,rwork,info)
     IF(info/=0)THEN
        print*,'ZHEEV err:info',info
        STOP
     ENDIF
     !---destroy array
     DEALLOCATE(rwork)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE diagM
  !----------------------------------------------------------
  !###################################################################!
  ! CALL ZHEGV(itype,JOBZ,uplo,N,A,lDA,B,LDB,w,work,lwork,rwork,INFO) !
  ! itype = 1 : Ax=eBx                                                !
  !         2 : ABx=ex                                                !
  !         3 : BAx=ex                                                !
  !###################################################################!
  SUBROUTINE GeneralizeEigen(dime,matA,matB,evec,eval)
    use m_time_evaluate, only: memory_sum,memory_free
     IMPLICIT NONE
     INTEGER(I4B),INTENT(IN) :: dime
     COMPLEX(DCP),INTENT(IN)  :: matA(:,:),matB(:,:)
     COMPLEX(DCP),INTENT(OUT) :: evec(:,:)
     REAL(DP),INTENT(OUT) :: eval(:)
     !LOCAL
     INTEGER(I4B) :: lda,ldb,info,lwork
     COMPLEX(DCP) :: U(dime,dime)
     COMPLEX(DCP) :: work(1000)
     REAL(DP)    :: rwork(3*dime-2)
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     call memory_sum("CeneralizeEigen",(real(size(U),DP)+size(work)*DP+size(rwork)*DP))
     IF(dime>MMAX)THEN
        WRITE(6,*) 'diagM:Check the dimension of MMAX'
        STOP
     ENDIF
     lda=dime
     ldb=dime
     !
     evec(:,:)=matA(:,:)
     U(:,:)=matB(:,:)
     !query the optimal workspace
     lwork=-1
     info=0
     CALL ZHEGV(1,'V','U',dime,evec,lda,U,ldb,eval,work,lwork,rwork,info)
     lwork=MIN(2*dime-1,INT(work(1)))
     !solve the eigen problem
     CALL ZHEGV(1,'V','U',dime,evec,lda,U,ldb,eval,work,lwork,rwork,info)
     !test
     IF(info/=0)THEN
        print*,'ZHEGV err:info',info
        STOP
     ENDIF
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
     call memory_free("CeneralizeEigen",(real(size(U),DP)+size(work)*DP+size(rwork)*DP))
  ENDSUBROUTINE GeneralizeEigen
  !#########################################################!
  !      CALL ZGEQRF(M,N,A,LDA,TAU,WORK,LWORK,INFO)         !
  !#########################################################!
   SUBROUTINE OrthNorm(mat)
      !
      IMPLICIT NONE
      COMPLEX(DCP) :: mat(:,:)
      !LOCAL
      INTEGER(I4B) :: m,n,lda,lwork,info
      COMPLEX(DCP) :: tau(SIZE(mat,2))
      COMPLEX(DCP) :: work(5000)!work(SIZE(mat,2))
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      m=SIZE(mat,1)
      n=SIZE(mat,2)
      lda=MAX(1,m)
      !
      !!!xlt!IF(m>MMAX .OR. n>NMAX .OR. m<n) THEN
      !!!xlt!   print*,'m,n',m,n
      !!!xlt!   WRITE(6,*) 'OrthNorm:Check the dimension of MMAX,NMAX,M,N'
      !!!xlt!   STOP
      !!!xlt!ENDIF
      !query the optimal workspace
      lwork=-1
      CALL ZGEQRF(m,n,mat,lda,tau,work,lwork,info)
      lwork=MAX(1,INT(work(1)))
      !ALLOCATE(work(lwork+2))
      !QR
      CALL ZGEQRF(m,n,mat,lda,tau,work,lwork,info)
      IF(info/=0)THEN
        print*,'ZGEQRF err:info',info
        STOP
      ENDIF
      !orth-normal
      CALL ZUNGQR(m,n,n,mat,lda,tau,work,lwork,info)
      IF(info/=0)THEN
        print*,'ZGEQRF err:info',info
        STOP
      ENDIF
      !
      !DEALLOCATE(work)
      !CALL X04DBF('General',' ',M,N,mat,LDA,'Bracketed',F7.4,'heelo','')
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE OrthNorm
  !-------------------------Norm2----------------------------
  FUNCTION Norm_2(mat,k)
     !
     IMPLICIT NONE
     COMPLEX(DCP),INTENT(IN) :: mat(:,:)
     INTEGER(I4B),INTENT(IN) :: k
     REAL(DP) :: Norm_2
     !
     COMPLEX(DCP) :: MM(k,k),evec(k,k)
     REAL(DP) :: eval(k)
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     MM(:,:)=MATMUL(TRANSPOSE(CONJG(mat)),mat)
     CALL diagM(MM,evec,eval)
     !
     Norm_2=SQRT((eval(k)))
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDFUNCTION NorM_2
  !#############################################################!
  ! matrix: C = op(A) * op(B)                                   !
  ! CALL ZGEMM(TRANSA,TRANSB,M,N,K,alpha,A,LDA,B,LDB,beta,c,LDC)!
  ! op*='N','T','C' for op(X)=(X), (X)' and (X*)'               !
  !#############################################################!
  !-------------------------matmat---------------------------
  SUBROUTINE matmat(matA,matB,opA,opB,matC)
     IMPLICIT NONE
     COMPLEX(DCP),INTENT(IN)  :: matA(:,:),matB(:,:)
     COMPLEX(DCP),INTENT(OUT) :: matC(:,:)
     CHARACTER(1),INTENT(IN) :: opA,opB
     !
     INTEGER(I4B) :: LDA,LDB,LDC,M,N,K
     COMPLEX(DCP)  :: alpha=1.d0  &
             &    ,  beta=0.d0
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     !C
     LDC=SIZE(matC,1)
     M=LDC
     N=SIZE(matC,2)
     !A
     IF(opA=='N'.OR.opA=='n')THEN
        K=SIZE(matA,2)
        LDA=MAX(1,M)
     ELSE
        K=SIZE(matA,1)
        LDA=MAX(1,K)
     ENDIF
     !B
     IF(opB=='N'.OR.opB=='n')THEN
        LDB=MAX(1,K)
     ELSE
        LDB=MAX(1,N)
     ENDIF
     !
     CALL ZGEMM(opA,opB,M,N,K,alpha,matA,LDA,matB,LDB,beta,matC,LDC)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE matmat
  !-------------------------invmat---------------------------
  SUBROUTINE invmat(mat)
     IMPLICIT NONE
     COMPLEX(DCP),INTENT(INOUT) :: mat(:,:)
     !LOCAL
     INTEGER(I4B) :: info,lda,m,n
     INTEGER(I4B) :: IPIV(SIZE(mat,1))
     !
     INTEGER(I4B) :: lwork
     COMPLEX(DCP) :: wkt(10)
     COMPLEX(DCP),ALLOCATABLE :: work(:)
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     n=SIZE(mat,1)
     m=n
     lda=n
     !LU decomposed
     CALL ZGETRF(m,n,mat,lda,IPIV,info)
     IF(info/=0)THEN
        print*,'DGETRF err:info',info
     ENDIF
     !inverse
     !inquire the work matrix
     lwork=-1
     CALL ZGETRI(n,mat,lda,IPIV,wkt,lwork,info)
     lwork=INT(REAL(wkt(1),DP))
     !allocate
     ALLOCATE(work(lwork))
     CALL ZGETRI(n,mat,lda,IPIV,work,lwork,info)
     IF(info/=0)THEN
        print*,'ZGETRF err:info',info
     ENDIF
     !deallocate
     DEALLOCATE(work) 
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE invmat

   !#############################################################!
   !for real matrix (NAG Library)                                !
   !#############################################################!
   !SUBROUTINE OrthNorm_nag(mat)
   !   IMPLICIT NONE
   !   REAL(DP),INTENT(INOUT) :: mat(:,:)
   !   !LOCAL
   !   INTEGER(I4B) :: lda,m,n1,n2,ICOL,IFAIL
   !   REAL(DP) :: S(SIZE(mat,2)),CC
   !   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !   m=SIZE(mat,1)
   !   n1=1
   !   n2=SIZE(mat,2)
   !   lda=MAX(1,m)
   !   !
   !   IFAIL=0
   !   CALL F05AAF(mat,lda,m,n1,n2,S,CC,ICOL,IFAIL)
   !   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   !ENDSUBROUTINE OrthNorm_nag
  !#########################################################!
  !      CALL DGEQRF(M,N,A,LDA,TAU,WORK,LWORK,INFO)         !
  !#########################################################!
   SUBROUTINE OrthNorm_real(mat)
      !
      IMPLICIT NONE
      REAL(DP) :: mat(:,:)
      !LOCAL
      INTEGER(I4B) :: m,n,lda,lwork,info
      REAL(DP) :: tau(SIZE(mat,2))
      !
      REAL(DP) :: wkt(10)
      REAL(DP),ALLOCATABLE :: work(:)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      m=SIZE(mat,1)
      n=SIZE(mat,2)
      lda=MAX(1,m)
      !
      !!!!xlt!IF(m>MMAX .OR. n>NMAX .OR. m<n) THEN
      !!!!xlt!   print*,'m,n',m,n
      !!!!xlt!   WRITE(6,*) 'OrthNorm:Check the dimension of MMAX,NMAX,M,N'
      !!!!xlt!   STOP
      !!!!xlt!ENDIF
      !query the optimal workspace
      lwork=-1
      CALL DGEQRF(m,n,mat,lda,tau,wkt,lwork,info)
      lwork=INT(wkt(1))
      !ALLOCATE
      ALLOCATE(work(lwork))
      !QR
      CALL DGEQRF(m,n,mat,lda,tau,work,lwork,info)
      IF(info/=0)THEN
        print*,'DGEQRF err:info',info
        STOP
      ENDIF
      !orth-normal
      CALL DORGQR(m,n,n,mat,lda,tau,work,lwork,info)
      IF(info/=0)THEN
        print*,'ZGEQRF err:info',info
        STOP
      ENDIF
      !Destory
      DEALLOCATE(work)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE OrthNorm_real
  !#########################################################!
  !  CALL DSYEV(jobs,uplo,n,a,lda,w,work,lwork,info)  !
  !#########################################################!
  SUBROUTINE diagM_real(mat,evec,eval)
     !
    USE m_time_evaluate, ONLY: memory_sum,memory_free
     IMPLICIT NONE
     REAL(DP),INTENT(IN)  :: mat(:,:) !IN:martix
     REAL(DP),INTENT(OUT) :: evec(:,:) !OUT:engivector
     REAL(DP),INTENT(OUT) :: eval(:)     !engivalue
     !
     INTEGER(I4B) :: lda            !see Lapack's ZHEEV
     INTEGER(I4B) :: lwork !,lwmax=1000
     REAL(DP) :: wkt(10)
     REAL(DP),ALLOCATABLE :: work(:)
     INTEGER(I4B) :: info
     INTEGER(I4B) :: dime
     !
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     dime=SIZE(mat,1)
     IF(dime>MMAX)THEN
        WRITE(6,*) 'diagM:Check the dimension of MMAX'
        STOP
     ENDIF
     lda=dime
     !use mat for lapack
     evec(:,:)=mat(:,:)

     !query the optimal workspace
     lwork=-1
     CALL DSYEV('V','U',dime,evec,lda,eval,wkt,lwork,info)
     !lwork=MIN(lwmax,INT(work(1)))
     lwork=INT(wkt(1))
     !allocate
     ALLOCATE(work(lwork))
     call memory_sum('diagM_real',real(size(work),DP)*DP)

     !solve engivalue
     CALL DSYEV('V','U',dime,evec,lda,eval,work,lwork,info)
     IF(info/=0)THEN
        print*,'DSYEV err:info',info
        STOP
     ENDIF

     !---destroy array
     call memory_free('diagM_real',real(size(work),DP)*DP)
     DEALLOCATE(work)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE diagM_real
  !###################################################################!
  ! CALL DSYGV(itype,JOBZ,uplo,N,A,lDA,B,LDB,W,work,lwork,INFO)       !
  ! itype = 1 : Ax=eBx                                                !
  !         2 : ABx=ex                                                !
  !         3 : BAx=ex                                                !
  !###################################################################!
  SUBROUTINE GeneralizeEigen_real(dime,matA,matB,evec,eval)
     IMPLICIT NONE
     INTEGER(I4B),INTENT(IN) :: dime
     REAL(DP),INTENT(IN)  :: matA(:,:),matB(:,:)
     REAL(DP),INTENT(OUT) :: evec(:,:)
     REAL(DP),INTENT(OUT) :: eval(:)
     !LOCAL
     INTEGER(I4B) :: lda,ldb,info,lwork
     REAL(DP) :: U(dime,dime)
     REAL(DP) :: wkt(10)
     REAL(DP),ALLOCATABLE :: work(:)
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     IF(dime>MMAX)THEN
        WRITE(6,*) 'diagM:Check the dimension of MMAX'
        STOP
     ENDIF
     lda=dime
     ldb=dime
     !
     evec(:,:)=matA(:,:)
     U(:,:)=matB(:,:)
     !query the optimal workspace
     lwork=-1
     info=0
     CALL DSYGV(1,'V','U',dime,evec,lda,U,ldb,eval,wkt,lwork,info)
     lwork=INT(wkt(1))
     !allocate
     ALLOCATE(work(lwork))
     !solve the eigen problem
     CALL DSYGV(1,'V','U',dime,evec,lda,U,ldb,eval,work,lwork,info)
     !test
     IF(info/=0)THEN
        print*,'DSYGV err:info',info
        STOP
     ENDIF
     !destroy
     DEALLOCATE(work)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE GeneralizeEigen_real
  !------------------selected eigenvalues--------------------
  SUBROUTINE diagMX_real(mat,dime,num,Il,Iu,evec,eval)
     IMPLICIT NONE
     !INOUT
     INTEGER(I4B),INTENT(IN) :: Il & !the lower index of eigen-pair
                            & , Iu,num & !the upper index of eigen-pair
                            & , dime !the dimension of mat
     REAL(DP),INTENT(IN) :: mat(dime,dime) !matrix
     REAL(DP),INTENT(OUT) :: evec(dime,num) &!eigen vector
                         &,  eval(num)   !eigen value
     !LOCAL
     REAL(DP) :: A(dime,dime)
     REAL(DP) :: vl,vu,wkt(10),w(dime) !,z(dime,num)
     INTEGER(I4B) :: lda,lwork,iwork(5*dime),m
     REAL(DP) :: abstol=1e-6
     INTEGER(I4B) :: ifail(dime),info,ldz
     REAL(DP),ALLOCATABLE :: work(:)
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     A(:,:)=mat(:,:)
     lda=dime
     ldz=dime
     m=num
     !query the optimal workspace
     lwork=-1
     CALL DSYEVX('V','I','U',dime,A,lda,vl,vu,Il,Iu,abstol &
         & ,m,w,evec,ldz,wkt,lwork,iwork,ifail,info)
     lwork=INT(wkt(1))
     !work array
     ALLOCATE(work(lwork))
     CALL DSYEVX('V','I','U',dime,A,lda,vl,vu,Il,Iu,abstol &
         & ,m,w,evec,ldz,work,lwork,iwork,ifail,info)
     IF(info/=0)THEN
        print*,'DSYEVX err:info',info
        STOP
     ENDIF
     !destory
     !DEALLOCATE(work)
     !-------------------OUT PUT-----------------------------
     !eigen-values
     eval(1:num)=w(1:num)
     DEALLOCATE(work)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE diagMX_real
  !-------------------------matmat---------------------------
  SUBROUTINE matmat_real(matA,matB,opA,opB,matC)
     IMPLICIT NONE
     REAL(DP),INTENT(IN)  :: matA(:,:),matB(:,:)
     REAL(DP),INTENT(OUT) :: matC(:,:)
     CHARACTER(1),INTENT(IN) :: opA,opB
     !
     INTEGER(I4B) :: LDA,LDB,LDC,M,N,K
     REAL(DP)  :: alpha=1.d0  &
             &    ,  beta=0.d0
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     !C
     matC(:,:)=0.d0
     LDC=SIZE(matC,1)
     M=LDC
     N=SIZE(matC,2)
     !A
     IF(opA=='N'.OR.opA=='n')THEN
        K=SIZE(matA,2)
        LDA=MAX(1,M)
     ELSE
        K=SIZE(matA,1)
        LDA=MAX(1,K)
     ENDIF
     !B
     IF(opB=='N'.OR.opB=='n')THEN
        LDB=MAX(1,K)
     ELSE
        LDB=MAX(1,N)
     ENDIF
     !
     CALL DGEMM(opA,opB,M,N,K,alpha,matA,LDA,matB,LDB,beta,matC,LDC)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE matmat_real
  !#########################################################!
  !          CALL DPOTRF(UPLO,N,A,LDA,INFO)                 !
  !#########################################################!
  !-------------------Cholesky factorization-----------------
  SUBROUTINE Cholesky_factor_real(mat)
     IMPLICIT NONE
     !IN/OUT
     REAL(DP),INTENT(INOUT) :: mat(:,:)
     !LOCAL
     INTEGER(I4B) :: info,lda,dimen
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     dimen=SIZE(mat,1)
     lda=dimen
     !
     CALL DPOTRF('U',dimen,mat,lda,info)
     !test
     IF(info/=0)THEN
        print*,'Cholesky err:info',info
        STOP
     ENDIF
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE Cholesky_factor_real
  !-----------------------inv_mat----------------------------
  SUBROUTINE invmat_real(mat)
     IMPLICIT NONE
     REAL(DP),INTENT(INOUT) :: mat(:,:)
     !LOCAL
     INTEGER(I4B) :: info,lda,m,n
     INTEGER(I4B) :: IPIV(SIZE(mat,1))
     !
     INTEGER(I4B) :: lwork
     REAL(DP) :: wkt(10)
     REAL(DP),ALLOCATABLE :: work(:)
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     n=SIZE(mat,1)
     m=n
     lda=n
     !LU decomposed
     CALL DGETRF(m,n,mat,lda,IPIV,info)
     IF(info/=0)THEN
        print*,'DGETRF err:info',info
     ENDIF
     !inverse
     !inquire the work matrix
     lwork=-1
     CALL DGETRI(n,mat,lda,IPIV,wkt,lwork,info)
     lwork=wkt(1)
     !allocate
     ALLOCATE(work(lwork))
     CALL DGETRI(n,mat,lda,IPIV,work,lwork,info)
     IF(info/=0)THEN
        print*,'DGETRF err:info',info
     ENDIF
     !deallocate
     DEALLOCATE(work) 
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE invmat_real
  !----------------------------------------------------------
ENDMODULE Lapack_module

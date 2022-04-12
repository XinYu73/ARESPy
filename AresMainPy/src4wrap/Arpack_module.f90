MODULE Arpack_module
!###########################################################!
!*For: Arpack                                               !
!*Author:Qiang Xu                                           !
!*Date:2017-7-20                                            !
!###########################################################!
  USE constants
  USE grid_module , ONLY : n,n1,n2,n3
  IMPLICIT NONE
  INTEGER(I4B),PARAMETER :: maxn=1000000,maxnev=500,maxncv=10000
CONTAINS
  !----------------------------------------------------------
  !#########################################################!
  !For complex hamiltonial
  !CALL znaupd
  !CALL zneupd
  !#########################################################!
  !SUBROUTINE diagH_arpack(mat,nev,evec,eval,resid_restart,nec,info)
  SUBROUTINE diagH_arpack(veff,Ik,nev,evec,eval,resid_restart,nec,info,maxmvs,TOL)
     USE matvec_module , ONLY : cmatvec
     !
     IMPLICIT NONE

     ! INCLUDE 'debug.h'
     !  integer  logfil, ndigit, mgetv0,                                     &
     ! &         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,     &
     ! &         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,     &
     ! &         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
     !  common /debug/                                                       &
     ! &         logfil, ndigit, mgetv0,                                     &
     ! &         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,     &
     ! &         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,     &
     ! &         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
     !
     REAL(DP),INTENT(IN)  :: veff(:,:,:)   !IN:veff
     INTEGER(I4B),INTENT(IN) :: Ik
     INTEGER(I4B),INTENT(IN)  :: nev      !number of state we need
     REAL(DP),INTENT(IN) :: TOL
     COMPLEX(DP),INTENT(OUT) :: evec(:,:) !OUT:engivector
     REAL(DP),INTENT(OUT) :: eval(nev)     !engivalue
     COMPLEX(DP),INTENT(INOUT) :: resid_restart(:)
     INTEGER(I4B),INTENT(INOUT) :: nec     !number of convergence
     INTEGER(I4B),INTENT(INOUT) :: info
     INTEGER(I4B),INTENT(INOUT) :: maxmvs
     !Arpack
     INTEGER(I4B)     :: ido,ncv,lworkl,mode1,         &
                    &    ishfts,maxitr,dimen,ldv,ierr
     INTEGER(I4B) :: iparam(11),ipntr(14),jj
     COMPLEX(DP) :: d(2*nev+5,2)
     COMPLEX(DP),ALLOCATABLE :: vtmp(:,:)
     COMPLEX(DP),ALLOCATABLE :: resid(:),workd(:),workl(:),workev(:)
     REAL(DP),ALLOCATABLE :: rwork(:)
     COMPLEX(DP) :: sigma
     !
     LOGICAL :: rvec,select(2*nev+5)
     !
     CHARACTER(LEN=1) :: bmat
     CHARACTER(LEN=2) :: which
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     !>>>
     !dimen=SIZE(,1)
     dimen=n
     ncv=2*nev+5
!     ncv=4
     ldv=dimen
     IF( dimen > maxn )THEN
        WRITE(6,*) 'ARPACK:ERROR with _NSIMP: N is greater than MAXN',dimen
        STOP
     ELSEIF ( nev > maxnev ) then
        WRITE(6,*) 'ARPACK:ERROR with _NSIMP: NEV is greater than MAXNEV',nev
        STOP
     ELSEIF ( ncv > maxncv ) then
        WRITE(6,*) 'ARPACK:ERROR with _NSIMP: NCV is greater than MAXNCV',ncv
        STOP
     ENDIF
     !<<<

     !debug variables>>>>
     !ndigit = -3
     !logfil = 6
     !mcaitr = 0
     !mcapps = 0
     !mcaupd = 1
     !mcaup2 = 0
     !mceigh = 0
     !mceupd = 0
     !<<<<debug variables
     !Initialize the arpack
     iparam(:)=0
     ipntr(:)=0
!!!

     bmat='I'
     which='SR'  !'LR'/'SR' for largest/smallest real eigenvaule

     lworkl=3*ncv**2+5*ncv

     ido=0
     mode1=1
     ishfts=1
     maxitr=nev*50
!maxitr=100
     iparam(1)=ishfts
     iparam(3)=maxitr
     iparam(7)=mode1

     !allocate>>>
     ALLOCATE(vtmp(dimen,2*nev+5))
     ALLOCATE(resid(dimen))
     ALLOCATE(workd(3*dimen))
     ALLOCATE(workl(lworkl))
     ALLOCATE(rwork(ncv))
     ALLOCATE(workev(3*ncv))
     vtmp(:,:)=(0.d0,0.d0)
     !<<<allocate
     resid(:)=resid_restart(:)
     !Main loop>>>(Reverse Communication)
!check
!
10   CONTINUE
        CALL znaupd (ido, bmat, dimen, which, nev, tol, resid &
          &  ,ncv, vtmp,ldv, iparam, ipntr, workd, workl,  &
          &  lworkl,rwork, info)
        !need to be change--------------------->
!print*,'ido',ido
        IF(ido==1.OR.ido==-1)THEN
           !CALL matvec (mat,workd(ipntr(1):ipntr(1)+dimen-1),workd(ipntr(2):ipntr(2)+dimen-1),1,dimen)
           CALL cmatvec(veff,Ik,workd(ipntr(1)),workd(ipntr(2)),dimen)
           IF(iparam(9)<maxmvs) GOTO 10
        ENDIF
     !<<<Main loop
     !debug
     IF ( info == 1) THEN
        WRITE(6,'(/,a,/)') 'ARPACK:Maximum number of iterations reached.'
        !STOP
     ELSEIF ( info .eq. 3) THEN
        WRITE(6,'(/,a,a,/)') 'ARPACK:No shifts could be applied during',  &
              &  'implicit Arnoldi update, try increasing NCV.'
        STOP
     ENDIF

     IF ( info .lt. 0 ) THEN
        WRITE(6,'(/,a,i5)') 'ARPACK:Error with _saupd, info = ', info
        WRITE(6,'(a,/)') 'ARPACK:Check documentation in _saupd '
        STOP
     ELSE
        rvec=.TRUE.
        CALL zneupd (rvec, 'A',select, d, vtmp, ldv,    &
       &       sigma,workev,bmat, dimen, which, nev, tol, resid, ncv,  &
       &       vtmp, ldv,iparam, ipntr, workd, workl, lworkl  &
       &       ,rwork, ierr)

        IF( ierr /= 0)THEN
           WRITE(6,*) ' Error with _seupd, ierr = ', ierr
           WRITE(6,*) ' Check the documentation of _seupd. '
           STOP
        ENDIF

        resid_restart(:) = (0.d0,0.d0)

        DO jj= 1,nev
           resid_restart(:) = resid_restart(:)+vtmp(:,jj)
        ENDDO
        resid_restart(:) = resid_restart(:) / REAL(nev,DP)

!     Print additional convergence information
!        print *, '-------------------ARPACK--------------------'
!        print *, '_NSIMP '
!        print *, '====== '
!        print *, ' '
!        print *, ' Size of the matrix is ', dimen
!        print *, ' The number of Ritz values requested is ', nev
!        print *, ' The number of Arnoldi vectors generated',  &
!    &            ' (NCV) is ', ncv
!        print *, ' What portion of the spectrum: ', which
!        print *, ' The number of converged Ritz values is ',  &
!    &              nec
!        print *, ' The number of Implicit Arnoldi update',    &
!    &            ' iterations taken is ', iparam(3)
!        print *, ' The number of OP*x is ', iparam(9)
!        print *, ' The convergence criterion is ', tol
!        print *, '---------------------------------------------'


     ENDIF
     !output value
     !
     eval(:)=REAL(d(1:nev,1),8)
     evec(:,:)=vtmp(:,1:nev)
     !
     !next
     nec=iparam(5)
     maxmvs=iparam(9)
     !------
     IF(ALLOCATED(vtmp))   DEALLOCATE(vtmp)
     IF(ALLOCATED(resid))  DEALLOCATE(resid)
     IF(ALLOCATED(workd))  DEALLOCATE(workd)
     IF(ALLOCATED(workl))  DEALLOCATE(workl)
     IF(ALLOCATED(rwork))  DEALLOCATE(rwork)
     IF(ALLOCATED(workev)) DEALLOCATE(workev)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE diagH_arpack
  !----------------------------------------------------------
  !###################################!
  !*For: Arpack for real Hamiltonian  !
  !*Author:Qiang Xu                   !
  !*Date:2017-8-28                    !
  !*CALL dsaupd                       !
  !*CALL dseupd                       !
  !###################################!
  SUBROUTINE real_diagH_arpack(veff,nev,evec,eval,resid_restart,nec,info,maxmvs,TOL)
     !
     USE matvec_module , ONLY : rmatvec
     IMPLICIT NONE
     ! INCLUDE 'debug.h'
     !  integer  logfil, ndigit, mgetv0,                                     &
     ! &         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,     &
     ! &         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,     &
     ! &         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
     !  common /debug/                                                       &
     ! &         logfil, ndigit, mgetv0,                                     &
     ! &         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,     &
     ! &         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,     &
     ! &         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
     !
     REAL(DP),INTENT(IN)  :: veff(:,:,:)   !IN:veff
     INTEGER(I4B),INTENT(IN)  :: nev      !number of state we need
     REAL(DP),INTENT(OUT) :: evec(:,:) !OUT:engivector
     REAL(DP),INTENT(OUT) :: eval(nev)     !engivalue
     REAL(DP),INTENT(INOUT) :: resid_restart(:)
     INTEGER(I4B),INTENT(INOUT) :: nec     !number of convergence
     INTEGER(I4B),INTENT(INOUT) :: info
     INTEGER(I4B),INTENT(INOUT) :: maxmvs
     REAL(DP),INTENT(IN) :: TOL
     !Arpack
     INTEGER(I4B)     :: ido,ncv,lworkl,mode1,         &
                    &    ishfts,maxitr,dimen,ldv,ierr
     INTEGER(I4B) :: iparam(11),ipntr(11),jj
     REAL(DP) :: d(2*nev+5,2)
     REAL(DP),ALLOCATABLE :: vtmp(:,:)
     REAL(DP),ALLOCATABLE :: resid(:),workd(:),workl(:)
     REAL(DP) :: sigma
     !
     LOGICAL :: rvec,select(2*nev+5)
     !
     CHARACTER(LEN=1) :: bmat
     CHARACTER(LEN=2) :: which
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     !>>>
     !dimen=SIZE(,1)
     dimen=n
     ncv=2*nev+5
!     ncv=4
     ldv=dimen
     IF( dimen > maxn )THEN
        WRITE(6,*) 'ARPACK:ERROR with _NSIMP: N is greater than MAXN',dimen
        STOP
     ELSEIF ( nev > maxnev ) then
        WRITE(6,*) 'ARPACK:ERROR with _NSIMP: NEV is greater than MAXNEV',nev
        STOP
     ELSEIF ( ncv > maxncv ) then
        WRITE(6,*) 'ARPACK:ERROR with _NSIMP: NCV is greater than MAXNCV',ncv
        STOP
     ENDIF
     !<<<

     !debug variables>>>>
     !ndigit = -3
     !logfil = 6
     !mcaitr = 0
     !mcapps = 0
     !mcaupd = 1
     !mcaup2 = 0
     !mceigh = 0
     !mceupd = 0
     !<<<<debug variables
     !Initialize the arpack
     iparam(:)=0
     ipntr(:)=0
!!!

     bmat='I'
     which='SA'  !'LA'/'SA' for largest/smallest real eigenvaule

     lworkl=ncv*(ncv+8)

     ido=0
     mode1=1
     ishfts=1
     maxitr=nev*50
!maxitr=100
     iparam(1)=ishfts
     iparam(3)=maxitr
     iparam(4)=1
     iparam(7)=mode1

     !allocate>>>
     ALLOCATE(vtmp(dimen,ncv))
     ALLOCATE(resid(dimen))
     ALLOCATE(workd(3*dimen))
     ALLOCATE(workl(lworkl))
     vtmp(:,:)=(0.d0,0.d0)
     !<<<allocate
     resid(:)=resid_restart(:)
     !Main loop>>>(Reverse Communication)
!
10   CONTINUE
        CALL dsaupd (ido, bmat, dimen, which, nev, tol, resid &
          &  ,ncv, vtmp,ldv, iparam, ipntr, workd, workl,  &
          &  lworkl, info)
        !need to be change--------------------->
        IF(ido==1.OR.ido==-1)THEN
           CALL rmatvec(veff,workd(ipntr(1)),workd(ipntr(2)),dimen)
           IF(iparam(9)<maxmvs) GOTO 10
        ENDIF
     !<<<Main loop
     !debug
     IF ( info == 1) THEN
        WRITE(6,'(/,a,/)') 'ARPACK:Maximum number of iterations reached.'
        !STOP
     ELSEIF ( info .eq. 3) THEN
        WRITE(6,'(/,a,a,/)') 'ARPACK:No shifts could be applied during',  &
              &  'implicit Arnoldi update, try increasing NCV.'
        STOP
     ENDIF

     IF ( info .lt. 0 ) THEN
        WRITE(6,'(/,a,i5)') 'ARPACK:Error with _saupd, info = ', info
        WRITE(6,'(a,/)') 'ARPACK:Check documentation in _saupd '
        STOP
     ELSE
        rvec=.TRUE.
        CALL dseupd (rvec, 'A',select, d, vtmp, ldv,    &
       &       sigma,bmat, dimen, which, nev, tol, resid, ncv,  &
       &       vtmp, ldv,iparam, ipntr, workd, workl, lworkl  &
       &       , ierr)

        IF( ierr /= 0)THEN
           WRITE(6,*) ' Error with _seupd, ierr = ', ierr
           WRITE(6,*) ' Check the documentation of _seupd. '
           STOP
        ENDIF

        resid_restart(:) = (0.d0,0.d0)

        DO jj= 1,nev
           resid_restart(:) = resid_restart(:)+vtmp(:,jj)
        ENDDO
        resid_restart(:) = resid_restart(:) / REAL(nev,DP)

!     Print additional convergence information
!        print *, '-------------------ARPACK--------------------'
!        print *, '_NSIMP '
!        print *, '====== '
!        print *, ' '
!        print *, ' Size of the matrix is ', dimen
!        print *, ' The number of Ritz values requested is ', nev
!        print *, ' The number of Arnoldi vectors generated',  &
!    &            ' (NCV) is ', ncv
!        print *, ' What portion of the spectrum: ', which
!        print *, ' The number of converged Ritz values is ',  &
!    &              nec
!        print *, ' The number of Implicit Arnoldi update',    &
!    &            ' iterations taken is ', iparam(3)
!        print *, ' The number of OP*x is ', iparam(9)
!        print *, ' The convergence criterion is ', tol
!        print *, '---------------------------------------------'


     ENDIF
     !output value
     !
     eval(:)=REAL(d(1:nev,1),DP)
     evec(:,:)=vtmp(:,1:nev)
     !
     !next
     nec=iparam(5)
     maxmvs=iparam(9)
     !------
     IF(ALLOCATED(vtmp))   DEALLOCATE(vtmp)
     IF(ALLOCATED(resid))  DEALLOCATE(resid)
     IF(ALLOCATED(workd))  DEALLOCATE(workd)
     IF(ALLOCATED(workl))  DEALLOCATE(workl)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE real_diagH_arpack
  !----------------------------------------------------------
  !###################################!
  !*For: Arpack for real matrix       !
  !*Author:Qiang Xu                   !
  !*Date:2018-3-06                    !
  !*CALL dsaupd                       !
  !*CALL dseupd                       !
  !###################################!
  SUBROUTINE real_diagM_arpk(mat,nev,evec,eval,resid_restart,nec,info,maxmvs,TOL)
     !
     USE matvec_module , ONLY : rmatvec_new
     IMPLICIT NONE
     ! INCLUDE 'debug.h'
     !  integer  logfil, ndigit, mgetv0,                                     &
     ! &         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,     &
     ! &         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,     &
     ! &         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
     !  common /debug/                                                       &
     ! &         logfil, ndigit, mgetv0,                                     &
     ! &         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,     &
     ! &         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,     &
     ! &         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
     !
     REAL(DP),INTENT(IN)  :: mat(:,:)   !IN:Martix
     INTEGER(I4B),INTENT(IN)  :: nev      !number of state we need
     REAL(DP),INTENT(OUT) :: evec(:,:) !OUT:engivector
     REAL(DP),INTENT(OUT) :: eval(nev)     !engivalue
     REAL(DP),INTENT(INOUT) :: resid_restart(:)
     INTEGER(I4B),INTENT(INOUT) :: nec     !number of convergence
     INTEGER(I4B),INTENT(INOUT) :: info
     INTEGER(I4B),INTENT(INOUT) :: maxmvs
     REAL(DP),INTENT(IN) :: TOL
     !Arpack
     INTEGER(I4B)     :: ido,ncv,lworkl,mode1,         &
                    &    ishfts,maxitr,dimen,ldv,ierr
     INTEGER(I4B) :: iparam(11),ipntr(11),jj
     REAL(DP) :: d(2*nev+5,2)
     REAL(DP),ALLOCATABLE :: vtmp(:,:)
     REAL(DP),ALLOCATABLE :: resid(:),workd(:),workl(:)
     REAL(DP) :: sigma
     !
     LOGICAL :: rvec,select(2*nev+5)
     !
     CHARACTER(LEN=1) :: bmat
     CHARACTER(LEN=2) :: which
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     !>>>
     dimen=SIZE(mat,1)
     ncv=MIN(2*nev+5,dimen)
!     ncv=4
     ldv=dimen
     IF( dimen > maxn )THEN
        WRITE(6,*) 'ARPACK:ERROR with _NSIMP: N is greater than MAXN',dimen
        STOP
     ELSEIF ( nev > maxnev ) then
        WRITE(6,*) 'ARPACK:ERROR with _NSIMP: NEV is greater than MAXNEV',nev
        STOP
     ELSEIF ( ncv > maxncv ) then
        WRITE(6,*) 'ARPACK:ERROR with _NSIMP: NCV is greater than MAXNCV',ncv
        STOP
     ENDIF
     !<<<
     !Initialize the arpack
     iparam(:)=0
     ipntr(:)=0
!!!

     bmat='I'
     which='LA'  !'LA'/'SA' for largest/smallest real eigenvaule

     lworkl=ncv*(ncv+8)

     ido=0
     mode1=1
     ishfts=1
     maxitr=nev*50
!maxitr=100
     iparam(1)=ishfts
     iparam(3)=maxitr
     iparam(4)=1
     iparam(7)=mode1

     !allocate>>>
     ALLOCATE(vtmp(dimen,ncv))
     ALLOCATE(resid(dimen))
     ALLOCATE(workd(3*dimen))
     ALLOCATE(workl(lworkl))
     vtmp(:,:)=(0.d0,0.d0)
     !<<<allocate
     resid(:)=resid_restart(:)
     !Main loop>>>(Reverse Communication)
!
10   CONTINUE
        CALL dsaupd (ido, bmat, dimen, which, nev, tol, resid &
          &  ,ncv, vtmp,ldv, iparam, ipntr, workd, workl,  &
          &  lworkl, info)
        !need to be change--------------------->
        IF(ido==1.OR.ido==-1)THEN
           CALL rmatvec_new(mat,workd(ipntr(1)),workd(ipntr(2)),dimen)
           IF(iparam(9)<maxmvs) GOTO 10
        ENDIF
     !<<<Main loop
     !debug
     IF ( info == 1) THEN
        WRITE(6,'(/,a,/)') 'ARPACK:Maximum number of iterations reached.'
        !STOP
     ELSEIF ( info .eq. 3) THEN
        WRITE(6,'(/,a,a,/)') 'ARPACK:No shifts could be applied during',  &
              &  'implicit Arnoldi update, try increasing NCV.'
        STOP
     ENDIF

     IF ( info .lt. 0 ) THEN
        WRITE(6,'(/,a,i5)') 'ARPACK:Error with _saupd, info = ', info
        WRITE(6,'(a,/)') 'ARPACK:Check documentation in _saupd '
        STOP
     ELSE
        rvec=.TRUE.
        CALL dseupd (rvec, 'A',select, d, vtmp, ldv,    &
       &       sigma,bmat, dimen, which, nev, tol, resid, ncv,  &
       &       vtmp, ldv,iparam, ipntr, workd, workl, lworkl  &
       &       , ierr)

        IF( ierr /= 0)THEN
           WRITE(6,*) ' Error with _seupd, ierr = ', ierr
           WRITE(6,*) ' Check the documentation of _seupd. '
           STOP
        ENDIF

        resid_restart(:) = (0.d0,0.d0)

        DO jj= 1,nev
           resid_restart(:) = resid_restart(:)+vtmp(:,jj)
        ENDDO
        resid_restart(:) = resid_restart(:) / REAL(nev,DP)

     ENDIF
     !output value
     !
     eval(:)=REAL(d(1:nev,1),DP)
     evec(:,:)=vtmp(:,1:nev)
     !
     !next
     nec=iparam(5)
     maxmvs=iparam(9)
     !------
     IF(ALLOCATED(vtmp))   DEALLOCATE(vtmp)
     IF(ALLOCATED(resid))  DEALLOCATE(resid)
     IF(ALLOCATED(workd))  DEALLOCATE(workd)
     IF(ALLOCATED(workl))  DEALLOCATE(workl)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE real_diagM_arpk
  !Interface
  SUBROUTINE rdiagM_arpk(mat,dimen,nev,evec,eval)
     USE math, ONLY : rsort_eigen
     IMPLICIT NONE
     !IN/OUT
     INTEGER(I4B),INTENT(IN) :: dimen & !dimension of mat
                      & ,nev !number of eigen-pairs we need
     REAL(DP),INTENT(IN)  :: mat(:,:) !matrix
     REAL(DP),INTENT(OUT) :: evec(:,:),eval(:) !eigenpairs
     !LOCAL
     INTEGER(I4B) :: nec,isc,info
     REAL(DP)     :: resid_int(dimen),diagTOL=1.0D-6
     INTEGER(I4B) :: maxmvs
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     !initialize arpk
     info=0
     resid_int(:)=0.d0
     diagH:DO isc=1,5
               maxmvs=isc*30000
               CALL real_diagM_arpk(mat,nev,evec,eval, &
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
        WRITE(6,*) 'diag(M) failed,STOP,increase NADST and try again'
        STOP
     ENDIF

     !IF(nev>1)THEN
     !    CALL rsort_eigen(nev,eval,evec)
     !ENDIF
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE rdiagM_arpk
  !-----------------------DIVIDER-LINE--------------------------
  !###################################!
  !*For: Arpack for real Hamiltonian  !
  !*Date:2018-4-10                    !
  !*CALL dsaupd                       !
  !*CALL dseupd                       !
  !###################################!
  SUBROUTINE ISO_diagH_arpack(veff,nev,evec,eval,resid_restart,nec,info,maxmvs,TOL)
     !
     USE matvec_module , ONLY : ISO_rmatvec,ISO_sphere2grid
     IMPLICIT NONE
     ! INCLUDE 'debug.h'
     !  integer  logfil, ndigit, mgetv0,                                     &
     ! &         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,     &
     ! &         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,     &
     ! &         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
     !  common /debug/                                                       &
     ! &         logfil, ndigit, mgetv0,                                     &
     ! &         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,     &
     ! &         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,     &
     ! &         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
     !
     REAL(DP),INTENT(IN)  :: veff(:)   !IN:veff
     INTEGER(I4B),INTENT(IN)  :: nev      !number of state we need
     REAL(DP),INTENT(OUT) :: evec(:,:) !OUT:engivector
     REAL(DP),INTENT(OUT) :: eval(nev)     !engivalue
     REAL(DP),INTENT(INOUT) :: resid_restart(:)
     INTEGER(I4B),INTENT(INOUT) :: nec     !number of convergence
     INTEGER(I4B),INTENT(INOUT) :: info
     INTEGER(I4B),INTENT(INOUT) :: maxmvs
     REAL(DP),INTENT(IN) :: TOL
     !Arpack
     INTEGER(I4B)     :: ido,ncv,lworkl,mode1,         &
                    &    ishfts,maxitr,dimen,ldv,ierr
     INTEGER(I4B) :: iparam(11),ipntr(11),jj
     REAL(DP) :: d(2*nev+5,2)
     REAL(DP),ALLOCATABLE :: vtmp(:,:)
     REAL(DP),ALLOCATABLE :: resid(:),workd(:),workl(:)
     REAL(DP) :: sigma
     REAL(DP) :: veff_grid(n1,n2,n3)
     !
     LOGICAL :: rvec,select(2*nev+5)
     !
     CHARACTER(LEN=1) :: bmat
     CHARACTER(LEN=2) :: which
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     !>>>
     dimen=SIZE(veff,1)
     !dimen=n
     ncv=2*nev+5
!     ncv=4
     ldv=dimen
     IF( dimen > maxn )THEN
        WRITE(6,*) 'ARPACK:ERROR with _NSIMP: N is greater than MAXN',dimen
        STOP
     ELSEIF ( nev > maxnev ) then
        WRITE(6,*) 'ARPACK:ERROR with _NSIMP: NEV is greater than MAXNEV',nev
        STOP
     ELSEIF ( ncv > maxncv ) then
        WRITE(6,*) 'ARPACK:ERROR with _NSIMP: NCV is greater than MAXNCV',ncv
        STOP
     ENDIF
     !<<<

     !debug variables>>>>
     !ndigit = -3
     !logfil = 6
     !mcaitr = 0
     !mcapps = 0
     !mcaupd = 1
     !mcaup2 = 0
     !mceigh = 0
     !mceupd = 0
     !<<<<debug variables
     !Initialize the arpack
     iparam(:)=0
     ipntr(:)=0
!!!

     bmat='I'
     which='SA'  !'LA'/'SA' for largest/smallest real eigenvaule

     lworkl=ncv*(ncv+8)

     ido=0
     mode1=1
     ishfts=1
     maxitr=nev*50
!maxitr=100
     iparam(1)=ishfts
     iparam(3)=maxitr
     iparam(4)=1
     iparam(7)=mode1

     !allocate>>>
     ALLOCATE(vtmp(dimen,ncv))
     ALLOCATE(resid(dimen))
     ALLOCATE(workd(3*dimen))
     ALLOCATE(workl(lworkl))
     vtmp(:,:)=(0.d0,0.d0)
     !<<<allocate
     resid(:)=resid_restart(:)
     !Main loop>>>(Reverse Communication)
!
10   CONTINUE
        CALL dsaupd (ido, bmat, dimen, which, nev, tol, resid &
          &  ,ncv, vtmp,ldv, iparam, ipntr, workd, workl,  &
          &  lworkl, info)
        !need to be change--------------------->
        IF(ido==1.OR.ido==-1)THEN
           CALL ISO_sphere2grid(veff,veff_grid)
           CALL ISO_rmatvec(veff_grid,workd(ipntr(1)),workd(ipntr(2)),dimen)
           IF(iparam(9)<maxmvs) GOTO 10
        ENDIF
     !<<<Main loop
     !debug
     IF ( info == 1) THEN
        WRITE(6,'(/,a,/)') 'ARPACK:Maximum number of iterations reached.'
        !STOP
     ELSEIF ( info .eq. 3) THEN
        WRITE(6,'(/,a,a,/)') 'ARPACK:No shifts could be applied during',  &
              &  'implicit Arnoldi update, try increasing NCV.'
        STOP
     ENDIF

     IF ( info .lt. 0 ) THEN
        WRITE(6,'(/,a,i5)') 'ARPACK:Error with _saupd, info = ', info
        WRITE(6,'(a,/)') 'ARPACK:Check documentation in _saupd '
        STOP
     ELSE
        rvec=.TRUE.
        CALL dseupd (rvec, 'A',select, d, vtmp, ldv,    &
       &       sigma,bmat, dimen, which, nev, tol, resid, ncv,  &
       &       vtmp, ldv,iparam, ipntr, workd, workl, lworkl  &
       &       , ierr)

        IF( ierr /= 0)THEN
           WRITE(6,*) ' Error with _seupd, ierr = ', ierr
           WRITE(6,*) ' Check the documentation of _seupd. '
           STOP
        ENDIF

        resid_restart(:) = (0.d0,0.d0)

        DO jj= 1,nev
           resid_restart(:) = resid_restart(:)+vtmp(:,jj)
        ENDDO
        resid_restart(:) = resid_restart(:) / REAL(nev,DP)

        !     Print additional convergence information
        !        print *, '-------------------ARPACK--------------------'
        !        print *, '_NSIMP '
        !        print *, '====== '
        !        print *, ' '
        !        print *, ' Size of the matrix is ', dimen
        !        print *, ' The number of Ritz values requested is ', nev
        !        print *, ' The number of Arnoldi vectors generated',  &
        !    &            ' (NCV) is ', ncv
        !        print *, ' What portion of the spectrum: ', which
        !        print *, ' The number of converged Ritz values is ',  &
        !    &              nec
        !        print *, ' The number of Implicit Arnoldi update',    &
        !    &            ' iterations taken is ', iparam(3)
        !        print *, ' The number of OP*x is ', iparam(9)
        !        print *, ' The convergence criterion is ', tol
        !        print *, '---------------------------------------------'


     ENDIF
     !output value
     !
     eval(:)=REAL(d(1:nev,1),DP)
     evec(:,:)=vtmp(:,1:nev)
     !
     !next
     nec=iparam(5)
     maxmvs=iparam(9)
     !------
     IF(ALLOCATED(vtmp))   DEALLOCATE(vtmp)
     IF(ALLOCATED(resid))  DEALLOCATE(resid)
     IF(ALLOCATED(workd))  DEALLOCATE(workd)
     IF(ALLOCATED(workl))  DEALLOCATE(workl)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ENDSUBROUTINE ISO_diagH_arpack
  !-----------------------PARTING-LINE--------------------------
#ifdef MPI
  SUBROUTINE diagH_arpack_band(veff,Ik,nev,evec,eval,resid_restart,nec,info,maxmvs,TOL)
    USE matvec_module , ONLY : cmatvec_band
    USE grid_module, ONLY: global_n
! #ifdef MPI
!      USE smpi_math_module, ONLY:parallel
! #endif
     !
     IMPLICIT NONE

     ! INCLUDE 'debug.h'
     !  integer  logfil, ndigit, mgetv0,                                     &
     ! &         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,     &
     ! &         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,     &
     ! &         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
     !  common /debug/                                                       &
     ! &         logfil, ndigit, mgetv0,                                     &
     ! &         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,     &
     ! &         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,     &
     ! &         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd

     REAL(DP),INTENT(IN)  :: veff(:,:,:)   !IN:veff
     INTEGER(I4B),INTENT(IN) :: Ik
     INTEGER(I4B),INTENT(IN)  :: nev      !number of state we need
     REAL(DP),INTENT(IN) :: TOL
     COMPLEX(DP),INTENT(OUT) :: evec(:,:) !OUT:engivector
     REAL(DP),INTENT(OUT) :: eval(nev)     !engivalue
     COMPLEX(DP),INTENT(INOUT) :: resid_restart(:)
     INTEGER(I4B),INTENT(INOUT) :: nec     !number of convergence
     INTEGER(I4B),INTENT(INOUT) :: info
     INTEGER(I4B),INTENT(INOUT) :: maxmvs
     !Arpack
     INTEGER(I4B)     :: ido,ncv,lworkl,mode1,         &
                    &    ishfts,maxitr,dimen,ldv,ierr
     INTEGER(I4B) :: iparam(11),ipntr(14),jj
     COMPLEX(DP) :: d(2*nev+5,2)
     COMPLEX(DP),ALLOCATABLE :: vtmp(:,:)
     COMPLEX(DP),ALLOCATABLE :: resid(:),workd(:),workl(:),workev(:)
     REAL(DP),ALLOCATABLE :: rwork(:)
     COMPLEX(DP) :: sigma
     !
     LOGICAL :: rvec,select(2*nev+5)
     !
     CHARACTER(LEN=1) :: bmat
     CHARACTER(LEN=2) :: which
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     !>>>
     !dimen=SIZE(,1)
     dimen=global_n
     ncv=2*nev+5
!     ncv=4
     ldv=dimen
     IF( dimen > maxn )THEN
        WRITE(6,*) 'ARPACK:ERROR with _NSIMP: N is greater than MAXN',dimen
        STOP
     ELSEIF ( nev > maxnev ) then
        WRITE(6,*) 'ARPACK:ERROR with _NSIMP: NEV is greater than MAXNEV',nev
        STOP
     ELSEIF ( ncv > maxncv ) then
        WRITE(6,*) 'ARPACK:ERROR with _NSIMP: NCV is greater than MAXNCV',ncv
        STOP
     ENDIF
     !<<<

     !debug variables>>>>
     ! ndigit = -3
     ! logfil = 6
     ! mcaitr = 0
     ! mcapps = 0
     ! mcaupd = 1
     ! mcaup2 = 0
     ! mceigh = 0
     ! mceupd = 0
     !<<<<debug variables
     !Initialize the arpack
     iparam(:)=0
     ipntr(:)=0
!!!

     bmat='I'
     which='SR'  !'LR'/'SR' for largest/smallest real eigenvaule

     lworkl=3*ncv**2+5*ncv

     ido=0
     mode1=1
     ishfts=1
     maxitr=nev*50
!maxitr=100
     iparam(1)=ishfts
     iparam(3)=maxitr
     iparam(7)=mode1

     !allocate>>>
     ALLOCATE(vtmp(dimen,2*nev+5))
     ALLOCATE(resid(dimen))
     ALLOCATE(workd(3*dimen))
     ALLOCATE(workl(lworkl))
     ALLOCATE(rwork(ncv))
     ALLOCATE(workev(3*ncv))
     vtmp(:,:)=(0.d0,0.d0)
     !<<<allocate
     resid(:)=resid_restart(:)
     !Main loop>>>(Reverse Communication)
!check
!
10   CONTINUE
! #ifdef MPI
!         CALL pznaupd (parallel%commx,ido, bmat, dimen, which, nev, tol, resid &
!           &  ,ncv, vtmp,ldv, iparam, ipntr, workd, workl,  &
!           &  lworkl,rwork, info)
! #else
        CALL znaupd (ido, bmat, dimen, which, nev, tol, resid &
          &  ,ncv, vtmp,ldv, iparam, ipntr, workd, workl,  &
          &  lworkl,rwork, info)
! #endif
        !need to be change--------------------->
!print*,'ido',ido
        IF(ido==1.OR.ido==-1)THEN
           !CALL matvec (mat,workd(ipntr(1):ipntr(1)+dimen-1),workd(ipntr(2):ipntr(2)+dimen-1),1,dimen)
           CALL cmatvec_band(veff,Ik,workd(ipntr(1)),workd(ipntr(2)),dimen)
           IF(iparam(9)<maxmvs) GOTO 10
        ENDIF
     !<<<Main loop
     !debug
     IF ( info == 1) THEN
        WRITE(6,'(/,a,/)') 'ARPACK:Maximum number of iterations reached.'
        !STOP
     ELSEIF ( info .eq. 3) THEN
        WRITE(6,'(/,a,a,/)') 'ARPACK:No shifts could be applied during',  &
              &  'implicit Arnoldi update, try increasing NCV.'
        STOP
     ENDIF

     IF ( info .lt. 0 ) THEN
! #ifdef MPI
!            if(parallel%isroot)then
! #endif
        WRITE(6,'(/,a,i5)') 'ARPACK:Error with _saupd, info = ', info
        WRITE(6,'(a,/)') 'ARPACK:Check documentation in _saupd '
! #ifdef MPI
!            endif
! #endif
        STOP
     ELSE
        rvec=.TRUE.
! #ifdef MPI
!         CALL pzneupd (parallel%commx,rvec, 'A',select, d, vtmp, ldv,    &
!        &       sigma,workev,bmat, dimen, which, nev, tol, resid, ncv,  &
!        &       vtmp, ldv,iparam, ipntr, workd, workl, lworkl  &
!        &       ,rwork, ierr)
! #else
        CALL zneupd (rvec, 'A',select, d, vtmp, ldv,    &
       &       sigma,workev,bmat, dimen, which, nev, tol, resid, ncv,  &
       &       vtmp, ldv,iparam, ipntr, workd, workl, lworkl  &
       &       ,rwork, ierr)
! #endif

        IF( ierr /= 0)THEN
! #ifdef MPI
!            if(parallel%isroot)then
! #endif
           WRITE(6,*) ' Error with _seupd, ierr = ', ierr
           WRITE(6,*) ' Check the documentation of _seupd. '
! #ifdef MPI
!            endif
! #endif
           STOP
        ENDIF

        resid_restart(:) = (0.d0,0.d0)

        DO jj= 1,nev
           ! CALL cmatvec(veff,Ik,vtmp(:,jj),resid(:),dimen)
           ! resid_restart(:) = resid_restart+resid-d(jj,1)*vtmp(:,jj)
           resid_restart(:) = resid_restart(:)+vtmp(:,jj)
        ENDDO
        resid_restart(:) = resid_restart(:) / REAL(nev,DP)

!     Print additional convergence information
   !     print *, '-------------------ARPACK--------------------'
   !     print *, '_NSIMP '
   !     print *, '====== '
   !     print *, ' '
   !     print *, ' Size of the matrix is ', dimen
   !     print *, ' The number of Ritz values requested is ', nev
   !     print *, ' The number of Arnoldi vectors generated',  &
   ! &            ' (NCV) is ', ncv
   !     print *, ' What portion of the spectrum: ', which
   !     print *, ' The number of converged Ritz values is ',  &
   ! &              nec
   !     print *, ' The number of Implicit Arnoldi update',    &
   ! &            ' iterations taken is ', iparam(3)
   !     print *, ' The number of OP*x is ', iparam(9)
   !     print *, ' The convergence criterion is ', tol
   !     print *, '---------------------------------------------'


     ENDIF
     !output value
     !
     eval(:)=REAL(d(1:nev,1),8)
     evec(:,:)=vtmp(:,1:nev)
     !
     !next
     nec=iparam(5)
     maxmvs=iparam(9)
     !------
     IF(ALLOCATED(vtmp))   DEALLOCATE(vtmp)
     IF(ALLOCATED(resid))  DEALLOCATE(resid)
     IF(ALLOCATED(workd))  DEALLOCATE(workd)
     IF(ALLOCATED(workl))  DEALLOCATE(workl)
     IF(ALLOCATED(rwork))  DEALLOCATE(rwork)
     IF(ALLOCATED(workev)) DEALLOCATE(workev)
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   ENDSUBROUTINE diagH_arpack_band
#endif
!==========================================================================

ENDMODULE Arpack_module

MODULE Arpack_module
!###########################################################!
!*For: Arpack                                               !
!*Author:Qiang Xu                                           !
!*Date:2017-7-20                                            !
!###########################################################!
    USE constants
    IMPLICIT NONE
    INTEGER(I4B),PARAMETER :: maxn=1000000, maxnev=500, maxncv=10000
!c
!c\SCCS Information: @(#)
!c FILE: debug.h   SID: 2.3   DATE OF SID: 11/16/95   RELEASE: 2
!c
!c     %---------------------------------%
!c     | See debug.doc for documentation |
!c     %---------------------------------%
    integer logfil, ndigit, mgetv0,                                     &
   &         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,     &
   &         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,     &
   &         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
    common/debug/                                                       &
   &         logfil, ndigit, mgetv0,                                     &
   &         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,     &
   &         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,     &
   &         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
CONTAINS
!----------------------------------------------------------
!###################################!
!*For: Arpack for real Hamiltonian  !
!*Author:Qiang Xu                   !
!*Date:2017-8-28                    !
!*CALL dsaupd                       !
!*CALL dseupd                       !
!###################################!
    SUBROUTINE real_diagH_arpack(n, veff, nev, evec, eval, resid_restart, nec, info, maxmvs, TOL)
!
        USE matvec_module, ONLY: real_matvec
        IMPLICIT NONE
!INCLUDE 'debug.h'
!
        INTEGER(I4B), INTENT(IN)  :: n, nev      !number of points/state we need
        REAL(DP), INTENT(IN)  :: veff(n)   !IN:veff
        REAL(DP), INTENT(OUT) :: evec(:, :) !OUT:engivector
        REAL(DP), INTENT(OUT) :: eval(nev)     !engivalue
        REAL(DP), INTENT(INOUT) :: resid_restart(:)
        INTEGER(I4B), INTENT(INOUT) :: nec     !number of convergence
        INTEGER(I4B), INTENT(INOUT) :: info
        INTEGER(I4B), INTENT(INOUT) :: maxmvs
        REAL(DP), INTENT(IN) :: TOL
!Arpack
        INTEGER(I4B)     :: ido, ncv, lworkl, mode1,         &
                       &    ishfts, maxitr, dimen, ldv, ierr
        INTEGER(I4B) :: iparam(11), ipntr(11), jj
        REAL(DP) :: d(2*nev + 5, 2)
        REAL(DP), ALLOCATABLE :: vtmp(:, :)
        REAL(DP), ALLOCATABLE :: resid(:), workd(:), workl(:)
        REAL(DP) :: sigma
!
        LOGICAL :: rvec, select(2*nev + 5)
!
        CHARACTER(LEN=1) :: bmat
        CHARACTER(LEN=2) :: which
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>
!dimen=SIZE(,1)
        dimen = n
        ncv = 2*nev + 5
!     ncv=4
        ldv = dimen
        IF (dimen > maxn) THEN
            WRITE (*, *) 'ARPACK:ERROR with _NSIMP: N is greater than MAXN', dimen
            STOP
        ELSEIF (nev > maxnev) then
            WRITE (*, *) 'ARPACK:ERROR with _NSIMP: NEV is greater than MAXNEV', nev
            STOP
        ELSEIF (ncv > maxncv) then
            WRITE (*, *) 'ARPACK:ERROR with _NSIMP: NCV is greater than MAXNCV', ncv
            STOP
        END IF
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
        iparam(:) = 0
        ipntr(:) = 0
!!!

        bmat = 'I'
        which = 'SA'  !'LA'/'SA' for largest/smallest real eigenvaule

        lworkl = ncv*(ncv + 8)

        ido = 0
        mode1 = 1
        ishfts = 1
        maxitr = nev*50
!maxitr=100
        iparam(1) = ishfts
        iparam(3) = maxitr
        iparam(4) = 1
        iparam(7) = mode1

!allocate>>>
        ALLOCATE (vtmp(dimen, ncv))
        ALLOCATE (resid(dimen))
        ALLOCATE (workd(3*dimen))
        ALLOCATE (workl(lworkl))
        vtmp(:, :) = (0.d0, 0.d0)
!<<<allocate
        resid(:) = resid_restart(:)
!Main loop>>>(Reverse Communication)
!
10      CONTINUE
        CALL dsaupd(ido, bmat, dimen, which, nev, tol, resid &
          &  , ncv, vtmp, ldv, iparam, ipntr, workd, workl,  &
          &  lworkl, info)
!need to be change--------------------->
        IF (ido == 1 .OR. ido == -1) THEN
            CALL real_matvec(veff, workd(ipntr(1)), workd(ipntr(2)), dimen)
            IF (iparam(9) < maxmvs) GOTO 10
        END IF
!<<<Main loop
!debug
        IF (info == 1) THEN
            WRITE (*, '(/,a,/)') 'ARPACK:Maximum number of iterations reached.'
!STOP
        ELSEIF (info .eq. 3) THEN
            write (*, '(/,a,a,/)') 'ARPACK:No shifts could be applied during',  &
                  &  'implicit Arnoldi update, try increasing NCV.'
            STOP
        END IF

        IF (info .lt. 0) THEN
            write (*, '(/,a,i5)') 'ARPACK:Error with _saupd, info = ', info
            write (*, '(a,/)') 'ARPACK:Check documentation in _saupd '
            STOP
        ELSE
            rvec = .TRUE.
            CALL dseupd(rvec, 'A', select, d, vtmp, ldv,    &
           &       sigma, bmat, dimen, which, nev, tol, resid, ncv,  &
           &       vtmp, ldv, iparam, ipntr, workd, workl, lworkl  &
           &       , ierr)

            IF (ierr /= 0) THEN
                write (*, *) ' Error with _seupd, ierr = ', ierr
                write (*, *) ' Check the documentation of _seupd. '
                STOP
            END IF

            resid_restart(:) = (0.d0, 0.d0)

            DO jj = 1, nev
                resid_restart(:) = resid_restart(:) + vtmp(:, jj)
            END DO
            resid_restart(:) = resid_restart(:)/REAL(nev, DP)

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

        END IF
!output value
!
        eval(:) = REAL(d(1:nev, 1), DP)
        evec(:, :) = vtmp(:, 1:nev)
!
!next
        nec = iparam(5)
        maxmvs = iparam(9)
!------
        IF (ALLOCATED(vtmp)) DEALLOCATE (vtmp)
        IF (ALLOCATED(resid)) DEALLOCATE (resid)
        IF (ALLOCATED(workd)) DEALLOCATE (workd)
        IF (ALLOCATED(workl)) DEALLOCATE (workl)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    END SUBROUTINE real_diagH_arpack
!----------------------------------------------------------

END MODULE Arpack_module

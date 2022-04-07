!#############################################################!
!        Chebyshev filter Method for diag(H)                  !
!*Author : Qiang Xu                                           !
!*Date   : 2017/11/28                                         !
!#############################################################!
MODULE chebyshev_module
    USE constants
    USE parameters, ONLY: CheM

    USE smpi_math_module
    USE parameters, ONLY: BLOCK_MBNB

    IMPLICIT NONE
    REAL(DP), PARAMETER :: LARGED = 1e4
    REAL(DP), ALLOCATABLE :: ad(:)
CONTAINS
!##############################################################!
!*First :            Pseudo Subspace Method                    !
!##############################################################!
    SUBROUTINE BuildSubspace(nps, nev, veff, eig)
        USE parameters, ONLY: nspin, CheM0, nwf0, Igamma &
                 & , nev_tot => Nstates_global
        USE potential_module, ONLY: calVeff
        USE pspot_module, ONLY: psp, max_nwfa
        USE grid_module, ONLY: eigen_type, nk
        USE struct_module, ONLY: naty, natom, struct
        IMPLICIT NONE
        !INOUT
        INTEGER(I4B), INTENT(IN) :: nps, nev !points and states
        REAL(DP), INTENT(IN) :: veff(nps, nspin)
        TYPE(eigen_type), INTENT(INOUT) :: eig
        !LOCAL
        INTEGER(I4B) :: Is, Ik, Ii, Ip  &
              &, Ity, Ib, nlm, Ia
        REAL(DP) :: randt, radius = 1._DP
        INTEGER(I4B) :: Nwfa  &! total orbital-like states we have
               &, Nrand, irand1  !# of random states in subspace
        REAL(DP) :: X0(nps, nev)

        INTEGER(I4B) :: Isl

        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        IF (max_nwfa > 0) THEN !revised later
            Nwfa = 0
            Nrand = 0
            DO Ity = 1, naty
                DO Ia = struct%eleid(Ity), struct%eleid(Ity + 1) - 1
                    DO Ib = 1, psp(Ity)%nwfa
                        nlm = 2*psp(Ity)%wfal(Ib) + 1
                        Nwfa = Nwfa + nlm
                    END DO
                END DO
            END DO

            IF (parallel%isroot) THEN

                PRINT *, 'CheFSI: Dim. of Pseudo Subspace:', Nwfa

            END IF

            !setting

            Nrand = nev_tot - Nwfa

            irand1 = Nwfa + 1
            !setting the atomic wvfs
            CALL real_PseudoSubspace(nps, nev, X0)
        ELSE
            Nrand = nev
            irand1 = 1
        END IF
        !random states
        IF (Nrand > 0) THEN

            IF (parallel%isroot) THEN

                PRINT *, 'CheFSI: # of random samping states:', Nrand

            END IF
            ! at least proceed here
            !initialize the temp subspace
            CALL random_seed()
            write (*, *) 'in Cheby ,1'
            DO Ii = irand1, nev_tot !for global states
                write (*, *) 'in Cheby ,2'
                DO Isl = 1, nev, 1  !for local states
                    write (*, *) 'in Cheby 3',Isl
                    IF (Ii == parallel%sub2sum(Isl, parallel%ranky + 1)) THEN
                        DO Ip = 1, nps
                            write (*, *) 'in Cheby ', Ip, 'of', shape(X0)
                            CALL random_number(randt)
                            randt = -radius + randt*2.d0*radius
                            X0(Ip, Isl) = randt
                        END DO
                    END IF
                END DO
            END DO
        END IF
!Raleigh-Ritz step
!Gamma point
        DO Is = 1, Nspin !spiner

            DO Ik = 1, nk !k-points

                IF (Ik /= IGamma) THEN
!for non-Gamma k-points
                    STOP 'Waitting for non-Gamma'
                ELSE
!for Gamma k-points
                    eig%wvfG(:, :, Is) = X0(:, :)
                    IF (CheM0 > 0) THEN
!filter
                        CALL real_first_filter(nps, nev, veff(:, Is), &
                          & eig%wvfG(:, :, Is), eig%val(:, IGamma, Is))
                    ELSE
!RR
                        CALL real_first_RRstep(nps, nev, veff(:, Is), &
                           &   eig%wvfG(:, :, Is), eig%val(:, Ik, Is))
                    END IF
                END IF
            END DO
        END DO
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    END SUBROUTINE BuildSubspace
!----------------------------------------------------------
    SUBROUTINE real_PseudoSubspace(nps, nev, initX)
        USE parameters, ONLY: nev_tot => Nstates_global
        USE grid_module, ONLY: n1, n2, n3, grid
        USE math, ONLY: interp
        USE struct_module, ONLY: struct, naty, lat_mat
        USE pspot_module, ONLY: psp
        USE nlpot_module, ONLY: apply_ylm
        IMPLICIT NONE
        INTEGER(I4B), INTENT(IN) :: nps, nev
        REAL(DP), INTENT(OUT) :: initX(nps, nev)
!
        INTEGER(I4B) :: Ity, Ia, Ib, Ip, Il, Im, Ii, icx, icy, icz
        REAL(DP) :: ra(3), rat(4), wfa, x, y, z
        REAL(DP) :: fac

        INTEGER(I4B) :: Isl, Ithis

        INTEGER(I4B) :: nrep = 1
        REAL(DP) :: wrcut = 4*ang2bohr
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        initX = 0._DP
        fac = 1.0_DP
        Ii = 0 !count global states

        Ithis = 0 !count this states

        DO Ity = 1, naty !all type

            IF (psp(Ity)%nwfa == 0) CYCLE

            DO Ia = struct%eleid(Ity), struct%eleid(Ity + 1) - 1 !all atoms
!store positions
                ra(:) = struct%poscar(:, Ia)
!
                DO Ib = 1, psp(Ity)%nwfa
                    Il = psp(Ity)%wfal(Ib)
                    DO Im = -Il, Il

                        Ii = Ii + 1
                        IF (Ii > nev_tot) Ii = 1

                        DO Isl = 1, nev, 1

                            IF (Ii == parallel%sub2sum(Isl, parallel%ranky + 1)) THEN
                                Ithis = Ithis + 1
                                IF (Ithis .gt. nev) Ithis = 1

!call states
                                DO ip = 1, nps
!replicas
                                    DO icz = -nrep, nrep
                                    DO icy = -nrep, nrep
                                    DO icx = -nrep, nrep

                                        rat(1:3) = icx*lat_mat(:, 1)  &
                                          & + icy*lat_mat(:, 2)  &
                                          & + icz*lat_mat(:, 3)  &
                                          & + grid%rVec(1:3, Ip)  &
                                          & - ra(:)

                                        rat(4) = SQRT(DOT_PRODUCT(rat(1:3), rat(1:3)))
!for efficient
                                        IF (rat(4) > wrcut) CYCLE

                                        IF (rat(4) > 0.d0) THEN
                                            x = rat(1)/rat(4)
                                            y = rat(2)/rat(4)
                                            z = rat(3)/rat(4)
                                        ELSE
                                            x = 0.d0
                                            y = 0.d0
                                            z = 0.d0
                                        END IF

                                        wfa = interp(psp(Ity)%numps, psp(Ity)%wfar(:, Ib), psp(Ity)%r, rat(4))
                                        CALL apply_ylm(Il, Im, fac, x, y, z, wfa)
!initilize the subspace

                                        initX(ip, Ithis) = initX(ip, Ithis) + wfa

                                    END DO
                                    END DO
                                    END DO
                                END DO
!

                            END IF
                        END DO !> Is

                    END DO
                END DO
            END DO
        END DO
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    END SUBROUTINE real_PseudoSubspace
!----------------------------------------------------------
    SUBROUTINE real_first_RRstep(nps, nev, veff, X, D)
        USE Lapack_module, ONLY: matmat_real, OrthNorm_real, diagM_real&
                         &, GeneralizeEigen_real

        USE grid_module, ONLY: global_n
        USE ScaLapack_module, ONLY: SL_GeneralizeEigen_real &
               &, SL_matmat_real_tn, SL_matmat_real_nn&
               &, SL_OrthNorm_real, SL_diagM_real, twoD_map
        USE parameters, ONLY: CheM0, LRROrthNorm &
                        &, nev_tot => Nstates_global, BLOCK_MBNB

        IMPLICIT NONE
!IN/OUT
        INTEGER(I4B), INTENT(IN) :: nps, nev
        REAL(DP), INTENT(IN) :: veff(nps)
        REAL(DP), INTENT(INOUT) :: X(nps, nev) !subspace
        REAL(DP), INTENT(OUT) :: D(:)      !eigenvalue
!LOCAL
!REAL(DP),DIMENSION(nev,nev) :: Shat,Hhat,Qs
        REAL(DP) :: Xnew(nps, nev)
        REAL(DP) :: a, b, al

        INTEGER(I4B) ::   m   & !dim. 1
                     &, n   & !dim. 2
                     &, mb  & !block of dim. 1
                     &, nb  &  !block of dim. 2
                     &, cm  &  !dim.1 of X^T*H*X
                     &, cn  &  !dim.2 of X^T*H*X
                     &, cmb & !block of dim.1 of X^T*H*X
                     &, cnb & !block of dim.2 of X^T*H*X
                     &, nstsp !number of sts of this proc
        REAL(DP), DIMENSION(twoD_map(1, parallel%rankx,&
             &parallel%ranky), nev) :: Hhat, Qs, Shat

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        m = global_n  !dim. 1
        n = nev_tot   !dim. 2
        mb = 1  !block of dim. 1
        nb = 1  !block of dim. 2
        cm = nev_tot  !dim.1 of X^T*H*X
        cn = nev_tot   !dim.2 of X^T*H*X
        cmb = 1  !block of dim.1 of X^T*H*X
        cnb = 1  !block of dim.2 of X^T*H*X
        nstsp = parallel%nstate_proc

!Raleigh-Ritz step
        IF (LRROrthNorm) THEN
!OrthNorm

            CALL SL_OrthNorm_real(X, m, n, mb, nb)

!RR
            CALL Rayleigh_quotient_real(nps, nev, veff, X, Hhat)
!eigen-decomposion

!STOP 'SL_diagM_real is needed'
            CALL SL_diagM_real(n, Hhat, n, n, n, n, cmb, cnb, cmb, cnb, Qs, n, n, D, cmb, cnb)

        ELSE
!Overlap matrix

            CALL SL_matmat_real_tn('T', 'N', X, X, Shat, m, n, m, n, mb, nb, mb, nb, cmb, cnb)

!projected hamiltonian
            CALL Rayleigh_quotient_real(nps, nev, veff, X, Hhat)

!eigen-decomposion

            CALL SL_GeneralizeEigen_real(n, Hhat, Shat, n, n, n, n  &
                &, cmb, cnb, cmb, cnb, Qs, n, n, D, cmb, cnb)

        END IF
!-------------------
!rotation

        CALL SL_matmat_real_nn('N', 'N', X, Qs, Xnew, m, n, n, n, mb, nb, cmb, cmb, cmb, cnb)

!eigen-value
        X(:, :) = Xnew(:, :)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    END SUBROUTINE real_first_RRstep
!----------------------------------------------------------
    SUBROUTINE init_uplow_real(nps, k, veff, v, a, b, al)
!
        USE matvec_module, ONLY: real_matvec
        USE Lapack_module, ONLY: diagM_real
        IMPLICIT NONE
!INOUT
        INTEGER(I4B), INTENT(IN) :: nps, k
!hamiltonian
        REAL(DP), INTENT(IN) :: veff(nps)
        REAL(DP), INTENT(INOUT) :: v(nps)
        REAL(DP), INTENT(OUT) :: a, b, al
!LOCAL
        REAL(DP) :: v0(nps), f(nps), alpha
        REAL(DP) :: beta, &
                &   fbeta = 0.5d0
        REAL(DP) :: T(k, k), evec(k, k), eval(k)
        INTEGER(I4B) :: J, Nj

        REAL(DP) :: alp_loc, bet_loc

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        T(:, :) = 0.d0
        Nj = MIN(k, 10)
!
        CALL real_matvec(veff, v, f, nps)

        alp_loc = DOT_PRODUCT(f, v)
        CALL MPI_ALLREDUCE(alp_loc, alpha, 1, MPI_REAL8, MPI_SUM, parallel%commx, mpinfo)

        f = f - alpha*v
        T(1, 1) = alpha
        DO J = 2, Nj
!beta=SQRT(DOT_PRODUCT(f,f))

            bet_loc = REAL(DOT_PRODUCT(f, f), 8)
            CALL MPI_ALLREDUCE(bet_loc, beta, 1, MPI_REAL8, MPI_SUM, parallel%commx, mpinfo)

            beta = SQRT(beta)

            v0 = v
            v = f/beta
            CALL real_matvec(veff, v, f, nps)
            f = f - beta*v0

            alp_loc = DOT_PRODUCT(f, v)
            CALL MPI_ALLREDUCE(alp_loc, alpha, 1, MPI_REAL8, MPI_SUM, parallel%commx, mpinfo)

            f = f - alpha*v
            T(J, J - 1) = beta
            T(J - 1, J) = beta
            T(J, J) = alpha
        END DO
!Rayleigh-Ritz value
!CALL diagM_real(T,evec,eval)

        IF (parallel%isroot) THEN

            CALL diagM_real(T, evec, eval)

        END IF !root
        CALL MPI_BCAST(evec, k*k, MPI_REAL8, parallel%rootid, parallel%comm, mpinfo)
        CALL MPI_BCAST(eval, k, MPI_REAL8, parallel%rootid, parallel%comm, mpinfo)
!resue beta
        bet_loc = REAL(DOT_PRODUCT(f, f), 8)
        CALL MPI_ALLREDUCE(bet_loc, beta, 1, MPI_REAL8, MPI_SUM, parallel%commx, mpinfo)

        a = fbeta*eval(1) + (1.d0 - fbeta)*eval(k)
        al = eval(1)
        b = eval(k) + SQRT(beta)*ABS(evec(k, k)) !+1e-10
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    END SUBROUTINE init_uplow_real
!----------------------------------------------------------
    SUBROUTINE real_first_filter(nps, nst, veff, X, eval)
        USE parameters, ONLY: CF0 => CheM0, LRROrthNorm, nev_tot => Nstates_global
        USE Lapack_module, ONLY: diagM_real, OrthNorm_real, matmat_real &
                         &, GeneralizeEigen_real

        USE grid_module, ONLY: nsp => global_n
        USE ScaLapack_module, ONLY: SL_GeneralizeEigen_real, twoD_map &
               &, SL_matmat_real_tn, SL_matmat_real_nn, SL_OrthNorm_real &
               &, SL_matmat_real, SL_diagM_real

        IMPLICIT NONE
        INTEGER(I4B), INTENT(IN) :: nps, nst !number of points/states
        REAL(DP), INTENT(IN) :: veff(nps)
        REAL(DP), INTENT(INOUT) :: X(:, :)
        REAL(DP), INTENT(OUT) :: eval(:)
!LOCAL
        REAL(DP) :: a, b, al, t
        INTEGER(I4B) :: I, Niter = 4
        REAL(DP) :: deval, TOL = 1e-8

        REAL(DP) :: evald(nev_tot)
        INTEGER(I4B) ::   m   & !dim. 1
                     &, n   & !dim. 2
                     &, mb  & !block of dim. 1
                     &, nb  &  !block of dim. 2
                     &, cm  &  !dim.1 of X^T*H*X
                     &, cn  &  !dim.2 of X^T*H*X
                     &, cmb & !block of dim.1 of X^T*H*X
                     &, cnb & !block of dim.2 of X^T*H*X
                     &, nstsp !number of sts of this proc
        REAL(DP), DIMENSION(twoD_map(1, parallel%rankx,&
             &parallel%ranky), nst) :: Shat, Hhat, Qs

        REAL(DP)  :: Xnew(nps, nst), vec(nps)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        m = nsp   !dim. 1
        n = nev_tot  !dim. 2
        mb = BLOCK_MBNB  !block of dim. 1
        nb = BLOCK_MBNB   !block of dim. 2
        cm = nev_tot   !dim.1 of X^T*H*X
        cn = nev_tot   !dim.2 of X^T*H*X
        cmb = BLOCK_MBNB  !block of dim.1 of X^T*H*X
        cnb = BLOCK_MBNB  !block of dim.2 of X^T*H*X
        nstsp = parallel%nstate_proc

!low up bound
        vec(:) = X(:, nst)
        CALL init_uplow_real(nps, 7, veff, vec, a, b, al)
!
        evald(:) = 612509.d0
        DO I = 1, Niter
!
            CALL chebyshev_filter_scaled_real(nps, nst, veff, X, CF0, a, b, al)
            IF (LRROrthNorm) THEN
!OrthNorm

                CALL SL_OrthNorm_real(X, m, n, mb, nb)

!xHx

                CALL Rayleigh_quotient_real(nps, nst, veff, X, Hhat)

!eigen-decomposion

                CALL SL_diagM_real(n, Hhat, n, n, n, n, cmb, cnb, cmb, cnb, Qs, n, n, eval, cmb, cnb)

            ELSE
!Overlap matrix
!CALL matmat_real(X,X,'T','N',Shat)

!CALL SL_matmat_real('T','N',X,X,Shat,m,n,m,n,mb,nb,mb,nb,cmb,cnb)
                CALL SL_matmat_real_tn('T', 'N', X, X, Shat, m, n, m, n, mb, nb, mb, nb, cmb, cnb)

!projected hamiltonian
!CALL Rayleigh_quotient_real(nps,nst,veff,X,Hhat)

                CALL Rayleigh_quotient_real(nps, nst, veff, X, Hhat)

!eigen-decomposion

                CALL SL_GeneralizeEigen_real(n, Hhat, Shat, n, n, n, n  &
                 &, cmb, cnb, cmb, cnb, Qs, n, n, eval, cmb, cnb)

            END IF

            deval = SUM(ABS(eval - evald))/n

!-----------------
            IF (deval < TOL) EXIT
!-----------------
!store old eigenvalue
            evald(:) = eval(:)
!update the new bound

            a = eval(n)

            al = eval(1)
        END DO

!rotation

        CALL SL_matmat_real_nn('N', 'N', X, Qs, Xnew, m, n, n, n, mb, nb, cmb, cmb, cmb, cnb)

        X = Xnew
!PRINT*,'Free diag>>>iter:',I
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    END SUBROUTINE real_first_filter
!##############################################################!
!*For    :     real Filtering (Standard RR)                    !
!*Author : Qiang Xu                                            !
!*Date   : 2018-03-05                                          !
!##############################################################!
!---------------------Rayleigh-quotient-------------------------
    SUBROUTINE Rayleigh_quotient_real(nps, nst, veff, x, xhx)
!xhx=(X,H,X)
        USE Lapack_module, ONLY: matmat_real

        USE ScaLapack_module, ONLY: SL_matmat_real_tn, twoD_map
        USE grid_module, ONLY: nsp => global_n
        USE parameters, ONLY: nev_tot => Nstates_global

        IMPLICIT NONE
        INTEGER(I4B), INTENT(IN) :: nps, nst
        REAL(DP), INTENT(IN) :: veff(nps)
        REAL(DP) :: x(nps, nst)

        REAL(DP) :: xhx(twoD_map(1, parallel%rankx,&
             &parallel%ranky), nst)

!LOCAL
        REAL(DP) :: hx(nps, nst)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        CALL cal_HX_real(nps, nst, veff, x, hx)
!xhx(:,:)=MATMUL( TRANSPOSE(CONJG(x)) , hx )

        CALL SL_matmat_real_tn('T', 'N', x, hx, xhx, nsp, nev_tot, nsp, nev_tot &
             &  , 1, BLOCK_MBNB, 1, BLOCK_MBNB, BLOCK_MBNB, BLOCK_MBNB)

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    END SUBROUTINE Rayleigh_quotient_real
!-------------------------------HX------------------------------
    SUBROUTINE cal_HX_real(nps, nst, veff, V, HV)
        USE matvec_module, ONLY: real_matvec
        IMPLICIT NONE
!INOUT
        INTEGER(I4B), INTENT(IN) :: nps, nst
        REAL(DP), INTENT(IN)  :: veff(nps) !effective potential
        REAL(DP), INTENT(IN)  :: V(nps, nst)
        REAL(DP), INTENT(OUT) :: HV(nps, nst)
!LOCAL
        INTEGER(I4B) :: Is
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        DO Is = 1, nst, 1
            CALL real_matvec(veff, V(:, Is), HV(:, Is), nps)
        END DO
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    END SUBROUTINE cal_HX_real
!-------------------------upper bound---------------------------
    SUBROUTINE Estupb_real(nps, k, veff, vec, b)
!
        USE matvec_module, ONLY: real_matvec
        USE Lapack_module, ONLY: diagM_real
        IMPLICIT NONE
!INOUT
        INTEGER(I4B), INTENT(IN) :: nps, k
!hamiltonian
        REAL(DP), INTENT(IN) :: veff(nps)
        REAL(DP), INTENT(IN) :: vec(nps)
        REAL(DP), INTENT(OUT) :: b
!LOCAL
        REAL(DP) :: v0(nps), v(nps), f(nps), alpha
        REAL(DP) :: beta, eval(k), mz
        REAL(DP) :: T(k, k), evec(k, k)
        INTEGER(I4B) :: J, Nj

        REAL(DP) :: alp_loc, bet_loc

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        v = vec !should normalize
        T(:, :) = 0.d0
        Nj = MIN(k, 10)
!
        CALL real_matvec(veff, v, f, nps)

        alp_loc = DOT_PRODUCT(f, v)
        CALL MPI_ALLREDUCE(alp_loc, alpha, 1, MPI_REAL8, MPI_SUM, parallel%commx, mpinfo)

        f = f - alpha*v
        T(1, 1) = alpha
        DO J = 2, Nj

            bet_loc = REAL(DOT_PRODUCT(f, f), 8)
            CALL MPI_ALLREDUCE(bet_loc, beta, 1, MPI_REAL8, MPI_SUM, parallel%commx, mpinfo)

            beta = SQRT(beta)
            v0 = v
            v = f/beta
            CALL real_matvec(veff, v, f, nps)
            f = f - beta*v0

            alp_loc = DOT_PRODUCT(f, v)
            CALL MPI_ALLREDUCE(alp_loc, alpha, 1, MPI_REAL8, MPI_SUM, parallel%commx, mpinfo)

            f = f - alpha*v
            T(J, J - 1) = beta
            T(J - 1, J) = beta
            T(J, J) = alpha
        END DO
!NORM2(T)
!b=Norm_2(T,k) + SQRT(REAL(DOT_PRODUCT(f,f),8))

        IF (parallel%isroot) THEN

            CALL diagM_real(T, evec, eval)

        END IF
        CALL MPI_BCAST(evec, k*k, MPI_REAL8, parallel%rootid, parallel%comm, mpinfo)
        CALL MPI_BCAST(eval, k, MPI_REAL8, parallel%rootid, parallel%comm, mpinfo)
!resue beta
        bet_loc = REAL(DOT_PRODUCT(f, f), 8)
        CALL MPI_ALLREDUCE(bet_loc, beta, 1, MPI_REAL8, MPI_SUM, parallel%commx, mpinfo)

!mz=MAX( ABS(evec(k,k)), ABS(evec(k,k-1)) , ABS(evec(k,k-2)) )
        mz = ABS(evec(k, k))
        b = eval(k) + SQRT(beta)*mz
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    END SUBROUTINE Estupb_real
!-------------------chebyshev_filter---------------------
    SUBROUTINE chebyshev_filter_real(nps, nst, veff, X, m, a, b)
!
        IMPLICIT NONE
!INOUT
        INTEGER(I4B), INTENT(IN) :: nps, nst !number of states
        REAL(DP), INTENT(IN) :: veff(nps) !veff
        REAL(DP), INTENT(INOUT) :: X(nps, nst)
        INTEGER(I4B), INTENT(IN)   :: m ! the m degree Chebyshev polynomial we used
        REAL(DP), INTENT(IN) :: a, b  !interval [a,b] to be filter
!LOCAL
        REAL(DP) :: e &   !(b-a)/2
              & , c    !(b+a)/2
!
        REAL(DP), DIMENSION(nps, nst)  :: HV, Y, Ynew
        INTEGER(I4B) :: Ic
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        IF (m < 2) THEN
            WRITE (*, *) 'Chebyshev_filter: the degree m must larger than 1'
        END IF

        e = (b - a)/2
        c = (b + a)/2

        CALL cal_HX_real(nps, nst, veff, X, HV)
        Y(:, :) = (HV(:, :) - c*X(:, :))/e

!Chebyshev filtering
        DO Ic = 2, m

!CALL HY
            CALL cal_HX_real(nps, nst, veff, Y, HV)
            Ynew = 2.d0*(HV - c*Y)/e - X
!store the Y

            X = Y

            Y = Ynew

        END DO
!out put filted space
        X = Ynew
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    END SUBROUTINE chebyshev_filter_real
!-------------------chebyshev_filter_scaled---------------------
    SUBROUTINE chebyshev_filter_scaled_real(nps, nst, veff, X, m, a, b, al)
!
        IMPLICIT NONE
!INOUT
        INTEGER(I4B), INTENT(IN) :: nps, nst !number of states
        REAL(DP), INTENT(IN) :: veff(nps) !veff
        REAL(DP), INTENT(INOUT) :: X(nps, nst)
        INTEGER(I4B), INTENT(IN)   :: m ! the m degree Chebyshev polynomial we used
        REAL(DP), INTENT(IN) :: a, b, al  !interval [a,b] to be filter
!LOCAL
        REAL(DP) :: e &   !(b-a)/2
              & , c &   !(b+a)/2
              & , sigma, sigmanew, tau
!
        REAL(DP), DIMENSION(nps, nst)  :: HV, Y, Ynew
        INTEGER(I4B) :: Ic
        REAL(DP) :: temp1, temp2
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        IF (m < 2) THEN
            WRITE (*, *) 'Chebyshev_filter: the degree m must larger than 1'
        END IF

        e = (b - a)/2
        c = (b + a)/2
        sigma = e/(c - al)
        tau = 2.d0/sigma
        CALL cal_HX_real(nps, nst, veff, X, HV)
        temp1 = sigma/e
        Y = (HV - c*X)*temp1
        DO Ic = 2, m
            sigmanew = 1.d0/(tau - sigma)
            CALL cal_HX_real(nps, nst, veff, Y, HV)
            temp1 = 2.d0*sigmanew/e
            temp2 = sigma*sigmanew
            Ynew = (HV - c*Y)*temp1 - temp2*X
            X = Y
            Y = Ynew
            sigma = sigmanew
        END DO
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    END SUBROUTINE chebyshev_filter_scaled_real
!-----------------------GRayleigh_Ritz--------------------------
    SUBROUTINE GRayleigh_Ritz_real(nps, nev, veff, X, D)
        USE Lapack_module, ONLY: GeneralizeEigen_real, matmat_real

        USE ScaLapack_module, ONLY: SL_GeneralizeEigen_real &
               &, SL_matmat_real_tn, SL_matmat_real_nn, twoD_map
        USE grid_module, ONLY: nsp => global_n
        USE parameters, ONLY: nev_tot => Nstates_global

        IMPLICIT NONE
!IN/OUT
        INTEGER(I4B), INTENT(IN) :: nps, nev
        REAL(DP), INTENT(IN) :: veff(nps)
        REAL(DP), INTENT(INOUT) :: X(nps, nev)
        REAL(DP), INTENT(OUT) :: D(:)
!

        INTEGER(I4B) :: m, n, mb, nb, cm, cn, cmb, cnb, nstsp
        REAL(DP), DIMENSION(twoD_map(1, parallel%rankx,&
             &parallel%ranky), nev) :: S_hat, H_hat, Q

        REAL(DP) :: Xnew(nps, nev)
!xq
!INTEGER(I4B) :: t1,t2,t3,t4,t5,t6
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        m = nsp
        n = nev_tot
        mb = BLOCK_MBNB
        nb = BLOCK_MBNB
        cm = nev_tot
        cn = nev_tot
        cmb = BLOCK_MBNB
        cnb = BLOCK_MBNB
        nstsp = parallel%nstate_proc

!Calculate the overlap matrix
!call system_clock(t1)

        CALL SL_matmat_real_tn('T', 'N', X, X, S_hat, m, n, m, n, mb, nb, mb, nb, cmb, cnb)

!call system_clock(t2)
!Calculate the project hamiltion
!call system_clock(t3)
        CALL Rayleigh_quotient_real(nps, nev, veff, X, H_hat)
!call system_clock(t4)
!solve the generalized eigenvalue problem
!call system_clock(t5)

        CALL SL_GeneralizeEigen_real(nev_tot, H_hat, S_hat, n, n, n, n &
             &, cmb, cnb, cmb, cnb, Q, n, n, D, cmb, cnb)

!call system_clock(t6)
!X=XQ

        CALL SL_matmat_real_nn('N', 'N', X, Q, Xnew, m, n, n, n, mb, nb, cmb, cnb, cmb, cnb)

!print*,'Overlap',(t2-t1)/10000.d0
!print*,'Hhat',(t4-t3)/10000.d0
!print*,'Genera',(t6-t5)/10000.d0
!STOP
        X = Xnew
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    END SUBROUTINE GRayleigh_Ritz_real
!---------------------Rayleigh-Ritz step------------------------
    SUBROUTINE Rayleigh_Ritz_real(nps, sn, veff, X, D)
        USE Lapack_module, ONLY: OrthNorm_real, diagM_real, matmat_real

        USE grid_module, ONLY: nsp => global_n
        USE ScaLapack_module, ONLY: twoD_map, SL_diagM_real,&
               &SL_matmat_real_tn, SL_matmat_real_nn, SL_OrthNorm_real
        USE parameters, ONLY: &
                        & nev_tot => Nstates_global, BLOCK_MBNB

        IMPLICIT NONE
        INTEGER(I4B), INTENT(IN) :: nps, sn
        REAL(DP), INTENT(IN) :: veff(nps)
        REAL(DP), INTENT(INOUT) :: X(nps, sn)
        REAL(DP), INTENT(OUT) :: D(:)
!LOCAL
        REAL(DP) :: Xnew(nps, sn)

        INTEGER(I4B) ::   m   & !dim. 1
                     &, n   & !dim. 2
                     &, mb  & !block of dim. 1
                     &, nb  &  !block of dim. 2
                     &, cm  &  !dim.1 of X^T*H*X
                     &, cn  &  !dim.2 of X^T*H*X
                     &, cmb & !block of dim.1 of X^T*H*X
                     &, cnb & !block of dim.2 of X^T*H*X
                     &, nstsp !number of sts of this proc
        REAL(DP), DIMENSION(twoD_map(1, parallel%rankx,&
             &parallel%ranky), sn) :: Hhat, Qs

!test
        INTEGER(I4B) :: I
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        m = nsp  !dim. 1
        n = nev_tot   !dim. 2
        mb = 1  !block of dim. 1
        nb = 1  !block of dim. 2
        cm = nev_tot  !dim.1 of X^T*H*X
        cn = nev_tot   !dim.2 of X^T*H*X
        cmb = 1  !block of dim.1 of X^T*H*X
        cnb = 1  !block of dim.2 of X^T*H*X
        nstsp = parallel%nstate_proc
!F0:OrthNorm
        CALL SL_OrthNorm_real(X, m, n, mb, nb)
!F1:Rayleigh_quotient
        CALL Rayleigh_quotient_real(nps, sn, veff, X, Hhat)
!F2:eigen-decomposition Q,D
        CALL SL_diagM_real(n, Hhat, n, n, n, n, cmb, cnb, cmb, cnb, Qs, n, n, D, cmb, cnb)
!X=MATMUL( X , Q )
        CALL SL_matmat_real_nn('N', 'N', X, Qs, Xnew, m, n, n, n, mb, nb, cmb, cmb, cmb, cnb)

        X = Xnew
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    END SUBROUTINE Rayleigh_Ritz_real
!-----------------non-OrthNorm Chebyshev_filter ----------------
    SUBROUTINE cheby_filtering_GRRr(nps, nev, veff, X, D)
        USE parameters, ONLY: LRROrthNorm, nev_tot => Nstates_global
!To avoid OrthNorm
        IMPLICIT NONE
        INTEGER(I4B), INTENT(IN) :: nps, nev
        REAL(DP), INTENT(IN) :: veff(nps)
        REAL(DP), INTENT(INOUT) :: X(nps, nev)
        REAL(DP), INTENT(INOUT) :: D(:)  !rayleigh-ritz value
!LOCAL
        REAL(DP) :: a, b, al
        REAL(DP) :: mixa = 0.3_DP

        INTEGER(I4B) :: bcastid, i, j

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!call (al,a,b)
!a=MAXVAL(D) !+ABS(ad(IGamma)-MAXVAL(D))*mixa
!print*,a,ad(IGamma)
!ad(IGamma)=a
        a = MAXVAL(D)
        al = MINVAL(D)
!CALL ChebyshevM(al,a,b,CheM)
!up boundary

!> find the core deal with the highest state
        aa: DO j = parallel%dims(1), 1, -1 !> all core
        DO i = size(parallel%sub2sum, 1), 1, -1  !> all states per core
            IF (parallel%sub2sum(i, j) == nev_tot) THEN
                bcastid = j - 1
                CALL Estupb_real(nps, 7, veff, X(:, parallel%nstate_proc), b)
                exit aa
            END IF
        END DO
        END DO aa
        CALL MPI_BCAST(b, 1, MPI_REAL8, bcastid, parallel%commy, mpinfo)

!if(parallel%isroot) print*,'a,b,al',a,b,al
!filtering (a,b)
!CALL chebyshev_filter_scaled_real(nps,nev,veff,X,CheM,a,b,al)
        CALL chebyshev_filter_real(nps, nev, veff, X, CheM, a, b)
        IF (LRROrthNorm) THEN
!RR (Rayleigh-Ritz Step)
            CALL Rayleigh_Ritz_real(nps, nev, veff, X, D)
        ELSE
!GRR (Rayleigh-Ritz Step)
            CALL GRayleigh_Ritz_real(nps, nev, veff, X, D)
        END IF
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    END SUBROUTINE cheby_filtering_GRRr
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
END MODULE chebyshev_module

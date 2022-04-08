MODULE ScaLapack_module
!###########################################################!
!*For: Lapack                                               !
!*Author : Sheng Wang                                       !
!###########################################################!
    USE constants
    USE smpi_math_module, ONLY: parallel
! BLACS INFO
    IMPLICIT NONE
    INTEGER, EXTERNAL  :: blacs_pnum
    INTEGER(i4b)  :: blacs_contxt
    INTEGER(I4B), parameter  :: DLEN = 9
    INTEGER(I4B)    :: MYROW, MYCOL
    INTEGER(I4B)    :: NPCOL, NPROW
    INTEGER(I4B)    :: my_blacs_id
    INTEGER(I4B)    :: NP, NQ
    INTEGER(I4B)    :: IAM, NPROCS
    LOGICAL, SAVE    :: INIT_CALLED = .false.
    INTEGER, EXTERNAL  :: NUMROC
    INTEGER(I4B)    :: INFO
    LOGICAL :: L_useless = .false.
    INTEGER(I4B), allocatable :: twoD_map(:, :, :)
CONTAINS
!-----------------------divided line---------------------------
    Subroutine Init_scala()  !{{{
        USE parameters, ONLY: Nstates_global

        implicit none
        LOGICAL :: Lpbc = .TRUE.
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        blacs_contxt = parallel%comm
        if (Lpbc) then
            nprow = parallel%dims(2) !1!parallel%dims(2)
            npcol = parallel%dims(1) !parallel%numprocs
        else
            nprow = parallel%numprocs!1!parallel%dims(2)
            npcol = 1!parallel%numprocs
        end if
        CALL BLACS_GET(-1, 0, blacs_contxt)
        CALL BLACS_GRIDINIT(blacs_contxt, 'Col-major', NPROW, NPCOL)
! CALL SL_INIT(blacs_contxt,NPROW,NPCOL)
        CALL BLACS_GRIDINFO(blacs_contxt, NPROW, NPCOL, MYROW, MYCOL)
! CALL BLACS_PInfo(iam,nprocs)
        INIT_CALLED = .true.
        if ((NPROW < 0) .or. (NPCOL < 0)) L_useless = .true.

        allocate (twoD_map(3, 0:nprow - 1, 0:npcol - 1))
        call twoD_map_set(Nstates_global, nprow, npcol, twoD_map)

    End Subroutine Init_scala  !}}}

! Subroutine Init_scala_sub()  !{{{

! implicit none
! blacs_contxt = parallel%subcomm
! nprow = parallel%dims(2)
! CALL BLACS_GET( -1, 0, blacs_contxt )
! CALL BLACS_GRIDINIT( blacs_contxt, 'Col-major', NPROW, NPCOL )
! CALL BLACS_GRIDINFO( blacs_contxt, NPROW, NPCOL, MYROW, MYCOL )
! CALL BLACS_PInfo(iam,nprocs)

! End Subroutine Init_scala_sub  !}}}

    Subroutine SL_orthnorm(amat, m, n, mb, nb)  !{{{
        implicit none
        complex(dcp)   :: amat(:, :)
        complex(dcp)    :: tau(10000)
        complex(dcp), allocatable     :: work(:)
        integer(i4b)    :: m, n
        integer(i4b)    :: mb, nb
        integer(i4b)    :: lwork, info
        integer(i4b)    :: desca(dlen)
        integer(i4b)    :: rsrc, csrc
!integer,external :: numroc
        integer(i4b)    :: np, nq, i, j, k

        NP = NUMROC(M, MB, MYROW, 0, NPROW)
        NQ = NUMROC(N, NB, MYCOL, 0, NPCOL)

        CALL DESCINIT(DESCA, M, N, MB, NB, 0, 0, blacs_contxt, MAX(1, NP), INFO)

        lwork = 1000000
        allocate (work(lwork))
        call pzgeqrf(m, n, amat, 1, 1, desca, tau, work, lwork, info)
        print *, info

        call pzungqr(m, n, n, amat, 1, 1, desca, tau, work, lwork, info)

        deallocate (work)

    End Subroutine SL_orthnorm  !}}}

    Subroutine SL_orthnorm_real(amat, m, n, mb, nb)  !{{{
! use system_info
        implicit none
        real(dp)   :: amat(:, :)
        real(dp), allocatable     :: tau(:)
        real(dp), allocatable     :: work(:)
        integer(i4b)    :: m, n
        integer(i4b)    :: mb, nb
        integer(i4b)    :: lwork, info
        integer(i4b)    :: desca(dlen)
        integer(i4b)    :: rsrc, csrc
        integer(i4b)    :: np, nq, i, j, k

        NP = NUMROC(M, MB, MYROW, 0, NPROW)
        NQ = NUMROC(N, NB, MYCOL, 0, NPCOL)

        CALL DESCINIT(DESCA, M, N, MB, NB, 0, 0, blacs_contxt, MAX(1, NP), INFO)

!
        allocate (tau(min(M, N)))
! lwork = 1000000
! allocate(work(lwork))
        allocate (work(3))
        lwork = -1
        call pdgeqrf(M, N, amat, 1, 1, desca, tau, work, lwork, info)
        lwork = work(1)
        deallocate (work)
        allocate (work(lwork))

        call pdgeqrf(M, N, amat, 1, 1, desca, tau, work, lwork, info)
        if (info < 0) then
            print *, 'error in scalapack pdgeqrf, info', info
            stop
        end if

        call pdorgqr(m, n, min(M, N), amat, 1, 1, desca, tau, work, lwork, info)
        if (info < 0) then
            print *, 'error in scalapack pdorgqr, info', info
            stop
        end if

        deallocate (work)
        deallocate (tau)

    End Subroutine SL_orthnorm_real  !}}}

    Subroutine SL_matmat(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb, cmbin, cnbin)  !{{{
!use system_info
        implicit none
        complex(dcp), intent(in)  :: amat(:, :), bmat(:, :)
        complex(dcp), intent(out)  :: cmat(:, :)
        character(1), intent(in) :: opa, opb
        integer(i4b)    :: m, n, k
        integer(i4b)    :: am, an, amb, anb, bm, bn, bmb, bnb, cm, cn, cmb, cnb
        integer(i4b), optional :: cmbin, cnbin
        integer(i4b)    :: desca(dlen), descb(dlen), descc(dlen)
        complex(dcp)     :: alpha = 1.d0, &
                          & beta = 0.d0
        real(dp) :: t
!---------------------------------------------------------------------

        if (NPROW <= 0 .or. NPCOL <= 0) then
            print *, 'err,some core useless in scala'
            return
        end if

        if (opa == 'N' .or. opa == 'n') then
            m = am
            cm = am
!cn = bn
        else
            m = an
            cm = an
!cn = bn
        end if

        if (opb == 'N' .or. opb == 'n') then
            cn = bn
        else
            cn = bm
        end if

        if (opb == 'N' .or. opb == 'n') then
            n = bn
            k = bm
        else
            n = bm
            k = bn
        end if

        if ((.not. present(cmbin)) .and. (.not. present(cnbin))) then
            cmb = m
            cnb = n
            cmat = 0.d0
        elseif (present(cmbin) .and. present(cnbin)) then
            cmb = cmbin
            cnb = cnbin
        end if
        NP = NUMROC(AM, AMB, MYROW, 0, NPROW)
        NQ = NUMROC(AN, ANB, MYCOL, 0, NPCOL)
        CALL DESCINIT(DESCA, AM, AN, AMB, ANB, 0, 0, blacs_contxt, MAX(1, NP), INFO)
        NP = NUMROC(BM, BMB, MYROW, 0, NPROW)
        NQ = NUMROC(BN, BNB, MYCOL, 0, NPCOL)
        CALL DESCINIT(DESCB, BM, BN, BMB, BNB, 0, 0, blacs_contxt, MAX(1, NP), INFO)

        NP = NUMROC(CM, CMB, MYROW, 0, NPROW)
        NQ = NUMROC(CN, CNB, MYCOL, 0, NPCOL)
        CALL DESCINIT(DESCC, CM, CN, CMB, CNB, 0, 0, blacs_contxt, MAX(1, NP), INFO)

        call pzgemm(opa, opb, m, n, k, alpha, amat, 1, 1, desca,  &
                    & bmat, 1, 1, descb, beta, cmat, 1, 1, descc, info)

    End Subroutine SL_matmat  !}}}

    Subroutine SL_matmat_cmplx_cn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb, cmb, cnb)  !{{{
        USE smpi_math_module, ONLY: MPI_complex16, mpinfo, MPI_SUM
        implicit none
        complex(dcp), intent(in)  :: amat(:, :), bmat(:, :)
        complex(dcp), intent(out)  :: cmat(:, :)
        character(1), intent(in) :: opa, opb
        integer(i4b)    :: m, n, k
        integer(i4b)    :: am, an, amb, anb, bm, bn, bmb, bnb, cm, cn, cmb, cnb
        integer(i4b)    :: desca(dlen), descb(dlen), descc(dlen)
        complex(dcp)     :: alpha = 1.d0, &
             & beta = 0.d0
        complex(dcp), allocatable :: amat_y(:, :), bmat_y(:, :)&
             & , cmat_local(:, :), cmat_global(:, :)
        integer(i4b) :: y_recvcounts(parallel%dims(1))&
             & , y_displs(parallel%dims(1))
        integer(i4b) :: lda, ldb, ldc
        integer(i4b) :: i
        integer(i4B) :: n_local
        integer(i4b) :: x_i, x_j, y_i, y_j, n_block, pad_block, pad_data
        integer(i4b) :: j, i_x(an), i_y(bn), il, im
!---------------------------------------------------------------------
        m = an
        n = an
        n_local = size(bmat, 2)
        cm = m
        cn = n
        k = size(amat, 1)
        lda = max(size(amat, 1), 1)
        ldb = max(size(bmat, 1), 1)
        ldc = max(m, 1)

!> shift x
        n_block = cm/cmb/parallel%dims(2)
        pad_block = mod(cm/cmb, parallel%dims(2))
        pad_data = mod(cm, cmb)
        if (parallel%rankx < pad_block) then
            x_i = parallel%rankx*(n_block + 1)*cmb + 1
        elseif (parallel%rankx == pad_block) then
            x_i = parallel%rankx*n_block*cmb + pad_block + pad_data + 1
        else
            x_i = parallel%rankx*n_block*cmb + pad_block + 1
        end if
        x_j = x_i + twoD_map(1, parallel%rankx, parallel%ranky) - 1

!> shift y
        n_block = cn/cnb/parallel%dims(1)
        pad_block = mod(cn/cnb, parallel%dims(1))
        pad_data = mod(cn, cnb)
        if (parallel%ranky < pad_block) then
            y_i = parallel%ranky*(n_block + 1)*cnb + 1
        elseif (parallel%ranky == pad_block) then
            y_i = parallel%ranky*n_block*cnb + pad_block + pad_data + 1
        else
            y_i = parallel%ranky*n_block*cnb + pad_block + 1
        end if
        y_j = y_i + twoD_map(2, parallel%rankx, parallel%ranky) - 1

!> gater in row
        allocate (amat_y(size(amat, 1), an))
        allocate (bmat_y(size(amat, 1), an))
        y_recvcounts = twoD_map(2, parallel%rankx, :)*size(amat, 1)
        y_displs(1) = 0
        do i = 2, parallel%dims(1)
            y_displs(i) = sum(y_recvcounts(:i - 1))
        end do
        call MPI_ALLGATHERV(amat, size(amat)&
             &, MPI_COMPLEX16, amat_y, y_recvcounts&
             &, y_displs, MPI_COMPLEX16&
             & , parallel%commy, mpinfo)

!> local matmat
        allocate (cmat_local(an, an))
        allocate (cmat_global(an, an))
        cmat_local = 0.d0
        CALL ZGEMM(opA, opB, M, N_local, K, alpha, amat_y, LDA, bmat, LDB, beta, cmat_local(:, y_i:y_j), LDC)

!> allreduce the local result
        CALL MPI_ALLREDUCE(Cmat_local, cmat_global, size(cmat_local), MPI_COMPLEX16,&
             & MPI_SUM, parallel%comm, mpinfo)
! & MPI_SUM, parallel%commx, mpinfo)

!> remain the local part
! Cmat=cmat_global(x_i:x_j,y_i:y_j)
        il = 1; im = 1
        do j = 0, parallel%dims(2) - 1
            do i = 1, an, 1
                if (mod(i - 1, parallel%dims(2)) == j) then
                    i_x(il) = i
                    il = il + 1
                end if
            end do
        end do
        do j = 0, parallel%dims(1) - 1
            do i = 1, an, 1
                if (mod(i - 1, parallel%dims(1)) == j) then
                    i_y(i) = im
                    im = im + 1
                end if
            end do
        end do
        do i = x_i, x_j, 1
            cmat(i + 1 - x_i, :) = cmat_global(i_y(i_x(i)), y_i:y_j)
        end do

!> test
! do j=1,parallel%numprocs,1
!    write(str_id,'(I4)')j
!    if(parallel%myid==j-1)then
!       open(1111+j,file="global_data"//trim(adjustl(str_id)))
!       ! write(1111+j,*)hhat_local
!       do i=1,size(cmat_global,2),1
!          write(1111+j,'(1000(F6.2,2X))')cmat_global(:,i)
!       enddo
!       close(1111+j)
!    endif
! enddo
!>
        deallocate (cmat_local, cmat_global)

    End Subroutine SL_matmat_cmplx_cn  !}}}

    Subroutine SL_matmat_cmplx_nn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb, cmb, cnb)  !{{{
        USE smpi_math_module, ONLY: MPI_complex16, mpinfo, MPI_SUM
        implicit none
        complex(dcp), intent(in)  :: amat(:, :), bmat(:, :)
        complex(dcp), intent(out)  :: cmat(:, :)
        character(1), intent(in) :: opa, opb
        integer(i4b)    :: m, n, k
        integer(i4b)    :: am, an, amb, anb, bm, bn, bmb, bnb, cm, cn, cmb, cnb
        integer(i4b)    :: desca(dlen), descb(dlen), descc(dlen)
        complex(dcp)     :: alpha = 1.d0, &
             & beta = 0.d0
        complex(dcp), allocatable :: Bmat_local(:, :), Bmat_global(:, :), amat_y(:, :)
        integer(i4b) :: lda, ldb, ldc
        integer(i4b) :: i, j, il, y_recvcounts(parallel%dims(1))&
           & , y_displs(parallel%dims(1))

!---------------------------------------------------------------------
!> cmat_local
        allocate (Bmat_local(bm, size(Bmat, 2)))
        allocate (Bmat_global(bm, size(Bmat, 2)))
!> m,n,k
        m = size(amat, 1)
        k = size(bmat_local, 1)
        n = size(bmat_local, 2)
!> lda,ldb,ldc
        lda = max(m, 1)
        ldb = max(k, 1)
        ldc = max(m, 1)

!> assign to local bmat
        il = 0
        Bmat_local = 0.d0
        do i = 1, bm, 1
            if (mod(i - 1, parallel%dims(2)) /= parallel%rankx) cycle
            il = il + 1
            Bmat_local(i, :) = Bmat(il, :)
        end do
        CALL MPI_ALLREDUCE(Bmat_local, Bmat_global, size(Bmat_local), MPI_COMPLEX16,&
             & MPI_SUM, parallel%commx, mpinfo)
!>> reorder the bmat_global
        il = 0
        do j = 0, parallel%dims(1) - 1
            do i = 1, bm, 1
                if (mod(i - 1, parallel%dims(1)) /= j) cycle
                il = il + 1
                Bmat_local(il, :) = Bmat_global(i, :)
            end do
        end do

!> amat_global
        allocate (amat_y(m, k))
        y_recvcounts = twoD_map(2, parallel%rankx, :)*size(amat, 1)
        y_displs(1) = 0
        do i = 2, parallel%dims(1)
            y_displs(i) = sum(y_recvcounts(:i - 1))
        end do
        call MPI_ALLGATHERV(amat, size(amat)&
             &, MPI_COMPLEX16, amat_y, y_recvcounts&
             &, y_displs, MPI_COMPLEX16&
             & , parallel%commy, mpinfo)

!> lapack
        CALL ZGEMM(opA, opB, M, N, K, alpha, Amat_y, LDA, Bmat_local, LDB, beta, Cmat, LDC)

!>
        deallocate (Bmat_local, Bmat_global, amat_y)

    End Subroutine SL_matmat_cmplx_nn  !}}}

    Subroutine SL_matmat_sub(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb, cmbin, cnbin)  !{{{
!use system_info
        implicit none
        complex(dcp), intent(in)  :: amat(:, :), bmat(:, :)
        complex(dcp), intent(out)  :: cmat(:, :)
        character(1), intent(in) :: opa, opb
        integer(i4b)    :: m, n, k
        integer(i4b)    :: am, an, amb, anb, bm, bn, bmb, bnb, cm, cn, cmb, cnb
        integer(i4b), optional :: cmbin, cnbin
        integer(i4b)    :: desca(dlen), descb(dlen), descc(dlen)
        complex(dcp)     :: alpha = 1.d0, &
                          & beta = 0.d0
        real(dp) :: t
!---------------------------------------------------------------------
!blacs_contxt = parallel%comm
!nprow = parallel%numprocs
!npcol = 1
!CALL BLACS_GET( -1, 0, blacs_contxt )
!CALL BLACS_GRIDINIT( blacs_contxt, 'Row-major', NPROW, NPCOL )
!CALL BLACS_GRIDINFO( blacs_contxt, NPROW, NPCOL, MYROW, MYCOL )
!call blacs_pinfo(iam,nprocs)

!am = grid%tn
!bm = grid%tn
!if ( .not. init_called) then
! call init_scala_sub()
!end if

        if (opa == 'N' .or. opa == 'n') then
            m = am
            cm = am
!cn = bn
        else
            m = an
            cm = an
!cn = bn
        end if

        if (opb == 'N' .or. opb == 'n') then
            cn = bn
        else
            cn = bm
        end if

        if (opb == 'N' .or. opb == 'n') then
            n = bn
            k = bm
        else
            n = bm
            k = bn
        end if

!if (parallel%isroot .and. present(cmbin) ) then
!print * , 'debug sl_matmat I' , cmbin,cnbin
!end if

        if ((.not. present(cmbin)) .and. (.not. present(cnbin))) then
            cmb = m
            cnb = n
            cmat = 0.d0
        elseif (present(cmbin) .and. present(cnbin)) then
            cmb = cmbin
            cnb = cnbin
        end if

!if (parallel%isroot) then
!print * , 'debug sl_matmat II' , cmb,cnb
!end if

        NP = NUMROC(AM, AMB, MYROW, 0, NPROW)
        NQ = NUMROC(AN, ANB, MYCOL, 0, NPCOL)
        CALL DESCINIT(DESCA, AM, AN, AMB, ANB, 0, 0, blacs_contxt, MAX(1, NP), INFO)

        NP = NUMROC(BM, BMB, MYROW, 0, NPROW)
        NQ = NUMROC(BN, BNB, MYCOL, 0, NPCOL)
        CALL DESCINIT(DESCB, BM, BN, BMB, BNB, 0, 0, blacs_contxt, MAX(1, NP), INFO)

!if (parallel%isroot .and. present(cmbin)) then
!print * , 'debug sl_matmat III' , cmbin,cnbin
!end if

        NP = NUMROC(CM, CMB, MYROW, 0, NPROW)
        NQ = NUMROC(CN, CNB, MYCOL, 0, NPCOL)
        CALL DESCINIT(DESCC, CM, CN, CMB, CNB, 0, 0, blacs_contxt, MAX(1, NP), INFO)
! CALL DESCINIT( DESCC, CM, CN, CMB, CNB, 0, 0, blacs_contxt,11,INFO )
!if (parallel%isroot) then
!print * , 'debug sl_matmat VI', am,an,bm,bn,amb,anb,bmb,bnb,cm,cn,cmb,cnb,m,n,k
!!print * , 'debug descc' , descc
!end if
!call mpi_barrier(parallel%comm,mpinfo)

!call start_time('scala')
        call pzgemm(opa, opb, m, n, k, alpha, amat, 1, 1, desca,  &
                    & bmat, 1, 1, descb, beta, cmat, 1, 1, descc, info)
!call end_time('scala')
!call print_time('scala',t)

!if (parallel%isroot) then
!print * , 'debug ok matmat' , t
!endif

    End Subroutine SL_matmat_sub  !}}}

    Subroutine SL_matmat_sub_real(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb, cmbin, cnbin)  !{{{
!use system_info
        implicit none
        real(dp), intent(in)  :: amat(:, :), bmat(:, :)
        real(dp), intent(out)  :: cmat(:, :)
        character(1), intent(in) :: opa, opb
        integer(i4b)    :: m, n, k
        integer(i4b)    :: am, an, amb, anb, bm, bn, bmb, bnb, cm, cn, cmb, cnb
        integer(i4b), optional :: cmbin, cnbin
        integer(i4b)    :: desca(dlen), descb(dlen), descc(dlen)
        real(dp)     :: alpha = 1.d0, &
                          & beta = 0.d0
        real(dp) :: t
!---------------------------------------------------------------------
!blacs_contxt = parallel%comm
!nprow = parallel%numprocs
!npcol = 1
!CALL BLACS_GET( -1, 0, blacs_contxt )
!CALL BLACS_GRIDINIT( blacs_contxt, 'Row-major', NPROW, NPCOL )
!CALL BLACS_GRIDINFO( blacs_contxt, NPROW, NPCOL, MYROW, MYCOL )
!call blacs_pinfo(iam,nprocs)

!am = grid%tn
!bm = grid%tn
!if ( .not. init_called) then
! call init_scala_sub()
!end if

        if (opa == 'N' .or. opa == 'n') then
            m = am
            cm = am
!cn = bn
        else
            m = an
            cm = an
!cn = bn
        end if

        if (opb == 'N' .or. opb == 'n') then
            cn = bn
        else
            cn = bm
        end if

        if (opb == 'N' .or. opb == 'n') then
            n = bn
            k = bm
        else
            n = bm
            k = bn
        end if

!if (parallel%isroot .and. present(cmbin) ) then
!print * , 'debug sl_matmat I' , cmbin,cnbin
!end if

        if ((.not. present(cmbin)) .and. (.not. present(cnbin))) then
            cmb = m
            cnb = n
            cmat = 0.d0
        elseif (present(cmbin) .and. present(cnbin)) then
            cmb = cmbin
            cnb = cnbin
        end if

!if (parallel%isroot) then
!print * , 'debug sl_matmat II' , cmb,cnb
!end if

        NP = NUMROC(AM, AMB, MYROW, 0, NPROW)
        NQ = NUMROC(AN, ANB, MYCOL, 0, NPCOL)
        CALL DESCINIT(DESCA, AM, AN, AMB, ANB, 0, 0, blacs_contxt, MAX(1, NP), INFO)

        NP = NUMROC(BM, BMB, MYROW, 0, NPROW)
        NQ = NUMROC(BN, BNB, MYCOL, 0, NPCOL)
        CALL DESCINIT(DESCB, BM, BN, BMB, BNB, 0, 0, blacs_contxt, MAX(1, NP), INFO)

!if (parallel%isroot .and. present(cmbin)) then
!print * , 'debug sl_matmat III' , cmbin,cnbin
!end if

        NP = NUMROC(CM, CMB, MYROW, 0, NPROW)
        NQ = NUMROC(CN, CNB, MYCOL, 0, NPCOL)
        CALL DESCINIT(DESCC, CM, CN, CMB, CNB, 0, 0, blacs_contxt, MAX(1, NP), INFO)
! CALL DESCINIT( DESCC, CM, CN, CMB, CNB, 0, 0, blacs_contxt,11,INFO )
!if (parallel%isroot) then
!print * , 'debug sl_matmat VI', am,an,bm,bn,amb,anb,bmb,bnb,cm,cn,cmb,cnb,m,n,k
!!print * , 'debug descc' , descc
!end if
!call mpi_barrier(parallel%comm,mpinfo)

!call start_time('scala')
        call pdgemm(opa, opb, m, n, k, alpha, amat, 1, 1, desca,  &
                    & bmat, 1, 1, descb, beta, cmat, 1, 1, descc, info)
!call end_time('scala')
!call print_time('scala',t)

!if (parallel%isroot) then
!print * , 'debug ok matmat' , t
!endif

    End Subroutine SL_matmat_sub_real  !}}}

    Subroutine SL_matmat_real(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb, cmbin, cnbin)  !{{{
!use system_info
        implicit none
!real(dp),intent(in)  :: amat(:,:),bmat(:,:)
        real(dp)  :: amat(:, :), bmat(:, :)
        real(dp), intent(out)  :: cmat(:, :)
        character(1), intent(in) :: opa, opb
        integer(i4b)    :: m, n, k
        integer(i4b)    :: am, an, amb, anb, bm, bn, bmb, bnb, cm, cn, cmb, cnb
        integer(i4b), optional :: cmbin, cnbin
        integer(i4b)    :: desca(dlen), descb(dlen), descc(dlen)
        complex(dcp)     :: alpha = 1.d0, &
                          & beta = 0.d0
        real(dp) :: t
        integer(i4b), save  :: itest = 0
!---------------------------------------------------------------------
        itest = itest + 1
!blacs_contxt = parallel%comm
!nprow = parallel%numprocs
!npcol = 1
!CALL BLACS_GET( -1, 0, blacs_contxt )
!CALL BLACS_GRIDINIT( blacs_contxt, 'Row-major', NPROW, NPCOL )
!CALL BLACS_GRIDINFO( blacs_contxt, NPROW, NPCOL, MYROW, MYCOL )
!call blacs_pinfo(iam,nprocs)

!am = grid%tn
!bm = grid%tn
!    cmat(:,:) = 0.d0

        if (opa == 'N' .or. opa == 'n') then
            m = am
            cm = am
!cn = bn
        else
            m = an
            cm = an
!cn = bn
        end if

        if (opb == 'N' .or. opb == 'n') then
            cn = bn
        else
            cn = bm
        end if

        if (opb == 'N' .or. opb == 'n') then
            n = bn
            k = bm
        else
            n = bm
            k = bn
        end if
!if (parallel%isroot .and. present(cmbin) ) then
!print * , 'debug sl_matmat I' , cmbin,cnbin
!end if

        if ((.not. present(cmbin)) .and. (.not. present(cnbin))) then
            cmb = m
            cnb = n
            cmat = 0.d0
        elseif (present(cmbin) .and. present(cnbin)) then
            cmb = cmbin
            cnb = cnbin
        end if

!if (parallel%coords(2) == 0) then
!print * , 'debug sl_matmat II' , cmb,cnb,parallel%myid
!end if

        NP = NUMROC(AM, AMB, MYROW, 0, NPROW)
        NQ = NUMROC(AN, ANB, MYCOL, 0, NPCOL)
!print *,"amb,anb...",amb,anb,bmb,bnb,cmb,cnb
        CALL DESCINIT(DESCA, AM, AN, AMB, ANB, 0, 0, blacs_contxt, MAX(1, NP), INFO)

        NP = NUMROC(BM, BMB, MYROW, 0, NPROW)
        NQ = NUMROC(BN, BNB, MYCOL, 0, NPCOL)
        CALL DESCINIT(DESCB, BM, BN, BMB, BNB, 0, 0, blacs_contxt, MAX(1, NP), INFO)

!if (parallel%isroot .and. present(cmbin)) then
!print * , 'debug sl_matmat III' , cmbin,cnbin
!end if

        NP = NUMROC(CM, CMB, MYROW, 0, NPROW)
        NQ = NUMROC(CN, CNB, MYCOL, 0, NPCOL)
        CALL DESCINIT(DESCC, CM, CN, CMB, CNB, 0, 0, blacs_contxt, MAX(1, NP), INFO)
! CALL DESCINIT( DESCC, CM, CN, CMB, CNB, 0, 0, blacs_contxt,11,INFO )
!if (parallel%coords(2) == 0) then
!print * , 'debug sl_matmat VI', am,an,bm,bn,amb,anb,bmb,bnb,cm,cn,cmb,cnb,m,n,k
!print * , ' np = ', np
!print * , 'debug descc' , descc
!end if

!call start_time('scala')
!if (parallel%isroot) then
!print * , 'debug advmpi desca' , desca
!print * , 'debug advmpi descb' , descb
!print * , 'debug advmpi descc' , descc
!end if
        call pdgemm(opa, opb, m, n, k, alpha, amat, 1, 1, desca,  &
                    & bmat, 1, 1, descb, beta, cmat, 1, 1, descc, info)

!print *,'info = ' , info , parallel%myid
!if (parallel%isroot) then
!print * , 'MATMAT Info' , info
!end if
!call mpi_barrier(parallel%comm,mpinfo)
!call end_time('scala')
!call print_time('scala',t)

!if (parallel%isroot) then
!print * , 'debug ok matmat' , t
!endif

    End Subroutine SL_matmat_real  !}}}

    SUBROUTINE SL_GeneralizeEigen(dime, amat, bmat, am, an, bm, bn, amb, anb, bmb, bnb, evec, cm, cn, eval, cmbin, cnbin)  !{{{
        IMPLICIT NONE
        integer(i4b), intent(in) :: dime
        COMPLEX(DCP), INTENT(IN)  :: amat(:, :), bmat(:, :)
        COMPLEX(DCP), INTENT(OUT) :: evec(:, :)
        REAL(DP), INTENT(OUT) :: eval(:)
        COMPLEX(DCP), allocatable :: U(:, :)
!LOCAL
        integer(i4b) :: am, an, bm, bn, amb, anb, bmb, bnb, cm, cn
        integer(i4b), optional  :: cmbin, cnbin
        integer(i4b) :: cmb, cnb
        integer(i4b) :: m, nevec
!integer(i4b),parameter :: lwork = 100000,lrwork=100000,liwork=100000
        integer(i4b) :: lwork, liwork, lrwork
        complex(dcp), allocatable :: work(:)
        real(dp), allocatable    :: rwork(:)
        integer(i4b), allocatable :: iwork(:)
        integer(i4b) :: lda
        integer(i4b)    :: desca(dlen), descb(dlen), descz(dlen)
!real(dp)     :: mone = 1E-3

! local arrays
        integer(I4B) :: iclustr(NPCOL*NPROW*2)
        integer(I4B) :: ifail(dime)
        real(dp)     :: w(am)
        real(dp)     :: gap(NPCOL*NPROW)
        real(dp)     :: abstol
        REAL(DP), external :: PDLAMCH
! integer(i4b) :: neig ,np0,nq0,nclust,nnp
! integer(i4b) :: nsytrd_lwopt,nsygst_lwopt,a_nb,sqnpc,nps
! real(dp),external :: pdlamch
! integer(i4b),external :: iceil,pjlaenv
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        if (NPROW <= 0 .and. NPCOL <= 0) then
            print *, 'error,rank=', parallel%myid, 'is useless'
            return
        end if
        if (present(cmbin) .and. present(cnbin)) then
            cmb = cmbin
            cnb = cnbin
        end if
        abstol = 2*PDLAMCH(blacs_contxt, 'S')
!> Init the descriptor of the array
        NP = NUMROC(AM, AMB, MYROW, 0, NPROW)
        NQ = NUMROC(AN, ANB, MYCOL, 0, NPCOL)
        CALL DESCINIT(DESCA, AM, AN, AMB, AMB, 0, 0, blacs_contxt, MAX(1, NP), INFO)
        NP = NUMROC(BM, BMB, MYROW, 0, NPROW)
        NQ = NUMROC(BN, BNB, MYCOL, 0, NPCOL)
        CALL DESCINIT(DESCB, BM, BN, BMB, BMB, 0, 0, blacs_contxt, MAX(1, NP), INFO)
        NP = NUMROC(CM, CMB, MYROW, 0, NPROW)
        NQ = NUMROC(CN, CNB, MYCOL, 0, NPCOL)
        CALL DESCINIT(DESCZ, CM, CN, CMB, CMB, 0, 0, blacs_contxt, MAX(1, NP), INFO)
!> query the size of [ir]?work
        lwork = -1
        liwork = -1
        lrwork = -1
! if (allocated(work)) deallocate(work)
! if (allocated(rwork)) deallocate(rwork)
! if (allocated(iwork)) deallocate(iwork)
        allocate (work(3))
        allocate (iwork(3))
        allocate (rwork(3))! allocate(rwork(5000000))
        call pzhegvx(1, 'V', 'A', 'U', an, amat, 1, 1, desca, bmat, 1, 1,&
             &descb, 0.d0, 0.d0, 0, 0, abstol, m, nevec, eval, 0.d0, evec, 1, 1&
             & , descz, work, lwork, rwork, lrwork, iwork, liwork&
             & , ifail, iclustr, gap, info)
        lwork = nint(real(work(1), DP), I4B)
        liwork = iwork(1)
        lrwork = nint(rwork(1), I4B)
        if (allocated(work)) deallocate (work)
        if (allocated(iwork)) deallocate (iwork)
        if (allocated(rwork)) deallocate (rwork)
        allocate (work(lwork))
        allocate (iwork(liwork))
        allocate (rwork(lrwork))

!> calculate the eigenvalues and eigenvectors
        call pzhegvx(1, 'V', 'A', 'U', an, amat, 1, 1, desca, bmat, 1, 1&
             &, descb, 0.d0, 0.d0, 0, 0, abstol, m, nevec, eval, 0.d0,&
             &evec, 1, 1, descz, work, lwork, rwork, lrwork, iwork, liwork&
             &, ifail, iclustr, gap, info)

!> error info
        IF (info /= 0) THEN
            print *, 'PZHEGV err:info', info, ifail(1), parallel%myid
            STOP
        END IF
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    END SUBROUTINE SL_GeneralizeEigen  !}}}

    SUBROUTINE SL_GeneralizeEigen_real(dime, amat, bmat, am, an, bm, bn, amb, anb, bmb, bnb, evec, cm, cn, eval, cmbin, cnbin)  !{{{
        USE m_time_evaluate, ONLY: memory_sum, memory_free
        IMPLICIT NONE
        integer(i4b), intent(in) :: dime
        real(DP), INTENT(IN)  :: amat(:, :), bmat(:, :)
        real(DP), INTENT(OUT) :: evec(:, :)
        REAL(DP), INTENT(OUT) :: eval(:)
        real(DP), allocatable :: U(:, :)
!LOCAL
        integer(i4b) :: am, an, bm, bn, amb, anb, bmb, bnb, cm, cn
        integer(i4b), optional  :: cmbin, cnbin
        integer(i4b) :: cmb, cnb
        integer(i4b) :: m, nevec
        integer(i4b) :: lwork, liwork
        real(dp), allocatable :: work(:)
        integer(i4b), allocatable :: iwork(:)
!integer(i4b) :: iwkt(10)
!real(dp) :: wkt(10)
        real(dp) :: abtol
!integer(i4b) :: lrowrk
        integer(i4b) :: lda
        integer(i4b) :: maxprocs
        integer(i4b)    :: desca(dlen), descb(dlen), descz(dlen)
!real(dp)     :: mone = 1E-3

! local arrays
!integer(i4b),allocatable :: iclustr(parallel%numprocs*2),ifail(dime)
        integer(i4b) :: iclustr(parallel%numprocs*2), ifail(dime)
        real(dp)     :: w(am)
        real(dp)     :: gap(parallel%numprocs)
        real(dp), external :: pdlamch
        integer(i4b), external :: iceil, pjlaenv
        integer(i4b) :: neig, np0, nq0, nclust, nnp
        integer(i4b) :: nsytrd_lwopt, nsygst_lwopt, a_nb, sqnpc, nps
        INTEGER(I4B) :: i, ia, ja !> For Debug
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!if (parallel%isroot) then
!print *,'debug sl_generalize',am,an,bm,bn,amb,anb,bmb,bnb
!endif
        if (present(cmbin) .and. present(cnbin)) then
            cmb = cmbin
            cnb = cnbin
        end if
!evec(:,:)=amat(:,:)
!if (allocated(u)) deallocate(u)
!allocate(u(size(bmat,1),size(bmat,2)))
        abtol = 1e-6   !2*pdlamch(blacs_contxt,'s')
!if (parallel%isroot) print * ,'debug abtol ' , abtol
!U(:,:)=bmat(:,:)
!if (parallel%coords(1) ==0) print * ,'debug sca advmpi I'
!call mpi_barrier(parallel%comm,mpinfo)
!
! Init the descriptor of the array
! print *,"AMB,ANB",AMB,ANB
        NP = NUMROC(AM, AMB, MYROW, 0, NPROW)
        NQ = NUMROC(AN, ANB, MYCOL, 0, NPCOL)
! print *,"MYROW,NPROW,NP",MYROW,NPROW,NP,"process_id,process_nstate",parallel%myid,parallel%nstate_proc
! print *,"MYCOL,NPCOL,NQ",MYCOL,NPCOL,NQ,"process_id,process_nstate",parallel%myid,parallel%nstate_proc
        CALL DESCINIT(DESCA, AM, AN, AMB, ANB, 0, 0, blacs_contxt, MAX(1, NP), INFO)
!CALL DESCINIT( DESCA, AM, AN, ANB, ANB, 0, 0, blacs_contxt, MAX( 1, NP ),INFO )
!call mpi_barrier(parallel%comm,mpinfo)

        NP = NUMROC(BM, BMB, MYROW, 0, NPROW)
        NQ = NUMROC(BN, BNB, MYCOL, 0, NPCOL)
!CALL DESCINIT( DESCB, BM, BN, BMB, BMB, 0, 0, blacs_contxt, MAX( 1, NP ),INFO )
        CALL DESCINIT(DESCB, BM, BN, BMB, BNB, 0, 0, blacs_contxt, MAX(1, NP), INFO)
!if (parallel%coords(2) ==0) print * ,'descb iner gene',descb
!CALL DESCINIT( DESCB, BM, BN, BNB, BNB, 0, 0, blacs_contxt, MAX( 1, NP ),INFO )

        NP = NUMROC(CM, CMB, MYROW, 0, NPROW)
        NQ = NUMROC(CN, CNB, MYCOL, 0, NPCOL)
        CALL DESCINIT(DESCZ, CM, CN, CMB, CNB, 0, 0, blacs_contxt, MAX(1, NP), INFO)
!CALL DESCINIT( DESCZ, CM, CN, CNB, CNB, 0, 0, blacs_contxt, MAX( 1, NP ),INFO )
!
!stop
! call the pzhegvx
!    lwork = -1
!    liwork = -1
!if (parallel%isroot) print * ,'debug lwork' , lwork ,'liwork ',liwork
!    wkt = 0.d0
!    iwkt = 0
!call pdsygvx(1,'V','A','U',an,evec,1,1,desca,u,1,1,descb,0.d0,0.d0, &
!            & 0,0,1E-5,m,nevec,eval,-1E-9,evec,1,1,descz,wkt,lwork,    &
!            & iwkt,liwork,ifail,iclustr,gap,info )
!lwork = int(wkt(1))
!liwork = int(iwkt(1))
!if (parallel%isroot) print * , 'calc lwork 1',lwork
        nnp = max(am, parallel%numprocs + 1, 4)
        neig = an
        np0 = numroc(am, amb, 0, 0, nprow)
        nq0 = max(numroc(neig, anb, 0, 0, npcol), anb)
        lwork = 5*am + max(5*am, np0*nq0 + 2*anb*anb) + iceil(neig, parallel%numprocs)*am + 2*anb*anb + 4*anb*am
        liwork = 6*am
!if (parallel%isroot) print * , 'ter 1',5*am,'ter2',max(5*am,np0*nq0+2*anb*anb),'ter 3',iceil(neig,parallel%numprocs)*am &
!    & ,'ter 3',2*anb*anb,'ter 4',3*anb*am

!if (parallel%isroot) print * , 'calc lwork 2',lwork

        a_nb = pjlaenv(blacs_contxt, 3, 'PDSYTTRD', 'L', 0, 0, 0, 0)
        sqnpc = int(sqrt(dble(parallel%numprocs)))
        nps = max(numroc(am, 1, 0, 0, sqnpc), 2*a_nb)
        nsytrd_lwopt = am + 2*(a_nb + 1)*(4*nps + 2) + (nps + 3)*nps
        nsygst_lwopt = 2*np0*anb + nq0*anb + anb*anb
        lwork = max(lwork, 5*am + nsytrd_lwopt, nsygst_lwopt)

!if (parallel%isroot) print * , 'calc lwork 3',lwork
!lwork = 5000000
!lwork = 10
!liwork = 3000000
!if (parallel%myid == 0) print * ,' debug iner gene 3'
        if (allocated(work)) then
            call memory_free('SL_GeneralizeEigen_real_work', real(size(work), DP)*DP)
            deallocate (work)
        end if
        if (allocated(iwork)) then
            call memory_free('SL_GeneralizeEigen_real_iwork', real(size(iwork), DP)*DP)
            deallocate (iwork)
        end if
        allocate (work(lwork))
        allocate (iwork(liwork))
        call memory_sum('SL_GeneralizeEigen_real', real(size(work), DP)*DP + size(iwork)*DP)
!print *,'debug desca',np
        call blacs_gridinfo(blacs_contxt, nprow, npcol, myrow, mycol)
!ia=parallel%sub2sum(1,parallel%myid+1)
!ja=1
!print*,"ia, process id: ",ia,parallel%myid
!print *,"N,shape(A):",an,shape(bmat)
! CALL MPI_Barrier(parallel%comm,mpinfo)
!CALL PDPOTRF("U", an, bmat, 1, 1, DESCB, INFO )
!print *,"the unimplement Cholesky factorization, info",info

! print *, 'debug size desca' , shape(desca)
! print *, 'debug gridinfo' , (desca(i),i=1,9,1)
!if (parallel%isroot) print * ,'debug lwork' , lwork ,'liwork',liwork,'iceil',iceil(neig,parallel%numprocs)
        call pdsygvx(1, 'V', 'A', 'U', an, amat, 1, 1, desca, bmat, 1, 1, descb, 0.d0, 0.d0, &
                    & 1, an, abtol, m, nevec, eval, 1e-12, evec, 1, 1, descz, work, lwork,    &
                    & iwork, liwork, ifail, iclustr, gap, info)

        if (allocated(work)) then
            call memory_free('SL_GeneralizeEigen_real', real(size(work), DP)*DP)
            deallocate (work)
        end if
        if (allocated(iwork)) then
            call memory_free('SL_GeneralizeEigen_real', real(size(iwork), DP)*DP)
            deallocate (iwork)
        end if
!call mpi_barrier(parallel%comm,mpinfo)
!stop

!if (parallel%isroot) print *,'the number of eigenvector',m
        IF (info /= 0) THEN
            print *, 'PZHEGV err:info', info, ifail(1), parallel%myid
! print *,'you shi scalapack','amat',amat,'bmat',bmat
            STOP
        END IF
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    END SUBROUTINE SL_GeneralizeEigen_real  !}}}

    SUBROUTINE SL_diagM_real(dime, amat, am, an, bm, bn, amb, anb, bmb, bnb, evec, cm, cn, eval, cmb, cnb)  !{
        USE m_time_evaluate, ONLY: memory_sum, memory_free
        IMPLICIT NONE
        integer(i4b), intent(in) :: dime
        real(DP), INTENT(IN)  :: amat(:, :)
        real(DP), INTENT(OUT) :: evec(:, :)
        REAL(DP), INTENT(OUT) :: eval(:)
        real(DP), allocatable :: U(:, :)
!LOCAL
        integer(i4b) :: am, an, bm, bn, amb, anb, bmb, bnb, cm, cn
        integer(i4b) :: cmb, cnb
        integer(i4b) :: m, nevec
        integer(i4b) :: lwork, liwork
        real(dp), allocatable :: work(:)
        integer(i4b), allocatable :: iwork(:)
!integer(i4b) :: iwkt(10)
!real(dp) :: wkt(10)
        real(dp) :: abtol
!integer(i4b) :: lrowrk
        integer(i4b) :: lda
        integer(i4b) :: maxprocs
        integer(i4b)    :: desca(dlen), descb(dlen), descz(dlen)
! local arrays
!integer(i4b),allocatable :: iclustr(parallel%numprocs*2),ifail(dime)
        integer(i4b) :: iclustr(parallel%numprocs*2), ifail(dime)
        real(dp)     :: w(am)
        real(dp)     :: gap(parallel%numprocs)
        real(dp), external :: pdlamch
        integer(i4b), external :: iceil, pjlaenv
        integer(i4b) :: neig, np0, nq0, nclust, nnp
        integer(i4b) :: nsytrd_lwopt, nsygst_lwopt, a_nb, sqnpc, nps
        INTEGER(I4B) :: i, ia, ja !> For Debug
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        abtol = 1e-6 !2*pdlamch(blacs_contxt,'s')
!> Init the descriptor of the array
        NP = NUMROC(AM, AMB, MYROW, 0, NPROW)
        NQ = NUMROC(AN, ANB, MYCOL, 0, NPCOL)
        CALL DESCINIT(DESCA, AM, AN, AMB, ANB, 0, 0, blacs_contxt, MAX(1, NP), INFO)
        NP = NUMROC(BM, BMB, MYROW, 0, NPROW)
        NQ = NUMROC(BN, BNB, MYCOL, 0, NPCOL)
        CALL DESCINIT(DESCB, BM, BN, BMB, BNB, 0, 0, blacs_contxt, MAX(1, NP), INFO)

        NP = NUMROC(CM, CMB, MYROW, 0, NPROW)
        NQ = NUMROC(CN, CNB, MYCOL, 0, NPCOL)
        CALL DESCINIT(DESCZ, CM, CN, CMB, CNB, 0, 0, blacs_contxt, MAX(1, NP), INFO)
        nnp = max(am, parallel%numprocs + 1, 4)
        neig = an
        np0 = numroc(am, amb, 0, 0, nprow)
        nq0 = max(numroc(neig, anb, 0, 0, npcol), anb)
        lwork = 5*am + max(5*am, np0*nq0 + 2*anb*anb) &
               & + iceil(neig, parallel%numprocs)*am + 2*anb*anb + 4*anb*am
        liwork = 6*am

        a_nb = pjlaenv(blacs_contxt, 3, 'PDSYTTRD', 'L', 0, 0, 0, 0)
        sqnpc = int(sqrt(dble(parallel%numprocs)))
        nps = max(numroc(am, 1, 0, 0, sqnpc), 2*a_nb)
        nsytrd_lwopt = am + 2*(a_nb + 1)*(4*nps + 2) + (nps + 3)*nps
        nsygst_lwopt = 2*np0*anb + nq0*anb + anb*anb
        lwork = max(lwork, 5*am + nsytrd_lwopt, nsygst_lwopt)
        if (allocated(work)) then
            call memory_free('SL_GeneralizeEigen_real_work', real(size(work), DP)*DP)
            deallocate (work)
        end if
        if (allocated(iwork)) then
            call memory_free('SL_diagM_real_iwork', real(size(iwork), DP)*DP)
            deallocate (iwork)
        end if
        allocate (work(lwork))
        allocate (iwork(liwork))
        call memory_sum('SL_diagM_real', real(size(work), DP)*DP + size(iwork)*DP)
!>JOBZ, RANGE, UPLO, N, A, IA, JA, DESCA, VL, VU, IL, IU, ABSTOL, M, NZ, W, ORFAC, Z, IZ, J
        call pdsyevx('V', 'A', 'U', an, amat, 1, 1, desca, 0.d0, 0.d0, &
                    & 1, an, abtol, m, nevec, eval, 1e-12, evec, 1, 1, descz, work, lwork,    &
                    & iwork, liwork, ifail, iclustr, gap, info)

        if (allocated(work)) then
            call memory_free('SL_diagM_real', real(size(work), DP)*DP)
            deallocate (work)
        end if
        if (allocated(iwork)) then
            call memory_free('SL_diagM_real', real(size(iwork), DP)*DP)
            deallocate (iwork)
        end if

        IF (info /= 0) THEN
            print *, 'PZHEGV err:info', info, ifail(1), parallel%myid
            STOP
        END IF
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    END SUBROUTINE SL_diagM_real  !}}}
    SUBROUTINE twoD_map_set(nstates, nrow, ncol&
         &, twoD_map)
        implicit none
        INTEGER(I4B), intent(in) :: nstates, nrow, ncol
        INTEGER(I4B), intent(out) :: twoD_map(3, nrow, ncol)
!> local
        INTEGER(I4B) :: irow, icol
        INTEGER(I4B) :: counter

!>hint: twoD_map(local nrow, local ncol, global id)
        twoD_map = 0
        counter = 0
!> row major
        do irow = 1, nrow, 1
            do icol = 1, ncol, 1
                twoD_map(1, irow, icol) = NUMROC(nstates, 1, irow, 1, nrow)
                twoD_map(2, irow, icol) = NUMROC(nstates, 1, icol, 1, ncol)
                twoD_map(3, irow, icol) = counter
                counter = counter + 1
            end do
        end do

    END SUBROUTINE twoD_map_set
    Subroutine SL_matmat_real_tn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb, cmb, cnb)  !{{{
        USE smpi_math_module, ONLY: MPI_real8, mpinfo, MPI_SUM
        implicit none
        real(dp), intent(in)  :: amat(:, :), bmat(:, :)
        real(dp), intent(out)  :: cmat(:, :)
        character(1), intent(in) :: opa, opb
        integer(i4b)    :: m, n, k
        integer(i4b)    :: am, an, amb, anb, bm, bn, bmb, bnb, cm, cn, cmb, cnb
        integer(i4b)    :: desca(dlen), descb(dlen), descc(dlen)
        real(dp)     :: alpha = 1.d0, &
             & beta = 0.d0
        real(dp), allocatable :: amat_y(:, :), bmat_y(:, :)&
             & , cmat_local(:, :), cmat_global(:, :)
        integer(i4b) :: y_recvcounts(parallel%dims(1))&
             & , y_displs(parallel%dims(1))
        integer(i4b) :: lda, ldb, ldc
        integer(i4b) :: i
        integer(i4B) :: n_local
        integer(i4b) :: x_i, x_j, y_i, y_j, n_block, pad_block, pad_data
        integer(i4b) :: j, i_x(an), i_y(bn), il, im
        !---------------------------------------------------------------------
        !write (*, *) 'line1101'
        m = an
        n = an
        n_local = size(bmat, 2)
        cm = m
        cn = n
        k = size(amat, 1)
        lda = max(size(amat, 1), 1)
        ldb = max(size(bmat, 1), 1)
        ldc = max(m, 1)
        !> shift x
        n_block = cm/cmb/parallel%dims(2)
        pad_block = mod(cm/cmb, parallel%dims(2))
        pad_data = mod(cm, cmb)
        if (parallel%rankx < pad_block) then
            !write (*, *) 'line1116'
            x_i = parallel%rankx*(n_block + 1)*cmb + 1
            !write (*, *) 'line1118'
        elseif (parallel%rankx == pad_block) then
            !write (*, *) 'line1120'
            x_i = parallel%rankx*n_block*cmb + pad_block + pad_data + 1
            !write (*, *) 'line1122'
        else
            !write (*, *) 'line1124'
            x_i = parallel%rankx*n_block*cmb + pad_block + 1
            !write (*, *) 'line1126'
        end if
        !write (*, *) 'line1128'
        x_j = x_i + twoD_map(1, parallel%rankx, parallel%ranky) - 1
        !> shift y
        !write (*, *) 'line1131'
        n_block = cn/cnb/parallel%dims(1)
        pad_block = mod(cn/cnb, parallel%dims(1))
        pad_data = mod(cn, cnb)
        if (parallel%ranky < pad_block) then
            y_i = parallel%ranky*(n_block + 1)*cnb + 1
        elseif (parallel%ranky == pad_block) then
            y_i = parallel%ranky*n_block*cnb + pad_block + pad_data + 1
        else
            y_i = parallel%ranky*n_block*cnb + pad_block + 1
        end if
        y_j = y_i + twoD_map(2, parallel%rankx, parallel%ranky) - 1
        !write (*, *) 'line1143'
        !> gater in row
        allocate (amat_y(size(amat, 1), an))
        allocate (bmat_y(size(amat, 1), an))
        y_recvcounts = twoD_map(2, parallel%rankx, :)*size(amat, 1)
        y_displs(1) = 0
        do i = 2, parallel%dims(1)
            y_displs(i) = sum(y_recvcounts(:i - 1))
        end do
        !write (*, *) 'line1152'
        !xqtest
        !print*,'size amat',size(amat),'rankx',parallel%myid,'rankx',parallel%rankx, &
        ! &'ranky',parallel%ranky,'size in amat_y' &
        ! &,y_recvcounts(parallel%ranky+1),'shift',y_displs(parallel%ranky+1),'size',shape(amat_y)
        !print *,'start rank',parallel%myid,'xy,',parallel%rankx,parallel%ranky
        !write (*, *) 'line1158'
        call MPI_ALLGATHERV(amat, size(amat), MPI_REAL8, amat_y, y_recvcounts, y_displs, MPI_REAL8, parallel%commy, mpinfo)
        !print *,'end rank',parallel%myid,'xy,',parallel%rankx,parallel%ranky
        !      call MPI_BARRIER(parallel%comm,mpinfo)
        !      call MPI_FINALIZE(mpinfo)
        !      stop
        !> local matmat
        !write (*, *) 'line1164'
        allocate (cmat_local(an, an))
        !write (*, *) 'line1167'
        allocate (cmat_global(an, an))
        !write (*, *) 'line1169'
        cmat_local = 0.d0
        CALL DGEMM(opA, opB, m, N_local, K, alpha, amat_y, LDA, bmat, LDB, beta, cmat_local(:, y_i:y_j), LDC)
        !write (*, *) 'line1173'
        !> allreduce the local result
        CALL MPI_ALLREDUCE(Cmat_local, cmat_global, size(cmat_local), MPI_REAL8, MPI_SUM, parallel%comm, mpinfo)
        !write (*, *) 'line1177'
        ! & MPI_SUM, parallel%commx, mpinfo)
        !> remain the local part
        ! Cmat=cmat_global(x_i:x_j,y_i:y_j)
        il = 1; im = 1
        do j = 0, parallel%dims(2) - 1
            do i = 1, an, 1
                if (mod(i - 1, parallel%dims(2)) == j) then
                    i_x(il) = i
                    il = il + 1
                end if
            end do
        end do
        !write (*, *) 'line1191'
        do j = 0, parallel%dims(1) - 1
            do i = 1, an, 1
                if (mod(i - 1, parallel%dims(1)) == j) then
                    i_y(i) = im
                    im = im + 1
                end if
            end do
        end do
        !write (*, *) 'line1200'
        do i = x_i, x_j, 1
            cmat(i + 1 - x_i, :) = cmat_global(i_y(i_x(i)), y_i:y_j)
        end do
        !write (*, *) 'line1204'
        !> test
        ! do j=1,parallel%numprocs,1
        !     write(str_id,'(I4)')j
        !    if(parallel%myid==j-1)then
        !       open(1111+j,file="global_data"//trim(adjustl(str_id)))
        !       ! write(1111+j,*)hhat_local
        !       do i=1,size(cmat_global,2),1
        !          write(1111+j,'(1000(F6.2,2X))')cmat_global(:,i)
        !       enddo
        !       close(1111+j)
        !    endif
        ! enddo
        !>
        deallocate (cmat_local, cmat_global)
        !write (*, *) 'line1219'
    End Subroutine SL_matmat_real_tn  !}}}

    Subroutine SL_matmat_real_nn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb, cmb, cnb)  !{{{
        USE smpi_math_module, ONLY: MPI_real8, mpinfo, MPI_SUM
        implicit none
        real(dp), intent(in)  :: amat(:, :), bmat(:, :)
        real(dp), intent(out)  :: cmat(:, :)
        character(1), intent(in) :: opa, opb
        integer(i4b)    :: m, n, k
        integer(i4b)    :: am, an, amb, anb, bm, bn, bmb, bnb, cm, cn, cmb, cnb
        integer(i4b)    :: desca(dlen), descb(dlen), descc(dlen)
        real(dp)     :: alpha = 1.d0, &
             & beta = 0.d0
        real(dp), allocatable :: Bmat_local(:, :), Bmat_global(:, :), amat_y(:, :)
        integer(i4b) :: lda, ldb, ldc
        integer(i4b) :: i, j, il, y_recvcounts(parallel%dims(1))&
           & , y_displs(parallel%dims(1))

!---------------------------------------------------------------------
!> cmat_local
        allocate (Bmat_local(bm, size(Bmat, 2)))
        allocate (Bmat_global(bm, size(Bmat, 2)))
!> m,n,k
        m = size(amat, 1)
        k = size(bmat_local, 1)
        n = size(bmat_local, 2)
!> lda,ldb,ldc
        lda = max(m, 1)
        ldb = max(k, 1)
        ldc = max(m, 1)

!> assign to local bmat
        il = 0
        Bmat_local = 0.d0
        do i = 1, bm, 1
            if (mod(i - 1, parallel%dims(2)) /= parallel%rankx) cycle
            il = il + 1
            Bmat_local(i, :) = Bmat(il, :)
        end do
        CALL MPI_ALLREDUCE(Bmat_local, Bmat_global, size(Bmat_local), MPI_REAL8,&
             & MPI_SUM, parallel%commx, mpinfo)
!>> reorder the bmat_global
        il = 0
        do j = 0, parallel%dims(1) - 1
            do i = 1, bm, 1
                if (mod(i - 1, parallel%dims(1)) /= j) cycle
                il = il + 1
                Bmat_local(il, :) = Bmat_global(i, :)
            end do
        end do

!> amat_global
        allocate (amat_y(m, k))
        y_recvcounts = twoD_map(2, parallel%rankx, :)*size(amat, 1)
        y_displs(1) = 0
        do i = 2, parallel%dims(1)
            y_displs(i) = sum(y_recvcounts(:i - 1))
        end do
        call MPI_ALLGATHERV(amat, size(amat)&
             &, MPI_REAL8, amat_y, y_recvcounts&
             &, y_displs, MPI_REAL8&
             & , parallel%commy, mpinfo)

!> lapack
        CALL DGEMM(opA, opB, M, N, K, alpha, Amat_y, LDA, Bmat_local, LDB, beta, Cmat, LDC)

!>
        deallocate (Bmat_local, Bmat_global, amat_y)

    End Subroutine SL_matmat_real_nn  !}}}

END MODULE ScaLapack_module

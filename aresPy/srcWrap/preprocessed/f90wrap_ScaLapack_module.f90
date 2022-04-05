! Module scalapack_module defined in file ScaLapack_module.fpp

subroutine f90wrap_init_scala
    use scalapack_module, only: init_scala
    implicit none
    
    call init_scala()
end subroutine f90wrap_init_scala

subroutine f90wrap_sl_orthnorm(amat, m, n, mb, nb, n0, n1)
    use scalapack_module, only: sl_orthnorm
    implicit none
    
    complex(8), dimension(n0,n1) :: amat
    integer(4) :: m
    integer(4) :: n
    integer(4) :: mb
    integer(4) :: nb
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(amat) :: n1 = shape(amat,1)
    call sl_orthnorm(amat=amat, m=m, n=n, mb=mb, nb=nb)
end subroutine f90wrap_sl_orthnorm

subroutine f90wrap_sl_orthnorm_real(amat, m, n, mb, nb, n0, n1)
    use scalapack_module, only: sl_orthnorm_real
    implicit none
    
    real(8), dimension(n0,n1) :: amat
    integer(4) :: m
    integer(4) :: n
    integer(4) :: mb
    integer(4) :: nb
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(amat) :: n1 = shape(amat,1)
    call sl_orthnorm_real(amat=amat, m=m, n=n, mb=mb, nb=nb)
end subroutine f90wrap_sl_orthnorm_real

subroutine f90wrap_sl_matmat(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb, cmbin, cnbin, n0, n1, n2, &
    n3, n4, n5)
    use scalapack_module, only: sl_matmat
    implicit none
    
    character(1), intent(in) :: opa
    character(1), intent(in) :: opb
    complex(8), intent(in), dimension(n0,n1) :: amat
    complex(8), intent(in), dimension(n2,n3) :: bmat
    complex(8), intent(inout), dimension(n4,n5) :: cmat
    integer(4) :: am
    integer(4) :: an
    integer(4) :: bm
    integer(4) :: bn
    integer(4) :: amb
    integer(4) :: anb
    integer(4) :: bmb
    integer(4) :: bnb
    integer(4), optional :: cmbin
    integer(4), optional :: cnbin
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(amat) :: n1 = shape(amat,1)
    integer :: n2
    !f2py intent(hide), depend(bmat) :: n2 = shape(bmat,0)
    integer :: n3
    !f2py intent(hide), depend(bmat) :: n3 = shape(bmat,1)
    integer :: n4
    !f2py intent(hide), depend(cmat) :: n4 = shape(cmat,0)
    integer :: n5
    !f2py intent(hide), depend(cmat) :: n5 = shape(cmat,1)
    call sl_matmat(opa=opa, opb=opb, amat=amat, bmat=bmat, cmat=cmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, &
        bnb=bnb, cmbin=cmbin, cnbin=cnbin)
end subroutine f90wrap_sl_matmat

subroutine f90wrap_sl_matmat_cmplx_cn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb, cmb, cnb, n0, n1, &
    n2, n3, n4, n5)
    use scalapack_module, only: sl_matmat_cmplx_cn
    implicit none
    
    character(1), intent(in) :: opa
    character(1), intent(in) :: opb
    complex(8), intent(in), dimension(n0,n1) :: amat
    complex(8), intent(in), dimension(n2,n3) :: bmat
    complex(8), intent(inout), dimension(n4,n5) :: cmat
    integer(4) :: am
    integer(4) :: an
    integer(4) :: bm
    integer(4) :: bn
    integer(4) :: amb
    integer(4) :: anb
    integer(4) :: bmb
    integer(4) :: bnb
    integer(4) :: cmb
    integer(4) :: cnb
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(amat) :: n1 = shape(amat,1)
    integer :: n2
    !f2py intent(hide), depend(bmat) :: n2 = shape(bmat,0)
    integer :: n3
    !f2py intent(hide), depend(bmat) :: n3 = shape(bmat,1)
    integer :: n4
    !f2py intent(hide), depend(cmat) :: n4 = shape(cmat,0)
    integer :: n5
    !f2py intent(hide), depend(cmat) :: n5 = shape(cmat,1)
    call sl_matmat_cmplx_cn(opa=opa, opb=opb, amat=amat, bmat=bmat, cmat=cmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, &
        bmb=bmb, bnb=bnb, cmb=cmb, cnb=cnb)
end subroutine f90wrap_sl_matmat_cmplx_cn

subroutine f90wrap_sl_matmat_cmplx_nn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb, cmb, cnb, n0, n1, &
    n2, n3, n4, n5)
    use scalapack_module, only: sl_matmat_cmplx_nn
    implicit none
    
    character(1), intent(in) :: opa
    character(1), intent(in) :: opb
    complex(8), intent(in), dimension(n0,n1) :: amat
    complex(8), intent(in), dimension(n2,n3) :: bmat
    complex(8), intent(inout), dimension(n4,n5) :: cmat
    integer(4) :: am
    integer(4) :: an
    integer(4) :: bm
    integer(4) :: bn
    integer(4) :: amb
    integer(4) :: anb
    integer(4) :: bmb
    integer(4) :: bnb
    integer(4) :: cmb
    integer(4) :: cnb
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(amat) :: n1 = shape(amat,1)
    integer :: n2
    !f2py intent(hide), depend(bmat) :: n2 = shape(bmat,0)
    integer :: n3
    !f2py intent(hide), depend(bmat) :: n3 = shape(bmat,1)
    integer :: n4
    !f2py intent(hide), depend(cmat) :: n4 = shape(cmat,0)
    integer :: n5
    !f2py intent(hide), depend(cmat) :: n5 = shape(cmat,1)
    call sl_matmat_cmplx_nn(opa=opa, opb=opb, amat=amat, bmat=bmat, cmat=cmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, &
        bmb=bmb, bnb=bnb, cmb=cmb, cnb=cnb)
end subroutine f90wrap_sl_matmat_cmplx_nn

subroutine f90wrap_sl_matmat_sub(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb, cmbin, cnbin, n0, n1, &
    n2, n3, n4, n5)
    use scalapack_module, only: sl_matmat_sub
    implicit none
    
    character(1), intent(in) :: opa
    character(1), intent(in) :: opb
    complex(8), intent(in), dimension(n0,n1) :: amat
    complex(8), intent(in), dimension(n2,n3) :: bmat
    complex(8), intent(inout), dimension(n4,n5) :: cmat
    integer(4) :: am
    integer(4) :: an
    integer(4) :: bm
    integer(4) :: bn
    integer(4) :: amb
    integer(4) :: anb
    integer(4) :: bmb
    integer(4) :: bnb
    integer(4), optional :: cmbin
    integer(4), optional :: cnbin
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(amat) :: n1 = shape(amat,1)
    integer :: n2
    !f2py intent(hide), depend(bmat) :: n2 = shape(bmat,0)
    integer :: n3
    !f2py intent(hide), depend(bmat) :: n3 = shape(bmat,1)
    integer :: n4
    !f2py intent(hide), depend(cmat) :: n4 = shape(cmat,0)
    integer :: n5
    !f2py intent(hide), depend(cmat) :: n5 = shape(cmat,1)
    call sl_matmat_sub(opa=opa, opb=opb, amat=amat, bmat=bmat, cmat=cmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, &
        bmb=bmb, bnb=bnb, cmbin=cmbin, cnbin=cnbin)
end subroutine f90wrap_sl_matmat_sub

subroutine f90wrap_sl_matmat_sub_real(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb, cmbin, cnbin, n0, &
    n1, n2, n3, n4, n5)
    use scalapack_module, only: sl_matmat_sub_real
    implicit none
    
    character(1), intent(in) :: opa
    character(1), intent(in) :: opb
    real(8), intent(in), dimension(n0,n1) :: amat
    real(8), intent(in), dimension(n2,n3) :: bmat
    real(8), intent(inout), dimension(n4,n5) :: cmat
    integer(4) :: am
    integer(4) :: an
    integer(4) :: bm
    integer(4) :: bn
    integer(4) :: amb
    integer(4) :: anb
    integer(4) :: bmb
    integer(4) :: bnb
    integer(4), optional :: cmbin
    integer(4), optional :: cnbin
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(amat) :: n1 = shape(amat,1)
    integer :: n2
    !f2py intent(hide), depend(bmat) :: n2 = shape(bmat,0)
    integer :: n3
    !f2py intent(hide), depend(bmat) :: n3 = shape(bmat,1)
    integer :: n4
    !f2py intent(hide), depend(cmat) :: n4 = shape(cmat,0)
    integer :: n5
    !f2py intent(hide), depend(cmat) :: n5 = shape(cmat,1)
    call sl_matmat_sub_real(opa=opa, opb=opb, amat=amat, bmat=bmat, cmat=cmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, &
        bmb=bmb, bnb=bnb, cmbin=cmbin, cnbin=cnbin)
end subroutine f90wrap_sl_matmat_sub_real

subroutine f90wrap_sl_matmat_real(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb, cmbin, cnbin, n0, n1, &
    n2, n3, n4, n5)
    use scalapack_module, only: sl_matmat_real
    implicit none
    
    character(1), intent(in) :: opa
    character(1), intent(in) :: opb
    real(8), dimension(n0,n1) :: amat
    real(8), dimension(n2,n3) :: bmat
    real(8), intent(inout), dimension(n4,n5) :: cmat
    integer(4) :: am
    integer(4) :: an
    integer(4) :: bm
    integer(4) :: bn
    integer(4) :: amb
    integer(4) :: anb
    integer(4) :: bmb
    integer(4) :: bnb
    integer(4), optional :: cmbin
    integer(4), optional :: cnbin
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(amat) :: n1 = shape(amat,1)
    integer :: n2
    !f2py intent(hide), depend(bmat) :: n2 = shape(bmat,0)
    integer :: n3
    !f2py intent(hide), depend(bmat) :: n3 = shape(bmat,1)
    integer :: n4
    !f2py intent(hide), depend(cmat) :: n4 = shape(cmat,0)
    integer :: n5
    !f2py intent(hide), depend(cmat) :: n5 = shape(cmat,1)
    call sl_matmat_real(opa=opa, opb=opb, amat=amat, bmat=bmat, cmat=cmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, &
        bmb=bmb, bnb=bnb, cmbin=cmbin, cnbin=cnbin)
end subroutine f90wrap_sl_matmat_real

subroutine f90wrap_sl_generalizeeigen(dime, amat, bmat, am, an, bm, bn, amb, anb, bmb, bnb, evec, cm, cn, eval, cmbin, &
    cnbin, n0, n1, n2, n3, n4, n5, n6)
    use scalapack_module, only: sl_generalizeeigen
    implicit none
    
    integer(4), intent(in) :: dime
    complex(8), intent(in), dimension(n0,n1) :: amat
    complex(8), intent(in), dimension(n2,n3) :: bmat
    integer(4) :: am
    integer(4) :: an
    integer(4) :: bm
    integer(4) :: bn
    integer(4) :: amb
    integer(4) :: anb
    integer(4) :: bmb
    integer(4) :: bnb
    complex(8), intent(inout), dimension(n4,n5) :: evec
    integer(4) :: cm
    integer(4) :: cn
    real(8), intent(inout), dimension(n6) :: eval
    integer(4), optional :: cmbin
    integer(4), optional :: cnbin
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(amat) :: n1 = shape(amat,1)
    integer :: n2
    !f2py intent(hide), depend(bmat) :: n2 = shape(bmat,0)
    integer :: n3
    !f2py intent(hide), depend(bmat) :: n3 = shape(bmat,1)
    integer :: n4
    !f2py intent(hide), depend(evec) :: n4 = shape(evec,0)
    integer :: n5
    !f2py intent(hide), depend(evec) :: n5 = shape(evec,1)
    integer :: n6
    !f2py intent(hide), depend(eval) :: n6 = shape(eval,0)
    call sl_generalizeeigen(dime=dime, amat=amat, bmat=bmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, &
        evec=evec, cm=cm, cn=cn, eval=eval, cmbin=cmbin, cnbin=cnbin)
end subroutine f90wrap_sl_generalizeeigen

subroutine f90wrap_sl_generalizeeigen_real(dime, amat, bmat, am, an, bm, bn, amb, anb, bmb, bnb, evec, cm, cn, eval, &
    cmbin, cnbin, n0, n1, n2, n3, n4, n5, n6)
    use scalapack_module, only: sl_generalizeeigen_real
    implicit none
    
    integer(4), intent(in) :: dime
    real(8), intent(in), dimension(n0,n1) :: amat
    real(8), intent(in), dimension(n2,n3) :: bmat
    integer(4) :: am
    integer(4) :: an
    integer(4) :: bm
    integer(4) :: bn
    integer(4) :: amb
    integer(4) :: anb
    integer(4) :: bmb
    integer(4) :: bnb
    real(8), intent(inout), dimension(n4,n5) :: evec
    integer(4) :: cm
    integer(4) :: cn
    real(8), intent(inout), dimension(n6) :: eval
    integer(4), optional :: cmbin
    integer(4), optional :: cnbin
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(amat) :: n1 = shape(amat,1)
    integer :: n2
    !f2py intent(hide), depend(bmat) :: n2 = shape(bmat,0)
    integer :: n3
    !f2py intent(hide), depend(bmat) :: n3 = shape(bmat,1)
    integer :: n4
    !f2py intent(hide), depend(evec) :: n4 = shape(evec,0)
    integer :: n5
    !f2py intent(hide), depend(evec) :: n5 = shape(evec,1)
    integer :: n6
    !f2py intent(hide), depend(eval) :: n6 = shape(eval,0)
    call sl_generalizeeigen_real(dime=dime, amat=amat, bmat=bmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, &
        bnb=bnb, evec=evec, cm=cm, cn=cn, eval=eval, cmbin=cmbin, cnbin=cnbin)
end subroutine f90wrap_sl_generalizeeigen_real

subroutine f90wrap_sl_diagm_real(dime, amat, am, an, bm, bn, amb, anb, bmb, bnb, evec, cm, cn, eval, cmb, cnb, n0, n1, &
    n2, n3, n4)
    use scalapack_module, only: sl_diagm_real
    implicit none
    
    integer(4), intent(in) :: dime
    real(8), intent(in), dimension(n0,n1) :: amat
    integer(4) :: am
    integer(4) :: an
    integer(4) :: bm
    integer(4) :: bn
    integer(4) :: amb
    integer(4) :: anb
    integer(4) :: bmb
    integer(4) :: bnb
    real(8), intent(inout), dimension(n2,n3) :: evec
    integer(4) :: cm
    integer(4) :: cn
    real(8), intent(inout), dimension(n4) :: eval
    integer(4) :: cmb
    integer(4) :: cnb
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(amat) :: n1 = shape(amat,1)
    integer :: n2
    !f2py intent(hide), depend(evec) :: n2 = shape(evec,0)
    integer :: n3
    !f2py intent(hide), depend(evec) :: n3 = shape(evec,1)
    integer :: n4
    !f2py intent(hide), depend(eval) :: n4 = shape(eval,0)
    call sl_diagm_real(dime=dime, amat=amat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, evec=evec, &
        cm=cm, cn=cn, eval=eval, cmb=cmb, cnb=cnb)
end subroutine f90wrap_sl_diagm_real

subroutine f90wrap_twod_map_set(nstates, nrow, ncol, twod_map, n0, n1)
    use scalapack_module, only: twod_map_set
    implicit none
    
    integer(4), intent(in) :: nstates
    integer(4), intent(in) :: nrow
    integer(4), intent(in) :: ncol
    integer(4), intent(inout), dimension(3,n0,n1) :: twod_map
    integer :: n0
    !f2py intent(hide), depend(twod_map) :: n0 = shape(twod_map,1)
    integer :: n1
    !f2py intent(hide), depend(twod_map) :: n1 = shape(twod_map,2)
    call twod_map_set(nstates=nstates, nrow=nrow, ncol=ncol, twoD_map=twod_map)
end subroutine f90wrap_twod_map_set

subroutine f90wrap_sl_matmat_real_tn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb, cmb, cnb, n0, n1, &
    n2, n3, n4, n5)
    use scalapack_module, only: sl_matmat_real_tn
    implicit none
    
    character(1), intent(in) :: opa
    character(1), intent(in) :: opb
    real(8), intent(in), dimension(n0,n1) :: amat
    real(8), intent(in), dimension(n2,n3) :: bmat
    real(8), intent(inout), dimension(n4,n5) :: cmat
    integer(4) :: am
    integer(4) :: an
    integer(4) :: bm
    integer(4) :: bn
    integer(4) :: amb
    integer(4) :: anb
    integer(4) :: bmb
    integer(4) :: bnb
    integer(4) :: cmb
    integer(4) :: cnb
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(amat) :: n1 = shape(amat,1)
    integer :: n2
    !f2py intent(hide), depend(bmat) :: n2 = shape(bmat,0)
    integer :: n3
    !f2py intent(hide), depend(bmat) :: n3 = shape(bmat,1)
    integer :: n4
    !f2py intent(hide), depend(cmat) :: n4 = shape(cmat,0)
    integer :: n5
    !f2py intent(hide), depend(cmat) :: n5 = shape(cmat,1)
    call sl_matmat_real_tn(opa=opa, opb=opb, amat=amat, bmat=bmat, cmat=cmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, &
        bmb=bmb, bnb=bnb, cmb=cmb, cnb=cnb)
end subroutine f90wrap_sl_matmat_real_tn

subroutine f90wrap_sl_matmat_real_nn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb, cmb, cnb, n0, n1, &
    n2, n3, n4, n5)
    use scalapack_module, only: sl_matmat_real_nn
    implicit none
    
    character(1), intent(in) :: opa
    character(1), intent(in) :: opb
    real(8), intent(in), dimension(n0,n1) :: amat
    real(8), intent(in), dimension(n2,n3) :: bmat
    real(8), intent(inout), dimension(n4,n5) :: cmat
    integer(4) :: am
    integer(4) :: an
    integer(4) :: bm
    integer(4) :: bn
    integer(4) :: amb
    integer(4) :: anb
    integer(4) :: bmb
    integer(4) :: bnb
    integer(4) :: cmb
    integer(4) :: cnb
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(amat) :: n1 = shape(amat,1)
    integer :: n2
    !f2py intent(hide), depend(bmat) :: n2 = shape(bmat,0)
    integer :: n3
    !f2py intent(hide), depend(bmat) :: n3 = shape(bmat,1)
    integer :: n4
    !f2py intent(hide), depend(cmat) :: n4 = shape(cmat,0)
    integer :: n5
    !f2py intent(hide), depend(cmat) :: n5 = shape(cmat,1)
    call sl_matmat_real_nn(opa=opa, opb=opb, amat=amat, bmat=bmat, cmat=cmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, &
        bmb=bmb, bnb=bnb, cmb=cmb, cnb=cnb)
end subroutine f90wrap_sl_matmat_real_nn

subroutine f90wrap_scalapack_module__get__blacs_contxt(f90wrap_blacs_contxt)
    use scalapack_module, only: scalapack_module_blacs_contxt => blacs_contxt
    implicit none
    integer(4), intent(out) :: f90wrap_blacs_contxt
    
    f90wrap_blacs_contxt = scalapack_module_blacs_contxt
end subroutine f90wrap_scalapack_module__get__blacs_contxt

subroutine f90wrap_scalapack_module__set__blacs_contxt(f90wrap_blacs_contxt)
    use scalapack_module, only: scalapack_module_blacs_contxt => blacs_contxt
    implicit none
    integer(4), intent(in) :: f90wrap_blacs_contxt
    
    scalapack_module_blacs_contxt = f90wrap_blacs_contxt
end subroutine f90wrap_scalapack_module__set__blacs_contxt

subroutine f90wrap_scalapack_module__get__DLEN(f90wrap_DLEN)
    use scalapack_module, only: scalapack_module_DLEN => DLEN
    implicit none
    integer(4), intent(out) :: f90wrap_DLEN
    
    f90wrap_DLEN = scalapack_module_DLEN
end subroutine f90wrap_scalapack_module__get__DLEN

subroutine f90wrap_scalapack_module__get__MYROW(f90wrap_MYROW)
    use scalapack_module, only: scalapack_module_MYROW => MYROW
    implicit none
    integer(4), intent(out) :: f90wrap_MYROW
    
    f90wrap_MYROW = scalapack_module_MYROW
end subroutine f90wrap_scalapack_module__get__MYROW

subroutine f90wrap_scalapack_module__set__MYROW(f90wrap_MYROW)
    use scalapack_module, only: scalapack_module_MYROW => MYROW
    implicit none
    integer(4), intent(in) :: f90wrap_MYROW
    
    scalapack_module_MYROW = f90wrap_MYROW
end subroutine f90wrap_scalapack_module__set__MYROW

subroutine f90wrap_scalapack_module__get__MYCOL(f90wrap_MYCOL)
    use scalapack_module, only: scalapack_module_MYCOL => MYCOL
    implicit none
    integer(4), intent(out) :: f90wrap_MYCOL
    
    f90wrap_MYCOL = scalapack_module_MYCOL
end subroutine f90wrap_scalapack_module__get__MYCOL

subroutine f90wrap_scalapack_module__set__MYCOL(f90wrap_MYCOL)
    use scalapack_module, only: scalapack_module_MYCOL => MYCOL
    implicit none
    integer(4), intent(in) :: f90wrap_MYCOL
    
    scalapack_module_MYCOL = f90wrap_MYCOL
end subroutine f90wrap_scalapack_module__set__MYCOL

subroutine f90wrap_scalapack_module__get__NPCOL(f90wrap_NPCOL)
    use scalapack_module, only: scalapack_module_NPCOL => NPCOL
    implicit none
    integer(4), intent(out) :: f90wrap_NPCOL
    
    f90wrap_NPCOL = scalapack_module_NPCOL
end subroutine f90wrap_scalapack_module__get__NPCOL

subroutine f90wrap_scalapack_module__set__NPCOL(f90wrap_NPCOL)
    use scalapack_module, only: scalapack_module_NPCOL => NPCOL
    implicit none
    integer(4), intent(in) :: f90wrap_NPCOL
    
    scalapack_module_NPCOL = f90wrap_NPCOL
end subroutine f90wrap_scalapack_module__set__NPCOL

subroutine f90wrap_scalapack_module__get__NPROW(f90wrap_NPROW)
    use scalapack_module, only: scalapack_module_NPROW => NPROW
    implicit none
    integer(4), intent(out) :: f90wrap_NPROW
    
    f90wrap_NPROW = scalapack_module_NPROW
end subroutine f90wrap_scalapack_module__get__NPROW

subroutine f90wrap_scalapack_module__set__NPROW(f90wrap_NPROW)
    use scalapack_module, only: scalapack_module_NPROW => NPROW
    implicit none
    integer(4), intent(in) :: f90wrap_NPROW
    
    scalapack_module_NPROW = f90wrap_NPROW
end subroutine f90wrap_scalapack_module__set__NPROW

subroutine f90wrap_scalapack_module__get__my_blacs_id(f90wrap_my_blacs_id)
    use scalapack_module, only: scalapack_module_my_blacs_id => my_blacs_id
    implicit none
    integer(4), intent(out) :: f90wrap_my_blacs_id
    
    f90wrap_my_blacs_id = scalapack_module_my_blacs_id
end subroutine f90wrap_scalapack_module__get__my_blacs_id

subroutine f90wrap_scalapack_module__set__my_blacs_id(f90wrap_my_blacs_id)
    use scalapack_module, only: scalapack_module_my_blacs_id => my_blacs_id
    implicit none
    integer(4), intent(in) :: f90wrap_my_blacs_id
    
    scalapack_module_my_blacs_id = f90wrap_my_blacs_id
end subroutine f90wrap_scalapack_module__set__my_blacs_id

subroutine f90wrap_scalapack_module__get__NP(f90wrap_NP)
    use scalapack_module, only: scalapack_module_NP => NP
    implicit none
    integer(4), intent(out) :: f90wrap_NP
    
    f90wrap_NP = scalapack_module_NP
end subroutine f90wrap_scalapack_module__get__NP

subroutine f90wrap_scalapack_module__set__NP(f90wrap_NP)
    use scalapack_module, only: scalapack_module_NP => NP
    implicit none
    integer(4), intent(in) :: f90wrap_NP
    
    scalapack_module_NP = f90wrap_NP
end subroutine f90wrap_scalapack_module__set__NP

subroutine f90wrap_scalapack_module__get__NQ(f90wrap_NQ)
    use scalapack_module, only: scalapack_module_NQ => NQ
    implicit none
    integer(4), intent(out) :: f90wrap_NQ
    
    f90wrap_NQ = scalapack_module_NQ
end subroutine f90wrap_scalapack_module__get__NQ

subroutine f90wrap_scalapack_module__set__NQ(f90wrap_NQ)
    use scalapack_module, only: scalapack_module_NQ => NQ
    implicit none
    integer(4), intent(in) :: f90wrap_NQ
    
    scalapack_module_NQ = f90wrap_NQ
end subroutine f90wrap_scalapack_module__set__NQ

subroutine f90wrap_scalapack_module__get__IAM(f90wrap_IAM)
    use scalapack_module, only: scalapack_module_IAM => IAM
    implicit none
    integer(4), intent(out) :: f90wrap_IAM
    
    f90wrap_IAM = scalapack_module_IAM
end subroutine f90wrap_scalapack_module__get__IAM

subroutine f90wrap_scalapack_module__set__IAM(f90wrap_IAM)
    use scalapack_module, only: scalapack_module_IAM => IAM
    implicit none
    integer(4), intent(in) :: f90wrap_IAM
    
    scalapack_module_IAM = f90wrap_IAM
end subroutine f90wrap_scalapack_module__set__IAM

subroutine f90wrap_scalapack_module__get__NPROCS(f90wrap_NPROCS)
    use scalapack_module, only: scalapack_module_NPROCS => NPROCS
    implicit none
    integer(4), intent(out) :: f90wrap_NPROCS
    
    f90wrap_NPROCS = scalapack_module_NPROCS
end subroutine f90wrap_scalapack_module__get__NPROCS

subroutine f90wrap_scalapack_module__set__NPROCS(f90wrap_NPROCS)
    use scalapack_module, only: scalapack_module_NPROCS => NPROCS
    implicit none
    integer(4), intent(in) :: f90wrap_NPROCS
    
    scalapack_module_NPROCS = f90wrap_NPROCS
end subroutine f90wrap_scalapack_module__set__NPROCS

subroutine f90wrap_scalapack_module__get__INIT_CALLED(f90wrap_INIT_CALLED)
    use scalapack_module, only: scalapack_module_INIT_CALLED => INIT_CALLED
    implicit none
    logical, intent(out) :: f90wrap_INIT_CALLED
    
    f90wrap_INIT_CALLED = scalapack_module_INIT_CALLED
end subroutine f90wrap_scalapack_module__get__INIT_CALLED

subroutine f90wrap_scalapack_module__set__INIT_CALLED(f90wrap_INIT_CALLED)
    use scalapack_module, only: scalapack_module_INIT_CALLED => INIT_CALLED
    implicit none
    logical, intent(in) :: f90wrap_INIT_CALLED
    
    scalapack_module_INIT_CALLED = f90wrap_INIT_CALLED
end subroutine f90wrap_scalapack_module__set__INIT_CALLED

subroutine f90wrap_scalapack_module__get__INFO(f90wrap_INFO)
    use scalapack_module, only: scalapack_module_INFO => INFO
    implicit none
    integer(4), intent(out) :: f90wrap_INFO
    
    f90wrap_INFO = scalapack_module_INFO
end subroutine f90wrap_scalapack_module__get__INFO

subroutine f90wrap_scalapack_module__set__INFO(f90wrap_INFO)
    use scalapack_module, only: scalapack_module_INFO => INFO
    implicit none
    integer(4), intent(in) :: f90wrap_INFO
    
    scalapack_module_INFO = f90wrap_INFO
end subroutine f90wrap_scalapack_module__set__INFO

subroutine f90wrap_scalapack_module__get__L_useless(f90wrap_L_useless)
    use scalapack_module, only: scalapack_module_L_useless => L_useless
    implicit none
    logical, intent(out) :: f90wrap_L_useless
    
    f90wrap_L_useless = scalapack_module_L_useless
end subroutine f90wrap_scalapack_module__get__L_useless

subroutine f90wrap_scalapack_module__set__L_useless(f90wrap_L_useless)
    use scalapack_module, only: scalapack_module_L_useless => L_useless
    implicit none
    logical, intent(in) :: f90wrap_L_useless
    
    scalapack_module_L_useless = f90wrap_L_useless
end subroutine f90wrap_scalapack_module__set__L_useless

subroutine f90wrap_scalapack_module__array__twoD_map(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use smpi_math_module
    use scalapack_module, only: scalapack_module_twod_map => twod_map
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 5
    if (allocated(scalapack_module_twoD_map)) then
        dshape(1:3) = shape(scalapack_module_twoD_map)
        dloc = loc(scalapack_module_twoD_map)
    else
        dloc = 0
    end if
end subroutine f90wrap_scalapack_module__array__twoD_map

! End of module scalapack_module defined in file ScaLapack_module.fpp


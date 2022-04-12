! Module finite_module defined in file Finite_module.fpp

subroutine f90wrap_init_finite(norder, h)
    use finite_module, only: init_finite
    implicit none
    
    integer(4) :: norder
    real(8), dimension(3) :: h
    call init_finite(norder=norder, h=h)
end subroutine f90wrap_init_finite

subroutine f90wrap_destroy_finite
    use finite_module, only: destroy_finite
    implicit none
    
    call destroy_finite()
end subroutine f90wrap_destroy_finite

subroutine f90wrap_trans_mat_full(mat, factor, int_miu, lad_gap, err)
    use finite_module, only: trans_mat_full
    implicit none
    
    real(8), intent(in), dimension(3,3) :: mat
    real(8), dimension(6), intent(inout) :: factor
    integer(4), dimension(3,3), intent(inout) :: int_miu
    real(8), dimension(3), intent(inout) :: lad_gap
    real(8), dimension(6), intent(inout) :: err
    call trans_mat_full(mat=mat, factor=factor, int_miu=int_miu, lad_gap=lad_gap, err=err)
end subroutine f90wrap_trans_mat_full

subroutine f90wrap_cmplx_nabla2(ifun, norder, ofun, n0, n1, n2, n3, n4, n5)
    use finite_module, only: cmplx_nabla2
    implicit none
    
    complex(8), intent(in), dimension(n0,n1,n2) :: ifun
    integer(4), intent(in) :: norder
    complex(8), intent(inout), dimension(n3,n4,n5) :: ofun
    integer :: n0
    !f2py intent(hide), depend(ifun) :: n0 = shape(ifun,0)
    integer :: n1
    !f2py intent(hide), depend(ifun) :: n1 = shape(ifun,1)
    integer :: n2
    !f2py intent(hide), depend(ifun) :: n2 = shape(ifun,2)
    integer :: n3
    !f2py intent(hide), depend(ofun) :: n3 = shape(ofun,0)
    integer :: n4
    !f2py intent(hide), depend(ofun) :: n4 = shape(ofun,1)
    integer :: n5
    !f2py intent(hide), depend(ofun) :: n5 = shape(ofun,2)
    call cmplx_nabla2(ifun=ifun, norder=norder, ofun=ofun)
end subroutine f90wrap_cmplx_nabla2

subroutine f90wrap_cmplx_nabla1(ifun, norder, derf, n0, n1, n2, n3, n4, n5, n6)
    use finite_module, only: cmplx_nabla1
    implicit none
    
    complex(8), intent(in), dimension(n0,n1,n2) :: ifun
    integer(4), intent(in) :: norder
    complex(8), intent(inout), dimension(n3,n4,n5,n6) :: derf
    integer :: n0
    !f2py intent(hide), depend(ifun) :: n0 = shape(ifun,0)
    integer :: n1
    !f2py intent(hide), depend(ifun) :: n1 = shape(ifun,1)
    integer :: n2
    !f2py intent(hide), depend(ifun) :: n2 = shape(ifun,2)
    integer :: n3
    !f2py intent(hide), depend(derf) :: n3 = shape(derf,0)
    integer :: n4
    !f2py intent(hide), depend(derf) :: n4 = shape(derf,1)
    integer :: n5
    !f2py intent(hide), depend(derf) :: n5 = shape(derf,2)
    integer :: n6
    !f2py intent(hide), depend(derf) :: n6 = shape(derf,3)
    call cmplx_nabla1(ifun=ifun, norder=norder, derf=derf)
end subroutine f90wrap_cmplx_nabla1

subroutine f90wrap_nabla(ifun, norder, derf, n0, n1, n2, n3, n4, n5, n6)
    use finite_module, only: nabla
    implicit none
    
    complex(8), intent(in), dimension(n0,n1,n2) :: ifun
    integer(4), intent(in) :: norder
    complex(8), intent(inout), dimension(n3,n4,n5,n6) :: derf
    integer :: n0
    !f2py intent(hide), depend(ifun) :: n0 = shape(ifun,0)
    integer :: n1
    !f2py intent(hide), depend(ifun) :: n1 = shape(ifun,1)
    integer :: n2
    !f2py intent(hide), depend(ifun) :: n2 = shape(ifun,2)
    integer :: n3
    !f2py intent(hide), depend(derf) :: n3 = shape(derf,0)
    integer :: n4
    !f2py intent(hide), depend(derf) :: n4 = shape(derf,1)
    integer :: n5
    !f2py intent(hide), depend(derf) :: n5 = shape(derf,2)
    integer :: n6
    !f2py intent(hide), depend(derf) :: n6 = shape(derf,3)
    call nabla(ifun=ifun, norder=norder, derf=derf)
end subroutine f90wrap_nabla

subroutine f90wrap_kskeperiod_op(psi, ik, kepsi, n0, n1, n2, n3, n4, n5)
    use finite_module, only: kskeperiod_op
    implicit none
    
    complex(8), intent(in), dimension(n0,n1,n2) :: psi
    integer(4), intent(in) :: ik
    complex(8), intent(inout), dimension(n3,n4,n5) :: kepsi
    integer :: n0
    !f2py intent(hide), depend(psi) :: n0 = shape(psi,0)
    integer :: n1
    !f2py intent(hide), depend(psi) :: n1 = shape(psi,1)
    integer :: n2
    !f2py intent(hide), depend(psi) :: n2 = shape(psi,2)
    integer :: n3
    !f2py intent(hide), depend(kepsi) :: n3 = shape(kepsi,0)
    integer :: n4
    !f2py intent(hide), depend(kepsi) :: n4 = shape(kepsi,1)
    integer :: n5
    !f2py intent(hide), depend(kepsi) :: n5 = shape(kepsi,2)
    call kskeperiod_op(psi=psi, Ik=ik, KEpsi=kepsi)
end subroutine f90wrap_kskeperiod_op

subroutine f90wrap_ke1(psi, ts1, n0, n1, n2, n3, n4, n5)
    use finite_module, only: ke1
    implicit none
    
    complex(8), intent(in), dimension(n0,n1,n2) :: psi
    complex(8), intent(inout), dimension(n3,n4,n5) :: ts1
    integer :: n0
    !f2py intent(hide), depend(psi) :: n0 = shape(psi,0)
    integer :: n1
    !f2py intent(hide), depend(psi) :: n1 = shape(psi,1)
    integer :: n2
    !f2py intent(hide), depend(psi) :: n2 = shape(psi,2)
    integer :: n3
    !f2py intent(hide), depend(ts1) :: n3 = shape(ts1,0)
    integer :: n4
    !f2py intent(hide), depend(ts1) :: n4 = shape(ts1,1)
    integer :: n5
    !f2py intent(hide), depend(ts1) :: n5 = shape(ts1,2)
    call ke1(psi=psi, Ts1=ts1)
end subroutine f90wrap_ke1

subroutine f90wrap_ke2(psi, ik, ts2, n0, n1, n2, n3, n4, n5)
    use finite_module, only: ke2
    implicit none
    
    complex(8), intent(in), dimension(n0,n1,n2) :: psi
    integer(4), intent(in) :: ik
    complex(8), intent(inout), dimension(n3,n4,n5) :: ts2
    integer :: n0
    !f2py intent(hide), depend(psi) :: n0 = shape(psi,0)
    integer :: n1
    !f2py intent(hide), depend(psi) :: n1 = shape(psi,1)
    integer :: n2
    !f2py intent(hide), depend(psi) :: n2 = shape(psi,2)
    integer :: n3
    !f2py intent(hide), depend(ts2) :: n3 = shape(ts2,0)
    integer :: n4
    !f2py intent(hide), depend(ts2) :: n4 = shape(ts2,1)
    integer :: n5
    !f2py intent(hide), depend(ts2) :: n5 = shape(ts2,2)
    call ke2(psi=psi, Ik=ik, Ts2=ts2)
end subroutine f90wrap_ke2

subroutine f90wrap_realnabla2(ifun, norder, ofun, n0, n1, n2, n3, n4, n5)
    use finite_module, only: realnabla2
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: ifun
    integer(4), intent(in) :: norder
    real(8), intent(inout), dimension(n3,n4,n5) :: ofun
    integer :: n0
    !f2py intent(hide), depend(ifun) :: n0 = shape(ifun,0)
    integer :: n1
    !f2py intent(hide), depend(ifun) :: n1 = shape(ifun,1)
    integer :: n2
    !f2py intent(hide), depend(ifun) :: n2 = shape(ifun,2)
    integer :: n3
    !f2py intent(hide), depend(ofun) :: n3 = shape(ofun,0)
    integer :: n4
    !f2py intent(hide), depend(ofun) :: n4 = shape(ofun,1)
    integer :: n5
    !f2py intent(hide), depend(ofun) :: n5 = shape(ofun,2)
    call realnabla2(ifun=ifun, norder=norder, ofun=ofun)
end subroutine f90wrap_realnabla2

subroutine f90wrap_realke(psi, kepsi, n0, n1, n2, n3, n4, n5)
    use finite_module, only: realke
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: psi
    real(8), intent(inout), dimension(n3,n4,n5) :: kepsi
    integer :: n0
    !f2py intent(hide), depend(psi) :: n0 = shape(psi,0)
    integer :: n1
    !f2py intent(hide), depend(psi) :: n1 = shape(psi,1)
    integer :: n2
    !f2py intent(hide), depend(psi) :: n2 = shape(psi,2)
    integer :: n3
    !f2py intent(hide), depend(kepsi) :: n3 = shape(kepsi,0)
    integer :: n4
    !f2py intent(hide), depend(kepsi) :: n4 = shape(kepsi,1)
    integer :: n5
    !f2py intent(hide), depend(kepsi) :: n5 = shape(kepsi,2)
    call realke(psi=psi, KEpsi=kepsi)
end subroutine f90wrap_realke

subroutine f90wrap_iso_realnabla2(ifun, norder, ofun, n0, n1)
    use finite_module, only: iso_realnabla2
    implicit none
    
    real(8), intent(in), dimension(n0) :: ifun
    integer(4), intent(in) :: norder
    real(8), intent(inout), dimension(n1) :: ofun
    integer :: n0
    !f2py intent(hide), depend(ifun) :: n0 = shape(ifun,0)
    integer :: n1
    !f2py intent(hide), depend(ofun) :: n1 = shape(ofun,0)
    call iso_realnabla2(ifun=ifun, norder=norder, ofun=ofun)
end subroutine f90wrap_iso_realnabla2

subroutine f90wrap_iso_realke(psi, kepsi, n0, n1)
    use finite_module, only: iso_realke
    implicit none
    
    real(8), intent(in), dimension(n0) :: psi
    real(8), intent(inout), dimension(n1) :: kepsi
    integer :: n0
    !f2py intent(hide), depend(psi) :: n0 = shape(psi,0)
    integer :: n1
    !f2py intent(hide), depend(kepsi) :: n1 = shape(kepsi,0)
    call iso_realke(psi=psi, KEpsi=kepsi)
end subroutine f90wrap_iso_realke

subroutine f90wrap_kskeperiod_op_band(psi, ik, kepsi, n0, n1, n2, n3, n4, n5)
    use finite_module, only: kskeperiod_op_band
    implicit none
    
    complex(8), intent(in), dimension(n0,n1,n2) :: psi
    integer(4), intent(in) :: ik
    complex(8), intent(inout), dimension(n3,n4,n5) :: kepsi
    integer :: n0
    !f2py intent(hide), depend(psi) :: n0 = shape(psi,0)
    integer :: n1
    !f2py intent(hide), depend(psi) :: n1 = shape(psi,1)
    integer :: n2
    !f2py intent(hide), depend(psi) :: n2 = shape(psi,2)
    integer :: n3
    !f2py intent(hide), depend(kepsi) :: n3 = shape(kepsi,0)
    integer :: n4
    !f2py intent(hide), depend(kepsi) :: n4 = shape(kepsi,1)
    integer :: n5
    !f2py intent(hide), depend(kepsi) :: n5 = shape(kepsi,2)
    call kskeperiod_op_band(psi=psi, Ik=ik, KEpsi=kepsi)
end subroutine f90wrap_kskeperiod_op_band

subroutine f90wrap_nabla2_np(ifun, norder, ofun, n0, n1, n2, n3, n4, n5)
    use finite_module, only: nabla2_np
    implicit none
    
    complex(8), intent(in), dimension(n0,n1,n2) :: ifun
    integer(4), intent(in) :: norder
    complex(8), intent(inout), dimension(n3,n4,n5) :: ofun
    integer :: n0
    !f2py intent(hide), depend(ifun) :: n0 = shape(ifun,0)
    integer :: n1
    !f2py intent(hide), depend(ifun) :: n1 = shape(ifun,1)
    integer :: n2
    !f2py intent(hide), depend(ifun) :: n2 = shape(ifun,2)
    integer :: n3
    !f2py intent(hide), depend(ofun) :: n3 = shape(ofun,0)
    integer :: n4
    !f2py intent(hide), depend(ofun) :: n4 = shape(ofun,1)
    integer :: n5
    !f2py intent(hide), depend(ofun) :: n5 = shape(ofun,2)
    call nabla2_np(ifun=ifun, norder=norder, ofun=ofun)
end subroutine f90wrap_nabla2_np

subroutine f90wrap_nabla1_np(ifun, norder, derf, n0, n1, n2, n3, n4, n5, n6)
    use finite_module, only: nabla1_np
    implicit none
    
    complex(8), intent(in), dimension(n0,n1,n2) :: ifun
    integer(4), intent(in) :: norder
    complex(8), intent(inout), dimension(n3,n4,n5,n6) :: derf
    integer :: n0
    !f2py intent(hide), depend(ifun) :: n0 = shape(ifun,0)
    integer :: n1
    !f2py intent(hide), depend(ifun) :: n1 = shape(ifun,1)
    integer :: n2
    !f2py intent(hide), depend(ifun) :: n2 = shape(ifun,2)
    integer :: n3
    !f2py intent(hide), depend(derf) :: n3 = shape(derf,0)
    integer :: n4
    !f2py intent(hide), depend(derf) :: n4 = shape(derf,1)
    integer :: n5
    !f2py intent(hide), depend(derf) :: n5 = shape(derf,2)
    integer :: n6
    !f2py intent(hide), depend(derf) :: n6 = shape(derf,3)
    call nabla1_np(ifun=ifun, norder=norder, derf=derf)
end subroutine f90wrap_nabla1_np

subroutine f90wrap_ke1_np(psi, ts1, n0, n1, n2, n3, n4, n5)
    use finite_module, only: ke1_np
    implicit none
    
    complex(8), intent(in), dimension(n0,n1,n2) :: psi
    complex(8), intent(inout), dimension(n3,n4,n5) :: ts1
    integer :: n0
    !f2py intent(hide), depend(psi) :: n0 = shape(psi,0)
    integer :: n1
    !f2py intent(hide), depend(psi) :: n1 = shape(psi,1)
    integer :: n2
    !f2py intent(hide), depend(psi) :: n2 = shape(psi,2)
    integer :: n3
    !f2py intent(hide), depend(ts1) :: n3 = shape(ts1,0)
    integer :: n4
    !f2py intent(hide), depend(ts1) :: n4 = shape(ts1,1)
    integer :: n5
    !f2py intent(hide), depend(ts1) :: n5 = shape(ts1,2)
    call ke1_np(psi=psi, Ts1=ts1)
end subroutine f90wrap_ke1_np

subroutine f90wrap_ke2_np(psi, ik, ts2, n0, n1, n2, n3, n4, n5)
    use finite_module, only: ke2_np
    implicit none
    
    complex(8), intent(in), dimension(n0,n1,n2) :: psi
    integer(4), intent(in) :: ik
    complex(8), intent(inout), dimension(n3,n4,n5) :: ts2
    integer :: n0
    !f2py intent(hide), depend(psi) :: n0 = shape(psi,0)
    integer :: n1
    !f2py intent(hide), depend(psi) :: n1 = shape(psi,1)
    integer :: n2
    !f2py intent(hide), depend(psi) :: n2 = shape(psi,2)
    integer :: n3
    !f2py intent(hide), depend(ts2) :: n3 = shape(ts2,0)
    integer :: n4
    !f2py intent(hide), depend(ts2) :: n4 = shape(ts2,1)
    integer :: n5
    !f2py intent(hide), depend(ts2) :: n5 = shape(ts2,2)
    call ke2_np(psi=psi, Ik=ik, Ts2=ts2)
end subroutine f90wrap_ke2_np

subroutine f90wrap_kskeperiod_op_gamma(psi, ik, kepsi, n0, n1, n2, n3, n4, n5)
    use finite_module, only: kskeperiod_op_gamma
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: psi
    integer(4), intent(in) :: ik
    real(8), intent(inout), dimension(n3,n4,n5) :: kepsi
    integer :: n0
    !f2py intent(hide), depend(psi) :: n0 = shape(psi,0)
    integer :: n1
    !f2py intent(hide), depend(psi) :: n1 = shape(psi,1)
    integer :: n2
    !f2py intent(hide), depend(psi) :: n2 = shape(psi,2)
    integer :: n3
    !f2py intent(hide), depend(kepsi) :: n3 = shape(kepsi,0)
    integer :: n4
    !f2py intent(hide), depend(kepsi) :: n4 = shape(kepsi,1)
    integer :: n5
    !f2py intent(hide), depend(kepsi) :: n5 = shape(kepsi,2)
    call kskeperiod_op_gamma(psi=psi, Ik=ik, KEpsi=kepsi)
end subroutine f90wrap_kskeperiod_op_gamma

subroutine f90wrap_ke_gamma(psi, ts1, n0, n1, n2, n3, n4, n5)
    use finite_module, only: ke_gamma
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: psi
    real(8), intent(inout), dimension(n3,n4,n5) :: ts1
    integer :: n0
    !f2py intent(hide), depend(psi) :: n0 = shape(psi,0)
    integer :: n1
    !f2py intent(hide), depend(psi) :: n1 = shape(psi,1)
    integer :: n2
    !f2py intent(hide), depend(psi) :: n2 = shape(psi,2)
    integer :: n3
    !f2py intent(hide), depend(ts1) :: n3 = shape(ts1,0)
    integer :: n4
    !f2py intent(hide), depend(ts1) :: n4 = shape(ts1,1)
    integer :: n5
    !f2py intent(hide), depend(ts1) :: n5 = shape(ts1,2)
    call ke_gamma(psi=psi, Ts1=ts1)
end subroutine f90wrap_ke_gamma

subroutine f90wrap_real_nabla2(ifun, norder, ofun, n0, n1, n2, n3, n4, n5)
    use finite_module, only: real_nabla2
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: ifun
    integer(4), intent(in) :: norder
    real(8), intent(inout), dimension(n3,n4,n5) :: ofun
    integer :: n0
    !f2py intent(hide), depend(ifun) :: n0 = shape(ifun,0)
    integer :: n1
    !f2py intent(hide), depend(ifun) :: n1 = shape(ifun,1)
    integer :: n2
    !f2py intent(hide), depend(ifun) :: n2 = shape(ifun,2)
    integer :: n3
    !f2py intent(hide), depend(ofun) :: n3 = shape(ofun,0)
    integer :: n4
    !f2py intent(hide), depend(ofun) :: n4 = shape(ofun,1)
    integer :: n5
    !f2py intent(hide), depend(ofun) :: n5 = shape(ofun,2)
    call real_nabla2(ifun=ifun, norder=norder, ofun=ofun)
end subroutine f90wrap_real_nabla2

subroutine f90wrap_real_nabla2_1d(ifun, ofun, n0, n1)
    use finite_module, only: real_nabla2_1d
    implicit none
    
    real(8), intent(in), dimension(n0) :: ifun
    real(8), intent(inout), dimension(n1) :: ofun
    integer :: n0
    !f2py intent(hide), depend(ifun) :: n0 = shape(ifun,0)
    integer :: n1
    !f2py intent(hide), depend(ofun) :: n1 = shape(ofun,0)
    call real_nabla2_1d(ifun=ifun, ofun=ofun)
end subroutine f90wrap_real_nabla2_1d

subroutine f90wrap_nabla_gamma(ifun, norder, derf, n0, n1, n2, n3, n4, n5, n6)
    use finite_module, only: nabla_gamma
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: ifun
    integer(4), intent(in) :: norder
    real(8), intent(inout), dimension(n3,n4,n5,n6) :: derf
    integer :: n0
    !f2py intent(hide), depend(ifun) :: n0 = shape(ifun,0)
    integer :: n1
    !f2py intent(hide), depend(ifun) :: n1 = shape(ifun,1)
    integer :: n2
    !f2py intent(hide), depend(ifun) :: n2 = shape(ifun,2)
    integer :: n3
    !f2py intent(hide), depend(derf) :: n3 = shape(derf,0)
    integer :: n4
    !f2py intent(hide), depend(derf) :: n4 = shape(derf,1)
    integer :: n5
    !f2py intent(hide), depend(derf) :: n5 = shape(derf,2)
    integer :: n6
    !f2py intent(hide), depend(derf) :: n6 = shape(derf,3)
    call nabla_gamma(ifun=ifun, norder=norder, derf=derf)
end subroutine f90wrap_nabla_gamma

subroutine f90wrap_real_nabla1(ifun, norder, derf, n0, n1, n2, n3, n4, n5, n6)
    use finite_module, only: real_nabla1
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: ifun
    integer(4), intent(in) :: norder
    real(8), intent(inout), dimension(n3,n4,n5,n6) :: derf
    integer :: n0
    !f2py intent(hide), depend(ifun) :: n0 = shape(ifun,0)
    integer :: n1
    !f2py intent(hide), depend(ifun) :: n1 = shape(ifun,1)
    integer :: n2
    !f2py intent(hide), depend(ifun) :: n2 = shape(ifun,2)
    integer :: n3
    !f2py intent(hide), depend(derf) :: n3 = shape(derf,0)
    integer :: n4
    !f2py intent(hide), depend(derf) :: n4 = shape(derf,1)
    integer :: n5
    !f2py intent(hide), depend(derf) :: n5 = shape(derf,2)
    integer :: n6
    !f2py intent(hide), depend(derf) :: n6 = shape(derf,3)
    call real_nabla1(ifun=ifun, norder=norder, derf=derf)
end subroutine f90wrap_real_nabla1

subroutine f90wrap_finite_module__array__Lapl(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use grid_module
    use finite_module, only: finite_module_lapl => lapl
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    if (allocated(finite_module_Lapl)) then
        dshape(1:2) = shape(finite_module_Lapl)
        dloc = loc(finite_module_Lapl)
    else
        dloc = 0
    end if
end subroutine f90wrap_finite_module__array__Lapl

subroutine f90wrap_finite_module__array__Grad(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use grid_module
    use finite_module, only: finite_module_grad => grad
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    if (allocated(finite_module_Grad)) then
        dshape(1:2) = shape(finite_module_Grad)
        dloc = loc(finite_module_Grad)
    else
        dloc = 0
    end if
end subroutine f90wrap_finite_module__array__Grad

subroutine f90wrap_finite_module__array__tBmat(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use grid_module
    use finite_module, only: finite_module_tbmat => tbmat
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    dshape(1:2) = shape(finite_module_tBmat)
    dloc = loc(finite_module_tBmat)
end subroutine f90wrap_finite_module__array__tBmat

subroutine f90wrap_finite_module__array__lap_add(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use grid_module
    use finite_module, only: finite_module_lap_add => lap_add
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    dshape(1:1) = shape(finite_module_lap_add)
    dloc = loc(finite_module_lap_add)
end subroutine f90wrap_finite_module__array__lap_add

subroutine f90wrap_finite_module__array__cell_mu(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use grid_module
    use finite_module, only: finite_module_cell_mu => cell_mu
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 5
    dshape(1:2) = shape(finite_module_cell_mu)
    dloc = loc(finite_module_cell_mu)
end subroutine f90wrap_finite_module__array__cell_mu

subroutine f90wrap_finite_module__array__cell_factor(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use grid_module
    use finite_module, only: finite_module_cell_factor => cell_factor
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    dshape(1:1) = shape(finite_module_cell_factor)
    dloc = loc(finite_module_cell_factor)
end subroutine f90wrap_finite_module__array__cell_factor

subroutine f90wrap_finite_module__array__wrap_box(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use grid_module
    use finite_module, only: finite_module_wrap_box => wrap_box
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 15
    if (allocated(finite_module_wrap_box)) then
        dshape(1:3) = shape(finite_module_wrap_box)
        dloc = loc(finite_module_wrap_box)
    else
        dloc = 0
    end if
end subroutine f90wrap_finite_module__array__wrap_box

subroutine f90wrap_finite_module__array__fun_global(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use grid_module
    use finite_module, only: finite_module_fun_global => fun_global
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 15
    if (allocated(finite_module_fun_global)) then
        dshape(1:3) = shape(finite_module_fun_global)
        dloc = loc(finite_module_fun_global)
    else
        dloc = 0
    end if
end subroutine f90wrap_finite_module__array__fun_global

subroutine f90wrap_finite_module__array__fun_1d(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use grid_module
    use finite_module, only: finite_module_fun_1d => fun_1d
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 15
    if (allocated(finite_module_fun_1d)) then
        dshape(1:1) = shape(finite_module_fun_1d)
        dloc = loc(finite_module_fun_1d)
    else
        dloc = 0
    end if
end subroutine f90wrap_finite_module__array__fun_1d

subroutine f90wrap_finite_module__array__wrap_box1d(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use grid_module
    use finite_module, only: finite_module_wrap_box1d => wrap_box1d
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 15
    if (allocated(finite_module_wrap_box1d)) then
        dshape(1:1) = shape(finite_module_wrap_box1d)
        dloc = loc(finite_module_wrap_box1d)
    else
        dloc = 0
    end if
end subroutine f90wrap_finite_module__array__wrap_box1d

subroutine f90wrap_finite_module__array__wrap_box_real(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use grid_module
    use finite_module, only: finite_module_wrap_box_real => wrap_box_real
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    if (allocated(finite_module_wrap_box_real)) then
        dshape(1:3) = shape(finite_module_wrap_box_real)
        dloc = loc(finite_module_wrap_box_real)
    else
        dloc = 0
    end if
end subroutine f90wrap_finite_module__array__wrap_box_real

subroutine f90wrap_finite_module__array__fun_global_real(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use grid_module
    use finite_module, only: finite_module_fun_global_real => fun_global_real
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    if (allocated(finite_module_fun_global_real)) then
        dshape(1:3) = shape(finite_module_fun_global_real)
        dloc = loc(finite_module_fun_global_real)
    else
        dloc = 0
    end if
end subroutine f90wrap_finite_module__array__fun_global_real

subroutine f90wrap_finite_module__array__fun_1d_real(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use grid_module
    use finite_module, only: finite_module_fun_1d_real => fun_1d_real
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(finite_module_fun_1d_real)) then
        dshape(1:1) = shape(finite_module_fun_1d_real)
        dloc = loc(finite_module_fun_1d_real)
    else
        dloc = 0
    end if
end subroutine f90wrap_finite_module__array__fun_1d_real

subroutine f90wrap_finite_module__array__wrap_box1d_real(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use grid_module
    use finite_module, only: finite_module_wrap_box1d_real => wrap_box1d_real
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(finite_module_wrap_box1d_real)) then
        dshape(1:1) = shape(finite_module_wrap_box1d_real)
        dloc = loc(finite_module_wrap_box1d_real)
    else
        dloc = 0
    end if
end subroutine f90wrap_finite_module__array__wrap_box1d_real

! End of module finite_module defined in file Finite_module.fpp


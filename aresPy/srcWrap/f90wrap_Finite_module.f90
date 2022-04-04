! Module finite_module defined in file Finite_module.f90

subroutine f90wrap_destroy_finite
    use finite_module, only: destroy_finite
    implicit none
    
    call destroy_finite()
end subroutine f90wrap_destroy_finite

subroutine f90wrap_init_finite(h)
    use finite_module, only: init_finite
    implicit none
    
    real(8), dimension(3) :: h
    call init_finite(h=h)
end subroutine f90wrap_init_finite

subroutine f90wrap_trans_mat_full(mat, factor, int_miu, lad_gap, err)
    use finite_module, only: trans_mat_full
    implicit none
    
    real(8), intent(in), dimension(3,3) :: mat
    real(8), dimension(6), intent(inout) :: factor
    integer(8), dimension(3,3), intent(inout) :: int_miu
    real(8), dimension(3), intent(inout) :: lad_gap
    real(8), dimension(6), intent(inout) :: err
    call trans_mat_full(mat=mat, factor=factor, int_miu=int_miu, lad_gap=lad_gap, err=err)
end subroutine f90wrap_trans_mat_full

subroutine f90wrap_cmplx_keop(uk, ik, ts, n0, n1, n2, n3, n4, n5)
    use finite_module, only: cmplx_keop
    implicit none
    
    complex(8), intent(in), dimension(n0,n1,n2) :: uk
    integer(8), intent(in) :: ik
    complex(8), intent(inout), dimension(n3,n4,n5) :: ts
    integer :: n0
    !f2py intent(hide), depend(uk) :: n0 = shape(uk,0)
    integer :: n1
    !f2py intent(hide), depend(uk) :: n1 = shape(uk,1)
    integer :: n2
    !f2py intent(hide), depend(uk) :: n2 = shape(uk,2)
    integer :: n3
    !f2py intent(hide), depend(ts) :: n3 = shape(ts,0)
    integer :: n4
    !f2py intent(hide), depend(ts) :: n4 = shape(ts,1)
    integer :: n5
    !f2py intent(hide), depend(ts) :: n5 = shape(ts,2)
    call cmplx_keop(Uk=uk, Ik=ik, Ts=ts)
end subroutine f90wrap_cmplx_keop

subroutine f90wrap_real_comm(ifun, n0)
    use finite_module, only: real_comm
    implicit none
    
    real(8), intent(in), dimension(n0) :: ifun
    integer :: n0
    !f2py intent(hide), depend(ifun) :: n0 = shape(ifun,0)
    call real_comm(ifun=ifun)
end subroutine f90wrap_real_comm

subroutine f90wrap_real_comm_clean
    use finite_module, only: real_comm_clean
    implicit none
    
    call real_comm_clean()
end subroutine f90wrap_real_comm_clean

subroutine f90wrap_real_nabla1_3d(ifun, norder, derf, n0, n1, n2, n3, n4)
    use finite_module, only: real_nabla1_3d
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: ifun
    integer(8), intent(in) :: norder
    real(8), intent(inout), dimension(n3,n4) :: derf
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
    call real_nabla1_3d(ifun=ifun, norder=norder, derf=derf)
end subroutine f90wrap_real_nabla1_3d

subroutine f90wrap_real_nabla2_3d(func, norder, ofun, n0, n1, n2, n3, n4, n5)
    use finite_module, only: real_nabla2_3d
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: func
    integer(8), intent(in) :: norder
    real(8), intent(inout), dimension(n3,n4,n5) :: ofun
    integer :: n0
    !f2py intent(hide), depend(func) :: n0 = shape(func,0)
    integer :: n1
    !f2py intent(hide), depend(func) :: n1 = shape(func,1)
    integer :: n2
    !f2py intent(hide), depend(func) :: n2 = shape(func,2)
    integer :: n3
    !f2py intent(hide), depend(ofun) :: n3 = shape(ofun,0)
    integer :: n4
    !f2py intent(hide), depend(ofun) :: n4 = shape(ofun,1)
    integer :: n5
    !f2py intent(hide), depend(ofun) :: n5 = shape(ofun,2)
    call real_nabla2_3d(func=func, norder=norder, ofun=ofun)
end subroutine f90wrap_real_nabla2_3d

subroutine f90wrap_real_nabla1_1d(ifun, derf, mgfun, n0, n1, n2)
    use finite_module, only: real_nabla1_1d
    implicit none
    
    real(8), intent(in), dimension(n0) :: ifun
    real(8), intent(inout), dimension(n1,3) :: derf
    real(8), optional, dimension(n2) :: mgfun
    integer :: n0
    !f2py intent(hide), depend(ifun) :: n0 = shape(ifun,0)
    integer :: n1
    !f2py intent(hide), depend(derf) :: n1 = shape(derf,0)
    integer :: n2
    !f2py intent(hide), depend(mgfun) :: n2 = shape(mgfun,0)
    call real_nabla1_1d(ifun=ifun, derf=derf, mgfun=mgfun)
end subroutine f90wrap_real_nabla1_1d

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
    dtype = 7
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

subroutine f90wrap_finite_module__array__wrap_real(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use grid_module
    use finite_module, only: finite_module_wrap_real => wrap_real
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(finite_module_wrap_real)) then
        dshape(1:1) = shape(finite_module_wrap_real)
        dloc = loc(finite_module_wrap_real)
    else
        dloc = 0
    end if
end subroutine f90wrap_finite_module__array__wrap_real

! End of module finite_module defined in file Finite_module.f90


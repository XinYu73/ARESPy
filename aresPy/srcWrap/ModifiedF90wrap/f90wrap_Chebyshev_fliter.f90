! Module chebyshev_module defined in file Chebyshev_fliter.fpp

subroutine f90wrap_buildsubspace(nps, nev, veff, eig, n0, n1)
    use chebyshev_module, only: buildsubspace
    use grid_module, only: eigen_type
    implicit none
    
    type eigen_type_ptr_type
        type(eigen_type), pointer :: p => NULL()
    end type eigen_type_ptr_type
    integer(4), intent(in) :: nps
    integer(4), intent(in) :: nev
    real(8), intent(in), dimension(n0,n1) :: veff
    type(eigen_type_ptr_type) :: eig_ptr
    integer, intent(in), dimension(2) :: eig
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    eig_ptr = transfer(eig, eig_ptr)
    call buildsubspace(nps=nps, nev=nev, veff=veff, eig=eig_ptr%p)
end subroutine f90wrap_buildsubspace

subroutine f90wrap_real_pseudosubspace(nps, nev, initx, n0, n1)
    use chebyshev_module, only: real_pseudosubspace
    implicit none
    
    integer(4), intent(in) :: nps
    integer(4), intent(in) :: nev
    real(8), intent(inout), dimension(n0,n1) :: initx
    integer :: n0
    !f2py intent(hide), depend(initx) :: n0 = shape(initx,0)
    integer :: n1
    !f2py intent(hide), depend(initx) :: n1 = shape(initx,1)
    call real_pseudosubspace(nps=nps, nev=nev, initX=initx)
end subroutine f90wrap_real_pseudosubspace

subroutine f90wrap_real_first_rrstep(nps, nev, veff, x, d, n0, n1, n2, n3)
    use chebyshev_module, only: real_first_rrstep
    implicit none
    
    integer(4), intent(in) :: nps
    integer(4), intent(in) :: nev
    real(8), intent(in), dimension(n0) :: veff
    real(8), intent(inout), dimension(n1,n2) :: x
    real(8), intent(inout), dimension(n3) :: d
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(x) :: n1 = shape(x,0)
    integer :: n2
    !f2py intent(hide), depend(x) :: n2 = shape(x,1)
    integer :: n3
    !f2py intent(hide), depend(d) :: n3 = shape(d,0)
    call real_first_rrstep(nps=nps, nev=nev, veff=veff, X=x, D=d)
end subroutine f90wrap_real_first_rrstep

subroutine f90wrap_init_uplow_real(nps, k, veff, v, a, b, al, n0, n1)
    use chebyshev_module, only: init_uplow_real
    implicit none
    
    integer(4), intent(in) :: nps
    integer(4), intent(in) :: k
    real(8), intent(in), dimension(n0) :: veff
    real(8), intent(inout), dimension(n1) :: v
    real(8), intent(out) :: a
    real(8), intent(out) :: b
    real(8), intent(out) :: al
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(v) :: n1 = shape(v,0)
    call init_uplow_real(nps=nps, k=k, veff=veff, v=v, a=a, b=b, al=al)
end subroutine f90wrap_init_uplow_real

subroutine f90wrap_real_first_filter(nps, nst, veff, x, eval, n0, n1, n2, n3)
    use chebyshev_module, only: real_first_filter
    implicit none
    
    integer(4), intent(in) :: nps
    integer(4), intent(in) :: nst
    real(8), intent(in), dimension(n0) :: veff
    real(8), intent(inout), dimension(n1,n2) :: x
    real(8), intent(inout), dimension(n3) :: eval
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(x) :: n1 = shape(x,0)
    integer :: n2
    !f2py intent(hide), depend(x) :: n2 = shape(x,1)
    integer :: n3
    !f2py intent(hide), depend(eval) :: n3 = shape(eval,0)
    call real_first_filter(nps=nps, nst=nst, veff=veff, X=x, eval=eval)
end subroutine f90wrap_real_first_filter

subroutine f90wrap_rayleigh_quotient_real(nps, nst, veff, x, xhx, n0, n1, n2, n3, n4)
    use chebyshev_module, only: rayleigh_quotient_real
    implicit none
    
    integer(4), intent(in) :: nps
    integer(4), intent(in) :: nst
    real(8), intent(in), dimension(n0) :: veff
    real(8), dimension(n1,n2) :: x
    real(8), dimension(n3,n4) :: xhx
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(x) :: n1 = shape(x,0)
    integer :: n2
    !f2py intent(hide), depend(x) :: n2 = shape(x,1)
    integer :: n3
    !f2py intent(hide), depend(xhx) :: n3 = shape(xhx,0)
    integer :: n4
    !f2py intent(hide), depend(xhx) :: n4 = shape(xhx,1)
    call rayleigh_quotient_real(nps=nps, nst=nst, veff=veff, x=x, xhx=xhx)
end subroutine f90wrap_rayleigh_quotient_real

subroutine f90wrap_cal_hx_real(nps, nst, veff, v, hv, n0, n1, n2, n3, n4)
    use chebyshev_module, only: cal_hx_real
    implicit none
    
    integer(4), intent(in) :: nps
    integer(4), intent(in) :: nst
    real(8), intent(in), dimension(n0) :: veff
    real(8), intent(in), dimension(n1,n2) :: v
    real(8), intent(inout), dimension(n3,n4) :: hv
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(v) :: n1 = shape(v,0)
    integer :: n2
    !f2py intent(hide), depend(v) :: n2 = shape(v,1)
    integer :: n3
    !f2py intent(hide), depend(hv) :: n3 = shape(hv,0)
    integer :: n4
    !f2py intent(hide), depend(hv) :: n4 = shape(hv,1)
    call cal_hx_real(nps=nps, nst=nst, veff=veff, V=v, HV=hv)
end subroutine f90wrap_cal_hx_real

subroutine f90wrap_estupb_real(nps, k, veff, vec, b, n0, n1)
    use chebyshev_module, only: estupb_real
    implicit none
    
    integer(4), intent(in) :: nps
    integer(4), intent(in) :: k
    real(8), intent(in), dimension(n0) :: veff
    real(8), intent(in), dimension(n1) :: vec
    real(8), intent(out) :: b
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(vec) :: n1 = shape(vec,0)
    call estupb_real(nps=nps, k=k, veff=veff, vec=vec, b=b)
end subroutine f90wrap_estupb_real

subroutine f90wrap_chebyshev_filter_real(nps, nst, veff, x, m, a, b, n0, n1, n2)
    use chebyshev_module, only: chebyshev_filter_real
    implicit none
    
    integer(4), intent(in) :: nps
    integer(4), intent(in) :: nst
    real(8), intent(in), dimension(n0) :: veff
    real(8), intent(inout), dimension(n1,n2) :: x
    integer(4), intent(in) :: m
    real(8), intent(in) :: a
    real(8), intent(in) :: b
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(x) :: n1 = shape(x,0)
    integer :: n2
    !f2py intent(hide), depend(x) :: n2 = shape(x,1)
    call chebyshev_filter_real(nps=nps, nst=nst, veff=veff, X=x, m=m, a=a, b=b)
end subroutine f90wrap_chebyshev_filter_real

subroutine f90wrap_chebyshev_filter_scaled_real(nps, nst, veff, x, m, a, b, al, n0, n1, n2)
    use chebyshev_module, only: chebyshev_filter_scaled_real
    implicit none
    
    integer(4), intent(in) :: nps
    integer(4), intent(in) :: nst
    real(8), intent(in), dimension(n0) :: veff
    real(8), intent(inout), dimension(n1,n2) :: x
    integer(4), intent(in) :: m
    real(8), intent(in) :: a
    real(8), intent(in) :: b
    real(8), intent(in) :: al
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(x) :: n1 = shape(x,0)
    integer :: n2
    !f2py intent(hide), depend(x) :: n2 = shape(x,1)
    call chebyshev_filter_scaled_real(nps=nps, nst=nst, veff=veff, X=x, m=m, a=a, b=b, al=al)
end subroutine f90wrap_chebyshev_filter_scaled_real

subroutine f90wrap_grayleigh_ritz_real(nps, nev, veff, x, d, n0, n1, n2, n3)
    use chebyshev_module, only: grayleigh_ritz_real
    implicit none
    
    integer(4), intent(in) :: nps
    integer(4), intent(in) :: nev
    real(8), intent(in), dimension(n0) :: veff
    real(8), intent(inout), dimension(n1,n2) :: x
    real(8), intent(inout), dimension(n3) :: d
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(x) :: n1 = shape(x,0)
    integer :: n2
    !f2py intent(hide), depend(x) :: n2 = shape(x,1)
    integer :: n3
    !f2py intent(hide), depend(d) :: n3 = shape(d,0)
    call grayleigh_ritz_real(nps=nps, nev=nev, veff=veff, X=x, D=d)
end subroutine f90wrap_grayleigh_ritz_real

subroutine f90wrap_rayleigh_ritz_real(nps, sn, veff, x, d, n0, n1, n2, n3)
    use chebyshev_module, only: rayleigh_ritz_real
    implicit none
    
    integer(4), intent(in) :: nps
    integer(4), intent(in) :: sn
    real(8), intent(in), dimension(n0) :: veff
    real(8), intent(inout), dimension(n1,n2) :: x
    real(8), intent(inout), dimension(n3) :: d
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(x) :: n1 = shape(x,0)
    integer :: n2
    !f2py intent(hide), depend(x) :: n2 = shape(x,1)
    integer :: n3
    !f2py intent(hide), depend(d) :: n3 = shape(d,0)
    call rayleigh_ritz_real(nps=nps, sn=sn, veff=veff, X=x, D=d)
end subroutine f90wrap_rayleigh_ritz_real

subroutine f90wrap_cheby_filtering_grrr(nps, nev, veff, x, d, n0, n1, n2, n3)
    use chebyshev_module, only: cheby_filtering_grrr
    implicit none
    
    integer(4), intent(in) :: nps
    integer(4), intent(in) :: nev
    real(8), intent(in), dimension(n0) :: veff
    real(8), intent(inout), dimension(n1,n2) :: x
    real(8), intent(inout), dimension(n3) :: d
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(x) :: n1 = shape(x,0)
    integer :: n2
    !f2py intent(hide), depend(x) :: n2 = shape(x,1)
    integer :: n3
    !f2py intent(hide), depend(d) :: n3 = shape(d,0)
    call cheby_filtering_grrr(nps=nps, nev=nev, veff=veff, X=x, D=d)
end subroutine f90wrap_cheby_filtering_grrr

subroutine f90wrap_chebyshev_module__get__LARGED(f90wrap_LARGED)
    use chebyshev_module, only: chebyshev_module_LARGED => LARGED
    implicit none
    real(8), intent(out) :: f90wrap_LARGED
    
    f90wrap_LARGED = chebyshev_module_LARGED
end subroutine f90wrap_chebyshev_module__get__LARGED

subroutine f90wrap_chebyshev_module__array__ad(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters
    use smpi_math_module
    use chebyshev_module, only: chebyshev_module_ad => ad
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(chebyshev_module_ad)) then
        dshape(1:1) = shape(chebyshev_module_ad)
        dloc = loc(chebyshev_module_ad)
    else
        dloc = 0
    end if
end subroutine f90wrap_chebyshev_module__array__ad

! End of module chebyshev_module defined in file Chebyshev_fliter.fpp


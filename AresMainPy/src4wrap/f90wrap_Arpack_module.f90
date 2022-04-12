! Module arpack_module defined in file Arpack_module.fpp

subroutine f90wrap_diagh_arpack(veff, ik, nev, evec, eval, resid_restart, nec, info, maxmvs, tol, n0, n1, n2, n3, n4, &
    n5, n6)
    use arpack_module, only: diagh_arpack
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    integer(4), intent(in) :: ik
    integer(4), intent(in) :: nev
    complex(8), intent(inout), dimension(n3,n4) :: evec
    real(8), intent(inout), dimension(n5) :: eval
    complex(8), intent(inout), dimension(n6) :: resid_restart
    integer(4), intent(inout) :: nec
    integer(4), intent(inout) :: info
    integer(4), intent(inout) :: maxmvs
    real(8), intent(in) :: tol
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(evec) :: n3 = shape(evec,0)
    integer :: n4
    !f2py intent(hide), depend(evec) :: n4 = shape(evec,1)
    integer :: n5
    !f2py intent(hide), depend(eval) :: n5 = shape(eval,0)
    integer :: n6
    !f2py intent(hide), depend(resid_restart) :: n6 = shape(resid_restart,0)
    call diagh_arpack(veff=veff, Ik=ik, nev=nev, evec=evec, eval=eval, resid_restart=resid_restart, nec=nec, info=info, &
        maxmvs=maxmvs, TOL=tol)
end subroutine f90wrap_diagh_arpack

subroutine f90wrap_real_diagh_arpack(veff, nev, evec, eval, resid_restart, nec, info, maxmvs, tol, n0, n1, n2, n3, n4, &
    n5, n6)
    use arpack_module, only: real_diagh_arpack
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    integer(4), intent(in) :: nev
    real(8), intent(inout), dimension(n3,n4) :: evec
    real(8), intent(inout), dimension(n5) :: eval
    real(8), intent(inout), dimension(n6) :: resid_restart
    integer(4), intent(inout) :: nec
    integer(4), intent(inout) :: info
    integer(4), intent(inout) :: maxmvs
    real(8), intent(in) :: tol
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(evec) :: n3 = shape(evec,0)
    integer :: n4
    !f2py intent(hide), depend(evec) :: n4 = shape(evec,1)
    integer :: n5
    !f2py intent(hide), depend(eval) :: n5 = shape(eval,0)
    integer :: n6
    !f2py intent(hide), depend(resid_restart) :: n6 = shape(resid_restart,0)
    call real_diagh_arpack(veff=veff, nev=nev, evec=evec, eval=eval, resid_restart=resid_restart, nec=nec, info=info, &
        maxmvs=maxmvs, TOL=tol)
end subroutine f90wrap_real_diagh_arpack

subroutine f90wrap_real_diagm_arpk(mat, nev, evec, eval, resid_restart, nec, info, maxmvs, tol, n0, n1, n2, n3, n4, n5)
    use arpack_module, only: real_diagm_arpk
    implicit none
    
    real(8), intent(in), dimension(n0,n1) :: mat
    integer(4), intent(in) :: nev
    real(8), intent(inout), dimension(n2,n3) :: evec
    real(8), intent(inout), dimension(n4) :: eval
    real(8), intent(inout), dimension(n5) :: resid_restart
    integer(4), intent(inout) :: nec
    integer(4), intent(inout) :: info
    integer(4), intent(inout) :: maxmvs
    real(8), intent(in) :: tol
    integer :: n0
    !f2py intent(hide), depend(mat) :: n0 = shape(mat,0)
    integer :: n1
    !f2py intent(hide), depend(mat) :: n1 = shape(mat,1)
    integer :: n2
    !f2py intent(hide), depend(evec) :: n2 = shape(evec,0)
    integer :: n3
    !f2py intent(hide), depend(evec) :: n3 = shape(evec,1)
    integer :: n4
    !f2py intent(hide), depend(eval) :: n4 = shape(eval,0)
    integer :: n5
    !f2py intent(hide), depend(resid_restart) :: n5 = shape(resid_restart,0)
    call real_diagm_arpk(mat=mat, nev=nev, evec=evec, eval=eval, resid_restart=resid_restart, nec=nec, info=info, &
        maxmvs=maxmvs, TOL=tol)
end subroutine f90wrap_real_diagm_arpk

subroutine f90wrap_rdiagm_arpk(mat, dimen, nev, evec, eval, n0, n1, n2, n3, n4)
    use arpack_module, only: rdiagm_arpk
    implicit none
    
    real(8), intent(in), dimension(n0,n1) :: mat
    integer(4), intent(in) :: dimen
    integer(4), intent(in) :: nev
    real(8), intent(inout), dimension(n2,n3) :: evec
    real(8), intent(inout), dimension(n4) :: eval
    integer :: n0
    !f2py intent(hide), depend(mat) :: n0 = shape(mat,0)
    integer :: n1
    !f2py intent(hide), depend(mat) :: n1 = shape(mat,1)
    integer :: n2
    !f2py intent(hide), depend(evec) :: n2 = shape(evec,0)
    integer :: n3
    !f2py intent(hide), depend(evec) :: n3 = shape(evec,1)
    integer :: n4
    !f2py intent(hide), depend(eval) :: n4 = shape(eval,0)
    call rdiagm_arpk(mat=mat, dimen=dimen, nev=nev, evec=evec, eval=eval)
end subroutine f90wrap_rdiagm_arpk

subroutine f90wrap_iso_diagh_arpack(veff, nev, evec, eval, resid_restart, nec, info, maxmvs, tol, n0, n1, n2, n3, n4)
    use arpack_module, only: iso_diagh_arpack
    implicit none
    
    real(8), intent(in), dimension(n0) :: veff
    integer(4), intent(in) :: nev
    real(8), intent(inout), dimension(n1,n2) :: evec
    real(8), intent(inout), dimension(n3) :: eval
    real(8), intent(inout), dimension(n4) :: resid_restart
    integer(4), intent(inout) :: nec
    integer(4), intent(inout) :: info
    integer(4), intent(inout) :: maxmvs
    real(8), intent(in) :: tol
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(evec) :: n1 = shape(evec,0)
    integer :: n2
    !f2py intent(hide), depend(evec) :: n2 = shape(evec,1)
    integer :: n3
    !f2py intent(hide), depend(eval) :: n3 = shape(eval,0)
    integer :: n4
    !f2py intent(hide), depend(resid_restart) :: n4 = shape(resid_restart,0)
    call iso_diagh_arpack(veff=veff, nev=nev, evec=evec, eval=eval, resid_restart=resid_restart, nec=nec, info=info, &
        maxmvs=maxmvs, TOL=tol)
end subroutine f90wrap_iso_diagh_arpack

subroutine f90wrap_diagh_arpack_band(veff, ik, nev, evec, eval, resid_restart, nec, info, maxmvs, tol, n0, n1, n2, n3, &
    n4, n5, n6)
    use arpack_module, only: diagh_arpack_band
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: veff
    integer(4), intent(in) :: ik
    integer(4), intent(in) :: nev
    complex(8), intent(inout), dimension(n3,n4) :: evec
    real(8), intent(inout), dimension(n5) :: eval
    complex(8), intent(inout), dimension(n6) :: resid_restart
    integer(4), intent(inout) :: nec
    integer(4), intent(inout) :: info
    integer(4), intent(inout) :: maxmvs
    real(8), intent(in) :: tol
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(veff) :: n2 = shape(veff,2)
    integer :: n3
    !f2py intent(hide), depend(evec) :: n3 = shape(evec,0)
    integer :: n4
    !f2py intent(hide), depend(evec) :: n4 = shape(evec,1)
    integer :: n5
    !f2py intent(hide), depend(eval) :: n5 = shape(eval,0)
    integer :: n6
    !f2py intent(hide), depend(resid_restart) :: n6 = shape(resid_restart,0)
    call diagh_arpack_band(veff=veff, Ik=ik, nev=nev, evec=evec, eval=eval, resid_restart=resid_restart, nec=nec, info=info, &
        maxmvs=maxmvs, TOL=tol)
end subroutine f90wrap_diagh_arpack_band

subroutine f90wrap_arpack_module__get__maxn(f90wrap_maxn)
    use arpack_module, only: arpack_module_maxn => maxn
    implicit none
    integer(4), intent(out) :: f90wrap_maxn
    
    f90wrap_maxn = arpack_module_maxn
end subroutine f90wrap_arpack_module__get__maxn

subroutine f90wrap_arpack_module__set__maxn(f90wrap_maxn)
    use arpack_module, only: arpack_module_maxn => maxn
    implicit none
    integer(4), intent(in) :: f90wrap_maxn
    
    arpack_module_maxn = f90wrap_maxn
end subroutine f90wrap_arpack_module__set__maxn

subroutine f90wrap_arpack_module__get__maxnev(f90wrap_maxnev)
    use arpack_module, only: arpack_module_maxnev => maxnev
    implicit none
    integer(4), intent(out) :: f90wrap_maxnev
    
    f90wrap_maxnev = arpack_module_maxnev
end subroutine f90wrap_arpack_module__get__maxnev

subroutine f90wrap_arpack_module__set__maxnev(f90wrap_maxnev)
    use arpack_module, only: arpack_module_maxnev => maxnev
    implicit none
    integer(4), intent(in) :: f90wrap_maxnev
    
    arpack_module_maxnev = f90wrap_maxnev
end subroutine f90wrap_arpack_module__set__maxnev

subroutine f90wrap_arpack_module__get__maxncv(f90wrap_maxncv)
    use arpack_module, only: arpack_module_maxncv => maxncv
    implicit none
    integer(4), intent(out) :: f90wrap_maxncv
    
    f90wrap_maxncv = arpack_module_maxncv
end subroutine f90wrap_arpack_module__get__maxncv

subroutine f90wrap_arpack_module__set__maxncv(f90wrap_maxncv)
    use arpack_module, only: arpack_module_maxncv => maxncv
    implicit none
    integer(4), intent(in) :: f90wrap_maxncv
    
    arpack_module_maxncv = f90wrap_maxncv
end subroutine f90wrap_arpack_module__set__maxncv

! End of module arpack_module defined in file Arpack_module.fpp


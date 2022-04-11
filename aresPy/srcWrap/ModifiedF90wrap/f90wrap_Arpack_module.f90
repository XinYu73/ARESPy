! Module arpack_module defined in file Arpack_module.fpp

subroutine f90wrap_real_diagh_arpack(n, veff, nev, evec, eval, resid_restart, nec, info, maxmvs, tol, n0, n1, n2, n3, &
    n4)
    use arpack_module, only: real_diagh_arpack
    implicit none
    
    integer(4), intent(in) :: n
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
    call real_diagh_arpack(n=n, veff=veff, nev=nev, evec=evec, eval=eval, resid_restart=resid_restart, nec=nec, info=info, &
        maxmvs=maxmvs, TOL=tol)
end subroutine f90wrap_real_diagh_arpack

subroutine f90wrap_arpack_module__get__maxn(f90wrap_maxn)
    use arpack_module, only: arpack_module_maxn => maxn
    implicit none
    integer(4), intent(out) :: f90wrap_maxn
    
    f90wrap_maxn = arpack_module_maxn
end subroutine f90wrap_arpack_module__get__maxn

subroutine f90wrap_arpack_module__get__maxnev(f90wrap_maxnev)
    use arpack_module, only: arpack_module_maxnev => maxnev
    implicit none
    integer(4), intent(out) :: f90wrap_maxnev
    
    f90wrap_maxnev = arpack_module_maxnev
end subroutine f90wrap_arpack_module__get__maxnev

subroutine f90wrap_arpack_module__get__maxncv(f90wrap_maxncv)
    use arpack_module, only: arpack_module_maxncv => maxncv
    implicit none
    integer(4), intent(out) :: f90wrap_maxncv
    
    f90wrap_maxncv = arpack_module_maxncv
end subroutine f90wrap_arpack_module__get__maxncv

subroutine f90wrap_arpack_module__get__logfil(f90wrap_logfil)
    use arpack_module, only: arpack_module_logfil => logfil
    implicit none
    integer, intent(out) :: f90wrap_logfil
    
    f90wrap_logfil = arpack_module_logfil
end subroutine f90wrap_arpack_module__get__logfil

subroutine f90wrap_arpack_module__set__logfil(f90wrap_logfil)
    use arpack_module, only: arpack_module_logfil => logfil
    implicit none
    integer, intent(in) :: f90wrap_logfil
    
    arpack_module_logfil = f90wrap_logfil
end subroutine f90wrap_arpack_module__set__logfil

subroutine f90wrap_arpack_module__get__ndigit(f90wrap_ndigit)
    use arpack_module, only: arpack_module_ndigit => ndigit
    implicit none
    integer, intent(out) :: f90wrap_ndigit
    
    f90wrap_ndigit = arpack_module_ndigit
end subroutine f90wrap_arpack_module__get__ndigit

subroutine f90wrap_arpack_module__set__ndigit(f90wrap_ndigit)
    use arpack_module, only: arpack_module_ndigit => ndigit
    implicit none
    integer, intent(in) :: f90wrap_ndigit
    
    arpack_module_ndigit = f90wrap_ndigit
end subroutine f90wrap_arpack_module__set__ndigit

subroutine f90wrap_arpack_module__get__mgetv0(f90wrap_mgetv0)
    use arpack_module, only: arpack_module_mgetv0 => mgetv0
    implicit none
    integer, intent(out) :: f90wrap_mgetv0
    
    f90wrap_mgetv0 = arpack_module_mgetv0
end subroutine f90wrap_arpack_module__get__mgetv0

subroutine f90wrap_arpack_module__set__mgetv0(f90wrap_mgetv0)
    use arpack_module, only: arpack_module_mgetv0 => mgetv0
    implicit none
    integer, intent(in) :: f90wrap_mgetv0
    
    arpack_module_mgetv0 = f90wrap_mgetv0
end subroutine f90wrap_arpack_module__set__mgetv0

subroutine f90wrap_arpack_module__get__msaupd(f90wrap_msaupd)
    use arpack_module, only: arpack_module_msaupd => msaupd
    implicit none
    integer, intent(out) :: f90wrap_msaupd
    
    f90wrap_msaupd = arpack_module_msaupd
end subroutine f90wrap_arpack_module__get__msaupd

subroutine f90wrap_arpack_module__set__msaupd(f90wrap_msaupd)
    use arpack_module, only: arpack_module_msaupd => msaupd
    implicit none
    integer, intent(in) :: f90wrap_msaupd
    
    arpack_module_msaupd = f90wrap_msaupd
end subroutine f90wrap_arpack_module__set__msaupd

subroutine f90wrap_arpack_module__get__msaup2(f90wrap_msaup2)
    use arpack_module, only: arpack_module_msaup2 => msaup2
    implicit none
    integer, intent(out) :: f90wrap_msaup2
    
    f90wrap_msaup2 = arpack_module_msaup2
end subroutine f90wrap_arpack_module__get__msaup2

subroutine f90wrap_arpack_module__set__msaup2(f90wrap_msaup2)
    use arpack_module, only: arpack_module_msaup2 => msaup2
    implicit none
    integer, intent(in) :: f90wrap_msaup2
    
    arpack_module_msaup2 = f90wrap_msaup2
end subroutine f90wrap_arpack_module__set__msaup2

subroutine f90wrap_arpack_module__get__msaitr(f90wrap_msaitr)
    use arpack_module, only: arpack_module_msaitr => msaitr
    implicit none
    integer, intent(out) :: f90wrap_msaitr
    
    f90wrap_msaitr = arpack_module_msaitr
end subroutine f90wrap_arpack_module__get__msaitr

subroutine f90wrap_arpack_module__set__msaitr(f90wrap_msaitr)
    use arpack_module, only: arpack_module_msaitr => msaitr
    implicit none
    integer, intent(in) :: f90wrap_msaitr
    
    arpack_module_msaitr = f90wrap_msaitr
end subroutine f90wrap_arpack_module__set__msaitr

subroutine f90wrap_arpack_module__get__mseigt(f90wrap_mseigt)
    use arpack_module, only: arpack_module_mseigt => mseigt
    implicit none
    integer, intent(out) :: f90wrap_mseigt
    
    f90wrap_mseigt = arpack_module_mseigt
end subroutine f90wrap_arpack_module__get__mseigt

subroutine f90wrap_arpack_module__set__mseigt(f90wrap_mseigt)
    use arpack_module, only: arpack_module_mseigt => mseigt
    implicit none
    integer, intent(in) :: f90wrap_mseigt
    
    arpack_module_mseigt = f90wrap_mseigt
end subroutine f90wrap_arpack_module__set__mseigt

subroutine f90wrap_arpack_module__get__msapps(f90wrap_msapps)
    use arpack_module, only: arpack_module_msapps => msapps
    implicit none
    integer, intent(out) :: f90wrap_msapps
    
    f90wrap_msapps = arpack_module_msapps
end subroutine f90wrap_arpack_module__get__msapps

subroutine f90wrap_arpack_module__set__msapps(f90wrap_msapps)
    use arpack_module, only: arpack_module_msapps => msapps
    implicit none
    integer, intent(in) :: f90wrap_msapps
    
    arpack_module_msapps = f90wrap_msapps
end subroutine f90wrap_arpack_module__set__msapps

subroutine f90wrap_arpack_module__get__msgets(f90wrap_msgets)
    use arpack_module, only: arpack_module_msgets => msgets
    implicit none
    integer, intent(out) :: f90wrap_msgets
    
    f90wrap_msgets = arpack_module_msgets
end subroutine f90wrap_arpack_module__get__msgets

subroutine f90wrap_arpack_module__set__msgets(f90wrap_msgets)
    use arpack_module, only: arpack_module_msgets => msgets
    implicit none
    integer, intent(in) :: f90wrap_msgets
    
    arpack_module_msgets = f90wrap_msgets
end subroutine f90wrap_arpack_module__set__msgets

subroutine f90wrap_arpack_module__get__mseupd(f90wrap_mseupd)
    use arpack_module, only: arpack_module_mseupd => mseupd
    implicit none
    integer, intent(out) :: f90wrap_mseupd
    
    f90wrap_mseupd = arpack_module_mseupd
end subroutine f90wrap_arpack_module__get__mseupd

subroutine f90wrap_arpack_module__set__mseupd(f90wrap_mseupd)
    use arpack_module, only: arpack_module_mseupd => mseupd
    implicit none
    integer, intent(in) :: f90wrap_mseupd
    
    arpack_module_mseupd = f90wrap_mseupd
end subroutine f90wrap_arpack_module__set__mseupd

subroutine f90wrap_arpack_module__get__mnaupd(f90wrap_mnaupd)
    use arpack_module, only: arpack_module_mnaupd => mnaupd
    implicit none
    integer, intent(out) :: f90wrap_mnaupd
    
    f90wrap_mnaupd = arpack_module_mnaupd
end subroutine f90wrap_arpack_module__get__mnaupd

subroutine f90wrap_arpack_module__set__mnaupd(f90wrap_mnaupd)
    use arpack_module, only: arpack_module_mnaupd => mnaupd
    implicit none
    integer, intent(in) :: f90wrap_mnaupd
    
    arpack_module_mnaupd = f90wrap_mnaupd
end subroutine f90wrap_arpack_module__set__mnaupd

subroutine f90wrap_arpack_module__get__mnaup2(f90wrap_mnaup2)
    use arpack_module, only: arpack_module_mnaup2 => mnaup2
    implicit none
    integer, intent(out) :: f90wrap_mnaup2
    
    f90wrap_mnaup2 = arpack_module_mnaup2
end subroutine f90wrap_arpack_module__get__mnaup2

subroutine f90wrap_arpack_module__set__mnaup2(f90wrap_mnaup2)
    use arpack_module, only: arpack_module_mnaup2 => mnaup2
    implicit none
    integer, intent(in) :: f90wrap_mnaup2
    
    arpack_module_mnaup2 = f90wrap_mnaup2
end subroutine f90wrap_arpack_module__set__mnaup2

subroutine f90wrap_arpack_module__get__mnaitr(f90wrap_mnaitr)
    use arpack_module, only: arpack_module_mnaitr => mnaitr
    implicit none
    integer, intent(out) :: f90wrap_mnaitr
    
    f90wrap_mnaitr = arpack_module_mnaitr
end subroutine f90wrap_arpack_module__get__mnaitr

subroutine f90wrap_arpack_module__set__mnaitr(f90wrap_mnaitr)
    use arpack_module, only: arpack_module_mnaitr => mnaitr
    implicit none
    integer, intent(in) :: f90wrap_mnaitr
    
    arpack_module_mnaitr = f90wrap_mnaitr
end subroutine f90wrap_arpack_module__set__mnaitr

subroutine f90wrap_arpack_module__get__mneigh(f90wrap_mneigh)
    use arpack_module, only: arpack_module_mneigh => mneigh
    implicit none
    integer, intent(out) :: f90wrap_mneigh
    
    f90wrap_mneigh = arpack_module_mneigh
end subroutine f90wrap_arpack_module__get__mneigh

subroutine f90wrap_arpack_module__set__mneigh(f90wrap_mneigh)
    use arpack_module, only: arpack_module_mneigh => mneigh
    implicit none
    integer, intent(in) :: f90wrap_mneigh
    
    arpack_module_mneigh = f90wrap_mneigh
end subroutine f90wrap_arpack_module__set__mneigh

subroutine f90wrap_arpack_module__get__mnapps(f90wrap_mnapps)
    use arpack_module, only: arpack_module_mnapps => mnapps
    implicit none
    integer, intent(out) :: f90wrap_mnapps
    
    f90wrap_mnapps = arpack_module_mnapps
end subroutine f90wrap_arpack_module__get__mnapps

subroutine f90wrap_arpack_module__set__mnapps(f90wrap_mnapps)
    use arpack_module, only: arpack_module_mnapps => mnapps
    implicit none
    integer, intent(in) :: f90wrap_mnapps
    
    arpack_module_mnapps = f90wrap_mnapps
end subroutine f90wrap_arpack_module__set__mnapps

subroutine f90wrap_arpack_module__get__mngets(f90wrap_mngets)
    use arpack_module, only: arpack_module_mngets => mngets
    implicit none
    integer, intent(out) :: f90wrap_mngets
    
    f90wrap_mngets = arpack_module_mngets
end subroutine f90wrap_arpack_module__get__mngets

subroutine f90wrap_arpack_module__set__mngets(f90wrap_mngets)
    use arpack_module, only: arpack_module_mngets => mngets
    implicit none
    integer, intent(in) :: f90wrap_mngets
    
    arpack_module_mngets = f90wrap_mngets
end subroutine f90wrap_arpack_module__set__mngets

subroutine f90wrap_arpack_module__get__mneupd(f90wrap_mneupd)
    use arpack_module, only: arpack_module_mneupd => mneupd
    implicit none
    integer, intent(out) :: f90wrap_mneupd
    
    f90wrap_mneupd = arpack_module_mneupd
end subroutine f90wrap_arpack_module__get__mneupd

subroutine f90wrap_arpack_module__set__mneupd(f90wrap_mneupd)
    use arpack_module, only: arpack_module_mneupd => mneupd
    implicit none
    integer, intent(in) :: f90wrap_mneupd
    
    arpack_module_mneupd = f90wrap_mneupd
end subroutine f90wrap_arpack_module__set__mneupd

subroutine f90wrap_arpack_module__get__mcaupd(f90wrap_mcaupd)
    use arpack_module, only: arpack_module_mcaupd => mcaupd
    implicit none
    integer, intent(out) :: f90wrap_mcaupd
    
    f90wrap_mcaupd = arpack_module_mcaupd
end subroutine f90wrap_arpack_module__get__mcaupd

subroutine f90wrap_arpack_module__set__mcaupd(f90wrap_mcaupd)
    use arpack_module, only: arpack_module_mcaupd => mcaupd
    implicit none
    integer, intent(in) :: f90wrap_mcaupd
    
    arpack_module_mcaupd = f90wrap_mcaupd
end subroutine f90wrap_arpack_module__set__mcaupd

subroutine f90wrap_arpack_module__get__mcaup2(f90wrap_mcaup2)
    use arpack_module, only: arpack_module_mcaup2 => mcaup2
    implicit none
    integer, intent(out) :: f90wrap_mcaup2
    
    f90wrap_mcaup2 = arpack_module_mcaup2
end subroutine f90wrap_arpack_module__get__mcaup2

subroutine f90wrap_arpack_module__set__mcaup2(f90wrap_mcaup2)
    use arpack_module, only: arpack_module_mcaup2 => mcaup2
    implicit none
    integer, intent(in) :: f90wrap_mcaup2
    
    arpack_module_mcaup2 = f90wrap_mcaup2
end subroutine f90wrap_arpack_module__set__mcaup2

subroutine f90wrap_arpack_module__get__mcaitr(f90wrap_mcaitr)
    use arpack_module, only: arpack_module_mcaitr => mcaitr
    implicit none
    integer, intent(out) :: f90wrap_mcaitr
    
    f90wrap_mcaitr = arpack_module_mcaitr
end subroutine f90wrap_arpack_module__get__mcaitr

subroutine f90wrap_arpack_module__set__mcaitr(f90wrap_mcaitr)
    use arpack_module, only: arpack_module_mcaitr => mcaitr
    implicit none
    integer, intent(in) :: f90wrap_mcaitr
    
    arpack_module_mcaitr = f90wrap_mcaitr
end subroutine f90wrap_arpack_module__set__mcaitr

subroutine f90wrap_arpack_module__get__mceigh(f90wrap_mceigh)
    use arpack_module, only: arpack_module_mceigh => mceigh
    implicit none
    integer, intent(out) :: f90wrap_mceigh
    
    f90wrap_mceigh = arpack_module_mceigh
end subroutine f90wrap_arpack_module__get__mceigh

subroutine f90wrap_arpack_module__set__mceigh(f90wrap_mceigh)
    use arpack_module, only: arpack_module_mceigh => mceigh
    implicit none
    integer, intent(in) :: f90wrap_mceigh
    
    arpack_module_mceigh = f90wrap_mceigh
end subroutine f90wrap_arpack_module__set__mceigh

subroutine f90wrap_arpack_module__get__mcapps(f90wrap_mcapps)
    use arpack_module, only: arpack_module_mcapps => mcapps
    implicit none
    integer, intent(out) :: f90wrap_mcapps
    
    f90wrap_mcapps = arpack_module_mcapps
end subroutine f90wrap_arpack_module__get__mcapps

subroutine f90wrap_arpack_module__set__mcapps(f90wrap_mcapps)
    use arpack_module, only: arpack_module_mcapps => mcapps
    implicit none
    integer, intent(in) :: f90wrap_mcapps
    
    arpack_module_mcapps = f90wrap_mcapps
end subroutine f90wrap_arpack_module__set__mcapps

subroutine f90wrap_arpack_module__get__mcgets(f90wrap_mcgets)
    use arpack_module, only: arpack_module_mcgets => mcgets
    implicit none
    integer, intent(out) :: f90wrap_mcgets
    
    f90wrap_mcgets = arpack_module_mcgets
end subroutine f90wrap_arpack_module__get__mcgets

subroutine f90wrap_arpack_module__set__mcgets(f90wrap_mcgets)
    use arpack_module, only: arpack_module_mcgets => mcgets
    implicit none
    integer, intent(in) :: f90wrap_mcgets
    
    arpack_module_mcgets = f90wrap_mcgets
end subroutine f90wrap_arpack_module__set__mcgets

subroutine f90wrap_arpack_module__get__mceupd(f90wrap_mceupd)
    use arpack_module, only: arpack_module_mceupd => mceupd
    implicit none
    integer, intent(out) :: f90wrap_mceupd
    
    f90wrap_mceupd = arpack_module_mceupd
end subroutine f90wrap_arpack_module__get__mceupd

subroutine f90wrap_arpack_module__set__mceupd(f90wrap_mceupd)
    use arpack_module, only: arpack_module_mceupd => mceupd
    implicit none
    integer, intent(in) :: f90wrap_mceupd
    
    arpack_module_mceupd = f90wrap_mceupd
end subroutine f90wrap_arpack_module__set__mceupd

! End of module arpack_module defined in file Arpack_module.fpp


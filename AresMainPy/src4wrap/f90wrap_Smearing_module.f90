! Module smearing_module defined in file Smearing_module.fpp

subroutine f90wrap_smear_init(nev)
    use smearing_module, only: smear_init
    implicit none
    
    integer(4) :: nev
    call smear_init(nev=nev)
end subroutine f90wrap_smear_init

subroutine f90wrap_destroy_smear
    use smearing_module, only: destroy_smear
    implicit none
    
    call destroy_smear()
end subroutine f90wrap_destroy_smear

subroutine f90wrap_hp(x, ret_hp, n)
    use smearing_module, only: hp
    implicit none
    
    real(8) :: x
    real(8), intent(out) :: ret_hp
    integer(4) :: n
    ret_hp = hp(x=x, n=n)
end subroutine f90wrap_hp

subroutine f90wrap_smearsn(y, y0, w, sn)
    use smearing_module, only: smearsn
    implicit none
    
    real(8), intent(in) :: y
    real(8), intent(in) :: y0
    real(8), intent(in) :: w
    real(8), intent(out) :: sn
    call smearsn(y=y, y0=y0, w=w, sn=sn)
end subroutine f90wrap_smearsn

subroutine f90wrap_fermilevel(ne, nev, nk, wk, eval, wke, fme, ets, n0, n1, n2, n3, n4, n5, n6)
    use smearing_module, only: fermilevel
    implicit none
    
    real(8), intent(in) :: ne
    integer(4), intent(in) :: nev
    integer(4), intent(in) :: nk
    real(8), intent(in), dimension(n0) :: wk
    real(8), intent(in), dimension(n1,n2,n3) :: eval
    real(8), dimension(n4,n5,n6) :: wke
    real(8), intent(out) :: fme
    real(8), intent(out) :: ets
    integer :: n0
    !f2py intent(hide), depend(wk) :: n0 = shape(wk,0)
    integer :: n1
    !f2py intent(hide), depend(eval) :: n1 = shape(eval,0)
    integer :: n2
    !f2py intent(hide), depend(eval) :: n2 = shape(eval,1)
    integer :: n3
    !f2py intent(hide), depend(eval) :: n3 = shape(eval,2)
    integer :: n4
    !f2py intent(hide), depend(wke) :: n4 = shape(wke,0)
    integer :: n5
    !f2py intent(hide), depend(wke) :: n5 = shape(wke,1)
    integer :: n6
    !f2py intent(hide), depend(wke) :: n6 = shape(wke,2)
    call fermilevel(ne=ne, nev=nev, nk=nk, wk=wk, eval=eval, wke=wke, fme=fme, ETS=ets)
end subroutine f90wrap_fermilevel

subroutine f90wrap_fermilevel_iso(ne, nev, nk, wk, eval, wke, fme, ets, n0, n1, n2, n3, n4, n5)
    use smearing_module, only: fermilevel_iso
    implicit none
    
    real(8), intent(in) :: ne
    integer(4), intent(in) :: nev
    integer(4), intent(in) :: nk
    real(8), intent(in), dimension(n0) :: wk
    real(8), intent(in), dimension(n1,n2) :: eval
    real(8), dimension(n3,n4,n5) :: wke
    real(8), intent(out) :: fme
    real(8), intent(out) :: ets
    integer :: n0
    !f2py intent(hide), depend(wk) :: n0 = shape(wk,0)
    integer :: n1
    !f2py intent(hide), depend(eval) :: n1 = shape(eval,0)
    integer :: n2
    !f2py intent(hide), depend(eval) :: n2 = shape(eval,1)
    integer :: n3
    !f2py intent(hide), depend(wke) :: n3 = shape(wke,0)
    integer :: n4
    !f2py intent(hide), depend(wke) :: n4 = shape(wke,1)
    integer :: n5
    !f2py intent(hide), depend(wke) :: n5 = shape(wke,2)
    call fermilevel_iso(ne=ne, nev=nev, nk=nk, wk=wk, eval=eval, wke=wke, fme=fme, ETS=ets)
end subroutine f90wrap_fermilevel_iso

subroutine f90wrap_enpy(e_mu, ret_enpy, fi)
    use smearing_module, only: enpy
    implicit none
    
    real(8), intent(in) :: e_mu
    real, intent(out) :: ret_enpy
    real(8), intent(in) :: fi
    ret_enpy = enpy(E_mu=e_mu, fi=fi)
end subroutine f90wrap_enpy

subroutine f90wrap_whg(x, ret_whg, n)
    use smearing_module, only: whg
    implicit none
    
    real(8), intent(in) :: x
    real(8), intent(out) :: ret_whg
    integer(4), intent(in) :: n
    ret_whg = whg(x=x, n=n)
end subroutine f90wrap_whg

subroutine f90wrap_smearing_module__array__wke(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters
    use smearing_module, only: smearing_module_wke => wke
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    if (allocated(smearing_module_wke)) then
        dshape(1:3) = shape(smearing_module_wke)
        dloc = loc(smearing_module_wke)
    else
        dloc = 0
    end if
end subroutine f90wrap_smearing_module__array__wke

subroutine f90wrap_smearing_module__array__subwke(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters
    use smearing_module, only: smearing_module_subwke => subwke
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 4
    dtype = 12
    if (allocated(smearing_module_subwke)) then
        dshape(1:4) = shape(smearing_module_subwke)
        dloc = loc(smearing_module_subwke)
    else
        dloc = 0
    end if
end subroutine f90wrap_smearing_module__array__subwke

! End of module smearing_module defined in file Smearing_module.fpp


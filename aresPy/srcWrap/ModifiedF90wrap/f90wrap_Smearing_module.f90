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

subroutine f90wrap_fermilevel(ne, nev, nk, wk, eval, sigma, n0, n1, n2, n3)
    use smearing_module, only: fermilevel
    implicit none
    
    real(8), intent(in) :: ne
    integer(4), intent(in) :: nev
    integer(4), intent(in) :: nk
    real(8), intent(in), dimension(n0) :: wk
    real(8), intent(in), dimension(n1,n2,n3) :: eval
    real(8), intent(in) :: sigma
    integer :: n0
    !f2py intent(hide), depend(wk) :: n0 = shape(wk,0)
    integer :: n1
    !f2py intent(hide), depend(eval) :: n1 = shape(eval,0)
    integer :: n2
    !f2py intent(hide), depend(eval) :: n2 = shape(eval,1)
    integer :: n3
    !f2py intent(hide), depend(eval) :: n3 = shape(eval,2)
    call fermilevel(ne=ne, nev=nev, nk=nk, wk=wk, eval=eval, sigma=sigma)
end subroutine f90wrap_fermilevel

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

subroutine f90wrap_updaterho_pbc(nps, nev, eig, wk, rhos, rho, n0, n1, n2, n3, n4, n5)
    use grid_module, only: eigen_type
    use smearing_module, only: updaterho_pbc
    implicit none
    
    type eigen_type_ptr_type
        type(eigen_type), pointer :: p => NULL()
    end type eigen_type_ptr_type
    integer(4), intent(in) :: nps
    integer(4), intent(in) :: nev
    type(eigen_type_ptr_type) :: eig_ptr
    integer, intent(in), dimension(2) :: eig
    real(8), intent(in), dimension(n0,n1,n2) :: wk
    real(8), intent(inout), dimension(n3,n4) :: rhos
    real(8), intent(inout), dimension(n5) :: rho
    integer :: n0
    !f2py intent(hide), depend(wk) :: n0 = shape(wk,0)
    integer :: n1
    !f2py intent(hide), depend(wk) :: n1 = shape(wk,1)
    integer :: n2
    !f2py intent(hide), depend(wk) :: n2 = shape(wk,2)
    integer :: n3
    !f2py intent(hide), depend(rhos) :: n3 = shape(rhos,0)
    integer :: n4
    !f2py intent(hide), depend(rhos) :: n4 = shape(rhos,1)
    integer :: n5
    !f2py intent(hide), depend(rho) :: n5 = shape(rho,0)
    eig_ptr = transfer(eig, eig_ptr)
    call updaterho_pbc(nps=nps, nev=nev, eig=eig_ptr%p, wk=wk, rhoS=rhos, rho=rho)
end subroutine f90wrap_updaterho_pbc

subroutine f90wrap_smear_updaterho(nps, nev, ne, eig, rhos, rho, n0, n1, n2)
    use grid_module, only: eigen_type
    use smearing_module, only: smear_updaterho
    implicit none
    
    type eigen_type_ptr_type
        type(eigen_type), pointer :: p => NULL()
    end type eigen_type_ptr_type
    integer(4), intent(in) :: nps
    integer(4), intent(in) :: nev
    real(8), intent(in) :: ne
    type(eigen_type_ptr_type) :: eig_ptr
    integer, intent(in), dimension(2) :: eig
    real(8), intent(inout), dimension(n0,n1) :: rhos
    real(8), intent(inout), dimension(n2) :: rho
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(rho) :: n2 = shape(rho,0)
    eig_ptr = transfer(eig, eig_ptr)
    call smear_updaterho(nps=nps, nev=nev, ne=ne, eig=eig_ptr%p, rhoS=rhos, rho=rho)
end subroutine f90wrap_smear_updaterho

subroutine f90wrap_smearing_module__array__wke(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters
    use smpi_math_module
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

subroutine f90wrap_smearing_module__get__fme(f90wrap_fme)
    use smearing_module, only: smearing_module_fme => fme
    implicit none
    real(8), intent(out) :: f90wrap_fme
    
    f90wrap_fme = smearing_module_fme
end subroutine f90wrap_smearing_module__get__fme

subroutine f90wrap_smearing_module__set__fme(f90wrap_fme)
    use smearing_module, only: smearing_module_fme => fme
    implicit none
    real(8), intent(in) :: f90wrap_fme
    
    smearing_module_fme = f90wrap_fme
end subroutine f90wrap_smearing_module__set__fme

subroutine f90wrap_smearing_module__get__ets(f90wrap_ets)
    use smearing_module, only: smearing_module_ets => ets
    implicit none
    real(8), intent(out) :: f90wrap_ets
    
    f90wrap_ets = smearing_module_ets
end subroutine f90wrap_smearing_module__get__ets

subroutine f90wrap_smearing_module__set__ets(f90wrap_ets)
    use smearing_module, only: smearing_module_ets => ets
    implicit none
    real(8), intent(in) :: f90wrap_ets
    
    smearing_module_ets = f90wrap_ets
end subroutine f90wrap_smearing_module__set__ets

! End of module smearing_module defined in file Smearing_module.fpp


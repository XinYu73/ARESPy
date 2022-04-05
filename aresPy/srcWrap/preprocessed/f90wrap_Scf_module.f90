! Module scf_module defined in file Scf_module.fpp

subroutine f90wrap_electronicscf
    use scf_module, only: electronicscf
    implicit none
    
    call electronicscf()
end subroutine f90wrap_electronicscf

subroutine f90wrap_arpackscf(nps, rhos, rho, eig, n0, n1, n2)
    use grid_module, only: eigen_type
    use scf_module, only: arpackscf
    implicit none
    
    type eigen_type_ptr_type
        type(eigen_type), pointer :: p => NULL()
    end type eigen_type_ptr_type
    integer(4), intent(in) :: nps
    real(8), intent(inout), dimension(n0,n1) :: rhos
    real(8), intent(inout), dimension(n2) :: rho
    type(eigen_type_ptr_type) :: eig_ptr
    integer, intent(in), dimension(2) :: eig
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(rho) :: n2 = shape(rho,0)
    eig_ptr = transfer(eig, eig_ptr)
    call arpackscf(nps=nps, rhoS=rhos, rho=rho, eig=eig_ptr%p)
end subroutine f90wrap_arpackscf

subroutine f90wrap_eigensolver_real(nps, nev, veff, psi, eval, diagtol, n0, n1, n2, n3, n4, n5, n6)
    use scf_module, only: eigensolver_real
    implicit none
    
    integer, intent(in) :: nps
    integer, intent(in) :: nev
    real(8), intent(in), dimension(n0,n1) :: veff
    real(8), intent(inout), dimension(n2,n3,n4) :: psi
    real(8), intent(inout), dimension(n5,n6) :: eval
    real(8), intent(in) :: diagtol
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(psi) :: n2 = shape(psi,0)
    integer :: n3
    !f2py intent(hide), depend(psi) :: n3 = shape(psi,1)
    integer :: n4
    !f2py intent(hide), depend(psi) :: n4 = shape(psi,2)
    integer :: n5
    !f2py intent(hide), depend(eval) :: n5 = shape(eval,0)
    integer :: n6
    !f2py intent(hide), depend(eval) :: n6 = shape(eval,1)
    call eigensolver_real(nps=nps, nev=nev, veff=veff, psi=psi, eval=eval, diagTOL=diagtol)
end subroutine f90wrap_eigensolver_real

subroutine f90wrap_chefsi(nps, rhos, rho, eig, n0, n1, n2)
    use grid_module, only: eigen_type
    use scf_module, only: chefsi
    implicit none
    
    type eigen_type_ptr_type
        type(eigen_type), pointer :: p => NULL()
    end type eigen_type_ptr_type
    integer(4), intent(in) :: nps
    real(8), intent(inout), dimension(n0,n1) :: rhos
    real(8), intent(inout), dimension(n2) :: rho
    type(eigen_type_ptr_type) :: eig_ptr
    integer, intent(in), dimension(2) :: eig
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(rho) :: n2 = shape(rho,0)
    eig_ptr = transfer(eig, eig_ptr)
    call chefsi(nps=nps, rhoS=rhos, rho=rho, eig=eig_ptr%p)
end subroutine f90wrap_chefsi

subroutine f90wrap_filter_spin_gamma(nps, nev, veff, x, d, n0, n1, n2, n3, n4, n5, n6)
    use scf_module, only: filter_spin_gamma
    implicit none
    
    integer(4), intent(in) :: nps
    integer(4), intent(in) :: nev
    real(8), intent(in), dimension(n0,n1) :: veff
    real(8), intent(inout), dimension(n2,n3,n4) :: x
    real(8), intent(inout), dimension(n5,n6) :: d
    integer :: n0
    !f2py intent(hide), depend(veff) :: n0 = shape(veff,0)
    integer :: n1
    !f2py intent(hide), depend(veff) :: n1 = shape(veff,1)
    integer :: n2
    !f2py intent(hide), depend(x) :: n2 = shape(x,0)
    integer :: n3
    !f2py intent(hide), depend(x) :: n3 = shape(x,1)
    integer :: n4
    !f2py intent(hide), depend(x) :: n4 = shape(x,2)
    integer :: n5
    !f2py intent(hide), depend(d) :: n5 = shape(d,0)
    integer :: n6
    !f2py intent(hide), depend(d) :: n6 = shape(d,1)
    call filter_spin_gamma(nps=nps, nev=nev, veff=veff, X=x, D=d)
end subroutine f90wrap_filter_spin_gamma

subroutine f90wrap_scf_module__get__IWD(f90wrap_IWD)
    use scf_module, only: scf_module_IWD => IWD
    implicit none
    integer(4), intent(out) :: f90wrap_IWD
    
    f90wrap_IWD = scf_module_IWD
end subroutine f90wrap_scf_module__get__IWD

subroutine f90wrap_scf_module__set__IWD(f90wrap_IWD)
    use scf_module, only: scf_module_IWD => IWD
    implicit none
    integer(4), intent(in) :: f90wrap_IWD
    
    scf_module_IWD = f90wrap_IWD
end subroutine f90wrap_scf_module__set__IWD

subroutine f90wrap_scf_module__get__LSCF(f90wrap_LSCF)
    use scf_module, only: scf_module_LSCF => LSCF
    implicit none
    logical, intent(out) :: f90wrap_LSCF
    
    f90wrap_LSCF = scf_module_LSCF
end subroutine f90wrap_scf_module__get__LSCF

subroutine f90wrap_scf_module__set__LSCF(f90wrap_LSCF)
    use scf_module, only: scf_module_LSCF => LSCF
    implicit none
    logical, intent(in) :: f90wrap_LSCF
    
    scf_module_LSCF = f90wrap_LSCF
end subroutine f90wrap_scf_module__set__LSCF

! End of module scf_module defined in file Scf_module.fpp


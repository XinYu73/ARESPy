! Module potential_module defined in file Potential_module.fpp

subroutine f90wrap_cal_veff(rhos, veffs, n0, n1, n2, n3, n4, n5, n6, n7)
    use potential_module, only: cal_veff
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2,n3) :: rhos
    real(8), intent(inout), dimension(n4,n5,n6,n7) :: veffs
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(rhos) :: n2 = shape(rhos,2)
    integer :: n3
    !f2py intent(hide), depend(rhos) :: n3 = shape(rhos,3)
    integer :: n4
    !f2py intent(hide), depend(veffs) :: n4 = shape(veffs,0)
    integer :: n5
    !f2py intent(hide), depend(veffs) :: n5 = shape(veffs,1)
    integer :: n6
    !f2py intent(hide), depend(veffs) :: n6 = shape(veffs,2)
    integer :: n7
    !f2py intent(hide), depend(veffs) :: n7 = shape(veffs,3)
    call cal_veff(rhoS=rhos, veffS=veffs)
end subroutine f90wrap_cal_veff

subroutine f90wrap_cal_veff_iso(rhos, veffs, n0, n1, n2, n3)
    use potential_module, only: cal_veff_iso
    implicit none
    
    real(8), intent(in), dimension(n0,n1) :: rhos
    real(8), intent(inout), dimension(n2,n3) :: veffs
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(veffs) :: n2 = shape(veffs,0)
    integer :: n3
    !f2py intent(hide), depend(veffs) :: n3 = shape(veffs,1)
    call cal_veff_iso(rhoS=rhos, veffS=veffs)
end subroutine f90wrap_cal_veff_iso

subroutine f90wrap_vhartree(rho, vhart, n0, n1, n2, n3, n4, n5)
    use potential_module, only: vhartree
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: rho
    real(8), intent(inout), dimension(n3,n4,n5) :: vhart
    integer :: n0
    !f2py intent(hide), depend(rho) :: n0 = shape(rho,0)
    integer :: n1
    !f2py intent(hide), depend(rho) :: n1 = shape(rho,1)
    integer :: n2
    !f2py intent(hide), depend(rho) :: n2 = shape(rho,2)
    integer :: n3
    !f2py intent(hide), depend(vhart) :: n3 = shape(vhart,0)
    integer :: n4
    !f2py intent(hide), depend(vhart) :: n4 = shape(vhart,1)
    integer :: n5
    !f2py intent(hide), depend(vhart) :: n5 = shape(vhart,2)
    call vhartree(rho=rho, vhart=vhart)
end subroutine f90wrap_vhartree

subroutine f90wrap_vlpp
    use potential_module, only: vlpp
    implicit none
    
    call vlpp()
end subroutine f90wrap_vlpp

subroutine f90wrap_vlpp_real
    use potential_module, only: vlpp_real
    implicit none
    
    call vlpp_real()
end subroutine f90wrap_vlpp_real

subroutine f90wrap_vlda(rho, ldapotential, n0, n1, n2, n3, n4, n5)
    use potential_module, only: vlda
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: rho
    real(8), intent(inout), dimension(n3,n4,n5) :: ldapotential
    integer :: n0
    !f2py intent(hide), depend(rho) :: n0 = shape(rho,0)
    integer :: n1
    !f2py intent(hide), depend(rho) :: n1 = shape(rho,1)
    integer :: n2
    !f2py intent(hide), depend(rho) :: n2 = shape(rho,2)
    integer :: n3
    !f2py intent(hide), depend(ldapotential) :: n3 = shape(ldapotential,0)
    integer :: n4
    !f2py intent(hide), depend(ldapotential) :: n4 = shape(ldapotential,1)
    integer :: n5
    !f2py intent(hide), depend(ldapotential) :: n5 = shape(ldapotential,2)
    call vlda(rho=rho, LDAPotential=ldapotential)
end subroutine f90wrap_vlda

subroutine f90wrap_sumrhos(rhos, rho, n0, n1, n2, n3, n4, n5, n6)
    use potential_module, only: sumrhos
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2,n3) :: rhos
    real(8), intent(inout), dimension(n4,n5,n6) :: rho
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(rhos) :: n2 = shape(rhos,2)
    integer :: n3
    !f2py intent(hide), depend(rhos) :: n3 = shape(rhos,3)
    integer :: n4
    !f2py intent(hide), depend(rho) :: n4 = shape(rho,0)
    integer :: n5
    !f2py intent(hide), depend(rho) :: n5 = shape(rho,1)
    integer :: n6
    !f2py intent(hide), depend(rho) :: n6 = shape(rho,2)
    call sumrhos(rhoS=rhos, rho=rho)
end subroutine f90wrap_sumrhos

subroutine f90wrap_sumrhos_iso(rhos, rho, n0, n1, n2)
    use potential_module, only: sumrhos_iso
    implicit none
    
    real(8), intent(in), dimension(n0,n1) :: rhos
    real(8), intent(inout), dimension(n2) :: rho
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(rho) :: n2 = shape(rho,0)
    call sumrhos_iso(rhoS=rhos, rho=rho)
end subroutine f90wrap_sumrhos_iso

subroutine f90wrap_potential_module__array__V_accelerate(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use potential_module, only: potential_module_v_accelerate => v_accelerate
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(potential_module_V_accelerate)) then
        dshape(1:1) = shape(potential_module_V_accelerate)
        dloc = loc(potential_module_V_accelerate)
    else
        dloc = 0
    end if
end subroutine f90wrap_potential_module__array__V_accelerate

subroutine f90wrap_potential_module__array__V_hxc_old(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use potential_module, only: potential_module_v_hxc_old => v_hxc_old
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    if (allocated(potential_module_V_hxc_old)) then
        dshape(1:3) = shape(potential_module_V_hxc_old)
        dloc = loc(potential_module_V_hxc_old)
    else
        dloc = 0
    end if
end subroutine f90wrap_potential_module__array__V_hxc_old

subroutine f90wrap_potential_module__array__V_hxc_new(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use potential_module, only: potential_module_v_hxc_new => v_hxc_new
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    if (allocated(potential_module_V_hxc_new)) then
        dshape(1:3) = shape(potential_module_V_hxc_new)
        dloc = loc(potential_module_V_hxc_new)
    else
        dloc = 0
    end if
end subroutine f90wrap_potential_module__array__V_hxc_new

subroutine f90wrap_potential_module__get__ACCELERATE(f90wrap_ACCELERATE)
    use potential_module, only: potential_module_ACCELERATE => ACCELERATE
    implicit none
    logical, intent(out) :: f90wrap_ACCELERATE
    
    f90wrap_ACCELERATE = potential_module_ACCELERATE
end subroutine f90wrap_potential_module__get__ACCELERATE

subroutine f90wrap_potential_module__set__ACCELERATE(f90wrap_ACCELERATE)
    use potential_module, only: potential_module_ACCELERATE => ACCELERATE
    implicit none
    logical, intent(in) :: f90wrap_ACCELERATE
    
    potential_module_ACCELERATE = f90wrap_ACCELERATE
end subroutine f90wrap_potential_module__set__ACCELERATE

! End of module potential_module defined in file Potential_module.fpp


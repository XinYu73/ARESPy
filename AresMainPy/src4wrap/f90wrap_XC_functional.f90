! Module libxc_module defined in file XC_functional.fpp

subroutine f90wrap_ldalib_energy(rhos, exc, n0, n1, n2, n3)
    use libxc_module, only: ldalib_energy
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2,n3) :: rhos
    real(8), intent(out) :: exc
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(rhos) :: n2 = shape(rhos,2)
    integer :: n3
    !f2py intent(hide), depend(rhos) :: n3 = shape(rhos,3)
    call ldalib_energy(rhoS=rhos, Exc=exc)
end subroutine f90wrap_ldalib_energy

subroutine f90wrap_ldalib_potential(rhos, vxcs, n0, n1, n2, n3, n4, n5, n6, n7)
    use libxc_module, only: ldalib_potential
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2,n3) :: rhos
    real(8), intent(inout), dimension(n4,n5,n6,n7) :: vxcs
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(rhos) :: n2 = shape(rhos,2)
    integer :: n3
    !f2py intent(hide), depend(rhos) :: n3 = shape(rhos,3)
    integer :: n4
    !f2py intent(hide), depend(vxcs) :: n4 = shape(vxcs,0)
    integer :: n5
    !f2py intent(hide), depend(vxcs) :: n5 = shape(vxcs,1)
    integer :: n6
    !f2py intent(hide), depend(vxcs) :: n6 = shape(vxcs,2)
    integer :: n7
    !f2py intent(hide), depend(vxcs) :: n7 = shape(vxcs,3)
    call ldalib_potential(rhoS=rhos, vxcS=vxcs)
end subroutine f90wrap_ldalib_potential

subroutine f90wrap_ldalib_energy_iso(rhos, exc, n0, n1)
    use libxc_module, only: ldalib_energy_iso
    implicit none
    
    real(8), intent(in), dimension(n0,n1) :: rhos
    real(8), intent(out) :: exc
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    call ldalib_energy_iso(rhoS=rhos, Exc=exc)
end subroutine f90wrap_ldalib_energy_iso

subroutine f90wrap_ldalib_potential_iso(rhos, vxcs, n0, n1, n2, n3)
    use libxc_module, only: ldalib_potential_iso
    implicit none
    
    real(8), intent(in), dimension(n0,n1) :: rhos
    real(8), intent(inout), dimension(n2,n3) :: vxcs
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(vxcs) :: n2 = shape(vxcs,0)
    integer :: n3
    !f2py intent(hide), depend(vxcs) :: n3 = shape(vxcs,1)
    call ldalib_potential_iso(rhoS=rhos, vxcS=vxcs)
end subroutine f90wrap_ldalib_potential_iso

! End of module libxc_module defined in file XC_functional.fpp


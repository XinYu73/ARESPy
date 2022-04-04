! Module potential_module defined in file Potential_module.f90

subroutine f90wrap_calveff(nps, rhos, rho, veffs, n0, n1, n2, n3, n4)
    use potential_module, only: calveff
    implicit none
    
    integer(8), intent(in) :: nps
    real(8), intent(in), dimension(n0,n1) :: rhos
    real(8), intent(in), dimension(n2) :: rho
    real(8), intent(inout), dimension(n3,n4) :: veffs
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(rho) :: n2 = shape(rho,0)
    integer :: n3
    !f2py intent(hide), depend(veffs) :: n3 = shape(veffs,0)
    integer :: n4
    !f2py intent(hide), depend(veffs) :: n4 = shape(veffs,1)
    call calveff(nps=nps, rhoS=rhos, rho=rho, veffS=veffs)
end subroutine f90wrap_calveff

subroutine f90wrap_vhartree(nps, rho, vhart, n0, n1)
    use potential_module, only: vhartree
    implicit none
    
    integer(8), intent(in) :: nps
    real(8), intent(in), dimension(n0) :: rho
    real(8), intent(inout), dimension(n1) :: vhart
    integer :: n0
    !f2py intent(hide), depend(rho) :: n0 = shape(rho,0)
    integer :: n1
    !f2py intent(hide), depend(vhart) :: n1 = shape(vhart,0)
    call vhartree(nps=nps, rho=rho, vhart=vhart)
end subroutine f90wrap_vhartree

subroutine f90wrap_vlpp
    use potential_module, only: vlpp
    implicit none
    
    call vlpp()
end subroutine f90wrap_vlpp

subroutine f90wrap_potential_module__array__V_accelerate(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use potential_module, only: potential_module_v_accelerate => v_accelerate
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    if (allocated(potential_module_V_accelerate)) then
        dshape(1:3) = shape(potential_module_V_accelerate)
        dloc = loc(potential_module_V_accelerate)
    else
        dloc = 0
    end if
end subroutine f90wrap_potential_module__array__V_accelerate

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

! End of module potential_module defined in file Potential_module.f90


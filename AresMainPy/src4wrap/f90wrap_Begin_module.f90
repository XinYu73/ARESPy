! Module begin_module defined in file Begin_module.fpp

subroutine f90wrap_initial_grid
    use begin_module, only: initial_grid
    implicit none
    
    call initial_grid()
end subroutine f90wrap_initial_grid

subroutine f90wrap_initial_grid_per
    use begin_module, only: initial_grid_per
    implicit none
    
    call initial_grid_per()
end subroutine f90wrap_initial_grid_per

subroutine f90wrap_initial_grid_iso
    use begin_module, only: initial_grid_iso
    implicit none
    
    call initial_grid_iso()
end subroutine f90wrap_initial_grid_iso

subroutine f90wrap_initial_density
    use begin_module, only: initial_density
    implicit none
    
    call initial_density()
end subroutine f90wrap_initial_density

subroutine f90wrap_inichrg
    use begin_module, only: inichrg
    implicit none
    
    call inichrg()
end subroutine f90wrap_inichrg

subroutine f90wrap_init_iso_density
    use begin_module, only: init_iso_density
    implicit none
    
    call init_iso_density()
end subroutine f90wrap_init_iso_density

subroutine f90wrap_inichrg_iso
    use begin_module, only: inichrg_iso
    implicit none
    
    call inichrg_iso()
end subroutine f90wrap_inichrg_iso

subroutine f90wrap_init_mo_density_iso
    use begin_module, only: init_mo_density_iso
    implicit none
    
    call init_mo_density_iso()
end subroutine f90wrap_init_mo_density_iso

subroutine f90wrap_iso_init_mo(initx_sto, n0, n1)
    use begin_module, only: iso_init_mo
    implicit none
    
    real(8), intent(inout), dimension(n0,n1) :: initx_sto
    integer :: n0
    !f2py intent(hide), depend(initx_sto) :: n0 = shape(initx_sto,0)
    integer :: n1
    !f2py intent(hide), depend(initx_sto) :: n1 = shape(initx_sto,1)
    call iso_init_mo(initX_sto=initx_sto)
end subroutine f90wrap_iso_init_mo

! End of module begin_module defined in file Begin_module.fpp


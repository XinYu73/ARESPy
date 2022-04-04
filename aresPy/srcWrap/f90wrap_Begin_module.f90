! Module begin_module defined in file Begin_module.f90

subroutine f90wrap_initial_grid_pbc
    use begin_module, only: initial_grid_pbc
    implicit none
    
    call initial_grid_pbc()
end subroutine f90wrap_initial_grid_pbc

subroutine f90wrap_initial_density
    use begin_module, only: initial_density
    implicit none
    
    call initial_density()
end subroutine f90wrap_initial_density

subroutine f90wrap_inichrg_sp
    use begin_module, only: inichrg_sp
    implicit none
    
    call inichrg_sp()
end subroutine f90wrap_inichrg_sp

! End of module begin_module defined in file Begin_module.f90


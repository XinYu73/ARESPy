! Module energy_module defined in file Energy_module.fpp

subroutine f90wrap_totalenergy(nps, eig, rhos, rho, n0, n1, n2)
    use grid_module, only: eigen_type
    use energy_module, only: totalenergy
    implicit none
    
    type eigen_type_ptr_type
        type(eigen_type), pointer :: p => NULL()
    end type eigen_type_ptr_type
    integer(4), intent(in) :: nps
    type(eigen_type_ptr_type) :: eig_ptr
    integer, intent(in), dimension(2) :: eig
    real(8), intent(in), dimension(n0,n1) :: rhos
    real(8), intent(in), dimension(n2) :: rho
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(rho) :: n2 = shape(rho,0)
    eig_ptr = transfer(eig, eig_ptr)
    call totalenergy(nps=nps, eig=eig_ptr%p, rhoS=rhos, rho=rho)
end subroutine f90wrap_totalenergy

subroutine f90wrap_ebands(eval, wke, eband, n0, n1, n2, n3, n4, n5)
    use energy_module, only: ebands
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: eval
    real(8), intent(in), dimension(n3,n4,n5) :: wke
    real(8), intent(out) :: eband
    integer :: n0
    !f2py intent(hide), depend(eval) :: n0 = shape(eval,0)
    integer :: n1
    !f2py intent(hide), depend(eval) :: n1 = shape(eval,1)
    integer :: n2
    !f2py intent(hide), depend(eval) :: n2 = shape(eval,2)
    integer :: n3
    !f2py intent(hide), depend(wke) :: n3 = shape(wke,0)
    integer :: n4
    !f2py intent(hide), depend(wke) :: n4 = shape(wke,1)
    integer :: n5
    !f2py intent(hide), depend(wke) :: n5 = shape(wke,2)
    call ebands(eval=eval, wke=wke, Eband=eband)
end subroutine f90wrap_ebands

subroutine f90wrap_energy_module__get__Etot(f90wrap_Etot)
    use energy_module, only: energy_module_Etot => Etot
    implicit none
    real(8), intent(out) :: f90wrap_Etot
    
    f90wrap_Etot = energy_module_Etot
end subroutine f90wrap_energy_module__get__Etot

subroutine f90wrap_energy_module__set__Etot(f90wrap_Etot)
    use energy_module, only: energy_module_Etot => Etot
    implicit none
    real(8), intent(in) :: f90wrap_Etot
    
    energy_module_Etot = f90wrap_Etot
end subroutine f90wrap_energy_module__set__Etot

subroutine f90wrap_energy_module__get__Eband(f90wrap_Eband)
    use energy_module, only: energy_module_Eband => Eband
    implicit none
    real(8), intent(out) :: f90wrap_Eband
    
    f90wrap_Eband = energy_module_Eband
end subroutine f90wrap_energy_module__get__Eband

subroutine f90wrap_energy_module__set__Eband(f90wrap_Eband)
    use energy_module, only: energy_module_Eband => Eband
    implicit none
    real(8), intent(in) :: f90wrap_Eband
    
    energy_module_Eband = f90wrap_Eband
end subroutine f90wrap_energy_module__set__Eband

subroutine f90wrap_energy_module__get__Eh(f90wrap_Eh)
    use energy_module, only: energy_module_Eh => Eh
    implicit none
    real(8), intent(out) :: f90wrap_Eh
    
    f90wrap_Eh = energy_module_Eh
end subroutine f90wrap_energy_module__get__Eh

subroutine f90wrap_energy_module__set__Eh(f90wrap_Eh)
    use energy_module, only: energy_module_Eh => Eh
    implicit none
    real(8), intent(in) :: f90wrap_Eh
    
    energy_module_Eh = f90wrap_Eh
end subroutine f90wrap_energy_module__set__Eh

subroutine f90wrap_energy_module__get__Eext(f90wrap_Eext)
    use energy_module, only: energy_module_Eext => Eext
    implicit none
    real(8), intent(out) :: f90wrap_Eext
    
    f90wrap_Eext = energy_module_Eext
end subroutine f90wrap_energy_module__get__Eext

subroutine f90wrap_energy_module__set__Eext(f90wrap_Eext)
    use energy_module, only: energy_module_Eext => Eext
    implicit none
    real(8), intent(in) :: f90wrap_Eext
    
    energy_module_Eext = f90wrap_Eext
end subroutine f90wrap_energy_module__set__Eext

subroutine f90wrap_energy_module__get__Exc(f90wrap_Exc)
    use energy_module, only: energy_module_Exc => Exc
    implicit none
    real(8), intent(out) :: f90wrap_Exc
    
    f90wrap_Exc = energy_module_Exc
end subroutine f90wrap_energy_module__get__Exc

subroutine f90wrap_energy_module__set__Exc(f90wrap_Exc)
    use energy_module, only: energy_module_Exc => Exc
    implicit none
    real(8), intent(in) :: f90wrap_Exc
    
    energy_module_Exc = f90wrap_Exc
end subroutine f90wrap_energy_module__set__Exc

subroutine f90wrap_energy_module__get__Eele(f90wrap_Eele)
    use energy_module, only: energy_module_Eele => Eele
    implicit none
    real(8), intent(out) :: f90wrap_Eele
    
    f90wrap_Eele = energy_module_Eele
end subroutine f90wrap_energy_module__get__Eele

subroutine f90wrap_energy_module__set__Eele(f90wrap_Eele)
    use energy_module, only: energy_module_Eele => Eele
    implicit none
    real(8), intent(in) :: f90wrap_Eele
    
    energy_module_Eele = f90wrap_Eele
end subroutine f90wrap_energy_module__set__Eele

subroutine f90wrap_energy_module__get__FE(f90wrap_FE)
    use energy_module, only: energy_module_FE => FE
    implicit none
    real(8), intent(out) :: f90wrap_FE
    
    f90wrap_FE = energy_module_FE
end subroutine f90wrap_energy_module__get__FE

subroutine f90wrap_energy_module__set__FE(f90wrap_FE)
    use energy_module, only: energy_module_FE => FE
    implicit none
    real(8), intent(in) :: f90wrap_FE
    
    energy_module_FE = f90wrap_FE
end subroutine f90wrap_energy_module__set__FE

subroutine f90wrap_energy_module__get__FE0(f90wrap_FE0)
    use energy_module, only: energy_module_FE0 => FE0
    implicit none
    real(8), intent(out) :: f90wrap_FE0
    
    f90wrap_FE0 = energy_module_FE0
end subroutine f90wrap_energy_module__get__FE0

subroutine f90wrap_energy_module__set__FE0(f90wrap_FE0)
    use energy_module, only: energy_module_FE0 => FE0
    implicit none
    real(8), intent(in) :: f90wrap_FE0
    
    energy_module_FE0 = f90wrap_FE0
end subroutine f90wrap_energy_module__set__FE0

! End of module energy_module defined in file Energy_module.fpp

! Module output_module defined in file Energy_module.fpp

subroutine f90wrap_output
    use output_module, only: output
    implicit none
    
    call output()
end subroutine f90wrap_output

subroutine f90wrap_write_density
    use output_module, only: write_density
    implicit none
    
    call write_density()
end subroutine f90wrap_write_density

subroutine f90wrap_write_band
    use output_module, only: write_band
    implicit none
    
    call write_band()
end subroutine f90wrap_write_band

subroutine f90wrap_output_module__get__time_total0(f90wrap_time_total0)
    use output_module, only: output_module_time_total0 => time_total0
    implicit none
    integer(4), intent(out) :: f90wrap_time_total0
    
    f90wrap_time_total0 = output_module_time_total0
end subroutine f90wrap_output_module__get__time_total0

subroutine f90wrap_output_module__set__time_total0(f90wrap_time_total0)
    use output_module, only: output_module_time_total0 => time_total0
    implicit none
    integer(4), intent(in) :: f90wrap_time_total0
    
    output_module_time_total0 = f90wrap_time_total0
end subroutine f90wrap_output_module__set__time_total0

subroutine f90wrap_output_module__get__time_total1(f90wrap_time_total1)
    use output_module, only: output_module_time_total1 => time_total1
    implicit none
    integer(4), intent(out) :: f90wrap_time_total1
    
    f90wrap_time_total1 = output_module_time_total1
end subroutine f90wrap_output_module__get__time_total1

subroutine f90wrap_output_module__set__time_total1(f90wrap_time_total1)
    use output_module, only: output_module_time_total1 => time_total1
    implicit none
    integer(4), intent(in) :: f90wrap_time_total1
    
    output_module_time_total1 = f90wrap_time_total1
end subroutine f90wrap_output_module__set__time_total1

subroutine f90wrap_output_module__get__time_scf0(f90wrap_time_scf0)
    use output_module, only: output_module_time_scf0 => time_scf0
    implicit none
    integer(4), intent(out) :: f90wrap_time_scf0
    
    f90wrap_time_scf0 = output_module_time_scf0
end subroutine f90wrap_output_module__get__time_scf0

subroutine f90wrap_output_module__set__time_scf0(f90wrap_time_scf0)
    use output_module, only: output_module_time_scf0 => time_scf0
    implicit none
    integer(4), intent(in) :: f90wrap_time_scf0
    
    output_module_time_scf0 = f90wrap_time_scf0
end subroutine f90wrap_output_module__set__time_scf0

subroutine f90wrap_output_module__get__time_scf1(f90wrap_time_scf1)
    use output_module, only: output_module_time_scf1 => time_scf1
    implicit none
    integer(4), intent(out) :: f90wrap_time_scf1
    
    f90wrap_time_scf1 = output_module_time_scf1
end subroutine f90wrap_output_module__get__time_scf1

subroutine f90wrap_output_module__set__time_scf1(f90wrap_time_scf1)
    use output_module, only: output_module_time_scf1 => time_scf1
    implicit none
    integer(4), intent(in) :: f90wrap_time_scf1
    
    output_module_time_scf1 = f90wrap_time_scf1
end subroutine f90wrap_output_module__set__time_scf1

! End of module output_module defined in file Energy_module.fpp


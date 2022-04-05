! Module parameters defined in file parameters.fpp

subroutine f90wrap_parameters__get__idp(f90wrap_idp)
    use parameters, only: parameters_idp => idp
    implicit none
    integer, intent(out) :: f90wrap_idp
    
    f90wrap_idp = parameters_idp
end subroutine f90wrap_parameters__get__idp

subroutine f90wrap_parameters__get__isp(f90wrap_isp)
    use parameters, only: parameters_isp => isp
    implicit none
    integer, intent(out) :: f90wrap_isp
    
    f90wrap_isp = parameters_isp
end subroutine f90wrap_parameters__get__isp

! End of module parameters defined in file parameters.fpp


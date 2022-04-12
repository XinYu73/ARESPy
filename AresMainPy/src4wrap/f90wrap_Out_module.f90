! Module out_module defined in file Out_module.fpp

subroutine f90wrap_ksout_list
    use out_module, only: ksout_list
    implicit none
    
    call ksout_list()
end subroutine f90wrap_ksout_list

subroutine f90wrap_ksout_list_per
    use out_module, only: ksout_list_per
    implicit none
    
    call ksout_list_per()
end subroutine f90wrap_ksout_list_per

subroutine f90wrap_ksout_list_iso
    use out_module, only: ksout_list_iso
    implicit none
    
    call ksout_list_iso()
end subroutine f90wrap_ksout_list_iso

! End of module out_module defined in file Out_module.fpp


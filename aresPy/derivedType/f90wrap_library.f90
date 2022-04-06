! Module library defined in file library.fpp

subroutine f90wrap_return_value_func(ret_val_out, val_in)
    use library, only: return_value_func
    implicit none
    
    integer(4), intent(out) :: ret_val_out
    integer(4) :: val_in
    ret_val_out = return_value_func(val_in=val_in)
end subroutine f90wrap_return_value_func

subroutine f90wrap_subfunc(ret_subfunc, a)
    use library, only: subfunc
    implicit none
    
    integer, intent(out) :: ret_subfunc
    integer, intent(in) :: a
    ret_subfunc = subfunc(a=a)
end subroutine f90wrap_subfunc

subroutine f90wrap_return_value_sub(val_in, val_out)
    use library, only: return_value_sub
    implicit none
    
    integer(4), intent(in) :: val_in
    integer(4), intent(out) :: val_out
    call return_value_sub(val_in=val_in, val_out=val_out)
end subroutine f90wrap_return_value_sub

subroutine f90wrap_return_a_dt_func(ret_dt)
    use library, only: return_a_dt_func
    use datatypes, only: different_types
    implicit none
    
    type different_types_ptr_type
        type(different_types), pointer :: p => NULL()
    end type different_types_ptr_type
    type(different_types_ptr_type) :: ret_dt_ptr
    integer, intent(out), dimension(2) :: ret_dt
    allocate(ret_dt_ptr%p)
    ret_dt_ptr%p = return_a_dt_func()
    ret_dt = transfer(ret_dt_ptr, ret_dt)
end subroutine f90wrap_return_a_dt_func

subroutine f90wrap_do_array_stuff(n, x, y, br, co, n0, n1, n2, n3)
    use library, only: do_array_stuff
    implicit none
    
    integer, intent(in) :: n
    real(8), intent(in), dimension(n0) :: x
    real(8), intent(in), dimension(n1) :: y
    real(8), intent(inout), dimension(n2) :: br
    real(8), intent(inout), dimension(4,n3) :: co
    integer :: n0
    !f2py intent(hide), depend(x) :: n0 = shape(x,0)
    integer :: n1
    !f2py intent(hide), depend(y) :: n1 = shape(y,0)
    integer :: n2
    !f2py intent(hide), depend(br) :: n2 = shape(br,0)
    integer :: n3
    !f2py intent(hide), depend(co) :: n3 = shape(co,1)
    call do_array_stuff(n=n, x=x, y=y, br=br, co=co)
end subroutine f90wrap_do_array_stuff

subroutine f90wrap_only_manipulate(n, array, n0)
    use library, only: only_manipulate
    implicit none
    
    integer, intent(in) :: n
    real(8), intent(inout), dimension(4,n0) :: array
    integer :: n0
    !f2py intent(hide), depend(array) :: n0 = shape(array,1)
    call only_manipulate(n=n, array=array)
end subroutine f90wrap_only_manipulate

subroutine f90wrap_set_derived_type(dt, dt_beta, dt_delta)
    use library, only: set_derived_type
    use datatypes, only: different_types
    implicit none
    
    type different_types_ptr_type
        type(different_types), pointer :: p => NULL()
    end type different_types_ptr_type
    type(different_types_ptr_type) :: dt_ptr
    integer, intent(out), dimension(2) :: dt
    integer(4), intent(in) :: dt_beta
    real(8), intent(in) :: dt_delta
    allocate(dt_ptr%p)
    call set_derived_type(dt=dt_ptr%p, dt_beta=dt_beta, dt_delta=dt_delta)
    dt = transfer(dt_ptr, dt)
end subroutine f90wrap_set_derived_type

subroutine f90wrap_modify_derived_types(dt1, dt2, dt3)
    use library, only: modify_derived_types
    use datatypes, only: different_types
    implicit none
    
    type different_types_ptr_type
        type(different_types), pointer :: p => NULL()
    end type different_types_ptr_type
    type(different_types_ptr_type) :: dt1_ptr
    integer, intent(in), dimension(2) :: dt1
    type(different_types_ptr_type) :: dt2_ptr
    integer, intent(in), dimension(2) :: dt2
    type(different_types_ptr_type) :: dt3_ptr
    integer, intent(in), dimension(2) :: dt3
    dt1_ptr = transfer(dt1, dt1_ptr)
    dt2_ptr = transfer(dt2, dt2_ptr)
    dt3_ptr = transfer(dt3, dt3_ptr)
    call modify_derived_types(dt1=dt1_ptr%p, dt2=dt2_ptr%p, dt3=dt3_ptr%p)
end subroutine f90wrap_modify_derived_types

subroutine f90wrap_modify_dertype_fixed_shape_arrays(dertype)
    use library, only: modify_dertype_fixed_shape_arrays
    use datatypes, only: fixed_shape_arrays
    implicit none
    
    type fixed_shape_arrays_ptr_type
        type(fixed_shape_arrays), pointer :: p => NULL()
    end type fixed_shape_arrays_ptr_type
    type(fixed_shape_arrays_ptr_type) :: dertype_ptr
    integer, intent(out), dimension(2) :: dertype
    allocate(dertype_ptr%p)
    call modify_dertype_fixed_shape_arrays(dertype=dertype_ptr%p)
    dertype = transfer(dertype_ptr, dertype)
end subroutine f90wrap_modify_dertype_fixed_shape_arrays

subroutine f90wrap_return_dertype_pointer_arrays(m, n, dertype)
    use library, only: return_dertype_pointer_arrays
    use datatypes, only: pointer_arrays
    implicit none
    
    type pointer_arrays_ptr_type
        type(pointer_arrays), pointer :: p => NULL()
    end type pointer_arrays_ptr_type
    integer, intent(in) :: m
    integer, intent(in) :: n
    type(pointer_arrays_ptr_type) :: dertype_ptr
    integer, intent(out), dimension(2) :: dertype
    allocate(dertype_ptr%p)
    call return_dertype_pointer_arrays(m=m, n=n, dertype=dertype_ptr%p)
    dertype = transfer(dertype_ptr, dertype)
end subroutine f90wrap_return_dertype_pointer_arrays

subroutine f90wrap_modify_dertype_pointer_arrays(dertype)
    use library, only: modify_dertype_pointer_arrays
    use datatypes, only: pointer_arrays
    implicit none
    
    type pointer_arrays_ptr_type
        type(pointer_arrays), pointer :: p => NULL()
    end type pointer_arrays_ptr_type
    type(pointer_arrays_ptr_type) :: dertype_ptr
    integer, intent(in), dimension(2) :: dertype
    dertype_ptr = transfer(dertype, dertype_ptr)
    call modify_dertype_pointer_arrays(dertype=dertype_ptr%p)
end subroutine f90wrap_modify_dertype_pointer_arrays

subroutine f90wrap_return_dertype_alloc_arrays(m, n, dertype)
    use library, only: return_dertype_alloc_arrays
    use datatypes_allocatable, only: alloc_arrays
    implicit none
    
    type alloc_arrays_ptr_type
        type(alloc_arrays), pointer :: p => NULL()
    end type alloc_arrays_ptr_type
    integer, intent(in) :: m
    integer, intent(in) :: n
    type(alloc_arrays_ptr_type) :: dertype_ptr
    integer, intent(out), dimension(2) :: dertype
    allocate(dertype_ptr%p)
    call return_dertype_alloc_arrays(m=m, n=n, dertype=dertype_ptr%p)
    dertype = transfer(dertype_ptr, dertype)
end subroutine f90wrap_return_dertype_alloc_arrays

subroutine f90wrap_modify_dertype_alloc_arrays(dertype)
    use library, only: modify_dertype_alloc_arrays
    use datatypes_allocatable, only: alloc_arrays
    implicit none
    
    type alloc_arrays_ptr_type
        type(alloc_arrays), pointer :: p => NULL()
    end type alloc_arrays_ptr_type
    type(alloc_arrays_ptr_type) :: dertype_ptr
    integer, intent(in), dimension(2) :: dertype
    dertype_ptr = transfer(dertype, dertype_ptr)
    call modify_dertype_alloc_arrays(dertype=dertype_ptr%p)
end subroutine f90wrap_modify_dertype_alloc_arrays

! End of module library defined in file library.fpp


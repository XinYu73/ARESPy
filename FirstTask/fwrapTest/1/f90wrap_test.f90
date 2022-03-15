! Module module_calcul defined in file test.fpp

subroutine f90wrap_type_ptmes__get__y(this, f90wrap_y)
    use module_calcul, only: type_ptmes
    implicit none
    type type_ptmes_ptr_type
        type(type_ptmes), pointer :: p => NULL()
    end type type_ptmes_ptr_type
    integer, intent(in)   :: this(2)
    type(type_ptmes_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_y
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_y = this_ptr%p%y
end subroutine f90wrap_type_ptmes__get__y

subroutine f90wrap_type_ptmes__set__y(this, f90wrap_y)
    use module_calcul, only: type_ptmes
    implicit none
    type type_ptmes_ptr_type
        type(type_ptmes), pointer :: p => NULL()
    end type type_ptmes_ptr_type
    integer, intent(in)   :: this(2)
    type(type_ptmes_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_y
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%y = f90wrap_y
end subroutine f90wrap_type_ptmes__set__y

subroutine f90wrap_type_ptmes_initialise(this)
    use module_calcul, only: type_ptmes
    implicit none
    
    type type_ptmes_ptr_type
        type(type_ptmes), pointer :: p => NULL()
    end type type_ptmes_ptr_type
    type(type_ptmes_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_type_ptmes_initialise

subroutine f90wrap_type_ptmes_finalise(this)
    use module_calcul, only: type_ptmes
    implicit none
    
    type type_ptmes_ptr_type
        type(type_ptmes), pointer :: p => NULL()
    end type type_ptmes_ptr_type
    type(type_ptmes_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_type_ptmes_finalise

subroutine f90wrap_array_type__array_getitem__x(f90wrap_this, f90wrap_i, xitem)
    
    use module_calcul, only: array_type, type_ptmes
    implicit none
    
    type array_type_ptr_type
        type(array_type), pointer :: p => NULL()
    end type array_type_ptr_type
    type type_ptmes_ptr_type
        type(type_ptmes), pointer :: p => NULL()
    end type type_ptmes_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(array_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: xitem(2)
    type(type_ptmes_ptr_type) :: x_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%x)) then
        call f90wrap_abort("array index out of range")
    else
        x_ptr%p => this_ptr%p%x(f90wrap_i)
        xitem = transfer(x_ptr,xitem)
    endif
end subroutine f90wrap_array_type__array_getitem__x

subroutine f90wrap_array_type__array_setitem__x(f90wrap_this, f90wrap_i, xitem)
    
    use module_calcul, only: array_type, type_ptmes
    implicit none
    
    type array_type_ptr_type
        type(array_type), pointer :: p => NULL()
    end type array_type_ptr_type
    type type_ptmes_ptr_type
        type(type_ptmes), pointer :: p => NULL()
    end type type_ptmes_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(array_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: xitem(2)
    type(type_ptmes_ptr_type) :: x_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%x)) then
        call f90wrap_abort("array index out of range")
    else
        x_ptr = transfer(xitem,x_ptr)
        this_ptr%p%x(f90wrap_i) = x_ptr%p
    endif
end subroutine f90wrap_array_type__array_setitem__x

subroutine f90wrap_array_type__array_len__x(f90wrap_this, f90wrap_n)
    
    use module_calcul, only: array_type, type_ptmes
    implicit none
    
    type array_type_ptr_type
        type(array_type), pointer :: p => NULL()
    end type array_type_ptr_type
    type type_ptmes_ptr_type
        type(type_ptmes), pointer :: p => NULL()
    end type type_ptmes_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(array_type_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    f90wrap_n = size(this_ptr%p%x)
end subroutine f90wrap_array_type__array_len__x

subroutine f90wrap_array_type_initialise(this)
    use module_calcul, only: array_type
    implicit none
    
    type array_type_ptr_type
        type(array_type), pointer :: p => NULL()
    end type array_type_ptr_type
    type(array_type_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_array_type_initialise

subroutine f90wrap_array_type_finalise(this)
    use module_calcul, only: array_type
    implicit none
    
    type array_type_ptr_type
        type(array_type), pointer :: p => NULL()
    end type array_type_ptr_type
    type(array_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_array_type_finalise

subroutine f90wrap_recup_point(x)
    use module_calcul, only: array_type, recup_point
    implicit none
    
    type array_type_ptr_type
        type(array_type), pointer :: p => NULL()
    end type array_type_ptr_type
    type(array_type_ptr_type) :: x_ptr
    integer, intent(in), dimension(2) :: x
    x_ptr = transfer(x, x_ptr)
    call recup_point(x=x_ptr%p)
end subroutine f90wrap_recup_point

subroutine f90wrap_xytest
    use module_calcul, only: xytest
    implicit none
    
    call xytest()
end subroutine f90wrap_xytest

subroutine f90wrap_module_calcul__array_getitem__xarr(dummy_this, f90wrap_i, xarritem)
    
    use module_calcul, only: array_type, module_calcul_xarr => xarr
    implicit none
    
    type array_type_ptr_type
        type(array_type), pointer :: p => NULL()
    end type array_type_ptr_type
    integer, intent(in) :: dummy_this(2)
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: xarritem(2)
    type(array_type_ptr_type) :: xarr_ptr
    
    if (f90wrap_i < 1 .or. f90wrap_i > size(module_calcul_xarr)) then
        call f90wrap_abort("array index out of range")
    else
        xarr_ptr%p => module_calcul_xarr(f90wrap_i)
        xarritem = transfer(xarr_ptr,xarritem)
    endif
end subroutine f90wrap_module_calcul__array_getitem__xarr

subroutine f90wrap_module_calcul__array_setitem__xarr(dummy_this, f90wrap_i, xarritem)
    
    use module_calcul, only: array_type, module_calcul_xarr => xarr
    implicit none
    
    type array_type_ptr_type
        type(array_type), pointer :: p => NULL()
    end type array_type_ptr_type
    integer, intent(in) :: dummy_this(2)
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: xarritem(2)
    type(array_type_ptr_type) :: xarr_ptr
    
    if (f90wrap_i < 1 .or. f90wrap_i > size(module_calcul_xarr)) then
        call f90wrap_abort("array index out of range")
    else
        xarr_ptr = transfer(xarritem,xarr_ptr)
        module_calcul_xarr(f90wrap_i) = xarr_ptr%p
    endif
end subroutine f90wrap_module_calcul__array_setitem__xarr

subroutine f90wrap_module_calcul__array_len__xarr(dummy_this, f90wrap_n)
    
    use module_calcul, only: array_type, module_calcul_xarr => xarr
    implicit none
    
    type array_type_ptr_type
        type(array_type), pointer :: p => NULL()
    end type array_type_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: dummy_this(2)
    f90wrap_n = size(module_calcul_xarr)
end subroutine f90wrap_module_calcul__array_len__xarr

! End of module module_calcul defined in file test.fpp


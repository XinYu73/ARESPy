! Module aresapi defined in file aresAPI.fpp

subroutine f90wrap_aresout__array__forces(this, nd, dtype, dshape, dloc)
    use aresapi, only: aresout
    implicit none
    type aresout_ptr_type
        type(aresout), pointer :: p => NULL()
    end type aresout_ptr_type
    integer, intent(in) :: this(2)
    type(aresout_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%forces)) then
        dshape(1:2) = shape(this_ptr%p%forces)
        dloc = loc(this_ptr%p%forces)
    else
        dloc = 0
    end if
end subroutine f90wrap_aresout__array__forces

subroutine f90wrap_aresout__array__stress(this, nd, dtype, dshape, dloc)
    use aresapi, only: aresout
    implicit none
    type aresout_ptr_type
        type(aresout), pointer :: p => NULL()
    end type aresout_ptr_type
    integer, intent(in) :: this(2)
    type(aresout_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%stress)) then
        dshape(1:2) = shape(this_ptr%p%stress)
        dloc = loc(this_ptr%p%stress)
    else
        dloc = 0
    end if
end subroutine f90wrap_aresout__array__stress

subroutine f90wrap_aresout__array__poscar(this, nd, dtype, dshape, dloc)
    use aresapi, only: aresout
    implicit none
    type aresout_ptr_type
        type(aresout), pointer :: p => NULL()
    end type aresout_ptr_type
    integer, intent(in) :: this(2)
    type(aresout_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%poscar)) then
        dshape(1:2) = shape(this_ptr%p%poscar)
        dloc = loc(this_ptr%p%poscar)
    else
        dloc = 0
    end if
end subroutine f90wrap_aresout__array__poscar

subroutine f90wrap_aresout__array__pos(this, nd, dtype, dshape, dloc)
    use aresapi, only: aresout
    implicit none
    type aresout_ptr_type
        type(aresout), pointer :: p => NULL()
    end type aresout_ptr_type
    integer, intent(in) :: this(2)
    type(aresout_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%pos)) then
        dshape(1:2) = shape(this_ptr%p%pos)
        dloc = loc(this_ptr%p%pos)
    else
        dloc = 0
    end if
end subroutine f90wrap_aresout__array__pos

subroutine f90wrap_aresout__array__chargeRho(this, nd, dtype, dshape, dloc)
    use aresapi, only: aresout
    implicit none
    type aresout_ptr_type
        type(aresout), pointer :: p => NULL()
    end type aresout_ptr_type
    integer, intent(in) :: this(2)
    type(aresout_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 4
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%chargeRho)) then
        dshape(1:4) = shape(this_ptr%p%chargeRho)
        dloc = loc(this_ptr%p%chargeRho)
    else
        dloc = 0
    end if
end subroutine f90wrap_aresout__array__chargeRho

subroutine f90wrap_aresout_initialise(this)
    use aresapi, only: aresout
    implicit none
    
    type aresout_ptr_type
        type(aresout), pointer :: p => NULL()
    end type aresout_ptr_type
    type(aresout_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_aresout_initialise

subroutine f90wrap_aresout_finalise(this)
    use aresapi, only: aresout
    implicit none
    
    type aresout_ptr_type
        type(aresout), pointer :: p => NULL()
    end type aresout_ptr_type
    type(aresout_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_aresout_finalise

subroutine f90wrap_init_alloc_arrays(dertype, nnatom)
    use aresapi, only: init_alloc_arrays, aresout
    implicit none
    
    type aresout_ptr_type
        type(aresout), pointer :: p => NULL()
    end type aresout_ptr_type
    type(aresout_ptr_type) :: dertype_ptr
    integer, intent(in), dimension(2) :: dertype
    integer(4), intent(in) :: nnatom
    dertype_ptr = transfer(dertype, dertype_ptr)
    call init_alloc_arrays(dertype=dertype_ptr%p, nnatom=nnatom)
end subroutine f90wrap_init_alloc_arrays

subroutine f90wrap_destroy_alloc_arrays(dertype)
    use aresapi, only: aresout, destroy_alloc_arrays
    implicit none
    
    type aresout_ptr_type
        type(aresout), pointer :: p => NULL()
    end type aresout_ptr_type
    type(aresout_ptr_type) :: dertype_ptr
    integer, intent(in), dimension(2) :: dertype
    dertype_ptr = transfer(dertype, dertype_ptr)
    call destroy_alloc_arrays(dertype=dertype_ptr%p)
end subroutine f90wrap_destroy_alloc_arrays

! End of module aresapi defined in file aresAPI.fpp


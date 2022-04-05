! Module datatypes_allocatable defined in file datatypes.fpp

subroutine f90wrap_alloc_arrays__array__chi(this, nd, dtype, dshape, dloc)
    use datatypes_allocatable, only: alloc_arrays
    implicit none
    type alloc_arrays_ptr_type
        type(alloc_arrays), pointer :: p => NULL()
    end type alloc_arrays_ptr_type
    integer, intent(in) :: this(2)
    type(alloc_arrays_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%chi)) then
        dshape(1:2) = shape(this_ptr%p%chi)
        dloc = loc(this_ptr%p%chi)
    else
        dloc = 0
    end if
end subroutine f90wrap_alloc_arrays__array__chi

subroutine f90wrap_alloc_arrays__array__psi(this, nd, dtype, dshape, dloc)
    use datatypes_allocatable, only: alloc_arrays
    implicit none
    type alloc_arrays_ptr_type
        type(alloc_arrays), pointer :: p => NULL()
    end type alloc_arrays_ptr_type
    integer, intent(in) :: this(2)
    type(alloc_arrays_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%psi)) then
        dshape(1:2) = shape(this_ptr%p%psi)
        dloc = loc(this_ptr%p%psi)
    else
        dloc = 0
    end if
end subroutine f90wrap_alloc_arrays__array__psi

subroutine f90wrap_alloc_arrays__array__chi_shape(this, nd, dtype, dshape, dloc)
    use datatypes_allocatable, only: alloc_arrays
    implicit none
    type alloc_arrays_ptr_type
        type(alloc_arrays), pointer :: p => NULL()
    end type alloc_arrays_ptr_type
    integer, intent(in) :: this(2)
    type(alloc_arrays_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%chi_shape)
    dloc = loc(this_ptr%p%chi_shape)
end subroutine f90wrap_alloc_arrays__array__chi_shape

subroutine f90wrap_alloc_arrays__array__psi_shape(this, nd, dtype, dshape, dloc)
    use datatypes_allocatable, only: alloc_arrays
    implicit none
    type alloc_arrays_ptr_type
        type(alloc_arrays), pointer :: p => NULL()
    end type alloc_arrays_ptr_type
    integer, intent(in) :: this(2)
    type(alloc_arrays_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%psi_shape)
    dloc = loc(this_ptr%p%psi_shape)
end subroutine f90wrap_alloc_arrays__array__psi_shape

subroutine f90wrap_alloc_arrays_initialise(this)
    use datatypes_allocatable, only: alloc_arrays
    implicit none
    
    type alloc_arrays_ptr_type
        type(alloc_arrays), pointer :: p => NULL()
    end type alloc_arrays_ptr_type
    type(alloc_arrays_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_alloc_arrays_initialise

subroutine f90wrap_alloc_arrays_finalise(this)
    use datatypes_allocatable, only: alloc_arrays
    implicit none
    
    type alloc_arrays_ptr_type
        type(alloc_arrays), pointer :: p => NULL()
    end type alloc_arrays_ptr_type
    type(alloc_arrays_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_alloc_arrays_finalise

subroutine f90wrap_init_alloc_arrays(dertype, m, n)
    use datatypes_allocatable, only: alloc_arrays, init_alloc_arrays
    implicit none
    
    type alloc_arrays_ptr_type
        type(alloc_arrays), pointer :: p => NULL()
    end type alloc_arrays_ptr_type
    type(alloc_arrays_ptr_type) :: dertype_ptr
    integer, intent(in), dimension(2) :: dertype
    integer(4), intent(in) :: m
    integer(4), intent(in) :: n
    dertype_ptr = transfer(dertype, dertype_ptr)
    call init_alloc_arrays(dertype=dertype_ptr%p, m=m, n=n)
end subroutine f90wrap_init_alloc_arrays

subroutine f90wrap_destroy_alloc_arrays(dertype)
    use datatypes_allocatable, only: alloc_arrays, destroy_alloc_arrays
    implicit none
    
    type alloc_arrays_ptr_type
        type(alloc_arrays), pointer :: p => NULL()
    end type alloc_arrays_ptr_type
    type(alloc_arrays_ptr_type) :: dertype_ptr
    integer, intent(in), dimension(2) :: dertype
    dertype_ptr = transfer(dertype, dertype_ptr)
    call destroy_alloc_arrays(dertype=dertype_ptr%p)
end subroutine f90wrap_destroy_alloc_arrays

! End of module datatypes_allocatable defined in file datatypes.fpp

! Module datatypes defined in file datatypes.fpp

subroutine f90wrap_different_types__get__alpha(this, f90wrap_alpha)
    use datatypes, only: different_types
    implicit none
    type different_types_ptr_type
        type(different_types), pointer :: p => NULL()
    end type different_types_ptr_type
    integer, intent(in)   :: this(2)
    type(different_types_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_alpha
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_alpha = this_ptr%p%alpha
end subroutine f90wrap_different_types__get__alpha

subroutine f90wrap_different_types__set__alpha(this, f90wrap_alpha)
    use datatypes, only: different_types
    implicit none
    type different_types_ptr_type
        type(different_types), pointer :: p => NULL()
    end type different_types_ptr_type
    integer, intent(in)   :: this(2)
    type(different_types_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_alpha
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%alpha = f90wrap_alpha
end subroutine f90wrap_different_types__set__alpha

subroutine f90wrap_different_types__get__beta(this, f90wrap_beta)
    use datatypes, only: different_types
    implicit none
    type different_types_ptr_type
        type(different_types), pointer :: p => NULL()
    end type different_types_ptr_type
    integer, intent(in)   :: this(2)
    type(different_types_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_beta
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_beta = this_ptr%p%beta
end subroutine f90wrap_different_types__get__beta

subroutine f90wrap_different_types__set__beta(this, f90wrap_beta)
    use datatypes, only: different_types
    implicit none
    type different_types_ptr_type
        type(different_types), pointer :: p => NULL()
    end type different_types_ptr_type
    integer, intent(in)   :: this(2)
    type(different_types_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_beta
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%beta = f90wrap_beta
end subroutine f90wrap_different_types__set__beta

subroutine f90wrap_different_types__get__delta(this, f90wrap_delta)
    use datatypes, only: different_types
    implicit none
    type different_types_ptr_type
        type(different_types), pointer :: p => NULL()
    end type different_types_ptr_type
    integer, intent(in)   :: this(2)
    type(different_types_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_delta
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_delta = this_ptr%p%delta
end subroutine f90wrap_different_types__get__delta

subroutine f90wrap_different_types__set__delta(this, f90wrap_delta)
    use datatypes, only: different_types
    implicit none
    type different_types_ptr_type
        type(different_types), pointer :: p => NULL()
    end type different_types_ptr_type
    integer, intent(in)   :: this(2)
    type(different_types_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_delta
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%delta = f90wrap_delta
end subroutine f90wrap_different_types__set__delta

subroutine f90wrap_different_types_initialise(this)
    use datatypes, only: different_types
    implicit none
    
    type different_types_ptr_type
        type(different_types), pointer :: p => NULL()
    end type different_types_ptr_type
    type(different_types_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_different_types_initialise

subroutine f90wrap_different_types_finalise(this)
    use datatypes, only: different_types
    implicit none
    
    type different_types_ptr_type
        type(different_types), pointer :: p => NULL()
    end type different_types_ptr_type
    type(different_types_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_different_types_finalise

subroutine f90wrap_fixed_shape_arrays__array__eta(this, nd, dtype, dshape, dloc)
    use datatypes, only: fixed_shape_arrays
    implicit none
    type fixed_shape_arrays_ptr_type
        type(fixed_shape_arrays), pointer :: p => NULL()
    end type fixed_shape_arrays_ptr_type
    integer, intent(in) :: this(2)
    type(fixed_shape_arrays_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:2) = shape(this_ptr%p%eta)
    dloc = loc(this_ptr%p%eta)
end subroutine f90wrap_fixed_shape_arrays__array__eta

subroutine f90wrap_fixed_shape_arrays__array__theta(this, nd, dtype, dshape, dloc)
    use datatypes, only: fixed_shape_arrays
    implicit none
    type fixed_shape_arrays_ptr_type
        type(fixed_shape_arrays), pointer :: p => NULL()
    end type fixed_shape_arrays_ptr_type
    integer, intent(in) :: this(2)
    type(fixed_shape_arrays_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    dshape(1:2) = shape(this_ptr%p%theta)
    dloc = loc(this_ptr%p%theta)
end subroutine f90wrap_fixed_shape_arrays__array__theta

subroutine f90wrap_fixed_shape_arrays__array__iota(this, nd, dtype, dshape, dloc)
    use datatypes, only: fixed_shape_arrays
    implicit none
    type fixed_shape_arrays_ptr_type
        type(fixed_shape_arrays), pointer :: p => NULL()
    end type fixed_shape_arrays_ptr_type
    integer, intent(in) :: this(2)
    type(fixed_shape_arrays_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    dshape(1:2) = shape(this_ptr%p%iota)
    dloc = loc(this_ptr%p%iota)
end subroutine f90wrap_fixed_shape_arrays__array__iota

subroutine f90wrap_fixed_shape_arrays_initialise(this)
    use datatypes, only: fixed_shape_arrays
    implicit none
    
    type fixed_shape_arrays_ptr_type
        type(fixed_shape_arrays), pointer :: p => NULL()
    end type fixed_shape_arrays_ptr_type
    type(fixed_shape_arrays_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_fixed_shape_arrays_initialise

subroutine f90wrap_fixed_shape_arrays_finalise(this)
    use datatypes, only: fixed_shape_arrays
    implicit none
    
    type fixed_shape_arrays_ptr_type
        type(fixed_shape_arrays), pointer :: p => NULL()
    end type fixed_shape_arrays_ptr_type
    type(fixed_shape_arrays_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_fixed_shape_arrays_finalise

subroutine f90wrap_nested__get__mu(this, f90wrap_mu)
    use datatypes, only: different_types, nested
    implicit none
    type nested_ptr_type
        type(nested), pointer :: p => NULL()
    end type nested_ptr_type
    type different_types_ptr_type
        type(different_types), pointer :: p => NULL()
    end type different_types_ptr_type
    integer, intent(in)   :: this(2)
    type(nested_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_mu(2)
    type(different_types_ptr_type) :: mu_ptr
    
    this_ptr = transfer(this, this_ptr)
    mu_ptr%p => this_ptr%p%mu
    f90wrap_mu = transfer(mu_ptr,f90wrap_mu)
end subroutine f90wrap_nested__get__mu

subroutine f90wrap_nested__set__mu(this, f90wrap_mu)
    use datatypes, only: different_types, nested
    implicit none
    type nested_ptr_type
        type(nested), pointer :: p => NULL()
    end type nested_ptr_type
    type different_types_ptr_type
        type(different_types), pointer :: p => NULL()
    end type different_types_ptr_type
    integer, intent(in)   :: this(2)
    type(nested_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_mu(2)
    type(different_types_ptr_type) :: mu_ptr
    
    this_ptr = transfer(this, this_ptr)
    mu_ptr = transfer(f90wrap_mu,mu_ptr)
    this_ptr%p%mu = mu_ptr%p
end subroutine f90wrap_nested__set__mu

subroutine f90wrap_nested__get__nu(this, f90wrap_nu)
    use datatypes, only: fixed_shape_arrays, nested
    implicit none
    type nested_ptr_type
        type(nested), pointer :: p => NULL()
    end type nested_ptr_type
    type fixed_shape_arrays_ptr_type
        type(fixed_shape_arrays), pointer :: p => NULL()
    end type fixed_shape_arrays_ptr_type
    integer, intent(in)   :: this(2)
    type(nested_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_nu(2)
    type(fixed_shape_arrays_ptr_type) :: nu_ptr
    
    this_ptr = transfer(this, this_ptr)
    nu_ptr%p => this_ptr%p%nu
    f90wrap_nu = transfer(nu_ptr,f90wrap_nu)
end subroutine f90wrap_nested__get__nu

subroutine f90wrap_nested__set__nu(this, f90wrap_nu)
    use datatypes, only: fixed_shape_arrays, nested
    implicit none
    type nested_ptr_type
        type(nested), pointer :: p => NULL()
    end type nested_ptr_type
    type fixed_shape_arrays_ptr_type
        type(fixed_shape_arrays), pointer :: p => NULL()
    end type fixed_shape_arrays_ptr_type
    integer, intent(in)   :: this(2)
    type(nested_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_nu(2)
    type(fixed_shape_arrays_ptr_type) :: nu_ptr
    
    this_ptr = transfer(this, this_ptr)
    nu_ptr = transfer(f90wrap_nu,nu_ptr)
    this_ptr%p%nu = nu_ptr%p
end subroutine f90wrap_nested__set__nu

subroutine f90wrap_nested_initialise(this)
    use datatypes, only: nested
    implicit none
    
    type nested_ptr_type
        type(nested), pointer :: p => NULL()
    end type nested_ptr_type
    type(nested_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_nested_initialise

subroutine f90wrap_nested_finalise(this)
    use datatypes, only: nested
    implicit none
    
    type nested_ptr_type
        type(nested), pointer :: p => NULL()
    end type nested_ptr_type
    type(nested_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_nested_finalise

subroutine f90wrap_pointer_arrays__array__chi(this, nd, dtype, dshape, dloc)
    use datatypes, only: pointer_arrays
    implicit none
    type pointer_arrays_ptr_type
        type(pointer_arrays), pointer :: p => NULL()
    end type pointer_arrays_ptr_type
    integer, intent(in) :: this(2)
    type(pointer_arrays_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    dshape(1:2) = shape(this_ptr%p%chi)
    dloc = loc(this_ptr%p%chi)
end subroutine f90wrap_pointer_arrays__array__chi

subroutine f90wrap_pointer_arrays__array__psi(this, nd, dtype, dshape, dloc)
    use datatypes, only: pointer_arrays
    implicit none
    type pointer_arrays_ptr_type
        type(pointer_arrays), pointer :: p => NULL()
    end type pointer_arrays_ptr_type
    integer, intent(in) :: this(2)
    type(pointer_arrays_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    dshape(1:2) = shape(this_ptr%p%psi)
    dloc = loc(this_ptr%p%psi)
end subroutine f90wrap_pointer_arrays__array__psi

subroutine f90wrap_pointer_arrays__array__chi_shape(this, nd, dtype, dshape, dloc)
    use datatypes, only: pointer_arrays
    implicit none
    type pointer_arrays_ptr_type
        type(pointer_arrays), pointer :: p => NULL()
    end type pointer_arrays_ptr_type
    integer, intent(in) :: this(2)
    type(pointer_arrays_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%chi_shape)
    dloc = loc(this_ptr%p%chi_shape)
end subroutine f90wrap_pointer_arrays__array__chi_shape

subroutine f90wrap_pointer_arrays__array__psi_shape(this, nd, dtype, dshape, dloc)
    use datatypes, only: pointer_arrays
    implicit none
    type pointer_arrays_ptr_type
        type(pointer_arrays), pointer :: p => NULL()
    end type pointer_arrays_ptr_type
    integer, intent(in) :: this(2)
    type(pointer_arrays_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%psi_shape)
    dloc = loc(this_ptr%p%psi_shape)
end subroutine f90wrap_pointer_arrays__array__psi_shape

subroutine f90wrap_pointer_arrays_initialise(this)
    use datatypes, only: pointer_arrays
    implicit none
    
    type pointer_arrays_ptr_type
        type(pointer_arrays), pointer :: p => NULL()
    end type pointer_arrays_ptr_type
    type(pointer_arrays_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_pointer_arrays_initialise

subroutine f90wrap_pointer_arrays_finalise(this)
    use datatypes, only: pointer_arrays
    implicit none
    
    type pointer_arrays_ptr_type
        type(pointer_arrays), pointer :: p => NULL()
    end type pointer_arrays_ptr_type
    type(pointer_arrays_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_pointer_arrays_finalise

subroutine f90wrap_alloc_arrays_2__array__chi(this, nd, dtype, dshape, dloc)
    use datatypes, only: alloc_arrays_2
    implicit none
    type alloc_arrays_2_ptr_type
        type(alloc_arrays_2), pointer :: p => NULL()
    end type alloc_arrays_2_ptr_type
    integer, intent(in) :: this(2)
    type(alloc_arrays_2_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%chi)) then
        dshape(1:2) = shape(this_ptr%p%chi)
        dloc = loc(this_ptr%p%chi)
    else
        dloc = 0
    end if
end subroutine f90wrap_alloc_arrays_2__array__chi

subroutine f90wrap_alloc_arrays_2__array__psi(this, nd, dtype, dshape, dloc)
    use datatypes, only: alloc_arrays_2
    implicit none
    type alloc_arrays_2_ptr_type
        type(alloc_arrays_2), pointer :: p => NULL()
    end type alloc_arrays_2_ptr_type
    integer, intent(in) :: this(2)
    type(alloc_arrays_2_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%psi)) then
        dshape(1:2) = shape(this_ptr%p%psi)
        dloc = loc(this_ptr%p%psi)
    else
        dloc = 0
    end if
end subroutine f90wrap_alloc_arrays_2__array__psi

subroutine f90wrap_alloc_arrays_2__array__chi_shape(this, nd, dtype, dshape, dloc)
    use datatypes, only: alloc_arrays_2
    implicit none
    type alloc_arrays_2_ptr_type
        type(alloc_arrays_2), pointer :: p => NULL()
    end type alloc_arrays_2_ptr_type
    integer, intent(in) :: this(2)
    type(alloc_arrays_2_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%chi_shape)
    dloc = loc(this_ptr%p%chi_shape)
end subroutine f90wrap_alloc_arrays_2__array__chi_shape

subroutine f90wrap_alloc_arrays_2__array__psi_shape(this, nd, dtype, dshape, dloc)
    use datatypes, only: alloc_arrays_2
    implicit none
    type alloc_arrays_2_ptr_type
        type(alloc_arrays_2), pointer :: p => NULL()
    end type alloc_arrays_2_ptr_type
    integer, intent(in) :: this(2)
    type(alloc_arrays_2_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%psi_shape)
    dloc = loc(this_ptr%p%psi_shape)
end subroutine f90wrap_alloc_arrays_2__array__psi_shape

subroutine f90wrap_alloc_arrays_2_initialise(this)
    use datatypes, only: alloc_arrays_2
    implicit none
    
    type alloc_arrays_2_ptr_type
        type(alloc_arrays_2), pointer :: p => NULL()
    end type alloc_arrays_2_ptr_type
    type(alloc_arrays_2_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_alloc_arrays_2_initialise

subroutine f90wrap_alloc_arrays_2_finalise(this)
    use datatypes, only: alloc_arrays_2
    implicit none
    
    type alloc_arrays_2_ptr_type
        type(alloc_arrays_2), pointer :: p => NULL()
    end type alloc_arrays_2_ptr_type
    type(alloc_arrays_2_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_alloc_arrays_2_finalise

subroutine f90wrap_array_nested__array_getitem__xi(f90wrap_this, f90wrap_i, xiitem)
    
    use datatypes, only: array_nested, different_types
    implicit none
    
    type array_nested_ptr_type
        type(array_nested), pointer :: p => NULL()
    end type array_nested_ptr_type
    type different_types_ptr_type
        type(different_types), pointer :: p => NULL()
    end type different_types_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(array_nested_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: xiitem(2)
    type(different_types_ptr_type) :: xi_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%xi)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%xi)) then
            call f90wrap_abort("array index out of range")
        else
            xi_ptr%p => this_ptr%p%xi(f90wrap_i)
            xiitem = transfer(xi_ptr,xiitem)
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_array_nested__array_getitem__xi

subroutine f90wrap_array_nested__array_setitem__xi(f90wrap_this, f90wrap_i, xiitem)
    
    use datatypes, only: array_nested, different_types
    implicit none
    
    type array_nested_ptr_type
        type(array_nested), pointer :: p => NULL()
    end type array_nested_ptr_type
    type different_types_ptr_type
        type(different_types), pointer :: p => NULL()
    end type different_types_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(array_nested_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: xiitem(2)
    type(different_types_ptr_type) :: xi_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%xi)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%xi)) then
            call f90wrap_abort("array index out of range")
        else
            xi_ptr = transfer(xiitem,xi_ptr)
            this_ptr%p%xi(f90wrap_i) = xi_ptr%p
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_array_nested__array_setitem__xi

subroutine f90wrap_array_nested__array_len__xi(f90wrap_this, f90wrap_n)
    
    use datatypes, only: array_nested, different_types
    implicit none
    
    type array_nested_ptr_type
        type(array_nested), pointer :: p => NULL()
    end type array_nested_ptr_type
    type different_types_ptr_type
        type(different_types), pointer :: p => NULL()
    end type different_types_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(array_nested_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%xi)) then
        f90wrap_n = size(this_ptr%p%xi)
    else
        f90wrap_n = 0
    end if
end subroutine f90wrap_array_nested__array_len__xi

subroutine f90wrap_array_nested__array_getitem__omicron(f90wrap_this, f90wrap_i, omicronitem)
    
    use datatypes, only: array_nested, fixed_shape_arrays
    implicit none
    
    type array_nested_ptr_type
        type(array_nested), pointer :: p => NULL()
    end type array_nested_ptr_type
    type fixed_shape_arrays_ptr_type
        type(fixed_shape_arrays), pointer :: p => NULL()
    end type fixed_shape_arrays_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(array_nested_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: omicronitem(2)
    type(fixed_shape_arrays_ptr_type) :: omicron_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%omicron)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%omicron)) then
            call f90wrap_abort("array index out of range")
        else
            omicron_ptr%p => this_ptr%p%omicron(f90wrap_i)
            omicronitem = transfer(omicron_ptr,omicronitem)
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_array_nested__array_getitem__omicron

subroutine f90wrap_array_nested__array_setitem__omicron(f90wrap_this, f90wrap_i, omicronitem)
    
    use datatypes, only: array_nested, fixed_shape_arrays
    implicit none
    
    type array_nested_ptr_type
        type(array_nested), pointer :: p => NULL()
    end type array_nested_ptr_type
    type fixed_shape_arrays_ptr_type
        type(fixed_shape_arrays), pointer :: p => NULL()
    end type fixed_shape_arrays_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(array_nested_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: omicronitem(2)
    type(fixed_shape_arrays_ptr_type) :: omicron_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%omicron)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%omicron)) then
            call f90wrap_abort("array index out of range")
        else
            omicron_ptr = transfer(omicronitem,omicron_ptr)
            this_ptr%p%omicron(f90wrap_i) = omicron_ptr%p
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_array_nested__array_setitem__omicron

subroutine f90wrap_array_nested__array_len__omicron(f90wrap_this, f90wrap_n)
    
    use datatypes, only: array_nested, fixed_shape_arrays
    implicit none
    
    type array_nested_ptr_type
        type(array_nested), pointer :: p => NULL()
    end type array_nested_ptr_type
    type fixed_shape_arrays_ptr_type
        type(fixed_shape_arrays), pointer :: p => NULL()
    end type fixed_shape_arrays_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(array_nested_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%omicron)) then
        f90wrap_n = size(this_ptr%p%omicron)
    else
        f90wrap_n = 0
    end if
end subroutine f90wrap_array_nested__array_len__omicron

subroutine f90wrap_array_nested__array_getitem__pi(f90wrap_this, f90wrap_i, piitem)
    
    use datatypes, only: array_nested
    use datatypes_allocatable, only: alloc_arrays
    implicit none
    
    type array_nested_ptr_type
        type(array_nested), pointer :: p => NULL()
    end type array_nested_ptr_type
    type alloc_arrays_ptr_type
        type(alloc_arrays), pointer :: p => NULL()
    end type alloc_arrays_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(array_nested_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: piitem(2)
    type(alloc_arrays_ptr_type) :: pi_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%pi)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%pi)) then
            call f90wrap_abort("array index out of range")
        else
            pi_ptr%p => this_ptr%p%pi(f90wrap_i)
            piitem = transfer(pi_ptr,piitem)
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_array_nested__array_getitem__pi

subroutine f90wrap_array_nested__array_setitem__pi(f90wrap_this, f90wrap_i, piitem)
    
    use datatypes, only: array_nested
    use datatypes_allocatable, only: alloc_arrays
    implicit none
    
    type array_nested_ptr_type
        type(array_nested), pointer :: p => NULL()
    end type array_nested_ptr_type
    type alloc_arrays_ptr_type
        type(alloc_arrays), pointer :: p => NULL()
    end type alloc_arrays_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(array_nested_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: piitem(2)
    type(alloc_arrays_ptr_type) :: pi_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%pi)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%pi)) then
            call f90wrap_abort("array index out of range")
        else
            pi_ptr = transfer(piitem,pi_ptr)
            this_ptr%p%pi(f90wrap_i) = pi_ptr%p
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_array_nested__array_setitem__pi

subroutine f90wrap_array_nested__array_len__pi(f90wrap_this, f90wrap_n)
    
    use datatypes, only: array_nested
    use datatypes_allocatable, only: alloc_arrays
    implicit none
    
    type array_nested_ptr_type
        type(array_nested), pointer :: p => NULL()
    end type array_nested_ptr_type
    type alloc_arrays_ptr_type
        type(alloc_arrays), pointer :: p => NULL()
    end type alloc_arrays_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(array_nested_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%pi)) then
        f90wrap_n = size(this_ptr%p%pi)
    else
        f90wrap_n = 0
    end if
end subroutine f90wrap_array_nested__array_len__pi

subroutine f90wrap_array_nested_initialise(this)
    use datatypes, only: array_nested
    implicit none
    
    type array_nested_ptr_type
        type(array_nested), pointer :: p => NULL()
    end type array_nested_ptr_type
    type(array_nested_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_array_nested_initialise

subroutine f90wrap_array_nested_finalise(this)
    use datatypes, only: array_nested
    implicit none
    
    type array_nested_ptr_type
        type(array_nested), pointer :: p => NULL()
    end type array_nested_ptr_type
    type(array_nested_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_array_nested_finalise

subroutine f90wrap_init_array_nested(dertype, size_bn)
    use datatypes, only: array_nested, init_array_nested
    implicit none
    
    type array_nested_ptr_type
        type(array_nested), pointer :: p => NULL()
    end type array_nested_ptr_type
    type(array_nested_ptr_type) :: dertype_ptr
    integer, intent(in), dimension(2) :: dertype
    integer(4), intent(in) :: size_bn
    dertype_ptr = transfer(dertype, dertype_ptr)
    call init_array_nested(dertype=dertype_ptr%p, size=size_bn)
end subroutine f90wrap_init_array_nested

subroutine f90wrap_destroy_array_nested(dertype)
    use datatypes, only: array_nested, destroy_array_nested
    implicit none
    
    type array_nested_ptr_type
        type(array_nested), pointer :: p => NULL()
    end type array_nested_ptr_type
    type(array_nested_ptr_type) :: dertype_ptr
    integer, intent(in), dimension(2) :: dertype
    dertype_ptr = transfer(dertype, dertype_ptr)
    call destroy_array_nested(dertype=dertype_ptr%p)
end subroutine f90wrap_destroy_array_nested

! End of module datatypes defined in file datatypes.fpp


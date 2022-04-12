! Module grid_module defined in file Grid_module.fpp

subroutine f90wrap_grid_type__array__rhoS(this, nd, dtype, dshape, dloc)
    use grid_module, only: grid_type
    implicit none
    type grid_type_ptr_type
        type(grid_type), pointer :: p => NULL()
    end type grid_type_ptr_type
    integer, intent(in) :: this(2)
    type(grid_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 4
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%rhoS)) then
        dshape(1:4) = shape(this_ptr%p%rhoS)
        dloc = loc(this_ptr%p%rhoS)
    else
        dloc = 0
    end if
end subroutine f90wrap_grid_type__array__rhoS

subroutine f90wrap_grid_type__array__vlpp(this, nd, dtype, dshape, dloc)
    use grid_module, only: grid_type
    implicit none
    type grid_type_ptr_type
        type(grid_type), pointer :: p => NULL()
    end type grid_type_ptr_type
    integer, intent(in) :: this(2)
    type(grid_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%vlpp)) then
        dshape(1:3) = shape(this_ptr%p%vlpp)
        dloc = loc(this_ptr%p%vlpp)
    else
        dloc = 0
    end if
end subroutine f90wrap_grid_type__array__vlpp

subroutine f90wrap_grid_type__array__eval(this, nd, dtype, dshape, dloc)
    use grid_module, only: grid_type
    implicit none
    type grid_type_ptr_type
        type(grid_type), pointer :: p => NULL()
    end type grid_type_ptr_type
    integer, intent(in) :: this(2)
    type(grid_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%eval)) then
        dshape(1:3) = shape(this_ptr%p%eval)
        dloc = loc(this_ptr%p%eval)
    else
        dloc = 0
    end if
end subroutine f90wrap_grid_type__array__eval

subroutine f90wrap_grid_type__array__gVec(this, nd, dtype, dshape, dloc)
    use grid_module, only: grid_type
    implicit none
    type grid_type_ptr_type
        type(grid_type), pointer :: p => NULL()
    end type grid_type_ptr_type
    integer, intent(in) :: this(2)
    type(grid_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%gVec)) then
        dshape(1:2) = shape(this_ptr%p%gVec)
        dloc = loc(this_ptr%p%gVec)
    else
        dloc = 0
    end if
end subroutine f90wrap_grid_type__array__gVec

subroutine f90wrap_grid_type__array__gMask(this, nd, dtype, dshape, dloc)
    use grid_module, only: grid_type
    implicit none
    type grid_type_ptr_type
        type(grid_type), pointer :: p => NULL()
    end type grid_type_ptr_type
    integer, intent(in) :: this(2)
    type(grid_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%gMask)) then
        dshape(1:1) = shape(this_ptr%p%gMask)
        dloc = loc(this_ptr%p%gMask)
    else
        dloc = 0
    end if
end subroutine f90wrap_grid_type__array__gMask

subroutine f90wrap_grid_type__array__rVec(this, nd, dtype, dshape, dloc)
    use grid_module, only: grid_type
    implicit none
    type grid_type_ptr_type
        type(grid_type), pointer :: p => NULL()
    end type grid_type_ptr_type
    integer, intent(in) :: this(2)
    type(grid_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%rVec)) then
        dshape(1:2) = shape(this_ptr%p%rVec)
        dloc = loc(this_ptr%p%rVec)
    else
        dloc = 0
    end if
end subroutine f90wrap_grid_type__array__rVec

subroutine f90wrap_grid_type_initialise(this)
    use grid_module, only: grid_type
    implicit none
    
    type grid_type_ptr_type
        type(grid_type), pointer :: p => NULL()
    end type grid_type_ptr_type
    type(grid_type_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_grid_type_initialise

subroutine f90wrap_grid_type_finalise(this)
    use grid_module, only: grid_type
    implicit none
    
    type grid_type_ptr_type
        type(grid_type), pointer :: p => NULL()
    end type grid_type_ptr_type
    type(grid_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_grid_type_finalise

subroutine f90wrap_kgrid_type__array__vec(this, nd, dtype, dshape, dloc)
    use grid_module, only: kgrid_type
    implicit none
    type kgrid_type_ptr_type
        type(kgrid_type), pointer :: p => NULL()
    end type kgrid_type_ptr_type
    integer, intent(in) :: this(2)
    type(kgrid_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%vec)) then
        dshape(1:2) = shape(this_ptr%p%vec)
        dloc = loc(this_ptr%p%vec)
    else
        dloc = 0
    end if
end subroutine f90wrap_kgrid_type__array__vec

subroutine f90wrap_kgrid_type__array__vcar(this, nd, dtype, dshape, dloc)
    use grid_module, only: kgrid_type
    implicit none
    type kgrid_type_ptr_type
        type(kgrid_type), pointer :: p => NULL()
    end type kgrid_type_ptr_type
    integer, intent(in) :: this(2)
    type(kgrid_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%vcar)) then
        dshape(1:2) = shape(this_ptr%p%vcar)
        dloc = loc(this_ptr%p%vcar)
    else
        dloc = 0
    end if
end subroutine f90wrap_kgrid_type__array__vcar

subroutine f90wrap_kgrid_type__array__wk(this, nd, dtype, dshape, dloc)
    use grid_module, only: kgrid_type
    implicit none
    type kgrid_type_ptr_type
        type(kgrid_type), pointer :: p => NULL()
    end type kgrid_type_ptr_type
    integer, intent(in) :: this(2)
    type(kgrid_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%wk)) then
        dshape(1:1) = shape(this_ptr%p%wk)
        dloc = loc(this_ptr%p%wk)
    else
        dloc = 0
    end if
end subroutine f90wrap_kgrid_type__array__wk

subroutine f90wrap_kgrid_type_initialise(this)
    use grid_module, only: kgrid_type
    implicit none
    
    type kgrid_type_ptr_type
        type(kgrid_type), pointer :: p => NULL()
    end type kgrid_type_ptr_type
    type(kgrid_type_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_kgrid_type_initialise

subroutine f90wrap_kgrid_type_finalise(this)
    use grid_module, only: kgrid_type
    implicit none
    
    type kgrid_type_ptr_type
        type(kgrid_type), pointer :: p => NULL()
    end type kgrid_type_ptr_type
    type(kgrid_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_kgrid_type_finalise

subroutine f90wrap_eigen_type__array__wvf(this, nd, dtype, dshape, dloc)
    use grid_module, only: eigen_type
    implicit none
    type eigen_type_ptr_type
        type(eigen_type), pointer :: p => NULL()
    end type eigen_type_ptr_type
    integer, intent(in) :: this(2)
    type(eigen_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 4
    dtype = 15
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%wvf)) then
        dshape(1:4) = shape(this_ptr%p%wvf)
        dloc = loc(this_ptr%p%wvf)
    else
        dloc = 0
    end if
end subroutine f90wrap_eigen_type__array__wvf

subroutine f90wrap_eigen_type__array__val(this, nd, dtype, dshape, dloc)
    use grid_module, only: eigen_type
    implicit none
    type eigen_type_ptr_type
        type(eigen_type), pointer :: p => NULL()
    end type eigen_type_ptr_type
    integer, intent(in) :: this(2)
    type(eigen_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%val)) then
        dshape(1:3) = shape(this_ptr%p%val)
        dloc = loc(this_ptr%p%val)
    else
        dloc = 0
    end if
end subroutine f90wrap_eigen_type__array__val

subroutine f90wrap_eigen_type_initialise(this)
    use grid_module, only: eigen_type
    implicit none
    
    type eigen_type_ptr_type
        type(eigen_type), pointer :: p => NULL()
    end type eigen_type_ptr_type
    type(eigen_type_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_eigen_type_initialise

subroutine f90wrap_eigen_type_finalise(this)
    use grid_module, only: eigen_type
    implicit none
    
    type eigen_type_ptr_type
        type(eigen_type), pointer :: p => NULL()
    end type eigen_type_ptr_type
    type(eigen_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_eigen_type_finalise

subroutine f90wrap_eigen_type_r__array__wvf(this, nd, dtype, dshape, dloc)
    use grid_module, only: eigen_type_r
    implicit none
    type eigen_type_r_ptr_type
        type(eigen_type_r), pointer :: p => NULL()
    end type eigen_type_r_ptr_type
    integer, intent(in) :: this(2)
    type(eigen_type_r_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 4
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%wvf)) then
        dshape(1:4) = shape(this_ptr%p%wvf)
        dloc = loc(this_ptr%p%wvf)
    else
        dloc = 0
    end if
end subroutine f90wrap_eigen_type_r__array__wvf

subroutine f90wrap_eigen_type_r__array__val(this, nd, dtype, dshape, dloc)
    use grid_module, only: eigen_type_r
    implicit none
    type eigen_type_r_ptr_type
        type(eigen_type_r), pointer :: p => NULL()
    end type eigen_type_r_ptr_type
    integer, intent(in) :: this(2)
    type(eigen_type_r_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%val)) then
        dshape(1:3) = shape(this_ptr%p%val)
        dloc = loc(this_ptr%p%val)
    else
        dloc = 0
    end if
end subroutine f90wrap_eigen_type_r__array__val

subroutine f90wrap_eigen_type_r_initialise(this)
    use grid_module, only: eigen_type_r
    implicit none
    
    type eigen_type_r_ptr_type
        type(eigen_type_r), pointer :: p => NULL()
    end type eigen_type_r_ptr_type
    type(eigen_type_r_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_eigen_type_r_initialise

subroutine f90wrap_eigen_type_r_finalise(this)
    use grid_module, only: eigen_type_r
    implicit none
    
    type eigen_type_r_ptr_type
        type(eigen_type_r), pointer :: p => NULL()
    end type eigen_type_r_ptr_type
    type(eigen_type_r_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_eigen_type_r_finalise

subroutine f90wrap_charge_sphere__get__OneDLength(this, f90wrap_OneDLength)
    use grid_module, only: charge_sphere
    implicit none
    type charge_sphere_ptr_type
        type(charge_sphere), pointer :: p => NULL()
    end type charge_sphere_ptr_type
    integer, intent(in)   :: this(2)
    type(charge_sphere_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_OneDLength
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_OneDLength = this_ptr%p%OneDLength
end subroutine f90wrap_charge_sphere__get__OneDLength

subroutine f90wrap_charge_sphere__set__OneDLength(this, f90wrap_OneDLength)
    use grid_module, only: charge_sphere
    implicit none
    type charge_sphere_ptr_type
        type(charge_sphere), pointer :: p => NULL()
    end type charge_sphere_ptr_type
    integer, intent(in)   :: this(2)
    type(charge_sphere_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_OneDLength
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%OneDLength = f90wrap_OneDLength
end subroutine f90wrap_charge_sphere__set__OneDLength

subroutine f90wrap_charge_sphere__get__volume(this, f90wrap_volume)
    use grid_module, only: charge_sphere
    implicit none
    type charge_sphere_ptr_type
        type(charge_sphere), pointer :: p => NULL()
    end type charge_sphere_ptr_type
    integer, intent(in)   :: this(2)
    type(charge_sphere_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_volume
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_volume = this_ptr%p%volume
end subroutine f90wrap_charge_sphere__get__volume

subroutine f90wrap_charge_sphere__set__volume(this, f90wrap_volume)
    use grid_module, only: charge_sphere
    implicit none
    type charge_sphere_ptr_type
        type(charge_sphere), pointer :: p => NULL()
    end type charge_sphere_ptr_type
    integer, intent(in)   :: this(2)
    type(charge_sphere_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_volume
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%volume = f90wrap_volume
end subroutine f90wrap_charge_sphere__set__volume

subroutine f90wrap_charge_sphere__array__x(this, nd, dtype, dshape, dloc)
    use grid_module, only: charge_sphere
    implicit none
    type charge_sphere_ptr_type
        type(charge_sphere), pointer :: p => NULL()
    end type charge_sphere_ptr_type
    integer, intent(in) :: this(2)
    type(charge_sphere_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%x)) then
        dshape(1:1) = shape(this_ptr%p%x)
        dloc = loc(this_ptr%p%x)
    else
        dloc = 0
    end if
end subroutine f90wrap_charge_sphere__array__x

subroutine f90wrap_charge_sphere__array__y(this, nd, dtype, dshape, dloc)
    use grid_module, only: charge_sphere
    implicit none
    type charge_sphere_ptr_type
        type(charge_sphere), pointer :: p => NULL()
    end type charge_sphere_ptr_type
    integer, intent(in) :: this(2)
    type(charge_sphere_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%y)) then
        dshape(1:1) = shape(this_ptr%p%y)
        dloc = loc(this_ptr%p%y)
    else
        dloc = 0
    end if
end subroutine f90wrap_charge_sphere__array__y

subroutine f90wrap_charge_sphere__array__z(this, nd, dtype, dshape, dloc)
    use grid_module, only: charge_sphere
    implicit none
    type charge_sphere_ptr_type
        type(charge_sphere), pointer :: p => NULL()
    end type charge_sphere_ptr_type
    integer, intent(in) :: this(2)
    type(charge_sphere_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%z)) then
        dshape(1:1) = shape(this_ptr%p%z)
        dloc = loc(this_ptr%p%z)
    else
        dloc = 0
    end if
end subroutine f90wrap_charge_sphere__array__z

subroutine f90wrap_charge_sphere__array__n(this, nd, dtype, dshape, dloc)
    use grid_module, only: charge_sphere
    implicit none
    type charge_sphere_ptr_type
        type(charge_sphere), pointer :: p => NULL()
    end type charge_sphere_ptr_type
    integer, intent(in) :: this(2)
    type(charge_sphere_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%n)) then
        dshape(1:1) = shape(this_ptr%p%n)
        dloc = loc(this_ptr%p%n)
    else
        dloc = 0
    end if
end subroutine f90wrap_charge_sphere__array__n

subroutine f90wrap_charge_sphere__array__OneDSphere(this, nd, dtype, dshape, dloc)
    use grid_module, only: charge_sphere
    implicit none
    type charge_sphere_ptr_type
        type(charge_sphere), pointer :: p => NULL()
    end type charge_sphere_ptr_type
    integer, intent(in) :: this(2)
    type(charge_sphere_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%OneDSphere)) then
        dshape(1:2) = shape(this_ptr%p%OneDSphere)
        dloc = loc(this_ptr%p%OneDSphere)
    else
        dloc = 0
    end if
end subroutine f90wrap_charge_sphere__array__OneDSphere

subroutine f90wrap_charge_sphere__array__OneDVeff(this, nd, dtype, dshape, dloc)
    use grid_module, only: charge_sphere
    implicit none
    type charge_sphere_ptr_type
        type(charge_sphere), pointer :: p => NULL()
    end type charge_sphere_ptr_type
    integer, intent(in) :: this(2)
    type(charge_sphere_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%OneDVeff)) then
        dshape(1:1) = shape(this_ptr%p%OneDVeff)
        dloc = loc(this_ptr%p%OneDVeff)
    else
        dloc = 0
    end if
end subroutine f90wrap_charge_sphere__array__OneDVeff

subroutine f90wrap_charge_sphere__array__OneDwvf(this, nd, dtype, dshape, dloc)
    use grid_module, only: charge_sphere
    implicit none
    type charge_sphere_ptr_type
        type(charge_sphere), pointer :: p => NULL()
    end type charge_sphere_ptr_type
    integer, intent(in) :: this(2)
    type(charge_sphere_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%OneDwvf)) then
        dshape(1:3) = shape(this_ptr%p%OneDwvf)
        dloc = loc(this_ptr%p%OneDwvf)
    else
        dloc = 0
    end if
end subroutine f90wrap_charge_sphere__array__OneDwvf

subroutine f90wrap_charge_sphere__array__OneDval(this, nd, dtype, dshape, dloc)
    use grid_module, only: charge_sphere
    implicit none
    type charge_sphere_ptr_type
        type(charge_sphere), pointer :: p => NULL()
    end type charge_sphere_ptr_type
    integer, intent(in) :: this(2)
    type(charge_sphere_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%OneDval)) then
        dshape(1:2) = shape(this_ptr%p%OneDval)
        dloc = loc(this_ptr%p%OneDval)
    else
        dloc = 0
    end if
end subroutine f90wrap_charge_sphere__array__OneDval

subroutine f90wrap_charge_sphere__array__initMO(this, nd, dtype, dshape, dloc)
    use grid_module, only: charge_sphere
    implicit none
    type charge_sphere_ptr_type
        type(charge_sphere), pointer :: p => NULL()
    end type charge_sphere_ptr_type
    integer, intent(in) :: this(2)
    type(charge_sphere_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%initMO)) then
        dshape(1:2) = shape(this_ptr%p%initMO)
        dloc = loc(this_ptr%p%initMO)
    else
        dloc = 0
    end if
end subroutine f90wrap_charge_sphere__array__initMO

subroutine f90wrap_charge_sphere__array__vlpp(this, nd, dtype, dshape, dloc)
    use grid_module, only: charge_sphere
    implicit none
    type charge_sphere_ptr_type
        type(charge_sphere), pointer :: p => NULL()
    end type charge_sphere_ptr_type
    integer, intent(in) :: this(2)
    type(charge_sphere_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%vlpp)) then
        dshape(1:1) = shape(this_ptr%p%vlpp)
        dloc = loc(this_ptr%p%vlpp)
    else
        dloc = 0
    end if
end subroutine f90wrap_charge_sphere__array__vlpp

subroutine f90wrap_charge_sphere_initialise(this)
    use grid_module, only: charge_sphere
    implicit none
    
    type charge_sphere_ptr_type
        type(charge_sphere), pointer :: p => NULL()
    end type charge_sphere_ptr_type
    type(charge_sphere_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_charge_sphere_initialise

subroutine f90wrap_charge_sphere_finalise(this)
    use grid_module, only: charge_sphere
    implicit none
    
    type charge_sphere_ptr_type
        type(charge_sphere), pointer :: p => NULL()
    end type charge_sphere_ptr_type
    type(charge_sphere_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_charge_sphere_finalise

subroutine f90wrap_grid_diff_map_type__array__nz_map(this, nd, dtype, dshape, dloc)
    use grid_module, only: grid_diff_map_type
    implicit none
    type grid_diff_map_type_ptr_type
        type(grid_diff_map_type), pointer :: p => NULL()
    end type grid_diff_map_type_ptr_type
    integer, intent(in) :: this(2)
    type(grid_diff_map_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%nz_map)) then
        dshape(1:2) = shape(this_ptr%p%nz_map)
        dloc = loc(this_ptr%p%nz_map)
    else
        dloc = 0
    end if
end subroutine f90wrap_grid_diff_map_type__array__nz_map

subroutine f90wrap_grid_diff_map_type__array__mycomm_cores(this, nd, dtype, dshape, dloc)
    use grid_module, only: grid_diff_map_type
    implicit none
    type grid_diff_map_type_ptr_type
        type(grid_diff_map_type), pointer :: p => NULL()
    end type grid_diff_map_type_ptr_type
    integer, intent(in) :: this(2)
    type(grid_diff_map_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%mycomm_cores)
    dloc = loc(this_ptr%p%mycomm_cores)
end subroutine f90wrap_grid_diff_map_type__array__mycomm_cores

subroutine f90wrap_grid_diff_map_type__array__mycomm_size(this, nd, dtype, dshape, dloc)
    use grid_module, only: grid_diff_map_type
    implicit none
    type grid_diff_map_type_ptr_type
        type(grid_diff_map_type), pointer :: p => NULL()
    end type grid_diff_map_type_ptr_type
    integer, intent(in) :: this(2)
    type(grid_diff_map_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%mycomm_size)) then
        dshape(1:2) = shape(this_ptr%p%mycomm_size)
        dloc = loc(this_ptr%p%mycomm_size)
    else
        dloc = 0
    end if
end subroutine f90wrap_grid_diff_map_type__array__mycomm_size

subroutine f90wrap_grid_diff_map_type__array__mysend_size(this, nd, dtype, dshape, dloc)
    use grid_module, only: grid_diff_map_type
    implicit none
    type grid_diff_map_type_ptr_type
        type(grid_diff_map_type), pointer :: p => NULL()
    end type grid_diff_map_type_ptr_type
    integer, intent(in) :: this(2)
    type(grid_diff_map_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%mysend_size)) then
        dshape(1:2) = shape(this_ptr%p%mysend_size)
        dloc = loc(this_ptr%p%mysend_size)
    else
        dloc = 0
    end if
end subroutine f90wrap_grid_diff_map_type__array__mysend_size

subroutine f90wrap_grid_diff_map_type__array__local_map(this, nd, dtype, dshape, dloc)
    use grid_module, only: grid_diff_map_type
    implicit none
    type grid_diff_map_type_ptr_type
        type(grid_diff_map_type), pointer :: p => NULL()
    end type grid_diff_map_type_ptr_type
    integer, intent(in) :: this(2)
    type(grid_diff_map_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%local_map)) then
        dshape(1:2) = shape(this_ptr%p%local_map)
        dloc = loc(this_ptr%p%local_map)
    else
        dloc = 0
    end if
end subroutine f90wrap_grid_diff_map_type__array__local_map

subroutine f90wrap_grid_diff_map_type__array__boundary(this, nd, dtype, dshape, dloc)
    use grid_module, only: grid_diff_map_type
    implicit none
    type grid_diff_map_type_ptr_type
        type(grid_diff_map_type), pointer :: p => NULL()
    end type grid_diff_map_type_ptr_type
    integer, intent(in) :: this(2)
    type(grid_diff_map_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:2) = shape(this_ptr%p%boundary)
    dloc = loc(this_ptr%p%boundary)
end subroutine f90wrap_grid_diff_map_type__array__boundary

subroutine f90wrap_grid_diff_map_type_initialise(this)
    use grid_module, only: grid_diff_map_type
    implicit none
    
    type grid_diff_map_type_ptr_type
        type(grid_diff_map_type), pointer :: p => NULL()
    end type grid_diff_map_type_ptr_type
    type(grid_diff_map_type_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_grid_diff_map_type_initialise

subroutine f90wrap_grid_diff_map_type_finalise(this)
    use grid_module, only: grid_diff_map_type
    implicit none
    
    type grid_diff_map_type_ptr_type
        type(grid_diff_map_type), pointer :: p => NULL()
    end type grid_diff_map_type_ptr_type
    type(grid_diff_map_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_grid_diff_map_type_finalise

subroutine f90wrap_build_rgrid
    use grid_module, only: build_rgrid
    implicit none
    
    call build_rgrid()
end subroutine f90wrap_build_rgrid

subroutine f90wrap_build_rgrid_iso
    use grid_module, only: build_rgrid_iso
    implicit none
    
    call build_rgrid_iso()
end subroutine f90wrap_build_rgrid_iso

subroutine f90wrap_destroy_rgrid
    use grid_module, only: destroy_rgrid
    implicit none
    
    call destroy_rgrid()
end subroutine f90wrap_destroy_rgrid

subroutine f90wrap_build_kgrid
    use grid_module, only: build_kgrid
    implicit none
    
    call build_kgrid()
end subroutine f90wrap_build_kgrid

subroutine f90wrap_destroy_kpt
    use grid_module, only: destroy_kpt
    implicit none
    
    call destroy_kpt()
end subroutine f90wrap_destroy_kpt

subroutine f90wrap_build_eigen
    use grid_module, only: build_eigen
    implicit none
    
    call build_eigen()
end subroutine f90wrap_build_eigen

subroutine f90wrap_destroy_eigen
    use grid_module, only: destroy_eigen
    implicit none
    
    call destroy_eigen()
end subroutine f90wrap_destroy_eigen

subroutine f90wrap_fillqtable
    use grid_module, only: fillqtable
    implicit none
    
    call fillqtable()
end subroutine f90wrap_fillqtable

subroutine f90wrap_fillrtable
    use grid_module, only: fillrtable
    implicit none
    
    call fillrtable()
end subroutine f90wrap_fillrtable

subroutine f90wrap_fillrtable_iso
    use grid_module, only: fillrtable_iso
    implicit none
    
    call fillrtable_iso()
end subroutine f90wrap_fillrtable_iso

subroutine f90wrap_build_iso_sphere_grid
    use grid_module, only: build_iso_sphere_grid
    implicit none
    
    call build_iso_sphere_grid()
end subroutine f90wrap_build_iso_sphere_grid

subroutine f90wrap_matchposcar_iso
    use grid_module, only: matchposcar_iso
    implicit none
    
    call matchposcar_iso()
end subroutine f90wrap_matchposcar_iso

subroutine f90wrap_set_car2spe(orig, x, y, z, r, cost, sint, cosp, sinp)
    use grid_module, only: set_car2spe
    implicit none
    
    real(8), intent(in), dimension(3) :: orig
    integer, intent(in) :: x
    integer, intent(in) :: y
    integer, intent(in) :: z
    real(8), intent(out) :: r
    real(8), intent(out) :: cost
    real(8), intent(out) :: sint
    real(8), intent(out) :: cosp
    real(8), intent(out) :: sinp
    call set_car2spe(ORIG=orig, x=x, y=y, z=z, r=r, cost=cost, sint=sint, cosp=cosp, sinp=sinp)
end subroutine f90wrap_set_car2spe

subroutine f90wrap_destroy_iso_sphere_grid
    use grid_module, only: destroy_iso_sphere_grid
    implicit none
    
    call destroy_iso_sphere_grid()
end subroutine f90wrap_destroy_iso_sphere_grid

subroutine f90wrap_iso_vsphere(vgrid, vsphere, n0, n1, n2, n3, n4, n5)
    use grid_module, only: iso_vsphere
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2,n3) :: vgrid
    real(8), intent(inout), dimension(n4,n5) :: vsphere
    integer :: n0
    !f2py intent(hide), depend(vgrid) :: n0 = shape(vgrid,0)
    integer :: n1
    !f2py intent(hide), depend(vgrid) :: n1 = shape(vgrid,1)
    integer :: n2
    !f2py intent(hide), depend(vgrid) :: n2 = shape(vgrid,2)
    integer :: n3
    !f2py intent(hide), depend(vgrid) :: n3 = shape(vgrid,3)
    integer :: n4
    !f2py intent(hide), depend(vsphere) :: n4 = shape(vsphere,0)
    integer :: n5
    !f2py intent(hide), depend(vsphere) :: n5 = shape(vsphere,1)
    call iso_vsphere(Vgrid=vgrid, Vsphere=vsphere)
end subroutine f90wrap_iso_vsphere

subroutine f90wrap_iso_rho2grid(rhosphere, rhogrid, n0, n1, n2, n3, n4, n5)
    use grid_module, only: iso_rho2grid
    implicit none
    
    real(8), intent(in), dimension(n0,n1) :: rhosphere
    real(8), intent(inout), dimension(n2,n3,n4,n5) :: rhogrid
    integer :: n0
    !f2py intent(hide), depend(rhosphere) :: n0 = shape(rhosphere,0)
    integer :: n1
    !f2py intent(hide), depend(rhosphere) :: n1 = shape(rhosphere,1)
    integer :: n2
    !f2py intent(hide), depend(rhogrid) :: n2 = shape(rhogrid,0)
    integer :: n3
    !f2py intent(hide), depend(rhogrid) :: n3 = shape(rhogrid,1)
    integer :: n4
    !f2py intent(hide), depend(rhogrid) :: n4 = shape(rhogrid,2)
    integer :: n5
    !f2py intent(hide), depend(rhogrid) :: n5 = shape(rhogrid,3)
    call iso_rho2grid(Rhosphere=rhosphere, Rhogrid=rhogrid)
end subroutine f90wrap_iso_rho2grid

subroutine f90wrap_parallel_s2g(p, thrq, n0, n1, n2, n3)
    use grid_module, only: parallel_s2g
    implicit none
    
    real(8), dimension(n0) :: p
    real(8), dimension(n1,n2,n3) :: thrq
    integer :: n0
    !f2py intent(hide), depend(p) :: n0 = shape(p,0)
    integer :: n1
    !f2py intent(hide), depend(thrq) :: n1 = shape(thrq,0)
    integer :: n2
    !f2py intent(hide), depend(thrq) :: n2 = shape(thrq,1)
    integer :: n3
    !f2py intent(hide), depend(thrq) :: n3 = shape(thrq,2)
    call parallel_s2g(p=p, thrq=thrq)
end subroutine f90wrap_parallel_s2g

subroutine f90wrap_parallel_g2s(p, thrq, n0, n1, n2, n3)
    use grid_module, only: parallel_g2s
    implicit none
    
    real(8), dimension(n0) :: p
    real(8), dimension(n1,n2,n3) :: thrq
    integer :: n0
    !f2py intent(hide), depend(p) :: n0 = shape(p,0)
    integer :: n1
    !f2py intent(hide), depend(thrq) :: n1 = shape(thrq,0)
    integer :: n2
    !f2py intent(hide), depend(thrq) :: n2 = shape(thrq,1)
    integer :: n3
    !f2py intent(hide), depend(thrq) :: n3 = shape(thrq,2)
    call parallel_g2s(p=p, thrq=thrq)
end subroutine f90wrap_parallel_g2s

subroutine f90wrap_get_z_range(my_z_range, delta_z, ngrid_z, cell_thick, cell_shape)
    use grid_module, only: get_z_range
    implicit none
    
    integer(4), dimension(2) :: my_z_range
    real(8) :: delta_z
    integer(4) :: ngrid_z
    real(8) :: cell_thick
    integer(4) :: cell_shape
    call get_z_range(my_z_range=my_z_range, delta_z=delta_z, ngrid_z=ngrid_z, cell_thick=cell_thick, cell_shape=cell_shape)
end subroutine f90wrap_get_z_range

subroutine f90wrap_confirm_iso_radius
    use grid_module, only: confirm_iso_radius
    implicit none
    
    call confirm_iso_radius()
end subroutine f90wrap_confirm_iso_radius

subroutine f90wrap_reshape_center
    use grid_module, only: reshape_center
    implicit none
    
    call reshape_center()
end subroutine f90wrap_reshape_center

subroutine f90wrap_reset_poscar(new_gap, new_ni)
    use grid_module, only: reset_poscar
    implicit none
    
    real(8) :: new_gap
    integer(4) :: new_ni
    call reset_poscar(new_gap=new_gap, new_ni=new_ni)
end subroutine f90wrap_reset_poscar

subroutine f90wrap_build_parallel_cubic_grid
    use grid_module, only: build_parallel_cubic_grid
    implicit none
    
    call build_parallel_cubic_grid()
end subroutine f90wrap_build_parallel_cubic_grid

subroutine f90wrap_rho_trans1d(rho3d, rho1d, n0, n1, n2, n3)
    use grid_module, only: rho_trans1d
    implicit none
    
    real(8), intent(in), dimension(n0,n1,n2) :: rho3d
    real(8), intent(inout), dimension(n3) :: rho1d
    integer :: n0
    !f2py intent(hide), depend(rho3d) :: n0 = shape(rho3d,0)
    integer :: n1
    !f2py intent(hide), depend(rho3d) :: n1 = shape(rho3d,1)
    integer :: n2
    !f2py intent(hide), depend(rho3d) :: n2 = shape(rho3d,2)
    integer :: n3
    !f2py intent(hide), depend(rho1d) :: n3 = shape(rho1d,0)
    call rho_trans1d(rho3d=rho3d, rho1d=rho1d)
end subroutine f90wrap_rho_trans1d

subroutine f90wrap_rho_trans3d(rho1d, rho3d, n0, n1, n2, n3)
    use grid_module, only: rho_trans3d
    implicit none
    
    real(8), intent(in), dimension(n0) :: rho1d
    real(8), intent(inout), dimension(n1,n2,n3) :: rho3d
    integer :: n0
    !f2py intent(hide), depend(rho1d) :: n0 = shape(rho1d,0)
    integer :: n1
    !f2py intent(hide), depend(rho3d) :: n1 = shape(rho3d,0)
    integer :: n2
    !f2py intent(hide), depend(rho3d) :: n2 = shape(rho3d,1)
    integer :: n3
    !f2py intent(hide), depend(rho3d) :: n3 = shape(rho3d,2)
    call rho_trans3d(rho1d=rho1d, rho3d=rho3d)
end subroutine f90wrap_rho_trans3d

subroutine f90wrap_fft_sph_r2c(ret_array_c, array_r, n0, n1, n2, n3, n4, n5)
    use grid_module, only: fft_sph
    implicit none
    
    complex(8), intent(out), dimension(n0,n1,n2) :: ret_array_c
    real(8), intent(in), dimension(n3,n4,n5) :: array_r
    integer :: n0
    integer :: n1
    integer :: n2
    integer :: n3
    !f2py intent(hide), depend(array_r) :: n3 = shape(array_r,0)
    integer :: n4
    !f2py intent(hide), depend(array_r) :: n4 = shape(array_r,1)
    integer :: n5
    !f2py intent(hide), depend(array_r) :: n5 = shape(array_r,2)
    ret_array_c = fft_sph(array_r=array_r)
end subroutine f90wrap_fft_sph_r2c

subroutine f90wrap_fft_sph_c2r(ret_array_r, array_c, n0, n1, n2, n3, n4, n5)
    use grid_module, only: fft_sph
    implicit none
    
    real(8), intent(out), dimension(n0,n1,n2) :: ret_array_r
    complex(8), intent(in), dimension(n3,n4,n5) :: array_c
    integer :: n0
    integer :: n1
    integer :: n2
    integer :: n3
    !f2py intent(hide), depend(array_c) :: n3 = shape(array_c,0)
    integer :: n4
    !f2py intent(hide), depend(array_c) :: n4 = shape(array_c,1)
    integer :: n5
    !f2py intent(hide), depend(array_c) :: n5 = shape(array_c,2)
    ret_array_r = fft_sph(array_c=array_c)
end subroutine f90wrap_fft_sph_c2r

subroutine f90wrap_grid_module__get__n1(f90wrap_n1)
    use grid_module, only: grid_module_n1 => n1
    implicit none
    integer(4), intent(out) :: f90wrap_n1
    
    f90wrap_n1 = grid_module_n1
end subroutine f90wrap_grid_module__get__n1

subroutine f90wrap_grid_module__set__n1(f90wrap_n1)
    use grid_module, only: grid_module_n1 => n1
    implicit none
    integer(4), intent(in) :: f90wrap_n1
    
    grid_module_n1 = f90wrap_n1
end subroutine f90wrap_grid_module__set__n1

subroutine f90wrap_grid_module__get__n2(f90wrap_n2)
    use grid_module, only: grid_module_n2 => n2
    implicit none
    integer(4), intent(out) :: f90wrap_n2
    
    f90wrap_n2 = grid_module_n2
end subroutine f90wrap_grid_module__get__n2

subroutine f90wrap_grid_module__set__n2(f90wrap_n2)
    use grid_module, only: grid_module_n2 => n2
    implicit none
    integer(4), intent(in) :: f90wrap_n2
    
    grid_module_n2 = f90wrap_n2
end subroutine f90wrap_grid_module__set__n2

subroutine f90wrap_grid_module__get__n3(f90wrap_n3)
    use grid_module, only: grid_module_n3 => n3
    implicit none
    integer(4), intent(out) :: f90wrap_n3
    
    f90wrap_n3 = grid_module_n3
end subroutine f90wrap_grid_module__get__n3

subroutine f90wrap_grid_module__set__n3(f90wrap_n3)
    use grid_module, only: grid_module_n3 => n3
    implicit none
    integer(4), intent(in) :: f90wrap_n3
    
    grid_module_n3 = f90wrap_n3
end subroutine f90wrap_grid_module__set__n3

subroutine f90wrap_grid_module__get__n(f90wrap_n)
    use grid_module, only: grid_module_n => n
    implicit none
    integer(4), intent(out) :: f90wrap_n
    
    f90wrap_n = grid_module_n
end subroutine f90wrap_grid_module__get__n

subroutine f90wrap_grid_module__set__n(f90wrap_n)
    use grid_module, only: grid_module_n => n
    implicit none
    integer(4), intent(in) :: f90wrap_n
    
    grid_module_n = f90wrap_n
end subroutine f90wrap_grid_module__set__n

subroutine f90wrap_grid_module__get__nsn(f90wrap_nsn)
    use grid_module, only: grid_module_nsn => nsn
    implicit none
    integer(4), intent(out) :: f90wrap_nsn
    
    f90wrap_nsn = grid_module_nsn
end subroutine f90wrap_grid_module__get__nsn

subroutine f90wrap_grid_module__set__nsn(f90wrap_nsn)
    use grid_module, only: grid_module_nsn => nsn
    implicit none
    integer(4), intent(in) :: f90wrap_nsn
    
    grid_module_nsn = f90wrap_nsn
end subroutine f90wrap_grid_module__set__nsn

subroutine f90wrap_grid_module__get__global_n1(f90wrap_global_n1)
    use grid_module, only: grid_module_global_n1 => global_n1
    implicit none
    integer(4), intent(out) :: f90wrap_global_n1
    
    f90wrap_global_n1 = grid_module_global_n1
end subroutine f90wrap_grid_module__get__global_n1

subroutine f90wrap_grid_module__set__global_n1(f90wrap_global_n1)
    use grid_module, only: grid_module_global_n1 => global_n1
    implicit none
    integer(4), intent(in) :: f90wrap_global_n1
    
    grid_module_global_n1 = f90wrap_global_n1
end subroutine f90wrap_grid_module__set__global_n1

subroutine f90wrap_grid_module__get__global_n2(f90wrap_global_n2)
    use grid_module, only: grid_module_global_n2 => global_n2
    implicit none
    integer(4), intent(out) :: f90wrap_global_n2
    
    f90wrap_global_n2 = grid_module_global_n2
end subroutine f90wrap_grid_module__get__global_n2

subroutine f90wrap_grid_module__set__global_n2(f90wrap_global_n2)
    use grid_module, only: grid_module_global_n2 => global_n2
    implicit none
    integer(4), intent(in) :: f90wrap_global_n2
    
    grid_module_global_n2 = f90wrap_global_n2
end subroutine f90wrap_grid_module__set__global_n2

subroutine f90wrap_grid_module__get__global_n3(f90wrap_global_n3)
    use grid_module, only: grid_module_global_n3 => global_n3
    implicit none
    integer(4), intent(out) :: f90wrap_global_n3
    
    f90wrap_global_n3 = grid_module_global_n3
end subroutine f90wrap_grid_module__get__global_n3

subroutine f90wrap_grid_module__set__global_n3(f90wrap_global_n3)
    use grid_module, only: grid_module_global_n3 => global_n3
    implicit none
    integer(4), intent(in) :: f90wrap_global_n3
    
    grid_module_global_n3 = f90wrap_global_n3
end subroutine f90wrap_grid_module__set__global_n3

subroutine f90wrap_grid_module__get__global_n(f90wrap_global_n)
    use grid_module, only: grid_module_global_n => global_n
    implicit none
    integer(4), intent(out) :: f90wrap_global_n
    
    f90wrap_global_n = grid_module_global_n
end subroutine f90wrap_grid_module__get__global_n

subroutine f90wrap_grid_module__set__global_n(f90wrap_global_n)
    use grid_module, only: grid_module_global_n => global_n
    implicit none
    integer(4), intent(in) :: f90wrap_global_n
    
    grid_module_global_n = f90wrap_global_n
end subroutine f90wrap_grid_module__set__global_n

subroutine f90wrap_grid_module__get__ng1(f90wrap_ng1)
    use grid_module, only: grid_module_ng1 => ng1
    implicit none
    integer(4), intent(out) :: f90wrap_ng1
    
    f90wrap_ng1 = grid_module_ng1
end subroutine f90wrap_grid_module__get__ng1

subroutine f90wrap_grid_module__set__ng1(f90wrap_ng1)
    use grid_module, only: grid_module_ng1 => ng1
    implicit none
    integer(4), intent(in) :: f90wrap_ng1
    
    grid_module_ng1 = f90wrap_ng1
end subroutine f90wrap_grid_module__set__ng1

subroutine f90wrap_grid_module__get__ng2(f90wrap_ng2)
    use grid_module, only: grid_module_ng2 => ng2
    implicit none
    integer(4), intent(out) :: f90wrap_ng2
    
    f90wrap_ng2 = grid_module_ng2
end subroutine f90wrap_grid_module__get__ng2

subroutine f90wrap_grid_module__set__ng2(f90wrap_ng2)
    use grid_module, only: grid_module_ng2 => ng2
    implicit none
    integer(4), intent(in) :: f90wrap_ng2
    
    grid_module_ng2 = f90wrap_ng2
end subroutine f90wrap_grid_module__set__ng2

subroutine f90wrap_grid_module__get__ng3(f90wrap_ng3)
    use grid_module, only: grid_module_ng3 => ng3
    implicit none
    integer(4), intent(out) :: f90wrap_ng3
    
    f90wrap_ng3 = grid_module_ng3
end subroutine f90wrap_grid_module__get__ng3

subroutine f90wrap_grid_module__set__ng3(f90wrap_ng3)
    use grid_module, only: grid_module_ng3 => ng3
    implicit none
    integer(4), intent(in) :: f90wrap_ng3
    
    grid_module_ng3 = f90wrap_ng3
end subroutine f90wrap_grid_module__set__ng3

subroutine f90wrap_grid_module__get__ng(f90wrap_ng)
    use grid_module, only: grid_module_ng => ng
    implicit none
    integer(4), intent(out) :: f90wrap_ng
    
    f90wrap_ng = grid_module_ng
end subroutine f90wrap_grid_module__get__ng

subroutine f90wrap_grid_module__set__ng(f90wrap_ng)
    use grid_module, only: grid_module_ng => ng
    implicit none
    integer(4), intent(in) :: f90wrap_ng
    
    grid_module_ng = f90wrap_ng
end subroutine f90wrap_grid_module__set__ng

subroutine f90wrap_grid_module__array__gap(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use grid_module, only: grid_module_gap => gap
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    dshape(1:1) = shape(grid_module_gap)
    dloc = loc(grid_module_gap)
end subroutine f90wrap_grid_module__array__gap

subroutine f90wrap_grid_module__get__dvol(f90wrap_dvol)
    use grid_module, only: grid_module_dvol => dvol
    implicit none
    real(8), intent(out) :: f90wrap_dvol
    
    f90wrap_dvol = grid_module_dvol
end subroutine f90wrap_grid_module__get__dvol

subroutine f90wrap_grid_module__set__dvol(f90wrap_dvol)
    use grid_module, only: grid_module_dvol => dvol
    implicit none
    real(8), intent(in) :: f90wrap_dvol
    
    grid_module_dvol = f90wrap_dvol
end subroutine f90wrap_grid_module__set__dvol

subroutine f90wrap_grid_module__get__nk1(f90wrap_nk1)
    use grid_module, only: grid_module_nk1 => nk1
    implicit none
    integer(4), intent(out) :: f90wrap_nk1
    
    f90wrap_nk1 = grid_module_nk1
end subroutine f90wrap_grid_module__get__nk1

subroutine f90wrap_grid_module__set__nk1(f90wrap_nk1)
    use grid_module, only: grid_module_nk1 => nk1
    implicit none
    integer(4), intent(in) :: f90wrap_nk1
    
    grid_module_nk1 = f90wrap_nk1
end subroutine f90wrap_grid_module__set__nk1

subroutine f90wrap_grid_module__get__nk2(f90wrap_nk2)
    use grid_module, only: grid_module_nk2 => nk2
    implicit none
    integer(4), intent(out) :: f90wrap_nk2
    
    f90wrap_nk2 = grid_module_nk2
end subroutine f90wrap_grid_module__get__nk2

subroutine f90wrap_grid_module__set__nk2(f90wrap_nk2)
    use grid_module, only: grid_module_nk2 => nk2
    implicit none
    integer(4), intent(in) :: f90wrap_nk2
    
    grid_module_nk2 = f90wrap_nk2
end subroutine f90wrap_grid_module__set__nk2

subroutine f90wrap_grid_module__get__nk3(f90wrap_nk3)
    use grid_module, only: grid_module_nk3 => nk3
    implicit none
    integer(4), intent(out) :: f90wrap_nk3
    
    f90wrap_nk3 = grid_module_nk3
end subroutine f90wrap_grid_module__get__nk3

subroutine f90wrap_grid_module__set__nk3(f90wrap_nk3)
    use grid_module, only: grid_module_nk3 => nk3
    implicit none
    integer(4), intent(in) :: f90wrap_nk3
    
    grid_module_nk3 = f90wrap_nk3
end subroutine f90wrap_grid_module__set__nk3

subroutine f90wrap_grid_module__get__nk(f90wrap_nk)
    use grid_module, only: grid_module_nk => nk
    implicit none
    integer(4), intent(out) :: f90wrap_nk
    
    f90wrap_nk = grid_module_nk
end subroutine f90wrap_grid_module__get__nk

subroutine f90wrap_grid_module__set__nk(f90wrap_nk)
    use grid_module, only: grid_module_nk => nk
    implicit none
    integer(4), intent(in) :: f90wrap_nk
    
    grid_module_nk = f90wrap_nk
end subroutine f90wrap_grid_module__set__nk

subroutine f90wrap_grid_module__array__kdispl(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use grid_module, only: grid_module_kdispl => kdispl
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    dshape(1:1) = shape(grid_module_kdispl)
    dloc = loc(grid_module_kdispl)
end subroutine f90wrap_grid_module__array__kdispl

subroutine f90wrap_grid_module__array__dr(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use grid_module, only: grid_module_dr => dr
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    if (allocated(grid_module_dr)) then
        dshape(1:3) = shape(grid_module_dr)
        dloc = loc(grid_module_dr)
    else
        dloc = 0
    end if
end subroutine f90wrap_grid_module__array__dr

subroutine f90wrap_grid_module__array__cost(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use grid_module, only: grid_module_cost => cost
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    if (allocated(grid_module_cost)) then
        dshape(1:3) = shape(grid_module_cost)
        dloc = loc(grid_module_cost)
    else
        dloc = 0
    end if
end subroutine f90wrap_grid_module__array__cost

subroutine f90wrap_grid_module__array__sint(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use grid_module, only: grid_module_sint => sint
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    if (allocated(grid_module_sint)) then
        dshape(1:3) = shape(grid_module_sint)
        dloc = loc(grid_module_sint)
    else
        dloc = 0
    end if
end subroutine f90wrap_grid_module__array__sint

subroutine f90wrap_grid_module__array__cns(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use grid_module, only: grid_module_cns => cns
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 15
    if (allocated(grid_module_cns)) then
        dshape(1:3) = shape(grid_module_cns)
        dloc = loc(grid_module_cns)
    else
        dloc = 0
    end if
end subroutine f90wrap_grid_module__array__cns

subroutine f90wrap_grid_module__array__indx(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use grid_module, only: grid_module_indx => indx
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 5
    if (allocated(grid_module_indx)) then
        dshape(1:3) = shape(grid_module_indx)
        dloc = loc(grid_module_indx)
    else
        dloc = 0
    end if
end subroutine f90wrap_grid_module__array__indx

! End of module grid_module defined in file Grid_module.fpp


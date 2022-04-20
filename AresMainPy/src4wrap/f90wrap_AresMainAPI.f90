! Module aresmainapi defined in file AresMainAPI.fpp

subroutine f90wrap_aresout__array__forces(this, nd, dtype, dshape, dloc)
    use aresmainapi, only: aresout
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

subroutine f90wrap_aresout__array__poscar(this, nd, dtype, dshape, dloc)
    use aresmainapi, only: aresout
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
    use aresmainapi, only: aresout
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
    use aresmainapi, only: aresout
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

subroutine f90wrap_aresout__array__stress(this, nd, dtype, dshape, dloc)
    use aresmainapi, only: aresout
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
    dshape(1:2) = shape(this_ptr%p%stress)
    dloc = loc(this_ptr%p%stress)
end subroutine f90wrap_aresout__array__stress

subroutine f90wrap_aresout__array__apilat_mat(this, nd, dtype, dshape, dloc)
    use aresmainapi, only: aresout
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
    dshape(1:2) = shape(this_ptr%p%apilat_mat)
    dloc = loc(this_ptr%p%apilat_mat)
end subroutine f90wrap_aresout__array__apilat_mat

subroutine f90wrap_aresout__array__apilat_para(this, nd, dtype, dshape, dloc)
    use aresmainapi, only: aresout
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
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%apilat_para)
    dloc = loc(this_ptr%p%apilat_para)
end subroutine f90wrap_aresout__array__apilat_para

subroutine f90wrap_aresout__get__comm(this, f90wrap_comm)
    use aresmainapi, only: aresout
    implicit none
    type aresout_ptr_type
        type(aresout), pointer :: p => NULL()
    end type aresout_ptr_type
    integer, intent(in)   :: this(2)
    type(aresout_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_comm
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_comm = this_ptr%p%comm
end subroutine f90wrap_aresout__get__comm

subroutine f90wrap_aresout__set__comm(this, f90wrap_comm)
    use aresmainapi, only: aresout
    implicit none
    type aresout_ptr_type
        type(aresout), pointer :: p => NULL()
    end type aresout_ptr_type
    integer, intent(in)   :: this(2)
    type(aresout_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_comm
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%comm = f90wrap_comm
end subroutine f90wrap_aresout__set__comm

subroutine f90wrap_aresout__get__myid(this, f90wrap_myid)
    use aresmainapi, only: aresout
    implicit none
    type aresout_ptr_type
        type(aresout), pointer :: p => NULL()
    end type aresout_ptr_type
    integer, intent(in)   :: this(2)
    type(aresout_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_myid
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_myid = this_ptr%p%myid
end subroutine f90wrap_aresout__get__myid

subroutine f90wrap_aresout__set__myid(this, f90wrap_myid)
    use aresmainapi, only: aresout
    implicit none
    type aresout_ptr_type
        type(aresout), pointer :: p => NULL()
    end type aresout_ptr_type
    integer, intent(in)   :: this(2)
    type(aresout_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_myid
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%myid = f90wrap_myid
end subroutine f90wrap_aresout__set__myid

subroutine f90wrap_aresout__get__numprocs(this, f90wrap_numprocs)
    use aresmainapi, only: aresout
    implicit none
    type aresout_ptr_type
        type(aresout), pointer :: p => NULL()
    end type aresout_ptr_type
    integer, intent(in)   :: this(2)
    type(aresout_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_numprocs
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_numprocs = this_ptr%p%numprocs
end subroutine f90wrap_aresout__get__numprocs

subroutine f90wrap_aresout__set__numprocs(this, f90wrap_numprocs)
    use aresmainapi, only: aresout
    implicit none
    type aresout_ptr_type
        type(aresout), pointer :: p => NULL()
    end type aresout_ptr_type
    integer, intent(in)   :: this(2)
    type(aresout_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_numprocs
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%numprocs = f90wrap_numprocs
end subroutine f90wrap_aresout__set__numprocs

subroutine f90wrap_aresout__get__rootid(this, f90wrap_rootid)
    use aresmainapi, only: aresout
    implicit none
    type aresout_ptr_type
        type(aresout), pointer :: p => NULL()
    end type aresout_ptr_type
    integer, intent(in)   :: this(2)
    type(aresout_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_rootid
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_rootid = this_ptr%p%rootid
end subroutine f90wrap_aresout__get__rootid

subroutine f90wrap_aresout__set__rootid(this, f90wrap_rootid)
    use aresmainapi, only: aresout
    implicit none
    type aresout_ptr_type
        type(aresout), pointer :: p => NULL()
    end type aresout_ptr_type
    integer, intent(in)   :: this(2)
    type(aresout_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_rootid
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%rootid = f90wrap_rootid
end subroutine f90wrap_aresout__set__rootid

subroutine f90wrap_aresout_initialise(this)
    use aresmainapi, only: aresout
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
    use aresmainapi, only: aresout
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
    use aresmainapi, only: init_alloc_arrays, aresout
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

subroutine f90wrap_assignment(dertype)
    use aresmainapi, only: assignment, aresout
    implicit none
    
    type aresout_ptr_type
        type(aresout), pointer :: p => NULL()
    end type aresout_ptr_type
    type(aresout_ptr_type) :: dertype_ptr
    integer, intent(in), dimension(2) :: dertype
    dertype_ptr = transfer(dertype, dertype_ptr)
    call assignment(dertype=dertype_ptr%p)
end subroutine f90wrap_assignment

subroutine f90wrap_destroy_alloc_arrays(dertype)
    use aresmainapi, only: aresout, destroy_alloc_arrays
    implicit none
    
    type aresout_ptr_type
        type(aresout), pointer :: p => NULL()
    end type aresout_ptr_type
    type(aresout_ptr_type) :: dertype_ptr
    integer, intent(in), dimension(2) :: dertype
    dertype_ptr = transfer(dertype, dertype_ptr)
    call destroy_alloc_arrays(dertype=dertype_ptr%p)
end subroutine f90wrap_destroy_alloc_arrays

subroutine f90wrap_updateions(pos, lattice, n0, n1, n2)
    use aresmainapi, only: updateions
    implicit none
    
    real(8), intent(in), dimension(n0,n1) :: pos
    real(8), intent(in), optional, dimension(3,n2) :: lattice
    integer :: n0
    !f2py intent(hide), depend(pos) :: n0 = shape(pos,0)
    integer :: n1
    !f2py intent(hide), depend(pos) :: n1 = shape(pos,1)
    integer :: n2
    !f2py intent(hide), depend(lattice) :: n2 = shape(lattice,1)
    call updateions(pos=pos, lattice=lattice)
end subroutine f90wrap_updateions

! End of module aresmainapi defined in file AresMainAPI.fpp


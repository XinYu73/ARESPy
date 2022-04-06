! Module smpi_math_module defined in file Smpi_math_module.fpp

subroutine f90wrap_parallel_type__get__comm(this, f90wrap_comm)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in)   :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_comm

    this_ptr = transfer(this, this_ptr)
    f90wrap_comm = this_ptr%p%comm
end subroutine f90wrap_parallel_type__get__comm

subroutine f90wrap_parallel_type__set__comm(this, f90wrap_comm)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in)   :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_comm

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%comm = f90wrap_comm
end subroutine f90wrap_parallel_type__set__comm

subroutine f90wrap_parallel_type__get__myid(this, f90wrap_myid)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in)   :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_myid

    this_ptr = transfer(this, this_ptr)
    f90wrap_myid = this_ptr%p%myid
end subroutine f90wrap_parallel_type__get__myid

subroutine f90wrap_parallel_type__set__myid(this, f90wrap_myid)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in)   :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_myid

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%myid = f90wrap_myid
end subroutine f90wrap_parallel_type__set__myid

subroutine f90wrap_parallel_type__get__numprocs(this, f90wrap_numprocs)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in)   :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_numprocs

    this_ptr = transfer(this, this_ptr)
    f90wrap_numprocs = this_ptr%p%numprocs
end subroutine f90wrap_parallel_type__get__numprocs

subroutine f90wrap_parallel_type__set__numprocs(this, f90wrap_numprocs)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in)   :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_numprocs

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%numprocs = f90wrap_numprocs
end subroutine f90wrap_parallel_type__set__numprocs

subroutine f90wrap_parallel_type__get__rootid(this, f90wrap_rootid)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in)   :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_rootid

    this_ptr = transfer(this, this_ptr)
    f90wrap_rootid = this_ptr%p%rootid
end subroutine f90wrap_parallel_type__get__rootid

subroutine f90wrap_parallel_type__set__rootid(this, f90wrap_rootid)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in)   :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_rootid

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%rootid = f90wrap_rootid
end subroutine f90wrap_parallel_type__set__rootid

subroutine f90wrap_parallel_type__get__isroot(this, f90wrap_isroot)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in)   :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_isroot

    this_ptr = transfer(this, this_ptr)
    f90wrap_isroot = this_ptr%p%isroot
end subroutine f90wrap_parallel_type__get__isroot

subroutine f90wrap_parallel_type__set__isroot(this, f90wrap_isroot)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in)   :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_isroot

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%isroot = f90wrap_isroot
end subroutine f90wrap_parallel_type__set__isroot

subroutine f90wrap_parallel_type__get__nstate_proc(this, f90wrap_nstate_proc)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in)   :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_nstate_proc

    this_ptr = transfer(this, this_ptr)
    f90wrap_nstate_proc = this_ptr%p%nstate_proc
end subroutine f90wrap_parallel_type__get__nstate_proc

subroutine f90wrap_parallel_type__set__nstate_proc(this, f90wrap_nstate_proc)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in)   :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_nstate_proc

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nstate_proc = f90wrap_nstate_proc
end subroutine f90wrap_parallel_type__set__nstate_proc

subroutine f90wrap_parallel_type__array__sub2sum(this, nd, dtype, dshape, dloc)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in) :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 2
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%sub2sum)) then
        dshape(1:2) = shape(this_ptr%p%sub2sum)
        dloc = loc(this_ptr%p%sub2sum)
    else
        dloc = 0
    end if
end subroutine f90wrap_parallel_type__array__sub2sum

subroutine f90wrap_parallel_type__array__mygrid_range(this, nd, dtype, dshape, dloc)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in) :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%mygrid_range)
    dloc = loc(this_ptr%p%mygrid_range)
end subroutine f90wrap_parallel_type__array__mygrid_range

subroutine f90wrap_parallel_type__array__recvcounts(this, nd, dtype, dshape, dloc)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in) :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%recvcounts)) then
        dshape(1:1) = shape(this_ptr%p%recvcounts)
        dloc = loc(this_ptr%p%recvcounts)
    else
        dloc = 0
    end if
end subroutine f90wrap_parallel_type__array__recvcounts

subroutine f90wrap_parallel_type__array__displs(this, nd, dtype, dshape, dloc)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in) :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%displs)) then
        dshape(1:1) = shape(this_ptr%p%displs)
        dloc = loc(this_ptr%p%displs)
    else
        dloc = 0
    end if
end subroutine f90wrap_parallel_type__array__displs

subroutine f90wrap_parallel_type__array__global_gridrange(this, nd, dtype, dshape, dloc)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in) :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 2
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%global_gridrange)) then
        dshape(1:2) = shape(this_ptr%p%global_gridrange)
        dloc = loc(this_ptr%p%global_gridrange)
    else
        dloc = 0
    end if
end subroutine f90wrap_parallel_type__array__global_gridrange

subroutine f90wrap_parallel_type__get__comm2d(this, f90wrap_comm2d)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in)   :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_comm2d

    this_ptr = transfer(this, this_ptr)
    f90wrap_comm2d = this_ptr%p%comm2d
end subroutine f90wrap_parallel_type__get__comm2d

subroutine f90wrap_parallel_type__set__comm2d(this, f90wrap_comm2d)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in)   :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_comm2d

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%comm2d = f90wrap_comm2d
end subroutine f90wrap_parallel_type__set__comm2d

subroutine f90wrap_parallel_type__get__commx(this, f90wrap_commx)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in)   :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_commx

    this_ptr = transfer(this, this_ptr)
    f90wrap_commx = this_ptr%p%commx
end subroutine f90wrap_parallel_type__get__commx

subroutine f90wrap_parallel_type__set__commx(this, f90wrap_commx)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in)   :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_commx

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%commx = f90wrap_commx
end subroutine f90wrap_parallel_type__set__commx

subroutine f90wrap_parallel_type__get__commy(this, f90wrap_commy)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in)   :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_commy

    this_ptr = transfer(this, this_ptr)
    f90wrap_commy = this_ptr%p%commy
end subroutine f90wrap_parallel_type__get__commy

subroutine f90wrap_parallel_type__set__commy(this, f90wrap_commy)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in)   :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_commy

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%commy = f90wrap_commy
end subroutine f90wrap_parallel_type__set__commy

subroutine f90wrap_parallel_type__get__rankx(this, f90wrap_rankx)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in)   :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_rankx

    this_ptr = transfer(this, this_ptr)
    f90wrap_rankx = this_ptr%p%rankx
end subroutine f90wrap_parallel_type__get__rankx

subroutine f90wrap_parallel_type__set__rankx(this, f90wrap_rankx)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in)   :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_rankx

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%rankx = f90wrap_rankx
end subroutine f90wrap_parallel_type__set__rankx

subroutine f90wrap_parallel_type__get__ranky(this, f90wrap_ranky)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in)   :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_ranky

    this_ptr = transfer(this, this_ptr)
    f90wrap_ranky = this_ptr%p%ranky
end subroutine f90wrap_parallel_type__get__ranky

subroutine f90wrap_parallel_type__set__ranky(this, f90wrap_ranky)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in)   :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_ranky

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ranky = f90wrap_ranky
end subroutine f90wrap_parallel_type__set__ranky

subroutine f90wrap_parallel_type__array__periods(this, nd, dtype, dshape, dloc)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in) :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%periods)
    dloc = loc(this_ptr%p%periods)
end subroutine f90wrap_parallel_type__array__periods

subroutine f90wrap_parallel_type__get__reorder(this, f90wrap_reorder)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in)   :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_reorder

    this_ptr = transfer(this, this_ptr)
    f90wrap_reorder = this_ptr%p%reorder
end subroutine f90wrap_parallel_type__get__reorder

subroutine f90wrap_parallel_type__set__reorder(this, f90wrap_reorder)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in)   :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_reorder

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%reorder = f90wrap_reorder
end subroutine f90wrap_parallel_type__set__reorder

subroutine f90wrap_parallel_type__array__remainX(this, nd, dtype, dshape, dloc)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in) :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%remainX)
    dloc = loc(this_ptr%p%remainX)
end subroutine f90wrap_parallel_type__array__remainX

subroutine f90wrap_parallel_type__array__remainY(this, nd, dtype, dshape, dloc)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in) :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%remainY)
    dloc = loc(this_ptr%p%remainY)
end subroutine f90wrap_parallel_type__array__remainY

subroutine f90wrap_parallel_type__get__ndims(this, f90wrap_ndims)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in)   :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_ndims

    this_ptr = transfer(this, this_ptr)
    f90wrap_ndims = this_ptr%p%ndims
end subroutine f90wrap_parallel_type__get__ndims

subroutine f90wrap_parallel_type__set__ndims(this, f90wrap_ndims)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in)   :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_ndims

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ndims = f90wrap_ndims
end subroutine f90wrap_parallel_type__set__ndims

subroutine f90wrap_parallel_type__array__dims(this, nd, dtype, dshape, dloc)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in) :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%dims)
    dloc = loc(this_ptr%p%dims)
end subroutine f90wrap_parallel_type__array__dims

subroutine f90wrap_parallel_type__get__commfft(this, f90wrap_commfft)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in)   :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_commfft

    this_ptr = transfer(this, this_ptr)
    f90wrap_commfft = this_ptr%p%commfft
end subroutine f90wrap_parallel_type__get__commfft

subroutine f90wrap_parallel_type__set__commfft(this, f90wrap_commfft)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in)   :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_commfft

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%commfft = f90wrap_commfft
end subroutine f90wrap_parallel_type__set__commfft

subroutine f90wrap_parallel_type__get__local_z(this, f90wrap_local_z)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in)   :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_local_z

    this_ptr = transfer(this, this_ptr)
    f90wrap_local_z = this_ptr%p%local_z
end subroutine f90wrap_parallel_type__get__local_z

subroutine f90wrap_parallel_type__set__local_z(this, f90wrap_local_z)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in)   :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_local_z

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%local_z = f90wrap_local_z
end subroutine f90wrap_parallel_type__set__local_z

subroutine f90wrap_parallel_type__get__local_z_start(this, f90wrap_local_z_start)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in)   :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_local_z_start

    this_ptr = transfer(this, this_ptr)
    f90wrap_local_z_start = this_ptr%p%local_z_start
end subroutine f90wrap_parallel_type__get__local_z_start

subroutine f90wrap_parallel_type__set__local_z_start(this, f90wrap_local_z_start)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in)   :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_local_z_start

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%local_z_start = f90wrap_local_z_start
end subroutine f90wrap_parallel_type__set__local_z_start

subroutine f90wrap_parallel_type__array__fft_grid_range(this, nd, dtype, dshape, dloc)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in) :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 2
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%fft_grid_range)) then
        dshape(1:2) = shape(this_ptr%p%fft_grid_range)
        dloc = loc(this_ptr%p%fft_grid_range)
    else
        dloc = 0
    end if
end subroutine f90wrap_parallel_type__array__fft_grid_range

subroutine f90wrap_parallel_type__array__fft_rcount(this, nd, dtype, dshape, dloc)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in) :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%fft_rcount)) then
        dshape(1:1) = shape(this_ptr%p%fft_rcount)
        dloc = loc(this_ptr%p%fft_rcount)
    else
        dloc = 0
    end if
end subroutine f90wrap_parallel_type__array__fft_rcount

subroutine f90wrap_parallel_type__array__fft_rdispls(this, nd, dtype, dshape, dloc)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in) :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%fft_rdispls)) then
        dshape(1:1) = shape(this_ptr%p%fft_rdispls)
        dloc = loc(this_ptr%p%fft_rdispls)
    else
        dloc = 0
    end if
end subroutine f90wrap_parallel_type__array__fft_rdispls

subroutine f90wrap_parallel_type__array__fft_scount(this, nd, dtype, dshape, dloc)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in) :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%fft_scount)) then
        dshape(1:1) = shape(this_ptr%p%fft_scount)
        dloc = loc(this_ptr%p%fft_scount)
    else
        dloc = 0
    end if
end subroutine f90wrap_parallel_type__array__fft_scount

subroutine f90wrap_parallel_type__array__fft_sdispls(this, nd, dtype, dshape, dloc)
    use smpi_math_module, only: parallel_type
    implicit none
    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    integer, intent(in) :: this(2)
    type(parallel_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%fft_sdispls)) then
        dshape(1:1) = shape(this_ptr%p%fft_sdispls)
        dloc = loc(this_ptr%p%fft_sdispls)
    else
        dloc = 0
    end if
end subroutine f90wrap_parallel_type__array__fft_sdispls

subroutine f90wrap_parallel_type_initialise(this)
    use smpi_math_module, only: parallel_type
    implicit none

    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    type(parallel_type_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate (this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_parallel_type_initialise

subroutine f90wrap_parallel_type_finalise(this)
    use smpi_math_module, only: parallel_type
    implicit none

    type parallel_type_ptr_type
        type(parallel_type), pointer :: p => NULL()
    end type parallel_type_ptr_type
    type(parallel_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate (this_ptr%p)
end subroutine f90wrap_parallel_type_finalise

subroutine f90wrap_smpi_root_type__array__natom_group(this, nd, dtype, dshape, dloc)
    use smpi_math_module, only: smpi_root_type
    implicit none
    type smpi_root_type_ptr_type
        type(smpi_root_type), pointer :: p => NULL()
    end type smpi_root_type_ptr_type
    integer, intent(in) :: this(2)
    type(smpi_root_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 2
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%natom_group)) then
        dshape(1:2) = shape(this_ptr%p%natom_group)
        dloc = loc(this_ptr%p%natom_group)
    else
        dloc = 0
    end if
end subroutine f90wrap_smpi_root_type__array__natom_group

subroutine f90wrap_smpi_root_type_initialise(this)
    use smpi_math_module, only: smpi_root_type
    implicit none

    type smpi_root_type_ptr_type
        type(smpi_root_type), pointer :: p => NULL()
    end type smpi_root_type_ptr_type
    type(smpi_root_type_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate (this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_smpi_root_type_initialise

subroutine f90wrap_smpi_root_type_finalise(this)
    use smpi_math_module, only: smpi_root_type
    implicit none

    type smpi_root_type_ptr_type
        type(smpi_root_type), pointer :: p => NULL()
    end type smpi_root_type_ptr_type
    type(smpi_root_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate (this_ptr%p)
end subroutine f90wrap_smpi_root_type_finalise

subroutine f90wrap_smpi_comm_type__array__atoms(this, nd, dtype, dshape, dloc)
    use smpi_math_module, only: smpi_comm_type
    implicit none
    type smpi_comm_type_ptr_type
        type(smpi_comm_type), pointer :: p => NULL()
    end type smpi_comm_type_ptr_type
    integer, intent(in) :: this(2)
    type(smpi_comm_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%atoms)) then
        dshape(1:1) = shape(this_ptr%p%atoms)
        dloc = loc(this_ptr%p%atoms)
    else
        dloc = 0
    end if
end subroutine f90wrap_smpi_comm_type__array__atoms

subroutine f90wrap_smpi_comm_type__array__displs(this, nd, dtype, dshape, dloc)
    use smpi_math_module, only: smpi_comm_type
    implicit none
    type smpi_comm_type_ptr_type
        type(smpi_comm_type), pointer :: p => NULL()
    end type smpi_comm_type_ptr_type
    integer, intent(in) :: this(2)
    type(smpi_comm_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%displs)) then
        dshape(1:1) = shape(this_ptr%p%displs)
        dloc = loc(this_ptr%p%displs)
    else
        dloc = 0
    end if
end subroutine f90wrap_smpi_comm_type__array__displs

subroutine f90wrap_smpi_comm_type_initialise(this)
    use smpi_math_module, only: smpi_comm_type
    implicit none

    type smpi_comm_type_ptr_type
        type(smpi_comm_type), pointer :: p => NULL()
    end type smpi_comm_type_ptr_type
    type(smpi_comm_type_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate (this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_smpi_comm_type_initialise

subroutine f90wrap_smpi_comm_type_finalise(this)
    use smpi_math_module, only: smpi_comm_type
    implicit none

    type smpi_comm_type_ptr_type
        type(smpi_comm_type), pointer :: p => NULL()
    end type smpi_comm_type_ptr_type
    type(smpi_comm_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate (this_ptr%p)
end subroutine f90wrap_smpi_comm_type_finalise

subroutine f90wrap_time_type__get__label(this, f90wrap_label)
    use smpi_math_module, only: time_type
    implicit none
    type time_type_ptr_type
        type(time_type), pointer :: p => NULL()
    end type time_type_ptr_type
    integer, intent(in)   :: this(2)
    type(time_type_ptr_type) :: this_ptr
    character(100), intent(out) :: f90wrap_label

    this_ptr = transfer(this, this_ptr)
    f90wrap_label = this_ptr%p%label
end subroutine f90wrap_time_type__get__label

subroutine f90wrap_time_type__set__label(this, f90wrap_label)
    use smpi_math_module, only: time_type
    implicit none
    type time_type_ptr_type
        type(time_type), pointer :: p => NULL()
    end type time_type_ptr_type
    integer, intent(in)   :: this(2)
    type(time_type_ptr_type) :: this_ptr
    character(100), intent(in) :: f90wrap_label

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%label = f90wrap_label
end subroutine f90wrap_time_type__set__label

subroutine f90wrap_time_type__get__tic(this, f90wrap_tic)
    use smpi_math_module, only: time_type
    implicit none
    type time_type_ptr_type
        type(time_type), pointer :: p => NULL()
    end type time_type_ptr_type
    integer, intent(in)   :: this(2)
    type(time_type_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_tic

    this_ptr = transfer(this, this_ptr)
    f90wrap_tic = this_ptr%p%tic
end subroutine f90wrap_time_type__get__tic

subroutine f90wrap_time_type__set__tic(this, f90wrap_tic)
    use smpi_math_module, only: time_type
    implicit none
    type time_type_ptr_type
        type(time_type), pointer :: p => NULL()
    end type time_type_ptr_type
    integer, intent(in)   :: this(2)
    type(time_type_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_tic

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%tic = f90wrap_tic
end subroutine f90wrap_time_type__set__tic

subroutine f90wrap_time_type__get__toc(this, f90wrap_toc)
    use smpi_math_module, only: time_type
    implicit none
    type time_type_ptr_type
        type(time_type), pointer :: p => NULL()
    end type time_type_ptr_type
    integer, intent(in)   :: this(2)
    type(time_type_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_toc

    this_ptr = transfer(this, this_ptr)
    f90wrap_toc = this_ptr%p%toc
end subroutine f90wrap_time_type__get__toc

subroutine f90wrap_time_type__set__toc(this, f90wrap_toc)
    use smpi_math_module, only: time_type
    implicit none
    type time_type_ptr_type
        type(time_type), pointer :: p => NULL()
    end type time_type_ptr_type
    integer, intent(in)   :: this(2)
    type(time_type_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_toc

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%toc = f90wrap_toc
end subroutine f90wrap_time_type__set__toc

subroutine f90wrap_time_type__get__total(this, f90wrap_total)
    use smpi_math_module, only: time_type
    implicit none
    type time_type_ptr_type
        type(time_type), pointer :: p => NULL()
    end type time_type_ptr_type
    integer, intent(in)   :: this(2)
    type(time_type_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_total

    this_ptr = transfer(this, this_ptr)
    f90wrap_total = this_ptr%p%total
end subroutine f90wrap_time_type__get__total

subroutine f90wrap_time_type__set__total(this, f90wrap_total)
    use smpi_math_module, only: time_type
    implicit none
    type time_type_ptr_type
        type(time_type), pointer :: p => NULL()
    end type time_type_ptr_type
    integer, intent(in)   :: this(2)
    type(time_type_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_total

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%total = f90wrap_total
end subroutine f90wrap_time_type__set__total

subroutine f90wrap_time_type__get__sum_total(this, f90wrap_sum_total)
    use smpi_math_module, only: time_type
    implicit none
    type time_type_ptr_type
        type(time_type), pointer :: p => NULL()
    end type time_type_ptr_type
    integer, intent(in)   :: this(2)
    type(time_type_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_sum_total

    this_ptr = transfer(this, this_ptr)
    f90wrap_sum_total = this_ptr%p%sum_total
end subroutine f90wrap_time_type__get__sum_total

subroutine f90wrap_time_type__set__sum_total(this, f90wrap_sum_total)
    use smpi_math_module, only: time_type
    implicit none
    type time_type_ptr_type
        type(time_type), pointer :: p => NULL()
    end type time_type_ptr_type
    integer, intent(in)   :: this(2)
    type(time_type_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_sum_total

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%sum_total = f90wrap_sum_total
end subroutine f90wrap_time_type__set__sum_total

subroutine f90wrap_time_type__get__num(this, f90wrap_num)
    use smpi_math_module, only: time_type
    implicit none
    type time_type_ptr_type
        type(time_type), pointer :: p => NULL()
    end type time_type_ptr_type
    integer, intent(in)   :: this(2)
    type(time_type_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_num

    this_ptr = transfer(this, this_ptr)
    f90wrap_num = this_ptr%p%num
end subroutine f90wrap_time_type__get__num

subroutine f90wrap_time_type__set__num(this, f90wrap_num)
    use smpi_math_module, only: time_type
    implicit none
    type time_type_ptr_type
        type(time_type), pointer :: p => NULL()
    end type time_type_ptr_type
    integer, intent(in)   :: this(2)
    type(time_type_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_num

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%num = f90wrap_num
end subroutine f90wrap_time_type__set__num

subroutine f90wrap_time_type_initialise(this)
    use smpi_math_module, only: time_type
    implicit none

    type time_type_ptr_type
        type(time_type), pointer :: p => NULL()
    end type time_type_ptr_type
    type(time_type_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate (this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_time_type_initialise

subroutine f90wrap_time_type_finalise(this)
    use smpi_math_module, only: time_type
    implicit none

    type time_type_ptr_type
        type(time_type), pointer :: p => NULL()
    end type time_type_ptr_type
    type(time_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate (this_ptr%p)
end subroutine f90wrap_time_type_finalise

subroutine f90wrap_mem_type__get__label(this, f90wrap_label)
    use smpi_math_module, only: mem_type
    implicit none
    type mem_type_ptr_type
        type(mem_type), pointer :: p => NULL()
    end type mem_type_ptr_type
    integer, intent(in)   :: this(2)
    type(mem_type_ptr_type) :: this_ptr
    character(100), intent(out) :: f90wrap_label

    this_ptr = transfer(this, this_ptr)
    f90wrap_label = this_ptr%p%label
end subroutine f90wrap_mem_type__get__label

subroutine f90wrap_mem_type__set__label(this, f90wrap_label)
    use smpi_math_module, only: mem_type
    implicit none
    type mem_type_ptr_type
        type(mem_type), pointer :: p => NULL()
    end type mem_type_ptr_type
    integer, intent(in)   :: this(2)
    type(mem_type_ptr_type) :: this_ptr
    character(100), intent(in) :: f90wrap_label

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%label = f90wrap_label
end subroutine f90wrap_mem_type__set__label

subroutine f90wrap_mem_type__get__memic(this, f90wrap_memic)
    use smpi_math_module, only: mem_type
    implicit none
    type mem_type_ptr_type
        type(mem_type), pointer :: p => NULL()
    end type mem_type_ptr_type
    integer, intent(in)   :: this(2)
    type(mem_type_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_memic

    this_ptr = transfer(this, this_ptr)
    f90wrap_memic = this_ptr%p%memic
end subroutine f90wrap_mem_type__get__memic

subroutine f90wrap_mem_type__set__memic(this, f90wrap_memic)
    use smpi_math_module, only: mem_type
    implicit none
    type mem_type_ptr_type
        type(mem_type), pointer :: p => NULL()
    end type mem_type_ptr_type
    integer, intent(in)   :: this(2)
    type(mem_type_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_memic

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%memic = f90wrap_memic
end subroutine f90wrap_mem_type__set__memic

subroutine f90wrap_mem_type__get__total(this, f90wrap_total)
    use smpi_math_module, only: mem_type
    implicit none
    type mem_type_ptr_type
        type(mem_type), pointer :: p => NULL()
    end type mem_type_ptr_type
    integer, intent(in)   :: this(2)
    type(mem_type_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_total

    this_ptr = transfer(this, this_ptr)
    f90wrap_total = this_ptr%p%total
end subroutine f90wrap_mem_type__get__total

subroutine f90wrap_mem_type__set__total(this, f90wrap_total)
    use smpi_math_module, only: mem_type
    implicit none
    type mem_type_ptr_type
        type(mem_type), pointer :: p => NULL()
    end type mem_type_ptr_type
    integer, intent(in)   :: this(2)
    type(mem_type_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_total

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%total = f90wrap_total
end subroutine f90wrap_mem_type__set__total

subroutine f90wrap_mem_type__get__num(this, f90wrap_num)
    use smpi_math_module, only: mem_type
    implicit none
    type mem_type_ptr_type
        type(mem_type), pointer :: p => NULL()
    end type mem_type_ptr_type
    integer, intent(in)   :: this(2)
    type(mem_type_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_num

    this_ptr = transfer(this, this_ptr)
    f90wrap_num = this_ptr%p%num
end subroutine f90wrap_mem_type__get__num

subroutine f90wrap_mem_type__set__num(this, f90wrap_num)
    use smpi_math_module, only: mem_type
    implicit none
    type mem_type_ptr_type
        type(mem_type), pointer :: p => NULL()
    end type mem_type_ptr_type
    integer, intent(in)   :: this(2)
    type(mem_type_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_num

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%num = f90wrap_num
end subroutine f90wrap_mem_type__set__num

subroutine f90wrap_mem_type_initialise(this)
    use smpi_math_module, only: mem_type
    implicit none

    type mem_type_ptr_type
        type(mem_type), pointer :: p => NULL()
    end type mem_type_ptr_type
    type(mem_type_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate (this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_mem_type_initialise

subroutine f90wrap_mem_type_finalise(this)
    use smpi_math_module, only: mem_type
    implicit none

    type mem_type_ptr_type
        type(mem_type), pointer :: p => NULL()
    end type mem_type_ptr_type
    type(mem_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate (this_ptr%p)
end subroutine f90wrap_mem_type_finalise

subroutine f90wrap_grid_diff_map_type__array__nz_map(this, nd, dtype, dshape, dloc)
    use smpi_math_module, only: grid_diff_map_type
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
    use smpi_math_module, only: grid_diff_map_type
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
    use smpi_math_module, only: grid_diff_map_type
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
    use smpi_math_module, only: grid_diff_map_type
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
    use smpi_math_module, only: grid_diff_map_type
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

subroutine f90wrap_grid_diff_map_type__array__local_map1d(this, nd, dtype, dshape, dloc)
    use smpi_math_module, only: grid_diff_map_type
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

    nd = 3
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%local_map1d)) then
        dshape(1:3) = shape(this_ptr%p%local_map1d)
        dloc = loc(this_ptr%p%local_map1d)
    else
        dloc = 0
    end if
end subroutine f90wrap_grid_diff_map_type__array__local_map1d

subroutine f90wrap_grid_diff_map_type__array__boundary(this, nd, dtype, dshape, dloc)
    use smpi_math_module, only: grid_diff_map_type
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

subroutine f90wrap_grid_diff_map_type__array__boundary1d(this, nd, dtype, dshape, dloc)
    use smpi_math_module, only: grid_diff_map_type
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
    dshape(1:1) = shape(this_ptr%p%boundary1d)
    dloc = loc(this_ptr%p%boundary1d)
end subroutine f90wrap_grid_diff_map_type__array__boundary1d

subroutine f90wrap_grid_diff_map_type__array__rcount(this, nd, dtype, dshape, dloc)
    use smpi_math_module, only: grid_diff_map_type
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
    if (allocated(this_ptr%p%rcount)) then
        dshape(1:1) = shape(this_ptr%p%rcount)
        dloc = loc(this_ptr%p%rcount)
    else
        dloc = 0
    end if
end subroutine f90wrap_grid_diff_map_type__array__rcount

subroutine f90wrap_grid_diff_map_type__array__rdispls(this, nd, dtype, dshape, dloc)
    use smpi_math_module, only: grid_diff_map_type
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
    if (allocated(this_ptr%p%rdispls)) then
        dshape(1:1) = shape(this_ptr%p%rdispls)
        dloc = loc(this_ptr%p%rdispls)
    else
        dloc = 0
    end if
end subroutine f90wrap_grid_diff_map_type__array__rdispls

subroutine f90wrap_grid_diff_map_type__array__scount(this, nd, dtype, dshape, dloc)
    use smpi_math_module, only: grid_diff_map_type
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
    if (allocated(this_ptr%p%scount)) then
        dshape(1:1) = shape(this_ptr%p%scount)
        dloc = loc(this_ptr%p%scount)
    else
        dloc = 0
    end if
end subroutine f90wrap_grid_diff_map_type__array__scount

subroutine f90wrap_grid_diff_map_type__array__sdispls(this, nd, dtype, dshape, dloc)
    use smpi_math_module, only: grid_diff_map_type
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
    if (allocated(this_ptr%p%sdispls)) then
        dshape(1:1) = shape(this_ptr%p%sdispls)
        dloc = loc(this_ptr%p%sdispls)
    else
        dloc = 0
    end if
end subroutine f90wrap_grid_diff_map_type__array__sdispls

subroutine f90wrap_grid_diff_map_type_initialise(this)
    use smpi_math_module, only: grid_diff_map_type
    implicit none

    type grid_diff_map_type_ptr_type
        type(grid_diff_map_type), pointer :: p => NULL()
    end type grid_diff_map_type_ptr_type
    type(grid_diff_map_type_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate (this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_grid_diff_map_type_initialise

subroutine f90wrap_grid_diff_map_type_finalise(this)
    use smpi_math_module, only: grid_diff_map_type
    implicit none

    type grid_diff_map_type_ptr_type
        type(grid_diff_map_type), pointer :: p => NULL()
    end type grid_diff_map_type_ptr_type
    type(grid_diff_map_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate (this_ptr%p)
end subroutine f90wrap_grid_diff_map_type_finalise

subroutine f90wrap_sphere_type__get__Length(this, f90wrap_Length)
    use smpi_math_module, only: sphere_type
    implicit none
    type sphere_type_ptr_type
        type(sphere_type), pointer :: p => NULL()
    end type sphere_type_ptr_type
    integer, intent(in)   :: this(2)
    type(sphere_type_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_Length

    this_ptr = transfer(this, this_ptr)
    f90wrap_Length = this_ptr%p%Length
end subroutine f90wrap_sphere_type__get__Length

subroutine f90wrap_sphere_type__set__Length(this, f90wrap_Length)
    use smpi_math_module, only: sphere_type
    implicit none
    type sphere_type_ptr_type
        type(sphere_type), pointer :: p => NULL()
    end type sphere_type_ptr_type
    integer, intent(in)   :: this(2)
    type(sphere_type_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_Length

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%Length = f90wrap_Length
end subroutine f90wrap_sphere_type__set__Length

subroutine f90wrap_sphere_type__array__map3d(this, nd, dtype, dshape, dloc)
    use smpi_math_module, only: sphere_type
    implicit none
    type sphere_type_ptr_type
        type(sphere_type), pointer :: p => NULL()
    end type sphere_type_ptr_type
    integer, intent(in) :: this(2)
    type(sphere_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 2
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%map3d)) then
        dshape(1:2) = shape(this_ptr%p%map3d)
        dloc = loc(this_ptr%p%map3d)
    else
        dloc = 0
    end if
end subroutine f90wrap_sphere_type__array__map3d

subroutine f90wrap_sphere_type_initialise(this)
    use smpi_math_module, only: sphere_type
    implicit none

    type sphere_type_ptr_type
        type(sphere_type), pointer :: p => NULL()
    end type sphere_type_ptr_type
    type(sphere_type_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate (this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_sphere_type_initialise

subroutine f90wrap_sphere_type_finalise(this)
    use smpi_math_module, only: sphere_type
    implicit none

    type sphere_type_ptr_type
        type(sphere_type), pointer :: p => NULL()
    end type sphere_type_ptr_type
    type(sphere_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate (this_ptr%p)
end subroutine f90wrap_sphere_type_finalise

subroutine f90wrap_smpi_init
    use smpi_math_module, only: smpi_init
    implicit none

    call smpi_init()
end subroutine f90wrap_smpi_init

subroutine f90wrap_smpi_init_pbc
    use smpi_math_module, only: smpi_init_pbc
    implicit none

    call smpi_init_pbc()
end subroutine f90wrap_smpi_init_pbc

subroutine f90wrap_smpi_exit
    use smpi_math_module, only: smpi_exit
    implicit none

    call smpi_exit()
end subroutine f90wrap_smpi_exit

subroutine f90wrap_smpi_stop(message)
    use smpi_math_module, only: smpi_stop
    implicit none

    character*(*) :: message
    call smpi_stop(message=message)
end subroutine f90wrap_smpi_stop

subroutine f90wrap_smpi_stop_info(message)
    use smpi_math_module, only: smpi_stop_info
    implicit none

    character*(*) :: message
    call smpi_stop_info(message=message)
end subroutine f90wrap_smpi_stop_info

subroutine f90wrap_nstates_split(m, np)
    use smpi_math_module, only: nstates_split
    implicit none

    integer(4), intent(inout) :: m
    integer(4) :: np
    call nstates_split(m=m, np=np)
end subroutine f90wrap_nstates_split

subroutine f90wrap_nstates_split_2(m, np)
    use smpi_math_module, only: nstates_split_2
    implicit none

    integer(4), intent(inout) :: m
    integer(4) :: np
    call nstates_split_2(m=m, np=np)
end subroutine f90wrap_nstates_split_2

subroutine f90wrap_smpi_reduce_sum_real(amat, na, ramat, n0, n1)
    use smpi_math_module, only: smpi_reduce_sum_real
    implicit none

    real(8), dimension(n0) :: amat
    integer(4) :: na
    real(8), optional, dimension(n1) :: ramat
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(ramat) :: n1 = shape(ramat,0)
    call smpi_reduce_sum_real(amat=amat, na=na, ramat=ramat)
end subroutine f90wrap_smpi_reduce_sum_real

subroutine f90wrap_start_time(inlabel, flag, tic)
    use smpi_math_module, only: start_time
    implicit none

    character*(*) :: inlabel
    logical :: flag
    real(8), optional :: tic
    call start_time(inlabel=inlabel, flag=flag, tic=tic)
end subroutine f90wrap_start_time

subroutine f90wrap_end_time(inlabel, flag, toc)
    use smpi_math_module, only: end_time
    implicit none

    character*(*) :: inlabel
    logical :: flag
    real(8), optional :: toc
    call end_time(inlabel=inlabel, flag=flag, toc=toc)
end subroutine f90wrap_end_time

subroutine f90wrap_write_time(inlabel, flag)
    use smpi_math_module, only: write_time
    implicit none

    character*(*) :: inlabel
    logical :: flag
    call write_time(inlabel=inlabel, flag=flag)
end subroutine f90wrap_write_time

subroutine f90wrap_write_sum_time(inlabel, flag)
    use smpi_math_module, only: write_sum_time
    implicit none

    character*(*) :: inlabel
    logical :: flag
    call write_sum_time(inlabel=inlabel, flag=flag)
end subroutine f90wrap_write_sum_time

subroutine f90wrap_print_time(inlabel, t)
    use smpi_math_module, only: print_time
    implicit none

    character*(*) :: inlabel
    real(8) :: t
    call print_time(inlabel=inlabel, t=t)
end subroutine f90wrap_print_time

subroutine f90wrap_states_split(nev)
    use smpi_math_module, only: states_split
    implicit none

    integer(4), intent(inout) :: nev
    call states_split(nev=nev)
end subroutine f90wrap_states_split

subroutine f90wrap_array_split(nev)
    use smpi_math_module, only: array_split
    implicit none

    integer(4), intent(in) :: nev
    call array_split(nev=nev)
end subroutine f90wrap_array_split

subroutine f90wrap_grid_split(ngrid, ncore, comm, id, grid_range, recvcounts, displs, gridrange_sum, n1, n2, n3, n, n0)
    use smpi_math_module, only: grid_split
    implicit none

    integer(4), intent(in) :: ngrid
    integer(4), intent(in) :: ncore
    integer(4), intent(in) :: comm
    integer(4), intent(in) :: id
    integer(4), dimension(3), intent(inout) :: grid_range
    integer(4), intent(inout), dimension(n0) :: recvcounts
    integer(4), intent(inout), dimension(n1) :: displs
    integer(4), optional, intent(inout), dimension(3, n2) :: gridrange_sum
    integer(4), intent(inout), optional :: n1
    integer(4), intent(inout), optional :: n2
    integer(4), intent(inout), optional :: n3
    integer(4), intent(inout), optional :: n
    integer :: n0
    !f2py intent(hide), depend(recvcounts) :: n0 = shape(recvcounts,0)
    call grid_split(ngrid=ngrid, ncore=ncore, comm=comm, id=id, grid_range=grid_range, recvcounts=recvcounts, displs=displs, &
                    gridrange_sum=gridrange_sum, n1=n1, n2=n2, n3=n3, n=n)
end subroutine f90wrap_grid_split

subroutine f90wrap_atom_split(mysize, natom, atom_index, n0)
    use smpi_math_module, only: atom_split
    implicit none

    integer(4) :: mysize
    integer(4) :: natom
    integer(4), dimension(n0) :: atom_index
    integer :: n0
    !f2py intent(hide), depend(atom_index) :: n0 = shape(atom_index,0)
    call atom_split(mysize=mysize, natom=natom, atom_index=atom_index)
end subroutine f90wrap_atom_split

subroutine f90wrap_grid_sphere_init(n1, n2, n3, norder)
    use smpi_math_module, only: grid_sphere_init
    implicit none

    integer(4) :: n1
    integer(4) :: n2
    integer(4) :: n3
    integer(4) :: norder
    call grid_sphere_init(n1=n1, n2=n2, n3=n3, norder=norder)
end subroutine f90wrap_grid_sphere_init

subroutine f90wrap_set_wrap_grid_iso(myrho, wrap_box, n0, n1, n2, n3)
    use smpi_math_module, only: set_wrap_grid_iso
    implicit none

    real(8), dimension(n0) :: myrho
    real(8), dimension(n1, n2, n3) :: wrap_box
    integer :: n0
    !f2py intent(hide), depend(myrho) :: n0 = shape(myrho,0)
    integer :: n1
    !f2py intent(hide), depend(wrap_box) :: n1 = shape(wrap_box,0)
    integer :: n2
    !f2py intent(hide), depend(wrap_box) :: n2 = shape(wrap_box,1)
    integer :: n3
    !f2py intent(hide), depend(wrap_box) :: n3 = shape(wrap_box,2)
    call set_wrap_grid_iso(myrho=myrho, wrap_box=wrap_box)
end subroutine f90wrap_set_wrap_grid_iso

subroutine f90wrap_destroy_diff_map
    use smpi_math_module, only: destroy_diff_map
    implicit none

    call destroy_diff_map()
end subroutine f90wrap_destroy_diff_map

subroutine f90wrap_smpi_diff_init_sph(n1, n2, n3, n, norder, cell_mu, lsp, n0)
    use smpi_math_module, only: smpi_diff_init_sph
    implicit none

    integer(4) :: n1
    integer(4) :: n2
    integer(4) :: n3
    integer(4) :: n
    integer(4) :: norder
    integer(4), dimension(3, 3) :: cell_mu
    logical, dimension(n0, n1, n2) :: lsp
    integer :: n0
    !f2py intent(hide), depend(lsp) :: n0 = shape(lsp,0)
    !f2py intent(hide), depend(lsp) :: n2 = shape(lsp,2)
    call smpi_diff_init_sph(n1=n1, n2=n2, n3=n3, n=n, norder=norder, cell_mu=cell_mu, lsp=lsp)
end subroutine f90wrap_smpi_diff_init_sph

subroutine f90wrap_smpi_diff_init(n1, n2, n3, n, norder, cell_mu)
    use smpi_math_module, only: smpi_diff_init
    implicit none

    integer(4) :: n1
    integer(4) :: n2
    integer(4) :: n3
    integer(4) :: n
    integer(4) :: norder
    integer(4), dimension(3, 3) :: cell_mu
    call smpi_diff_init(n1=n1, n2=n2, n3=n3, n=n, norder=norder, cell_mu=cell_mu)
end subroutine f90wrap_smpi_diff_init

subroutine f90wrap_set_wrap_sph_pbc_ata_real(myrho, wrap_box1d, n0, n1)
    use smpi_math_module, only: set_wrap_sph_pbc_ata_real
    implicit none

    real(8), intent(in), dimension(n0) :: myrho
    real(8), intent(inout), dimension(n1) :: wrap_box1d
    integer :: n0
    !f2py intent(hide), depend(myrho) :: n0 = shape(myrho,0)
    integer :: n1
    !f2py intent(hide), depend(wrap_box1d) :: n1 = shape(wrap_box1d,0)
    call set_wrap_sph_pbc_ata_real(myrho=myrho, wrap_box1d=wrap_box1d)
end subroutine f90wrap_set_wrap_sph_pbc_ata_real

subroutine f90wrap_set_fft_alltoallv(fft_grid_range_temp, n0)
    use smpi_math_module, only: set_fft_alltoallv
    implicit none

    integer(4), dimension(n0) :: fft_grid_range_temp
    integer :: n0
    !f2py intent(hide), depend(fft_grid_range_temp) :: n0 = shape(fft_grid_range_temp,0)
    call set_fft_alltoallv(fft_grid_range_temp=fft_grid_range_temp)
end subroutine f90wrap_set_fft_alltoallv

subroutine f90wrap_destroy_fft_alltoallv
    use smpi_math_module, only: destroy_fft_alltoallv
    implicit none

    call destroy_fft_alltoallv()
end subroutine f90wrap_destroy_fft_alltoallv

subroutine f90wrap_sum_real_1d(ret_totals, amat, n0, n1, n2)
    use smpi_math_module, only: sompsum
    implicit none

    real(8), intent(out) :: ret_totals
    real(8), dimension(n0, n1, n2) :: amat
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(amat) :: n1 = shape(amat,1)
    integer :: n2
    !f2py intent(hide), depend(amat) :: n2 = shape(amat,2)
    ret_totals = sompsum(amat=amat)
end subroutine f90wrap_sum_real_1d

subroutine f90wrap_sum_real_2d(amat, ret_totals, bmat, n0, n1, n2, n3, n4, n5)
    use smpi_math_module, only: sompsum
    implicit none

    real(8), dimension(n0, n1, n2) :: amat
    real(8), intent(out) :: ret_totals
    real(8), dimension(n3, n4, n5) :: bmat
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(amat) :: n1 = shape(amat,1)
    integer :: n2
    !f2py intent(hide), depend(amat) :: n2 = shape(amat,2)
    integer :: n3
    !f2py intent(hide), depend(bmat) :: n3 = shape(bmat,0)
    integer :: n4
    !f2py intent(hide), depend(bmat) :: n4 = shape(bmat,1)
    integer :: n5
    !f2py intent(hide), depend(bmat) :: n5 = shape(bmat,2)
    ret_totals = sompsum(amat=amat, bmat=bmat)
end subroutine f90wrap_sum_real_2d

subroutine f90wrap_sum_real_3d(amat, bmat, ret_totals, cmat, n0, n1, n2, n3, n4, n5, n6, n7, n8)
    use smpi_math_module, only: sompsum
    implicit none

    real(8), dimension(n0, n1, n2) :: amat
    real(8), dimension(n3, n4, n5) :: bmat
    real(8), intent(out) :: ret_totals
    real(8), dimension(n6, n7, n8) :: cmat
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(amat) :: n1 = shape(amat,1)
    integer :: n2
    !f2py intent(hide), depend(amat) :: n2 = shape(amat,2)
    integer :: n3
    !f2py intent(hide), depend(bmat) :: n3 = shape(bmat,0)
    integer :: n4
    !f2py intent(hide), depend(bmat) :: n4 = shape(bmat,1)
    integer :: n5
    !f2py intent(hide), depend(bmat) :: n5 = shape(bmat,2)
    integer :: n6
    !f2py intent(hide), depend(cmat) :: n6 = shape(cmat,0)
    integer :: n7
    !f2py intent(hide), depend(cmat) :: n7 = shape(cmat,1)
    integer :: n8
    !f2py intent(hide), depend(cmat) :: n8 = shape(cmat,2)
    ret_totals = sompsum(amat=amat, bmat=bmat, cmat=cmat)
end subroutine f90wrap_sum_real_3d

subroutine f90wrap_sum_cplx_1d(ret_totals, amat, n0, n1, n2)
    use smpi_math_module, only: sompsum
    implicit none

    complex(8), intent(out) :: ret_totals
    complex(8), dimension(n0, n1, n2) :: amat
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(amat) :: n1 = shape(amat,1)
    integer :: n2
    !f2py intent(hide), depend(amat) :: n2 = shape(amat,2)
    ret_totals = sompsum(amat=amat)
end subroutine f90wrap_sum_cplx_1d

subroutine f90wrap_sum_cplx_2d(amat, ret_totals, bmat, n0, n1, n2, n3, n4, n5)
    use smpi_math_module, only: sompsum
    implicit none

    complex(8), dimension(n0, n1, n2) :: amat
    complex(8), intent(out) :: ret_totals
    complex(8), dimension(n3, n4, n5) :: bmat
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(amat) :: n1 = shape(amat,1)
    integer :: n2
    !f2py intent(hide), depend(amat) :: n2 = shape(amat,2)
    integer :: n3
    !f2py intent(hide), depend(bmat) :: n3 = shape(bmat,0)
    integer :: n4
    !f2py intent(hide), depend(bmat) :: n4 = shape(bmat,1)
    integer :: n5
    !f2py intent(hide), depend(bmat) :: n5 = shape(bmat,2)
    ret_totals = sompsum(amat=amat, bmat=bmat)
end subroutine f90wrap_sum_cplx_2d

subroutine f90wrap_sum_cplx_3d(amat, bmat, ret_totals, cmat, n0, n1, n2, n3, n4, n5, n6, n7, n8)
    use smpi_math_module, only: sompsum
    implicit none

    complex(8), dimension(n0, n1, n2) :: amat
    complex(8), dimension(n3, n4, n5) :: bmat
    complex(8), intent(out) :: ret_totals
    complex(8), dimension(n6, n7, n8) :: cmat
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(amat) :: n1 = shape(amat,1)
    integer :: n2
    !f2py intent(hide), depend(amat) :: n2 = shape(amat,2)
    integer :: n3
    !f2py intent(hide), depend(bmat) :: n3 = shape(bmat,0)
    integer :: n4
    !f2py intent(hide), depend(bmat) :: n4 = shape(bmat,1)
    integer :: n5
    !f2py intent(hide), depend(bmat) :: n5 = shape(bmat,2)
    integer :: n6
    !f2py intent(hide), depend(cmat) :: n6 = shape(cmat,0)
    integer :: n7
    !f2py intent(hide), depend(cmat) :: n7 = shape(cmat,1)
    integer :: n8
    !f2py intent(hide), depend(cmat) :: n8 = shape(cmat,2)
    ret_totals = sompsum(amat=amat, bmat=bmat, cmat=cmat)
end subroutine f90wrap_sum_cplx_3d

subroutine f90wrap_smpi_sum_int_1s(ret_sumx, x)
    use smpi_math_module, only: smpisum
    implicit none

    integer(4), intent(out) :: ret_sumx
    integer(4) :: x
    ret_sumx = smpisum(x=x)
end subroutine f90wrap_smpi_sum_int_1s

subroutine f90wrap_smpi_sum_cplx_1s(ret_sumx, x)
    use smpi_math_module, only: smpisum
    implicit none

    complex(8), intent(out) :: ret_sumx
    complex(8) :: x
    ret_sumx = smpisum(x=x)
end subroutine f90wrap_smpi_sum_cplx_1s

subroutine f90wrap_smpi_sum_real_1s(ret_sumx, x)
    use smpi_math_module, only: smpisum
    implicit none

    real(8), intent(out) :: ret_sumx
    real(8) :: x
    ret_sumx = smpisum(x=x)
end subroutine f90wrap_smpi_sum_real_1s

subroutine f90wrap_smpi_sum_real_1d(ret_suma, amat, n0, n1, n2)
    use smpi_math_module, only: smpisum
    implicit none

    real(8), intent(out) :: ret_suma
    real(8), dimension(n0, n1, n2) :: amat
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(amat) :: n1 = shape(amat,1)
    integer :: n2
    !f2py intent(hide), depend(amat) :: n2 = shape(amat,2)
    ret_suma = smpisum(amat=amat)
end subroutine f90wrap_smpi_sum_real_1d

subroutine f90wrap_smpi_sum_real_2d(amat, ret_suma, bmat, n0, n1, n2, n3, n4, n5)
    use smpi_math_module, only: smpisum
    implicit none

    real(8), dimension(n0, n1, n2) :: amat
    real(8), intent(out) :: ret_suma
    real(8), dimension(n3, n4, n5) :: bmat
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(amat) :: n1 = shape(amat,1)
    integer :: n2
    !f2py intent(hide), depend(amat) :: n2 = shape(amat,2)
    integer :: n3
    !f2py intent(hide), depend(bmat) :: n3 = shape(bmat,0)
    integer :: n4
    !f2py intent(hide), depend(bmat) :: n4 = shape(bmat,1)
    integer :: n5
    !f2py intent(hide), depend(bmat) :: n5 = shape(bmat,2)
    ret_suma = smpisum(amat=amat, bmat=bmat)
end subroutine f90wrap_smpi_sum_real_2d

subroutine f90wrap_smpi_sum_real_3d(amat, bmat, ret_suma, cmat, n0, n1, n2, n3, n4, n5, n6, n7, n8)
    use smpi_math_module, only: smpisum
    implicit none

    real(8), dimension(n0, n1, n2) :: amat
    real(8), dimension(n3, n4, n5) :: bmat
    real(8), intent(out) :: ret_suma
    real(8), dimension(n6, n7, n8) :: cmat
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(amat) :: n1 = shape(amat,1)
    integer :: n2
    !f2py intent(hide), depend(amat) :: n2 = shape(amat,2)
    integer :: n3
    !f2py intent(hide), depend(bmat) :: n3 = shape(bmat,0)
    integer :: n4
    !f2py intent(hide), depend(bmat) :: n4 = shape(bmat,1)
    integer :: n5
    !f2py intent(hide), depend(bmat) :: n5 = shape(bmat,2)
    integer :: n6
    !f2py intent(hide), depend(cmat) :: n6 = shape(cmat,0)
    integer :: n7
    !f2py intent(hide), depend(cmat) :: n7 = shape(cmat,1)
    integer :: n8
    !f2py intent(hide), depend(cmat) :: n8 = shape(cmat,2)
    ret_suma = smpisum(amat=amat, bmat=bmat, cmat=cmat)
end subroutine f90wrap_smpi_sum_real_3d

subroutine f90wrap_smpi_sum_mem_1d(munit, ret_summem, amat, n0)
    use smpi_math_module, only: smpisummem
    implicit none

    character(8) :: munit
    real(8), intent(out) :: ret_summem
    real(8), dimension(n0) :: amat
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    ret_summem = smpisummem(munit=munit, amat=amat)
end subroutine f90wrap_smpi_sum_mem_1d

subroutine f90wrap_smpi_sum_mem_2d(munit, ret_summem, amat, n0, n1)
    use smpi_math_module, only: smpisummem
    implicit none

    character(8) :: munit
    real(8), intent(out) :: ret_summem
    real(8), dimension(n0, n1) :: amat
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(amat) :: n1 = shape(amat,1)
    ret_summem = smpisummem(munit=munit, amat=amat)
end subroutine f90wrap_smpi_sum_mem_2d

subroutine f90wrap_smpi_sum_mem_3d(munit, ret_summem, amat, n0, n1, n2)
    use smpi_math_module, only: smpisummem
    implicit none

    character(8) :: munit
    real(8), intent(out) :: ret_summem
    real(8), dimension(n0, n1, n2) :: amat
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(amat) :: n1 = shape(amat,1)
    integer :: n2
    !f2py intent(hide), depend(amat) :: n2 = shape(amat,2)
    ret_summem = smpisummem(munit=munit, amat=amat)
end subroutine f90wrap_smpi_sum_mem_3d

subroutine f90wrap_smpi_reduce_sum_real_1d(amat, ramat, n0, n1)
    use smpi_math_module, only: smpireducesum
    implicit none

    real(8), dimension(n0) :: amat
    real(8), optional, dimension(n1) :: ramat
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(ramat) :: n1 = shape(ramat,0)
    call smpireducesum(amat=amat, ramat=ramat)
end subroutine f90wrap_smpi_reduce_sum_real_1d

subroutine f90wrap_smpi_reduce_sum_int_1d(amat, ramat, n0, n1)
    use smpi_math_module, only: smpireducesum
    implicit none

    integer(4), dimension(n0) :: amat
    integer(4), optional, dimension(n1) :: ramat
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(ramat) :: n1 = shape(ramat,0)
    call smpireducesum(amat=amat, ramat=ramat)
end subroutine f90wrap_smpi_reduce_sum_int_1d

subroutine f90wrap_smpi_reduce_sum_cplx_1d(amat, ramat, n0, n1)
    use smpi_math_module, only: smpireducesum
    implicit none

    complex(8), dimension(n0) :: amat
    complex(8), optional, dimension(n1) :: ramat
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(ramat) :: n1 = shape(ramat,0)
    call smpireducesum(amat=amat, ramat=ramat)
end subroutine f90wrap_smpi_reduce_sum_cplx_1d

subroutine f90wrap_smpi_reduce_sum_real_2d(amat, ramat, n0, n1, n2, n3)
    use smpi_math_module, only: smpireducesum
    implicit none

    real(8), dimension(n0, n1) :: amat
    real(8), optional, dimension(n2, n3) :: ramat
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(amat) :: n1 = shape(amat,1)
    integer :: n2
    !f2py intent(hide), depend(ramat) :: n2 = shape(ramat,0)
    integer :: n3
    !f2py intent(hide), depend(ramat) :: n3 = shape(ramat,1)
    call smpireducesum(amat=amat, ramat=ramat)
end subroutine f90wrap_smpi_reduce_sum_real_2d

subroutine f90wrap_sum_pow_int(amat, ret_totals, pow, n0, n1, n2)
    use smpi_math_module, only: sompsumpow
    implicit none

    real(8), dimension(n0, n1, n2) :: amat
    real(8), intent(out) :: ret_totals
    integer(4) :: pow
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(amat) :: n1 = shape(amat,1)
    integer :: n2
    !f2py intent(hide), depend(amat) :: n2 = shape(amat,2)
    ret_totals = sompsumpow(amat=amat, pow=pow)
end subroutine f90wrap_sum_pow_int

subroutine f90wrap_sum_pow_real(amat, ret_totals, pow, n0, n1, n2)
    use smpi_math_module, only: sompsumpow
    implicit none

    real(8), dimension(n0, n1, n2) :: amat
    real(8), intent(out) :: ret_totals
    real(8) :: pow
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(amat) :: n1 = shape(amat,1)
    integer :: n2
    !f2py intent(hide), depend(amat) :: n2 = shape(amat,2)
    ret_totals = sompsumpow(amat=amat, pow=pow)
end subroutine f90wrap_sum_pow_real

subroutine f90wrap_smpi_sum_pow_int(amat, ret_suma, pow, n0, n1, n2)
    use smpi_math_module, only: smpisumpow
    implicit none

    real(8), dimension(n0, n1, n2) :: amat
    real(8), intent(out) :: ret_suma
    integer(4) :: pow
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(amat) :: n1 = shape(amat,1)
    integer :: n2
    !f2py intent(hide), depend(amat) :: n2 = shape(amat,2)
    ret_suma = smpisumpow(amat=amat, pow=pow)
end subroutine f90wrap_smpi_sum_pow_int

subroutine f90wrap_smpi_sum_pow_real(amat, ret_suma, pow, n0, n1, n2)
    use smpi_math_module, only: smpisumpow
    implicit none

    real(8), dimension(n0, n1, n2) :: amat
    real(8), intent(out) :: ret_suma
    real(8) :: pow
    integer :: n0
    !f2py intent(hide), depend(amat) :: n0 = shape(amat,0)
    integer :: n1
    !f2py intent(hide), depend(amat) :: n1 = shape(amat,1)
    integer :: n2
    !f2py intent(hide), depend(amat) :: n2 = shape(amat,2)
    ret_suma = smpisumpow(amat=amat, pow=pow)
end subroutine f90wrap_smpi_sum_pow_real

subroutine f90wrap_set_wrap_grid_pbc_ata_cmplx(myrho, wrap_box1d, global_n, global_n1, global_n2, n0, n1)
    use smpi_math_module, only: set_wrap_grid_pbc_ata
    implicit none

    complex(8), dimension(n0) :: myrho
    complex(8), dimension(n1) :: wrap_box1d
    integer(4) :: global_n
    integer(4) :: global_n1
    integer(4) :: global_n2
    integer :: n0
    !f2py intent(hide), depend(myrho) :: n0 = shape(myrho,0)
    integer :: n1
    !f2py intent(hide), depend(wrap_box1d) :: n1 = shape(wrap_box1d,0)
    call set_wrap_grid_pbc_ata(myrho=myrho, wrap_box1d=wrap_box1d, global_n=global_n, global_n1=global_n1, &
                               global_n2=global_n2)
end subroutine f90wrap_set_wrap_grid_pbc_ata_cmplx

subroutine f90wrap_set_wrap_grid_pbc_ata_real(myrho, wrap_box1d, global_n, global_n1, global_n2, n0, n1)
    use smpi_math_module, only: set_wrap_grid_pbc_ata
    implicit none

    real(8), intent(in), dimension(n0) :: myrho
    real(8), intent(inout), dimension(n1) :: wrap_box1d
    integer(4) :: global_n
    integer(4) :: global_n1
    integer(4) :: global_n2
    integer :: n0
    !f2py intent(hide), depend(myrho) :: n0 = shape(myrho,0)
    integer :: n1
    !f2py intent(hide), depend(wrap_box1d) :: n1 = shape(wrap_box1d,0)
    call set_wrap_grid_pbc_ata(myrho=myrho, wrap_box1d=wrap_box1d, global_n=global_n, global_n1=global_n1, &
                               global_n2=global_n2)
end subroutine f90wrap_set_wrap_grid_pbc_ata_real

subroutine f90wrap_smpi_math_module__get__rtic(f90wrap_rtic)
    use smpi_math_module, only: smpi_math_module_rtic => rtic
    implicit none
    real(8), intent(out) :: f90wrap_rtic

    f90wrap_rtic = smpi_math_module_rtic
end subroutine f90wrap_smpi_math_module__get__rtic

subroutine f90wrap_smpi_math_module__set__rtic(f90wrap_rtic)
    use smpi_math_module, only: smpi_math_module_rtic => rtic
    implicit none
    real(8), intent(in) :: f90wrap_rtic

    smpi_math_module_rtic = f90wrap_rtic
end subroutine f90wrap_smpi_math_module__set__rtic

subroutine f90wrap_smpi_math_module__get__rtoc(f90wrap_rtoc)
    use smpi_math_module, only: smpi_math_module_rtoc => rtoc
    implicit none
    real(8), intent(out) :: f90wrap_rtoc

    f90wrap_rtoc = smpi_math_module_rtoc
end subroutine f90wrap_smpi_math_module__get__rtoc

subroutine f90wrap_smpi_math_module__set__rtoc(f90wrap_rtoc)
    use smpi_math_module, only: smpi_math_module_rtoc => rtoc
    implicit none
    real(8), intent(in) :: f90wrap_rtoc

    smpi_math_module_rtoc = f90wrap_rtoc
end subroutine f90wrap_smpi_math_module__set__rtoc

subroutine f90wrap_smpi_math_module__get__mpinfo(f90wrap_mpinfo)
    use smpi_math_module, only: smpi_math_module_mpinfo => mpinfo
    implicit none
    integer(4), intent(out) :: f90wrap_mpinfo

    f90wrap_mpinfo = smpi_math_module_mpinfo
end subroutine f90wrap_smpi_math_module__get__mpinfo

subroutine f90wrap_smpi_math_module__set__mpinfo(f90wrap_mpinfo)
    use smpi_math_module, only: smpi_math_module_mpinfo => mpinfo
    implicit none
    integer(4), intent(in) :: f90wrap_mpinfo

    smpi_math_module_mpinfo = f90wrap_mpinfo
end subroutine f90wrap_smpi_math_module__set__mpinfo

subroutine f90wrap_smpi_math_module__array__smpi_status(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use smpi_math_module, only: smpi_math_module_smpi_status => smpi_status
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 5
    dshape(1:1) = shape(smpi_math_module_smpi_status)
    dloc = loc(smpi_math_module_smpi_status)
end subroutine f90wrap_smpi_math_module__array__smpi_status

subroutine f90wrap_smpi_math_module__get__Lall_grid(f90wrap_Lall_grid)
    use smpi_math_module, only: smpi_math_module_Lall_grid => Lall_grid
    implicit none
    logical, intent(out) :: f90wrap_Lall_grid

    f90wrap_Lall_grid = smpi_math_module_Lall_grid
end subroutine f90wrap_smpi_math_module__get__Lall_grid

subroutine f90wrap_smpi_math_module__set__Lall_grid(f90wrap_Lall_grid)
    use smpi_math_module, only: smpi_math_module_Lall_grid => Lall_grid
    implicit none
    logical, intent(in) :: f90wrap_Lall_grid

    smpi_math_module_Lall_grid = f90wrap_Lall_grid
end subroutine f90wrap_smpi_math_module__set__Lall_grid

! End of module smpi_math_module defined in file Smpi_math_module.fpp


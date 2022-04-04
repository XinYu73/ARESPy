! Module m_time_evaluate defined in file MPI_time_evaluate.f90

subroutine f90wrap_time_record__array__str(this, nd, dtype, dshape, dloc)
    use m_time_evaluate, only: time_record
    implicit none
    type time_record_ptr_type
        type(time_record), pointer :: p => NULL()
    end type time_record_ptr_type
    integer, intent(in) :: this(2)
    type(time_record_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 2
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%str)) then
        dshape(1:2) = (/len(this_ptr%p%str(1)), shape(this_ptr%p%str)/)
        dloc = loc(this_ptr%p%str)
    else
        dloc = 0
    end if
end subroutine f90wrap_time_record__array__str

subroutine f90wrap_time_record__array__t1(this, nd, dtype, dshape, dloc)
    use m_time_evaluate, only: time_record
    implicit none
    type time_record_ptr_type
        type(time_record), pointer :: p => NULL()
    end type time_record_ptr_type
    integer, intent(in) :: this(2)
    type(time_record_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 7
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%t1)) then
        dshape(1:1) = shape(this_ptr%p%t1)
        dloc = loc(this_ptr%p%t1)
    else
        dloc = 0
    end if
end subroutine f90wrap_time_record__array__t1

subroutine f90wrap_time_record__array__t2(this, nd, dtype, dshape, dloc)
    use m_time_evaluate, only: time_record
    implicit none
    type time_record_ptr_type
        type(time_record), pointer :: p => NULL()
    end type time_record_ptr_type
    integer, intent(in) :: this(2)
    type(time_record_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 7
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%t2)) then
        dshape(1:1) = shape(this_ptr%p%t2)
        dloc = loc(this_ptr%p%t2)
    else
        dloc = 0
    end if
end subroutine f90wrap_time_record__array__t2

subroutine f90wrap_time_record__get__length(this, f90wrap_length)
    use m_time_evaluate, only: time_record
    implicit none
    type time_record_ptr_type
        type(time_record), pointer :: p => NULL()
    end type time_record_ptr_type
    integer, intent(in)   :: this(2)
    type(time_record_ptr_type) :: this_ptr
    integer(8), intent(out) :: f90wrap_length
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_length = this_ptr%p%length
end subroutine f90wrap_time_record__get__length

subroutine f90wrap_time_record__set__length(this, f90wrap_length)
    use m_time_evaluate, only: time_record
    implicit none
    type time_record_ptr_type
        type(time_record), pointer :: p => NULL()
    end type time_record_ptr_type
    integer, intent(in)   :: this(2)
    type(time_record_ptr_type) :: this_ptr
    integer(8), intent(in) :: f90wrap_length
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%length = f90wrap_length
end subroutine f90wrap_time_record__set__length

subroutine f90wrap_time_record_initialise(this)
    use m_time_evaluate, only: time_record
    implicit none
    
    type time_record_ptr_type
        type(time_record), pointer :: p => NULL()
    end type time_record_ptr_type
    type(time_record_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_time_record_initialise

subroutine f90wrap_time_record_finalise(this)
    use m_time_evaluate, only: time_record
    implicit none
    
    type time_record_ptr_type
        type(time_record), pointer :: p => NULL()
    end type time_record_ptr_type
    type(time_record_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_time_record_finalise

subroutine f90wrap_record(str, t)
    use m_time_evaluate, only: record
    implicit none
    
    character(100), intent(in) :: str
    integer(8), intent(in) :: t
    call record(str=str, t=t)
end subroutine f90wrap_record

subroutine f90wrap_reset_data(myflag, str, t)
    use m_time_evaluate, only: reset_data
    implicit none
    
    logical, intent(in) :: myflag
    character(100) :: str
    integer(8) :: t
    call reset_data(myflag=myflag, str=str, t=t)
end subroutine f90wrap_reset_data

subroutine f90wrap_init_recrod
    use m_time_evaluate, only: init_recrod
    implicit none
    
    call init_recrod()
end subroutine f90wrap_init_recrod

subroutine f90wrap_time_start(str, flag)
    use m_time_evaluate, only: time_start
    implicit none
    
    character*(*) :: str
    logical :: flag
    call time_start(str=str, flag=flag)
end subroutine f90wrap_time_start

subroutine f90wrap_time_end(str, flag)
    use m_time_evaluate, only: time_end
    implicit none
    
    character*(*) :: str
    logical :: flag
    call time_end(str=str, flag=flag)
end subroutine f90wrap_time_end

subroutine f90wrap_time_output(str, flag)
    use m_time_evaluate, only: time_output
    implicit none
    
    character*(*) :: str
    logical :: flag
    call time_output(str=str, flag=flag)
end subroutine f90wrap_time_output

subroutine f90wrap_memory_sum(label, memory_size)
    use m_time_evaluate, only: memory_sum
    implicit none
    
    character*(*) :: label
    real(8) :: memory_size
    call memory_sum(label=label, memory_size=memory_size)
end subroutine f90wrap_memory_sum

subroutine f90wrap_memory_free(label, memory_size)
    use m_time_evaluate, only: memory_free
    implicit none
    
    character*(*) :: label
    real(8) :: memory_size
    call memory_free(label=label, memory_size=memory_size)
end subroutine f90wrap_memory_free

subroutine f90wrap_m_time_evaluate__get__total_memory_consum(f90wrap_total_memory_consum)
    use m_time_evaluate, only: m_time_evaluate_total_memory_consum => total_memory_consum
    implicit none
    real(8), intent(out) :: f90wrap_total_memory_consum
    
    f90wrap_total_memory_consum = m_time_evaluate_total_memory_consum
end subroutine f90wrap_m_time_evaluate__get__total_memory_consum

subroutine f90wrap_m_time_evaluate__set__total_memory_consum(f90wrap_total_memory_consum)
    use m_time_evaluate, only: m_time_evaluate_total_memory_consum => total_memory_consum
    implicit none
    real(8), intent(in) :: f90wrap_total_memory_consum
    
    m_time_evaluate_total_memory_consum = f90wrap_total_memory_consum
end subroutine f90wrap_m_time_evaluate__set__total_memory_consum

subroutine f90wrap_m_time_evaluate__get__filename(f90wrap_filename)
    use m_time_evaluate, only: m_time_evaluate_filename => filename
    implicit none
    character(20), intent(out) :: f90wrap_filename
    
    f90wrap_filename = m_time_evaluate_filename
end subroutine f90wrap_m_time_evaluate__get__filename

subroutine f90wrap_m_time_evaluate__set__filename(f90wrap_filename)
    use m_time_evaluate, only: m_time_evaluate_filename => filename
    implicit none
    character(20), intent(in) :: f90wrap_filename
    
    m_time_evaluate_filename = f90wrap_filename
end subroutine f90wrap_m_time_evaluate__set__filename

subroutine f90wrap_m_time_evaluate__get__file_unit(f90wrap_file_unit)
    use m_time_evaluate, only: m_time_evaluate_file_unit => file_unit
    implicit none
    integer(8), intent(out) :: f90wrap_file_unit
    
    f90wrap_file_unit = m_time_evaluate_file_unit
end subroutine f90wrap_m_time_evaluate__get__file_unit

subroutine f90wrap_m_time_evaluate__set__file_unit(f90wrap_file_unit)
    use m_time_evaluate, only: m_time_evaluate_file_unit => file_unit
    implicit none
    integer(8), intent(in) :: f90wrap_file_unit
    
    m_time_evaluate_file_unit = f90wrap_file_unit
end subroutine f90wrap_m_time_evaluate__set__file_unit

! End of module m_time_evaluate defined in file MPI_time_evaluate.f90


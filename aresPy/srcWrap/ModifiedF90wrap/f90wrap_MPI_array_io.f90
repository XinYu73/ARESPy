! Module array_io defined in file MPI_array_io.fpp

subroutine f90wrap_out_label__get__label(this, f90wrap_label)
    use array_io, only: out_label
    implicit none
    type out_label_ptr_type
        type(out_label), pointer :: p => NULL()
    end type out_label_ptr_type
    integer, intent(in)   :: this(2)
    type(out_label_ptr_type) :: this_ptr
    character(100), intent(out) :: f90wrap_label
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_label = this_ptr%p%label
end subroutine f90wrap_out_label__get__label

subroutine f90wrap_out_label__set__label(this, f90wrap_label)
    use array_io, only: out_label
    implicit none
    type out_label_ptr_type
        type(out_label), pointer :: p => NULL()
    end type out_label_ptr_type
    integer, intent(in)   :: this(2)
    type(out_label_ptr_type) :: this_ptr
    character(100), intent(in) :: f90wrap_label
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%label = f90wrap_label
end subroutine f90wrap_out_label__set__label

subroutine f90wrap_out_label__get__num(this, f90wrap_num)
    use array_io, only: out_label
    implicit none
    type out_label_ptr_type
        type(out_label), pointer :: p => NULL()
    end type out_label_ptr_type
    integer, intent(in)   :: this(2)
    type(out_label_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_num
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_num = this_ptr%p%num
end subroutine f90wrap_out_label__get__num

subroutine f90wrap_out_label__set__num(this, f90wrap_num)
    use array_io, only: out_label
    implicit none
    type out_label_ptr_type
        type(out_label), pointer :: p => NULL()
    end type out_label_ptr_type
    integer, intent(in)   :: this(2)
    type(out_label_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_num
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%num = f90wrap_num
end subroutine f90wrap_out_label__set__num

subroutine f90wrap_out_label_initialise(this)
    use array_io, only: out_label
    implicit none
    
    type out_label_ptr_type
        type(out_label), pointer :: p => NULL()
    end type out_label_ptr_type
    type(out_label_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_out_label_initialise

subroutine f90wrap_out_label_finalise(this)
    use array_io, only: out_label
    implicit none
    
    type out_label_ptr_type
        type(out_label), pointer :: p => NULL()
    end type out_label_ptr_type
    type(out_label_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_out_label_finalise

subroutine f90wrap_init_outfile
    use array_io, only: init_outfile
    implicit none
    
    call init_outfile()
end subroutine f90wrap_init_outfile

subroutine f90wrap_output_r(size_bn, array, name, n0)
    use array_io, only: output
    implicit none
    
    integer(4) :: size_bn
    real(8), intent(in), dimension(n0) :: array
    character*(*) :: name
    integer :: n0
    !f2py intent(hide), depend(array) :: n0 = shape(array,0)
    call output(size=size_bn, array=array, name=name)
end subroutine f90wrap_output_r

subroutine f90wrap_output_i(size_bn, array, name, n0)
    use array_io, only: output
    implicit none
    
    integer(4) :: size_bn
    integer(4), intent(in), dimension(n0) :: array
    character*(*) :: name
    integer :: n0
    !f2py intent(hide), depend(array) :: n0 = shape(array,0)
    call output(size=size_bn, array=array, name=name)
end subroutine f90wrap_output_i

subroutine f90wrap_input_r(size_bn, array, name, n0)
    use array_io, only: input
    implicit none
    
    integer(4) :: size_bn
    real(8), intent(inout), dimension(n0) :: array
    character*(*) :: name
    integer :: n0
    !f2py intent(hide), depend(array) :: n0 = shape(array,0)
    call input(size=size_bn, array=array, name=name)
end subroutine f90wrap_input_r

subroutine f90wrap_input_i(size_bn, array, name, n0)
    use array_io, only: input
    implicit none
    
    integer(4) :: size_bn
    integer(4), intent(inout), dimension(n0) :: array
    character*(*) :: name
    integer :: n0
    !f2py intent(hide), depend(array) :: n0 = shape(array,0)
    call input(size=size_bn, array=array, name=name)
end subroutine f90wrap_input_i

! End of module array_io defined in file MPI_array_io.fpp


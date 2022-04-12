! Module dyn_module defined in file dyn_module.fpp

subroutine f90wrap_dynamics__get__lstop(this, f90wrap_lstop)
    use dyn_module, only: dynamics
    implicit none
    type dynamics_ptr_type
        type(dynamics), pointer :: p => NULL()
    end type dynamics_ptr_type
    integer, intent(in)   :: this(2)
    type(dynamics_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_lstop
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_lstop = this_ptr%p%lstop
end subroutine f90wrap_dynamics__get__lstop

subroutine f90wrap_dynamics__set__lstop(this, f90wrap_lstop)
    use dyn_module, only: dynamics
    implicit none
    type dynamics_ptr_type
        type(dynamics), pointer :: p => NULL()
    end type dynamics_ptr_type
    integer, intent(in)   :: this(2)
    type(dynamics_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_lstop
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%lstop = f90wrap_lstop
end subroutine f90wrap_dynamics__set__lstop

subroutine f90wrap_dynamics__array__posion(this, nd, dtype, dshape, dloc)
    use dyn_module, only: dynamics
    implicit none
    type dynamics_ptr_type
        type(dynamics), pointer :: p => NULL()
    end type dynamics_ptr_type
    integer, intent(in) :: this(2)
    type(dynamics_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%posion)) then
        dshape(1:2) = shape(this_ptr%p%posion)
        dloc = loc(this_ptr%p%posion)
    else
        dloc = 0
    end if
end subroutine f90wrap_dynamics__array__posion

subroutine f90wrap_dynamics__array__posioc(this, nd, dtype, dshape, dloc)
    use dyn_module, only: dynamics
    implicit none
    type dynamics_ptr_type
        type(dynamics), pointer :: p => NULL()
    end type dynamics_ptr_type
    integer, intent(in) :: this(2)
    type(dynamics_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%posioc)) then
        dshape(1:2) = shape(this_ptr%p%posioc)
        dloc = loc(this_ptr%p%posioc)
    else
        dloc = 0
    end if
end subroutine f90wrap_dynamics__array__posioc

subroutine f90wrap_dynamics__array__d2(this, nd, dtype, dshape, dloc)
    use dyn_module, only: dynamics
    implicit none
    type dynamics_ptr_type
        type(dynamics), pointer :: p => NULL()
    end type dynamics_ptr_type
    integer, intent(in) :: this(2)
    type(dynamics_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%d2)) then
        dshape(1:2) = shape(this_ptr%p%d2)
        dloc = loc(this_ptr%p%d2)
    else
        dloc = 0
    end if
end subroutine f90wrap_dynamics__array__d2

subroutine f90wrap_dynamics__array__d2c(this, nd, dtype, dshape, dloc)
    use dyn_module, only: dynamics
    implicit none
    type dynamics_ptr_type
        type(dynamics), pointer :: p => NULL()
    end type dynamics_ptr_type
    integer, intent(in) :: this(2)
    type(dynamics_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%d2c)) then
        dshape(1:2) = shape(this_ptr%p%d2c)
        dloc = loc(this_ptr%p%d2c)
    else
        dloc = 0
    end if
end subroutine f90wrap_dynamics__array__d2c

subroutine f90wrap_dynamics__array__d3(this, nd, dtype, dshape, dloc)
    use dyn_module, only: dynamics
    implicit none
    type dynamics_ptr_type
        type(dynamics), pointer :: p => NULL()
    end type dynamics_ptr_type
    integer, intent(in) :: this(2)
    type(dynamics_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%d3)) then
        dshape(1:2) = shape(this_ptr%p%d3)
        dloc = loc(this_ptr%p%d3)
    else
        dloc = 0
    end if
end subroutine f90wrap_dynamics__array__d3

subroutine f90wrap_dynamics__array__A(this, nd, dtype, dshape, dloc)
    use dyn_module, only: dynamics
    implicit none
    type dynamics_ptr_type
        type(dynamics), pointer :: p => NULL()
    end type dynamics_ptr_type
    integer, intent(in) :: this(2)
    type(dynamics_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    dshape(1:2) = shape(this_ptr%p%A)
    dloc = loc(this_ptr%p%A)
end subroutine f90wrap_dynamics__array__A

subroutine f90wrap_dynamics__array__B(this, nd, dtype, dshape, dloc)
    use dyn_module, only: dynamics
    implicit none
    type dynamics_ptr_type
        type(dynamics), pointer :: p => NULL()
    end type dynamics_ptr_type
    integer, intent(in) :: this(2)
    type(dynamics_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    dshape(1:2) = shape(this_ptr%p%B)
    dloc = loc(this_ptr%p%B)
end subroutine f90wrap_dynamics__array__B

subroutine f90wrap_dynamics__array__AC(this, nd, dtype, dshape, dloc)
    use dyn_module, only: dynamics
    implicit none
    type dynamics_ptr_type
        type(dynamics), pointer :: p => NULL()
    end type dynamics_ptr_type
    integer, intent(in) :: this(2)
    type(dynamics_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    dshape(1:2) = shape(this_ptr%p%AC)
    dloc = loc(this_ptr%p%AC)
end subroutine f90wrap_dynamics__array__AC

subroutine f90wrap_dynamics__get__POTIM(this, f90wrap_POTIM)
    use dyn_module, only: dynamics
    implicit none
    type dynamics_ptr_type
        type(dynamics), pointer :: p => NULL()
    end type dynamics_ptr_type
    integer, intent(in)   :: this(2)
    type(dynamics_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_POTIM
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_POTIM = this_ptr%p%POTIM
end subroutine f90wrap_dynamics__get__POTIM

subroutine f90wrap_dynamics__set__POTIM(this, f90wrap_POTIM)
    use dyn_module, only: dynamics
    implicit none
    type dynamics_ptr_type
        type(dynamics), pointer :: p => NULL()
    end type dynamics_ptr_type
    integer, intent(in)   :: this(2)
    type(dynamics_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_POTIM
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%POTIM = f90wrap_POTIM
end subroutine f90wrap_dynamics__set__POTIM

subroutine f90wrap_dynamics__get__EDIFFG(this, f90wrap_EDIFFG)
    use dyn_module, only: dynamics
    implicit none
    type dynamics_ptr_type
        type(dynamics), pointer :: p => NULL()
    end type dynamics_ptr_type
    integer, intent(in)   :: this(2)
    type(dynamics_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_EDIFFG
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_EDIFFG = this_ptr%p%EDIFFG
end subroutine f90wrap_dynamics__get__EDIFFG

subroutine f90wrap_dynamics__set__EDIFFG(this, f90wrap_EDIFFG)
    use dyn_module, only: dynamics
    implicit none
    type dynamics_ptr_type
        type(dynamics), pointer :: p => NULL()
    end type dynamics_ptr_type
    integer, intent(in)   :: this(2)
    type(dynamics_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_EDIFFG
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%EDIFFG = f90wrap_EDIFFG
end subroutine f90wrap_dynamics__set__EDIFFG

subroutine f90wrap_dynamics__get__PSTRESS(this, f90wrap_PSTRESS)
    use dyn_module, only: dynamics
    implicit none
    type dynamics_ptr_type
        type(dynamics), pointer :: p => NULL()
    end type dynamics_ptr_type
    integer, intent(in)   :: this(2)
    type(dynamics_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_PSTRESS
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_PSTRESS = this_ptr%p%PSTRESS
end subroutine f90wrap_dynamics__get__PSTRESS

subroutine f90wrap_dynamics__set__PSTRESS(this, f90wrap_PSTRESS)
    use dyn_module, only: dynamics
    implicit none
    type dynamics_ptr_type
        type(dynamics), pointer :: p => NULL()
    end type dynamics_ptr_type
    integer, intent(in)   :: this(2)
    type(dynamics_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_PSTRESS
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%PSTRESS = f90wrap_PSTRESS
end subroutine f90wrap_dynamics__set__PSTRESS

subroutine f90wrap_dynamics__get__IBRION(this, f90wrap_IBRION)
    use dyn_module, only: dynamics
    implicit none
    type dynamics_ptr_type
        type(dynamics), pointer :: p => NULL()
    end type dynamics_ptr_type
    integer, intent(in)   :: this(2)
    type(dynamics_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_IBRION
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_IBRION = this_ptr%p%IBRION
end subroutine f90wrap_dynamics__get__IBRION

subroutine f90wrap_dynamics__set__IBRION(this, f90wrap_IBRION)
    use dyn_module, only: dynamics
    implicit none
    type dynamics_ptr_type
        type(dynamics), pointer :: p => NULL()
    end type dynamics_ptr_type
    integer, intent(in)   :: this(2)
    type(dynamics_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_IBRION
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%IBRION = f90wrap_IBRION
end subroutine f90wrap_dynamics__set__IBRION

subroutine f90wrap_dynamics__get__ISIF(this, f90wrap_ISIF)
    use dyn_module, only: dynamics
    implicit none
    type dynamics_ptr_type
        type(dynamics), pointer :: p => NULL()
    end type dynamics_ptr_type
    integer, intent(in)   :: this(2)
    type(dynamics_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_ISIF
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ISIF = this_ptr%p%ISIF
end subroutine f90wrap_dynamics__get__ISIF

subroutine f90wrap_dynamics__set__ISIF(this, f90wrap_ISIF)
    use dyn_module, only: dynamics
    implicit none
    type dynamics_ptr_type
        type(dynamics), pointer :: p => NULL()
    end type dynamics_ptr_type
    integer, intent(in)   :: this(2)
    type(dynamics_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_ISIF
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ISIF = f90wrap_ISIF
end subroutine f90wrap_dynamics__set__ISIF

subroutine f90wrap_dynamics__get__NFREE(this, f90wrap_NFREE)
    use dyn_module, only: dynamics
    implicit none
    type dynamics_ptr_type
        type(dynamics), pointer :: p => NULL()
    end type dynamics_ptr_type
    integer, intent(in)   :: this(2)
    type(dynamics_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_NFREE
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_NFREE = this_ptr%p%NFREE
end subroutine f90wrap_dynamics__get__NFREE

subroutine f90wrap_dynamics__set__NFREE(this, f90wrap_NFREE)
    use dyn_module, only: dynamics
    implicit none
    type dynamics_ptr_type
        type(dynamics), pointer :: p => NULL()
    end type dynamics_ptr_type
    integer, intent(in)   :: this(2)
    type(dynamics_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_NFREE
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%NFREE = f90wrap_NFREE
end subroutine f90wrap_dynamics__set__NFREE

subroutine f90wrap_dynamics_initialise(this)
    use dyn_module, only: dynamics
    implicit none
    
    type dynamics_ptr_type
        type(dynamics), pointer :: p => NULL()
    end type dynamics_ptr_type
    type(dynamics_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_dynamics_initialise

subroutine f90wrap_dynamics_finalise(this)
    use dyn_module, only: dynamics
    implicit none
    
    type dynamics_ptr_type
        type(dynamics), pointer :: p => NULL()
    end type dynamics_ptr_type
    type(dynamics_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_dynamics_finalise

subroutine f90wrap_create_dyn(na, dyn)
    use dyn_module, only: dynamics, create_dyn
    implicit none
    
    type dynamics_ptr_type
        type(dynamics), pointer :: p => NULL()
    end type dynamics_ptr_type
    integer(4), intent(in) :: na
    type(dynamics_ptr_type) :: dyn_ptr
    integer, intent(in), dimension(2) :: dyn
    dyn_ptr = transfer(dyn, dyn_ptr)
    call create_dyn(na=na, dyn=dyn_ptr%p)
end subroutine f90wrap_create_dyn

subroutine f90wrap_destroy_dyn(dyn)
    use dyn_module, only: dynamics, destroy_dyn
    implicit none
    
    type dynamics_ptr_type
        type(dynamics), pointer :: p => NULL()
    end type dynamics_ptr_type
    type(dynamics_ptr_type) :: dyn_ptr
    integer, intent(in), dimension(2) :: dyn
    dyn_ptr = transfer(dyn, dyn_ptr)
    call destroy_dyn(dyn=dyn_ptr%p)
end subroutine f90wrap_destroy_dyn

! End of module dyn_module defined in file dyn_module.fpp


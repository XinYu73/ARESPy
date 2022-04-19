! Module cg_relax defined in file cg_module.fpp

subroutine f90wrap_lattice__array__A(this, nd, dtype, dshape, dloc)
    use cg_relax, only: lattice
    implicit none
    type lattice_ptr_type
        type(lattice), pointer :: p => NULL()
    end type lattice_ptr_type
    integer, intent(in) :: this(2)
    type(lattice_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    dshape(1:2) = shape(this_ptr%p%A)
    dloc = loc(this_ptr%p%A)
end subroutine f90wrap_lattice__array__A

subroutine f90wrap_lattice__array__B(this, nd, dtype, dshape, dloc)
    use cg_relax, only: lattice
    implicit none
    type lattice_ptr_type
        type(lattice), pointer :: p => NULL()
    end type lattice_ptr_type
    integer, intent(in) :: this(2)
    type(lattice_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    dshape(1:2) = shape(this_ptr%p%B)
    dloc = loc(this_ptr%p%B)
end subroutine f90wrap_lattice__array__B

subroutine f90wrap_lattice__array__anorm(this, nd, dtype, dshape, dloc)
    use cg_relax, only: lattice
    implicit none
    type lattice_ptr_type
        type(lattice), pointer :: p => NULL()
    end type lattice_ptr_type
    integer, intent(in) :: this(2)
    type(lattice_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%anorm)
    dloc = loc(this_ptr%p%anorm)
end subroutine f90wrap_lattice__array__anorm

subroutine f90wrap_lattice__array__bnorm(this, nd, dtype, dshape, dloc)
    use cg_relax, only: lattice
    implicit none
    type lattice_ptr_type
        type(lattice), pointer :: p => NULL()
    end type lattice_ptr_type
    integer, intent(in) :: this(2)
    type(lattice_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%bnorm)
    dloc = loc(this_ptr%p%bnorm)
end subroutine f90wrap_lattice__array__bnorm

subroutine f90wrap_lattice__get__omega(this, f90wrap_omega)
    use cg_relax, only: lattice
    implicit none
    type lattice_ptr_type
        type(lattice), pointer :: p => NULL()
    end type lattice_ptr_type
    integer, intent(in)   :: this(2)
    type(lattice_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_omega
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_omega = this_ptr%p%omega
end subroutine f90wrap_lattice__get__omega

subroutine f90wrap_lattice__set__omega(this, f90wrap_omega)
    use cg_relax, only: lattice
    implicit none
    type lattice_ptr_type
        type(lattice), pointer :: p => NULL()
    end type lattice_ptr_type
    integer, intent(in)   :: this(2)
    type(lattice_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_omega
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%omega = f90wrap_omega
end subroutine f90wrap_lattice__set__omega

subroutine f90wrap_lattice_initialise(this)
    use cg_relax, only: lattice
    implicit none
    
    type lattice_ptr_type
        type(lattice), pointer :: p => NULL()
    end type lattice_ptr_type
    type(lattice_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_lattice_initialise

subroutine f90wrap_lattice_finalise(this)
    use cg_relax, only: lattice
    implicit none
    
    type lattice_ptr_type
        type(lattice), pointer :: p => NULL()
    end type lattice_ptr_type
    type(lattice_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_lattice_finalise

subroutine f90wrap_cg_relax_vasp_interface(nstep, lfopt, lopt, lexceed)
    use cg_relax, only: cg_relax_vasp_interface
    implicit none
    
    integer(4), intent(in) :: nstep
    logical, intent(inout) :: lfopt
    logical, intent(inout) :: lopt
    logical :: lexceed
    call cg_relax_vasp_interface(nstep=nstep, lfopt=lfopt, lopt=lopt, lexceed=lexceed)
end subroutine f90wrap_cg_relax_vasp_interface

subroutine f90wrap_check_opt(na, flag)
    use cg_relax, only: check_opt
    implicit none
    
    integer(4), intent(in) :: na
    logical :: flag
    call check_opt(na=na, flag=flag)
end subroutine f90wrap_check_opt

subroutine f90wrap_lattic(mylatt)
    use cg_relax, only: lattic, lattice
    implicit none
    
    type lattice_ptr_type
        type(lattice), pointer :: p => NULL()
    end type lattice_ptr_type
    type(lattice_ptr_type) :: mylatt_ptr
    integer, intent(in), dimension(2) :: mylatt
    mylatt_ptr = transfer(mylatt, mylatt_ptr)
    call lattic(mylatt=mylatt_ptr%p)
end subroutine f90wrap_lattic

subroutine f90wrap_expro(h, u1, u2)
    use cg_relax, only: expro
    implicit none
    
    real(8), dimension(3) :: h
    real(8), dimension(3) :: u1
    real(8), dimension(3) :: u2
    call expro(H=h, U1=u1, U2=u2)
end subroutine f90wrap_expro

subroutine f90wrap_check_distance(l, pos, lwrong, n0, n1)
    use cg_relax, only: check_distance
    implicit none
    
    real(8), intent(in) :: l
    real(8), intent(in), dimension(n0,n1) :: pos
    logical :: lwrong
    integer :: n0
    !f2py intent(hide), depend(pos) :: n0 = shape(pos,0)
    integer :: n1
    !f2py intent(hide), depend(pos) :: n1 = shape(pos,1)
    call check_distance(l=l, pos=pos, lwrong=lwrong)
end subroutine f90wrap_check_distance

subroutine f90wrap_cg_relax__get__MAXFORCE(f90wrap_MAXFORCE)
    use cg_relax, only: cg_relax_MAXFORCE => MAXFORCE
    implicit none
    real(8), intent(out) :: f90wrap_MAXFORCE
    
    f90wrap_MAXFORCE = cg_relax_MAXFORCE
end subroutine f90wrap_cg_relax__get__MAXFORCE

subroutine f90wrap_cg_relax__set__MAXFORCE(f90wrap_MAXFORCE)
    use cg_relax, only: cg_relax_MAXFORCE => MAXFORCE
    implicit none
    real(8), intent(in) :: f90wrap_MAXFORCE
    
    cg_relax_MAXFORCE = f90wrap_MAXFORCE
end subroutine f90wrap_cg_relax__set__MAXFORCE

subroutine f90wrap_cg_relax__get__MAXSTRESS(f90wrap_MAXSTRESS)
    use cg_relax, only: cg_relax_MAXSTRESS => MAXSTRESS
    implicit none
    real(8), intent(out) :: f90wrap_MAXSTRESS
    
    f90wrap_MAXSTRESS = cg_relax_MAXSTRESS
end subroutine f90wrap_cg_relax__get__MAXSTRESS

subroutine f90wrap_cg_relax__set__MAXSTRESS(f90wrap_MAXSTRESS)
    use cg_relax, only: cg_relax_MAXSTRESS => MAXSTRESS
    implicit none
    real(8), intent(in) :: f90wrap_MAXSTRESS
    
    cg_relax_MAXSTRESS = f90wrap_MAXSTRESS
end subroutine f90wrap_cg_relax__set__MAXSTRESS

subroutine f90wrap_cg_relax__array__TSTRESS(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use cg_relax, only: cg_relax_tstress => tstress
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    dshape(1:2) = shape(cg_relax_TSTRESS)
    dloc = loc(cg_relax_TSTRESS)
end subroutine f90wrap_cg_relax__array__TSTRESS

subroutine f90wrap_cg_relax__array__Latmark(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use cg_relax, only: cg_relax_latmark => latmark
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    dshape(1:2) = shape(cg_relax_Latmark)
    dloc = loc(cg_relax_Latmark)
end subroutine f90wrap_cg_relax__array__Latmark

! End of module cg_relax defined in file cg_module.fpp


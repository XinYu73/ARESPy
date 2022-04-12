! Module pspot_module defined in file Psp_module.fpp

subroutine f90wrap_pspot__get__Zion(this, f90wrap_Zion)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in)   :: this(2)
    type(pspot_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_Zion
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_Zion = this_ptr%p%Zion
end subroutine f90wrap_pspot__get__Zion

subroutine f90wrap_pspot__set__Zion(this, f90wrap_Zion)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in)   :: this(2)
    type(pspot_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_Zion
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%Zion = f90wrap_Zion
end subroutine f90wrap_pspot__set__Zion

subroutine f90wrap_pspot__get__numps(this, f90wrap_numps)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in)   :: this(2)
    type(pspot_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_numps
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_numps = this_ptr%p%numps
end subroutine f90wrap_pspot__get__numps

subroutine f90wrap_pspot__set__numps(this, f90wrap_numps)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in)   :: this(2)
    type(pspot_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_numps
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%numps = f90wrap_numps
end subroutine f90wrap_pspot__set__numps

subroutine f90wrap_pspot__get__qnumps(this, f90wrap_qnumps)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in)   :: this(2)
    type(pspot_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_qnumps
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_qnumps = this_ptr%p%qnumps
end subroutine f90wrap_pspot__get__qnumps

subroutine f90wrap_pspot__set__qnumps(this, f90wrap_qnumps)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in)   :: this(2)
    type(pspot_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_qnumps
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%qnumps = f90wrap_qnumps
end subroutine f90wrap_pspot__set__qnumps

subroutine f90wrap_pspot__get__numps_den(this, f90wrap_numps_den)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in)   :: this(2)
    type(pspot_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_numps_den
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_numps_den = this_ptr%p%numps_den
end subroutine f90wrap_pspot__get__numps_den

subroutine f90wrap_pspot__set__numps_den(this, f90wrap_numps_den)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in)   :: this(2)
    type(pspot_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_numps_den
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%numps_den = f90wrap_numps_den
end subroutine f90wrap_pspot__set__numps_den

subroutine f90wrap_pspot__get__qmax(this, f90wrap_qmax)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in)   :: this(2)
    type(pspot_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_qmax
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_qmax = this_ptr%p%qmax
end subroutine f90wrap_pspot__get__qmax

subroutine f90wrap_pspot__set__qmax(this, f90wrap_qmax)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in)   :: this(2)
    type(pspot_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_qmax
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%qmax = f90wrap_qmax
end subroutine f90wrap_pspot__set__qmax

subroutine f90wrap_pspot__get__qspacing(this, f90wrap_qspacing)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in)   :: this(2)
    type(pspot_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_qspacing
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_qspacing = this_ptr%p%qspacing
end subroutine f90wrap_pspot__get__qspacing

subroutine f90wrap_pspot__set__qspacing(this, f90wrap_qspacing)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in)   :: this(2)
    type(pspot_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_qspacing
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%qspacing = f90wrap_qspacing
end subroutine f90wrap_pspot__set__qspacing

subroutine f90wrap_pspot__array__qmesh(this, nd, dtype, dshape, dloc)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in) :: this(2)
    type(pspot_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%qmesh)) then
        dshape(1:1) = shape(this_ptr%p%qmesh)
        dloc = loc(this_ptr%p%qmesh)
    else
        dloc = 0
    end if
end subroutine f90wrap_pspot__array__qmesh

subroutine f90wrap_pspot__array__Vlocq(this, nd, dtype, dshape, dloc)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in) :: this(2)
    type(pspot_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%Vlocq)) then
        dshape(1:1) = shape(this_ptr%p%Vlocq)
        dloc = loc(this_ptr%p%Vlocq)
    else
        dloc = 0
    end if
end subroutine f90wrap_pspot__array__Vlocq

subroutine f90wrap_pspot__array__VlocqS(this, nd, dtype, dshape, dloc)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in) :: this(2)
    type(pspot_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%VlocqS)) then
        dshape(1:1) = shape(this_ptr%p%VlocqS)
        dloc = loc(this_ptr%p%VlocqS)
    else
        dloc = 0
    end if
end subroutine f90wrap_pspot__array__VlocqS

subroutine f90wrap_pspot__array__ddVl_dq2(this, nd, dtype, dshape, dloc)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in) :: this(2)
    type(pspot_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%ddVl_dq2)) then
        dshape(1:1) = shape(this_ptr%p%ddVl_dq2)
        dloc = loc(this_ptr%p%ddVl_dq2)
    else
        dloc = 0
    end if
end subroutine f90wrap_pspot__array__ddVl_dq2

subroutine f90wrap_pspot__get__nproj(this, f90wrap_nproj)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in)   :: this(2)
    type(pspot_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_nproj
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nproj = this_ptr%p%nproj
end subroutine f90wrap_pspot__get__nproj

subroutine f90wrap_pspot__set__nproj(this, f90wrap_nproj)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in)   :: this(2)
    type(pspot_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_nproj
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nproj = f90wrap_nproj
end subroutine f90wrap_pspot__set__nproj

subroutine f90wrap_pspot__array__proj_l(this, nd, dtype, dshape, dloc)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in) :: this(2)
    type(pspot_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%proj_l)) then
        dshape(1:1) = shape(this_ptr%p%proj_l)
        dloc = loc(this_ptr%p%proj_l)
    else
        dloc = 0
    end if
end subroutine f90wrap_pspot__array__proj_l

subroutine f90wrap_pspot__array__proj_m(this, nd, dtype, dshape, dloc)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in) :: this(2)
    type(pspot_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%proj_m)) then
        dshape(1:1) = shape(this_ptr%p%proj_m)
        dloc = loc(this_ptr%p%proj_m)
    else
        dloc = 0
    end if
end subroutine f90wrap_pspot__array__proj_m

subroutine f90wrap_pspot__get__rcut(this, f90wrap_rcut)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in)   :: this(2)
    type(pspot_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_rcut
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_rcut = this_ptr%p%rcut
end subroutine f90wrap_pspot__get__rcut

subroutine f90wrap_pspot__set__rcut(this, f90wrap_rcut)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in)   :: this(2)
    type(pspot_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_rcut
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%rcut = f90wrap_rcut
end subroutine f90wrap_pspot__set__rcut

subroutine f90wrap_pspot__get__rmax(this, f90wrap_rmax)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in)   :: this(2)
    type(pspot_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_rmax
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_rmax = this_ptr%p%rmax
end subroutine f90wrap_pspot__get__rmax

subroutine f90wrap_pspot__set__rmax(this, f90wrap_rmax)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in)   :: this(2)
    type(pspot_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_rmax
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%rmax = f90wrap_rmax
end subroutine f90wrap_pspot__set__rmax

subroutine f90wrap_pspot__get__rspacing(this, f90wrap_rspacing)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in)   :: this(2)
    type(pspot_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_rspacing
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_rspacing = this_ptr%p%rspacing
end subroutine f90wrap_pspot__get__rspacing

subroutine f90wrap_pspot__set__rspacing(this, f90wrap_rspacing)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in)   :: this(2)
    type(pspot_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_rspacing
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%rspacing = f90wrap_rspacing
end subroutine f90wrap_pspot__set__rspacing

subroutine f90wrap_pspot__array__D0(this, nd, dtype, dshape, dloc)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in) :: this(2)
    type(pspot_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%D0)) then
        dshape(1:2) = shape(this_ptr%p%D0)
        dloc = loc(this_ptr%p%D0)
    else
        dloc = 0
    end if
end subroutine f90wrap_pspot__array__D0

subroutine f90wrap_pspot__array__beta_r(this, nd, dtype, dshape, dloc)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in) :: this(2)
    type(pspot_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%beta_r)) then
        dshape(1:2) = shape(this_ptr%p%beta_r)
        dloc = loc(this_ptr%p%beta_r)
    else
        dloc = 0
    end if
end subroutine f90wrap_pspot__array__beta_r

subroutine f90wrap_pspot__array__dbeta_dr(this, nd, dtype, dshape, dloc)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in) :: this(2)
    type(pspot_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%dbeta_dr)) then
        dshape(1:2) = shape(this_ptr%p%dbeta_dr)
        dloc = loc(this_ptr%p%dbeta_dr)
    else
        dloc = 0
    end if
end subroutine f90wrap_pspot__array__dbeta_dr

subroutine f90wrap_pspot__array__denr(this, nd, dtype, dshape, dloc)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in) :: this(2)
    type(pspot_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%denr)) then
        dshape(1:1) = shape(this_ptr%p%denr)
        dloc = loc(this_ptr%p%denr)
    else
        dloc = 0
    end if
end subroutine f90wrap_pspot__array__denr

subroutine f90wrap_pspot__array__ddden_dr2(this, nd, dtype, dshape, dloc)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in) :: this(2)
    type(pspot_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%ddden_dr2)) then
        dshape(1:1) = shape(this_ptr%p%ddden_dr2)
        dloc = loc(this_ptr%p%ddden_dr2)
    else
        dloc = 0
    end if
end subroutine f90wrap_pspot__array__ddden_dr2

subroutine f90wrap_pspot__array__r_real(this, nd, dtype, dshape, dloc)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in) :: this(2)
    type(pspot_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%r_real)) then
        dshape(1:1) = shape(this_ptr%p%r_real)
        dloc = loc(this_ptr%p%r_real)
    else
        dloc = 0
    end if
end subroutine f90wrap_pspot__array__r_real

subroutine f90wrap_pspot__array__V_loc(this, nd, dtype, dshape, dloc)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in) :: this(2)
    type(pspot_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%V_loc)) then
        dshape(1:1) = shape(this_ptr%p%V_loc)
        dloc = loc(this_ptr%p%V_loc)
    else
        dloc = 0
    end if
end subroutine f90wrap_pspot__array__V_loc

subroutine f90wrap_pspot_initialise(this)
    use pspot_module, only: pspot
    implicit none
    
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    type(pspot_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_pspot_initialise

subroutine f90wrap_pspot_finalise(this)
    use pspot_module, only: pspot
    implicit none
    
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    type(pspot_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_pspot_finalise

subroutine f90wrap_pspot_module__get__max_nproj(f90wrap_max_nproj)
    use pspot_module, only: pspot_module_max_nproj => max_nproj
    implicit none
    integer(4), intent(out) :: f90wrap_max_nproj
    
    f90wrap_max_nproj = pspot_module_max_nproj
end subroutine f90wrap_pspot_module__get__max_nproj

subroutine f90wrap_pspot_module__set__max_nproj(f90wrap_max_nproj)
    use pspot_module, only: pspot_module_max_nproj => max_nproj
    implicit none
    integer(4), intent(in) :: f90wrap_max_nproj
    
    pspot_module_max_nproj = f90wrap_max_nproj
end subroutine f90wrap_pspot_module__set__max_nproj

subroutine f90wrap_pspot_module__get__max_rcut(f90wrap_max_rcut)
    use pspot_module, only: pspot_module_max_rcut => max_rcut
    implicit none
    real(8), intent(out) :: f90wrap_max_rcut
    
    f90wrap_max_rcut = pspot_module_max_rcut
end subroutine f90wrap_pspot_module__get__max_rcut

subroutine f90wrap_pspot_module__set__max_rcut(f90wrap_max_rcut)
    use pspot_module, only: pspot_module_max_rcut => max_rcut
    implicit none
    real(8), intent(in) :: f90wrap_max_rcut
    
    pspot_module_max_rcut = f90wrap_max_rcut
end subroutine f90wrap_pspot_module__set__max_rcut

subroutine f90wrap_pspot_module__array__tknots(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use pspot_module, only: pspot_module_tknots => tknots
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(pspot_module_tknots)) then
        dshape(1:1) = shape(pspot_module_tknots)
        dloc = loc(pspot_module_tknots)
    else
        dloc = 0
    end if
end subroutine f90wrap_pspot_module__array__tknots

! End of module pspot_module defined in file Psp_module.fpp


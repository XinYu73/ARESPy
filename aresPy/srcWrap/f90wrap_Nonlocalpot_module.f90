! Module nlpot_module defined in file Nonlocalpot_module.f90

subroutine f90wrap_nol_type__get__npts(this, f90wrap_npts)
    use nlpot_module, only: nol_type
    implicit none
    type nol_type_ptr_type
        type(nol_type), pointer :: p => NULL()
    end type nol_type_ptr_type
    integer, intent(in)   :: this(2)
    type(nol_type_ptr_type) :: this_ptr
    integer(8), intent(out) :: f90wrap_npts
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_npts = this_ptr%p%npts
end subroutine f90wrap_nol_type__get__npts

subroutine f90wrap_nol_type__set__npts(this, f90wrap_npts)
    use nlpot_module, only: nol_type
    implicit none
    type nol_type_ptr_type
        type(nol_type), pointer :: p => NULL()
    end type nol_type_ptr_type
    integer, intent(in)   :: this(2)
    type(nol_type_ptr_type) :: this_ptr
    integer(8), intent(in) :: f90wrap_npts
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%npts = f90wrap_npts
end subroutine f90wrap_nol_type__set__npts

subroutine f90wrap_nol_type__array__Id(this, nd, dtype, dshape, dloc)
    use nlpot_module, only: nol_type
    implicit none
    type nol_type_ptr_type
        type(nol_type), pointer :: p => NULL()
    end type nol_type_ptr_type
    integer, intent(in) :: this(2)
    type(nol_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 7
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%Id)) then
        dshape(1:1) = shape(this_ptr%p%Id)
        dloc = loc(this_ptr%p%Id)
    else
        dloc = 0
    end if
end subroutine f90wrap_nol_type__array__Id

subroutine f90wrap_nol_type__array__rRvec(this, nd, dtype, dshape, dloc)
    use nlpot_module, only: nol_type
    implicit none
    type nol_type_ptr_type
        type(nol_type), pointer :: p => NULL()
    end type nol_type_ptr_type
    integer, intent(in) :: this(2)
    type(nol_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%rRvec)) then
        dshape(1:2) = shape(this_ptr%p%rRvec)
        dloc = loc(this_ptr%p%rRvec)
    else
        dloc = 0
    end if
end subroutine f90wrap_nol_type__array__rRvec

subroutine f90wrap_nol_type__array__proj0(this, nd, dtype, dshape, dloc)
    use nlpot_module, only: nol_type
    implicit none
    type nol_type_ptr_type
        type(nol_type), pointer :: p => NULL()
    end type nol_type_ptr_type
    integer, intent(in) :: this(2)
    type(nol_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%proj0)) then
        dshape(1:2) = shape(this_ptr%p%proj0)
        dloc = loc(this_ptr%p%proj0)
    else
        dloc = 0
    end if
end subroutine f90wrap_nol_type__array__proj0

subroutine f90wrap_nol_type__array__proj(this, nd, dtype, dshape, dloc)
    use nlpot_module, only: nol_type
    implicit none
    type nol_type_ptr_type
        type(nol_type), pointer :: p => NULL()
    end type nol_type_ptr_type
    integer, intent(in) :: this(2)
    type(nol_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%proj)) then
        dshape(1:2) = shape(this_ptr%p%proj)
        dloc = loc(this_ptr%p%proj)
    else
        dloc = 0
    end if
end subroutine f90wrap_nol_type__array__proj

subroutine f90wrap_nol_type__array__proj_phs(this, nd, dtype, dshape, dloc)
    use nlpot_module, only: nol_type
    implicit none
    type nol_type_ptr_type
        type(nol_type), pointer :: p => NULL()
    end type nol_type_ptr_type
    integer, intent(in) :: this(2)
    type(nol_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 15
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%proj_phs)) then
        dshape(1:3) = shape(this_ptr%p%proj_phs)
        dloc = loc(this_ptr%p%proj_phs)
    else
        dloc = 0
    end if
end subroutine f90wrap_nol_type__array__proj_phs

subroutine f90wrap_nol_type__array__proj0_dg(this, nd, dtype, dshape, dloc)
    use nlpot_module, only: nol_type
    implicit none
    type nol_type_ptr_type
        type(nol_type), pointer :: p => NULL()
    end type nol_type_ptr_type
    integer, intent(in) :: this(2)
    type(nol_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%proj0_dg)) then
        dshape(1:2) = shape(this_ptr%p%proj0_dg)
        dloc = loc(this_ptr%p%proj0_dg)
    else
        dloc = 0
    end if
end subroutine f90wrap_nol_type__array__proj0_dg

subroutine f90wrap_nol_type__array__proj_dg(this, nd, dtype, dshape, dloc)
    use nlpot_module, only: nol_type
    implicit none
    type nol_type_ptr_type
        type(nol_type), pointer :: p => NULL()
    end type nol_type_ptr_type
    integer, intent(in) :: this(2)
    type(nol_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%proj_dg)) then
        dshape(1:2) = shape(this_ptr%p%proj_dg)
        dloc = loc(this_ptr%p%proj_dg)
    else
        dloc = 0
    end if
end subroutine f90wrap_nol_type__array__proj_dg

subroutine f90wrap_nol_type__array__proj_phs_dg(this, nd, dtype, dshape, dloc)
    use nlpot_module, only: nol_type
    implicit none
    type nol_type_ptr_type
        type(nol_type), pointer :: p => NULL()
    end type nol_type_ptr_type
    integer, intent(in) :: this(2)
    type(nol_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 15
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%proj_phs_dg)) then
        dshape(1:3) = shape(this_ptr%p%proj_phs_dg)
        dloc = loc(this_ptr%p%proj_phs_dg)
    else
        dloc = 0
    end if
end subroutine f90wrap_nol_type__array__proj_phs_dg

subroutine f90wrap_nol_type_initialise(this)
    use nlpot_module, only: nol_type
    implicit none
    
    type nol_type_ptr_type
        type(nol_type), pointer :: p => NULL()
    end type nol_type_ptr_type
    type(nol_type_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_nol_type_initialise

subroutine f90wrap_nol_type_finalise(this)
    use nlpot_module, only: nol_type
    implicit none
    
    type nol_type_ptr_type
        type(nol_type), pointer :: p => NULL()
    end type nol_type_ptr_type
    type(nol_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_nol_type_finalise

subroutine f90wrap_initialize_nlpot
    use nlpot_module, only: initialize_nlpot
    implicit none
    
    call initialize_nlpot()
end subroutine f90wrap_initialize_nlpot

subroutine f90wrap_destroy_nlpot
    use nlpot_module, only: destroy_nlpot
    implicit none
    
    call destroy_nlpot()
end subroutine f90wrap_destroy_nlpot

subroutine f90wrap_set_beta_real
    use nlpot_module, only: set_beta_real
    implicit none
    
    call set_beta_real()
end subroutine f90wrap_set_beta_real

subroutine f90wrap_nlp_beta_interp_r(ity, ia, beta_init, n0, n1)
    use nlpot_module, only: nlp_beta_interp_r
    implicit none
    
    integer(8), intent(in) :: ity
    integer(8), intent(in) :: ia
    real(8), intent(inout), dimension(n0,n1) :: beta_init
    integer :: n0
    !f2py intent(hide), depend(beta_init) :: n0 = shape(beta_init,0)
    integer :: n1
    !f2py intent(hide), depend(beta_init) :: n1 = shape(beta_init,1)
    call nlp_beta_interp_r(Ity=ity, Ia=ia, beta_init=beta_init)
end subroutine f90wrap_nlp_beta_interp_r

subroutine f90wrap_nlp_beta_ylm_r(ity, ia, beta_ylm, n0, n1)
    use nlpot_module, only: nlp_beta_ylm_r
    implicit none
    
    integer(8), intent(in) :: ity
    integer(8), intent(in) :: ia
    real(8), intent(inout), dimension(n0,n1) :: beta_ylm
    integer :: n0
    !f2py intent(hide), depend(beta_ylm) :: n0 = shape(beta_ylm,0)
    integer :: n1
    !f2py intent(hide), depend(beta_ylm) :: n1 = shape(beta_ylm,1)
    call nlp_beta_ylm_r(Ity=ity, Ia=ia, beta_Ylm=beta_ylm)
end subroutine f90wrap_nlp_beta_ylm_r

subroutine f90wrap_nlp_beta_phase_r(ity, ia, beta_ylm, beta_phase, n0, n1, n2, n3, n4)
    use nlpot_module, only: nlp_beta_phase_r
    implicit none
    
    integer(8), intent(in) :: ity
    integer(8), intent(in) :: ia
    real(8), intent(in), dimension(n0,n1) :: beta_ylm
    complex(8), intent(inout), dimension(n2,n3,n4) :: beta_phase
    integer :: n0
    !f2py intent(hide), depend(beta_ylm) :: n0 = shape(beta_ylm,0)
    integer :: n1
    !f2py intent(hide), depend(beta_ylm) :: n1 = shape(beta_ylm,1)
    integer :: n2
    !f2py intent(hide), depend(beta_phase) :: n2 = shape(beta_phase,0)
    integer :: n3
    !f2py intent(hide), depend(beta_phase) :: n3 = shape(beta_phase,1)
    integer :: n4
    !f2py intent(hide), depend(beta_phase) :: n4 = shape(beta_phase,2)
    call nlp_beta_phase_r(Ity=ity, Ia=ia, beta_ylm=beta_ylm, beta_phase=beta_phase)
end subroutine f90wrap_nlp_beta_phase_r

subroutine f90wrap_apply_ylm(l, m, fac, x, y, z, f)
    use nlpot_module, only: apply_ylm
    implicit none
    
    integer(8), intent(in) :: l
    integer(8), intent(in) :: m
    real(8), intent(in) :: fac
    real(8), intent(in) :: x
    real(8), intent(in) :: y
    real(8), intent(in) :: z
    real(8), intent(inout) :: f
    call apply_ylm(l=l, m=m, fac=fac, x=x, y=y, z=z, f=f)
end subroutine f90wrap_apply_ylm

subroutine f90wrap_nlp_init_partialcore(rhoc, n0)
    use nlpot_module, only: nlp_init_partialcore
    implicit none
    
    real(8), intent(inout), dimension(n0) :: rhoc
    integer :: n0
    !f2py intent(hide), depend(rhoc) :: n0 = shape(rhoc,0)
    call nlp_init_partialcore(rhoc=rhoc)
end subroutine f90wrap_nlp_init_partialcore

subroutine f90wrap_nlpot_module__get__max_nlnpts(f90wrap_max_nlnpts)
    use nlpot_module, only: nlpot_module_max_nlnpts => max_nlnpts
    implicit none
    integer(8), intent(out) :: f90wrap_max_nlnpts
    
    f90wrap_max_nlnpts = nlpot_module_max_nlnpts
end subroutine f90wrap_nlpot_module__get__max_nlnpts

subroutine f90wrap_nlpot_module__set__max_nlnpts(f90wrap_max_nlnpts)
    use nlpot_module, only: nlpot_module_max_nlnpts => max_nlnpts
    implicit none
    integer(8), intent(in) :: f90wrap_max_nlnpts
    
    nlpot_module_max_nlnpts = f90wrap_max_nlnpts
end subroutine f90wrap_nlpot_module__set__max_nlnpts

! End of module nlpot_module defined in file Nonlocalpot_module.f90


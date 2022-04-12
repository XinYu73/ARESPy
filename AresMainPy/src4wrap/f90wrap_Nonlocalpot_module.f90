! Module nlpot_module defined in file Nonlocalpot_module.fpp

subroutine f90wrap_nol_type__get__npts(this, f90wrap_npts)
    use nlpot_module, only: nol_type
    implicit none
    type nol_type_ptr_type
        type(nol_type), pointer :: p => NULL()
    end type nol_type_ptr_type
    integer, intent(in)   :: this(2)
    type(nol_type_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_npts
    
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
    integer(4), intent(in) :: f90wrap_npts
    
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
    dtype = 5
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

subroutine f90wrap_nol_type__array__Id_iso(this, nd, dtype, dshape, dloc)
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
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%Id_iso)) then
        dshape(1:1) = shape(this_ptr%p%Id_iso)
        dloc = loc(this_ptr%p%Id_iso)
    else
        dloc = 0
    end if
end subroutine f90wrap_nol_type__array__Id_iso

subroutine f90wrap_nol_type__array__rRvec_iso(this, nd, dtype, dshape, dloc)
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
    if (allocated(this_ptr%p%rRvec_iso)) then
        dshape(1:2) = shape(this_ptr%p%rRvec_iso)
        dloc = loc(this_ptr%p%rRvec_iso)
    else
        dloc = 0
    end if
end subroutine f90wrap_nol_type__array__rRvec_iso

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

subroutine f90wrap_initialize_nlpot_per
    use nlpot_module, only: initialize_nlpot_per
    implicit none
    
    call initialize_nlpot_per()
end subroutine f90wrap_initialize_nlpot_per

subroutine f90wrap_initialize_nlpot_iso
    use nlpot_module, only: initialize_nlpot_iso
    implicit none
    
    call initialize_nlpot_iso()
end subroutine f90wrap_initialize_nlpot_iso

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
    
    integer(4), intent(in) :: ity
    integer(4), intent(in) :: ia
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
    
    integer(4), intent(in) :: ity
    integer(4), intent(in) :: ia
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
    
    integer(4), intent(in) :: ity
    integer(4), intent(in) :: ia
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
    
    integer(4), intent(in) :: l
    integer(4), intent(in) :: m
    real(8), intent(in) :: fac
    real(8), intent(in) :: x
    real(8), intent(in) :: y
    real(8), intent(in) :: z
    real(8), intent(inout) :: f
    call apply_ylm(l=l, m=m, fac=fac, x=x, y=y, z=z, f=f)
end subroutine f90wrap_apply_ylm

subroutine f90wrap_nlp_wijk_dg(ndg, inpol, numdg, lft, rit, wijk, drijk, n0, n1)
    use nlpot_module, only: nlp_wijk_dg
    implicit none
    
    integer(4), intent(in) :: ndg
    integer(4), intent(in) :: inpol
    integer(4), intent(in) :: numdg
    integer(4), intent(in) :: lft
    integer(4), intent(in) :: rit
    real(8), intent(inout), dimension(n0) :: wijk
    real(8), intent(inout), dimension(3,n1) :: drijk
    integer :: n0
    !f2py intent(hide), depend(wijk) :: n0 = shape(wijk,0)
    integer :: n1
    !f2py intent(hide), depend(drijk) :: n1 = shape(drijk,1)
    call nlp_wijk_dg(ndg=ndg, inpol=inpol, numdg=numdg, lft=lft, rit=rit, wijk=wijk, drijk=drijk)
end subroutine f90wrap_nlp_wijk_dg

subroutine f90wrap_nlp_wijk_dg_old(ndg, n_near, numdg, lft, rit, wijk, drijk, n0, n1)
    use nlpot_module, only: nlp_wijk_dg_old
    implicit none
    
    integer(4), intent(in) :: ndg
    integer(4), intent(in) :: n_near
    integer(4), intent(in) :: numdg
    integer(4), intent(in) :: lft
    integer(4), intent(in) :: rit
    real(8), intent(inout), dimension(n0) :: wijk
    real(8), intent(inout), dimension(3,n1) :: drijk
    integer :: n0
    !f2py intent(hide), depend(wijk) :: n0 = shape(wijk,0)
    integer :: n1
    !f2py intent(hide), depend(drijk) :: n1 = shape(drijk,1)
    call nlp_wijk_dg_old(ndg=ndg, n_near=n_near, numdg=numdg, lft=lft, rit=rit, wijk=wijk, drijk=drijk)
end subroutine f90wrap_nlp_wijk_dg_old

subroutine f90wrap_nlp_interp_betalm_dg(ity, ia, ip, numdg, lft, rit, wijk, drijk, beta0_ijk, betalm_ijk, n0, n1, n2, &
    n3)
    use nlpot_module, only: nlp_interp_betalm_dg
    implicit none
    
    integer(4), intent(in) :: ity
    integer(4), intent(in) :: ia
    integer(4), intent(in) :: ip
    integer(4), intent(in) :: numdg
    integer(4), intent(in) :: lft
    integer(4), intent(in) :: rit
    real(8), intent(in), dimension(n0) :: wijk
    real(8), intent(in), dimension(3,n1) :: drijk
    real(8), intent(inout), dimension(n2) :: beta0_ijk
    real(8), intent(inout), dimension(n3) :: betalm_ijk
    integer :: n0
    !f2py intent(hide), depend(wijk) :: n0 = shape(wijk,0)
    integer :: n1
    !f2py intent(hide), depend(drijk) :: n1 = shape(drijk,1)
    integer :: n2
    !f2py intent(hide), depend(beta0_ijk) :: n2 = shape(beta0_ijk,0)
    integer :: n3
    !f2py intent(hide), depend(betalm_ijk) :: n3 = shape(betalm_ijk,0)
    call nlp_interp_betalm_dg(Ity=ity, Ia=ia, Ip=ip, numdg=numdg, lft=lft, rit=rit, wijk=wijk, drijk=drijk, &
        beta0_IJK=beta0_ijk, betalm_IJK=betalm_ijk)
end subroutine f90wrap_nlp_interp_betalm_dg

subroutine f90wrap_nlp_interp_betalm_dg_iso(ity, ia, ip, numdg, lft, rit, wijk, drijk, beta0_ijk, betalm_ijk, n0, n1, &
    n2, n3)
    use nlpot_module, only: nlp_interp_betalm_dg_iso
    implicit none
    
    integer(4), intent(in) :: ity
    integer(4), intent(in) :: ia
    integer(4), intent(in) :: ip
    integer(4), intent(in) :: numdg
    integer(4), intent(in) :: lft
    integer(4), intent(in) :: rit
    real(8), intent(in), dimension(n0) :: wijk
    real(8), intent(in), dimension(3,n1) :: drijk
    real(8), intent(inout), dimension(n2) :: beta0_ijk
    real(8), intent(inout), dimension(n3) :: betalm_ijk
    integer :: n0
    !f2py intent(hide), depend(wijk) :: n0 = shape(wijk,0)
    integer :: n1
    !f2py intent(hide), depend(drijk) :: n1 = shape(drijk,1)
    integer :: n2
    !f2py intent(hide), depend(beta0_ijk) :: n2 = shape(beta0_ijk,0)
    integer :: n3
    !f2py intent(hide), depend(betalm_ijk) :: n3 = shape(betalm_ijk,0)
    call nlp_interp_betalm_dg_iso(Ity=ity, Ia=ia, Ip=ip, numdg=numdg, lft=lft, rit=rit, wijk=wijk, drijk=drijk, &
        beta0_IJK=beta0_ijk, betalm_IJK=betalm_ijk)
end subroutine f90wrap_nlp_interp_betalm_dg_iso

subroutine f90wrap_nlp_beta_ylm_r_dg(ity, ia, numdg, lft, rit, wijk, drijk, beta0, beta_ylm, n0, n1, n2, n3, n4, n5, n6)
    use nlpot_module, only: nlp_beta_ylm_r_dg
    implicit none
    
    integer(4), intent(in) :: ity
    integer(4), intent(in) :: ia
    integer(4), intent(in) :: numdg
    integer(4), intent(in) :: lft
    integer(4), intent(in) :: rit
    real(8), intent(in), dimension(n0) :: wijk
    real(8), intent(in), dimension(n1,n2) :: drijk
    real(8), intent(inout), dimension(n3,n4) :: beta0
    real(8), intent(inout), dimension(n5,n6) :: beta_ylm
    integer :: n0
    !f2py intent(hide), depend(wijk) :: n0 = shape(wijk,0)
    integer :: n1
    !f2py intent(hide), depend(drijk) :: n1 = shape(drijk,0)
    integer :: n2
    !f2py intent(hide), depend(drijk) :: n2 = shape(drijk,1)
    integer :: n3
    !f2py intent(hide), depend(beta0) :: n3 = shape(beta0,0)
    integer :: n4
    !f2py intent(hide), depend(beta0) :: n4 = shape(beta0,1)
    integer :: n5
    !f2py intent(hide), depend(beta_ylm) :: n5 = shape(beta_ylm,0)
    integer :: n6
    !f2py intent(hide), depend(beta_ylm) :: n6 = shape(beta_ylm,1)
    call nlp_beta_ylm_r_dg(Ity=ity, Ia=ia, numdg=numdg, lft=lft, rit=rit, wijk=wijk, drijk=drijk, beta0=beta0, &
        beta_Ylm=beta_ylm)
end subroutine f90wrap_nlp_beta_ylm_r_dg

subroutine f90wrap_nlp_beta_ylm_r_dg_iso(ity, ia, numdg, lft, rit, wijk, drijk, beta0, beta_ylm, n0, n1, n2, n3, n4, n5, &
    n6)
    use nlpot_module, only: nlp_beta_ylm_r_dg_iso
    implicit none
    
    integer(4), intent(in) :: ity
    integer(4), intent(in) :: ia
    integer(4), intent(in) :: numdg
    integer(4), intent(in) :: lft
    integer(4), intent(in) :: rit
    real(8), intent(in), dimension(n0) :: wijk
    real(8), intent(in), dimension(n1,n2) :: drijk
    real(8), intent(inout), dimension(n3,n4) :: beta0
    real(8), intent(inout), dimension(n5,n6) :: beta_ylm
    integer :: n0
    !f2py intent(hide), depend(wijk) :: n0 = shape(wijk,0)
    integer :: n1
    !f2py intent(hide), depend(drijk) :: n1 = shape(drijk,0)
    integer :: n2
    !f2py intent(hide), depend(drijk) :: n2 = shape(drijk,1)
    integer :: n3
    !f2py intent(hide), depend(beta0) :: n3 = shape(beta0,0)
    integer :: n4
    !f2py intent(hide), depend(beta0) :: n4 = shape(beta0,1)
    integer :: n5
    !f2py intent(hide), depend(beta_ylm) :: n5 = shape(beta_ylm,0)
    integer :: n6
    !f2py intent(hide), depend(beta_ylm) :: n6 = shape(beta_ylm,1)
    call nlp_beta_ylm_r_dg_iso(Ity=ity, Ia=ia, numdg=numdg, lft=lft, rit=rit, wijk=wijk, drijk=drijk, beta0=beta0, &
        beta_Ylm=beta_ylm)
end subroutine f90wrap_nlp_beta_ylm_r_dg_iso

subroutine f90wrap_set_beta_real_dg
    use nlpot_module, only: set_beta_real_dg
    implicit none
    
    call set_beta_real_dg()
end subroutine f90wrap_set_beta_real_dg

subroutine f90wrap_initialize_nlpot_band
    use nlpot_module, only: initialize_nlpot_band
    implicit none
    
    call initialize_nlpot_band()
end subroutine f90wrap_initialize_nlpot_band

subroutine f90wrap_nlpot_module__get__max_nlnpts(f90wrap_max_nlnpts)
    use nlpot_module, only: nlpot_module_max_nlnpts => max_nlnpts
    implicit none
    integer(4), intent(out) :: f90wrap_max_nlnpts
    
    f90wrap_max_nlnpts = nlpot_module_max_nlnpts
end subroutine f90wrap_nlpot_module__get__max_nlnpts

subroutine f90wrap_nlpot_module__set__max_nlnpts(f90wrap_max_nlnpts)
    use nlpot_module, only: nlpot_module_max_nlnpts => max_nlnpts
    implicit none
    integer(4), intent(in) :: f90wrap_max_nlnpts
    
    nlpot_module_max_nlnpts = f90wrap_max_nlnpts
end subroutine f90wrap_nlpot_module__set__max_nlnpts

! End of module nlpot_module defined in file Nonlocalpot_module.fpp


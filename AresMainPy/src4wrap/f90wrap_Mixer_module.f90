! Module mixer_module defined in file Mixer_module.fpp

subroutine f90wrap_mixer_data__get__NHMIX(this, f90wrap_NHMIX)
    use mixer_module, only: mixer_data
    implicit none
    type mixer_data_ptr_type
        type(mixer_data), pointer :: p => NULL()
    end type mixer_data_ptr_type
    integer, intent(in)   :: this(2)
    type(mixer_data_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_NHMIX
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_NHMIX = this_ptr%p%NHMIX
end subroutine f90wrap_mixer_data__get__NHMIX

subroutine f90wrap_mixer_data__set__NHMIX(this, f90wrap_NHMIX)
    use mixer_module, only: mixer_data
    implicit none
    type mixer_data_ptr_type
        type(mixer_data), pointer :: p => NULL()
    end type mixer_data_ptr_type
    integer, intent(in)   :: this(2)
    type(mixer_data_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_NHMIX
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%NHMIX = f90wrap_NHMIX
end subroutine f90wrap_mixer_data__set__NHMIX

subroutine f90wrap_mixer_data__get__NHMIN(this, f90wrap_NHMIN)
    use mixer_module, only: mixer_data
    implicit none
    type mixer_data_ptr_type
        type(mixer_data), pointer :: p => NULL()
    end type mixer_data_ptr_type
    integer, intent(in)   :: this(2)
    type(mixer_data_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_NHMIN
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_NHMIN = this_ptr%p%NHMIN
end subroutine f90wrap_mixer_data__get__NHMIN

subroutine f90wrap_mixer_data__set__NHMIN(this, f90wrap_NHMIN)
    use mixer_module, only: mixer_data
    implicit none
    type mixer_data_ptr_type
        type(mixer_data), pointer :: p => NULL()
    end type mixer_data_ptr_type
    integer, intent(in)   :: this(2)
    type(mixer_data_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_NHMIN
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%NHMIN = f90wrap_NHMIN
end subroutine f90wrap_mixer_data__set__NHMIN

subroutine f90wrap_mixer_data__get__NAM(this, f90wrap_NAM)
    use mixer_module, only: mixer_data
    implicit none
    type mixer_data_ptr_type
        type(mixer_data), pointer :: p => NULL()
    end type mixer_data_ptr_type
    integer, intent(in)   :: this(2)
    type(mixer_data_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_NAM
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_NAM = this_ptr%p%NAM
end subroutine f90wrap_mixer_data__get__NAM

subroutine f90wrap_mixer_data__set__NAM(this, f90wrap_NAM)
    use mixer_module, only: mixer_data
    implicit none
    type mixer_data_ptr_type
        type(mixer_data), pointer :: p => NULL()
    end type mixer_data_ptr_type
    integer, intent(in)   :: this(2)
    type(mixer_data_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_NAM
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%NAM = f90wrap_NAM
end subroutine f90wrap_mixer_data__set__NAM

subroutine f90wrap_mixer_data__get__NITRA(this, f90wrap_NITRA)
    use mixer_module, only: mixer_data
    implicit none
    type mixer_data_ptr_type
        type(mixer_data), pointer :: p => NULL()
    end type mixer_data_ptr_type
    integer, intent(in)   :: this(2)
    type(mixer_data_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_NITRA
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_NITRA = this_ptr%p%NITRA
end subroutine f90wrap_mixer_data__get__NITRA

subroutine f90wrap_mixer_data__set__NITRA(this, f90wrap_NITRA)
    use mixer_module, only: mixer_data
    implicit none
    type mixer_data_ptr_type
        type(mixer_data), pointer :: p => NULL()
    end type mixer_data_ptr_type
    integer, intent(in)   :: this(2)
    type(mixer_data_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_NITRA
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%NITRA = f90wrap_NITRA
end subroutine f90wrap_mixer_data__set__NITRA

subroutine f90wrap_mixer_data__get__alpha(this, f90wrap_alpha)
    use mixer_module, only: mixer_data
    implicit none
    type mixer_data_ptr_type
        type(mixer_data), pointer :: p => NULL()
    end type mixer_data_ptr_type
    integer, intent(in)   :: this(2)
    type(mixer_data_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_alpha
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_alpha = this_ptr%p%alpha
end subroutine f90wrap_mixer_data__get__alpha

subroutine f90wrap_mixer_data__set__alpha(this, f90wrap_alpha)
    use mixer_module, only: mixer_data
    implicit none
    type mixer_data_ptr_type
        type(mixer_data), pointer :: p => NULL()
    end type mixer_data_ptr_type
    integer, intent(in)   :: this(2)
    type(mixer_data_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_alpha
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%alpha = f90wrap_alpha
end subroutine f90wrap_mixer_data__set__alpha

subroutine f90wrap_mixer_data__get__beta(this, f90wrap_beta)
    use mixer_module, only: mixer_data
    implicit none
    type mixer_data_ptr_type
        type(mixer_data), pointer :: p => NULL()
    end type mixer_data_ptr_type
    integer, intent(in)   :: this(2)
    type(mixer_data_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_beta
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_beta = this_ptr%p%beta
end subroutine f90wrap_mixer_data__get__beta

subroutine f90wrap_mixer_data__set__beta(this, f90wrap_beta)
    use mixer_module, only: mixer_data
    implicit none
    type mixer_data_ptr_type
        type(mixer_data), pointer :: p => NULL()
    end type mixer_data_ptr_type
    integer, intent(in)   :: this(2)
    type(mixer_data_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_beta
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%beta = f90wrap_beta
end subroutine f90wrap_mixer_data__set__beta

subroutine f90wrap_mixer_data__get__w0(this, f90wrap_w0)
    use mixer_module, only: mixer_data
    implicit none
    type mixer_data_ptr_type
        type(mixer_data), pointer :: p => NULL()
    end type mixer_data_ptr_type
    integer, intent(in)   :: this(2)
    type(mixer_data_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_w0
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_w0 = this_ptr%p%w0
end subroutine f90wrap_mixer_data__get__w0

subroutine f90wrap_mixer_data__set__w0(this, f90wrap_w0)
    use mixer_module, only: mixer_data
    implicit none
    type mixer_data_ptr_type
        type(mixer_data), pointer :: p => NULL()
    end type mixer_data_ptr_type
    integer, intent(in)   :: this(2)
    type(mixer_data_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_w0
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%w0 = f90wrap_w0
end subroutine f90wrap_mixer_data__set__w0

subroutine f90wrap_mixer_data__array__DXL(this, nd, dtype, dshape, dloc)
    use mixer_module, only: mixer_data
    implicit none
    type mixer_data_ptr_type
        type(mixer_data), pointer :: p => NULL()
    end type mixer_data_ptr_type
    integer, intent(in) :: this(2)
    type(mixer_data_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%DXL)) then
        dshape(1:2) = shape(this_ptr%p%DXL)
        dloc = loc(this_ptr%p%DXL)
    else
        dloc = 0
    end if
end subroutine f90wrap_mixer_data__array__DXL

subroutine f90wrap_mixer_data__array__DFL(this, nd, dtype, dshape, dloc)
    use mixer_module, only: mixer_data
    implicit none
    type mixer_data_ptr_type
        type(mixer_data), pointer :: p => NULL()
    end type mixer_data_ptr_type
    integer, intent(in) :: this(2)
    type(mixer_data_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%DFL)) then
        dshape(1:2) = shape(this_ptr%p%DFL)
        dloc = loc(this_ptr%p%DFL)
    else
        dloc = 0
    end if
end subroutine f90wrap_mixer_data__array__DFL

subroutine f90wrap_mixer_data__array__VOMA(this, nd, dtype, dshape, dloc)
    use mixer_module, only: mixer_data
    implicit none
    type mixer_data_ptr_type
        type(mixer_data), pointer :: p => NULL()
    end type mixer_data_ptr_type
    integer, intent(in) :: this(2)
    type(mixer_data_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%VOMA)) then
        dshape(1:2) = shape(this_ptr%p%VOMA)
        dloc = loc(this_ptr%p%VOMA)
    else
        dloc = 0
    end if
end subroutine f90wrap_mixer_data__array__VOMA

subroutine f90wrap_mixer_data__get__SP(this, f90wrap_SP)
    use mixer_module, only: mixer_data
    implicit none
    type mixer_data_ptr_type
        type(mixer_data), pointer :: p => NULL()
    end type mixer_data_ptr_type
    integer, intent(in)   :: this(2)
    type(mixer_data_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_SP
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_SP = this_ptr%p%SP
end subroutine f90wrap_mixer_data__get__SP

subroutine f90wrap_mixer_data__set__SP(this, f90wrap_SP)
    use mixer_module, only: mixer_data
    implicit none
    type mixer_data_ptr_type
        type(mixer_data), pointer :: p => NULL()
    end type mixer_data_ptr_type
    integer, intent(in)   :: this(2)
    type(mixer_data_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_SP
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%SP = f90wrap_SP
end subroutine f90wrap_mixer_data__set__SP

subroutine f90wrap_mixer_data__array__kerker(this, nd, dtype, dshape, dloc)
    use mixer_module, only: mixer_data
    implicit none
    type mixer_data_ptr_type
        type(mixer_data), pointer :: p => NULL()
    end type mixer_data_ptr_type
    integer, intent(in) :: this(2)
    type(mixer_data_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%kerker)) then
        dshape(1:1) = shape(this_ptr%p%kerker)
        dloc = loc(this_ptr%p%kerker)
    else
        dloc = 0
    end if
end subroutine f90wrap_mixer_data__array__kerker

subroutine f90wrap_mixer_data__array__DXGL(this, nd, dtype, dshape, dloc)
    use mixer_module, only: mixer_data
    implicit none
    type mixer_data_ptr_type
        type(mixer_data), pointer :: p => NULL()
    end type mixer_data_ptr_type
    integer, intent(in) :: this(2)
    type(mixer_data_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 15
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%DXGL)) then
        dshape(1:2) = shape(this_ptr%p%DXGL)
        dloc = loc(this_ptr%p%DXGL)
    else
        dloc = 0
    end if
end subroutine f90wrap_mixer_data__array__DXGL

subroutine f90wrap_mixer_data__array__DFGL(this, nd, dtype, dshape, dloc)
    use mixer_module, only: mixer_data
    implicit none
    type mixer_data_ptr_type
        type(mixer_data), pointer :: p => NULL()
    end type mixer_data_ptr_type
    integer, intent(in) :: this(2)
    type(mixer_data_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 15
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%DFGL)) then
        dshape(1:2) = shape(this_ptr%p%DFGL)
        dloc = loc(this_ptr%p%DFGL)
    else
        dloc = 0
    end if
end subroutine f90wrap_mixer_data__array__DFGL

subroutine f90wrap_mixer_data__get__BMIX(this, f90wrap_BMIX)
    use mixer_module, only: mixer_data
    implicit none
    type mixer_data_ptr_type
        type(mixer_data), pointer :: p => NULL()
    end type mixer_data_ptr_type
    integer, intent(in)   :: this(2)
    type(mixer_data_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_BMIX
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_BMIX = this_ptr%p%BMIX
end subroutine f90wrap_mixer_data__get__BMIX

subroutine f90wrap_mixer_data__set__BMIX(this, f90wrap_BMIX)
    use mixer_module, only: mixer_data
    implicit none
    type mixer_data_ptr_type
        type(mixer_data), pointer :: p => NULL()
    end type mixer_data_ptr_type
    integer, intent(in)   :: this(2)
    type(mixer_data_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_BMIX
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%BMIX = f90wrap_BMIX
end subroutine f90wrap_mixer_data__set__BMIX

subroutine f90wrap_mixer_data__get__AMIX(this, f90wrap_AMIX)
    use mixer_module, only: mixer_data
    implicit none
    type mixer_data_ptr_type
        type(mixer_data), pointer :: p => NULL()
    end type mixer_data_ptr_type
    integer, intent(in)   :: this(2)
    type(mixer_data_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_AMIX
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_AMIX = this_ptr%p%AMIX
end subroutine f90wrap_mixer_data__get__AMIX

subroutine f90wrap_mixer_data__set__AMIX(this, f90wrap_AMIX)
    use mixer_module, only: mixer_data
    implicit none
    type mixer_data_ptr_type
        type(mixer_data), pointer :: p => NULL()
    end type mixer_data_ptr_type
    integer, intent(in)   :: this(2)
    type(mixer_data_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_AMIX
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%AMIX = f90wrap_AMIX
end subroutine f90wrap_mixer_data__set__AMIX

subroutine f90wrap_mixer_data_initialise(this)
    use mixer_module, only: mixer_data
    implicit none
    
    type mixer_data_ptr_type
        type(mixer_data), pointer :: p => NULL()
    end type mixer_data_ptr_type
    type(mixer_data_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_mixer_data_initialise

subroutine f90wrap_mixer_data_finalise(this)
    use mixer_module, only: mixer_data
    implicit none
    
    type mixer_data_ptr_type
        type(mixer_data), pointer :: p => NULL()
    end type mixer_data_ptr_type
    type(mixer_data_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_mixer_data_finalise

subroutine f90wrap_init_mixer_data
    use mixer_module, only: init_mixer_data
    implicit none
    
    call init_mixer_data()
end subroutine f90wrap_init_mixer_data

subroutine f90wrap_init_mixer_data_per
    use mixer_module, only: init_mixer_data_per
    implicit none
    
    call init_mixer_data_per()
end subroutine f90wrap_init_mixer_data_per

subroutine f90wrap_init_mixer_data_iso
    use mixer_module, only: init_mixer_data_iso
    implicit none
    
    call init_mixer_data_iso()
end subroutine f90wrap_init_mixer_data_iso

subroutine f90wrap_destroy_mixer
    use mixer_module, only: destroy_mixer
    implicit none
    
    call destroy_mixer()
end subroutine f90wrap_destroy_mixer

subroutine f90wrap_mixing(iter, xout, xin, res, n0, n1, n2, n3, n4, n5, n6, n7)
    use mixer_module, only: mixing
    implicit none
    
    integer, intent(in) :: iter
    real(8), intent(inout), dimension(n0,n1,n2,n3) :: xout
    real(8), intent(inout), dimension(n4,n5,n6,n7) :: xin
    real(8), intent(out) :: res
    integer :: n0
    !f2py intent(hide), depend(xout) :: n0 = shape(xout,0)
    integer :: n1
    !f2py intent(hide), depend(xout) :: n1 = shape(xout,1)
    integer :: n2
    !f2py intent(hide), depend(xout) :: n2 = shape(xout,2)
    integer :: n3
    !f2py intent(hide), depend(xout) :: n3 = shape(xout,3)
    integer :: n4
    !f2py intent(hide), depend(xin) :: n4 = shape(xin,0)
    integer :: n5
    !f2py intent(hide), depend(xin) :: n5 = shape(xin,1)
    integer :: n6
    !f2py intent(hide), depend(xin) :: n6 = shape(xin,2)
    integer :: n7
    !f2py intent(hide), depend(xin) :: n7 = shape(xin,3)
    call mixing(Iter=iter, XOUT=xout, XIN=xin, Res=res)
end subroutine f90wrap_mixing

subroutine f90wrap_mixing_iso(iter, xout, xin, res, n0, n1, n2, n3)
    use mixer_module, only: mixing_iso
    implicit none
    
    integer, intent(in) :: iter
    real(8), intent(inout), dimension(n0,n1) :: xout
    real(8), intent(inout), dimension(n2,n3) :: xin
    real(8), intent(out) :: res
    integer :: n0
    !f2py intent(hide), depend(xout) :: n0 = shape(xout,0)
    integer :: n1
    !f2py intent(hide), depend(xout) :: n1 = shape(xout,1)
    integer :: n2
    !f2py intent(hide), depend(xin) :: n2 = shape(xin,0)
    integer :: n3
    !f2py intent(hide), depend(xin) :: n3 = shape(xin,1)
    call mixing_iso(Iter=iter, XOUT=xout, XIN=xin, Res=res)
end subroutine f90wrap_mixing_iso

subroutine f90wrap_anderson_mixing(iiter, xout, xin, err, n0, n1)
    use mixer_module, only: anderson_mixing
    implicit none
    
    integer, intent(in) :: iiter
    real(8), intent(inout), dimension(n0) :: xout
    real(8), intent(inout), dimension(n1) :: xin
    real(8), intent(out) :: err
    integer :: n0
    !f2py intent(hide), depend(xout) :: n0 = shape(xout,0)
    integer :: n1
    !f2py intent(hide), depend(xin) :: n1 = shape(xin,0)
    call anderson_mixing(IITER=iiter, XOUT=xout, XIN=xin, err=err)
end subroutine f90wrap_anderson_mixing

subroutine f90wrap_om1c(nam, nuh, sp, dfp, voma, n0, n1, n2, n3)
    use mixer_module, only: om1c
    implicit none
    
    integer, intent(in) :: nam
    integer, intent(in) :: nuh
    real(8), intent(in) :: sp
    real(8), intent(in), dimension(n0,n1) :: dfp
    real(8), intent(inout), dimension(n2,n3) :: voma
    integer :: n0
    !f2py intent(hide), depend(dfp) :: n0 = shape(dfp,0)
    integer :: n1
    !f2py intent(hide), depend(dfp) :: n1 = shape(dfp,1)
    integer :: n2
    !f2py intent(hide), depend(voma) :: n2 = shape(voma,0)
    integer :: n3
    !f2py intent(hide), depend(voma) :: n3 = shape(voma,1)
    call om1c(NAM=nam, NUH=nuh, SP=sp, DFP=dfp, VOMA=voma)
end subroutine f90wrap_om1c

subroutine f90wrap_amst(beta, w0, nam, nuh, dxp, dfp, sp, xl, fl, voma, xn, n0, n1, n2, n3, n4, n5, n6, n7, n8)
    use mixer_module, only: amst
    implicit none
    
    real(8), intent(in) :: beta
    real(8), intent(in) :: w0
    integer, intent(in) :: nam
    integer, intent(in) :: nuh
    real(8), intent(in), dimension(n0,n1) :: dxp
    real(8), intent(in), dimension(n2,n3) :: dfp
    real(8), intent(in) :: sp
    real(8), intent(in), dimension(n4) :: xl
    real(8), intent(in), dimension(n5) :: fl
    real(8), intent(in), dimension(n6,n7) :: voma
    real(8), intent(inout), dimension(n8) :: xn
    integer :: n0
    !f2py intent(hide), depend(dxp) :: n0 = shape(dxp,0)
    integer :: n1
    !f2py intent(hide), depend(dxp) :: n1 = shape(dxp,1)
    integer :: n2
    !f2py intent(hide), depend(dfp) :: n2 = shape(dfp,0)
    integer :: n3
    !f2py intent(hide), depend(dfp) :: n3 = shape(dfp,1)
    integer :: n4
    !f2py intent(hide), depend(xl) :: n4 = shape(xl,0)
    integer :: n5
    !f2py intent(hide), depend(fl) :: n5 = shape(fl,0)
    integer :: n6
    !f2py intent(hide), depend(voma) :: n6 = shape(voma,0)
    integer :: n7
    !f2py intent(hide), depend(voma) :: n7 = shape(voma,1)
    integer :: n8
    !f2py intent(hide), depend(xn) :: n8 = shape(xn,0)
    call amst(BETA=beta, W0=w0, NAM=nam, NUH=nuh, DXP=dxp, DFP=dfp, SP=sp, XL=xl, FL=fl, VOMA=voma, XN=xn)
end subroutine f90wrap_amst

subroutine f90wrap_rpulay_mixing(iter, xout, xin, err, n0, n1)
    use mixer_module, only: rpulay_mixing
    implicit none
    
    integer(4), intent(in) :: iter
    real(8), intent(inout), dimension(n0) :: xout
    real(8), intent(inout), dimension(n1) :: xin
    real(8), intent(out) :: err
    integer :: n0
    !f2py intent(hide), depend(xout) :: n0 = shape(xout,0)
    integer :: n1
    !f2py intent(hide), depend(xin) :: n1 = shape(xin,0)
    call rpulay_mixing(Iter=iter, XOUT=xout, XIN=xin, err=err)
end subroutine f90wrap_rpulay_mixing

subroutine f90wrap_rpulay_mix(beta, w0, dime, nh, dxl, drl, xl, rl, xn, n0, n1, n2, n3, n4, n5, n6)
    use mixer_module, only: rpulay_mix
    implicit none
    
    real(8), intent(in) :: beta
    real(8), intent(in) :: w0
    integer(4), intent(in) :: dime
    integer(4), intent(in) :: nh
    real(8), intent(in), dimension(n0,n1) :: dxl
    real(8), intent(in), dimension(n2,n3) :: drl
    real(8), intent(in), dimension(n4) :: xl
    real(8), intent(in), dimension(n5) :: rl
    real(8), intent(inout), dimension(n6) :: xn
    integer :: n0
    !f2py intent(hide), depend(dxl) :: n0 = shape(dxl,0)
    integer :: n1
    !f2py intent(hide), depend(dxl) :: n1 = shape(dxl,1)
    integer :: n2
    !f2py intent(hide), depend(drl) :: n2 = shape(drl,0)
    integer :: n3
    !f2py intent(hide), depend(drl) :: n3 = shape(drl,1)
    integer :: n4
    !f2py intent(hide), depend(xl) :: n4 = shape(xl,0)
    integer :: n5
    !f2py intent(hide), depend(rl) :: n5 = shape(rl,0)
    integer :: n6
    !f2py intent(hide), depend(xn) :: n6 = shape(xn,0)
    call rpulay_mix(beta=beta, w0=w0, dime=dime, NH=nh, DXL=dxl, DRL=drl, XL=xl, RL=rl, XN=xn)
end subroutine f90wrap_rpulay_mix

subroutine f90wrap_init_kerker
    use mixer_module, only: init_kerker
    implicit none
    
    call init_kerker()
end subroutine f90wrap_init_kerker

subroutine f90wrap_rpulayk_mixing(iter, rlg, xing, n0, n1)
    use mixer_module, only: rpulayk_mixing
    implicit none
    
    integer(4), intent(in) :: iter
    complex(8), intent(inout), dimension(n0) :: rlg
    complex(8), intent(inout), dimension(n1) :: xing
    integer :: n0
    !f2py intent(hide), depend(rlg) :: n0 = shape(rlg,0)
    integer :: n1
    !f2py intent(hide), depend(xing) :: n1 = shape(xing,0)
    call rpulayk_mixing(Iter=iter, RLg=rlg, XINg=xing)
end subroutine f90wrap_rpulayk_mixing

subroutine f90wrap_rpulayk_mix(beta, w0, dime, nh, dxl, drl, xl, rl, xn, n0, n1, n2, n3, n4, n5, n6)
    use mixer_module, only: rpulayk_mix
    implicit none
    
    real(8), intent(in) :: beta
    real(8), intent(in) :: w0
    integer(4), intent(in) :: dime
    integer(4), intent(in) :: nh
    complex(8), intent(in), dimension(n0,n1) :: dxl
    complex(8), intent(in), dimension(n2,n3) :: drl
    complex(8), intent(in), dimension(n4) :: xl
    complex(8), intent(in), dimension(n5) :: rl
    complex(8), intent(inout), dimension(n6) :: xn
    integer :: n0
    !f2py intent(hide), depend(dxl) :: n0 = shape(dxl,0)
    integer :: n1
    !f2py intent(hide), depend(dxl) :: n1 = shape(dxl,1)
    integer :: n2
    !f2py intent(hide), depend(drl) :: n2 = shape(drl,0)
    integer :: n3
    !f2py intent(hide), depend(drl) :: n3 = shape(drl,1)
    integer :: n4
    !f2py intent(hide), depend(xl) :: n4 = shape(xl,0)
    integer :: n5
    !f2py intent(hide), depend(rl) :: n5 = shape(rl,0)
    integer :: n6
    !f2py intent(hide), depend(xn) :: n6 = shape(xn,0)
    call rpulayk_mix(beta=beta, w0=w0, dime=dime, NH=nh, DXL=dxl, DRL=drl, XL=xl, RL=rl, XN=xn)
end subroutine f90wrap_rpulayk_mix

subroutine f90wrap_resta_mixing(iter, rlg, xing, n0, n1)
    use mixer_module, only: resta_mixing
    implicit none
    
    integer(4), intent(in) :: iter
    complex(8), intent(inout), dimension(n0) :: rlg
    complex(8), intent(inout), dimension(n1) :: xing
    integer :: n0
    !f2py intent(hide), depend(rlg) :: n0 = shape(rlg,0)
    integer :: n1
    !f2py intent(hide), depend(xing) :: n1 = shape(xing,0)
    call resta_mixing(Iter=iter, RLg=rlg, XINg=xing)
end subroutine f90wrap_resta_mixing

subroutine f90wrap_init_resta
    use mixer_module, only: init_resta
    implicit none
    
    call init_resta()
end subroutine f90wrap_init_resta

! End of module mixer_module defined in file Mixer_module.fpp


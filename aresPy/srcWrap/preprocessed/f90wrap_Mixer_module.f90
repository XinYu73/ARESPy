! Module mixer_module defined in file Mixer_module.fpp

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

subroutine f90wrap_mixer_data__array__DRGL(this, nd, dtype, dshape, dloc)
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
    if (allocated(this_ptr%p%DRGL)) then
        dshape(1:2) = shape(this_ptr%p%DRGL)
        dloc = loc(this_ptr%p%DRGL)
    else
        dloc = 0
    end if
end subroutine f90wrap_mixer_data__array__DRGL

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

subroutine f90wrap_init_mixer(nps)
    use mixer_module, only: init_mixer
    implicit none
    
    integer(4), intent(in) :: nps
    call init_mixer(nps=nps)
end subroutine f90wrap_init_mixer

subroutine f90wrap_destroy_mixer
    use mixer_module, only: destroy_mixer
    implicit none
    
    call destroy_mixer()
end subroutine f90wrap_destroy_mixer

subroutine f90wrap_mixing(iter, xout, xin, res, n0, n1, n2, n3)
    use mixer_module, only: mixing
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
    call mixing(Iter=iter, XOUT=xout, XIN=xin, Res=res)
end subroutine f90wrap_mixing

subroutine f90wrap_anderson_mixing(iiter, xin, xout, err, n0, n1)
    use mixer_module, only: anderson_mixing
    implicit none
    
    integer, intent(in) :: iiter
    real(8), intent(inout), dimension(n0) :: xin
    real(8), intent(inout), dimension(n1) :: xout
    real(8), intent(out) :: err
    integer :: n0
    !f2py intent(hide), depend(xin) :: n0 = shape(xin,0)
    integer :: n1
    !f2py intent(hide), depend(xout) :: n1 = shape(xout,0)
    call anderson_mixing(IITER=iiter, XIN=xin, XOUT=xout, err=err)
end subroutine f90wrap_anderson_mixing

subroutine f90wrap_om1c(nam, nuh, sp, dfp, voma)
    use mixer_module, only: om1c
    implicit none
    
    integer, intent(in) :: nam
    integer, intent(in) :: nuh
    real :: sp
    real :: dfp
    real :: voma
    call om1c(NAM=nam, NUH=nuh, sp=sp, dfp=dfp, voma=voma)
end subroutine f90wrap_om1c

subroutine f90wrap_amst(beta, w0, nam, nuh, dxp, dfp, sp, xl, fl, voma, xn)
    use mixer_module, only: amst
    implicit none
    
    real :: beta
    real :: w0
    integer, intent(in) :: nam
    integer, intent(in) :: nuh
    real :: dxp
    real :: dfp
    real :: sp
    real :: xl
    real :: fl
    real :: voma
    real :: xn
    call amst(beta=beta, w0=w0, NAM=nam, NUH=nuh, dxp=dxp, dfp=dfp, sp=sp, xl=xl, fl=fl, voma=voma, xn=xn)
end subroutine f90wrap_amst

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

subroutine f90wrap_mixer_module__get__NAM(f90wrap_NAM)
    use mixer_module, only: mixer_module_NAM => NAM
    implicit none
    integer(4), intent(out) :: f90wrap_NAM
    
    f90wrap_NAM = mixer_module_NAM
end subroutine f90wrap_mixer_module__get__NAM

subroutine f90wrap_mixer_module__set__NAM(f90wrap_NAM)
    use mixer_module, only: mixer_module_NAM => NAM
    implicit none
    integer(4), intent(in) :: f90wrap_NAM
    
    mixer_module_NAM = f90wrap_NAM
end subroutine f90wrap_mixer_module__set__NAM

subroutine f90wrap_mixer_module__get__NUH(f90wrap_NUH)
    use mixer_module, only: mixer_module_NUH => NUH
    implicit none
    integer(4), intent(out) :: f90wrap_NUH
    
    f90wrap_NUH = mixer_module_NUH
end subroutine f90wrap_mixer_module__get__NUH

subroutine f90wrap_mixer_module__set__NUH(f90wrap_NUH)
    use mixer_module, only: mixer_module_NUH => NUH
    implicit none
    integer(4), intent(in) :: f90wrap_NUH
    
    mixer_module_NUH = f90wrap_NUH
end subroutine f90wrap_mixer_module__set__NUH

! End of module mixer_module defined in file Mixer_module.fpp


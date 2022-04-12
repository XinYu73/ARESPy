! Module relax_module defined in file Relax_module.fpp

subroutine f90wrap_relax_type__get__TIMSmax(this, f90wrap_TIMSmax)
    use relax_module, only: relax_type
    implicit none
    type relax_type_ptr_type
        type(relax_type), pointer :: p => NULL()
    end type relax_type_ptr_type
    integer, intent(in)   :: this(2)
    type(relax_type_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_TIMSmax
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_TIMSmax = this_ptr%p%TIMSmax
end subroutine f90wrap_relax_type__get__TIMSmax

subroutine f90wrap_relax_type__set__TIMSmax(this, f90wrap_TIMSmax)
    use relax_module, only: relax_type
    implicit none
    type relax_type_ptr_type
        type(relax_type), pointer :: p => NULL()
    end type relax_type_ptr_type
    integer, intent(in)   :: this(2)
    type(relax_type_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_TIMSmax
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%TIMSmax = f90wrap_TIMSmax
end subroutine f90wrap_relax_type__set__TIMSmax

subroutine f90wrap_relax_type__get__TIMS(this, f90wrap_TIMS)
    use relax_module, only: relax_type
    implicit none
    type relax_type_ptr_type
        type(relax_type), pointer :: p => NULL()
    end type relax_type_ptr_type
    integer, intent(in)   :: this(2)
    type(relax_type_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_TIMS
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_TIMS = this_ptr%p%TIMS
end subroutine f90wrap_relax_type__get__TIMS

subroutine f90wrap_relax_type__set__TIMS(this, f90wrap_TIMS)
    use relax_module, only: relax_type
    implicit none
    type relax_type_ptr_type
        type(relax_type), pointer :: p => NULL()
    end type relax_type_ptr_type
    integer, intent(in)   :: this(2)
    type(relax_type_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_TIMS
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%TIMS = f90wrap_TIMS
end subroutine f90wrap_relax_type__set__TIMS

subroutine f90wrap_relax_type__get__alpha(this, f90wrap_alpha)
    use relax_module, only: relax_type
    implicit none
    type relax_type_ptr_type
        type(relax_type), pointer :: p => NULL()
    end type relax_type_ptr_type
    integer, intent(in)   :: this(2)
    type(relax_type_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_alpha
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_alpha = this_ptr%p%alpha
end subroutine f90wrap_relax_type__get__alpha

subroutine f90wrap_relax_type__set__alpha(this, f90wrap_alpha)
    use relax_module, only: relax_type
    implicit none
    type relax_type_ptr_type
        type(relax_type), pointer :: p => NULL()
    end type relax_type_ptr_type
    integer, intent(in)   :: this(2)
    type(relax_type_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_alpha
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%alpha = f90wrap_alpha
end subroutine f90wrap_relax_type__set__alpha

subroutine f90wrap_relax_type__array__FACTion(this, nd, dtype, dshape, dloc)
    use relax_module, only: relax_type
    implicit none
    type relax_type_ptr_type
        type(relax_type), pointer :: p => NULL()
    end type relax_type_ptr_type
    integer, intent(in) :: this(2)
    type(relax_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%FACTion)) then
        dshape(1:1) = shape(this_ptr%p%FACTion)
        dloc = loc(this_ptr%p%FACTion)
    else
        dloc = 0
    end if
end subroutine f90wrap_relax_type__array__FACTion

subroutine f90wrap_relax_type__array__Fiond(this, nd, dtype, dshape, dloc)
    use relax_module, only: relax_type
    implicit none
    type relax_type_ptr_type
        type(relax_type), pointer :: p => NULL()
    end type relax_type_ptr_type
    integer, intent(in) :: this(2)
    type(relax_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%Fiond)) then
        dshape(1:2) = shape(this_ptr%p%Fiond)
        dloc = loc(this_ptr%p%Fiond)
    else
        dloc = 0
    end if
end subroutine f90wrap_relax_type__array__Fiond

subroutine f90wrap_relax_type__array__Fceld(this, nd, dtype, dshape, dloc)
    use relax_module, only: relax_type
    implicit none
    type relax_type_ptr_type
        type(relax_type), pointer :: p => NULL()
    end type relax_type_ptr_type
    integer, intent(in) :: this(2)
    type(relax_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%Fceld)) then
        dshape(1:2) = shape(this_ptr%p%Fceld)
        dloc = loc(this_ptr%p%Fceld)
    else
        dloc = 0
    end if
end subroutine f90wrap_relax_type__array__Fceld

subroutine f90wrap_relax_type__array__VELion(this, nd, dtype, dshape, dloc)
    use relax_module, only: relax_type
    implicit none
    type relax_type_ptr_type
        type(relax_type), pointer :: p => NULL()
    end type relax_type_ptr_type
    integer, intent(in) :: this(2)
    type(relax_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%VELion)) then
        dshape(1:2) = shape(this_ptr%p%VELion)
        dloc = loc(this_ptr%p%VELion)
    else
        dloc = 0
    end if
end subroutine f90wrap_relax_type__array__VELion

subroutine f90wrap_relax_type__array__VELcel(this, nd, dtype, dshape, dloc)
    use relax_module, only: relax_type
    implicit none
    type relax_type_ptr_type
        type(relax_type), pointer :: p => NULL()
    end type relax_type_ptr_type
    integer, intent(in) :: this(2)
    type(relax_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%VELcel)) then
        dshape(1:2) = shape(this_ptr%p%VELcel)
        dloc = loc(this_ptr%p%VELcel)
    else
        dloc = 0
    end if
end subroutine f90wrap_relax_type__array__VELcel

subroutine f90wrap_relax_type__get__FACTcel(this, f90wrap_FACTcel)
    use relax_module, only: relax_type
    implicit none
    type relax_type_ptr_type
        type(relax_type), pointer :: p => NULL()
    end type relax_type_ptr_type
    integer, intent(in)   :: this(2)
    type(relax_type_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_FACTcel
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_FACTcel = this_ptr%p%FACTcel
end subroutine f90wrap_relax_type__get__FACTcel

subroutine f90wrap_relax_type__set__FACTcel(this, f90wrap_FACTcel)
    use relax_module, only: relax_type
    implicit none
    type relax_type_ptr_type
        type(relax_type), pointer :: p => NULL()
    end type relax_type_ptr_type
    integer, intent(in)   :: this(2)
    type(relax_type_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_FACTcel
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%FACTcel = f90wrap_FACTcel
end subroutine f90wrap_relax_type__set__FACTcel

subroutine f90wrap_relax_type__get__NEG(this, f90wrap_NEG)
    use relax_module, only: relax_type
    implicit none
    type relax_type_ptr_type
        type(relax_type), pointer :: p => NULL()
    end type relax_type_ptr_type
    integer, intent(in)   :: this(2)
    type(relax_type_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_NEG
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_NEG = this_ptr%p%NEG
end subroutine f90wrap_relax_type__get__NEG

subroutine f90wrap_relax_type__set__NEG(this, f90wrap_NEG)
    use relax_module, only: relax_type
    implicit none
    type relax_type_ptr_type
        type(relax_type), pointer :: p => NULL()
    end type relax_type_ptr_type
    integer, intent(in)   :: this(2)
    type(relax_type_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_NEG
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%NEG = f90wrap_NEG
end subroutine f90wrap_relax_type__set__NEG

subroutine f90wrap_relax_type__get__lNEG(this, f90wrap_lNEG)
    use relax_module, only: relax_type
    implicit none
    type relax_type_ptr_type
        type(relax_type), pointer :: p => NULL()
    end type relax_type_ptr_type
    integer, intent(in)   :: this(2)
    type(relax_type_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_lNEG
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_lNEG = this_ptr%p%lNEG
end subroutine f90wrap_relax_type__get__lNEG

subroutine f90wrap_relax_type__set__lNEG(this, f90wrap_lNEG)
    use relax_module, only: relax_type
    implicit none
    type relax_type_ptr_type
        type(relax_type), pointer :: p => NULL()
    end type relax_type_ptr_type
    integer, intent(in)   :: this(2)
    type(relax_type_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_lNEG
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%lNEG = f90wrap_lNEG
end subroutine f90wrap_relax_type__set__lNEG

subroutine f90wrap_relax_type_initialise(this)
    use relax_module, only: relax_type
    implicit none
    
    type relax_type_ptr_type
        type(relax_type), pointer :: p => NULL()
    end type relax_type_ptr_type
    type(relax_type_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_relax_type_initialise

subroutine f90wrap_relax_type_finalise(this)
    use relax_module, only: relax_type
    implicit none
    
    type relax_type_ptr_type
        type(relax_type), pointer :: p => NULL()
    end type relax_type_ptr_type
    type(relax_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_relax_type_finalise

subroutine f90wrap_initialize_relax
    use relax_module, only: initialize_relax
    implicit none
    
    call initialize_relax()
end subroutine f90wrap_initialize_relax

subroutine f90wrap_destroy_relax
    use relax_module, only: destroy_relax
    implicit none
    
    call destroy_relax()
end subroutine f90wrap_destroy_relax

subroutine f90wrap_relaxer(nstep)
    use relax_module, only: relaxer
    implicit none
    
    integer(4), intent(in) :: nstep
    call relaxer(nstep=nstep)
end subroutine f90wrap_relaxer

subroutine f90wrap_fire_relax(lat, pos, fion, fcel, n0, n1, n2, n3, n4, n5, n6, n7)
    use relax_module, only: fire_relax
    implicit none
    
    real(8), intent(inout), dimension(n0,n1) :: lat
    real(8), intent(inout), dimension(n2,n3) :: pos
    real(8), intent(in), dimension(n4,n5) :: fion
    real(8), intent(in), dimension(n6,n7) :: fcel
    integer :: n0
    !f2py intent(hide), depend(lat) :: n0 = shape(lat,0)
    integer :: n1
    !f2py intent(hide), depend(lat) :: n1 = shape(lat,1)
    integer :: n2
    !f2py intent(hide), depend(pos) :: n2 = shape(pos,0)
    integer :: n3
    !f2py intent(hide), depend(pos) :: n3 = shape(pos,1)
    integer :: n4
    !f2py intent(hide), depend(fion) :: n4 = shape(fion,0)
    integer :: n5
    !f2py intent(hide), depend(fion) :: n5 = shape(fion,1)
    integer :: n6
    !f2py intent(hide), depend(fcel) :: n6 = shape(fcel,0)
    integer :: n7
    !f2py intent(hide), depend(fcel) :: n7 = shape(fcel,1)
    call fire_relax(LAT=lat, POS=pos, Fion=fion, Fcel=fcel)
end subroutine f90wrap_fire_relax

subroutine f90wrap_force_on_cell(stress, lat, fcel, n0, n1, n2, n3, n4, n5)
    use relax_module, only: force_on_cell
    implicit none
    
    real(8), intent(in), dimension(n0,n1) :: stress
    real(8), intent(in), dimension(n2,n3) :: lat
    real(8), intent(inout), dimension(n4,n5) :: fcel
    integer :: n0
    !f2py intent(hide), depend(stress) :: n0 = shape(stress,0)
    integer :: n1
    !f2py intent(hide), depend(stress) :: n1 = shape(stress,1)
    integer :: n2
    !f2py intent(hide), depend(lat) :: n2 = shape(lat,0)
    integer :: n3
    !f2py intent(hide), depend(lat) :: n3 = shape(lat,1)
    integer :: n4
    !f2py intent(hide), depend(fcel) :: n4 = shape(fcel,0)
    integer :: n5
    !f2py intent(hide), depend(fcel) :: n5 = shape(fcel,1)
    call force_on_cell(stress=stress, LAT=lat, Fcel=fcel)
end subroutine f90wrap_force_on_cell

subroutine f90wrap_fix_direction(stress, n0, n1)
    use relax_module, only: fix_direction
    implicit none
    
    real(8), dimension(n0,n1) :: stress
    integer :: n0
    !f2py intent(hide), depend(stress) :: n0 = shape(stress,0)
    integer :: n1
    !f2py intent(hide), depend(stress) :: n1 = shape(stress,1)
    call fix_direction(stress=stress)
end subroutine f90wrap_fix_direction

subroutine f90wrap_relax_module__get__pstress(f90wrap_pstress)
    use relax_module, only: relax_module_pstress => pstress
    implicit none
    real(8), intent(out) :: f90wrap_pstress
    
    f90wrap_pstress = relax_module_pstress
end subroutine f90wrap_relax_module__get__pstress

subroutine f90wrap_relax_module__set__pstress(f90wrap_pstress)
    use relax_module, only: relax_module_pstress => pstress
    implicit none
    real(8), intent(in) :: f90wrap_pstress
    
    relax_module_pstress = f90wrap_pstress
end subroutine f90wrap_relax_module__set__pstress

subroutine f90wrap_relax_module__get__Ldone(f90wrap_Ldone)
    use relax_module, only: relax_module_Ldone => Ldone
    implicit none
    logical, intent(out) :: f90wrap_Ldone
    
    f90wrap_Ldone = relax_module_Ldone
end subroutine f90wrap_relax_module__get__Ldone

subroutine f90wrap_relax_module__set__Ldone(f90wrap_Ldone)
    use relax_module, only: relax_module_Ldone => Ldone
    implicit none
    logical, intent(in) :: f90wrap_Ldone
    
    relax_module_Ldone = f90wrap_Ldone
end subroutine f90wrap_relax_module__set__Ldone

subroutine f90wrap_relax_module__get__Lfirst(f90wrap_Lfirst)
    use relax_module, only: relax_module_Lfirst => Lfirst
    implicit none
    logical, intent(out) :: f90wrap_Lfirst
    
    f90wrap_Lfirst = relax_module_Lfirst
end subroutine f90wrap_relax_module__get__Lfirst

subroutine f90wrap_relax_module__set__Lfirst(f90wrap_Lfirst)
    use relax_module, only: relax_module_Lfirst => Lfirst
    implicit none
    logical, intent(in) :: f90wrap_Lfirst
    
    relax_module_Lfirst = f90wrap_Lfirst
end subroutine f90wrap_relax_module__set__Lfirst

! End of module relax_module defined in file Relax_module.fpp


! Module struct_module defined in file Struct_module.fpp

subroutine f90wrap_struct_type__array__Zion(this, nd, dtype, dshape, dloc)
    use struct_module, only: struct_type
    implicit none
    type struct_type_ptr_type
        type(struct_type), pointer :: p => NULL()
    end type struct_type_ptr_type
    integer, intent(in) :: this(2)
    type(struct_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%Zion)) then
        dshape(1:1) = shape(this_ptr%p%Zion)
        dloc = loc(this_ptr%p%Zion)
    else
        dloc = 0
    end if
end subroutine f90wrap_struct_type__array__Zion

subroutine f90wrap_struct_type__array__nati(this, nd, dtype, dshape, dloc)
    use struct_module, only: struct_type
    implicit none
    type struct_type_ptr_type
        type(struct_type), pointer :: p => NULL()
    end type struct_type_ptr_type
    integer, intent(in) :: this(2)
    type(struct_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%nati)) then
        dshape(1:1) = shape(this_ptr%p%nati)
        dloc = loc(this_ptr%p%nati)
    else
        dloc = 0
    end if
end subroutine f90wrap_struct_type__array__nati

subroutine f90wrap_struct_type__array__eleid(this, nd, dtype, dshape, dloc)
    use struct_module, only: struct_type
    implicit none
    type struct_type_ptr_type
        type(struct_type), pointer :: p => NULL()
    end type struct_type_ptr_type
    integer, intent(in) :: this(2)
    type(struct_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%eleid)) then
        dshape(1:1) = shape(this_ptr%p%eleid)
        dloc = loc(this_ptr%p%eleid)
    else
        dloc = 0
    end if
end subroutine f90wrap_struct_type__array__eleid

subroutine f90wrap_struct_type__array__pos(this, nd, dtype, dshape, dloc)
    use struct_module, only: struct_type
    implicit none
    type struct_type_ptr_type
        type(struct_type), pointer :: p => NULL()
    end type struct_type_ptr_type
    integer, intent(in) :: this(2)
    type(struct_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%pos)) then
        dshape(1:2) = shape(this_ptr%p%pos)
        dloc = loc(this_ptr%p%pos)
    else
        dloc = 0
    end if
end subroutine f90wrap_struct_type__array__pos

subroutine f90wrap_struct_type__array__poscar(this, nd, dtype, dshape, dloc)
    use struct_module, only: struct_type
    implicit none
    type struct_type_ptr_type
        type(struct_type), pointer :: p => NULL()
    end type struct_type_ptr_type
    integer, intent(in) :: this(2)
    type(struct_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%poscar)) then
        dshape(1:2) = shape(this_ptr%p%poscar)
        dloc = loc(this_ptr%p%poscar)
    else
        dloc = 0
    end if
end subroutine f90wrap_struct_type__array__poscar

subroutine f90wrap_struct_type__array__stress(this, nd, dtype, dshape, dloc)
    use struct_module, only: struct_type
    implicit none
    type struct_type_ptr_type
        type(struct_type), pointer :: p => NULL()
    end type struct_type_ptr_type
    integer, intent(in) :: this(2)
    type(struct_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    dshape(1:2) = shape(this_ptr%p%stress)
    dloc = loc(this_ptr%p%stress)
end subroutine f90wrap_struct_type__array__stress

subroutine f90wrap_struct_type__array__forces(this, nd, dtype, dshape, dloc)
    use struct_module, only: struct_type
    implicit none
    type struct_type_ptr_type
        type(struct_type), pointer :: p => NULL()
    end type struct_type_ptr_type
    integer, intent(in) :: this(2)
    type(struct_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%forces)) then
        dshape(1:2) = shape(this_ptr%p%forces)
        dloc = loc(this_ptr%p%forces)
    else
        dloc = 0
    end if
end subroutine f90wrap_struct_type__array__forces

subroutine f90wrap_struct_type__array__mass(this, nd, dtype, dshape, dloc)
    use struct_module, only: struct_type
    implicit none
    type struct_type_ptr_type
        type(struct_type), pointer :: p => NULL()
    end type struct_type_ptr_type
    integer, intent(in) :: this(2)
    type(struct_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%mass)) then
        dshape(1:1) = shape(this_ptr%p%mass)
        dloc = loc(this_ptr%p%mass)
    else
        dloc = 0
    end if
end subroutine f90wrap_struct_type__array__mass

subroutine f90wrap_struct_type__array__zeta(this, nd, dtype, dshape, dloc)
    use struct_module, only: struct_type
    implicit none
    type struct_type_ptr_type
        type(struct_type), pointer :: p => NULL()
    end type struct_type_ptr_type
    integer, intent(in) :: this(2)
    type(struct_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%zeta)) then
        dshape(1:2) = shape(this_ptr%p%zeta)
        dloc = loc(this_ptr%p%zeta)
    else
        dloc = 0
    end if
end subroutine f90wrap_struct_type__array__zeta

subroutine f90wrap_struct_type__array__prinq(this, nd, dtype, dshape, dloc)
    use struct_module, only: struct_type
    implicit none
    type struct_type_ptr_type
        type(struct_type), pointer :: p => NULL()
    end type struct_type_ptr_type
    integer, intent(in) :: this(2)
    type(struct_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%prinq)) then
        dshape(1:2) = shape(this_ptr%p%prinq)
        dloc = loc(this_ptr%p%prinq)
    else
        dloc = 0
    end if
end subroutine f90wrap_struct_type__array__prinq

subroutine f90wrap_struct_type__array__Lmax(this, nd, dtype, dshape, dloc)
    use struct_module, only: struct_type
    implicit none
    type struct_type_ptr_type
        type(struct_type), pointer :: p => NULL()
    end type struct_type_ptr_type
    integer, intent(in) :: this(2)
    type(struct_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%Lmax)) then
        dshape(1:1) = shape(this_ptr%p%Lmax)
        dloc = loc(this_ptr%p%Lmax)
    else
        dloc = 0
    end if
end subroutine f90wrap_struct_type__array__Lmax

subroutine f90wrap_struct_type__array__elements(this, nd, dtype, dshape, dloc)
    use struct_module, only: struct_type
    implicit none
    type struct_type_ptr_type
        type(struct_type), pointer :: p => NULL()
    end type struct_type_ptr_type
    integer, intent(in) :: this(2)
    type(struct_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 2
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%elements)) then
        dshape(1:2) = (/len(this_ptr%p%elements(1)), shape(this_ptr%p%elements)/)
        dloc = loc(this_ptr%p%elements)
    else
        dloc = 0
    end if
end subroutine f90wrap_struct_type__array__elements

subroutine f90wrap_struct_type__array__coeff(this, nd, dtype, dshape, dloc)
    use struct_module, only: struct_type
    implicit none
    type struct_type_ptr_type
        type(struct_type), pointer :: p => NULL()
    end type struct_type_ptr_type
    integer, intent(in) :: this(2)
    type(struct_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 4
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%coeff)) then
        dshape(1:4) = shape(this_ptr%p%coeff)
        dloc = loc(this_ptr%p%coeff)
    else
        dloc = 0
    end if
end subroutine f90wrap_struct_type__array__coeff

subroutine f90wrap_struct_type__array__occupy(this, nd, dtype, dshape, dloc)
    use struct_module, only: struct_type
    implicit none
    type struct_type_ptr_type
        type(struct_type), pointer :: p => NULL()
    end type struct_type_ptr_type
    integer, intent(in) :: this(2)
    type(struct_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%occupy)) then
        dshape(1:1) = shape(this_ptr%p%occupy)
        dloc = loc(this_ptr%p%occupy)
    else
        dloc = 0
    end if
end subroutine f90wrap_struct_type__array__occupy

subroutine f90wrap_struct_type__get__Noccupy(this, f90wrap_Noccupy)
    use struct_module, only: struct_type
    implicit none
    type struct_type_ptr_type
        type(struct_type), pointer :: p => NULL()
    end type struct_type_ptr_type
    integer, intent(in)   :: this(2)
    type(struct_type_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_Noccupy
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_Noccupy = this_ptr%p%Noccupy
end subroutine f90wrap_struct_type__get__Noccupy

subroutine f90wrap_struct_type__set__Noccupy(this, f90wrap_Noccupy)
    use struct_module, only: struct_type
    implicit none
    type struct_type_ptr_type
        type(struct_type), pointer :: p => NULL()
    end type struct_type_ptr_type
    integer, intent(in)   :: this(2)
    type(struct_type_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_Noccupy
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%Noccupy = f90wrap_Noccupy
end subroutine f90wrap_struct_type__set__Noccupy

subroutine f90wrap_struct_type_initialise(this)
    use struct_module, only: struct_type
    implicit none
    
    type struct_type_ptr_type
        type(struct_type), pointer :: p => NULL()
    end type struct_type_ptr_type
    type(struct_type_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_struct_type_initialise

subroutine f90wrap_struct_type_finalise(this)
    use struct_module, only: struct_type
    implicit none
    
    type struct_type_ptr_type
        type(struct_type), pointer :: p => NULL()
    end type struct_type_ptr_type
    type(struct_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_struct_type_finalise

subroutine f90wrap_creat_struct(numtyp, numatom)
    use struct_module, only: creat_struct
    implicit none
    
    integer(4), intent(in) :: numtyp
    integer(4), intent(in) :: numatom
    call creat_struct(numtyp=numtyp, numatom=numatom)
end subroutine f90wrap_creat_struct

subroutine f90wrap_destroy_struct
    use struct_module, only: destroy_struct
    implicit none
    
    call destroy_struct()
end subroutine f90wrap_destroy_struct

subroutine f90wrap_struct_module__get__natom(f90wrap_natom)
    use struct_module, only: struct_module_natom => natom
    implicit none
    integer(4), intent(out) :: f90wrap_natom
    
    f90wrap_natom = struct_module_natom
end subroutine f90wrap_struct_module__get__natom

subroutine f90wrap_struct_module__set__natom(f90wrap_natom)
    use struct_module, only: struct_module_natom => natom
    implicit none
    integer(4), intent(in) :: f90wrap_natom
    
    struct_module_natom = f90wrap_natom
end subroutine f90wrap_struct_module__set__natom

subroutine f90wrap_struct_module__get__naty(f90wrap_naty)
    use struct_module, only: struct_module_naty => naty
    implicit none
    integer(4), intent(out) :: f90wrap_naty
    
    f90wrap_naty = struct_module_naty
end subroutine f90wrap_struct_module__get__naty

subroutine f90wrap_struct_module__set__naty(f90wrap_naty)
    use struct_module, only: struct_module_naty => naty
    implicit none
    integer(4), intent(in) :: f90wrap_naty
    
    struct_module_naty = f90wrap_naty
end subroutine f90wrap_struct_module__set__naty

subroutine f90wrap_struct_module__get__ncharge(f90wrap_ncharge)
    use struct_module, only: struct_module_ncharge => ncharge
    implicit none
    integer(4), intent(out) :: f90wrap_ncharge
    
    f90wrap_ncharge = struct_module_ncharge
end subroutine f90wrap_struct_module__get__ncharge

subroutine f90wrap_struct_module__set__ncharge(f90wrap_ncharge)
    use struct_module, only: struct_module_ncharge => ncharge
    implicit none
    integer(4), intent(in) :: f90wrap_ncharge
    
    struct_module_ncharge = f90wrap_ncharge
end subroutine f90wrap_struct_module__set__ncharge

subroutine f90wrap_struct_module__get__charge_ave(f90wrap_charge_ave)
    use struct_module, only: struct_module_charge_ave => charge_ave
    implicit none
    real(8), intent(out) :: f90wrap_charge_ave
    
    f90wrap_charge_ave = struct_module_charge_ave
end subroutine f90wrap_struct_module__get__charge_ave

subroutine f90wrap_struct_module__set__charge_ave(f90wrap_charge_ave)
    use struct_module, only: struct_module_charge_ave => charge_ave
    implicit none
    real(8), intent(in) :: f90wrap_charge_ave
    
    struct_module_charge_ave = f90wrap_charge_ave
end subroutine f90wrap_struct_module__set__charge_ave

subroutine f90wrap_struct_module__get__volume(f90wrap_volume)
    use struct_module, only: struct_module_volume => volume
    implicit none
    real(8), intent(out) :: f90wrap_volume
    
    f90wrap_volume = struct_module_volume
end subroutine f90wrap_struct_module__get__volume

subroutine f90wrap_struct_module__set__volume(f90wrap_volume)
    use struct_module, only: struct_module_volume => volume
    implicit none
    real(8), intent(in) :: f90wrap_volume
    
    struct_module_volume = f90wrap_volume
end subroutine f90wrap_struct_module__set__volume

subroutine f90wrap_struct_module__array__lat_mat(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use struct_module, only: struct_module_lat_mat => lat_mat
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    dshape(1:2) = shape(struct_module_lat_mat)
    dloc = loc(struct_module_lat_mat)
end subroutine f90wrap_struct_module__array__lat_mat

subroutine f90wrap_struct_module__array__lat_para(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use struct_module, only: struct_module_lat_para => lat_para
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    dshape(1:1) = shape(struct_module_lat_para)
    dloc = loc(struct_module_lat_para)
end subroutine f90wrap_struct_module__array__lat_para

subroutine f90wrap_struct_module__array__recip_lat(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use struct_module, only: struct_module_recip_lat => recip_lat
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    dshape(1:2) = shape(struct_module_recip_lat)
    dloc = loc(struct_module_recip_lat)
end subroutine f90wrap_struct_module__array__recip_lat

subroutine f90wrap_struct_module__array__reclat_para(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use struct_module, only: struct_module_reclat_para => reclat_para
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    dshape(1:1) = shape(struct_module_reclat_para)
    dloc = loc(struct_module_reclat_para)
end subroutine f90wrap_struct_module__array__reclat_para

subroutine f90wrap_struct_module__array__energy(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use struct_module, only: struct_module_energy => energy
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    dshape(1:1) = shape(struct_module_energy)
    dloc = loc(struct_module_energy)
end subroutine f90wrap_struct_module__array__energy

! End of module struct_module defined in file Struct_module.fpp


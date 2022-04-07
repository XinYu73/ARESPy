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

subroutine f90wrap_struct_type_initialise(this)
    use struct_module, only: struct_type
    implicit none

    type struct_type_ptr_type
        type(struct_type), pointer :: p => NULL()
    end type struct_type_ptr_type
    type(struct_type_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate (this_ptr%p)
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
    deallocate (this_ptr%p)
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

subroutine f90wrap_struct_module__get__nzion(f90wrap_nzion)
    use struct_module, only: struct_module_nzion => nzion
    implicit none
    integer(4), intent(out) :: f90wrap_nzion

    f90wrap_nzion = struct_module_nzion
end subroutine f90wrap_struct_module__get__nzion

subroutine f90wrap_struct_module__set__nzion(f90wrap_nzion)
    use struct_module, only: struct_module_nzion => nzion
    implicit none
    integer(4), intent(in) :: f90wrap_nzion

    struct_module_nzion = f90wrap_nzion
end subroutine f90wrap_struct_module__set__nzion

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
    real(8), intent(out) :: f90wrap_ncharge

    f90wrap_ncharge = struct_module_ncharge
end subroutine f90wrap_struct_module__get__ncharge

subroutine f90wrap_struct_module__set__ncharge(f90wrap_ncharge)
    use struct_module, only: struct_module_ncharge => ncharge
    implicit none
    real(8), intent(in) :: f90wrap_ncharge

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

subroutine f90wrap_struct_module__get__volsp(f90wrap_volsp)
    use struct_module, only: struct_module_volsp => volsp
    implicit none
    real(8), intent(out) :: f90wrap_volsp

    f90wrap_volsp = struct_module_volsp
end subroutine f90wrap_struct_module__get__volsp

subroutine f90wrap_struct_module__set__volsp(f90wrap_volsp)
    use struct_module, only: struct_module_volsp => volsp
    implicit none
    real(8), intent(in) :: f90wrap_volsp

    struct_module_volsp = f90wrap_volsp
end subroutine f90wrap_struct_module__set__volsp

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

subroutine f90wrap_struct_module__get__Eionion(f90wrap_Eionion)
    use struct_module, only: struct_module_Eionion => Eionion
    implicit none
    real(8), intent(out) :: f90wrap_Eionion

    f90wrap_Eionion = struct_module_Eionion
end subroutine f90wrap_struct_module__get__Eionion

subroutine f90wrap_struct_module__set__Eionion(f90wrap_Eionion)
    use struct_module, only: struct_module_Eionion => Eionion
    implicit none
    real(8), intent(in) :: f90wrap_Eionion

    struct_module_Eionion = f90wrap_Eionion
end subroutine f90wrap_struct_module__set__Eionion

subroutine f90wrap_struct_module__get__Eshift_ps(f90wrap_Eshift_ps)
    use struct_module, only: struct_module_Eshift_ps => Eshift_ps
    implicit none
    real(8), intent(out) :: f90wrap_Eshift_ps

    f90wrap_Eshift_ps = struct_module_Eshift_ps
end subroutine f90wrap_struct_module__get__Eshift_ps

subroutine f90wrap_struct_module__set__Eshift_ps(f90wrap_Eshift_ps)
    use struct_module, only: struct_module_Eshift_ps => Eshift_ps
    implicit none
    real(8), intent(in) :: f90wrap_Eshift_ps

    struct_module_Eshift_ps = f90wrap_Eshift_ps
end subroutine f90wrap_struct_module__set__Eshift_ps

subroutine f90wrap_struct_module__get__Eshift_tot(f90wrap_Eshift_tot)
    use struct_module, only: struct_module_Eshift_tot => Eshift_tot
    implicit none
    real(8), intent(out) :: f90wrap_Eshift_tot

    f90wrap_Eshift_tot = struct_module_Eshift_tot
end subroutine f90wrap_struct_module__get__Eshift_tot

subroutine f90wrap_struct_module__set__Eshift_tot(f90wrap_Eshift_tot)
    use struct_module, only: struct_module_Eshift_tot => Eshift_tot
    implicit none
    real(8), intent(in) :: f90wrap_Eshift_tot

    struct_module_Eshift_tot = f90wrap_Eshift_tot
end subroutine f90wrap_struct_module__set__Eshift_tot

subroutine f90wrap_struct_module__array__Opsym(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use struct_module, only: struct_module_opsym => opsym
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 3
    dtype = 12
    dshape(1:3) = shape(struct_module_Opsym)
    dloc = loc(struct_module_Opsym)
end subroutine f90wrap_struct_module__array__Opsym

subroutine f90wrap_struct_module__array__Otrans(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use struct_module, only: struct_module_otrans => otrans
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 2
    dtype = 12
    dshape(1:2) = shape(struct_module_Otrans)
    dloc = loc(struct_module_Otrans)
end subroutine f90wrap_struct_module__array__Otrans

subroutine f90wrap_struct_module__get__nsym(f90wrap_nsym)
    use struct_module, only: struct_module_nsym => nsym
    implicit none
    integer(4), intent(out) :: f90wrap_nsym

    f90wrap_nsym = struct_module_nsym
end subroutine f90wrap_struct_module__get__nsym

subroutine f90wrap_struct_module__set__nsym(f90wrap_nsym)
    use struct_module, only: struct_module_nsym => nsym
    implicit none
    integer(4), intent(in) :: f90wrap_nsym

    struct_module_nsym = f90wrap_nsym
end subroutine f90wrap_struct_module__set__nsym

subroutine f90wrap_struct_module__get__num_t(f90wrap_num_t)
    use struct_module, only: struct_module_num_t => num_t
    implicit none
    integer(4), intent(out) :: f90wrap_num_t

    f90wrap_num_t = struct_module_num_t
end subroutine f90wrap_struct_module__get__num_t

subroutine f90wrap_struct_module__set__num_t(f90wrap_num_t)
    use struct_module, only: struct_module_num_t => num_t
    implicit none
    integer(4), intent(in) :: f90wrap_num_t

    struct_module_num_t = f90wrap_num_t
end subroutine f90wrap_struct_module__set__num_t

subroutine f90wrap_struct_module__array__c_i(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use struct_module, only: struct_module_c_i => c_i
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 5
    dshape(1:1) = shape(struct_module_c_i)
    dloc = loc(struct_module_c_i)
end subroutine f90wrap_struct_module__array__c_i

subroutine f90wrap_struct_module__array__Odet(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use struct_module, only: struct_module_odet => odet
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 5
    dshape(1:1) = shape(struct_module_Odet)
    dloc = loc(struct_module_Odet)
end subroutine f90wrap_struct_module__array__Odet

! End of module struct_module defined in file Struct_module.fpp

! Module pspot_module defined in file Struct_module.fpp

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

subroutine f90wrap_pspot__get__elename(this, f90wrap_elename)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in)   :: this(2)
    type(pspot_ptr_type) :: this_ptr
    character(3), intent(out) :: f90wrap_elename

    this_ptr = transfer(this, this_ptr)
    f90wrap_elename = this_ptr%p%elename
end subroutine f90wrap_pspot__get__elename

subroutine f90wrap_pspot__set__elename(this, f90wrap_elename)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in)   :: this(2)
    type(pspot_ptr_type) :: this_ptr
    character(3), intent(in) :: f90wrap_elename

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%elename = f90wrap_elename
end subroutine f90wrap_pspot__set__elename

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

subroutine f90wrap_pspot__array__r(this, nd, dtype, dshape, dloc)
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
    if (allocated(this_ptr%p%r)) then
        dshape(1:1) = shape(this_ptr%p%r)
        dloc = loc(this_ptr%p%r)
    else
        dloc = 0
    end if
end subroutine f90wrap_pspot__array__r

subroutine f90wrap_pspot__array__rab(this, nd, dtype, dshape, dloc)
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
    if (allocated(this_ptr%p%rab)) then
        dshape(1:1) = shape(this_ptr%p%rab)
        dloc = loc(this_ptr%p%rab)
    else
        dloc = 0
    end if
end subroutine f90wrap_pspot__array__rab

subroutine f90wrap_pspot__array__Vlocr(this, nd, dtype, dshape, dloc)
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
    if (allocated(this_ptr%p%Vlocr)) then
        dshape(1:1) = shape(this_ptr%p%Vlocr)
        dloc = loc(this_ptr%p%Vlocr)
    else
        dloc = 0
    end if
end subroutine f90wrap_pspot__array__Vlocr

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

subroutine f90wrap_pspot__array__indx(this, nd, dtype, dshape, dloc)
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
    if (allocated(this_ptr%p%indx)) then
        dshape(1:1) = shape(this_ptr%p%indx)
        dloc = loc(this_ptr%p%indx)
    else
        dloc = 0
    end if
end subroutine f90wrap_pspot__array__indx

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

subroutine f90wrap_pspot__array__Dij(this, nd, dtype, dshape, dloc)
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
    if (allocated(this_ptr%p%Dij)) then
        dshape(1:2) = shape(this_ptr%p%Dij)
        dloc = loc(this_ptr%p%Dij)
    else
        dloc = 0
    end if
end subroutine f90wrap_pspot__array__Dij

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

subroutine f90wrap_pspot__get__lden(this, f90wrap_lden)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in)   :: this(2)
    type(pspot_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_lden

    this_ptr = transfer(this, this_ptr)
    f90wrap_lden = this_ptr%p%lden
end subroutine f90wrap_pspot__get__lden

subroutine f90wrap_pspot__set__lden(this, f90wrap_lden)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in)   :: this(2)
    type(pspot_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_lden

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%lden = f90wrap_lden
end subroutine f90wrap_pspot__set__lden

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

subroutine f90wrap_pspot__array__dden_dr(this, nd, dtype, dshape, dloc)
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
    if (allocated(this_ptr%p%dden_dr)) then
        dshape(1:1) = shape(this_ptr%p%dden_dr)
        dloc = loc(this_ptr%p%dden_dr)
    else
        dloc = 0
    end if
end subroutine f90wrap_pspot__array__dden_dr

subroutine f90wrap_pspot__get__lcore(this, f90wrap_lcore)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in)   :: this(2)
    type(pspot_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_lcore

    this_ptr = transfer(this, this_ptr)
    f90wrap_lcore = this_ptr%p%lcore
end subroutine f90wrap_pspot__get__lcore

subroutine f90wrap_pspot__set__lcore(this, f90wrap_lcore)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in)   :: this(2)
    type(pspot_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_lcore

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%lcore = f90wrap_lcore
end subroutine f90wrap_pspot__set__lcore

subroutine f90wrap_pspot__array__denc(this, nd, dtype, dshape, dloc)
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
    if (allocated(this_ptr%p%denc)) then
        dshape(1:1) = shape(this_ptr%p%denc)
        dloc = loc(this_ptr%p%denc)
    else
        dloc = 0
    end if
end subroutine f90wrap_pspot__array__denc

subroutine f90wrap_pspot__array__ddenc_dr(this, nd, dtype, dshape, dloc)
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
    if (allocated(this_ptr%p%ddenc_dr)) then
        dshape(1:1) = shape(this_ptr%p%ddenc_dr)
        dloc = loc(this_ptr%p%ddenc_dr)
    else
        dloc = 0
    end if
end subroutine f90wrap_pspot__array__ddenc_dr

subroutine f90wrap_pspot__get__rnoverlap(this, f90wrap_rnoverlap)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in)   :: this(2)
    type(pspot_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_rnoverlap

    this_ptr = transfer(this, this_ptr)
    f90wrap_rnoverlap = this_ptr%p%rnoverlap
end subroutine f90wrap_pspot__get__rnoverlap

subroutine f90wrap_pspot__set__rnoverlap(this, f90wrap_rnoverlap)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in)   :: this(2)
    type(pspot_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_rnoverlap

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%rnoverlap = f90wrap_rnoverlap
end subroutine f90wrap_pspot__set__rnoverlap

subroutine f90wrap_pspot__get__eps(this, f90wrap_eps)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in)   :: this(2)
    type(pspot_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_eps

    this_ptr = transfer(this, this_ptr)
    f90wrap_eps = this_ptr%p%eps
end subroutine f90wrap_pspot__get__eps

subroutine f90wrap_pspot__set__eps(this, f90wrap_eps)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in)   :: this(2)
    type(pspot_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_eps

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%eps = f90wrap_eps
end subroutine f90wrap_pspot__set__eps

subroutine f90wrap_pspot__get__eae(this, f90wrap_eae)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in)   :: this(2)
    type(pspot_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_eae

    this_ptr = transfer(this, this_ptr)
    f90wrap_eae = this_ptr%p%eae
end subroutine f90wrap_pspot__get__eae

subroutine f90wrap_pspot__set__eae(this, f90wrap_eae)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in)   :: this(2)
    type(pspot_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_eae

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%eae = f90wrap_eae
end subroutine f90wrap_pspot__set__eae

subroutine f90wrap_pspot__get__nwfa(this, f90wrap_nwfa)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in)   :: this(2)
    type(pspot_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_nwfa

    this_ptr = transfer(this, this_ptr)
    f90wrap_nwfa = this_ptr%p%nwfa
end subroutine f90wrap_pspot__get__nwfa

subroutine f90wrap_pspot__set__nwfa(this, f90wrap_nwfa)
    use pspot_module, only: pspot
    implicit none
    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer, intent(in)   :: this(2)
    type(pspot_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_nwfa

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nwfa = f90wrap_nwfa
end subroutine f90wrap_pspot__set__nwfa

subroutine f90wrap_pspot__array__wfal(this, nd, dtype, dshape, dloc)
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
    if (allocated(this_ptr%p%wfal)) then
        dshape(1:1) = shape(this_ptr%p%wfal)
        dloc = loc(this_ptr%p%wfal)
    else
        dloc = 0
    end if
end subroutine f90wrap_pspot__array__wfal

subroutine f90wrap_pspot__array__wfar(this, nd, dtype, dshape, dloc)
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
    if (allocated(this_ptr%p%wfar)) then
        dshape(1:2) = shape(this_ptr%p%wfar)
        dloc = loc(this_ptr%p%wfar)
    else
        dloc = 0
    end if
end subroutine f90wrap_pspot__array__wfar

subroutine f90wrap_pspot_initialise(this)
    use pspot_module, only: pspot
    implicit none

    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    type(pspot_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate (this_ptr%p)
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
    deallocate (this_ptr%p)
end subroutine f90wrap_pspot_finalise

subroutine f90wrap_attribute__get__value(this, f90wrap_value)
    use pspot_module, only: attribute
    implicit none
    type attribute_ptr_type
        type(attribute), pointer :: p => NULL()
    end type attribute_ptr_type
    integer, intent(in)   :: this(2)
    type(attribute_ptr_type) :: this_ptr
    character(150), intent(out) :: f90wrap_value

    this_ptr = transfer(this, this_ptr)
    f90wrap_value = this_ptr%p%value
end subroutine f90wrap_attribute__get__value

subroutine f90wrap_attribute__set__value(this, f90wrap_value)
    use pspot_module, only: attribute
    implicit none
    type attribute_ptr_type
        type(attribute), pointer :: p => NULL()
    end type attribute_ptr_type
    integer, intent(in)   :: this(2)
    type(attribute_ptr_type) :: this_ptr
    character(150), intent(in) :: f90wrap_value

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%value = f90wrap_value
end subroutine f90wrap_attribute__set__value

subroutine f90wrap_attribute_initialise(this)
    use pspot_module, only: attribute
    implicit none

    type attribute_ptr_type
        type(attribute), pointer :: p => NULL()
    end type attribute_ptr_type
    type(attribute_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate (this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_attribute_initialise

subroutine f90wrap_attribute_finalise(this)
    use pspot_module, only: attribute
    implicit none

    type attribute_ptr_type
        type(attribute), pointer :: p => NULL()
    end type attribute_ptr_type
    type(attribute_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate (this_ptr%p)
end subroutine f90wrap_attribute_finalise

subroutine f90wrap_read_pspot(nty, filenames, n0)
    use pspot_module, only: read_pspot
    implicit none

    integer(4), intent(in) :: nty
    character(30), intent(in), dimension(n0) :: filenames
    integer :: n0
    !f2py intent(hide), depend(filenames) :: n0 = shape(filenames,0)
    call read_pspot(nty=nty, filenames=filenames)
end subroutine f90wrap_read_pspot

subroutine f90wrap_read_psupf_atom(ity, filename, ps)
    use pspot_module, only: pspot, read_psupf_atom
    implicit none

    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer(4), intent(in) :: ity
    character(30), intent(in) :: filename
    type(pspot_ptr_type) :: ps_ptr
    integer, intent(in), dimension(2) :: ps
    ps_ptr = transfer(ps, ps_ptr)
    call read_psupf_atom(Ity=ity, filename=filename, ps=ps_ptr%p)
end subroutine f90wrap_read_psupf_atom

subroutine f90wrap_scan_head(file_unit, title, start_old)
    use pspot_module, only: scan_head
    implicit none

    integer(4) :: file_unit
    character*(*) :: title
    logical :: start_old
    call scan_head(file_unit=file_unit, title=title, start_old=start_old)
end subroutine f90wrap_scan_head

subroutine f90wrap_scan_tail(file_unit, title)
    use pspot_module, only: scan_tail
    implicit none

    integer(4) :: file_unit
    character*(*) :: title
    call scan_tail(file_unit=file_unit, title=title)
end subroutine f90wrap_scan_tail

subroutine f90wrap_read_pseudo_header(ps, nproj)
    use pspot_module, only: pspot, read_pseudo_header
    implicit none

    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    type(pspot_ptr_type) :: ps_ptr
    integer, intent(out), dimension(2) :: ps
    integer(4), intent(out) :: nproj
    allocate (ps_ptr%p)
    call read_pseudo_header(ps=ps_ptr%p, nproj=nproj)
    ps = transfer(ps_ptr, ps)
end subroutine f90wrap_read_pseudo_header

subroutine f90wrap_exist_in(string1, ret_exist_in, string2)
    use pspot_module, only: exist_in
    implicit none

    character*(*) :: string1
    logical, intent(out) :: ret_exist_in
    character*(*) :: string2
    ret_exist_in = exist_in(string1=string1, string2=string2)
end subroutine f90wrap_exist_in

subroutine f90wrap_exist_ibegin(string1, ret_exist_ibegin, string2)
    use pspot_module, only: exist_ibegin
    implicit none

    character*(*) :: string1
    integer(4), intent(out) :: ret_exist_ibegin
    character*(*) :: string2
    ret_exist_ibegin = exist_ibegin(string1=string1, string2=string2)
end subroutine f90wrap_exist_ibegin

subroutine f90wrap_read_pseudo_nonlocal(unit_upf, nl, beta_r, d0, rcut, proj_l, n0, n1, n2, n3, n4)
    use pspot_module, only: read_pseudo_nonlocal
    implicit none

    integer(4) :: unit_upf
    integer(4) :: nl
    real(8), dimension(n0, n1) :: beta_r
    real(8), dimension(n2, n3) :: d0
    real(8) :: rcut
    integer(4), dimension(n4) :: proj_l
    integer :: n0
    !f2py intent(hide), depend(beta_r) :: n0 = shape(beta_r,0)
    integer :: n1
    !f2py intent(hide), depend(beta_r) :: n1 = shape(beta_r,1)
    integer :: n2
    !f2py intent(hide), depend(d0) :: n2 = shape(d0,0)
    integer :: n3
    !f2py intent(hide), depend(d0) :: n3 = shape(d0,1)
    integer :: n4
    !f2py intent(hide), depend(proj_l) :: n4 = shape(proj_l,0)
    call read_pseudo_nonlocal(unit_UPF=unit_upf, nl=nl, beta_r=beta_r, D0=d0, rcut=rcut, proj_l=proj_l)
end subroutine f90wrap_read_pseudo_nonlocal

subroutine f90wrap_read_pseudo_pswfc(unit_upf, nwfc, nps, wfcl, wfcr, n0, n1, n2)
    use pspot_module, only: read_pseudo_pswfc
    implicit none

    integer(4), intent(in) :: unit_upf
    integer(4), intent(in) :: nwfc
    integer(4), intent(in) :: nps
    integer(4), intent(inout), dimension(n0) :: wfcl
    real(8), intent(inout), dimension(n1, n2) :: wfcr
    integer :: n0
    !f2py intent(hide), depend(wfcl) :: n0 = shape(wfcl,0)
    integer :: n1
    !f2py intent(hide), depend(wfcr) :: n1 = shape(wfcr,0)
    integer :: n2
    !f2py intent(hide), depend(wfcr) :: n2 = shape(wfcr,1)
    call read_pseudo_pswfc(unit_UPF=unit_upf, nwfc=nwfc, nps=nps, wfcl=wfcl, wfcr=wfcr)
end subroutine f90wrap_read_pseudo_pswfc

subroutine f90wrap_aep_generator(nz, ps)
    use pspot_module, only: aep_generator, pspot
    implicit none

    type pspot_ptr_type
        type(pspot), pointer :: p => NULL()
    end type pspot_ptr_type
    integer(4), intent(in) :: nz
    type(pspot_ptr_type) :: ps_ptr
    integer, intent(out), dimension(2) :: ps
    allocate (ps_ptr%p)
    call aep_generator(nz=nz, ps=ps_ptr%p)
    ps = transfer(ps_ptr, ps)
end subroutine f90wrap_aep_generator

subroutine f90wrap_get_value_int(char_in, char_find, variable, find_flag)
    use pspot_module, only: get_value
    implicit none

    character(120), intent(in) :: char_in
    character*(*), intent(in) :: char_find
    integer(4) :: variable
    logical :: find_flag
    call get_value(char_in=char_in, char_find=char_find, variable=variable, find_flag=find_flag)
end subroutine f90wrap_get_value_int

subroutine f90wrap_get_value_real(char_in, char_find, variable, find_flag)
    use pspot_module, only: get_value
    implicit none

    character(120), intent(in) :: char_in
    character*(*), intent(in) :: char_find
    real(8) :: variable
    logical :: find_flag
    call get_value(char_in=char_in, char_find=char_find, variable=variable, find_flag=find_flag)
end subroutine f90wrap_get_value_real

subroutine f90wrap_get_value_char(char_in, char_find, variable, find_flag)
    use pspot_module, only: get_value
    implicit none

    character(120), intent(in) :: char_in
    character*(*), intent(in) :: char_find
    character*(*) :: variable
    logical :: find_flag
    call get_value(char_in=char_in, char_find=char_find, variable=variable, find_flag=find_flag)
end subroutine f90wrap_get_value_char

subroutine f90wrap_get_value_logic(char_in, char_find, variable, find_flag)
    use pspot_module, only: get_value
    implicit none

    character(120), intent(in) :: char_in
    character*(*), intent(in) :: char_find
    logical :: variable
    logical :: find_flag
    call get_value(char_in=char_in, char_find=char_find, variable=variable, find_flag=find_flag)
end subroutine f90wrap_get_value_logic

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

subroutine f90wrap_pspot_module__get__max_nwfa(f90wrap_max_nwfa)
    use pspot_module, only: pspot_module_max_nwfa => max_nwfa
    implicit none
    integer(4), intent(out) :: f90wrap_max_nwfa

    f90wrap_max_nwfa = pspot_module_max_nwfa
end subroutine f90wrap_pspot_module__get__max_nwfa

subroutine f90wrap_pspot_module__set__max_nwfa(f90wrap_max_nwfa)
    use pspot_module, only: pspot_module_max_nwfa => max_nwfa
    implicit none
    integer(4), intent(in) :: f90wrap_max_nwfa

    pspot_module_max_nwfa = f90wrap_max_nwfa
end subroutine f90wrap_pspot_module__set__max_nwfa

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

! End of module pspot_module defined in file Struct_module.fpp

! Module grid_module defined in file Struct_module.fpp

subroutine f90wrap_grid_type__array__rhoS(this, nd, dtype, dshape, dloc)
    use grid_module, only: grid_type
    implicit none
    type grid_type_ptr_type
        type(grid_type), pointer :: p => NULL()
    end type grid_type_ptr_type
    integer, intent(in) :: this(2)
    type(grid_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%rhoS)) then
        dshape(1:2) = shape(this_ptr%p%rhoS)
        dloc = loc(this_ptr%p%rhoS)
    else
        dloc = 0
    end if
end subroutine f90wrap_grid_type__array__rhoS

subroutine f90wrap_grid_type__array__rho(this, nd, dtype, dshape, dloc)
    use grid_module, only: grid_type
    implicit none
    type grid_type_ptr_type
        type(grid_type), pointer :: p => NULL()
    end type grid_type_ptr_type
    integer, intent(in) :: this(2)
    type(grid_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%rho)) then
        dshape(1:1) = shape(this_ptr%p%rho)
        dloc = loc(this_ptr%p%rho)
    else
        dloc = 0
    end if
end subroutine f90wrap_grid_type__array__rho

subroutine f90wrap_grid_type__array__vxcS(this, nd, dtype, dshape, dloc)
    use grid_module, only: grid_type
    implicit none
    type grid_type_ptr_type
        type(grid_type), pointer :: p => NULL()
    end type grid_type_ptr_type
    integer, intent(in) :: this(2)
    type(grid_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%vxcS)) then
        dshape(1:2) = shape(this_ptr%p%vxcS)
        dloc = loc(this_ptr%p%vxcS)
    else
        dloc = 0
    end if
end subroutine f90wrap_grid_type__array__vxcS

subroutine f90wrap_grid_type__array__vhxcd(this, nd, dtype, dshape, dloc)
    use grid_module, only: grid_type
    implicit none
    type grid_type_ptr_type
        type(grid_type), pointer :: p => NULL()
    end type grid_type_ptr_type
    integer, intent(in) :: this(2)
    type(grid_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%vhxcd)) then
        dshape(1:2) = shape(this_ptr%p%vhxcd)
        dloc = loc(this_ptr%p%vhxcd)
    else
        dloc = 0
    end if
end subroutine f90wrap_grid_type__array__vhxcd

subroutine f90wrap_grid_type__array__vlpp(this, nd, dtype, dshape, dloc)
    use grid_module, only: grid_type
    implicit none
    type grid_type_ptr_type
        type(grid_type), pointer :: p => NULL()
    end type grid_type_ptr_type
    integer, intent(in) :: this(2)
    type(grid_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%vlpp)) then
        dshape(1:1) = shape(this_ptr%p%vlpp)
        dloc = loc(this_ptr%p%vlpp)
    else
        dloc = 0
    end if
end subroutine f90wrap_grid_type__array__vlpp

subroutine f90wrap_grid_type__array__vh(this, nd, dtype, dshape, dloc)
    use grid_module, only: grid_type
    implicit none
    type grid_type_ptr_type
        type(grid_type), pointer :: p => NULL()
    end type grid_type_ptr_type
    integer, intent(in) :: this(2)
    type(grid_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%vh)) then
        dshape(1:1) = shape(this_ptr%p%vh)
        dloc = loc(this_ptr%p%vh)
    else
        dloc = 0
    end if
end subroutine f90wrap_grid_type__array__vh

subroutine f90wrap_grid_type__array__eval(this, nd, dtype, dshape, dloc)
    use grid_module, only: grid_type
    implicit none
    type grid_type_ptr_type
        type(grid_type), pointer :: p => NULL()
    end type grid_type_ptr_type
    integer, intent(in) :: this(2)
    type(grid_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 3
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%eval)) then
        dshape(1:3) = shape(this_ptr%p%eval)
        dloc = loc(this_ptr%p%eval)
    else
        dloc = 0
    end if
end subroutine f90wrap_grid_type__array__eval

subroutine f90wrap_grid_type__array__rhoc(this, nd, dtype, dshape, dloc)
    use grid_module, only: grid_type
    implicit none
    type grid_type_ptr_type
        type(grid_type), pointer :: p => NULL()
    end type grid_type_ptr_type
    integer, intent(in) :: this(2)
    type(grid_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%rhoc)) then
        dshape(1:1) = shape(this_ptr%p%rhoc)
        dloc = loc(this_ptr%p%rhoc)
    else
        dloc = 0
    end if
end subroutine f90wrap_grid_type__array__rhoc

subroutine f90wrap_grid_type__array__gVec(this, nd, dtype, dshape, dloc)
    use grid_module, only: grid_type
    implicit none
    type grid_type_ptr_type
        type(grid_type), pointer :: p => NULL()
    end type grid_type_ptr_type
    integer, intent(in) :: this(2)
    type(grid_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%gVec)) then
        dshape(1:2) = shape(this_ptr%p%gVec)
        dloc = loc(this_ptr%p%gVec)
    else
        dloc = 0
    end if
end subroutine f90wrap_grid_type__array__gVec

subroutine f90wrap_grid_type__array__gMask(this, nd, dtype, dshape, dloc)
    use grid_module, only: grid_type
    implicit none
    type grid_type_ptr_type
        type(grid_type), pointer :: p => NULL()
    end type grid_type_ptr_type
    integer, intent(in) :: this(2)
    type(grid_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%gMask)) then
        dshape(1:1) = shape(this_ptr%p%gMask)
        dloc = loc(this_ptr%p%gMask)
    else
        dloc = 0
    end if
end subroutine f90wrap_grid_type__array__gMask

subroutine f90wrap_grid_type__array__rVec(this, nd, dtype, dshape, dloc)
    use grid_module, only: grid_type
    implicit none
    type grid_type_ptr_type
        type(grid_type), pointer :: p => NULL()
    end type grid_type_ptr_type
    integer, intent(in) :: this(2)
    type(grid_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%rVec)) then
        dshape(1:2) = shape(this_ptr%p%rVec)
        dloc = loc(this_ptr%p%rVec)
    else
        dloc = 0
    end if
end subroutine f90wrap_grid_type__array__rVec

subroutine f90wrap_grid_type__array__lsp(this, nd, dtype, dshape, dloc)
    use grid_module, only: grid_type
    implicit none
    type grid_type_ptr_type
        type(grid_type), pointer :: p => NULL()
    end type grid_type_ptr_type
    integer, intent(in) :: this(2)
    type(grid_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 3
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%lsp)) then
        dshape(1:3) = shape(this_ptr%p%lsp)
        dloc = loc(this_ptr%p%lsp)
    else
        dloc = 0
    end if
end subroutine f90wrap_grid_type__array__lsp

subroutine f90wrap_grid_type__get__oneDlength(this, f90wrap_oneDlength)
    use grid_module, only: grid_type
    implicit none
    type grid_type_ptr_type
        type(grid_type), pointer :: p => NULL()
    end type grid_type_ptr_type
    integer, intent(in)   :: this(2)
    type(grid_type_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_oneDlength

    this_ptr = transfer(this, this_ptr)
    f90wrap_oneDlength = this_ptr%p%oneDlength
end subroutine f90wrap_grid_type__get__oneDlength

subroutine f90wrap_grid_type__set__oneDlength(this, f90wrap_oneDlength)
    use grid_module, only: grid_type
    implicit none
    type grid_type_ptr_type
        type(grid_type), pointer :: p => NULL()
    end type grid_type_ptr_type
    integer, intent(in)   :: this(2)
    type(grid_type_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_oneDlength

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%oneDlength = f90wrap_oneDlength
end subroutine f90wrap_grid_type__set__oneDlength

subroutine f90wrap_grid_type_initialise(this)
    use grid_module, only: grid_type
    implicit none

    type grid_type_ptr_type
        type(grid_type), pointer :: p => NULL()
    end type grid_type_ptr_type
    type(grid_type_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate (this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_grid_type_initialise

subroutine f90wrap_grid_type_finalise(this)
    use grid_module, only: grid_type
    implicit none

    type grid_type_ptr_type
        type(grid_type), pointer :: p => NULL()
    end type grid_type_ptr_type
    type(grid_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate (this_ptr%p)
end subroutine f90wrap_grid_type_finalise

subroutine f90wrap_kgrid_type__array__vec(this, nd, dtype, dshape, dloc)
    use grid_module, only: kgrid_type
    implicit none
    type kgrid_type_ptr_type
        type(kgrid_type), pointer :: p => NULL()
    end type kgrid_type_ptr_type
    integer, intent(in) :: this(2)
    type(kgrid_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%vec)) then
        dshape(1:2) = shape(this_ptr%p%vec)
        dloc = loc(this_ptr%p%vec)
    else
        dloc = 0
    end if
end subroutine f90wrap_kgrid_type__array__vec

subroutine f90wrap_kgrid_type__array__vcar(this, nd, dtype, dshape, dloc)
    use grid_module, only: kgrid_type
    implicit none
    type kgrid_type_ptr_type
        type(kgrid_type), pointer :: p => NULL()
    end type kgrid_type_ptr_type
    integer, intent(in) :: this(2)
    type(kgrid_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%vcar)) then
        dshape(1:2) = shape(this_ptr%p%vcar)
        dloc = loc(this_ptr%p%vcar)
    else
        dloc = 0
    end if
end subroutine f90wrap_kgrid_type__array__vcar

subroutine f90wrap_kgrid_type__array__wk(this, nd, dtype, dshape, dloc)
    use grid_module, only: kgrid_type
    implicit none
    type kgrid_type_ptr_type
        type(kgrid_type), pointer :: p => NULL()
    end type kgrid_type_ptr_type
    integer, intent(in) :: this(2)
    type(kgrid_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%wk)) then
        dshape(1:1) = shape(this_ptr%p%wk)
        dloc = loc(this_ptr%p%wk)
    else
        dloc = 0
    end if
end subroutine f90wrap_kgrid_type__array__wk

subroutine f90wrap_kgrid_type_initialise(this)
    use grid_module, only: kgrid_type
    implicit none

    type kgrid_type_ptr_type
        type(kgrid_type), pointer :: p => NULL()
    end type kgrid_type_ptr_type
    type(kgrid_type_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate (this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_kgrid_type_initialise

subroutine f90wrap_kgrid_type_finalise(this)
    use grid_module, only: kgrid_type
    implicit none

    type kgrid_type_ptr_type
        type(kgrid_type), pointer :: p => NULL()
    end type kgrid_type_ptr_type
    type(kgrid_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate (this_ptr%p)
end subroutine f90wrap_kgrid_type_finalise

subroutine f90wrap_eigen_type__array__wvf(this, nd, dtype, dshape, dloc)
    use grid_module, only: eigen_type
    implicit none
    type eigen_type_ptr_type
        type(eigen_type), pointer :: p => NULL()
    end type eigen_type_ptr_type
    integer, intent(in) :: this(2)
    type(eigen_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 4
    dtype = 15
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%wvf)) then
        dshape(1:4) = shape(this_ptr%p%wvf)
        dloc = loc(this_ptr%p%wvf)
    else
        dloc = 0
    end if
end subroutine f90wrap_eigen_type__array__wvf

subroutine f90wrap_eigen_type__array__val(this, nd, dtype, dshape, dloc)
    use grid_module, only: eigen_type
    implicit none
    type eigen_type_ptr_type
        type(eigen_type), pointer :: p => NULL()
    end type eigen_type_ptr_type
    integer, intent(in) :: this(2)
    type(eigen_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 3
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%val)) then
        dshape(1:3) = shape(this_ptr%p%val)
        dloc = loc(this_ptr%p%val)
    else
        dloc = 0
    end if
end subroutine f90wrap_eigen_type__array__val

subroutine f90wrap_eigen_type__array__wvfG(this, nd, dtype, dshape, dloc)
    use grid_module, only: eigen_type
    implicit none
    type eigen_type_ptr_type
        type(eigen_type), pointer :: p => NULL()
    end type eigen_type_ptr_type
    integer, intent(in) :: this(2)
    type(eigen_type_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 3
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%wvfG)) then
        dshape(1:3) = shape(this_ptr%p%wvfG)
        dloc = loc(this_ptr%p%wvfG)
    else
        dloc = 0
    end if
end subroutine f90wrap_eigen_type__array__wvfG

subroutine f90wrap_eigen_type_initialise(this)
    use grid_module, only: eigen_type
    implicit none

    type eigen_type_ptr_type
        type(eigen_type), pointer :: p => NULL()
    end type eigen_type_ptr_type
    type(eigen_type_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate (this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_eigen_type_initialise

subroutine f90wrap_eigen_type_finalise(this)
    use grid_module, only: eigen_type
    implicit none

    type eigen_type_ptr_type
        type(eigen_type), pointer :: p => NULL()
    end type eigen_type_ptr_type
    type(eigen_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate (this_ptr%p)
end subroutine f90wrap_eigen_type_finalise

subroutine f90wrap_build_rgrid
    use grid_module, only: build_rgrid
    implicit none

    call build_rgrid()
end subroutine f90wrap_build_rgrid

subroutine f90wrap_destroy_rgrid
    use grid_module, only: destroy_rgrid
    implicit none

    call destroy_rgrid()
end subroutine f90wrap_destroy_rgrid

subroutine f90wrap_build_kgrid
    use grid_module, only: build_kgrid
    implicit none

    call build_kgrid()
end subroutine f90wrap_build_kgrid

subroutine f90wrap_destroy_kpt
    use grid_module, only: destroy_kpt
    implicit none

    call destroy_kpt()
end subroutine f90wrap_destroy_kpt

subroutine f90wrap_build_eigen
    use grid_module, only: build_eigen
    implicit none

    call build_eigen()
end subroutine f90wrap_build_eigen

subroutine f90wrap_destroy_eigen
    use grid_module, only: destroy_eigen
    implicit none

    call destroy_eigen()
end subroutine f90wrap_destroy_eigen

subroutine f90wrap_fillqtable
    use grid_module, only: fillqtable
    implicit none

    call fillqtable()
end subroutine f90wrap_fillqtable

subroutine f90wrap_fillrtable
    use grid_module, only: fillrtable
    implicit none

    call fillrtable()
end subroutine f90wrap_fillrtable

subroutine f90wrap_symm_kgrid
    use grid_module, only: symm_kgrid
    implicit none

    call symm_kgrid()
end subroutine f90wrap_symm_kgrid

subroutine f90wrap_symm_density(rho, n0, n1, n2)
    use grid_module, only: symm_density
    implicit none

    real(8), intent(inout), dimension(n0, n1, n2) :: rho
    integer :: n0
    !f2py intent(hide), depend(rho) :: n0 = shape(rho,0)
    integer :: n1
    !f2py intent(hide), depend(rho) :: n1 = shape(rho,1)
    integer :: n2
    !f2py intent(hide), depend(rho) :: n2 = shape(rho,2)
    call symm_density(rho=rho)
end subroutine f90wrap_symm_density

subroutine f90wrap_isymm_apply(isy, nsize, vin, vout, lsymm)
    use grid_module, only: isymm_apply
    implicit none

    integer(4), intent(in) :: isy
    integer(4), intent(in), dimension(3) :: nsize
    integer(4), intent(in), dimension(3) :: vin
    integer(4), dimension(3), intent(inout) :: vout
    logical, intent(out) :: lsymm
    call isymm_apply(Isy=isy, nsize=nsize, vin=vin, vout=vout, lsymm=lsymm)
end subroutine f90wrap_isymm_apply

subroutine f90wrap_build_parallel_3d_grid
    use grid_module, only: build_parallel_3d_grid
    implicit none

    call build_parallel_3d_grid()
end subroutine f90wrap_build_parallel_3d_grid

subroutine f90wrap_build_parallel_sph_grid
    use grid_module, only: build_parallel_sph_grid
    implicit none

    call build_parallel_sph_grid()
end subroutine f90wrap_build_parallel_sph_grid

subroutine f90wrap_grid_split_sph(ngrid, ncore, comm, id, grid_range, recvcounts, displs, gridrange_sum, n1, n2, n3, n, &
                                  n0, n11, n22)
    use grid_module, only: grid_split_sph
    implicit none

    integer(4), intent(in) :: ngrid
    integer(4), intent(in) :: ncore
    integer(4), intent(in) :: comm
    integer(4), intent(in) :: id
    integer(4), dimension(3), intent(inout) :: grid_range
    integer(4), intent(inout), dimension(n0) :: recvcounts
    integer(4), intent(inout), dimension(n1) :: displs
    integer(4), intent(inout), dimension(3, n2) :: gridrange_sum
    integer(4), intent(inout) :: n1
    integer(4), intent(inout) :: n2
    integer(4), intent(inout) :: n3
    integer(4), intent(inout) :: n
    integer :: n0
!f2py intent(hide), depend(recvcounts) :: n0 = shape(recvcounts,0)
    integer :: n11
!f2py intent(hide), depend(displs) :: n11 = shape(displs,0)
    integer :: n22
!f2py intent(hide), depend(gridrange_sum) :: n22 = shape(gridrange_sum,1)
    call grid_split_sph(ngrid=ngrid, ncore=ncore, comm=comm, id=id, grid_range=grid_range, recvcounts=recvcounts, &
                        displs=displs, gridrange_sum=gridrange_sum, n1=n1, n2=n2, n3=n3, n=n)
end subroutine f90wrap_grid_split_sph

subroutine f90wrap_sphere_region(n1, n2, n3, lsphere, nz_map, nsphere, n0, n11, n22, n33)
    use grid_module, only: sphere_region
    implicit none

    integer(4) :: n1
    integer(4) :: n2
    integer(4) :: n3
    logical, dimension(n0, n1, n2) :: lsphere
    integer(4), dimension(2, n3) :: nz_map
    integer(4) :: nsphere
    integer :: n0
    !f2py intent(hide), depend(lsphere) :: n0 = shape(lsphere,0)
    integer :: n11
    !f2py intent(hide), depend(lsphere) :: n11 = shape(lsphere,1)
    integer :: n22
    !f2py intent(hide), depend(lsphere) :: n22 = shape(lsphere,2)
    integer :: n33
    !f2py intent(hide), depend(nz_map) :: n33 = shape(nz_map,1)
    call sphere_region(n1=n1, n2=n2, n3=n3, Lsphere=lsphere, nz_map=nz_map, nsphere=nsphere)
end subroutine f90wrap_sphere_region

subroutine f90wrap_sphere2cubic(nps, f1d, f3d, rfill, n0, n1, n2, n3)
    use grid_module, only: sphere2cubic
    implicit none

    integer(4), intent(in) :: nps
    real(8), intent(in), dimension(n0) :: f1d
    real(8), intent(inout), dimension(n1, n2, n3) :: f3d
    real(8), optional :: rfill
    integer :: n0
    !f2py intent(hide), depend(f1d) :: n0 = shape(f1d,0)
    integer :: n1
    !f2py intent(hide), depend(f3d) :: n1 = shape(f3d,0)
    integer :: n2
    !f2py intent(hide), depend(f3d) :: n2 = shape(f3d,1)
    integer :: n3
    !f2py intent(hide), depend(f3d) :: n3 = shape(f3d,2)
    call sphere2cubic(nps=nps, f1d=f1d, f3d=f3d, rfill=rfill)
end subroutine f90wrap_sphere2cubic

subroutine f90wrap_cubic2sphere(nps, f3d, f1d, n0, n1, n2, n3)
    use grid_module, only: cubic2sphere
    implicit none

    integer(4), intent(in) :: nps
    real(8), intent(in), dimension(n0, n1, n2) :: f3d
    real(8), intent(inout), dimension(n3) :: f1d
    integer :: n0
    !f2py intent(hide), depend(f3d) :: n0 = shape(f3d,0)
    integer :: n1
    !f2py intent(hide), depend(f3d) :: n1 = shape(f3d,1)
    integer :: n2
    !f2py intent(hide), depend(f3d) :: n2 = shape(f3d,2)
    integer :: n3
    !f2py intent(hide), depend(f1d) :: n3 = shape(f1d,0)
    call cubic2sphere(nps=nps, f3d=f3d, f1d=f1d)
end subroutine f90wrap_cubic2sphere

subroutine f90wrap_cubic2sphere_fft(nps, f3d, f1d, shiftn, shiftz, n0, n1, n2, n3)
    use grid_module, only: cubic2sphere_fft
    implicit none

    integer(4), intent(in) :: nps
    real(8), intent(in), dimension(n0, n1, n2) :: f3d
    real(8), intent(inout), dimension(n3) :: f1d
    integer(4) :: shiftn
    integer(4) :: shiftz
    integer :: n0
    !f2py intent(hide), depend(f3d) :: n0 = shape(f3d,0)
    integer :: n1
    !f2py intent(hide), depend(f3d) :: n1 = shape(f3d,1)
    integer :: n2
    !f2py intent(hide), depend(f3d) :: n2 = shape(f3d,2)
    integer :: n3
    !f2py intent(hide), depend(f1d) :: n3 = shape(f1d,0)
    call cubic2sphere_fft(nps=nps, f3d=f3d, f1d=f1d, shiftn=shiftn, shiftz=shiftz)
end subroutine f90wrap_cubic2sphere_fft

subroutine f90wrap_sphere2cubic_fft(nps, f1d, f3d, shiftn, shiftz, n0, n1, n2, n3)
    use grid_module, only: sphere2cubic_fft
    implicit none

    integer(4), intent(in) :: nps
    real(8), intent(in), dimension(n0) :: f1d
    real(8), intent(inout), dimension(n1, n2, n3) :: f3d
    integer(4) :: shiftn
    integer(4) :: shiftz
    integer :: n0
    !f2py intent(hide), depend(f1d) :: n0 = shape(f1d,0)
    integer :: n1
    !f2py intent(hide), depend(f3d) :: n1 = shape(f3d,0)
    integer :: n2
    !f2py intent(hide), depend(f3d) :: n2 = shape(f3d,1)
    integer :: n3
    !f2py intent(hide), depend(f3d) :: n3 = shape(f3d,2)
    call sphere2cubic_fft(nps=nps, f1d=f1d, f3d=f3d, shiftn=shiftn, shiftz=shiftz)
end subroutine f90wrap_sphere2cubic_fft

subroutine f90wrap_sumrhos(nps, rhos, rho, n0, n1, n2)
    use grid_module, only: sumrhos
    implicit none

    integer(4) :: nps
    real(8), intent(in), dimension(n0, n1) :: rhos
    real(8), intent(inout), dimension(n2) :: rho
    integer :: n0
    !f2py intent(hide), depend(rhos) :: n0 = shape(rhos,0)
    integer :: n1
    !f2py intent(hide), depend(rhos) :: n1 = shape(rhos,1)
    integer :: n2
    !f2py intent(hide), depend(rho) :: n2 = shape(rho,0)
    call sumrhos(nps=nps, rhoS=rhos, rho=rho)
end subroutine f90wrap_sumrhos

subroutine f90wrap_fft_sph_r2c(ret_array_c, array_r, n0, n1, n2, n3)
    use grid_module, only: fft_sph
    implicit none

    complex(8), intent(out), dimension(n0, n1, n2) :: ret_array_c
    real(8), intent(in), dimension(n3) :: array_r
    integer :: n0
    integer :: n1
    integer :: n2
    integer :: n3
    !f2py intent(hide), depend(array_r) :: n3 = shape(array_r,0)
    ret_array_c = fft_sph(array_r=array_r)
end subroutine f90wrap_fft_sph_r2c

subroutine f90wrap_fft_sph_c2r(ret_array_r, array_c, n0, n1, n2, n3)
    use grid_module, only: fft_sph
    implicit none

    real(8), intent(out), dimension(n0) :: ret_array_r
    complex(8), intent(in), dimension(n1, n2, n3) :: array_c
    integer :: n0
    integer :: n1
    !f2py intent(hide), depend(array_c) :: n1 = shape(array_c,0)
    integer :: n2
    !f2py intent(hide), depend(array_c) :: n2 = shape(array_c,1)
    integer :: n3
    !f2py intent(hide), depend(array_c) :: n3 = shape(array_c,2)
    ret_array_r = fft_sph(array_c=array_c)
end subroutine f90wrap_fft_sph_c2r

subroutine f90wrap_grid_module__get__n1(f90wrap_n1)
    use grid_module, only: grid_module_n1 => n1
    implicit none
    integer(4), intent(out) :: f90wrap_n1

    f90wrap_n1 = grid_module_n1
end subroutine f90wrap_grid_module__get__n1

subroutine f90wrap_grid_module__set__n1(f90wrap_n1)
    use grid_module, only: grid_module_n1 => n1
    implicit none
    integer(4), intent(in) :: f90wrap_n1

    grid_module_n1 = f90wrap_n1
end subroutine f90wrap_grid_module__set__n1

subroutine f90wrap_grid_module__get__n2(f90wrap_n2)
    use grid_module, only: grid_module_n2 => n2
    implicit none
    integer(4), intent(out) :: f90wrap_n2

    f90wrap_n2 = grid_module_n2
end subroutine f90wrap_grid_module__get__n2

subroutine f90wrap_grid_module__set__n2(f90wrap_n2)
    use grid_module, only: grid_module_n2 => n2
    implicit none
    integer(4), intent(in) :: f90wrap_n2

    grid_module_n2 = f90wrap_n2
end subroutine f90wrap_grid_module__set__n2

subroutine f90wrap_grid_module__get__n3(f90wrap_n3)
    use grid_module, only: grid_module_n3 => n3
    implicit none
    integer(4), intent(out) :: f90wrap_n3

    f90wrap_n3 = grid_module_n3
end subroutine f90wrap_grid_module__get__n3

subroutine f90wrap_grid_module__set__n3(f90wrap_n3)
    use grid_module, only: grid_module_n3 => n3
    implicit none
    integer(4), intent(in) :: f90wrap_n3

    grid_module_n3 = f90wrap_n3
end subroutine f90wrap_grid_module__set__n3

subroutine f90wrap_grid_module__get__n(f90wrap_n)
    use grid_module, only: grid_module_n => n
    implicit none
    integer(4), intent(out) :: f90wrap_n

    f90wrap_n = grid_module_n
end subroutine f90wrap_grid_module__get__n

subroutine f90wrap_grid_module__set__n(f90wrap_n)
    use grid_module, only: grid_module_n => n
    implicit none
    integer(4), intent(in) :: f90wrap_n

    grid_module_n = f90wrap_n
end subroutine f90wrap_grid_module__set__n

subroutine f90wrap_grid_module__get__ng1(f90wrap_ng1)
    use grid_module, only: grid_module_ng1 => ng1
    implicit none
    integer(4), intent(out) :: f90wrap_ng1

    f90wrap_ng1 = grid_module_ng1
end subroutine f90wrap_grid_module__get__ng1

subroutine f90wrap_grid_module__set__ng1(f90wrap_ng1)
    use grid_module, only: grid_module_ng1 => ng1
    implicit none
    integer(4), intent(in) :: f90wrap_ng1

    grid_module_ng1 = f90wrap_ng1
end subroutine f90wrap_grid_module__set__ng1

subroutine f90wrap_grid_module__get__ng2(f90wrap_ng2)
    use grid_module, only: grid_module_ng2 => ng2
    implicit none
    integer(4), intent(out) :: f90wrap_ng2

    f90wrap_ng2 = grid_module_ng2
end subroutine f90wrap_grid_module__get__ng2

subroutine f90wrap_grid_module__set__ng2(f90wrap_ng2)
    use grid_module, only: grid_module_ng2 => ng2
    implicit none
    integer(4), intent(in) :: f90wrap_ng2

    grid_module_ng2 = f90wrap_ng2
end subroutine f90wrap_grid_module__set__ng2

subroutine f90wrap_grid_module__get__ng3(f90wrap_ng3)
    use grid_module, only: grid_module_ng3 => ng3
    implicit none
    integer(4), intent(out) :: f90wrap_ng3

    f90wrap_ng3 = grid_module_ng3
end subroutine f90wrap_grid_module__get__ng3

subroutine f90wrap_grid_module__set__ng3(f90wrap_ng3)
    use grid_module, only: grid_module_ng3 => ng3
    implicit none
    integer(4), intent(in) :: f90wrap_ng3

    grid_module_ng3 = f90wrap_ng3
end subroutine f90wrap_grid_module__set__ng3

subroutine f90wrap_grid_module__get__ng(f90wrap_ng)
    use grid_module, only: grid_module_ng => ng
    implicit none
    integer(4), intent(out) :: f90wrap_ng

    f90wrap_ng = grid_module_ng
end subroutine f90wrap_grid_module__get__ng

subroutine f90wrap_grid_module__set__ng(f90wrap_ng)
    use grid_module, only: grid_module_ng => ng
    implicit none
    integer(4), intent(in) :: f90wrap_ng

    grid_module_ng = f90wrap_ng
end subroutine f90wrap_grid_module__set__ng

subroutine f90wrap_grid_module__array__gap(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use grid_module, only: grid_module_gap => gap
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 12
    dshape(1:1) = shape(grid_module_gap)
    dloc = loc(grid_module_gap)
end subroutine f90wrap_grid_module__array__gap

subroutine f90wrap_grid_module__get__dvol(f90wrap_dvol)
    use grid_module, only: grid_module_dvol => dvol
    implicit none
    real(8), intent(out) :: f90wrap_dvol

    f90wrap_dvol = grid_module_dvol
end subroutine f90wrap_grid_module__get__dvol

subroutine f90wrap_grid_module__set__dvol(f90wrap_dvol)
    use grid_module, only: grid_module_dvol => dvol
    implicit none
    real(8), intent(in) :: f90wrap_dvol

    grid_module_dvol = f90wrap_dvol
end subroutine f90wrap_grid_module__set__dvol

subroutine f90wrap_grid_module__get__nk1(f90wrap_nk1)
    use grid_module, only: grid_module_nk1 => nk1
    implicit none
    integer(4), intent(out) :: f90wrap_nk1

    f90wrap_nk1 = grid_module_nk1
end subroutine f90wrap_grid_module__get__nk1

subroutine f90wrap_grid_module__set__nk1(f90wrap_nk1)
    use grid_module, only: grid_module_nk1 => nk1
    implicit none
    integer(4), intent(in) :: f90wrap_nk1

    grid_module_nk1 = f90wrap_nk1
end subroutine f90wrap_grid_module__set__nk1

subroutine f90wrap_grid_module__get__nk2(f90wrap_nk2)
    use grid_module, only: grid_module_nk2 => nk2
    implicit none
    integer(4), intent(out) :: f90wrap_nk2

    f90wrap_nk2 = grid_module_nk2
end subroutine f90wrap_grid_module__get__nk2

subroutine f90wrap_grid_module__set__nk2(f90wrap_nk2)
    use grid_module, only: grid_module_nk2 => nk2
    implicit none
    integer(4), intent(in) :: f90wrap_nk2

    grid_module_nk2 = f90wrap_nk2
end subroutine f90wrap_grid_module__set__nk2

subroutine f90wrap_grid_module__get__nk3(f90wrap_nk3)
    use grid_module, only: grid_module_nk3 => nk3
    implicit none
    integer(4), intent(out) :: f90wrap_nk3

    f90wrap_nk3 = grid_module_nk3
end subroutine f90wrap_grid_module__get__nk3

subroutine f90wrap_grid_module__set__nk3(f90wrap_nk3)
    use grid_module, only: grid_module_nk3 => nk3
    implicit none
    integer(4), intent(in) :: f90wrap_nk3

    grid_module_nk3 = f90wrap_nk3
end subroutine f90wrap_grid_module__set__nk3

subroutine f90wrap_grid_module__get__nk(f90wrap_nk)
    use grid_module, only: grid_module_nk => nk
    implicit none
    integer(4), intent(out) :: f90wrap_nk

    f90wrap_nk = grid_module_nk
end subroutine f90wrap_grid_module__get__nk

subroutine f90wrap_grid_module__set__nk(f90wrap_nk)
    use grid_module, only: grid_module_nk => nk
    implicit none
    integer(4), intent(in) :: f90wrap_nk

    grid_module_nk = f90wrap_nk
end subroutine f90wrap_grid_module__set__nk

subroutine f90wrap_grid_module__array__kdispl(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use grid_module, only: grid_module_kdispl => kdispl
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 12
    dshape(1:1) = shape(grid_module_kdispl)
    dloc = loc(grid_module_kdispl)
end subroutine f90wrap_grid_module__array__kdispl

subroutine f90wrap_grid_module__get__global_n1(f90wrap_global_n1)
    use grid_module, only: grid_module_global_n1 => global_n1
    implicit none
    integer(4), intent(out) :: f90wrap_global_n1

    f90wrap_global_n1 = grid_module_global_n1
end subroutine f90wrap_grid_module__get__global_n1

subroutine f90wrap_grid_module__set__global_n1(f90wrap_global_n1)
    use grid_module, only: grid_module_global_n1 => global_n1
    implicit none
    integer(4), intent(in) :: f90wrap_global_n1

    grid_module_global_n1 = f90wrap_global_n1
end subroutine f90wrap_grid_module__set__global_n1

subroutine f90wrap_grid_module__get__global_n2(f90wrap_global_n2)
    use grid_module, only: grid_module_global_n2 => global_n2
    implicit none
    integer(4), intent(out) :: f90wrap_global_n2

    f90wrap_global_n2 = grid_module_global_n2
end subroutine f90wrap_grid_module__get__global_n2

subroutine f90wrap_grid_module__set__global_n2(f90wrap_global_n2)
    use grid_module, only: grid_module_global_n2 => global_n2
    implicit none
    integer(4), intent(in) :: f90wrap_global_n2

    grid_module_global_n2 = f90wrap_global_n2
end subroutine f90wrap_grid_module__set__global_n2

subroutine f90wrap_grid_module__get__global_n3(f90wrap_global_n3)
    use grid_module, only: grid_module_global_n3 => global_n3
    implicit none
    integer(4), intent(out) :: f90wrap_global_n3

    f90wrap_global_n3 = grid_module_global_n3
end subroutine f90wrap_grid_module__get__global_n3

subroutine f90wrap_grid_module__set__global_n3(f90wrap_global_n3)
    use grid_module, only: grid_module_global_n3 => global_n3
    implicit none
    integer(4), intent(in) :: f90wrap_global_n3

    grid_module_global_n3 = f90wrap_global_n3
end subroutine f90wrap_grid_module__set__global_n3

subroutine f90wrap_grid_module__get__global_n(f90wrap_global_n)
    use grid_module, only: grid_module_global_n => global_n
    implicit none
    integer(4), intent(out) :: f90wrap_global_n

    f90wrap_global_n = grid_module_global_n
end subroutine f90wrap_grid_module__get__global_n

subroutine f90wrap_grid_module__set__global_n(f90wrap_global_n)
    use grid_module, only: grid_module_global_n => global_n
    implicit none
    integer(4), intent(in) :: f90wrap_global_n

    grid_module_global_n = f90wrap_global_n
end subroutine f90wrap_grid_module__set__global_n

! End of module grid_module defined in file Struct_module.fpp


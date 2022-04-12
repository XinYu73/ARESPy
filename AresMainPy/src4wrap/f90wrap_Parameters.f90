! Module parameters defined in file Parameters.fpp

subroutine f90wrap_io_info_type__get__IU0(this, f90wrap_IU0)
    use parameters, only: io_info_type
    implicit none
    type io_info_type_ptr_type
        type(io_info_type), pointer :: p => NULL()
    end type io_info_type_ptr_type
    integer, intent(in)   :: this(2)
    type(io_info_type_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_IU0
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_IU0 = this_ptr%p%IU0
end subroutine f90wrap_io_info_type__get__IU0

subroutine f90wrap_io_info_type__set__IU0(this, f90wrap_IU0)
    use parameters, only: io_info_type
    implicit none
    type io_info_type_ptr_type
        type(io_info_type), pointer :: p => NULL()
    end type io_info_type_ptr_type
    integer, intent(in)   :: this(2)
    type(io_info_type_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_IU0
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%IU0 = f90wrap_IU0
end subroutine f90wrap_io_info_type__set__IU0

subroutine f90wrap_io_info_type__get__IU1(this, f90wrap_IU1)
    use parameters, only: io_info_type
    implicit none
    type io_info_type_ptr_type
        type(io_info_type), pointer :: p => NULL()
    end type io_info_type_ptr_type
    integer, intent(in)   :: this(2)
    type(io_info_type_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_IU1
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_IU1 = this_ptr%p%IU1
end subroutine f90wrap_io_info_type__get__IU1

subroutine f90wrap_io_info_type__set__IU1(this, f90wrap_IU1)
    use parameters, only: io_info_type
    implicit none
    type io_info_type_ptr_type
        type(io_info_type), pointer :: p => NULL()
    end type io_info_type_ptr_type
    integer, intent(in)   :: this(2)
    type(io_info_type_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_IU1
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%IU1 = f90wrap_IU1
end subroutine f90wrap_io_info_type__set__IU1

subroutine f90wrap_io_info_type__get__IU2(this, f90wrap_IU2)
    use parameters, only: io_info_type
    implicit none
    type io_info_type_ptr_type
        type(io_info_type), pointer :: p => NULL()
    end type io_info_type_ptr_type
    integer, intent(in)   :: this(2)
    type(io_info_type_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_IU2
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_IU2 = this_ptr%p%IU2
end subroutine f90wrap_io_info_type__get__IU2

subroutine f90wrap_io_info_type__set__IU2(this, f90wrap_IU2)
    use parameters, only: io_info_type
    implicit none
    type io_info_type_ptr_type
        type(io_info_type), pointer :: p => NULL()
    end type io_info_type_ptr_type
    integer, intent(in)   :: this(2)
    type(io_info_type_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_IU2
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%IU2 = f90wrap_IU2
end subroutine f90wrap_io_info_type__set__IU2

subroutine f90wrap_io_info_type__get__IU3(this, f90wrap_IU3)
    use parameters, only: io_info_type
    implicit none
    type io_info_type_ptr_type
        type(io_info_type), pointer :: p => NULL()
    end type io_info_type_ptr_type
    integer, intent(in)   :: this(2)
    type(io_info_type_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_IU3
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_IU3 = this_ptr%p%IU3
end subroutine f90wrap_io_info_type__get__IU3

subroutine f90wrap_io_info_type__set__IU3(this, f90wrap_IU3)
    use parameters, only: io_info_type
    implicit none
    type io_info_type_ptr_type
        type(io_info_type), pointer :: p => NULL()
    end type io_info_type_ptr_type
    integer, intent(in)   :: this(2)
    type(io_info_type_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_IU3
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%IU3 = f90wrap_IU3
end subroutine f90wrap_io_info_type__set__IU3

subroutine f90wrap_io_info_type__get__IU4(this, f90wrap_IU4)
    use parameters, only: io_info_type
    implicit none
    type io_info_type_ptr_type
        type(io_info_type), pointer :: p => NULL()
    end type io_info_type_ptr_type
    integer, intent(in)   :: this(2)
    type(io_info_type_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_IU4
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_IU4 = this_ptr%p%IU4
end subroutine f90wrap_io_info_type__get__IU4

subroutine f90wrap_io_info_type__set__IU4(this, f90wrap_IU4)
    use parameters, only: io_info_type
    implicit none
    type io_info_type_ptr_type
        type(io_info_type), pointer :: p => NULL()
    end type io_info_type_ptr_type
    integer, intent(in)   :: this(2)
    type(io_info_type_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_IU4
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%IU4 = f90wrap_IU4
end subroutine f90wrap_io_info_type__set__IU4

subroutine f90wrap_io_info_type__get__IU5(this, f90wrap_IU5)
    use parameters, only: io_info_type
    implicit none
    type io_info_type_ptr_type
        type(io_info_type), pointer :: p => NULL()
    end type io_info_type_ptr_type
    integer, intent(in)   :: this(2)
    type(io_info_type_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_IU5
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_IU5 = this_ptr%p%IU5
end subroutine f90wrap_io_info_type__get__IU5

subroutine f90wrap_io_info_type__set__IU5(this, f90wrap_IU5)
    use parameters, only: io_info_type
    implicit none
    type io_info_type_ptr_type
        type(io_info_type), pointer :: p => NULL()
    end type io_info_type_ptr_type
    integer, intent(in)   :: this(2)
    type(io_info_type_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_IU5
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%IU5 = f90wrap_IU5
end subroutine f90wrap_io_info_type__set__IU5

subroutine f90wrap_io_info_type__get__IU6(this, f90wrap_IU6)
    use parameters, only: io_info_type
    implicit none
    type io_info_type_ptr_type
        type(io_info_type), pointer :: p => NULL()
    end type io_info_type_ptr_type
    integer, intent(in)   :: this(2)
    type(io_info_type_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_IU6
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_IU6 = this_ptr%p%IU6
end subroutine f90wrap_io_info_type__get__IU6

subroutine f90wrap_io_info_type__set__IU6(this, f90wrap_IU6)
    use parameters, only: io_info_type
    implicit none
    type io_info_type_ptr_type
        type(io_info_type), pointer :: p => NULL()
    end type io_info_type_ptr_type
    integer, intent(in)   :: this(2)
    type(io_info_type_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_IU6
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%IU6 = f90wrap_IU6
end subroutine f90wrap_io_info_type__set__IU6

subroutine f90wrap_io_info_type__get__IU8(this, f90wrap_IU8)
    use parameters, only: io_info_type
    implicit none
    type io_info_type_ptr_type
        type(io_info_type), pointer :: p => NULL()
    end type io_info_type_ptr_type
    integer, intent(in)   :: this(2)
    type(io_info_type_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_IU8
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_IU8 = this_ptr%p%IU8
end subroutine f90wrap_io_info_type__get__IU8

subroutine f90wrap_io_info_type__set__IU8(this, f90wrap_IU8)
    use parameters, only: io_info_type
    implicit none
    type io_info_type_ptr_type
        type(io_info_type), pointer :: p => NULL()
    end type io_info_type_ptr_type
    integer, intent(in)   :: this(2)
    type(io_info_type_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_IU8
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%IU8 = f90wrap_IU8
end subroutine f90wrap_io_info_type__set__IU8

subroutine f90wrap_io_info_type_initialise(this)
    use parameters, only: io_info_type
    implicit none
    
    type io_info_type_ptr_type
        type(io_info_type), pointer :: p => NULL()
    end type io_info_type_ptr_type
    type(io_info_type_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_io_info_type_initialise

subroutine f90wrap_io_info_type_finalise(this)
    use parameters, only: io_info_type
    implicit none
    
    type io_info_type_ptr_type
        type(io_info_type), pointer :: p => NULL()
    end type io_info_type_ptr_type
    type(io_info_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_io_info_type_finalise

subroutine f90wrap_generate_infile(infile_name)
    use parameters, only: generate_infile
    implicit none
    
    character*(*) :: infile_name
    call generate_infile(infile_name=infile_name)
end subroutine f90wrap_generate_infile

subroutine f90wrap_parameters__get__iKEDF(f90wrap_iKEDF)
    use parameters, only: parameters_iKEDF => iKEDF
    implicit none
    integer(4), intent(out) :: f90wrap_iKEDF
    
    f90wrap_iKEDF = parameters_iKEDF
end subroutine f90wrap_parameters__get__iKEDF

subroutine f90wrap_parameters__set__iKEDF(f90wrap_iKEDF)
    use parameters, only: parameters_iKEDF => iKEDF
    implicit none
    integer(4), intent(in) :: f90wrap_iKEDF
    
    parameters_iKEDF = f90wrap_iKEDF
end subroutine f90wrap_parameters__set__iKEDF

subroutine f90wrap_parameters__get__iXCDF(f90wrap_iXCDF)
    use parameters, only: parameters_iXCDF => iXCDF
    implicit none
    integer(4), intent(out) :: f90wrap_iXCDF
    
    f90wrap_iXCDF = parameters_iXCDF
end subroutine f90wrap_parameters__get__iXCDF

subroutine f90wrap_parameters__set__iXCDF(f90wrap_iXCDF)
    use parameters, only: parameters_iXCDF => iXCDF
    implicit none
    integer(4), intent(in) :: f90wrap_iXCDF
    
    parameters_iXCDF = f90wrap_iXCDF
end subroutine f90wrap_parameters__set__iXCDF

subroutine f90wrap_parameters__get__finite_order(f90wrap_finite_order)
    use parameters, only: parameters_finite_order => finite_order
    implicit none
    integer(4), intent(out) :: f90wrap_finite_order
    
    f90wrap_finite_order = parameters_finite_order
end subroutine f90wrap_parameters__get__finite_order

subroutine f90wrap_parameters__set__finite_order(f90wrap_finite_order)
    use parameters, only: parameters_finite_order => finite_order
    implicit none
    integer(4), intent(in) :: f90wrap_finite_order
    
    parameters_finite_order = f90wrap_finite_order
end subroutine f90wrap_parameters__set__finite_order

subroutine f90wrap_parameters__get__ntype(f90wrap_ntype)
    use parameters, only: parameters_ntype => ntype
    implicit none
    integer(4), intent(out) :: f90wrap_ntype
    
    f90wrap_ntype = parameters_ntype
end subroutine f90wrap_parameters__get__ntype

subroutine f90wrap_parameters__set__ntype(f90wrap_ntype)
    use parameters, only: parameters_ntype => ntype
    implicit none
    integer(4), intent(in) :: f90wrap_ntype
    
    parameters_ntype = f90wrap_ntype
end subroutine f90wrap_parameters__set__ntype

subroutine f90wrap_parameters__get__NADDSTATES(f90wrap_NADDSTATES)
    use parameters, only: parameters_NADDSTATES => NADDSTATES
    implicit none
    integer(4), intent(out) :: f90wrap_NADDSTATES
    
    f90wrap_NADDSTATES = parameters_NADDSTATES
end subroutine f90wrap_parameters__get__NADDSTATES

subroutine f90wrap_parameters__set__NADDSTATES(f90wrap_NADDSTATES)
    use parameters, only: parameters_NADDSTATES => NADDSTATES
    implicit none
    integer(4), intent(in) :: f90wrap_NADDSTATES
    
    parameters_NADDSTATES = f90wrap_NADDSTATES
end subroutine f90wrap_parameters__set__NADDSTATES

subroutine f90wrap_parameters__get__ISTART(f90wrap_ISTART)
    use parameters, only: parameters_ISTART => ISTART
    implicit none
    integer(4), intent(out) :: f90wrap_ISTART
    
    f90wrap_ISTART = parameters_ISTART
end subroutine f90wrap_parameters__get__ISTART

subroutine f90wrap_parameters__set__ISTART(f90wrap_ISTART)
    use parameters, only: parameters_ISTART => ISTART
    implicit none
    integer(4), intent(in) :: f90wrap_ISTART
    
    parameters_ISTART = f90wrap_ISTART
end subroutine f90wrap_parameters__set__ISTART

subroutine f90wrap_parameters__get__Idiag(f90wrap_Idiag)
    use parameters, only: parameters_Idiag => Idiag
    implicit none
    integer(4), intent(out) :: f90wrap_Idiag
    
    f90wrap_Idiag = parameters_Idiag
end subroutine f90wrap_parameters__get__Idiag

subroutine f90wrap_parameters__set__Idiag(f90wrap_Idiag)
    use parameters, only: parameters_Idiag => Idiag
    implicit none
    integer(4), intent(in) :: f90wrap_Idiag
    
    parameters_Idiag = f90wrap_Idiag
end subroutine f90wrap_parameters__set__Idiag

subroutine f90wrap_parameters__get__CheM(f90wrap_CheM)
    use parameters, only: parameters_CheM => CheM
    implicit none
    integer(4), intent(out) :: f90wrap_CheM
    
    f90wrap_CheM = parameters_CheM
end subroutine f90wrap_parameters__get__CheM

subroutine f90wrap_parameters__set__CheM(f90wrap_CheM)
    use parameters, only: parameters_CheM => CheM
    implicit none
    integer(4), intent(in) :: f90wrap_CheM
    
    parameters_CheM = f90wrap_CheM
end subroutine f90wrap_parameters__set__CheM

subroutine f90wrap_parameters__get__CheM0(f90wrap_CheM0)
    use parameters, only: parameters_CheM0 => CheM0
    implicit none
    integer(4), intent(out) :: f90wrap_CheM0
    
    f90wrap_CheM0 = parameters_CheM0
end subroutine f90wrap_parameters__get__CheM0

subroutine f90wrap_parameters__set__CheM0(f90wrap_CheM0)
    use parameters, only: parameters_CheM0 => CheM0
    implicit none
    integer(4), intent(in) :: f90wrap_CheM0
    
    parameters_CheM0 = f90wrap_CheM0
end subroutine f90wrap_parameters__set__CheM0

subroutine f90wrap_parameters__get__Nstates(f90wrap_Nstates)
    use parameters, only: parameters_Nstates => Nstates
    implicit none
    integer(4), intent(out) :: f90wrap_Nstates
    
    f90wrap_Nstates = parameters_Nstates
end subroutine f90wrap_parameters__get__Nstates

subroutine f90wrap_parameters__set__Nstates(f90wrap_Nstates)
    use parameters, only: parameters_Nstates => Nstates
    implicit none
    integer(4), intent(in) :: f90wrap_Nstates
    
    parameters_Nstates = f90wrap_Nstates
end subroutine f90wrap_parameters__set__Nstates

subroutine f90wrap_parameters__get__Nstates_global(f90wrap_Nstates_global)
    use parameters, only: parameters_Nstates_global => Nstates_global
    implicit none
    integer(4), intent(out) :: f90wrap_Nstates_global
    
    f90wrap_Nstates_global = parameters_Nstates_global
end subroutine f90wrap_parameters__get__Nstates_global

subroutine f90wrap_parameters__set__Nstates_global(f90wrap_Nstates_global)
    use parameters, only: parameters_Nstates_global => Nstates_global
    implicit none
    integer(4), intent(in) :: f90wrap_Nstates_global
    
    parameters_Nstates_global = f90wrap_Nstates_global
end subroutine f90wrap_parameters__set__Nstates_global

subroutine f90wrap_parameters__get__NPRR(f90wrap_NPRR)
    use parameters, only: parameters_NPRR => NPRR
    implicit none
    integer(4), intent(out) :: f90wrap_NPRR
    
    f90wrap_NPRR = parameters_NPRR
end subroutine f90wrap_parameters__get__NPRR

subroutine f90wrap_parameters__set__NPRR(f90wrap_NPRR)
    use parameters, only: parameters_NPRR => NPRR
    implicit none
    integer(4), intent(in) :: f90wrap_NPRR
    
    parameters_NPRR = f90wrap_NPRR
end subroutine f90wrap_parameters__set__NPRR

subroutine f90wrap_parameters__get__Nssp(f90wrap_Nssp)
    use parameters, only: parameters_Nssp => Nssp
    implicit none
    integer(4), intent(out) :: f90wrap_Nssp
    
    f90wrap_Nssp = parameters_Nssp
end subroutine f90wrap_parameters__get__Nssp

subroutine f90wrap_parameters__set__Nssp(f90wrap_Nssp)
    use parameters, only: parameters_Nssp => Nssp
    implicit none
    integer(4), intent(in) :: f90wrap_Nssp
    
    parameters_Nssp = f90wrap_Nssp
end subroutine f90wrap_parameters__set__Nssp

subroutine f90wrap_parameters__get__KSPACING(f90wrap_KSPACING)
    use parameters, only: parameters_KSPACING => KSPACING
    implicit none
    real(8), intent(out) :: f90wrap_KSPACING
    
    f90wrap_KSPACING = parameters_KSPACING
end subroutine f90wrap_parameters__get__KSPACING

subroutine f90wrap_parameters__set__KSPACING(f90wrap_KSPACING)
    use parameters, only: parameters_KSPACING => KSPACING
    implicit none
    real(8), intent(in) :: f90wrap_KSPACING
    
    parameters_KSPACING = f90wrap_KSPACING
end subroutine f90wrap_parameters__set__KSPACING

subroutine f90wrap_parameters__array__KSHIFT(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters, only: parameters_kshift => kshift
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    dshape(1:1) = shape(parameters_KSHIFT)
    dloc = loc(parameters_KSHIFT)
end subroutine f90wrap_parameters__array__KSHIFT

subroutine f90wrap_parameters__get__init_gap(f90wrap_init_gap)
    use parameters, only: parameters_init_gap => init_gap
    implicit none
    real(8), intent(out) :: f90wrap_init_gap
    
    f90wrap_init_gap = parameters_init_gap
end subroutine f90wrap_parameters__get__init_gap

subroutine f90wrap_parameters__set__init_gap(f90wrap_init_gap)
    use parameters, only: parameters_init_gap => init_gap
    implicit none
    real(8), intent(in) :: f90wrap_init_gap
    
    parameters_init_gap = f90wrap_init_gap
end subroutine f90wrap_parameters__set__init_gap

subroutine f90wrap_parameters__get__Ecut(f90wrap_Ecut)
    use parameters, only: parameters_Ecut => Ecut
    implicit none
    real(8), intent(out) :: f90wrap_Ecut
    
    f90wrap_Ecut = parameters_Ecut
end subroutine f90wrap_parameters__get__Ecut

subroutine f90wrap_parameters__set__Ecut(f90wrap_Ecut)
    use parameters, only: parameters_Ecut => Ecut
    implicit none
    real(8), intent(in) :: f90wrap_Ecut
    
    parameters_Ecut = f90wrap_Ecut
end subroutine f90wrap_parameters__set__Ecut

subroutine f90wrap_parameters__get__LKP(f90wrap_LKP)
    use parameters, only: parameters_LKP => LKP
    implicit none
    logical, intent(out) :: f90wrap_LKP
    
    f90wrap_LKP = parameters_LKP
end subroutine f90wrap_parameters__get__LKP

subroutine f90wrap_parameters__set__LKP(f90wrap_LKP)
    use parameters, only: parameters_LKP => LKP
    implicit none
    logical, intent(in) :: f90wrap_LKP
    
    parameters_LKP = f90wrap_LKP
end subroutine f90wrap_parameters__set__LKP

subroutine f90wrap_parameters__get__PP_identifer(f90wrap_PP_identifer)
    use parameters, only: parameters_PP_identifer => PP_identifer
    implicit none
    integer(4), intent(out) :: f90wrap_PP_identifer
    
    f90wrap_PP_identifer = parameters_PP_identifer
end subroutine f90wrap_parameters__get__PP_identifer

subroutine f90wrap_parameters__set__PP_identifer(f90wrap_PP_identifer)
    use parameters, only: parameters_PP_identifer => PP_identifer
    implicit none
    integer(4), intent(in) :: f90wrap_PP_identifer
    
    parameters_PP_identifer = f90wrap_PP_identifer
end subroutine f90wrap_parameters__set__PP_identifer

subroutine f90wrap_parameters__get__Lpbc(f90wrap_Lpbc)
    use parameters, only: parameters_Lpbc => Lpbc
    implicit none
    logical, intent(out) :: f90wrap_Lpbc
    
    f90wrap_Lpbc = parameters_Lpbc
end subroutine f90wrap_parameters__get__Lpbc

subroutine f90wrap_parameters__set__Lpbc(f90wrap_Lpbc)
    use parameters, only: parameters_Lpbc => Lpbc
    implicit none
    logical, intent(in) :: f90wrap_Lpbc
    
    parameters_Lpbc = f90wrap_Lpbc
end subroutine f90wrap_parameters__set__Lpbc

subroutine f90wrap_parameters__get__CalForce(f90wrap_CalForce)
    use parameters, only: parameters_CalForce => CalForce
    implicit none
    logical, intent(out) :: f90wrap_CalForce
    
    f90wrap_CalForce = parameters_CalForce
end subroutine f90wrap_parameters__get__CalForce

subroutine f90wrap_parameters__set__CalForce(f90wrap_CalForce)
    use parameters, only: parameters_CalForce => CalForce
    implicit none
    logical, intent(in) :: f90wrap_CalForce
    
    parameters_CalForce = f90wrap_CalForce
end subroutine f90wrap_parameters__set__CalForce

subroutine f90wrap_parameters__get__hartree_method(f90wrap_hartree_method)
    use parameters, only: parameters_hartree_method => hartree_method
    implicit none
    integer(4), intent(out) :: f90wrap_hartree_method
    
    f90wrap_hartree_method = parameters_hartree_method
end subroutine f90wrap_parameters__get__hartree_method

subroutine f90wrap_parameters__set__hartree_method(f90wrap_hartree_method)
    use parameters, only: parameters_hartree_method => hartree_method
    implicit none
    integer(4), intent(in) :: f90wrap_hartree_method
    
    parameters_hartree_method = f90wrap_hartree_method
end subroutine f90wrap_parameters__set__hartree_method

subroutine f90wrap_parameters__get__Lcell(f90wrap_Lcell)
    use parameters, only: parameters_Lcell => Lcell
    implicit none
    real(8), intent(out) :: f90wrap_Lcell
    
    f90wrap_Lcell = parameters_Lcell
end subroutine f90wrap_parameters__get__Lcell

subroutine f90wrap_parameters__set__Lcell(f90wrap_Lcell)
    use parameters, only: parameters_Lcell => Lcell
    implicit none
    real(8), intent(in) :: f90wrap_Lcell
    
    parameters_Lcell = f90wrap_Lcell
end subroutine f90wrap_parameters__set__Lcell

subroutine f90wrap_parameters__get__Lcellbohr(f90wrap_Lcellbohr)
    use parameters, only: parameters_Lcellbohr => Lcellbohr
    implicit none
    real(8), intent(out) :: f90wrap_Lcellbohr
    
    f90wrap_Lcellbohr = parameters_Lcellbohr
end subroutine f90wrap_parameters__get__Lcellbohr

subroutine f90wrap_parameters__set__Lcellbohr(f90wrap_Lcellbohr)
    use parameters, only: parameters_Lcellbohr => Lcellbohr
    implicit none
    real(8), intent(in) :: f90wrap_Lcellbohr
    
    parameters_Lcellbohr = f90wrap_Lcellbohr
end subroutine f90wrap_parameters__set__Lcellbohr

subroutine f90wrap_parameters__get__IsoRmax(f90wrap_IsoRmax)
    use parameters, only: parameters_IsoRmax => IsoRmax
    implicit none
    real(8), intent(out) :: f90wrap_IsoRmax
    
    f90wrap_IsoRmax = parameters_IsoRmax
end subroutine f90wrap_parameters__get__IsoRmax

subroutine f90wrap_parameters__set__IsoRmax(f90wrap_IsoRmax)
    use parameters, only: parameters_IsoRmax => IsoRmax
    implicit none
    real(8), intent(in) :: f90wrap_IsoRmax
    
    parameters_IsoRmax = f90wrap_IsoRmax
end subroutine f90wrap_parameters__set__IsoRmax

subroutine f90wrap_parameters__get__IsoRmaxbohr(f90wrap_IsoRmaxbohr)
    use parameters, only: parameters_IsoRmaxbohr => IsoRmaxbohr
    implicit none
    real(8), intent(out) :: f90wrap_IsoRmaxbohr
    
    f90wrap_IsoRmaxbohr = parameters_IsoRmaxbohr
end subroutine f90wrap_parameters__get__IsoRmaxbohr

subroutine f90wrap_parameters__set__IsoRmaxbohr(f90wrap_IsoRmaxbohr)
    use parameters, only: parameters_IsoRmaxbohr => IsoRmaxbohr
    implicit none
    real(8), intent(in) :: f90wrap_IsoRmaxbohr
    
    parameters_IsoRmaxbohr = f90wrap_IsoRmaxbohr
end subroutine f90wrap_parameters__set__IsoRmaxbohr

subroutine f90wrap_parameters__get__RadiusMax(f90wrap_RadiusMax)
    use parameters, only: parameters_RadiusMax => RadiusMax
    implicit none
    real(8), intent(out) :: f90wrap_RadiusMax
    
    f90wrap_RadiusMax = parameters_RadiusMax
end subroutine f90wrap_parameters__get__RadiusMax

subroutine f90wrap_parameters__set__RadiusMax(f90wrap_RadiusMax)
    use parameters, only: parameters_RadiusMax => RadiusMax
    implicit none
    real(8), intent(in) :: f90wrap_RadiusMax
    
    parameters_RadiusMax = f90wrap_RadiusMax
end subroutine f90wrap_parameters__set__RadiusMax

subroutine f90wrap_parameters__get__NVC(f90wrap_NVC)
    use parameters, only: parameters_NVC => NVC
    implicit none
    integer(4), intent(out) :: f90wrap_NVC
    
    f90wrap_NVC = parameters_NVC
end subroutine f90wrap_parameters__get__NVC

subroutine f90wrap_parameters__set__NVC(f90wrap_NVC)
    use parameters, only: parameters_NVC => NVC
    implicit none
    integer(4), intent(in) :: f90wrap_NVC
    
    parameters_NVC = f90wrap_NVC
end subroutine f90wrap_parameters__set__NVC

subroutine f90wrap_parameters__get__ISOnorder(f90wrap_ISOnorder)
    use parameters, only: parameters_ISOnorder => ISOnorder
    implicit none
    integer(4), intent(out) :: f90wrap_ISOnorder
    
    f90wrap_ISOnorder = parameters_ISOnorder
end subroutine f90wrap_parameters__get__ISOnorder

subroutine f90wrap_parameters__set__ISOnorder(f90wrap_ISOnorder)
    use parameters, only: parameters_ISOnorder => ISOnorder
    implicit none
    integer(4), intent(in) :: f90wrap_ISOnorder
    
    parameters_ISOnorder = f90wrap_ISOnorder
end subroutine f90wrap_parameters__set__ISOnorder

subroutine f90wrap_parameters__get__TOLCG(f90wrap_TOLCG)
    use parameters, only: parameters_TOLCG => TOLCG
    implicit none
    real(8), intent(out) :: f90wrap_TOLCG
    
    f90wrap_TOLCG = parameters_TOLCG
end subroutine f90wrap_parameters__get__TOLCG

subroutine f90wrap_parameters__set__TOLCG(f90wrap_TOLCG)
    use parameters, only: parameters_TOLCG => TOLCG
    implicit none
    real(8), intent(in) :: f90wrap_TOLCG
    
    parameters_TOLCG = f90wrap_TOLCG
end subroutine f90wrap_parameters__set__TOLCG

subroutine f90wrap_parameters__get__NFCD(f90wrap_NFCD)
    use parameters, only: parameters_NFCD => NFCD
    implicit none
    integer(4), intent(out) :: f90wrap_NFCD
    
    f90wrap_NFCD = parameters_NFCD
end subroutine f90wrap_parameters__get__NFCD

subroutine f90wrap_parameters__set__NFCD(f90wrap_NFCD)
    use parameters, only: parameters_NFCD => NFCD
    implicit none
    integer(4), intent(in) :: f90wrap_NFCD
    
    parameters_NFCD = f90wrap_NFCD
end subroutine f90wrap_parameters__set__NFCD

subroutine f90wrap_parameters__get__ISOLmax(f90wrap_ISOLmax)
    use parameters, only: parameters_ISOLmax => ISOLmax
    implicit none
    integer(4), intent(out) :: f90wrap_ISOLmax
    
    f90wrap_ISOLmax = parameters_ISOLmax
end subroutine f90wrap_parameters__get__ISOLmax

subroutine f90wrap_parameters__set__ISOLmax(f90wrap_ISOLmax)
    use parameters, only: parameters_ISOLmax => ISOLmax
    implicit none
    integer(4), intent(in) :: f90wrap_ISOLmax
    
    parameters_ISOLmax = f90wrap_ISOLmax
end subroutine f90wrap_parameters__set__ISOLmax

subroutine f90wrap_parameters__get__iprec_fmm(f90wrap_iprec_fmm)
    use parameters, only: parameters_iprec_fmm => iprec_fmm
    implicit none
    integer(4), intent(out) :: f90wrap_iprec_fmm
    
    f90wrap_iprec_fmm = parameters_iprec_fmm
end subroutine f90wrap_parameters__get__iprec_fmm

subroutine f90wrap_parameters__set__iprec_fmm(f90wrap_iprec_fmm)
    use parameters, only: parameters_iprec_fmm => iprec_fmm
    implicit none
    integer(4), intent(in) :: f90wrap_iprec_fmm
    
    parameters_iprec_fmm = f90wrap_iprec_fmm
end subroutine f90wrap_parameters__set__iprec_fmm

subroutine f90wrap_parameters__get__ADDcharge(f90wrap_ADDcharge)
    use parameters, only: parameters_ADDcharge => ADDcharge
    implicit none
    integer(4), intent(out) :: f90wrap_ADDcharge
    
    f90wrap_ADDcharge = parameters_ADDcharge
end subroutine f90wrap_parameters__get__ADDcharge

subroutine f90wrap_parameters__set__ADDcharge(f90wrap_ADDcharge)
    use parameters, only: parameters_ADDcharge => ADDcharge
    implicit none
    integer(4), intent(in) :: f90wrap_ADDcharge
    
    parameters_ADDcharge = f90wrap_ADDcharge
end subroutine f90wrap_parameters__set__ADDcharge

subroutine f90wrap_parameters__get__cell_shape(f90wrap_cell_shape)
    use parameters, only: parameters_cell_shape => cell_shape
    implicit none
    integer(4), intent(out) :: f90wrap_cell_shape
    
    f90wrap_cell_shape = parameters_cell_shape
end subroutine f90wrap_parameters__get__cell_shape

subroutine f90wrap_parameters__set__cell_shape(f90wrap_cell_shape)
    use parameters, only: parameters_cell_shape => cell_shape
    implicit none
    integer(4), intent(in) :: f90wrap_cell_shape
    
    parameters_cell_shape = f90wrap_cell_shape
end subroutine f90wrap_parameters__set__cell_shape

subroutine f90wrap_parameters__get__cell_thick(f90wrap_cell_thick)
    use parameters, only: parameters_cell_thick => cell_thick
    implicit none
    real(8), intent(out) :: f90wrap_cell_thick
    
    f90wrap_cell_thick = parameters_cell_thick
end subroutine f90wrap_parameters__get__cell_thick

subroutine f90wrap_parameters__set__cell_thick(f90wrap_cell_thick)
    use parameters, only: parameters_cell_thick => cell_thick
    implicit none
    real(8), intent(in) :: f90wrap_cell_thick
    
    parameters_cell_thick = f90wrap_cell_thick
end subroutine f90wrap_parameters__set__cell_thick

subroutine f90wrap_parameters__get__Lpbc2iso(f90wrap_Lpbc2iso)
    use parameters, only: parameters_Lpbc2iso => Lpbc2iso
    implicit none
    logical, intent(out) :: f90wrap_Lpbc2iso
    
    f90wrap_Lpbc2iso = parameters_Lpbc2iso
end subroutine f90wrap_parameters__get__Lpbc2iso

subroutine f90wrap_parameters__set__Lpbc2iso(f90wrap_Lpbc2iso)
    use parameters, only: parameters_Lpbc2iso => Lpbc2iso
    implicit none
    logical, intent(in) :: f90wrap_Lpbc2iso
    
    parameters_Lpbc2iso = f90wrap_Lpbc2iso
end subroutine f90wrap_parameters__set__Lpbc2iso

subroutine f90wrap_parameters__get__Lradius_auto(f90wrap_Lradius_auto)
    use parameters, only: parameters_Lradius_auto => Lradius_auto
    implicit none
    logical, intent(out) :: f90wrap_Lradius_auto
    
    f90wrap_Lradius_auto = parameters_Lradius_auto
end subroutine f90wrap_parameters__get__Lradius_auto

subroutine f90wrap_parameters__set__Lradius_auto(f90wrap_Lradius_auto)
    use parameters, only: parameters_Lradius_auto => Lradius_auto
    implicit none
    logical, intent(in) :: f90wrap_Lradius_auto
    
    parameters_Lradius_auto = f90wrap_Lradius_auto
end subroutine f90wrap_parameters__set__Lradius_auto

subroutine f90wrap_parameters__get__BLOCK_MBNB(f90wrap_BLOCK_MBNB)
    use parameters, only: parameters_BLOCK_MBNB => BLOCK_MBNB
    implicit none
    integer(4), intent(out) :: f90wrap_BLOCK_MBNB
    
    f90wrap_BLOCK_MBNB = parameters_BLOCK_MBNB
end subroutine f90wrap_parameters__get__BLOCK_MBNB

subroutine f90wrap_parameters__set__BLOCK_MBNB(f90wrap_BLOCK_MBNB)
    use parameters, only: parameters_BLOCK_MBNB => BLOCK_MBNB
    implicit none
    integer(4), intent(in) :: f90wrap_BLOCK_MBNB
    
    parameters_BLOCK_MBNB = f90wrap_BLOCK_MBNB
end subroutine f90wrap_parameters__set__BLOCK_MBNB

subroutine f90wrap_parameters__get__debug_out(f90wrap_debug_out)
    use parameters, only: parameters_debug_out => debug_out
    implicit none
    integer(4), intent(out) :: f90wrap_debug_out
    
    f90wrap_debug_out = parameters_debug_out
end subroutine f90wrap_parameters__get__debug_out

subroutine f90wrap_parameters__set__debug_out(f90wrap_debug_out)
    use parameters, only: parameters_debug_out => debug_out
    implicit none
    integer(4), intent(in) :: f90wrap_debug_out
    
    parameters_debug_out = f90wrap_debug_out
end subroutine f90wrap_parameters__set__debug_out

subroutine f90wrap_parameters__get__Wexict(f90wrap_Wexict)
    use parameters, only: parameters_Wexict => Wexict
    implicit none
    real(8), intent(out) :: f90wrap_Wexict
    
    f90wrap_Wexict = parameters_Wexict
end subroutine f90wrap_parameters__get__Wexict

subroutine f90wrap_parameters__set__Wexict(f90wrap_Wexict)
    use parameters, only: parameters_Wexict => Wexict
    implicit none
    real(8), intent(in) :: f90wrap_Wexict
    
    parameters_Wexict = f90wrap_Wexict
end subroutine f90wrap_parameters__set__Wexict

subroutine f90wrap_parameters__get__LBvK(f90wrap_LBvK)
    use parameters, only: parameters_LBvK => LBvK
    implicit none
    logical, intent(out) :: f90wrap_LBvK
    
    f90wrap_LBvK = parameters_LBvK
end subroutine f90wrap_parameters__get__LBvK

subroutine f90wrap_parameters__set__LBvK(f90wrap_LBvK)
    use parameters, only: parameters_LBvK => LBvK
    implicit none
    logical, intent(in) :: f90wrap_LBvK
    
    parameters_LBvK = f90wrap_LBvK
end subroutine f90wrap_parameters__set__LBvK

subroutine f90wrap_parameters__get__LFIRST(f90wrap_LFIRST)
    use parameters, only: parameters_LFIRST => LFIRST
    implicit none
    logical, intent(out) :: f90wrap_LFIRST
    
    f90wrap_LFIRST = parameters_LFIRST
end subroutine f90wrap_parameters__get__LFIRST

subroutine f90wrap_parameters__set__LFIRST(f90wrap_LFIRST)
    use parameters, only: parameters_LFIRST => LFIRST
    implicit none
    logical, intent(in) :: f90wrap_LFIRST
    
    parameters_LFIRST = f90wrap_LFIRST
end subroutine f90wrap_parameters__set__LFIRST

subroutine f90wrap_parameters__get__LINRHO(f90wrap_LINRHO)
    use parameters, only: parameters_LINRHO => LINRHO
    implicit none
    logical, intent(out) :: f90wrap_LINRHO
    
    f90wrap_LINRHO = parameters_LINRHO
end subroutine f90wrap_parameters__get__LINRHO

subroutine f90wrap_parameters__set__LINRHO(f90wrap_LINRHO)
    use parameters, only: parameters_LINRHO => LINRHO
    implicit none
    logical, intent(in) :: f90wrap_LINRHO
    
    parameters_LINRHO = f90wrap_LINRHO
end subroutine f90wrap_parameters__set__LINRHO

subroutine f90wrap_parameters__get__LRadRho(f90wrap_LRadRho)
    use parameters, only: parameters_LRadRho => LRadRho
    implicit none
    logical, intent(out) :: f90wrap_LRadRho
    
    f90wrap_LRadRho = parameters_LRadRho
end subroutine f90wrap_parameters__get__LRadRho

subroutine f90wrap_parameters__set__LRadRho(f90wrap_LRadRho)
    use parameters, only: parameters_LRadRho => LRadRho
    implicit none
    logical, intent(in) :: f90wrap_LRadRho
    
    parameters_LRadRho = f90wrap_LRadRho
end subroutine f90wrap_parameters__set__LRadRho

subroutine f90wrap_parameters__get__LRROrthNorm(f90wrap_LRROrthNorm)
    use parameters, only: parameters_LRROrthNorm => LRROrthNorm
    implicit none
    logical, intent(out) :: f90wrap_LRROrthNorm
    
    f90wrap_LRROrthNorm = parameters_LRROrthNorm
end subroutine f90wrap_parameters__get__LRROrthNorm

subroutine f90wrap_parameters__set__LRROrthNorm(f90wrap_LRROrthNorm)
    use parameters, only: parameters_LRROrthNorm => LRROrthNorm
    implicit none
    logical, intent(in) :: f90wrap_LRROrthNorm
    
    parameters_LRROrthNorm = f90wrap_LRROrthNorm
end subroutine f90wrap_parameters__set__LRROrthNorm

subroutine f90wrap_parameters__get__Lrandom(f90wrap_Lrandom)
    use parameters, only: parameters_Lrandom => Lrandom
    implicit none
    logical, intent(out) :: f90wrap_Lrandom
    
    f90wrap_Lrandom = parameters_Lrandom
end subroutine f90wrap_parameters__get__Lrandom

subroutine f90wrap_parameters__set__Lrandom(f90wrap_Lrandom)
    use parameters, only: parameters_Lrandom => Lrandom
    implicit none
    logical, intent(in) :: f90wrap_Lrandom
    
    parameters_Lrandom = f90wrap_Lrandom
end subroutine f90wrap_parameters__set__Lrandom

subroutine f90wrap_parameters__get__system_name(f90wrap_system_name)
    use parameters, only: parameters_system_name => system_name
    implicit none
    character(30), intent(out) :: f90wrap_system_name
    
    f90wrap_system_name = parameters_system_name
end subroutine f90wrap_parameters__get__system_name

subroutine f90wrap_parameters__set__system_name(f90wrap_system_name)
    use parameters, only: parameters_system_name => system_name
    implicit none
    character(30), intent(in) :: f90wrap_system_name
    
    parameters_system_name = f90wrap_system_name
end subroutine f90wrap_parameters__set__system_name

subroutine f90wrap_parameters__get__cellfile_name(f90wrap_cellfile_name)
    use parameters, only: parameters_cellfile_name => cellfile_name
    implicit none
    character(30), intent(out) :: f90wrap_cellfile_name
    
    f90wrap_cellfile_name = parameters_cellfile_name
end subroutine f90wrap_parameters__get__cellfile_name

subroutine f90wrap_parameters__set__cellfile_name(f90wrap_cellfile_name)
    use parameters, only: parameters_cellfile_name => cellfile_name
    implicit none
    character(30), intent(in) :: f90wrap_cellfile_name
    
    parameters_cellfile_name = f90wrap_cellfile_name
end subroutine f90wrap_parameters__set__cellfile_name

subroutine f90wrap_parameters__get__outfile(f90wrap_outfile)
    use parameters, only: parameters_outfile => outfile
    implicit none
    character(30), intent(out) :: f90wrap_outfile
    
    f90wrap_outfile = parameters_outfile
end subroutine f90wrap_parameters__get__outfile

subroutine f90wrap_parameters__set__outfile(f90wrap_outfile)
    use parameters, only: parameters_outfile => outfile
    implicit none
    character(30), intent(in) :: f90wrap_outfile
    
    parameters_outfile = f90wrap_outfile
end subroutine f90wrap_parameters__set__outfile

subroutine f90wrap_parameters__array__ppfile_name(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters, only: parameters_ppfile_name => ppfile_name
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 2
    dshape(1:2) = (/len(parameters_ppfile_name(1)), shape(parameters_ppfile_name)/)
    dloc = loc(parameters_ppfile_name)
end subroutine f90wrap_parameters__array__ppfile_name

subroutine f90wrap_parameters__array__elements(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters, only: parameters_elements => elements
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 2
    dshape(1:2) = (/len(parameters_elements(1)), shape(parameters_elements)/)
    dloc = loc(parameters_elements)
end subroutine f90wrap_parameters__array__elements

subroutine f90wrap_parameters__get__Nspin(f90wrap_Nspin)
    use parameters, only: parameters_Nspin => Nspin
    implicit none
    integer(4), intent(out) :: f90wrap_Nspin
    
    f90wrap_Nspin = parameters_Nspin
end subroutine f90wrap_parameters__get__Nspin

subroutine f90wrap_parameters__set__Nspin(f90wrap_Nspin)
    use parameters, only: parameters_Nspin => Nspin
    implicit none
    integer(4), intent(in) :: f90wrap_Nspin
    
    parameters_Nspin = f90wrap_Nspin
end subroutine f90wrap_parameters__set__Nspin

subroutine f90wrap_parameters__get__Nsmear(f90wrap_Nsmear)
    use parameters, only: parameters_Nsmear => Nsmear
    implicit none
    integer(4), intent(out) :: f90wrap_Nsmear
    
    f90wrap_Nsmear = parameters_Nsmear
end subroutine f90wrap_parameters__get__Nsmear

subroutine f90wrap_parameters__set__Nsmear(f90wrap_Nsmear)
    use parameters, only: parameters_Nsmear => Nsmear
    implicit none
    integer(4), intent(in) :: f90wrap_Nsmear
    
    parameters_Nsmear = f90wrap_Nsmear
end subroutine f90wrap_parameters__set__Nsmear

subroutine f90wrap_parameters__get__Wsmear(f90wrap_Wsmear)
    use parameters, only: parameters_Wsmear => Wsmear
    implicit none
    real(8), intent(out) :: f90wrap_Wsmear
    
    f90wrap_Wsmear = parameters_Wsmear
end subroutine f90wrap_parameters__get__Wsmear

subroutine f90wrap_parameters__set__Wsmear(f90wrap_Wsmear)
    use parameters, only: parameters_Wsmear => Wsmear
    implicit none
    real(8), intent(in) :: f90wrap_Wsmear
    
    parameters_Wsmear = f90wrap_Wsmear
end subroutine f90wrap_parameters__set__Wsmear

subroutine f90wrap_parameters__get__IMIXER(f90wrap_IMIXER)
    use parameters, only: parameters_IMIXER => IMIXER
    implicit none
    integer(4), intent(out) :: f90wrap_IMIXER
    
    f90wrap_IMIXER = parameters_IMIXER
end subroutine f90wrap_parameters__get__IMIXER

subroutine f90wrap_parameters__set__IMIXER(f90wrap_IMIXER)
    use parameters, only: parameters_IMIXER => IMIXER
    implicit none
    integer(4), intent(in) :: f90wrap_IMIXER
    
    parameters_IMIXER = f90wrap_IMIXER
end subroutine f90wrap_parameters__set__IMIXER

subroutine f90wrap_parameters__get__NMITER(f90wrap_NMITER)
    use parameters, only: parameters_NMITER => NMITER
    implicit none
    integer(4), intent(out) :: f90wrap_NMITER
    
    f90wrap_NMITER = parameters_NMITER
end subroutine f90wrap_parameters__get__NMITER

subroutine f90wrap_parameters__set__NMITER(f90wrap_NMITER)
    use parameters, only: parameters_NMITER => NMITER
    implicit none
    integer(4), intent(in) :: f90wrap_NMITER
    
    parameters_NMITER = f90wrap_NMITER
end subroutine f90wrap_parameters__set__NMITER

subroutine f90wrap_parameters__get__NSMIX(f90wrap_NSMIX)
    use parameters, only: parameters_NSMIX => NSMIX
    implicit none
    integer(4), intent(out) :: f90wrap_NSMIX
    
    f90wrap_NSMIX = parameters_NSMIX
end subroutine f90wrap_parameters__get__NSMIX

subroutine f90wrap_parameters__set__NSMIX(f90wrap_NSMIX)
    use parameters, only: parameters_NSMIX => NSMIX
    implicit none
    integer(4), intent(in) :: f90wrap_NSMIX
    
    parameters_NSMIX = f90wrap_NSMIX
end subroutine f90wrap_parameters__set__NSMIX

subroutine f90wrap_parameters__get__NHMIX(f90wrap_NHMIX)
    use parameters, only: parameters_NHMIX => NHMIX
    implicit none
    integer(4), intent(out) :: f90wrap_NHMIX
    
    f90wrap_NHMIX = parameters_NHMIX
end subroutine f90wrap_parameters__get__NHMIX

subroutine f90wrap_parameters__set__NHMIX(f90wrap_NHMIX)
    use parameters, only: parameters_NHMIX => NHMIX
    implicit none
    integer(4), intent(in) :: f90wrap_NHMIX
    
    parameters_NHMIX = f90wrap_NHMIX
end subroutine f90wrap_parameters__set__NHMIX

subroutine f90wrap_parameters__get__NHMIN(f90wrap_NHMIN)
    use parameters, only: parameters_NHMIN => NHMIN
    implicit none
    integer(4), intent(out) :: f90wrap_NHMIN
    
    f90wrap_NHMIN = parameters_NHMIN
end subroutine f90wrap_parameters__get__NHMIN

subroutine f90wrap_parameters__set__NHMIN(f90wrap_NHMIN)
    use parameters, only: parameters_NHMIN => NHMIN
    implicit none
    integer(4), intent(in) :: f90wrap_NHMIN
    
    parameters_NHMIN = f90wrap_NHMIN
end subroutine f90wrap_parameters__set__NHMIN

subroutine f90wrap_parameters__get__MALPHA(f90wrap_MALPHA)
    use parameters, only: parameters_MALPHA => MALPHA
    implicit none
    real(8), intent(out) :: f90wrap_MALPHA
    
    f90wrap_MALPHA = parameters_MALPHA
end subroutine f90wrap_parameters__get__MALPHA

subroutine f90wrap_parameters__set__MALPHA(f90wrap_MALPHA)
    use parameters, only: parameters_MALPHA => MALPHA
    implicit none
    real(8), intent(in) :: f90wrap_MALPHA
    
    parameters_MALPHA = f90wrap_MALPHA
end subroutine f90wrap_parameters__set__MALPHA

subroutine f90wrap_parameters__get__MBETA(f90wrap_MBETA)
    use parameters, only: parameters_MBETA => MBETA
    implicit none
    real(8), intent(out) :: f90wrap_MBETA
    
    f90wrap_MBETA = parameters_MBETA
end subroutine f90wrap_parameters__get__MBETA

subroutine f90wrap_parameters__set__MBETA(f90wrap_MBETA)
    use parameters, only: parameters_MBETA => MBETA
    implicit none
    real(8), intent(in) :: f90wrap_MBETA
    
    parameters_MBETA = f90wrap_MBETA
end subroutine f90wrap_parameters__set__MBETA

subroutine f90wrap_parameters__get__AMIX(f90wrap_AMIX)
    use parameters, only: parameters_AMIX => AMIX
    implicit none
    real(8), intent(out) :: f90wrap_AMIX
    
    f90wrap_AMIX = parameters_AMIX
end subroutine f90wrap_parameters__get__AMIX

subroutine f90wrap_parameters__set__AMIX(f90wrap_AMIX)
    use parameters, only: parameters_AMIX => AMIX
    implicit none
    real(8), intent(in) :: f90wrap_AMIX
    
    parameters_AMIX = f90wrap_AMIX
end subroutine f90wrap_parameters__set__AMIX

subroutine f90wrap_parameters__get__BMIX(f90wrap_BMIX)
    use parameters, only: parameters_BMIX => BMIX
    implicit none
    real(8), intent(out) :: f90wrap_BMIX
    
    f90wrap_BMIX = parameters_BMIX
end subroutine f90wrap_parameters__get__BMIX

subroutine f90wrap_parameters__set__BMIX(f90wrap_BMIX)
    use parameters, only: parameters_BMIX => BMIX
    implicit none
    real(8), intent(in) :: f90wrap_BMIX
    
    parameters_BMIX = f90wrap_BMIX
end subroutine f90wrap_parameters__set__BMIX

subroutine f90wrap_parameters__array__RESTA(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters, only: parameters_resta => resta
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    dshape(1:1) = shape(parameters_RESTA)
    dloc = loc(parameters_RESTA)
end subroutine f90wrap_parameters__array__RESTA

subroutine f90wrap_parameters__get__W0AM(f90wrap_W0AM)
    use parameters, only: parameters_W0AM => W0AM
    implicit none
    real(8), intent(out) :: f90wrap_W0AM
    
    f90wrap_W0AM = parameters_W0AM
end subroutine f90wrap_parameters__get__W0AM

subroutine f90wrap_parameters__set__W0AM(f90wrap_W0AM)
    use parameters, only: parameters_W0AM => W0AM
    implicit none
    real(8), intent(in) :: f90wrap_W0AM
    
    parameters_W0AM = f90wrap_W0AM
end subroutine f90wrap_parameters__set__W0AM

subroutine f90wrap_parameters__get__RTOL(f90wrap_RTOL)
    use parameters, only: parameters_RTOL => RTOL
    implicit none
    real(8), intent(out) :: f90wrap_RTOL
    
    f90wrap_RTOL = parameters_RTOL
end subroutine f90wrap_parameters__get__RTOL

subroutine f90wrap_parameters__set__RTOL(f90wrap_RTOL)
    use parameters, only: parameters_RTOL => RTOL
    implicit none
    real(8), intent(in) :: f90wrap_RTOL
    
    parameters_RTOL = f90wrap_RTOL
end subroutine f90wrap_parameters__set__RTOL

subroutine f90wrap_parameters__get__ETOL(f90wrap_ETOL)
    use parameters, only: parameters_ETOL => ETOL
    implicit none
    real(8), intent(out) :: f90wrap_ETOL
    
    f90wrap_ETOL = parameters_ETOL
end subroutine f90wrap_parameters__get__ETOL

subroutine f90wrap_parameters__set__ETOL(f90wrap_ETOL)
    use parameters, only: parameters_ETOL => ETOL
    implicit none
    real(8), intent(in) :: f90wrap_ETOL
    
    parameters_ETOL = f90wrap_ETOL
end subroutine f90wrap_parameters__set__ETOL

subroutine f90wrap_parameters__get__LSUB(f90wrap_LSUB)
    use parameters, only: parameters_LSUB => LSUB
    implicit none
    logical, intent(out) :: f90wrap_LSUB
    
    f90wrap_LSUB = parameters_LSUB
end subroutine f90wrap_parameters__get__LSUB

subroutine f90wrap_parameters__set__LSUB(f90wrap_LSUB)
    use parameters, only: parameters_LSUB => LSUB
    implicit none
    logical, intent(in) :: f90wrap_LSUB
    
    parameters_LSUB = f90wrap_LSUB
end subroutine f90wrap_parameters__set__LSUB

subroutine f90wrap_parameters__get__LFAT(f90wrap_LFAT)
    use parameters, only: parameters_LFAT => LFAT
    implicit none
    logical, intent(out) :: f90wrap_LFAT
    
    f90wrap_LFAT = parameters_LFAT
end subroutine f90wrap_parameters__get__LFAT

subroutine f90wrap_parameters__set__LFAT(f90wrap_LFAT)
    use parameters, only: parameters_LFAT => LFAT
    implicit none
    logical, intent(in) :: f90wrap_LFAT
    
    parameters_LFAT = f90wrap_LFAT
end subroutine f90wrap_parameters__set__LFAT

subroutine f90wrap_parameters__get__LONE(f90wrap_LONE)
    use parameters, only: parameters_LONE => LONE
    implicit none
    logical, intent(out) :: f90wrap_LONE
    
    f90wrap_LONE = parameters_LONE
end subroutine f90wrap_parameters__get__LONE

subroutine f90wrap_parameters__set__LONE(f90wrap_LONE)
    use parameters, only: parameters_LONE => LONE
    implicit none
    logical, intent(in) :: f90wrap_LONE
    
    parameters_LONE = f90wrap_LONE
end subroutine f90wrap_parameters__set__LONE

subroutine f90wrap_parameters__get__Nsub(f90wrap_Nsub)
    use parameters, only: parameters_Nsub => Nsub
    implicit none
    integer(4), intent(out) :: f90wrap_Nsub
    
    f90wrap_Nsub = parameters_Nsub
end subroutine f90wrap_parameters__get__Nsub

subroutine f90wrap_parameters__set__Nsub(f90wrap_Nsub)
    use parameters, only: parameters_Nsub => Nsub
    implicit none
    integer(4), intent(in) :: f90wrap_Nsub
    
    parameters_Nsub = f90wrap_Nsub
end subroutine f90wrap_parameters__set__Nsub

subroutine f90wrap_parameters__get__NFAT(f90wrap_NFAT)
    use parameters, only: parameters_NFAT => NFAT
    implicit none
    integer(4), intent(out) :: f90wrap_NFAT
    
    f90wrap_NFAT = parameters_NFAT
end subroutine f90wrap_parameters__get__NFAT

subroutine f90wrap_parameters__set__NFAT(f90wrap_NFAT)
    use parameters, only: parameters_NFAT => NFAT
    implicit none
    integer(4), intent(in) :: f90wrap_NFAT
    
    parameters_NFAT = f90wrap_NFAT
end subroutine f90wrap_parameters__set__NFAT

subroutine f90wrap_parameters__array__TFVW(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters, only: parameters_tfvw => tfvw
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    dshape(1:1) = shape(parameters_TFVW)
    dloc = loc(parameters_TFVW)
end subroutine f90wrap_parameters__array__TFVW

subroutine f90wrap_parameters__get__LOFDFT(f90wrap_LOFDFT)
    use parameters, only: parameters_LOFDFT => LOFDFT
    implicit none
    logical, intent(out) :: f90wrap_LOFDFT
    
    f90wrap_LOFDFT = parameters_LOFDFT
end subroutine f90wrap_parameters__get__LOFDFT

subroutine f90wrap_parameters__set__LOFDFT(f90wrap_LOFDFT)
    use parameters, only: parameters_LOFDFT => LOFDFT
    implicit none
    logical, intent(in) :: f90wrap_LOFDFT
    
    parameters_LOFDFT = f90wrap_LOFDFT
end subroutine f90wrap_parameters__set__LOFDFT

subroutine f90wrap_parameters__get__Lband(f90wrap_Lband)
    use parameters, only: parameters_Lband => Lband
    implicit none
    logical, intent(out) :: f90wrap_Lband
    
    f90wrap_Lband = parameters_Lband
end subroutine f90wrap_parameters__get__Lband

subroutine f90wrap_parameters__set__Lband(f90wrap_Lband)
    use parameters, only: parameters_Lband => Lband
    implicit none
    logical, intent(in) :: f90wrap_Lband
    
    parameters_Lband = f90wrap_Lband
end subroutine f90wrap_parameters__set__Lband

subroutine f90wrap_parameters__get__IOPM(f90wrap_IOPM)
    use parameters, only: parameters_IOPM => IOPM
    implicit none
    integer(4), intent(out) :: f90wrap_IOPM
    
    f90wrap_IOPM = parameters_IOPM
end subroutine f90wrap_parameters__get__IOPM

subroutine f90wrap_parameters__set__IOPM(f90wrap_IOPM)
    use parameters, only: parameters_IOPM => IOPM
    implicit none
    integer(4), intent(in) :: f90wrap_IOPM
    
    parameters_IOPM = f90wrap_IOPM
end subroutine f90wrap_parameters__set__IOPM

subroutine f90wrap_parameters__get__IGOAL(f90wrap_IGOAL)
    use parameters, only: parameters_IGOAL => IGOAL
    implicit none
    integer(4), intent(out) :: f90wrap_IGOAL
    
    f90wrap_IGOAL = parameters_IGOAL
end subroutine f90wrap_parameters__get__IGOAL

subroutine f90wrap_parameters__set__IGOAL(f90wrap_IGOAL)
    use parameters, only: parameters_IGOAL => IGOAL
    implicit none
    integer(4), intent(in) :: f90wrap_IGOAL
    
    parameters_IGOAL = f90wrap_IGOAL
end subroutine f90wrap_parameters__set__IGOAL

subroutine f90wrap_parameters__get__PRESS(f90wrap_PRESS)
    use parameters, only: parameters_PRESS => PRESS
    implicit none
    real(8), intent(out) :: f90wrap_PRESS
    
    f90wrap_PRESS = parameters_PRESS
end subroutine f90wrap_parameters__get__PRESS

subroutine f90wrap_parameters__set__PRESS(f90wrap_PRESS)
    use parameters, only: parameters_PRESS => PRESS
    implicit none
    real(8), intent(in) :: f90wrap_PRESS
    
    parameters_PRESS = f90wrap_PRESS
end subroutine f90wrap_parameters__set__PRESS

subroutine f90wrap_parameters__get__TOLF(f90wrap_TOLF)
    use parameters, only: parameters_TOLF => TOLF
    implicit none
    real(8), intent(out) :: f90wrap_TOLF
    
    f90wrap_TOLF = parameters_TOLF
end subroutine f90wrap_parameters__get__TOLF

subroutine f90wrap_parameters__set__TOLF(f90wrap_TOLF)
    use parameters, only: parameters_TOLF => TOLF
    implicit none
    real(8), intent(in) :: f90wrap_TOLF
    
    parameters_TOLF = f90wrap_TOLF
end subroutine f90wrap_parameters__set__TOLF

subroutine f90wrap_parameters__get__TOLP(f90wrap_TOLP)
    use parameters, only: parameters_TOLP => TOLP
    implicit none
    real(8), intent(out) :: f90wrap_TOLP
    
    f90wrap_TOLP = parameters_TOLP
end subroutine f90wrap_parameters__get__TOLP

subroutine f90wrap_parameters__set__TOLP(f90wrap_TOLP)
    use parameters, only: parameters_TOLP => TOLP
    implicit none
    real(8), intent(in) :: f90wrap_TOLP
    
    parameters_TOLP = f90wrap_TOLP
end subroutine f90wrap_parameters__set__TOLP

subroutine f90wrap_parameters__get__TIMES(f90wrap_TIMES)
    use parameters, only: parameters_TIMES => TIMES
    implicit none
    real(8), intent(out) :: f90wrap_TIMES
    
    f90wrap_TIMES = parameters_TIMES
end subroutine f90wrap_parameters__get__TIMES

subroutine f90wrap_parameters__set__TIMES(f90wrap_TIMES)
    use parameters, only: parameters_TIMES => TIMES
    implicit none
    real(8), intent(in) :: f90wrap_TIMES
    
    parameters_TIMES = f90wrap_TIMES
end subroutine f90wrap_parameters__set__TIMES

subroutine f90wrap_parameters__get__ISP(f90wrap_ISP)
    use parameters, only: parameters_ISP => ISP
    implicit none
    integer(4), intent(out) :: f90wrap_ISP
    
    f90wrap_ISP = parameters_ISP
end subroutine f90wrap_parameters__get__ISP

subroutine f90wrap_parameters__set__ISP(f90wrap_ISP)
    use parameters, only: parameters_ISP => ISP
    implicit none
    integer(4), intent(in) :: f90wrap_ISP
    
    parameters_ISP = f90wrap_ISP
end subroutine f90wrap_parameters__set__ISP

subroutine f90wrap_parameters__get__LDG(f90wrap_LDG)
    use parameters, only: parameters_LDG => LDG
    implicit none
    logical, intent(out) :: f90wrap_LDG
    
    f90wrap_LDG = parameters_LDG
end subroutine f90wrap_parameters__get__LDG

subroutine f90wrap_parameters__set__LDG(f90wrap_LDG)
    use parameters, only: parameters_LDG => LDG
    implicit none
    logical, intent(in) :: f90wrap_LDG
    
    parameters_LDG = f90wrap_LDG
end subroutine f90wrap_parameters__set__LDG

subroutine f90wrap_parameters__get__NDG(f90wrap_NDG)
    use parameters, only: parameters_NDG => NDG
    implicit none
    integer(4), intent(out) :: f90wrap_NDG
    
    f90wrap_NDG = parameters_NDG
end subroutine f90wrap_parameters__get__NDG

subroutine f90wrap_parameters__set__NDG(f90wrap_NDG)
    use parameters, only: parameters_NDG => NDG
    implicit none
    integer(4), intent(in) :: f90wrap_NDG
    
    parameters_NDG = f90wrap_NDG
end subroutine f90wrap_parameters__set__NDG

subroutine f90wrap_parameters__get__n_near(f90wrap_n_near)
    use parameters, only: parameters_n_near => n_near
    implicit none
    integer(4), intent(out) :: f90wrap_n_near
    
    f90wrap_n_near = parameters_n_near
end subroutine f90wrap_parameters__get__n_near

subroutine f90wrap_parameters__set__n_near(f90wrap_n_near)
    use parameters, only: parameters_n_near => n_near
    implicit none
    integer(4), intent(in) :: f90wrap_n_near
    
    parameters_n_near = f90wrap_n_near
end subroutine f90wrap_parameters__set__n_near

subroutine f90wrap_parameters__get__inpol(f90wrap_inpol)
    use parameters, only: parameters_inpol => inpol
    implicit none
    integer(4), intent(out) :: f90wrap_inpol
    
    f90wrap_inpol = parameters_inpol
end subroutine f90wrap_parameters__get__inpol

subroutine f90wrap_parameters__set__inpol(f90wrap_inpol)
    use parameters, only: parameters_inpol => inpol
    implicit none
    integer(4), intent(in) :: f90wrap_inpol
    
    parameters_inpol = f90wrap_inpol
end subroutine f90wrap_parameters__set__inpol

subroutine f90wrap_parameters__get__IDinit(f90wrap_IDinit)
    use parameters, only: parameters_IDinit => IDinit
    implicit none
    integer(4), intent(out) :: f90wrap_IDinit
    
    f90wrap_IDinit = parameters_IDinit
end subroutine f90wrap_parameters__get__IDinit

subroutine f90wrap_parameters__set__IDinit(f90wrap_IDinit)
    use parameters, only: parameters_IDinit => IDinit
    implicit none
    integer(4), intent(in) :: f90wrap_IDinit
    
    parameters_IDinit = f90wrap_IDinit
end subroutine f90wrap_parameters__set__IDinit

subroutine f90wrap_parameters__get__MO_file(f90wrap_MO_file)
    use parameters, only: parameters_MO_file => MO_file
    implicit none
    character(30), intent(out) :: f90wrap_MO_file
    
    f90wrap_MO_file = parameters_MO_file
end subroutine f90wrap_parameters__get__MO_file

subroutine f90wrap_parameters__set__MO_file(f90wrap_MO_file)
    use parameters, only: parameters_MO_file => MO_file
    implicit none
    character(30), intent(in) :: f90wrap_MO_file
    
    parameters_MO_file = f90wrap_MO_file
end subroutine f90wrap_parameters__set__MO_file

subroutine f90wrap_parameters__get__POTIM(f90wrap_POTIM)
    use parameters, only: parameters_POTIM => POTIM
    implicit none
    real(8), intent(out) :: f90wrap_POTIM
    
    f90wrap_POTIM = parameters_POTIM
end subroutine f90wrap_parameters__get__POTIM

subroutine f90wrap_parameters__set__POTIM(f90wrap_POTIM)
    use parameters, only: parameters_POTIM => POTIM
    implicit none
    real(8), intent(in) :: f90wrap_POTIM
    
    parameters_POTIM = f90wrap_POTIM
end subroutine f90wrap_parameters__set__POTIM

subroutine f90wrap_parameters__get__EDIFFG(f90wrap_EDIFFG)
    use parameters, only: parameters_EDIFFG => EDIFFG
    implicit none
    real(8), intent(out) :: f90wrap_EDIFFG
    
    f90wrap_EDIFFG = parameters_EDIFFG
end subroutine f90wrap_parameters__get__EDIFFG

subroutine f90wrap_parameters__set__EDIFFG(f90wrap_EDIFFG)
    use parameters, only: parameters_EDIFFG => EDIFFG
    implicit none
    real(8), intent(in) :: f90wrap_EDIFFG
    
    parameters_EDIFFG = f90wrap_EDIFFG
end subroutine f90wrap_parameters__set__EDIFFG

subroutine f90wrap_parameters__get__step_fixrho(f90wrap_step_fixrho)
    use parameters, only: parameters_step_fixrho => step_fixrho
    implicit none
    integer(4), intent(out) :: f90wrap_step_fixrho
    
    f90wrap_step_fixrho = parameters_step_fixrho
end subroutine f90wrap_parameters__get__step_fixrho

subroutine f90wrap_parameters__set__step_fixrho(f90wrap_step_fixrho)
    use parameters, only: parameters_step_fixrho => step_fixrho
    implicit none
    integer(4), intent(in) :: f90wrap_step_fixrho
    
    parameters_step_fixrho = f90wrap_step_fixrho
end subroutine f90wrap_parameters__set__step_fixrho

subroutine f90wrap_parameters__array__fix_xyz(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters, only: parameters_fix_xyz => fix_xyz
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    dshape(1:1) = shape(parameters_fix_xyz)
    dloc = loc(parameters_fix_xyz)
end subroutine f90wrap_parameters__array__fix_xyz

subroutine f90wrap_parameters__get__pstress(f90wrap_pstress)
    use parameters, only: parameters_pstress => pstress
    implicit none
    real(8), intent(out) :: f90wrap_pstress
    
    f90wrap_pstress = parameters_pstress
end subroutine f90wrap_parameters__get__pstress

subroutine f90wrap_parameters__set__pstress(f90wrap_pstress)
    use parameters, only: parameters_pstress => pstress
    implicit none
    real(8), intent(in) :: f90wrap_pstress
    
    parameters_pstress = f90wrap_pstress
end subroutine f90wrap_parameters__set__pstress

subroutine f90wrap_parameters__get__ramax(f90wrap_ramax)
    use parameters, only: parameters_ramax => ramax
    implicit none
    real(8), intent(out) :: f90wrap_ramax
    
    f90wrap_ramax = parameters_ramax
end subroutine f90wrap_parameters__get__ramax

subroutine f90wrap_parameters__set__ramax(f90wrap_ramax)
    use parameters, only: parameters_ramax => ramax
    implicit none
    real(8), intent(in) :: f90wrap_ramax
    
    parameters_ramax = f90wrap_ramax
end subroutine f90wrap_parameters__set__ramax

subroutine f90wrap_parameters__get__Grad_order(f90wrap_Grad_order)
    use parameters, only: parameters_Grad_order => Grad_order
    implicit none
    integer(4), intent(out) :: f90wrap_Grad_order
    
    f90wrap_Grad_order = parameters_Grad_order
end subroutine f90wrap_parameters__get__Grad_order

subroutine f90wrap_parameters__set__Grad_order(f90wrap_Grad_order)
    use parameters, only: parameters_Grad_order => Grad_order
    implicit none
    integer(4), intent(in) :: f90wrap_Grad_order
    
    parameters_Grad_order = f90wrap_Grad_order
end subroutine f90wrap_parameters__set__Grad_order

subroutine f90wrap_parameters__get__Maxnpts(f90wrap_Maxnpts)
    use parameters, only: parameters_Maxnpts => Maxnpts
    implicit none
    integer(4), intent(out) :: f90wrap_Maxnpts
    
    f90wrap_Maxnpts = parameters_Maxnpts
end subroutine f90wrap_parameters__get__Maxnpts

subroutine f90wrap_parameters__set__Maxnpts(f90wrap_Maxnpts)
    use parameters, only: parameters_Maxnpts => Maxnpts
    implicit none
    integer(4), intent(in) :: f90wrap_Maxnpts
    
    parameters_Maxnpts = f90wrap_Maxnpts
end subroutine f90wrap_parameters__set__Maxnpts

subroutine f90wrap_parameters__get__nevshift(f90wrap_nevshift)
    use parameters, only: parameters_nevshift => nevshift
    implicit none
    integer(4), intent(out) :: f90wrap_nevshift
    
    f90wrap_nevshift = parameters_nevshift
end subroutine f90wrap_parameters__get__nevshift

subroutine f90wrap_parameters__set__nevshift(f90wrap_nevshift)
    use parameters, only: parameters_nevshift => nevshift
    implicit none
    integer(4), intent(in) :: f90wrap_nevshift
    
    parameters_nevshift = f90wrap_nevshift
end subroutine f90wrap_parameters__set__nevshift

subroutine f90wrap_parameters__array__gridn(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters, only: parameters_gridn => gridn
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    dshape(1:1) = shape(parameters_gridn)
    dloc = loc(parameters_gridn)
end subroutine f90wrap_parameters__array__gridn

subroutine f90wrap_parameters__get__LGamma(f90wrap_LGamma)
    use parameters, only: parameters_LGamma => LGamma
    implicit none
    logical, intent(out) :: f90wrap_LGamma
    
    f90wrap_LGamma = parameters_LGamma
end subroutine f90wrap_parameters__get__LGamma

subroutine f90wrap_parameters__set__LGamma(f90wrap_LGamma)
    use parameters, only: parameters_LGamma => LGamma
    implicit none
    logical, intent(in) :: f90wrap_LGamma
    
    parameters_LGamma = f90wrap_LGamma
end subroutine f90wrap_parameters__set__LGamma

subroutine f90wrap_parameters__get__Lcore_val(f90wrap_Lcore_val)
    use parameters, only: parameters_Lcore_val => Lcore_val
    implicit none
    logical, intent(out) :: f90wrap_Lcore_val
    
    f90wrap_Lcore_val = parameters_Lcore_val
end subroutine f90wrap_parameters__get__Lcore_val

subroutine f90wrap_parameters__set__Lcore_val(f90wrap_Lcore_val)
    use parameters, only: parameters_Lcore_val => Lcore_val
    implicit none
    logical, intent(in) :: f90wrap_Lcore_val
    
    parameters_Lcore_val = f90wrap_Lcore_val
end subroutine f90wrap_parameters__set__Lcore_val

subroutine f90wrap_parameters__array__IOUNITS(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters, only: parameters_iounits => iounits
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    dshape(1:1) = shape(parameters_IOUNITS)
    dloc = loc(parameters_IOUNITS)
end subroutine f90wrap_parameters__array__IOUNITS

subroutine f90wrap_parameters__get__Isym(f90wrap_Isym)
    use parameters, only: parameters_Isym => Isym
    implicit none
    integer(4), intent(out) :: f90wrap_Isym
    
    f90wrap_Isym = parameters_Isym
end subroutine f90wrap_parameters__get__Isym

subroutine f90wrap_parameters__set__Isym(f90wrap_Isym)
    use parameters, only: parameters_Isym => Isym
    implicit none
    integer(4), intent(in) :: f90wrap_Isym
    
    parameters_Isym = f90wrap_Isym
end subroutine f90wrap_parameters__set__Isym

subroutine f90wrap_parameters__get__Lfinite_full(f90wrap_Lfinite_full)
    use parameters, only: parameters_Lfinite_full => Lfinite_full
    implicit none
    logical, intent(out) :: f90wrap_Lfinite_full
    
    f90wrap_Lfinite_full = parameters_Lfinite_full
end subroutine f90wrap_parameters__get__Lfinite_full

subroutine f90wrap_parameters__set__Lfinite_full(f90wrap_Lfinite_full)
    use parameters, only: parameters_Lfinite_full => Lfinite_full
    implicit none
    logical, intent(in) :: f90wrap_Lfinite_full
    
    parameters_Lfinite_full = f90wrap_Lfinite_full
end subroutine f90wrap_parameters__set__Lfinite_full

subroutine f90wrap_parameters__get__lforce(f90wrap_lforce)
    use parameters, only: parameters_lforce => lforce
    implicit none
    logical, intent(out) :: f90wrap_lforce
    
    f90wrap_lforce = parameters_lforce
end subroutine f90wrap_parameters__get__lforce

subroutine f90wrap_parameters__set__lforce(f90wrap_lforce)
    use parameters, only: parameters_lforce => lforce
    implicit none
    logical, intent(in) :: f90wrap_lforce
    
    parameters_lforce = f90wrap_lforce
end subroutine f90wrap_parameters__set__lforce

subroutine f90wrap_parameters__get__lstress(f90wrap_lstress)
    use parameters, only: parameters_lstress => lstress
    implicit none
    logical, intent(out) :: f90wrap_lstress
    
    f90wrap_lstress = parameters_lstress
end subroutine f90wrap_parameters__get__lstress

subroutine f90wrap_parameters__set__lstress(f90wrap_lstress)
    use parameters, only: parameters_lstress => lstress
    implicit none
    logical, intent(in) :: f90wrap_lstress
    
    parameters_lstress = f90wrap_lstress
end subroutine f90wrap_parameters__set__lstress

subroutine f90wrap_parameters__get__temper_elec(f90wrap_temper_elec)
    use parameters, only: parameters_temper_elec => temper_elec
    implicit none
    real(8), intent(out) :: f90wrap_temper_elec
    
    f90wrap_temper_elec = parameters_temper_elec
end subroutine f90wrap_parameters__get__temper_elec

subroutine f90wrap_parameters__set__temper_elec(f90wrap_temper_elec)
    use parameters, only: parameters_temper_elec => temper_elec
    implicit none
    real(8), intent(in) :: f90wrap_temper_elec
    
    parameters_temper_elec = f90wrap_temper_elec
end subroutine f90wrap_parameters__set__temper_elec

subroutine f90wrap_parameters__get__lwstyle(f90wrap_lwstyle)
    use parameters, only: parameters_lwstyle => lwstyle
    implicit none
    character(20), intent(out) :: f90wrap_lwstyle
    
    f90wrap_lwstyle = parameters_lwstyle
end subroutine f90wrap_parameters__get__lwstyle

subroutine f90wrap_parameters__set__lwstyle(f90wrap_lwstyle)
    use parameters, only: parameters_lwstyle => lwstyle
    implicit none
    character(20), intent(in) :: f90wrap_lwstyle
    
    parameters_lwstyle = f90wrap_lwstyle
end subroutine f90wrap_parameters__set__lwstyle

subroutine f90wrap_parameters__get__LOPT(f90wrap_LOPT)
    use parameters, only: parameters_LOPT => LOPT
    implicit none
    logical, intent(out) :: f90wrap_LOPT
    
    f90wrap_LOPT = parameters_LOPT
end subroutine f90wrap_parameters__get__LOPT

subroutine f90wrap_parameters__set__LOPT(f90wrap_LOPT)
    use parameters, only: parameters_LOPT => LOPT
    implicit none
    logical, intent(in) :: f90wrap_LOPT
    
    parameters_LOPT = f90wrap_LOPT
end subroutine f90wrap_parameters__set__LOPT

subroutine f90wrap_parameters__get__IBRON(f90wrap_IBRON)
    use parameters, only: parameters_IBRON => IBRON
    implicit none
    integer(4), intent(out) :: f90wrap_IBRON
    
    f90wrap_IBRON = parameters_IBRON
end subroutine f90wrap_parameters__get__IBRON

subroutine f90wrap_parameters__set__IBRON(f90wrap_IBRON)
    use parameters, only: parameters_IBRON => IBRON
    implicit none
    integer(4), intent(in) :: f90wrap_IBRON
    
    parameters_IBRON = f90wrap_IBRON
end subroutine f90wrap_parameters__set__IBRON

subroutine f90wrap_parameters__get__maxsave(f90wrap_maxsave)
    use parameters, only: parameters_maxsave => maxsave
    implicit none
    integer(4), intent(out) :: f90wrap_maxsave
    
    f90wrap_maxsave = parameters_maxsave
end subroutine f90wrap_parameters__get__maxsave

subroutine f90wrap_parameters__set__maxsave(f90wrap_maxsave)
    use parameters, only: parameters_maxsave => maxsave
    implicit none
    integer(4), intent(in) :: f90wrap_maxsave
    
    parameters_maxsave = f90wrap_maxsave
end subroutine f90wrap_parameters__set__maxsave

subroutine f90wrap_parameters__get__IGamma(f90wrap_IGamma)
    use parameters, only: parameters_IGamma => IGamma
    implicit none
    integer(4), intent(out) :: f90wrap_IGamma
    
    f90wrap_IGamma = parameters_IGamma
end subroutine f90wrap_parameters__get__IGamma

subroutine f90wrap_parameters__set__IGamma(f90wrap_IGamma)
    use parameters, only: parameters_IGamma => IGamma
    implicit none
    integer(4), intent(in) :: f90wrap_IGamma
    
    parameters_IGamma = f90wrap_IGamma
end subroutine f90wrap_parameters__set__IGamma

subroutine f90wrap_parameters__array__kgrid(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters, only: parameters_kgrid => kgrid
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    dshape(1:1) = shape(parameters_kgrid)
    dloc = loc(parameters_kgrid)
end subroutine f90wrap_parameters__array__kgrid

subroutine f90wrap_parameters__get__LMD(f90wrap_LMD)
    use parameters, only: parameters_LMD => LMD
    implicit none
    logical, intent(out) :: f90wrap_LMD
    
    f90wrap_LMD = parameters_LMD
end subroutine f90wrap_parameters__get__LMD

subroutine f90wrap_parameters__set__LMD(f90wrap_LMD)
    use parameters, only: parameters_LMD => LMD
    implicit none
    logical, intent(in) :: f90wrap_LMD
    
    parameters_LMD = f90wrap_LMD
end subroutine f90wrap_parameters__set__LMD

subroutine f90wrap_parameters__get__nhis(f90wrap_nhis)
    use parameters, only: parameters_nhis => nhis
    implicit none
    integer(4), intent(out) :: f90wrap_nhis
    
    f90wrap_nhis = parameters_nhis
end subroutine f90wrap_parameters__get__nhis

subroutine f90wrap_parameters__set__nhis(f90wrap_nhis)
    use parameters, only: parameters_nhis => nhis
    implicit none
    integer(4), intent(in) :: f90wrap_nhis
    
    parameters_nhis = f90wrap_nhis
end subroutine f90wrap_parameters__set__nhis

subroutine f90wrap_parameters__get__rdfr(f90wrap_rdfr)
    use parameters, only: parameters_rdfr => rdfr
    implicit none
    real(8), intent(out) :: f90wrap_rdfr
    
    f90wrap_rdfr = parameters_rdfr
end subroutine f90wrap_parameters__get__rdfr

subroutine f90wrap_parameters__set__rdfr(f90wrap_rdfr)
    use parameters, only: parameters_rdfr => rdfr
    implicit none
    real(8), intent(in) :: f90wrap_rdfr
    
    parameters_rdfr = f90wrap_rdfr
end subroutine f90wrap_parameters__set__rdfr

subroutine f90wrap_parameters__get__thermostat(f90wrap_thermostat)
    use parameters, only: parameters_thermostat => thermostat
    implicit none
    character(30), intent(out) :: f90wrap_thermostat
    
    f90wrap_thermostat = parameters_thermostat
end subroutine f90wrap_parameters__get__thermostat

subroutine f90wrap_parameters__set__thermostat(f90wrap_thermostat)
    use parameters, only: parameters_thermostat => thermostat
    implicit none
    character(30), intent(in) :: f90wrap_thermostat
    
    parameters_thermostat = f90wrap_thermostat
end subroutine f90wrap_parameters__set__thermostat

subroutine f90wrap_parameters__get__integrator(f90wrap_integrator)
    use parameters, only: parameters_integrator => integrator
    implicit none
    character(30), intent(out) :: f90wrap_integrator
    
    f90wrap_integrator = parameters_integrator
end subroutine f90wrap_parameters__get__integrator

subroutine f90wrap_parameters__set__integrator(f90wrap_integrator)
    use parameters, only: parameters_integrator => integrator
    implicit none
    character(30), intent(in) :: f90wrap_integrator
    
    parameters_integrator = f90wrap_integrator
end subroutine f90wrap_parameters__set__integrator

subroutine f90wrap_parameters__get__sfreq(f90wrap_sfreq)
    use parameters, only: parameters_sfreq => sfreq
    implicit none
    integer(4), intent(out) :: f90wrap_sfreq
    
    f90wrap_sfreq = parameters_sfreq
end subroutine f90wrap_parameters__get__sfreq

subroutine f90wrap_parameters__set__sfreq(f90wrap_sfreq)
    use parameters, only: parameters_sfreq => sfreq
    implicit none
    integer(4), intent(in) :: f90wrap_sfreq
    
    parameters_sfreq = f90wrap_sfreq
end subroutine f90wrap_parameters__set__sfreq

subroutine f90wrap_parameters__get__cfreq(f90wrap_cfreq)
    use parameters, only: parameters_cfreq => cfreq
    implicit none
    real(8), intent(out) :: f90wrap_cfreq
    
    f90wrap_cfreq = parameters_cfreq
end subroutine f90wrap_parameters__get__cfreq

subroutine f90wrap_parameters__set__cfreq(f90wrap_cfreq)
    use parameters, only: parameters_cfreq => cfreq
    implicit none
    real(8), intent(in) :: f90wrap_cfreq
    
    parameters_cfreq = f90wrap_cfreq
end subroutine f90wrap_parameters__set__cfreq

subroutine f90wrap_parameters__get__temperature(f90wrap_temperature)
    use parameters, only: parameters_temperature => temperature
    implicit none
    real(8), intent(out) :: f90wrap_temperature
    
    f90wrap_temperature = parameters_temperature
end subroutine f90wrap_parameters__get__temperature

subroutine f90wrap_parameters__set__temperature(f90wrap_temperature)
    use parameters, only: parameters_temperature => temperature
    implicit none
    real(8), intent(in) :: f90wrap_temperature
    
    parameters_temperature = f90wrap_temperature
end subroutine f90wrap_parameters__set__temperature

subroutine f90wrap_parameters__get__mste(f90wrap_mste)
    use parameters, only: parameters_mste => mste
    implicit none
    real(8), intent(out) :: f90wrap_mste
    
    f90wrap_mste = parameters_mste
end subroutine f90wrap_parameters__get__mste

subroutine f90wrap_parameters__set__mste(f90wrap_mste)
    use parameters, only: parameters_mste => mste
    implicit none
    real(8), intent(in) :: f90wrap_mste
    
    parameters_mste = f90wrap_mste
end subroutine f90wrap_parameters__set__mste

subroutine f90wrap_parameters__get__delt(f90wrap_delt)
    use parameters, only: parameters_delt => delt
    implicit none
    real(8), intent(out) :: f90wrap_delt
    
    f90wrap_delt = parameters_delt
end subroutine f90wrap_parameters__get__delt

subroutine f90wrap_parameters__set__delt(f90wrap_delt)
    use parameters, only: parameters_delt => delt
    implicit none
    real(8), intent(in) :: f90wrap_delt
    
    parameters_delt = f90wrap_delt
end subroutine f90wrap_parameters__set__delt

subroutine f90wrap_parameters__get__relaxt(f90wrap_relaxt)
    use parameters, only: parameters_relaxt => relaxt
    implicit none
    real(8), intent(out) :: f90wrap_relaxt
    
    f90wrap_relaxt = parameters_relaxt
end subroutine f90wrap_parameters__get__relaxt

subroutine f90wrap_parameters__set__relaxt(f90wrap_relaxt)
    use parameters, only: parameters_relaxt => relaxt
    implicit none
    real(8), intent(in) :: f90wrap_relaxt
    
    parameters_relaxt = f90wrap_relaxt
end subroutine f90wrap_parameters__set__relaxt

subroutine f90wrap_parameters__get__ensemble(f90wrap_ensemble)
    use parameters, only: parameters_ensemble => ensemble
    implicit none
    character(10), intent(out) :: f90wrap_ensemble
    
    f90wrap_ensemble = parameters_ensemble
end subroutine f90wrap_parameters__get__ensemble

subroutine f90wrap_parameters__set__ensemble(f90wrap_ensemble)
    use parameters, only: parameters_ensemble => ensemble
    implicit none
    character(10), intent(in) :: f90wrap_ensemble
    
    parameters_ensemble = f90wrap_ensemble
end subroutine f90wrap_parameters__set__ensemble

subroutine f90wrap_parameters__get__dof(f90wrap_dof)
    use parameters, only: parameters_dof => dof
    implicit none
    integer(4), intent(out) :: f90wrap_dof
    
    f90wrap_dof = parameters_dof
end subroutine f90wrap_parameters__get__dof

subroutine f90wrap_parameters__set__dof(f90wrap_dof)
    use parameters, only: parameters_dof => dof
    implicit none
    integer(4), intent(in) :: f90wrap_dof
    
    parameters_dof = f90wrap_dof
end subroutine f90wrap_parameters__set__dof

subroutine f90wrap_parameters__get__pext(f90wrap_pext)
    use parameters, only: parameters_pext => pext
    implicit none
    real(8), intent(out) :: f90wrap_pext
    
    f90wrap_pext = parameters_pext
end subroutine f90wrap_parameters__get__pext

subroutine f90wrap_parameters__set__pext(f90wrap_pext)
    use parameters, only: parameters_pext => pext
    implicit none
    real(8), intent(in) :: f90wrap_pext
    
    parameters_pext = f90wrap_pext
end subroutine f90wrap_parameters__set__pext

subroutine f90wrap_parameters__get__iresmd(f90wrap_iresmd)
    use parameters, only: parameters_iresmd => iresmd
    implicit none
    integer(4), intent(out) :: f90wrap_iresmd
    
    f90wrap_iresmd = parameters_iresmd
end subroutine f90wrap_parameters__get__iresmd

subroutine f90wrap_parameters__set__iresmd(f90wrap_iresmd)
    use parameters, only: parameters_iresmd => iresmd
    implicit none
    integer(4), intent(in) :: f90wrap_iresmd
    
    parameters_iresmd = f90wrap_iresmd
end subroutine f90wrap_parameters__set__iresmd

subroutine f90wrap_parameters__get__nresn(f90wrap_nresn)
    use parameters, only: parameters_nresn => nresn
    implicit none
    integer(4), intent(out) :: f90wrap_nresn
    
    f90wrap_nresn = parameters_nresn
end subroutine f90wrap_parameters__get__nresn

subroutine f90wrap_parameters__set__nresn(f90wrap_nresn)
    use parameters, only: parameters_nresn => nresn
    implicit none
    integer(4), intent(in) :: f90wrap_nresn
    
    parameters_nresn = f90wrap_nresn
end subroutine f90wrap_parameters__set__nresn

subroutine f90wrap_parameters__get__nyosh(f90wrap_nyosh)
    use parameters, only: parameters_nyosh => nyosh
    implicit none
    integer(4), intent(out) :: f90wrap_nyosh
    
    f90wrap_nyosh = parameters_nyosh
end subroutine f90wrap_parameters__get__nyosh

subroutine f90wrap_parameters__set__nyosh(f90wrap_nyosh)
    use parameters, only: parameters_nyosh => nyosh
    implicit none
    integer(4), intent(in) :: f90wrap_nyosh
    
    parameters_nyosh = f90wrap_nyosh
end subroutine f90wrap_parameters__set__nyosh

subroutine f90wrap_parameters__get__nnhc(f90wrap_nnhc)
    use parameters, only: parameters_nnhc => nnhc
    implicit none
    integer(4), intent(out) :: f90wrap_nnhc
    
    f90wrap_nnhc = parameters_nnhc
end subroutine f90wrap_parameters__get__nnhc

subroutine f90wrap_parameters__set__nnhc(f90wrap_nnhc)
    use parameters, only: parameters_nnhc => nnhc
    implicit none
    integer(4), intent(in) :: f90wrap_nnhc
    
    parameters_nnhc = f90wrap_nnhc
end subroutine f90wrap_parameters__set__nnhc

subroutine f90wrap_parameters__array__wdti2(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters, only: parameters_wdti2 => wdti2
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(parameters_wdti2)) then
        dshape(1:1) = shape(parameters_wdti2)
        dloc = loc(parameters_wdti2)
    else
        dloc = 0
    end if
end subroutine f90wrap_parameters__array__wdti2

subroutine f90wrap_parameters__array__wdti4(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters, only: parameters_wdti4 => wdti4
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(parameters_wdti4)) then
        dshape(1:1) = shape(parameters_wdti4)
        dloc = loc(parameters_wdti4)
    else
        dloc = 0
    end if
end subroutine f90wrap_parameters__array__wdti4

subroutine f90wrap_parameters__array__wdti8(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters, only: parameters_wdti8 => wdti8
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(parameters_wdti8)) then
        dshape(1:1) = shape(parameters_wdti8)
        dloc = loc(parameters_wdti8)
    else
        dloc = 0
    end if
end subroutine f90wrap_parameters__array__wdti8

subroutine f90wrap_parameters__array__qmass(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters, only: parameters_qmass => qmass
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(parameters_qmass)) then
        dshape(1:1) = shape(parameters_qmass)
        dloc = loc(parameters_qmass)
    else
        dloc = 0
    end if
end subroutine f90wrap_parameters__array__qmass

subroutine f90wrap_parameters__array__syin_coeff(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters, only: parameters_syin_coeff => syin_coeff
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(parameters_syin_coeff)) then
        dshape(1:1) = shape(parameters_syin_coeff)
        dloc = loc(parameters_syin_coeff)
    else
        dloc = 0
    end if
end subroutine f90wrap_parameters__array__syin_coeff

subroutine f90wrap_parameters__get__bmass(f90wrap_bmass)
    use parameters, only: parameters_bmass => bmass
    implicit none
    real(8), intent(out) :: f90wrap_bmass
    
    f90wrap_bmass = parameters_bmass
end subroutine f90wrap_parameters__get__bmass

subroutine f90wrap_parameters__set__bmass(f90wrap_bmass)
    use parameters, only: parameters_bmass => bmass
    implicit none
    real(8), intent(in) :: f90wrap_bmass
    
    parameters_bmass = f90wrap_bmass
end subroutine f90wrap_parameters__set__bmass

subroutine f90wrap_parameters__get__fthermo(f90wrap_fthermo)
    use parameters, only: parameters_fthermo => fthermo
    implicit none
    real(8), intent(out) :: f90wrap_fthermo
    
    f90wrap_fthermo = parameters_fthermo
end subroutine f90wrap_parameters__get__fthermo

subroutine f90wrap_parameters__set__fthermo(f90wrap_fthermo)
    use parameters, only: parameters_fthermo => fthermo
    implicit none
    real(8), intent(in) :: f90wrap_fthermo
    
    parameters_fthermo = f90wrap_fthermo
end subroutine f90wrap_parameters__set__fthermo

subroutine f90wrap_parameters__get__fbaro(f90wrap_fbaro)
    use parameters, only: parameters_fbaro => fbaro
    implicit none
    real(8), intent(out) :: f90wrap_fbaro
    
    f90wrap_fbaro = parameters_fbaro
end subroutine f90wrap_parameters__get__fbaro

subroutine f90wrap_parameters__set__fbaro(f90wrap_fbaro)
    use parameters, only: parameters_fbaro => fbaro
    implicit none
    real(8), intent(in) :: f90wrap_fbaro
    
    parameters_fbaro = f90wrap_fbaro
end subroutine f90wrap_parameters__set__fbaro

subroutine f90wrap_parameters__get__fthermown(f90wrap_fthermown)
    use parameters, only: parameters_fthermown => fthermown
    implicit none
    real(8), intent(out) :: f90wrap_fthermown
    
    f90wrap_fthermown = parameters_fthermown
end subroutine f90wrap_parameters__get__fthermown

subroutine f90wrap_parameters__set__fthermown(f90wrap_fthermown)
    use parameters, only: parameters_fthermown => fthermown
    implicit none
    real(8), intent(in) :: f90wrap_fthermown
    
    parameters_fthermown = f90wrap_fthermown
end subroutine f90wrap_parameters__set__fthermown

subroutine f90wrap_parameters__get__fbarown(f90wrap_fbarown)
    use parameters, only: parameters_fbarown => fbarown
    implicit none
    real(8), intent(out) :: f90wrap_fbarown
    
    f90wrap_fbarown = parameters_fbarown
end subroutine f90wrap_parameters__get__fbarown

subroutine f90wrap_parameters__set__fbarown(f90wrap_fbarown)
    use parameters, only: parameters_fbarown => fbarown
    implicit none
    real(8), intent(in) :: f90wrap_fbarown
    
    parameters_fbarown = f90wrap_fbarown
end subroutine f90wrap_parameters__set__fbarown

subroutine f90wrap_parameters__get__lbin(f90wrap_lbin)
    use parameters, only: parameters_lbin => lbin
    implicit none
    logical, intent(out) :: f90wrap_lbin
    
    f90wrap_lbin = parameters_lbin
end subroutine f90wrap_parameters__get__lbin

subroutine f90wrap_parameters__set__lbin(f90wrap_lbin)
    use parameters, only: parameters_lbin => lbin
    implicit none
    logical, intent(in) :: f90wrap_lbin
    
    parameters_lbin = f90wrap_lbin
end subroutine f90wrap_parameters__set__lbin

subroutine f90wrap_parameters__array__p_flag(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters, only: parameters_p_flag => p_flag
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    dshape(1:1) = shape(parameters_p_flag)
    dloc = loc(parameters_p_flag)
end subroutine f90wrap_parameters__array__p_flag

subroutine f90wrap_parameters__array__p_start(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters, only: parameters_p_start => p_start
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    dshape(1:1) = shape(parameters_p_start)
    dloc = loc(parameters_p_start)
end subroutine f90wrap_parameters__array__p_start

subroutine f90wrap_parameters__array__p_stop(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters, only: parameters_p_stop => p_stop
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    dshape(1:1) = shape(parameters_p_stop)
    dloc = loc(parameters_p_stop)
end subroutine f90wrap_parameters__array__p_stop

subroutine f90wrap_parameters__array__p_freq(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters, only: parameters_p_freq => p_freq
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    dshape(1:1) = shape(parameters_p_freq)
    dloc = loc(parameters_p_freq)
end subroutine f90wrap_parameters__array__p_freq

subroutine f90wrap_parameters__array__p_target(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters, only: parameters_p_target => p_target
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    dshape(1:1) = shape(parameters_p_target)
    dloc = loc(parameters_p_target)
end subroutine f90wrap_parameters__array__p_target

subroutine f90wrap_parameters__get__t_start(f90wrap_t_start)
    use parameters, only: parameters_t_start => t_start
    implicit none
    real(8), intent(out) :: f90wrap_t_start
    
    f90wrap_t_start = parameters_t_start
end subroutine f90wrap_parameters__get__t_start

subroutine f90wrap_parameters__set__t_start(f90wrap_t_start)
    use parameters, only: parameters_t_start => t_start
    implicit none
    real(8), intent(in) :: f90wrap_t_start
    
    parameters_t_start = f90wrap_t_start
end subroutine f90wrap_parameters__set__t_start

subroutine f90wrap_parameters__get__t_stop(f90wrap_t_stop)
    use parameters, only: parameters_t_stop => t_stop
    implicit none
    real(8), intent(out) :: f90wrap_t_stop
    
    f90wrap_t_stop = parameters_t_stop
end subroutine f90wrap_parameters__get__t_stop

subroutine f90wrap_parameters__set__t_stop(f90wrap_t_stop)
    use parameters, only: parameters_t_stop => t_stop
    implicit none
    real(8), intent(in) :: f90wrap_t_stop
    
    parameters_t_stop = f90wrap_t_stop
end subroutine f90wrap_parameters__set__t_stop

subroutine f90wrap_parameters__get__t_freq(f90wrap_t_freq)
    use parameters, only: parameters_t_freq => t_freq
    implicit none
    real(8), intent(out) :: f90wrap_t_freq
    
    f90wrap_t_freq = parameters_t_freq
end subroutine f90wrap_parameters__get__t_freq

subroutine f90wrap_parameters__set__t_freq(f90wrap_t_freq)
    use parameters, only: parameters_t_freq => t_freq
    implicit none
    real(8), intent(in) :: f90wrap_t_freq
    
    parameters_t_freq = f90wrap_t_freq
end subroutine f90wrap_parameters__set__t_freq

subroutine f90wrap_parameters__get__t_target(f90wrap_t_target)
    use parameters, only: parameters_t_target => t_target
    implicit none
    real(8), intent(out) :: f90wrap_t_target
    
    f90wrap_t_target = parameters_t_target
end subroutine f90wrap_parameters__get__t_target

subroutine f90wrap_parameters__set__t_target(f90wrap_t_target)
    use parameters, only: parameters_t_target => t_target
    implicit none
    real(8), intent(in) :: f90wrap_t_target
    
    parameters_t_target = f90wrap_t_target
end subroutine f90wrap_parameters__set__t_target

subroutine f90wrap_parameters__get__press_control(f90wrap_press_control)
    use parameters, only: parameters_press_control => press_control
    implicit none
    character(1024), intent(out) :: f90wrap_press_control
    
    f90wrap_press_control = parameters_press_control
end subroutine f90wrap_parameters__get__press_control

subroutine f90wrap_parameters__set__press_control(f90wrap_press_control)
    use parameters, only: parameters_press_control => press_control
    implicit none
    character(1024), intent(in) :: f90wrap_press_control
    
    parameters_press_control = f90wrap_press_control
end subroutine f90wrap_parameters__set__press_control

subroutine f90wrap_parameters__get__pcontrol(f90wrap_pcontrol)
    use parameters, only: parameters_pcontrol => pcontrol
    implicit none
    character(1024), intent(out) :: f90wrap_pcontrol
    
    f90wrap_pcontrol = parameters_pcontrol
end subroutine f90wrap_parameters__get__pcontrol

subroutine f90wrap_parameters__set__pcontrol(f90wrap_pcontrol)
    use parameters, only: parameters_pcontrol => pcontrol
    implicit none
    character(1024), intent(in) :: f90wrap_pcontrol
    
    parameters_pcontrol = f90wrap_pcontrol
end subroutine f90wrap_parameters__set__pcontrol

subroutine f90wrap_parameters__get__pmass(f90wrap_pmass)
    use parameters, only: parameters_pmass => pmass
    implicit none
    real(8), intent(out) :: f90wrap_pmass
    
    f90wrap_pmass = parameters_pmass
end subroutine f90wrap_parameters__get__pmass

subroutine f90wrap_parameters__set__pmass(f90wrap_pmass)
    use parameters, only: parameters_pmass => pmass
    implicit none
    real(8), intent(in) :: f90wrap_pmass
    
    parameters_pmass = f90wrap_pmass
end subroutine f90wrap_parameters__set__pmass

subroutine f90wrap_parameters__get__pdrag(f90wrap_pdrag)
    use parameters, only: parameters_pdrag => pdrag
    implicit none
    real(8), intent(out) :: f90wrap_pdrag
    
    f90wrap_pdrag = parameters_pdrag
end subroutine f90wrap_parameters__get__pdrag

subroutine f90wrap_parameters__set__pdrag(f90wrap_pdrag)
    use parameters, only: parameters_pdrag => pdrag
    implicit none
    real(8), intent(in) :: f90wrap_pdrag
    
    parameters_pdrag = f90wrap_pdrag
end subroutine f90wrap_parameters__set__pdrag

subroutine f90wrap_parameters__get__tdrag(f90wrap_tdrag)
    use parameters, only: parameters_tdrag => tdrag
    implicit none
    real(8), intent(out) :: f90wrap_tdrag
    
    f90wrap_tdrag = parameters_tdrag
end subroutine f90wrap_parameters__get__tdrag

subroutine f90wrap_parameters__set__tdrag(f90wrap_tdrag)
    use parameters, only: parameters_tdrag => tdrag
    implicit none
    real(8), intent(in) :: f90wrap_tdrag
    
    parameters_tdrag = f90wrap_tdrag
end subroutine f90wrap_parameters__set__tdrag

subroutine f90wrap_parameters__get__erate(f90wrap_erate)
    use parameters, only: parameters_erate => erate
    implicit none
    real(8), intent(out) :: f90wrap_erate
    
    f90wrap_erate = parameters_erate
end subroutine f90wrap_parameters__get__erate

subroutine f90wrap_parameters__set__erate(f90wrap_erate)
    use parameters, only: parameters_erate => erate
    implicit none
    real(8), intent(in) :: f90wrap_erate
    
    parameters_erate = f90wrap_erate
end subroutine f90wrap_parameters__set__erate

subroutine f90wrap_parameters__get__mstrain(f90wrap_mstrain)
    use parameters, only: parameters_mstrain => mstrain
    implicit none
    real(8), intent(out) :: f90wrap_mstrain
    
    f90wrap_mstrain = parameters_mstrain
end subroutine f90wrap_parameters__get__mstrain

subroutine f90wrap_parameters__set__mstrain(f90wrap_mstrain)
    use parameters, only: parameters_mstrain => mstrain
    implicit none
    real(8), intent(in) :: f90wrap_mstrain
    
    parameters_mstrain = f90wrap_mstrain
end subroutine f90wrap_parameters__set__mstrain

subroutine f90wrap_parameters__get__sdir(f90wrap_sdir)
    use parameters, only: parameters_sdir => sdir
    implicit none
    character(1024), intent(out) :: f90wrap_sdir
    
    f90wrap_sdir = parameters_sdir
end subroutine f90wrap_parameters__get__sdir

subroutine f90wrap_parameters__set__sdir(f90wrap_sdir)
    use parameters, only: parameters_sdir => sdir
    implicit none
    character(1024), intent(in) :: f90wrap_sdir
    
    parameters_sdir = f90wrap_sdir
end subroutine f90wrap_parameters__set__sdir

subroutine f90wrap_parameters__get__Lke_pot(f90wrap_Lke_pot)
    use parameters, only: parameters_Lke_pot => Lke_pot
    implicit none
    logical, intent(out) :: f90wrap_Lke_pot
    
    f90wrap_Lke_pot = parameters_Lke_pot
end subroutine f90wrap_parameters__get__Lke_pot

subroutine f90wrap_parameters__set__Lke_pot(f90wrap_Lke_pot)
    use parameters, only: parameters_Lke_pot => Lke_pot
    implicit none
    logical, intent(in) :: f90wrap_Lke_pot
    
    parameters_Lke_pot = f90wrap_Lke_pot
end subroutine f90wrap_parameters__set__Lke_pot

! End of module parameters defined in file Parameters.fpp


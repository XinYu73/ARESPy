! Module parameters defined in file Parameters.fpp

subroutine f90wrap_parameters__get__ixc(f90wrap_ixc)
    use parameters, only: parameters_ixc => ixc
    implicit none
    integer(4), intent(out) :: f90wrap_ixc
    
    f90wrap_ixc = parameters_ixc
end subroutine f90wrap_parameters__get__ixc

subroutine f90wrap_parameters__set__ixc(f90wrap_ixc)
    use parameters, only: parameters_ixc => ixc
    implicit none
    integer(4), intent(in) :: f90wrap_ixc
    
    parameters_ixc = f90wrap_ixc
end subroutine f90wrap_parameters__set__ixc

subroutine f90wrap_parameters__get__nspin(f90wrap_nspin)
    use parameters, only: parameters_nspin => nspin
    implicit none
    integer(4), intent(out) :: f90wrap_nspin
    
    f90wrap_nspin = parameters_nspin
end subroutine f90wrap_parameters__get__nspin

subroutine f90wrap_parameters__set__nspin(f90wrap_nspin)
    use parameters, only: parameters_nspin => nspin
    implicit none
    integer(4), intent(in) :: f90wrap_nspin
    
    parameters_nspin = f90wrap_nspin
end subroutine f90wrap_parameters__set__nspin

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

subroutine f90wrap_parameters__get__Naddstates(f90wrap_Naddstates)
    use parameters, only: parameters_Naddstates => Naddstates
    implicit none
    integer(4), intent(out) :: f90wrap_Naddstates
    
    f90wrap_Naddstates = parameters_Naddstates
end subroutine f90wrap_parameters__get__Naddstates

subroutine f90wrap_parameters__set__Naddstates(f90wrap_Naddstates)
    use parameters, only: parameters_Naddstates => Naddstates
    implicit none
    integer(4), intent(in) :: f90wrap_Naddstates
    
    parameters_Naddstates = f90wrap_Naddstates
end subroutine f90wrap_parameters__set__Naddstates

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

subroutine f90wrap_parameters__get__dcharge(f90wrap_dcharge)
    use parameters, only: parameters_dcharge => dcharge
    implicit none
    real(8), intent(out) :: f90wrap_dcharge
    
    f90wrap_dcharge = parameters_dcharge
end subroutine f90wrap_parameters__get__dcharge

subroutine f90wrap_parameters__set__dcharge(f90wrap_dcharge)
    use parameters, only: parameters_dcharge => dcharge
    implicit none
    real(8), intent(in) :: f90wrap_dcharge
    
    parameters_dcharge = f90wrap_dcharge
end subroutine f90wrap_parameters__set__dcharge

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

subroutine f90wrap_parameters__get__kspacing(f90wrap_kspacing)
    use parameters, only: parameters_kspacing => kspacing
    implicit none
    real(8), intent(out) :: f90wrap_kspacing
    
    f90wrap_kspacing = parameters_kspacing
end subroutine f90wrap_parameters__get__kspacing

subroutine f90wrap_parameters__set__kspacing(f90wrap_kspacing)
    use parameters, only: parameters_kspacing => kspacing
    implicit none
    real(8), intent(in) :: f90wrap_kspacing
    
    parameters_kspacing = f90wrap_kspacing
end subroutine f90wrap_parameters__set__kspacing

subroutine f90wrap_parameters__get__AtomRc(f90wrap_AtomRc)
    use parameters, only: parameters_AtomRc => AtomRc
    implicit none
    real(8), intent(out) :: f90wrap_AtomRc
    
    f90wrap_AtomRc = parameters_AtomRc
end subroutine f90wrap_parameters__get__AtomRc

subroutine f90wrap_parameters__set__AtomRc(f90wrap_AtomRc)
    use parameters, only: parameters_AtomRc => AtomRc
    implicit none
    real(8), intent(in) :: f90wrap_AtomRc
    
    parameters_AtomRc = f90wrap_AtomRc
end subroutine f90wrap_parameters__set__AtomRc

subroutine f90wrap_parameters__get__Snlcc(f90wrap_Snlcc)
    use parameters, only: parameters_Snlcc => Snlcc
    implicit none
    real(8), intent(out) :: f90wrap_Snlcc
    
    f90wrap_Snlcc = parameters_Snlcc
end subroutine f90wrap_parameters__get__Snlcc

subroutine f90wrap_parameters__set__Snlcc(f90wrap_Snlcc)
    use parameters, only: parameters_Snlcc => Snlcc
    implicit none
    real(8), intent(in) :: f90wrap_Snlcc
    
    parameters_Snlcc = f90wrap_Snlcc
end subroutine f90wrap_parameters__set__Snlcc

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

subroutine f90wrap_parameters__get__LAtomRho(f90wrap_LAtomRho)
    use parameters, only: parameters_LAtomRho => LAtomRho
    implicit none
    logical, intent(out) :: f90wrap_LAtomRho
    
    f90wrap_LAtomRho = parameters_LAtomRho
end subroutine f90wrap_parameters__get__LAtomRho

subroutine f90wrap_parameters__set__LAtomRho(f90wrap_LAtomRho)
    use parameters, only: parameters_LAtomRho => LAtomRho
    implicit none
    logical, intent(in) :: f90wrap_LAtomRho
    
    parameters_LAtomRho = f90wrap_LAtomRho
end subroutine f90wrap_parameters__set__LAtomRho

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

subroutine f90wrap_parameters__get__NHMAX(f90wrap_NHMAX)
    use parameters, only: parameters_NHMAX => NHMAX
    implicit none
    integer(4), intent(out) :: f90wrap_NHMAX
    
    f90wrap_NHMAX = parameters_NHMAX
end subroutine f90wrap_parameters__get__NHMAX

subroutine f90wrap_parameters__set__NHMAX(f90wrap_NHMAX)
    use parameters, only: parameters_NHMAX => NHMAX
    implicit none
    integer(4), intent(in) :: f90wrap_NHMAX
    
    parameters_NHMAX = f90wrap_NHMAX
end subroutine f90wrap_parameters__set__NHMAX

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

subroutine f90wrap_parameters__get__LWAVE(f90wrap_LWAVE)
    use parameters, only: parameters_LWAVE => LWAVE
    implicit none
    logical, intent(out) :: f90wrap_LWAVE
    
    f90wrap_LWAVE = parameters_LWAVE
end subroutine f90wrap_parameters__get__LWAVE

subroutine f90wrap_parameters__set__LWAVE(f90wrap_LWAVE)
    use parameters, only: parameters_LWAVE => LWAVE
    implicit none
    logical, intent(in) :: f90wrap_LWAVE
    
    parameters_LWAVE = f90wrap_LWAVE
end subroutine f90wrap_parameters__set__LWAVE

subroutine f90wrap_parameters__get__LCHARGE(f90wrap_LCHARGE)
    use parameters, only: parameters_LCHARGE => LCHARGE
    implicit none
    logical, intent(out) :: f90wrap_LCHARGE
    
    f90wrap_LCHARGE = parameters_LCHARGE
end subroutine f90wrap_parameters__get__LCHARGE

subroutine f90wrap_parameters__set__LCHARGE(f90wrap_LCHARGE)
    use parameters, only: parameters_LCHARGE => LCHARGE
    implicit none
    logical, intent(in) :: f90wrap_LCHARGE
    
    parameters_LCHARGE = f90wrap_LCHARGE
end subroutine f90wrap_parameters__set__LCHARGE

subroutine f90wrap_parameters__get__LMOM(f90wrap_LMOM)
    use parameters, only: parameters_LMOM => LMOM
    implicit none
    logical, intent(out) :: f90wrap_LMOM
    
    f90wrap_LMOM = parameters_LMOM
end subroutine f90wrap_parameters__get__LMOM

subroutine f90wrap_parameters__set__LMOM(f90wrap_LMOM)
    use parameters, only: parameters_LMOM => LMOM
    implicit none
    logical, intent(in) :: f90wrap_LMOM
    
    parameters_LMOM = f90wrap_LMOM
end subroutine f90wrap_parameters__set__LMOM

subroutine f90wrap_parameters__get__MOMsigma(f90wrap_MOMsigma)
    use parameters, only: parameters_MOMsigma => MOMsigma
    implicit none
    real(8), intent(out) :: f90wrap_MOMsigma
    
    f90wrap_MOMsigma = parameters_MOMsigma
end subroutine f90wrap_parameters__get__MOMsigma

subroutine f90wrap_parameters__set__MOMsigma(f90wrap_MOMsigma)
    use parameters, only: parameters_MOMsigma => MOMsigma
    implicit none
    real(8), intent(in) :: f90wrap_MOMsigma
    
    parameters_MOMsigma = f90wrap_MOMsigma
end subroutine f90wrap_parameters__set__MOMsigma

subroutine f90wrap_parameters__get__nwf0(f90wrap_nwf0)
    use parameters, only: parameters_nwf0 => nwf0
    implicit none
    integer(4), intent(out) :: f90wrap_nwf0
    
    f90wrap_nwf0 = parameters_nwf0
end subroutine f90wrap_parameters__get__nwf0

subroutine f90wrap_parameters__set__nwf0(f90wrap_nwf0)
    use parameters, only: parameters_nwf0 => nwf0
    implicit none
    integer(4), intent(in) :: f90wrap_nwf0
    
    parameters_nwf0 = f90wrap_nwf0
end subroutine f90wrap_parameters__set__nwf0

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

subroutine f90wrap_parameters__array__NPARA(dummy_this, nd, dtype, dshape, dloc)
    use constants
    use parameters, only: parameters_npara => npara
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    dshape(1:1) = shape(parameters_NPARA)
    dloc = loc(parameters_NPARA)
end subroutine f90wrap_parameters__array__NPARA

! End of module parameters defined in file Parameters.fpp


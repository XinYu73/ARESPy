"""
Module parameters


Defined at Parameters.fpp lines 5-83

"""
from __future__ import print_function, absolute_import, division
import _arespy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def get_ixc():
    """
    Element ixc ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 11
    
    """
    return _arespy_pkg.f90wrap_parameters__get__ixc()

def set_ixc(ixc):
    _arespy_pkg.f90wrap_parameters__set__ixc(ixc)

def get_nspin():
    """
    Element nspin ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 11
    
    """
    return _arespy_pkg.f90wrap_parameters__get__nspin()

def set_nspin(nspin):
    _arespy_pkg.f90wrap_parameters__set__nspin(nspin)

def get_finite_order():
    """
    Element finite_order ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 12
    
    """
    return _arespy_pkg.f90wrap_parameters__get__finite_order()

def set_finite_order(finite_order):
    _arespy_pkg.f90wrap_parameters__set__finite_order(finite_order)

def get_ntype():
    """
    Element ntype ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 13
    
    """
    return _arespy_pkg.f90wrap_parameters__get__ntype()

def set_ntype(ntype):
    _arespy_pkg.f90wrap_parameters__set__ntype(ntype)

def get_naddstates():
    """
    Element naddstates ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 24
    
    """
    return _arespy_pkg.f90wrap_parameters__get__naddstates()

def set_naddstates(naddstates):
    _arespy_pkg.f90wrap_parameters__set__naddstates(naddstates)

def get_igamma():
    """
    Element igamma ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 24
    
    """
    return _arespy_pkg.f90wrap_parameters__get__igamma()

def set_igamma(igamma):
    _arespy_pkg.f90wrap_parameters__set__igamma(igamma)

def get_istart():
    """
    Element istart ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 24
    
    """
    return _arespy_pkg.f90wrap_parameters__get__istart()

def set_istart(istart):
    _arespy_pkg.f90wrap_parameters__set__istart(istart)

def get_isym():
    """
    Element isym ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 24
    
    """
    return _arespy_pkg.f90wrap_parameters__get__isym()

def set_isym(isym):
    _arespy_pkg.f90wrap_parameters__set__isym(isym)

def get_idiag():
    """
    Element idiag ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 24
    
    """
    return _arespy_pkg.f90wrap_parameters__get__idiag()

def set_idiag(idiag):
    _arespy_pkg.f90wrap_parameters__set__idiag(idiag)

def get_chem():
    """
    Element chem ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 24
    
    """
    return _arespy_pkg.f90wrap_parameters__get__chem()

def set_chem(chem):
    _arespy_pkg.f90wrap_parameters__set__chem(chem)

def get_chem0():
    """
    Element chem0 ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 24
    
    """
    return _arespy_pkg.f90wrap_parameters__get__chem0()

def set_chem0(chem0):
    _arespy_pkg.f90wrap_parameters__set__chem0(chem0)

def get_nstates():
    """
    Element nstates ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 24
    
    """
    return _arespy_pkg.f90wrap_parameters__get__nstates()

def set_nstates(nstates):
    _arespy_pkg.f90wrap_parameters__set__nstates(nstates)

def get_array_gridn():
    """
    Element gridn ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 24
    
    """
    global gridn
    array_ndim, array_type, array_shape, array_handle = \
        _arespy_pkg.f90wrap_parameters__array__gridn(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        gridn = _arrays[array_handle]
    else:
        gridn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_parameters__array__gridn)
        _arrays[array_handle] = gridn
    return gridn

def set_array_gridn(gridn):
    gridn[...] = gridn

def get_array_kgrid():
    """
    Element kgrid ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 24
    
    """
    global kgrid
    array_ndim, array_type, array_shape, array_handle = \
        _arespy_pkg.f90wrap_parameters__array__kgrid(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        kgrid = _arrays[array_handle]
    else:
        kgrid = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_parameters__array__kgrid)
        _arrays[array_handle] = kgrid
    return kgrid

def set_array_kgrid(kgrid):
    kgrid[...] = kgrid

def get_nssp():
    """
    Element nssp ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 24
    
    """
    return _arespy_pkg.f90wrap_parameters__get__nssp()

def set_nssp(nssp):
    _arespy_pkg.f90wrap_parameters__set__nssp(nssp)

def get_dcharge():
    """
    Element dcharge ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 26
    
    """
    return _arespy_pkg.f90wrap_parameters__get__dcharge()

def set_dcharge(dcharge):
    _arespy_pkg.f90wrap_parameters__set__dcharge(dcharge)

def get_init_gap():
    """
    Element init_gap ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 31
    
    """
    return _arespy_pkg.f90wrap_parameters__get__init_gap()

def set_init_gap(init_gap):
    _arespy_pkg.f90wrap_parameters__set__init_gap(init_gap)

def get_ecut():
    """
    Element ecut ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 31
    
    """
    return _arespy_pkg.f90wrap_parameters__get__ecut()

def set_ecut(ecut):
    _arespy_pkg.f90wrap_parameters__set__ecut(ecut)

def get_kspacing():
    """
    Element kspacing ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 31
    
    """
    return _arespy_pkg.f90wrap_parameters__get__kspacing()

def set_kspacing(kspacing):
    _arespy_pkg.f90wrap_parameters__set__kspacing(kspacing)

def get_atomrc():
    """
    Element atomrc ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 31
    
    """
    return _arespy_pkg.f90wrap_parameters__get__atomrc()

def set_atomrc(atomrc):
    _arespy_pkg.f90wrap_parameters__set__atomrc(atomrc)

def get_snlcc():
    """
    Element snlcc ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 31
    
    """
    return _arespy_pkg.f90wrap_parameters__get__snlcc()

def set_snlcc(snlcc):
    _arespy_pkg.f90wrap_parameters__set__snlcc(snlcc)

def get_lfirst():
    """
    Element lfirst ftype=logical pytype=bool
    
    
    Defined at Parameters.fpp line 39
    
    """
    return _arespy_pkg.f90wrap_parameters__get__lfirst()

def set_lfirst(lfirst):
    _arespy_pkg.f90wrap_parameters__set__lfirst(lfirst)

def get_linrho():
    """
    Element linrho ftype=logical pytype=bool
    
    
    Defined at Parameters.fpp line 39
    
    """
    return _arespy_pkg.f90wrap_parameters__get__linrho()

def set_linrho(linrho):
    _arespy_pkg.f90wrap_parameters__set__linrho(linrho)

def get_latomrho():
    """
    Element latomrho ftype=logical pytype=bool
    
    
    Defined at Parameters.fpp line 39
    
    """
    return _arespy_pkg.f90wrap_parameters__get__latomrho()

def set_latomrho(latomrho):
    _arespy_pkg.f90wrap_parameters__set__latomrho(latomrho)

def get_lrrorthnorm():
    """
    Element lrrorthnorm ftype=logical pytype=bool
    
    
    Defined at Parameters.fpp line 39
    
    """
    return _arespy_pkg.f90wrap_parameters__get__lrrorthnorm()

def set_lrrorthnorm(lrrorthnorm):
    _arespy_pkg.f90wrap_parameters__set__lrrorthnorm(lrrorthnorm)

def get_lrandom():
    """
    Element lrandom ftype=logical pytype=bool
    
    
    Defined at Parameters.fpp line 39
    
    """
    return _arespy_pkg.f90wrap_parameters__get__lrandom()

def set_lrandom(lrandom):
    _arespy_pkg.f90wrap_parameters__set__lrandom(lrandom)

def get_lcore_val():
    """
    Element lcore_val ftype=logical pytype=bool
    
    
    Defined at Parameters.fpp line 39
    
    """
    return _arespy_pkg.f90wrap_parameters__get__lcore_val()

def set_lcore_val(lcore_val):
    _arespy_pkg.f90wrap_parameters__set__lcore_val(lcore_val)

def get_system_name():
    """
    Element system_name ftype=character(30) pytype=str
    
    
    Defined at Parameters.fpp line 41
    
    """
    return _arespy_pkg.f90wrap_parameters__get__system_name()

def set_system_name(system_name):
    _arespy_pkg.f90wrap_parameters__set__system_name(system_name)

def get_cellfile_name():
    """
    Element cellfile_name ftype=character(30) pytype=str
    
    
    Defined at Parameters.fpp line 42
    
    """
    return _arespy_pkg.f90wrap_parameters__get__cellfile_name()

def set_cellfile_name(cellfile_name):
    _arespy_pkg.f90wrap_parameters__set__cellfile_name(cellfile_name)

def get_array_ppfile_name():
    """
    Element ppfile_name ftype=character(30) pytype=str
    
    
    Defined at Parameters.fpp line 43
    
    """
    global ppfile_name
    array_ndim, array_type, array_shape, array_handle = \
        _arespy_pkg.f90wrap_parameters__array__ppfile_name(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        ppfile_name = _arrays[array_handle]
    else:
        ppfile_name = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_parameters__array__ppfile_name)
        _arrays[array_handle] = ppfile_name
    return ppfile_name

def set_array_ppfile_name(ppfile_name):
    ppfile_name[...] = ppfile_name

def get_array_elements():
    """
    Element elements ftype=character(30) pytype=str
    
    
    Defined at Parameters.fpp line 44
    
    """
    global elements
    array_ndim, array_type, array_shape, array_handle = \
        _arespy_pkg.f90wrap_parameters__array__elements(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        elements = _arrays[array_handle]
    else:
        elements = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_parameters__array__elements)
        _arrays[array_handle] = elements
    return elements

def set_array_elements(elements):
    elements[...] = elements

def get_nsmear():
    """
    Element nsmear ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 46
    
    """
    return _arespy_pkg.f90wrap_parameters__get__nsmear()

def set_nsmear(nsmear):
    _arespy_pkg.f90wrap_parameters__set__nsmear(nsmear)

def get_wsmear():
    """
    Element wsmear ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 47
    
    """
    return _arespy_pkg.f90wrap_parameters__get__wsmear()

def set_wsmear(wsmear):
    _arespy_pkg.f90wrap_parameters__set__wsmear(wsmear)

def get_imixer():
    """
    Element imixer ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 50
    
    """
    return _arespy_pkg.f90wrap_parameters__get__imixer()

def set_imixer(imixer):
    _arespy_pkg.f90wrap_parameters__set__imixer(imixer)

def get_nmiter():
    """
    Element nmiter ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 50
    
    """
    return _arespy_pkg.f90wrap_parameters__get__nmiter()

def set_nmiter(nmiter):
    _arespy_pkg.f90wrap_parameters__set__nmiter(nmiter)

def get_nsmix():
    """
    Element nsmix ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 53
    
    """
    return _arespy_pkg.f90wrap_parameters__get__nsmix()

def set_nsmix(nsmix):
    _arespy_pkg.f90wrap_parameters__set__nsmix(nsmix)

def get_nhmax():
    """
    Element nhmax ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 53
    
    """
    return _arespy_pkg.f90wrap_parameters__get__nhmax()

def set_nhmax(nhmax):
    _arespy_pkg.f90wrap_parameters__set__nhmax(nhmax)

def get_nhmin():
    """
    Element nhmin ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 53
    
    """
    return _arespy_pkg.f90wrap_parameters__get__nhmin()

def set_nhmin(nhmin):
    _arespy_pkg.f90wrap_parameters__set__nhmin(nhmin)

def get_malpha():
    """
    Element malpha ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 57
    
    """
    return _arespy_pkg.f90wrap_parameters__get__malpha()

def set_malpha(malpha):
    _arespy_pkg.f90wrap_parameters__set__malpha(malpha)

def get_mbeta():
    """
    Element mbeta ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 57
    
    """
    return _arespy_pkg.f90wrap_parameters__get__mbeta()

def set_mbeta(mbeta):
    _arespy_pkg.f90wrap_parameters__set__mbeta(mbeta)

def get_amix():
    """
    Element amix ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 57
    
    """
    return _arespy_pkg.f90wrap_parameters__get__amix()

def set_amix(amix):
    _arespy_pkg.f90wrap_parameters__set__amix(amix)

def get_bmix():
    """
    Element bmix ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 57
    
    """
    return _arespy_pkg.f90wrap_parameters__get__bmix()

def set_bmix(bmix):
    _arespy_pkg.f90wrap_parameters__set__bmix(bmix)

def get_array_resta():
    """
    Element resta ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 58
    
    """
    global resta
    array_ndim, array_type, array_shape, array_handle = \
        _arespy_pkg.f90wrap_parameters__array__resta(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        resta = _arrays[array_handle]
    else:
        resta = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_parameters__array__resta)
        _arrays[array_handle] = resta
    return resta

def set_array_resta(resta):
    resta[...] = resta

def get_w0am():
    """
    Element w0am ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 60
    
    """
    return _arespy_pkg.f90wrap_parameters__get__w0am()

def set_w0am(w0am):
    _arespy_pkg.f90wrap_parameters__set__w0am(w0am)

def get_rtol():
    """
    Element rtol ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 62
    
    """
    return _arespy_pkg.f90wrap_parameters__get__rtol()

def set_rtol(rtol):
    _arespy_pkg.f90wrap_parameters__set__rtol(rtol)

def get_etol():
    """
    Element etol ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 62
    
    """
    return _arespy_pkg.f90wrap_parameters__get__etol()

def set_etol(etol):
    _arespy_pkg.f90wrap_parameters__set__etol(etol)

def get_lband():
    """
    Element lband ftype=logical pytype=bool
    
    
    Defined at Parameters.fpp line 64
    
    """
    return _arespy_pkg.f90wrap_parameters__get__lband()

def set_lband(lband):
    _arespy_pkg.f90wrap_parameters__set__lband(lband)

def get_iopm():
    """
    Element iopm ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 68
    
    """
    return _arespy_pkg.f90wrap_parameters__get__iopm()

def set_iopm(iopm):
    _arespy_pkg.f90wrap_parameters__set__iopm(iopm)

def get_igoal():
    """
    Element igoal ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 68
    
    """
    return _arespy_pkg.f90wrap_parameters__get__igoal()

def set_igoal(igoal):
    _arespy_pkg.f90wrap_parameters__set__igoal(igoal)

def get_press():
    """
    Element press ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 72
    
    """
    return _arespy_pkg.f90wrap_parameters__get__press()

def set_press(press):
    _arespy_pkg.f90wrap_parameters__set__press(press)

def get_tolf():
    """
    Element tolf ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 72
    
    """
    return _arespy_pkg.f90wrap_parameters__get__tolf()

def set_tolf(tolf):
    _arespy_pkg.f90wrap_parameters__set__tolf(tolf)

def get_tolp():
    """
    Element tolp ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 72
    
    """
    return _arespy_pkg.f90wrap_parameters__get__tolp()

def set_tolp(tolp):
    _arespy_pkg.f90wrap_parameters__set__tolp(tolp)

def get_times():
    """
    Element times ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 72
    
    """
    return _arespy_pkg.f90wrap_parameters__get__times()

def set_times(times):
    _arespy_pkg.f90wrap_parameters__set__times(times)

def get_lwave():
    """
    Element lwave ftype=logical pytype=bool
    
    
    Defined at Parameters.fpp line 75
    
    """
    return _arespy_pkg.f90wrap_parameters__get__lwave()

def set_lwave(lwave):
    _arespy_pkg.f90wrap_parameters__set__lwave(lwave)

def get_lcharge():
    """
    Element lcharge ftype=logical pytype=bool
    
    
    Defined at Parameters.fpp line 75
    
    """
    return _arespy_pkg.f90wrap_parameters__get__lcharge()

def set_lcharge(lcharge):
    _arespy_pkg.f90wrap_parameters__set__lcharge(lcharge)

def get_lmom():
    """
    Element lmom ftype=logical pytype=bool
    
    
    Defined at Parameters.fpp line 77
    
    """
    return _arespy_pkg.f90wrap_parameters__get__lmom()

def set_lmom(lmom):
    _arespy_pkg.f90wrap_parameters__set__lmom(lmom)

def get_momsigma():
    """
    Element momsigma ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 78
    
    """
    return _arespy_pkg.f90wrap_parameters__get__momsigma()

def set_momsigma(momsigma):
    _arespy_pkg.f90wrap_parameters__set__momsigma(momsigma)

def get_nwf0():
    """
    Element nwf0 ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 79
    
    """
    return _arespy_pkg.f90wrap_parameters__get__nwf0()

def set_nwf0(nwf0):
    _arespy_pkg.f90wrap_parameters__set__nwf0(nwf0)

def get_block_mbnb():
    """
    Element block_mbnb ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 82
    
    """
    return _arespy_pkg.f90wrap_parameters__get__block_mbnb()

def set_block_mbnb(block_mbnb):
    _arespy_pkg.f90wrap_parameters__set__block_mbnb(block_mbnb)

def get_nstates_global():
    """
    Element nstates_global ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 82
    
    """
    return _arespy_pkg.f90wrap_parameters__get__nstates_global()

def set_nstates_global(nstates_global):
    _arespy_pkg.f90wrap_parameters__set__nstates_global(nstates_global)

def get_array_npara():
    """
    Element npara ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 83
    
    """
    global npara
    array_ndim, array_type, array_shape, array_handle = \
        _arespy_pkg.f90wrap_parameters__array__npara(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        npara = _arrays[array_handle]
    else:
        npara = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_parameters__array__npara)
        _arrays[array_handle] = npara
    return npara

def set_array_npara(npara):
    npara[...] = npara


_array_initialisers = [get_array_gridn, get_array_kgrid, get_array_ppfile_name, \
    get_array_elements, get_array_resta, get_array_npara]
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "parameters".')

for func in _dt_array_initialisers:
    func()

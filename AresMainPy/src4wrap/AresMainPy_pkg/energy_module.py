"""
Module energy_module


Defined at Energy_module.fpp lines 5-763

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def totalenergy(psi, rhos, eval, fmenergy, ets, llast):
    """
    totalenergy(psi, rhos, eval, fmenergy, ets, llast)
    
    
    Defined at Energy_module.fpp lines 32-127
    
    Parameters
    ----------
    psi : complex array
    rhos : float array
    eval : float array
    fmenergy : float
    ets : float
    llast : bool
    
    """
    _AresMainPy_pkg.f90wrap_totalenergy(psi=psi, rhos=rhos, eval=eval, \
        fmenergy=fmenergy, ets=ets, llast=llast)

def out_ke_potential(efm, vie, vh, vxc):
    """
    out_ke_potential(efm, vie, vh, vxc)
    
    
    Defined at Energy_module.fpp lines 129-156
    
    Parameters
    ----------
    efm : float
    vie : float array
    vh : float array
    vxc : float array
    
    """
    _AresMainPy_pkg.f90wrap_out_ke_potential(efm=efm, vie=vie, vh=vh, vxc=vxc)

def totalenergy_gamma(psi, rhos, eval, fmenergy, ets, llast):
    """
    totalenergy_gamma(psi, rhos, eval, fmenergy, ets, llast)
    
    
    Defined at Energy_module.fpp lines 159-250
    
    Parameters
    ----------
    psi : float array
    rhos : float array
    eval : float array
    fmenergy : float
    ets : float
    llast : bool
    
    """
    _AresMainPy_pkg.f90wrap_totalenergy_gamma(psi=psi, rhos=rhos, eval=eval, \
        fmenergy=fmenergy, ets=ets, llast=llast)

def ehartree(rho, vhart):
    """
    ehart = ehartree(rho, vhart)
    
    
    Defined at Energy_module.fpp lines 253-270
    
    Parameters
    ----------
    rho : float array
    vhart : float array
    
    Returns
    -------
    ehart : float
    
    """
    ehart = _AresMainPy_pkg.f90wrap_ehartree(rho=rho, vhart=vhart)
    return ehart

def ehartree_iso(rho, vhart):
    """
    ehart = ehartree_iso(rho, vhart)
    
    
    Defined at Energy_module.fpp lines 273-284
    
    Parameters
    ----------
    rho : float array
    vhart : float array
    
    Returns
    -------
    ehart : float
    
    """
    ehart = _AresMainPy_pkg.f90wrap_ehartree_iso(rho=rho, vhart=vhart)
    return ehart

def e_ie(rho, vion):
    """
    eie = e_ie(rho, vion)
    
    
    Defined at Energy_module.fpp lines 287-304
    
    Parameters
    ----------
    rho : float array
    vion : float array
    
    Returns
    -------
    eie : float
    
    """
    eie = _AresMainPy_pkg.f90wrap_e_ie(rho=rho, vion=vion)
    return eie

def ebands(eval, wke):
    """
    eband = ebands(eval, wke)
    
    
    Defined at Energy_module.fpp lines 307-338
    
    Parameters
    ----------
    eval : float array
    wke : float array
    
    Returns
    -------
    eband : float
    
    """
    eband = _AresMainPy_pkg.f90wrap_ebands(eval=eval, wke=wke)
    return eband

def bvk_totalenergy(rhos, eval, fmenergy, ets, llast):
    """
    bvk_totalenergy(rhos, eval, fmenergy, ets, llast)
    
    
    Defined at Energy_module.fpp lines 344-394
    
    Parameters
    ----------
    rhos : float array
    eval : float array
    fmenergy : float
    ets : float
    llast : bool
    
    """
    _AresMainPy_pkg.f90wrap_bvk_totalenergy(rhos=rhos, eval=eval, fmenergy=fmenergy, \
        ets=ets, llast=llast)

def bvk_ebands(eval, focc):
    """
    eband = bvk_ebands(eval, focc)
    
    
    Defined at Energy_module.fpp lines 397-420
    
    Parameters
    ----------
    eval : float array
    focc : float array
    
    Returns
    -------
    eband : float
    
    """
    eband = _AresMainPy_pkg.f90wrap_bvk_ebands(eval=eval, focc=focc)
    return eband

def prr_ebands(veff, phi, pbar):
    """
    eband = prr_ebands(veff, phi, pbar)
    
    
    Defined at Energy_module.fpp lines 426-456
    
    Parameters
    ----------
    veff : float array
    phi : float array
    pbar : float array
    
    Returns
    -------
    eband : float
    
    """
    eband = _AresMainPy_pkg.f90wrap_prr_ebands(veff=veff, phi=phi, pbar=pbar)
    return eband

def prr_totalenergy(rhos, phi, pbar, fmenergy, ets, llast):
    """
    prr_totalenergy(rhos, phi, pbar, fmenergy, ets, llast)
    
    
    Defined at Energy_module.fpp lines 459-516
    
    Parameters
    ----------
    rhos : float array
    phi : float array
    pbar : float array
    fmenergy : float
    ets : float
    llast : bool
    
    """
    _AresMainPy_pkg.f90wrap_prr_totalenergy(rhos=rhos, phi=phi, pbar=pbar, \
        fmenergy=fmenergy, ets=ets, llast=llast)

def iso_totalenergy(rhos, eval, fmenergy, ets, llast):
    """
    iso_totalenergy(rhos, eval, fmenergy, ets, llast)
    
    
    Defined at Energy_module.fpp lines 522-646
    
    Parameters
    ----------
    rhos : float array
    eval : float array
    fmenergy : float
    ets : float
    llast : bool
    
    ------------------total energy calculation------------------
    """
    _AresMainPy_pkg.f90wrap_iso_totalenergy(rhos=rhos, eval=eval, fmenergy=fmenergy, \
        ets=ets, llast=llast)

def iso_ebands(eval, focc):
    """
    eband = iso_ebands(eval, focc)
    
    
    Defined at Energy_module.fpp lines 650-673
    
    Parameters
    ----------
    eval : float array
    focc : float array
    
    Returns
    -------
    eband : float
    
    """
    eband = _AresMainPy_pkg.f90wrap_iso_ebands(eval=eval, focc=focc)
    return eband

def iso_e_ie(rho, vionlpp):
    """
    eie = iso_e_ie(rho, vionlpp)
    
    
    Defined at Energy_module.fpp lines 676-708
    
    Parameters
    ----------
    rho : float array
    vionlpp : float array
    
    Returns
    -------
    eie : float
    
    """
    eie = _AresMainPy_pkg.f90wrap_iso_e_ie(rho=rho, vionlpp=vionlpp)
    return eie

def lda_energy(rhoreal):
    """
    lda_energy = lda_energy(rhoreal)
    
    
    Defined at Energy_module.fpp lines 712-762
    
    Parameters
    ----------
    rhoreal : float array
    
    Returns
    -------
    lda_energy : float
    
    """
    lda_energy = _AresMainPy_pkg.f90wrap_lda_energy(rhoreal=rhoreal)
    return lda_energy

def get_etot():
    """
    Element etot ftype=real(dp) pytype=float
    
    
    Defined at Energy_module.fpp line 23
    
    """
    return _AresMainPy_pkg.f90wrap_energy_module__get__etot()

def set_etot(etot):
    _AresMainPy_pkg.f90wrap_energy_module__set__etot(etot)

def get_ekine():
    """
    Element ekine ftype=real(dp) pytype=float
    
    
    Defined at Energy_module.fpp line 23
    
    """
    return _AresMainPy_pkg.f90wrap_energy_module__get__ekine()

def set_ekine(ekine):
    _AresMainPy_pkg.f90wrap_energy_module__set__ekine(ekine)

def get_ehart():
    """
    Element ehart ftype=real(dp) pytype=float
    
    
    Defined at Energy_module.fpp line 23
    
    """
    return _AresMainPy_pkg.f90wrap_energy_module__get__ehart()

def set_ehart(ehart):
    _AresMainPy_pkg.f90wrap_energy_module__set__ehart(ehart)

def get_exc():
    """
    Element exc ftype=real(dp) pytype=float
    
    
    Defined at Energy_module.fpp line 23
    
    """
    return _AresMainPy_pkg.f90wrap_energy_module__get__exc()

def set_exc(exc):
    _AresMainPy_pkg.f90wrap_energy_module__set__exc(exc)

def get_eband():
    """
    Element eband ftype=real(dp) pytype=float
    
    
    Defined at Energy_module.fpp line 23
    
    """
    return _AresMainPy_pkg.f90wrap_energy_module__get__eband()

def set_eband(eband):
    _AresMainPy_pkg.f90wrap_energy_module__set__eband(eband)

def get_eie():
    """
    Element eie ftype=real(dp) pytype=float
    
    
    Defined at Energy_module.fpp line 23
    
    """
    return _AresMainPy_pkg.f90wrap_energy_module__get__eie()

def set_eie(eie):
    _AresMainPy_pkg.f90wrap_energy_module__set__eie(eie)

def get_eienl():
    """
    Element eienl ftype=real(dp) pytype=float
    
    
    Defined at Energy_module.fpp line 23
    
    """
    return _AresMainPy_pkg.f90wrap_energy_module__get__eienl()

def set_eienl(eienl):
    _AresMainPy_pkg.f90wrap_energy_module__set__eienl(eienl)

def get_eewald():
    """
    Element eewald ftype=real(dp) pytype=float
    
    
    Defined at Energy_module.fpp line 23
    
    """
    return _AresMainPy_pkg.f90wrap_energy_module__get__eewald()

def set_eewald(eewald):
    _AresMainPy_pkg.f90wrap_energy_module__set__eewald(eewald)

def get_efm():
    """
    Element efm ftype=real(dp) pytype=float
    
    
    Defined at Energy_module.fpp line 23
    
    """
    return _AresMainPy_pkg.f90wrap_energy_module__get__efm()

def set_efm(efm):
    _AresMainPy_pkg.f90wrap_energy_module__set__efm(efm)

def get_tnad():
    """
    Element tnad ftype=real(dp) pytype=float
    
    
    Defined at Energy_module.fpp line 23
    
    """
    return _AresMainPy_pkg.f90wrap_energy_module__get__tnad()

def set_tnad(tnad):
    _AresMainPy_pkg.f90wrap_energy_module__set__tnad(tnad)

def get_fe():
    """
    Element fe ftype=real(dp) pytype=float
    
    
    Defined at Energy_module.fpp line 23
    
    """
    return _AresMainPy_pkg.f90wrap_energy_module__get__fe()

def set_fe(fe):
    _AresMainPy_pkg.f90wrap_energy_module__set__fe(fe)

def get_fe0():
    """
    Element fe0 ftype=real(dp) pytype=float
    
    
    Defined at Energy_module.fpp line 23
    
    """
    return _AresMainPy_pkg.f90wrap_energy_module__get__fe0()

def set_fe0(fe0):
    _AresMainPy_pkg.f90wrap_energy_module__set__fe0(fe0)

def get_wmaxl():
    """
    Element wmaxl ftype=real(dp) pytype=float
    
    
    Defined at Energy_module.fpp line 25
    
    """
    return _AresMainPy_pkg.f90wrap_energy_module__get__wmaxl()

def set_wmaxl(wmaxl):
    _AresMainPy_pkg.f90wrap_energy_module__set__wmaxl(wmaxl)

def get_wminl():
    """
    Element wminl ftype=real(dp) pytype=float
    
    
    Defined at Energy_module.fpp line 25
    
    """
    return _AresMainPy_pkg.f90wrap_energy_module__get__wminl()

def set_wminl(wminl):
    _AresMainPy_pkg.f90wrap_energy_module__set__wminl(wminl)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "energy_module".')

for func in _dt_array_initialisers:
    func()

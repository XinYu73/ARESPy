"""
Module smearing_module


Defined at Smearing_module.fpp lines 5-358

"""
from __future__ import print_function, absolute_import, division
import _arespy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def smear_init(nev):
    """
    smear_init(nev)
    
    
    Defined at Smearing_module.fpp lines 22-35
    
    Parameters
    ----------
    nev : int
    
    """
    _arespy_pkg.f90wrap_smear_init(nev=nev)

def destroy_smear():
    """
    destroy_smear()
    
    
    Defined at Smearing_module.fpp lines 38-46
    
    
    """
    _arespy_pkg.f90wrap_destroy_smear()

def hp(x, n):
    """
    hp = hp(x, n)
    
    
    Defined at Smearing_module.fpp lines 49-79
    
    Parameters
    ----------
    x : float
    n : int
    
    Returns
    -------
    hp : float
    
    """
    hp = _arespy_pkg.f90wrap_hp(x=x, n=n)
    return hp

def smearsn(y, y0, w):
    """
    sn = smearsn(y, y0, w)
    
    
    Defined at Smearing_module.fpp lines 82-120
    
    Parameters
    ----------
    y : float
    y0 : float
    w : float
    
    Returns
    -------
    sn : float
    
    """
    sn = _arespy_pkg.f90wrap_smearsn(y=y, y0=y0, w=w)
    return sn

def fermilevel(ne, nev, nk, wk, eval, sigma):
    """
    fermilevel(ne, nev, nk, wk, eval, sigma)
    
    
    Defined at Smearing_module.fpp lines 123-222
    
    Parameters
    ----------
    ne : float
    nev : int
    nk : int
    wk : float array
    eval : float array
    sigma : float
    
    """
    _arespy_pkg.f90wrap_fermilevel(ne=ne, nev=nev, nk=nk, wk=wk, eval=eval, \
        sigma=sigma)

def enpy(e_mu, fi):
    """
    enpy = enpy(e_mu, fi)
    
    
    Defined at Smearing_module.fpp lines 225-244
    
    Parameters
    ----------
    e_mu : float
    fi : float
    
    Returns
    -------
    enpy : float
    
    """
    enpy = _arespy_pkg.f90wrap_enpy(e_mu=e_mu, fi=fi)
    return enpy

def whg(x, n):
    """
    whg = whg(x, n)
    
    
    Defined at Smearing_module.fpp lines 247-274
    
    Parameters
    ----------
    x : float
    n : int
    
    Returns
    -------
    whg : float
    
    """
    whg = _arespy_pkg.f90wrap_whg(x=x, n=n)
    return whg

def updaterho_pbc(nps, nev, eig, wk, rhos, rho):
    """
    updaterho_pbc(nps, nev, eig, wk, rhos, rho)
    
    
    Defined at Smearing_module.fpp lines 277-323
    
    Parameters
    ----------
    nps : int
    nev : int
    eig : Eigen_Type
    wk : float array
    rhos : float array
    rho : float array
    
    """
    _arespy_pkg.f90wrap_updaterho_pbc(nps=nps, nev=nev, eig=eig._handle, wk=wk, \
        rhos=rhos, rho=rho)

def smear_updaterho(nps, nev, ne, eig, rhos, rho):
    """
    smear_updaterho(nps, nev, ne, eig, rhos, rho)
    
    
    Defined at Smearing_module.fpp lines 326-357
    
    Parameters
    ----------
    nps : int
    nev : int
    ne : float
    eig : Eigen_Type
    rhos : float array
    rho : float array
    
    """
    _arespy_pkg.f90wrap_smear_updaterho(nps=nps, nev=nev, ne=ne, eig=eig._handle, \
        rhos=rhos, rho=rho)

def get_array_wke():
    """
    Element wke ftype=real(dp) pytype=float
    
    
    Defined at Smearing_module.fpp line 17
    
    """
    global wke
    array_ndim, array_type, array_shape, array_handle = \
        _arespy_pkg.f90wrap_smearing_module__array__wke(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        wke = _arrays[array_handle]
    else:
        wke = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_smearing_module__array__wke)
        _arrays[array_handle] = wke
    return wke

def set_array_wke(wke):
    wke[...] = wke

def get_fme():
    """
    Element fme ftype=real(dp) pytype=float
    
    
    Defined at Smearing_module.fpp line 19
    
    """
    return _arespy_pkg.f90wrap_smearing_module__get__fme()

def set_fme(fme):
    _arespy_pkg.f90wrap_smearing_module__set__fme(fme)

def get_ets():
    """
    Element ets ftype=real(dp) pytype=float
    
    
    Defined at Smearing_module.fpp line 19
    
    """
    return _arespy_pkg.f90wrap_smearing_module__get__ets()

def set_ets(ets):
    _arespy_pkg.f90wrap_smearing_module__set__ets(ets)


_array_initialisers = [get_array_wke]
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "smearing_module".')

for func in _dt_array_initialisers:
    func()

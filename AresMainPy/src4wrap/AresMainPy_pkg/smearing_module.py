"""
Module smearing_module


Defined at Smearing_module.fpp lines 5-380

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def smear_init(nev):
    """
    smear_init(nev)
    
    
    Defined at Smearing_module.fpp lines 21-48
    
    Parameters
    ----------
    nev : int
    
    """
    _AresMainPy_pkg.f90wrap_smear_init(nev=nev)

def destroy_smear():
    """
    destroy_smear()
    
    
    Defined at Smearing_module.fpp lines 51-63
    
    
    """
    _AresMainPy_pkg.f90wrap_destroy_smear()

def hp(x, n):
    """
    hp = hp(x, n)
    
    
    Defined at Smearing_module.fpp lines 66-96
    
    Parameters
    ----------
    x : float
    n : int
    
    Returns
    -------
    hp : float
    
    """
    hp = _AresMainPy_pkg.f90wrap_hp(x=x, n=n)
    return hp

def smearsn(y, y0, w):
    """
    sn = smearsn(y, y0, w)
    
    
    Defined at Smearing_module.fpp lines 99-137
    
    Parameters
    ----------
    y : float
    y0 : float
    w : float
    
    Returns
    -------
    sn : float
    
    """
    sn = _AresMainPy_pkg.f90wrap_smearsn(y=y, y0=y0, w=w)
    return sn

def fermilevel(ne, nev, nk, wk, eval, wke):
    """
    fme, ets = fermilevel(ne, nev, nk, wk, eval, wke)
    
    
    Defined at Smearing_module.fpp lines 140-229
    
    Parameters
    ----------
    ne : float
    nev : int
    nk : int
    wk : float array
    eval : float array
    wke : float array
    
    Returns
    -------
    fme : float
    ets : float
    
    """
    fme, ets = _AresMainPy_pkg.f90wrap_fermilevel(ne=ne, nev=nev, nk=nk, wk=wk, \
        eval=eval, wke=wke)
    return fme, ets

def fermilevel_iso(ne, nev, nk, wk, eval, wke):
    """
    fme, ets = fermilevel_iso(ne, nev, nk, wk, eval, wke)
    
    
    Defined at Smearing_module.fpp lines 232-327
    
    Parameters
    ----------
    ne : float
    nev : int
    nk : int
    wk : float array
    eval : float array
    wke : float array
    
    Returns
    -------
    fme : float
    ets : float
    
    """
    fme, ets = _AresMainPy_pkg.f90wrap_fermilevel_iso(ne=ne, nev=nev, nk=nk, wk=wk, \
        eval=eval, wke=wke)
    return fme, ets

def enpy(e_mu, fi):
    """
    enpy = enpy(e_mu, fi)
    
    
    Defined at Smearing_module.fpp lines 330-349
    
    Parameters
    ----------
    e_mu : float
    fi : float
    
    Returns
    -------
    enpy : float
    
    """
    enpy = _AresMainPy_pkg.f90wrap_enpy(e_mu=e_mu, fi=fi)
    return enpy

def whg(x, n):
    """
    whg = whg(x, n)
    
    
    Defined at Smearing_module.fpp lines 352-379
    
    Parameters
    ----------
    x : float
    n : int
    
    Returns
    -------
    whg : float
    
    """
    whg = _AresMainPy_pkg.f90wrap_whg(x=x, n=n)
    return whg

def get_array_wke():
    """
    Element wke ftype=real(dp) pytype=float
    
    
    Defined at Smearing_module.fpp line 16
    
    """
    global wke
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_smearing_module__array__wke(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        wke = _arrays[array_handle]
    else:
        wke = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_smearing_module__array__wke)
        _arrays[array_handle] = wke
    return wke

def set_array_wke(wke):
    wke[...] = wke

def get_array_subwke():
    """
    Element subwke ftype=real(dp) pytype=float
    
    
    Defined at Smearing_module.fpp line 17
    
    """
    global subwke
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_smearing_module__array__subwke(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        subwke = _arrays[array_handle]
    else:
        subwke = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_smearing_module__array__subwke)
        _arrays[array_handle] = subwke
    return subwke

def set_array_subwke(subwke):
    subwke[...] = subwke


_array_initialisers = [get_array_wke, get_array_subwke]
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "smearing_module".')

for func in _dt_array_initialisers:
    func()

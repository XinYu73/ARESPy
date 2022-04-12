"""
Module arpack_module


Defined at Arpack_module.fpp lines 5-919

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def diagh_arpack(veff, ik, nev, evec, eval, resid_restart, nec, info, maxmvs, \
    tol):
    """
    diagh_arpack(veff, ik, nev, evec, eval, resid_restart, nec, info, maxmvs, tol)
    
    
    Defined at Arpack_module.fpp lines 23-192
    
    Parameters
    ----------
    veff : float array
    ik : int
    nev : int
    evec : complex array
    eval : float array
    resid_restart : complex array
    nec : int
    info : int
    maxmvs : int
    tol : float
    
    ------
    """
    _AresMainPy_pkg.f90wrap_diagh_arpack(veff=veff, ik=ik, nev=nev, evec=evec, \
        eval=eval, resid_restart=resid_restart, nec=nec, info=info, maxmvs=maxmvs, \
        tol=tol)

def real_diagh_arpack(veff, nev, evec, eval, resid_restart, nec, info, maxmvs, \
    tol):
    """
    real_diagh_arpack(veff, nev, evec, eval, resid_restart, nec, info, maxmvs, tol)
    
    
    Defined at Arpack_module.fpp lines 202-363
    
    Parameters
    ----------
    veff : float array
    nev : int
    evec : float array
    eval : float array
    resid_restart : float array
    nec : int
    info : int
    maxmvs : int
    tol : float
    
    ------
    """
    _AresMainPy_pkg.f90wrap_real_diagh_arpack(veff=veff, nev=nev, evec=evec, \
        eval=eval, resid_restart=resid_restart, nec=nec, info=info, maxmvs=maxmvs, \
        tol=tol)

def real_diagm_arpk(mat, nev, evec, eval, resid_restart, nec, info, maxmvs, \
    tol):
    """
    real_diagm_arpk(mat, nev, evec, eval, resid_restart, nec, info, maxmvs, tol)
    
    
    Defined at Arpack_module.fpp lines 373-506
    
    Parameters
    ----------
    mat : float array
    nev : int
    evec : float array
    eval : float array
    resid_restart : float array
    nec : int
    info : int
    maxmvs : int
    tol : float
    
    ------
    """
    _AresMainPy_pkg.f90wrap_real_diagm_arpk(mat=mat, nev=nev, evec=evec, eval=eval, \
        resid_restart=resid_restart, nec=nec, info=info, maxmvs=maxmvs, tol=tol)

def rdiagm_arpk(mat, dimen, nev, evec, eval):
    """
    rdiagm_arpk(mat, dimen, nev, evec, eval)
    
    
    Defined at Arpack_module.fpp lines 509-544
    
    Parameters
    ----------
    mat : float array
    dimen : int
    nev : int
    evec : float array
    eval : float array
    
    """
    _AresMainPy_pkg.f90wrap_rdiagm_arpk(mat=mat, dimen=dimen, nev=nev, evec=evec, \
        eval=eval)

def iso_diagh_arpack(veff, nev, evec, eval, resid_restart, nec, info, maxmvs, \
    tol):
    """
    iso_diagh_arpack(veff, nev, evec, eval, resid_restart, nec, info, maxmvs, tol)
    
    
    Defined at Arpack_module.fpp lines 553-716
    
    Parameters
    ----------
    veff : float array
    nev : int
    evec : float array
    eval : float array
    resid_restart : float array
    nec : int
    info : int
    maxmvs : int
    tol : float
    
    ------
    """
    _AresMainPy_pkg.f90wrap_iso_diagh_arpack(veff=veff, nev=nev, evec=evec, \
        eval=eval, resid_restart=resid_restart, nec=nec, info=info, maxmvs=maxmvs, \
        tol=tol)

def diagh_arpack_band(veff, ik, nev, evec, eval, resid_restart, nec, info, \
    maxmvs, tol):
    """
    diagh_arpack_band(veff, ik, nev, evec, eval, resid_restart, nec, info, maxmvs, \
        tol)
    
    
    Defined at Arpack_module.fpp lines 719-918
    
    Parameters
    ----------
    veff : float array
    ik : int
    nev : int
    evec : complex array
    eval : float array
    resid_restart : complex array
    nec : int
    info : int
    maxmvs : int
    tol : float
    
    ------
    """
    _AresMainPy_pkg.f90wrap_diagh_arpack_band(veff=veff, ik=ik, nev=nev, evec=evec, \
        eval=eval, resid_restart=resid_restart, nec=nec, info=info, maxmvs=maxmvs, \
        tol=tol)

def get_maxn():
    """
    Element maxn ftype=integer(i4b) pytype=int
    
    
    Defined at Arpack_module.fpp line 14
    
    """
    return _AresMainPy_pkg.f90wrap_arpack_module__get__maxn()

maxn = get_maxn()

def get_maxnev():
    """
    Element maxnev ftype=integer(i4b) pytype=int
    
    
    Defined at Arpack_module.fpp line 14
    
    """
    return _AresMainPy_pkg.f90wrap_arpack_module__get__maxnev()

maxnev = get_maxnev()

def get_maxncv():
    """
    Element maxncv ftype=integer(i4b) pytype=int
    
    
    Defined at Arpack_module.fpp line 14
    
    """
    return _AresMainPy_pkg.f90wrap_arpack_module__get__maxncv()

maxncv = get_maxncv()


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "arpack_module".')

for func in _dt_array_initialisers:
    func()

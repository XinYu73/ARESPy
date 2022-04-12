"""
Module arpack_module


Defined at Arpack_module.fpp lines 5-920

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
    
    
    Defined at Arpack_module.fpp lines 24-193
    
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
    
    
    Defined at Arpack_module.fpp lines 203-364
    
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
    
    
    Defined at Arpack_module.fpp lines 374-507
    
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
    
    
    Defined at Arpack_module.fpp lines 510-545
    
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
    
    
    Defined at Arpack_module.fpp lines 554-717
    
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
    
    
    Defined at Arpack_module.fpp lines 720-919
    
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

def set_maxn(maxn):
    _AresMainPy_pkg.f90wrap_arpack_module__set__maxn(maxn)

def get_maxnev():
    """
    Element maxnev ftype=integer(i4b) pytype=int
    
    
    Defined at Arpack_module.fpp line 14
    
    """
    return _AresMainPy_pkg.f90wrap_arpack_module__get__maxnev()

def set_maxnev(maxnev):
    _AresMainPy_pkg.f90wrap_arpack_module__set__maxnev(maxnev)

def get_maxncv():
    """
    Element maxncv ftype=integer(i4b) pytype=int
    
    
    Defined at Arpack_module.fpp line 14
    
    """
    return _AresMainPy_pkg.f90wrap_arpack_module__get__maxncv()

def set_maxncv(maxncv):
    _AresMainPy_pkg.f90wrap_arpack_module__set__maxncv(maxncv)


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

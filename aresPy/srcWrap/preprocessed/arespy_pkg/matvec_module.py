"""
Module matvec_module


Defined at Matvec_module.fpp lines 5-195

"""
from __future__ import print_function, absolute_import, division
import _arespy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def cmplx_matvec(ik, veff1d, p, q, dimen):
    """
    cmplx_matvec(ik, veff1d, p, q, dimen)
    
    
    Defined at Matvec_module.fpp lines 10-28
    
    Parameters
    ----------
    ik : int
    veff1d : float array
    p : complex array
    q : complex array
    dimen : int
    
    """
    _arespy_pkg.f90wrap_cmplx_matvec(ik=ik, veff1d=veff1d, p=p, q=q, dimen=dimen)

def cmplx_nlocmatvec(ik, p, q):
    """
    cmplx_nlocmatvec(ik, p, q)
    
    
    Defined at Matvec_module.fpp lines 31-94
    
    Parameters
    ----------
    ik : int
    p : complex array
    q : complex array
    
    """
    _arespy_pkg.f90wrap_cmplx_nlocmatvec(ik=ik, p=p, q=q)

def real_matvec(veff1d, p, q, dimen):
    """
    real_matvec(veff1d, p, q, dimen)
    
    
    Defined at Matvec_module.fpp lines 97-117
    
    Parameters
    ----------
    veff1d : float array
    p : float array
    q : float array
    dimen : int
    
    """
    _arespy_pkg.f90wrap_real_matvec(veff1d=veff1d, p=p, q=q, dimen=dimen)

def real_nlocmatvec(p, q):
    """
    real_nlocmatvec(p, q)
    
    
    Defined at Matvec_module.fpp lines 120-183
    
    Parameters
    ----------
    p : float array
    q : float array
    
    """
    _arespy_pkg.f90wrap_real_nlocmatvec(p=p, q=q)

def real_matvec_m(mat, p, q, dimen):
    """
    real_matvec_m(mat, p, q, dimen)
    
    
    Defined at Matvec_module.fpp lines 186-194
    
    Parameters
    ----------
    mat : float array
    p : float array
    q : float array
    dimen : int
    
    """
    _arespy_pkg.f90wrap_real_matvec_m(mat=mat, p=p, q=q, dimen=dimen)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "matvec_module".')

for func in _dt_array_initialisers:
    func()

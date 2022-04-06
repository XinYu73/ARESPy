"""
Module potential_module


Defined at Potential_module.fpp lines 5-190

"""
from __future__ import print_function, absolute_import, division
import _arespy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def calveff(nps, rhos, rho, veffs):
    """
    calveff(nps, rhos, rho, veffs)
    
    
    Defined at Potential_module.fpp lines 19-56
    
    Parameters
    ----------
    nps : int
    rhos : float array
    rho : float array
    veffs : float array
    
    """
    _arespy_pkg.f90wrap_calveff(nps=nps, rhos=rhos, rho=rho, veffs=veffs)

def vhartree(nps, rho, vhart):
    """
    vhartree(nps, rho, vhart)
    
    
    Defined at Potential_module.fpp lines 59-124
    
    Parameters
    ----------
    nps : int
    rho : float array
    vhart : float array
    
    """
    _arespy_pkg.f90wrap_vhartree(nps=nps, rho=rho, vhart=vhart)

def vlpp():
    """
    vlpp()
    
    
    Defined at Potential_module.fpp lines 127-189
    
    
    """
    _arespy_pkg.f90wrap_vlpp()

def get_array_v_accelerate():
    """
    Element v_accelerate ftype=real(dp) pytype=float
    
    
    Defined at Potential_module.fpp line 13
    
    """
    global v_accelerate
    array_ndim, array_type, array_shape, array_handle = \
        _arespy_pkg.f90wrap_potential_module__array__v_accelerate(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        v_accelerate = _arrays[array_handle]
    else:
        v_accelerate = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_potential_module__array__v_accelerate)
        _arrays[array_handle] = v_accelerate
    return v_accelerate

def set_array_v_accelerate(v_accelerate):
    v_accelerate[...] = v_accelerate


_array_initialisers = [get_array_v_accelerate]
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "potential_module".')

for func in _dt_array_initialisers:
    func()

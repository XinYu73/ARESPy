"""
Module potential_module


Defined at Potential_module.fpp lines 5-436

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def cal_veff(rhos, veffs):
    """
    cal_veff(rhos, veffs)
    
    
    Defined at Potential_module.fpp lines 19-73
    
    Parameters
    ----------
    rhos : float array
    veffs : float array
    
    """
    _AresMainPy_pkg.f90wrap_cal_veff(rhos=rhos, veffs=veffs)

def cal_veff_iso(rhos, veffs):
    """
    cal_veff_iso(rhos, veffs)
    
    
    Defined at Potential_module.fpp lines 76-125
    
    Parameters
    ----------
    rhos : float array
    veffs : float array
    
    """
    _AresMainPy_pkg.f90wrap_cal_veff_iso(rhos=rhos, veffs=veffs)

def vhartree(rho, vhart):
    """
    vhartree(rho, vhart)
    
    
    Defined at Potential_module.fpp lines 128-205
    
    Parameters
    ----------
    rho : float array
    vhart : float array
    
    """
    _AresMainPy_pkg.f90wrap_vhartree(rho=rho, vhart=vhart)

def vlpp():
    """
    vlpp()
    
    
    Defined at Potential_module.fpp lines 208-276
    
    
    """
    _AresMainPy_pkg.f90wrap_vlpp()

def vlpp_real():
    """
    vlpp_real()
    
    
    Defined at Potential_module.fpp lines 279-324
    
    
    """
    _AresMainPy_pkg.f90wrap_vlpp_real()

def vlda(rho, ldapotential):
    """
    vlda(rho, ldapotential)
    
    
    Defined at Potential_module.fpp lines 327-403
    
    Parameters
    ----------
    rho : float array
    ldapotential : float array
    
    ------------------------------------------------------------------------------
     DESCRIPTION:
       This function computes the exchange-correlation potential in the Local
       Density Approximation(LDA) based on the real-space electron density.
     GLOBAL/MODULE VARIABLES CHANGED:
     CONDITIONS AND ASSUMPTIONS:
       We can use LDAPointPot for moudularity.
     FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
       Complete calculation to handle spin-polarized cases.
     REFERENCES:
       [1]  Perdew J.P. and Zunger A., Phys. Rev. B 23(10), 1981, 5048-79
    ------------------------------------------------------------------------------
    """
    _AresMainPy_pkg.f90wrap_vlda(rho=rho, ldapotential=ldapotential)

def sumrhos(rhos, rho):
    """
    sumrhos(rhos, rho)
    
    
    Defined at Potential_module.fpp lines 406-419
    
    Parameters
    ----------
    rhos : float array
    rho : float array
    
    """
    _AresMainPy_pkg.f90wrap_sumrhos(rhos=rhos, rho=rho)

def sumrhos_iso(rhos, rho):
    """
    sumrhos_iso(rhos, rho)
    
    
    Defined at Potential_module.fpp lines 422-435
    
    Parameters
    ----------
    rhos : float array
    rho : float array
    
    """
    _AresMainPy_pkg.f90wrap_sumrhos_iso(rhos=rhos, rho=rho)

def get_array_v_accelerate():
    """
    Element v_accelerate ftype=real(dp) pytype=float
    
    
    Defined at Potential_module.fpp line 13
    
    """
    global v_accelerate
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_potential_module__array__v_accelerate(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        v_accelerate = _arrays[array_handle]
    else:
        v_accelerate = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_potential_module__array__v_accelerate)
        _arrays[array_handle] = v_accelerate
    return v_accelerate

def set_array_v_accelerate(v_accelerate):
    v_accelerate[...] = v_accelerate

def get_array_v_hxc_old():
    """
    Element v_hxc_old ftype=real(dp) pytype=float
    
    
    Defined at Potential_module.fpp line 14
    
    """
    global v_hxc_old
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_potential_module__array__v_hxc_old(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        v_hxc_old = _arrays[array_handle]
    else:
        v_hxc_old = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_potential_module__array__v_hxc_old)
        _arrays[array_handle] = v_hxc_old
    return v_hxc_old

def set_array_v_hxc_old(v_hxc_old):
    v_hxc_old[...] = v_hxc_old

def get_array_v_hxc_new():
    """
    Element v_hxc_new ftype=real(dp) pytype=float
    
    
    Defined at Potential_module.fpp line 14
    
    """
    global v_hxc_new
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_potential_module__array__v_hxc_new(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        v_hxc_new = _arrays[array_handle]
    else:
        v_hxc_new = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_potential_module__array__v_hxc_new)
        _arrays[array_handle] = v_hxc_new
    return v_hxc_new

def set_array_v_hxc_new(v_hxc_new):
    v_hxc_new[...] = v_hxc_new

def get_accelerate():
    """
    Element accelerate ftype=logical pytype=bool
    
    
    Defined at Potential_module.fpp line 15
    
    """
    return _AresMainPy_pkg.f90wrap_potential_module__get__accelerate()

def set_accelerate(accelerate):
    _AresMainPy_pkg.f90wrap_potential_module__set__accelerate(accelerate)


_array_initialisers = [get_array_v_accelerate, get_array_v_hxc_old, \
    get_array_v_hxc_new]
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "potential_module".')

for func in _dt_array_initialisers:
    func()

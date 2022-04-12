"""
Module libxc_module


Defined at XC_functional.fpp lines 5-344

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def ldalib_energy(rhos):
    """
    exc = ldalib_energy(rhos)
    
    
    Defined at XC_functional.fpp lines 19-138
    
    Parameters
    ----------
    rhos : float array
    
    Returns
    -------
    exc : float
    
    """
    exc = _AresMainPy_pkg.f90wrap_ldalib_energy(rhos=rhos)
    return exc

def ldalib_potential(rhos, vxcs):
    """
    ldalib_potential(rhos, vxcs)
    
    
    Defined at XC_functional.fpp lines 141-249
    
    Parameters
    ----------
    rhos : float array
    vxcs : float array
    
    """
    _AresMainPy_pkg.f90wrap_ldalib_potential(rhos=rhos, vxcs=vxcs)

def ldalib_energy_iso(rhos):
    """
    exc = ldalib_energy_iso(rhos)
    
    
    Defined at XC_functional.fpp lines 252-298
    
    Parameters
    ----------
    rhos : float array
    
    Returns
    -------
    exc : float
    
    """
    exc = _AresMainPy_pkg.f90wrap_ldalib_energy_iso(rhos=rhos)
    return exc

def ldalib_potential_iso(rhos, vxcs):
    """
    ldalib_potential_iso(rhos, vxcs)
    
    
    Defined at XC_functional.fpp lines 301-343
    
    Parameters
    ----------
    rhos : float array
    vxcs : float array
    
    """
    _AresMainPy_pkg.f90wrap_ldalib_potential_iso(rhos=rhos, vxcs=vxcs)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "libxc_module".')

for func in _dt_array_initialisers:
    func()

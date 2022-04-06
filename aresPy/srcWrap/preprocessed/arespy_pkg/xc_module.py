"""
Module xc_module


Defined at XC_functional.fpp lines 9-428

"""
from __future__ import print_function, absolute_import, division
import _arespy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def xc_functional(nps, rhos, vxc, exc=None):
    """
    xc_functional(nps, rhos, vxc[, exc])
    
    
    Defined at XC_functional.fpp lines 21-46
    
    Parameters
    ----------
    nps : int
    rhos : float array
    vxc : float array
    exc : float
    
    """
    _arespy_pkg.f90wrap_xc_functional(nps=nps, rhos=rhos, vxc=vxc, exc=exc)

def libxc_lda_set(nps, rhos, vxcs, exc=None):
    """
    libxc_lda_set(nps, rhos, vxcs[, exc])
    
    
    Defined at XC_functional.fpp lines 51-112
    
    Parameters
    ----------
    nps : int
    rhos : float array
    vxcs : float array
    exc : float
    
    """
    _arespy_pkg.f90wrap_libxc_lda_set(nps=nps, rhos=rhos, vxcs=vxcs, exc=exc)

def libxc_lda_x(nps, rhos, vxs=None):
    """
    ex = libxc_lda_x(nps, rhos[, vxs])
    
    
    Defined at XC_functional.fpp lines 117-160
    
    Parameters
    ----------
    nps : int
    rhos : float array
    vxs : float array
    
    Returns
    -------
    ex : float
    
    """
    ex = _arespy_pkg.f90wrap_libxc_lda_x(nps=nps, rhos=rhos, vxs=vxs)
    return ex

def libxc_vwn1rpa_c(nps, rhos, vcs=None):
    """
    ec = libxc_vwn1rpa_c(nps, rhos[, vcs])
    
    
    Defined at XC_functional.fpp lines 165-208
    
    Parameters
    ----------
    nps : int
    rhos : float array
    vcs : float array
    
    Returns
    -------
    ec : float
    
    """
    ec = _arespy_pkg.f90wrap_libxc_vwn1rpa_c(nps=nps, rhos=rhos, vcs=vcs)
    return ec

def libxc_gga_set(nps, rhos, vxcs, exc=None):
    """
    libxc_gga_set(nps, rhos, vxcs[, exc])
    
    
    Defined at XC_functional.fpp lines 213-319
    
    Parameters
    ----------
    nps : int
    rhos : float array
    vxcs : float array
    exc : float
    
    """
    _arespy_pkg.f90wrap_libxc_gga_set(nps=nps, rhos=rhos, vxcs=vxcs, exc=exc)

def libxc_b88_x(nps, rhos, sigma, vxs=None):
    """
    ex = libxc_b88_x(nps, rhos, sigma[, vxs])
    
    
    Defined at XC_functional.fpp lines 324-371
    
    Parameters
    ----------
    nps : int
    rhos : float array
    sigma : float array
    vxs : float array
    
    Returns
    -------
    ex : float
    
    """
    ex = _arespy_pkg.f90wrap_libxc_b88_x(nps=nps, rhos=rhos, sigma=sigma, vxs=vxs)
    return ex

def libxc_lyp_c(nps, rhos, sigma, vcs=None):
    """
    ec = libxc_lyp_c(nps, rhos, sigma[, vcs])
    
    
    Defined at XC_functional.fpp lines 376-424
    
    Parameters
    ----------
    nps : int
    rhos : float array
    sigma : float array
    vcs : float array
    
    Returns
    -------
    ec : float
    
    """
    ec = _arespy_pkg.f90wrap_libxc_lyp_c(nps=nps, rhos=rhos, sigma=sigma, vcs=vcs)
    return ec


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "xc_module".')

for func in _dt_array_initialisers:
    func()

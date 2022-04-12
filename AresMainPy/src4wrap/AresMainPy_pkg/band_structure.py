"""
Module band_structure


Defined at Bands_module.fpp lines 10-405

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def init_bandstruct(numk, kvec):
    """
    init_bandstruct(numk, kvec)
    
    
    Defined at Bands_module.fpp lines 29-49
    
    Parameters
    ----------
    numk : int
    kvec : float array
    
    """
    _AresMainPy_pkg.f90wrap_init_bandstruct(numk=numk, kvec=kvec)

def band_begin(numk, kvec):
    """
    band_begin(numk, kvec)
    
    
    Defined at Bands_module.fpp lines 52-107
    
    Parameters
    ----------
    numk : int
    kvec : float array
    
    """
    _AresMainPy_pkg.f90wrap_band_begin(numk=numk, kvec=kvec)

def build_bands():
    """
    build_bands()
    
    
    Defined at Bands_module.fpp lines 110-171
    
    
    """
    _AresMainPy_pkg.f90wrap_build_bands()

def read_bandpath(infile):
    """
    read_bandpath(infile)
    
    
    Defined at Bands_module.fpp lines 174-268
    
    Parameters
    ----------
    infile : str
    
    -------start input.dat-----------
    """
    _AresMainPy_pkg.f90wrap_read_bandpath(infile=infile)

def band_pathread(infile):
    """
    band_pathread(infile)
    
    
    Defined at Bands_module.fpp lines 271-304
    
    Parameters
    ----------
    infile : str
    
    """
    _AresMainPy_pkg.f90wrap_band_pathread(infile=infile)

def read_density(infile, rho):
    """
    read_density(infile, rho)
    
    
    Defined at Bands_module.fpp lines 307-351
    
    Parameters
    ----------
    infile : str
    rho : float array
    
    """
    _AresMainPy_pkg.f90wrap_read_density(infile=infile, rho=rho)

def cal_band(rhos, nev, eigval):
    """
    cal_band(rhos, nev, eigval)
    
    
    Defined at Bands_module.fpp lines 354-404
    
    Parameters
    ----------
    rhos : float array
    nev : int
    eigval : float array
    
    """
    _AresMainPy_pkg.f90wrap_cal_band(rhos=rhos, nev=nev, eigval=eigval)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "band_structure".')

for func in _dt_array_initialisers:
    func()

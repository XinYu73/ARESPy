"""
Module read_module


Defined at Read_module.fpp lines 5-629

"""
from __future__ import print_function, absolute_import, division
import _arespy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def read_file(infile):
    """
    read_file(infile)
    
    
    Defined at Read_module.fpp lines 16-448
    
    Parameters
    ----------
    infile : str
    
    -----------------------------------------------------------
    """
    _arespy_pkg.f90wrap_read_file(infile=infile)

def read_poscar(nty, filename):
    """
    read_poscar(nty, filename)
    
    
    Defined at Read_module.fpp lines 451-566
    
    Parameters
    ----------
    nty : int
    filename : str
    
    """
    _arespy_pkg.f90wrap_read_poscar(nty=nty, filename=filename)

def resetpos(natom, lat, pos, poscar):
    """
    resetpos(natom, lat, pos, poscar)
    
    
    Defined at Read_module.fpp lines 569-628
    
    Parameters
    ----------
    natom : int
    lat : float array
    pos : float array
    poscar : float array
    
    """
    _arespy_pkg.f90wrap_resetpos(natom=natom, lat=lat, pos=pos, poscar=poscar)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "read_module".')

for func in _dt_array_initialisers:
    func()

"""
Module begin_module


Defined at Begin_module.fpp lines 5-203

"""
from __future__ import print_function, absolute_import, division
import _arespy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def initial_grid_pbc():
    """
    initial_grid_pbc()
    
    
    Defined at Begin_module.fpp lines 10-72
    
    
    """
    _arespy_pkg.f90wrap_initial_grid_pbc()

def initial_density():
    """
    initial_density()
    
    
    Defined at Begin_module.fpp lines 75-141
    
    
    """
    _arespy_pkg.f90wrap_initial_density()

def inichrg_sp():
    """
    inichrg_sp()
    
    
    Defined at Begin_module.fpp lines 144-202
    
    
    """
    _arespy_pkg.f90wrap_inichrg_sp()


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "begin_module".')

for func in _dt_array_initialisers:
    func()

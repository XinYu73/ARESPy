"""
Module end_module


Defined at End_module.fpp lines 5-29

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def destroy_beast():
    """
    destroy_beast()
    
    
    Defined at End_module.fpp lines 13-28
    
    
    """
    _AresMainPy_pkg.f90wrap_destroy_beast()


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "end_module".')

for func in _dt_array_initialisers:
    func()

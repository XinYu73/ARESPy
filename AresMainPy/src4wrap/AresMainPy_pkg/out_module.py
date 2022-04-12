"""
Module out_module


Defined at Out_module.fpp lines 5-130

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def ksout_list():
    """
    ksout_list()
    
    
    Defined at Out_module.fpp lines 9-16
    
    
    """
    _AresMainPy_pkg.f90wrap_ksout_list()

def ksout_list_per():
    """
    ksout_list_per()
    
    
    Defined at Out_module.fpp lines 19-69
    
    
    """
    _AresMainPy_pkg.f90wrap_ksout_list_per()

def ksout_list_iso():
    """
    ksout_list_iso()
    
    
    Defined at Out_module.fpp lines 72-130
    
    
    """
    _AresMainPy_pkg.f90wrap_ksout_list_iso()


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "out_module".')

for func in _dt_array_initialisers:
    func()

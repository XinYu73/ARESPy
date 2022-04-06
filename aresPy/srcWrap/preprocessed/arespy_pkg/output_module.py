"""
Module output_module


Defined at Energy_module.fpp lines 95-192

"""
from __future__ import print_function, absolute_import, division
import _arespy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def output():
    """
    output()
    
    
    Defined at Energy_module.fpp lines 105-142
    
    
    """
    _arespy_pkg.f90wrap_output()

def write_density():
    """
    write_density()
    
    
    Defined at Energy_module.fpp lines 145-168
    
    
    """
    _arespy_pkg.f90wrap_write_density()

def write_band():
    """
    write_band()
    
    
    Defined at Energy_module.fpp lines 172-191
    
    
    """
    _arespy_pkg.f90wrap_write_band()

def get_time_total0():
    """
    Element time_total0 ftype=integer(i4b) pytype=int
    
    
    Defined at Energy_module.fpp line 102
    
    """
    return _arespy_pkg.f90wrap_output_module__get__time_total0()

def set_time_total0(time_total0):
    _arespy_pkg.f90wrap_output_module__set__time_total0(time_total0)

def get_time_total1():
    """
    Element time_total1 ftype=integer(i4b) pytype=int
    
    
    Defined at Energy_module.fpp line 102
    
    """
    return _arespy_pkg.f90wrap_output_module__get__time_total1()

def set_time_total1(time_total1):
    _arespy_pkg.f90wrap_output_module__set__time_total1(time_total1)

def get_time_scf0():
    """
    Element time_scf0 ftype=integer(i4b) pytype=int
    
    
    Defined at Energy_module.fpp line 102
    
    """
    return _arespy_pkg.f90wrap_output_module__get__time_scf0()

def set_time_scf0(time_scf0):
    _arespy_pkg.f90wrap_output_module__set__time_scf0(time_scf0)

def get_time_scf1():
    """
    Element time_scf1 ftype=integer(i4b) pytype=int
    
    
    Defined at Energy_module.fpp line 102
    
    """
    return _arespy_pkg.f90wrap_output_module__get__time_scf1()

def set_time_scf1(time_scf1):
    _arespy_pkg.f90wrap_output_module__set__time_scf1(time_scf1)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "output_module".')

for func in _dt_array_initialisers:
    func()

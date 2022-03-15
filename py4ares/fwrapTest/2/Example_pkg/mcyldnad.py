"""
Module mcyldnad


Defined at cyldnad.fpp lines 5-14

"""
from __future__ import print_function, absolute_import, division
import _Example_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def cyldnad(self, height):
    """
    vol = cyldnad(self, height)
    
    
    Defined at cyldnad.fpp lines 8-13
    
    Parameters
    ----------
    radius : Dual_Num
    height : Dual_Num
    
    Returns
    -------
    vol : Dual_Num
    
    """
    vol = _Example_pkg.f90wrap_cyldnad(radius=self._handle, height=height._handle)
    vol = f90wrap.runtime.lookup_class("Example_pkg.DUAL_NUM").from_handle(vol, \
        alloc=True)
    return vol


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "mcyldnad".')

for func in _dt_array_initialisers:
    func()

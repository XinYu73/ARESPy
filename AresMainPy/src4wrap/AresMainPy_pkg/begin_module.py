"""
Module begin_module


Defined at Begin_module.fpp lines 5-650

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def initial_grid():
    """
    initial_grid()
    
    
    Defined at Begin_module.fpp lines 11-18
    
    
    """
    _AresMainPy_pkg.f90wrap_initial_grid()

def initial_grid_per():
    """
    initial_grid_per()
    
    
    Defined at Begin_module.fpp lines 21-80
    
    
    """
    _AresMainPy_pkg.f90wrap_initial_grid_per()

def initial_grid_iso():
    """
    initial_grid_iso()
    
    
    Defined at Begin_module.fpp lines 83-143
    
    
    """
    _AresMainPy_pkg.f90wrap_initial_grid_iso()

def initial_density():
    """
    initial_density()
    
    
    Defined at Begin_module.fpp lines 146-185
    
    
    =================================================
    SHIELD
    """
    _AresMainPy_pkg.f90wrap_initial_density()

def inichrg():
    """
    inichrg()
    
    
    Defined at Begin_module.fpp lines 188-375
    
    
    """
    _AresMainPy_pkg.f90wrap_inichrg()

def init_iso_density():
    """
    init_iso_density()
    
    
    Defined at Begin_module.fpp lines 378-405
    
    
    =====================================================================
    ##JUDGE THE INITIALIZE METHOD
    """
    _AresMainPy_pkg.f90wrap_init_iso_density()

def inichrg_iso():
    """
    inichrg_iso()
    
    
    Defined at Begin_module.fpp lines 408-503
    
    
    """
    _AresMainPy_pkg.f90wrap_inichrg_iso()

def init_mo_density_iso():
    """
    init_mo_density_iso()
    
    
    Defined at Begin_module.fpp lines 506-569
    
    
    ===============================================================
    ##assignment to grid
    """
    _AresMainPy_pkg.f90wrap_init_mo_density_iso()

def iso_init_mo(initx_sto):
    """
    iso_init_mo(initx_sto)
    
    
    Defined at Begin_module.fpp lines 572-648
    
    Parameters
    ----------
    initx_sto : float array
    
    """
    _AresMainPy_pkg.f90wrap_iso_init_mo(initx_sto=initx_sto)


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

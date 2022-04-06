"""
Module energy_module


Defined at Energy_module.fpp lines 11-80

"""
from __future__ import print_function, absolute_import, division
import _arespy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def totalenergy(nps, eig, rhos, rho):
    """
    totalenergy(nps, eig, rhos, rho)
    
    
    Defined at Energy_module.fpp lines 26-55
    
    Parameters
    ----------
    nps : int
    eig : Eigen_Type
    rhos : float array
    rho : float array
    
    """
    _arespy_pkg.f90wrap_totalenergy(nps=nps, eig=eig._handle, rhos=rhos, rho=rho)

def ebands(eval, wke):
    """
    eband = ebands(eval, wke)
    
    
    Defined at Energy_module.fpp lines 58-78
    
    Parameters
    ----------
    eval : float array
    wke : float array
    
    Returns
    -------
    eband : float
    
    """
    eband = _arespy_pkg.f90wrap_ebands(eval=eval, wke=wke)
    return eband

def get_etot():
    """
    Element etot ftype=real(dp) pytype=float
    
    
    Defined at Energy_module.fpp line 23
    
    """
    return _arespy_pkg.f90wrap_energy_module__get__etot()

def set_etot(etot):
    _arespy_pkg.f90wrap_energy_module__set__etot(etot)

def get_eband():
    """
    Element eband ftype=real(dp) pytype=float
    
    
    Defined at Energy_module.fpp line 23
    
    """
    return _arespy_pkg.f90wrap_energy_module__get__eband()

def set_eband(eband):
    _arespy_pkg.f90wrap_energy_module__set__eband(eband)

def get_eh():
    """
    Element eh ftype=real(dp) pytype=float
    
    
    Defined at Energy_module.fpp line 23
    
    """
    return _arespy_pkg.f90wrap_energy_module__get__eh()

def set_eh(eh):
    _arespy_pkg.f90wrap_energy_module__set__eh(eh)

def get_eext():
    """
    Element eext ftype=real(dp) pytype=float
    
    
    Defined at Energy_module.fpp line 23
    
    """
    return _arespy_pkg.f90wrap_energy_module__get__eext()

def set_eext(eext):
    _arespy_pkg.f90wrap_energy_module__set__eext(eext)

def get_exc():
    """
    Element exc ftype=real(dp) pytype=float
    
    
    Defined at Energy_module.fpp line 23
    
    """
    return _arespy_pkg.f90wrap_energy_module__get__exc()

def set_exc(exc):
    _arespy_pkg.f90wrap_energy_module__set__exc(exc)

def get_eele():
    """
    Element eele ftype=real(dp) pytype=float
    
    
    Defined at Energy_module.fpp line 23
    
    """
    return _arespy_pkg.f90wrap_energy_module__get__eele()

def set_eele(eele):
    _arespy_pkg.f90wrap_energy_module__set__eele(eele)

def get_fe():
    """
    Element fe ftype=real(dp) pytype=float
    
    
    Defined at Energy_module.fpp line 23
    
    """
    return _arespy_pkg.f90wrap_energy_module__get__fe()

def set_fe(fe):
    _arespy_pkg.f90wrap_energy_module__set__fe(fe)

def get_fe0():
    """
    Element fe0 ftype=real(dp) pytype=float
    
    
    Defined at Energy_module.fpp line 23
    
    """
    return _arespy_pkg.f90wrap_energy_module__get__fe0()

def set_fe0(fe0):
    _arespy_pkg.f90wrap_energy_module__set__fe0(fe0)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "energy_module".')

for func in _dt_array_initialisers:
    func()

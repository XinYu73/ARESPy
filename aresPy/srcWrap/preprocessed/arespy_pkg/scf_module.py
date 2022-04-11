"""
Module scf_module


Defined at Scf_module.fpp lines 5-299

"""
from __future__ import print_function, absolute_import, division
import _arespy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def electronicscf():
    """
    electronicscf()
    
    
    Defined at Scf_module.fpp lines 18-43
    
    
    """
    _arespy_pkg.f90wrap_electronicscf()

def arpackscf(nps, rhos, rho, eig):
    """
    arpackscf(nps, rhos, rho, eig)
    
    
    Defined at Scf_module.fpp lines 50-119
    
    Parameters
    ----------
    nps : int
    rhos : float array
    rho : float array
    eig : Eigen_Type
    
    ===============================================================
    """
    _arespy_pkg.f90wrap_arpackscf(nps=nps, rhos=rhos, rho=rho, eig=eig._handle)

def eigensolver_real(nps, nev, veff, psi, eval, diagtol):
    """
    eigensolver_real(nps, nev, veff, psi, eval, diagtol)
    
    
    Defined at Scf_module.fpp lines 122-165
    
    Parameters
    ----------
    nps : int
    nev : int
    veff : float array
    psi : float array
    eval : float array
    diagtol : float
    
    """
    _arespy_pkg.f90wrap_eigensolver_real(nps=nps, nev=nev, veff=veff, psi=psi, \
        eval=eval, diagtol=diagtol)

def chefsi(nps, rhos, rho, eig):
    """
    chefsi(nps, rhos, rho, eig)
    
    
    Defined at Scf_module.fpp lines 173-280
    
    Parameters
    ----------
    nps : int
    rhos : float array
    rho : float array
    eig : Eigen_Type
    
    ===============================================================
    """
    _arespy_pkg.f90wrap_chefsi(nps=nps, rhos=rhos, rho=rho, eig=eig._handle)

def filter_spin_gamma(nps, nev, veff, x, d):
    """
    filter_spin_gamma(nps, nev, veff, x, d)
    
    
    Defined at Scf_module.fpp lines 283-298
    
    Parameters
    ----------
    nps : int
    nev : int
    veff : float array
    x : float array
    d : float array
    
    """
    _arespy_pkg.f90wrap_filter_spin_gamma(nps=nps, nev=nev, veff=veff, x=x, d=d)

def get_iwd():
    """
    Element iwd ftype=integer(i4b) pytype=int
    
    
    Defined at Scf_module.fpp line 14
    
    """
    return _arespy_pkg.f90wrap_scf_module__get__iwd()

def set_iwd(iwd):
    _arespy_pkg.f90wrap_scf_module__set__iwd(iwd)

def get_lscf():
    """
    Element lscf ftype=logical pytype=bool
    
    
    Defined at Scf_module.fpp line 15
    
    """
    return _arespy_pkg.f90wrap_scf_module__get__lscf()

def set_lscf(lscf):
    _arespy_pkg.f90wrap_scf_module__set__lscf(lscf)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "scf_module".')

for func in _dt_array_initialisers:
    func()

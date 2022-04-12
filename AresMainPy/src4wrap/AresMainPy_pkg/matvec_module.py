"""
Module matvec_module


Defined at Matvec_module.fpp lines 5-945

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def cmatvec(veff, ik, p, q, dimen):
    """
    cmatvec(veff, ik, p, q, dimen)
    
    
    Defined at Matvec_module.fpp lines 9-163
    
    Parameters
    ----------
    veff : float array
    ik : int
    p : complex array
    q : complex array
    dimen : int
    
    """
    _AresMainPy_pkg.f90wrap_cmatvec(veff=veff, ik=ik, p=p, q=q, dimen=dimen)

def nlocmatvec(ik, p, q):
    """
    nlocmatvec(ik, p, q)
    
    
    Defined at Matvec_module.fpp lines 166-227
    
    Parameters
    ----------
    ik : int
    p : complex array
    q : complex array
    
    """
    _AresMainPy_pkg.f90wrap_nlocmatvec(ik=ik, p=p, q=q)

def nlocmatvec_dg(ik, p, q):
    """
    nlocmatvec_dg(ik, p, q)
    
    
    Defined at Matvec_module.fpp lines 230-284
    
    Parameters
    ----------
    ik : int
    p : complex array
    q : complex array
    
    """
    _AresMainPy_pkg.f90wrap_nlocmatvec_dg(ik=ik, p=p, q=q)

def rmatvec(veff, p, q, dimen):
    """
    rmatvec(veff, p, q, dimen)
    
    
    Defined at Matvec_module.fpp lines 287-347
    
    Parameters
    ----------
    veff : float array
    p : float array
    q : float array
    dimen : int
    
    =====================================================
    ##OUTPUT NONLOCAL TERM TO TEST ITS CORRECTNESS
    open(unit=521,file="nl_Vrho_period")
    write(521,*)n1,n2,n3
    write(521,*)q
    close(521)
    =====================================================
    """
    _AresMainPy_pkg.f90wrap_rmatvec(veff=veff, p=p, q=q, dimen=dimen)

def nlocmatvec_r(p, q):
    """
    nlocmatvec_r(p, q)
    
    
    Defined at Matvec_module.fpp lines 350-405
    
    Parameters
    ----------
    p : float array
    q : float array
    
    """
    _AresMainPy_pkg.f90wrap_nlocmatvec_r(p=p, q=q)

def rmatvec_new(mat, p, q, dimen):
    """
    rmatvec_new(mat, p, q, dimen)
    
    
    Defined at Matvec_module.fpp lines 408-416
    
    Parameters
    ----------
    mat : float array
    p : float array
    q : float array
    dimen : int
    
    """
    _AresMainPy_pkg.f90wrap_rmatvec_new(mat=mat, p=p, q=q, dimen=dimen)

def iso_rmatvec(veff_3d, p, q, dimen):
    """
    iso_rmatvec(veff_3d, p, q, dimen)
    
    
    Defined at Matvec_module.fpp lines 420-457
    
    Parameters
    ----------
    veff_3d : float array
    p : float array
    q : float array
    dimen : int
    
    """
    _AresMainPy_pkg.f90wrap_iso_rmatvec(veff_3d=veff_3d, p=p, q=q, dimen=dimen)

def iso_grid2sphere(grid_v, grid_v1, sphere_v):
    """
    iso_grid2sphere(grid_v, grid_v1, sphere_v)
    
    
    Defined at Matvec_module.fpp lines 460-473
    
    Parameters
    ----------
    grid_v : float array
    grid_v1 : float array
    sphere_v : float array
    
    """
    _AresMainPy_pkg.f90wrap_iso_grid2sphere(grid_v=grid_v, grid_v1=grid_v1, \
        sphere_v=sphere_v)

def iso_sphere2grid(sphere_v, grid_v):
    """
    iso_sphere2grid(sphere_v, grid_v)
    
    
    Defined at Matvec_module.fpp lines 476-488
    
    Parameters
    ----------
    sphere_v : float array
    grid_v : float array
    
    """
    _AresMainPy_pkg.f90wrap_iso_sphere2grid(sphere_v=sphere_v, grid_v=grid_v)

def nlocmatvec_iso(p, q):
    """
    nlocmatvec_iso(p, q)
    
    
    Defined at Matvec_module.fpp lines 491-549
    
    Parameters
    ----------
    p : float array
    q : float array
    
    """
    _AresMainPy_pkg.f90wrap_nlocmatvec_iso(p=p, q=q)

def nlocmatvec_iso_dg(p, q):
    """
    nlocmatvec_iso_dg(p, q)
    
    
    Defined at Matvec_module.fpp lines 552-609
    
    Parameters
    ----------
    p : float array
    q : float array
    
    """
    _AresMainPy_pkg.f90wrap_nlocmatvec_iso_dg(p=p, q=q)

def cmatvec_band(veff, ik, p, q, dimen):
    """
    cmatvec_band(veff, ik, p, q, dimen)
    
    
    Defined at Matvec_module.fpp lines 613-699
    
    Parameters
    ----------
    veff : float array
    ik : int
    p : complex array
    q : complex array
    dimen : int
    
    """
    _AresMainPy_pkg.f90wrap_cmatvec_band(veff=veff, ik=ik, p=p, q=q, dimen=dimen)

def nlocmatvec_band(ik, p, q):
    """
    nlocmatvec_band(ik, p, q)
    
    
    Defined at Matvec_module.fpp lines 702-760
    
    Parameters
    ----------
    ik : int
    p : complex array
    q : complex array
    
    """
    _AresMainPy_pkg.f90wrap_nlocmatvec_band(ik=ik, p=p, q=q)

def rmatvec_gamma(veff, ik, p, q, dimen):
    """
    rmatvec_gamma(veff, ik, p, q, dimen)
    
    
    Defined at Matvec_module.fpp lines 763-881
    
    Parameters
    ----------
    veff : float array
    ik : int
    p : float array
    q : float array
    dimen : int
    
    """
    _AresMainPy_pkg.f90wrap_rmatvec_gamma(veff=veff, ik=ik, p=p, q=q, dimen=dimen)

def nlocmatvec_gamma(ik, p, q):
    """
    nlocmatvec_gamma(ik, p, q)
    
    
    Defined at Matvec_module.fpp lines 884-945
    
    Parameters
    ----------
    ik : int
    p : float array
    q : float array
    
    """
    _AresMainPy_pkg.f90wrap_nlocmatvec_gamma(ik=ik, p=p, q=q)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "matvec_module".')

for func in _dt_array_initialisers:
    func()

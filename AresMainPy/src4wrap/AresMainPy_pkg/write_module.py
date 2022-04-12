"""
Module write_module


Defined at Write_module.fpp lines 5-203

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def write_poscar(filename, lat_mat, eleid, pos, atom_symbol=None, fixpos=None):
    """
    write_poscar(filename, lat_mat, eleid, pos[, atom_symbol, fixpos])
    
    
    Defined at Write_module.fpp lines 22-62
    
    Parameters
    ----------
    filename : str
    lat_mat : float array
    eleid : int array
    pos : float array
    atom_symbol : str array
    fixpos : bool array
    
    """
    _AresMainPy_pkg.f90wrap_write_poscar(filename=filename, lat_mat=lat_mat, \
        eleid=eleid, pos=pos, atom_symbol=atom_symbol, fixpos=fixpos)

def write_cif(filename, lat_para, nati, pos, atom_symbol=None):
    """
    write_cif(filename, lat_para, nati, pos[, atom_symbol])
    
    
    Defined at Write_module.fpp lines 66-131
    
    Parameters
    ----------
    filename : str
    lat_para : float array
    nati : int array
    pos : float array
    atom_symbol : str array
    
    """
    _AresMainPy_pkg.f90wrap_write_cif(filename=filename, lat_para=lat_para, \
        nati=nati, pos=pos, atom_symbol=atom_symbol)

def write3dat(filename, rho):
    """
    write3dat(filename, rho)
    
    
    Defined at Write_module.fpp lines 135-164
    
    Parameters
    ----------
    filename : str
    rho : float array
    
    -----------------------------------------------------------------------
    """
    _AresMainPy_pkg.f90wrap_write3dat(filename=filename, rho=rho)

def write3d(filename, rho):
    """
    write3d(filename, rho)
    
    
    Defined at Write_module.fpp lines 168-192
    
    Parameters
    ----------
    filename : str
    rho : float array
    
    -----------------------------------------------------------------------
    """
    _AresMainPy_pkg.f90wrap_write3d(filename=filename, rho=rho)

def writedensity(filename, rho):
    """
    writedensity(filename, rho)
    
    
    Defined at Write_module.fpp lines 198-201
    
    Parameters
    ----------
    filename : str
    rho : float array
    
    """
    _AresMainPy_pkg.f90wrap_writedensity(filename=filename, rho=rho)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "write_module".')

for func in _dt_array_initialisers:
    func()

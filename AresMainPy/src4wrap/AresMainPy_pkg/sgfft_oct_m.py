"""
Module sgfft_oct_m


Defined at sgfft.fpp lines 26-4425

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def fourier_dim(n):
    """
    n_next = fourier_dim(n)
    
    
    Defined at sgfft.fpp lines 51-81
    
    Parameters
    ----------
    n : int
    
    Returns
    -------
    n_next : int
    
    """
    n_next = _AresMainPy_pkg.f90wrap_fourier_dim(n=n)
    return n_next

def fft(n1, n2, n3, nd1, nd2, nd3, z, isign, inzee):
    """
    fft(n1, n2, n3, nd1, nd2, nd3, z, isign, inzee)
    
    
    Defined at sgfft.fpp lines 161-385
    
    Parameters
    ----------
    n1 : int
    n2 : int
    n3 : int
    nd1 : int
    nd2 : int
    nd3 : int
    z : float array
    isign : int
    inzee : int
    
    """
    _AresMainPy_pkg.f90wrap_fft(n1=n1, n2=n2, n3=n3, nd1=nd1, nd2=nd2, nd3=nd3, z=z, \
        isign=isign, inzee=inzee)

def convolxc_off(n1, n2, n3, nd1, nd2, nd3, md1, md2, md3, nproc, iproc, pot, \
    zf, scal, comm):
    """
    convolxc_off(n1, n2, n3, nd1, nd2, nd3, md1, md2, md3, nproc, iproc, pot, zf, \
        scal, comm)
    
    
    Defined at sgfft.fpp lines 3504-3796
    
    Parameters
    ----------
    n1 : int
    n2 : int
    n3 : int
    nd1 : int
    nd2 : int
    nd3 : int
    md1 : int
    md2 : int
    md3 : int
    nproc : int
    iproc : int
    pot : float array
    zf : float array
    scal : float
    comm : int
    
    """
    _AresMainPy_pkg.f90wrap_convolxc_off(n1=n1, n2=n2, n3=n3, nd1=nd1, nd2=nd2, \
        nd3=nd3, md1=md1, md2=md2, md3=md3, nproc=nproc, iproc=iproc, pot=pot, \
        zf=zf, scal=scal, comm=comm)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "sgfft_oct_m".')

for func in _dt_array_initialisers:
    func()

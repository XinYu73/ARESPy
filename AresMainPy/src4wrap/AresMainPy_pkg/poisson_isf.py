"""
Module poisson_isf


Defined at poisson_isf.fpp lines 5-1430

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def karray_set():
    """
    karray_set()
    
    
    Defined at poisson_isf.fpp lines 72-94
    
    
    """
    _AresMainPy_pkg.f90wrap_karray_set()

def poisson_isf_method(rho, vh):
    """
    poisson_isf_method(rho, vh)
    
    
    Defined at poisson_isf.fpp lines 97-110
    
    Parameters
    ----------
    rho : float array
    vh : float array
    
    """
    _AresMainPy_pkg.f90wrap_poisson_isf_method(rho=rho, vh=vh)

def build_kernel(n01, n02, n03, nfft1, nfft2, nfft3, hgrid, itype_scf, \
    karrayout):
    """
    build_kernel(n01, n02, n03, nfft1, nfft2, nfft3, hgrid, itype_scf, karrayout)
    
    
    Defined at poisson_isf.fpp lines 113-249
    
    Parameters
    ----------
    n01 : int
    n02 : int
    n03 : int
    nfft1 : int
    nfft2 : int
    nfft3 : int
    hgrid : float
    itype_scf : int
    karrayout : float array
    
    """
    _AresMainPy_pkg.f90wrap_build_kernel(n01=n01, n02=n02, n03=n03, nfft1=nfft1, \
        nfft2=nfft2, nfft3=nfft3, hgrid=hgrid, itype_scf=itype_scf, \
        karrayout=karrayout)

def scaling_function(itype, nd, a, x):
    """
    nrange = scaling_function(itype, nd, a, x)
    
    
    Defined at poisson_isf.fpp lines 252-290
    
    Parameters
    ----------
    itype : int
    nd : int
    a : float array
    x : float array
    
    Returns
    -------
    nrange : int
    
    """
    nrange = _AresMainPy_pkg.f90wrap_scaling_function(itype=itype, nd=nd, a=a, x=x)
    return nrange

def back_trans_8(nd, nt, x, y):
    """
    back_trans_8(nd, nt, x, y)
    
    
    Defined at poisson_isf.fpp lines 293-383
    
    Parameters
    ----------
    nd : int
    nt : int
    x : float array
    y : float array
    
    """
    _AresMainPy_pkg.f90wrap_back_trans_8(nd=nd, nt=nt, x=x, y=y)

def gequad(n_gauss, p_gauss, w_gauss):
    """
    ur_gauss, dr_gauss, acc_gauss = gequad(n_gauss, p_gauss, w_gauss)
    
    
    Defined at poisson_isf.fpp lines 386-582
    
    Parameters
    ----------
    n_gauss : int
    p_gauss : float array
    w_gauss : float array
    
    Returns
    -------
    ur_gauss : float
    dr_gauss : float
    acc_gauss : float
    
    """
    ur_gauss, dr_gauss, acc_gauss = _AresMainPy_pkg.f90wrap_gequad(n_gauss=n_gauss, \
        p_gauss=p_gauss, w_gauss=w_gauss)
    return ur_gauss, dr_gauss, acc_gauss

def calculate_dimensions(n01, n02, n03):
    """
    nfft1, nfft2, nfft3 = calculate_dimensions(n01, n02, n03)
    
    
    Defined at poisson_isf.fpp lines 585-619
    
    Parameters
    ----------
    n01 : int
    n02 : int
    n03 : int
    
    Returns
    -------
    nfft1 : int
    nfft2 : int
    nfft3 : int
    
    """
    nfft1, nfft2, nfft3 = _AresMainPy_pkg.f90wrap_calculate_dimensions(n01=n01, \
        n02=n02, n03=n03)
    return nfft1, nfft2, nfft3

def psolver_kernel(n01, n02, n03, nfft1, nfft2, nfft3, hgrid, karray, rhopot):
    """
    psolver_kernel(n01, n02, n03, nfft1, nfft2, nfft3, hgrid, karray, rhopot)
    
    
    Defined at poisson_isf.fpp lines 622-663
    
    Parameters
    ----------
    n01 : int
    n02 : int
    n03 : int
    nfft1 : int
    nfft2 : int
    nfft3 : int
    hgrid : float
    karray : float array
    rhopot : float array
    
    """
    _AresMainPy_pkg.f90wrap_psolver_kernel(n01=n01, n02=n02, n03=n03, nfft1=nfft1, \
        nfft2=nfft2, nfft3=nfft3, hgrid=hgrid, karray=karray, rhopot=rhopot)

def zarray_in(n01, n02, n03, nd1, nd2, nd3, density, zarray):
    """
    zarray_in(n01, n02, n03, nd1, nd2, nd3, density, zarray)
    
    
    Defined at poisson_isf.fpp lines 666-700
    
    Parameters
    ----------
    n01 : int
    n02 : int
    n03 : int
    nd1 : int
    nd2 : int
    nd3 : int
    density : float array
    zarray : float array
    
    """
    _AresMainPy_pkg.f90wrap_zarray_in(n01=n01, n02=n02, n03=n03, nd1=nd1, nd2=nd2, \
        nd3=nd3, density=density, zarray=zarray)

def zarray_out(n01, n02, n03, nd1, nd2, nd3, rhopot, zarray, factor):
    """
    zarray_out(n01, n02, n03, nd1, nd2, nd3, rhopot, zarray, factor)
    
    
    Defined at poisson_isf.fpp lines 703-716
    
    Parameters
    ----------
    n01 : int
    n02 : int
    n03 : int
    nd1 : int
    nd2 : int
    nd3 : int
    rhopot : float array
    zarray : float array
    factor : float
    
    """
    _AresMainPy_pkg.f90wrap_zarray_out(n01=n01, n02=n02, n03=n03, nd1=nd1, nd2=nd2, \
        nd3=nd3, rhopot=rhopot, zarray=zarray, factor=factor)

def scf_recursion(itype, n_iter, n_range, kernel_scf, kern_1_scf):
    """
    scf_recursion(itype, n_iter, n_range, kernel_scf, kern_1_scf)
    
    
    Defined at poisson_isf.fpp lines 729-745
    
    Parameters
    ----------
    itype : int
    n_iter : int
    n_range : int
    kernel_scf : float array
    kern_1_scf : float array
    
    """
    _AresMainPy_pkg.f90wrap_scf_recursion(itype=itype, n_iter=n_iter, \
        n_range=n_range, kernel_scf=kernel_scf, kern_1_scf=kern_1_scf)

def scf_recursion_8(n_iter, n_range, kernel_scf, kern_1_scf):
    """
    scf_recursion_8(n_iter, n_range, kernel_scf, kern_1_scf)
    
    
    Defined at poisson_isf.fpp lines 758-852
    
    Parameters
    ----------
    n_iter : int
    n_range : int
    kernel_scf : float array
    kern_1_scf : float array
    
    """
    _AresMainPy_pkg.f90wrap_scf_recursion_8(n_iter=n_iter, n_range=n_range, \
        kernel_scf=kernel_scf, kern_1_scf=kern_1_scf)

def karrayhalf_in(n01, n02, n03, n1k, n2k, n3k, nfft1, nfft2, nfft3, nd1, nd2, \
    nd3, kernel, karrayhalf):
    """
    karrayhalf_in(n01, n02, n03, n1k, n2k, n3k, nfft1, nfft2, nfft3, nd1, nd2, nd3, \
        kernel, karrayhalf)
    
    
    Defined at poisson_isf.fpp lines 864-905
    
    Parameters
    ----------
    n01 : int
    n02 : int
    n03 : int
    n1k : int
    n2k : int
    n3k : int
    nfft1 : int
    nfft2 : int
    nfft3 : int
    nd1 : int
    nd2 : int
    nd3 : int
    kernel : float array
    karrayhalf : float array
    
    """
    _AresMainPy_pkg.f90wrap_karrayhalf_in(n01=n01, n02=n02, n03=n03, n1k=n1k, \
        n2k=n2k, n3k=n3k, nfft1=nfft1, nfft2=nfft2, nfft3=nfft3, nd1=nd1, nd2=nd2, \
        nd3=nd3, kernel=kernel, karrayhalf=karrayhalf)

def kernel_recon(n1k, n2k, n3k, nfft1, nfft2, nfft3, nd1, nd2, nd3, zarray, \
    karray):
    """
    kernel_recon(n1k, n2k, n3k, nfft1, nfft2, nfft3, nd1, nd2, nd3, zarray, karray)
    
    
    Defined at poisson_isf.fpp lines 919-969
    
    Parameters
    ----------
    n1k : int
    n2k : int
    n3k : int
    nfft1 : int
    nfft2 : int
    nfft3 : int
    nd1 : int
    nd2 : int
    nd3 : int
    zarray : float array
    karray : float array
    
    """
    _AresMainPy_pkg.f90wrap_kernel_recon(n1k=n1k, n2k=n2k, n3k=n3k, nfft1=nfft1, \
        nfft2=nfft2, nfft3=nfft3, nd1=nd1, nd2=nd2, nd3=nd3, zarray=zarray, \
        karray=karray)

def norm_ind(nd1, nd2, nd3, i1, i2, i3, ind):
    """
    norm_ind(nd1, nd2, nd3, i1, i2, i3, ind)
    
    
    Defined at poisson_isf.fpp lines 981-1001
    
    Parameters
    ----------
    nd1 : int
    nd2 : int
    nd3 : int
    i1 : int
    i2 : int
    i3 : int
    ind : int
    
    """
    _AresMainPy_pkg.f90wrap_norm_ind(nd1=nd1, nd2=nd2, nd3=nd3, i1=i1, i2=i2, i3=i3, \
        ind=ind)

def symm_ind(nd1, nd2, nd3, i1, i2, i3, ind):
    """
    symm_ind(nd1, nd2, nd3, i1, i2, i3, ind)
    
    
    Defined at poisson_isf.fpp lines 1014-1033
    
    Parameters
    ----------
    nd1 : int
    nd2 : int
    nd3 : int
    i1 : int
    i2 : int
    i3 : int
    ind : int
    
    """
    _AresMainPy_pkg.f90wrap_symm_ind(nd1=nd1, nd2=nd2, nd3=nd3, i1=i1, i2=i2, i3=i3, \
        ind=ind)

def kernel_application(n1xy, n2xy, n3xy, nd1h, nd2, nd3, nfft1, nfft2, nfft3, \
    zarray, karray, inzee):
    """
    kernel_application(n1xy, n2xy, n3xy, nd1h, nd2, nd3, nfft1, nfft2, nfft3, \
        zarray, karray, inzee)
    
    
    Defined at poisson_isf.fpp lines 1037-1429
    
    Parameters
    ----------
    n1xy : int
    n2xy : int
    n3xy : int
    nd1h : int
    nd2 : int
    nd3 : int
    nfft1 : int
    nfft2 : int
    nfft3 : int
    zarray : float array
    karray : float array
    inzee : int
    
    --------------------------------------------
    --- Starting reconstruction half -> full ---
    --------------------------------------------
    -------------Case i3 = 1
    """
    _AresMainPy_pkg.f90wrap_kernel_application(n1xy=n1xy, n2xy=n2xy, n3xy=n3xy, \
        nd1h=nd1h, nd2=nd2, nd3=nd3, nfft1=nfft1, nfft2=nfft2, nfft3=nfft3, \
        zarray=zarray, karray=karray, inzee=inzee)

def get_array_karray():
    """
    Element karray ftype=real(dp) pytype=float
    
    
    Defined at poisson_isf.fpp line 39
    
    """
    global karray
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_poisson_isf__array__karray(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        karray = _arrays[array_handle]
    else:
        karray = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_poisson_isf__array__karray)
        _arrays[array_handle] = karray
    return karray

def set_array_karray(karray):
    karray[...] = karray


_array_initialisers = [get_array_karray]
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "poisson_isf".')

for func in _dt_array_initialisers:
    func()

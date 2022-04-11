"""
Module chebyshev_module


Defined at Chebyshev_fliter.fpp lines 10-747

"""
from __future__ import print_function, absolute_import, division
import _arespy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def buildsubspace(nps, nev, veff, eig):
    """
    buildsubspace(nps, nev, veff, eig)
    
    
    Defined at Chebyshev_fliter.fpp lines 22-119
    
    Parameters
    ----------
    nps : int
    nev : int
    veff : float array
    eig : Eigen_Type
    
    """
    _arespy_pkg.f90wrap_buildsubspace(nps=nps, nev=nev, veff=veff, eig=eig._handle)

def real_pseudosubspace(nps, nev, initx):
    """
    real_pseudosubspace(nps, nev, initx)
    
    
    Defined at Chebyshev_fliter.fpp lines 122-197
    
    Parameters
    ----------
    nps : int
    nev : int
    initx : float array
    
    """
    _arespy_pkg.f90wrap_real_pseudosubspace(nps=nps, nev=nev, initx=initx)

def real_first_rrstep(nps, nev, veff, x, d):
    """
    real_first_rrstep(nps, nev, veff, x, d)
    
    
    Defined at Chebyshev_fliter.fpp lines 200-276
    
    Parameters
    ----------
    nps : int
    nev : int
    veff : float array
    x : float array
    d : float array
    
    -------------------
    rotation
    write(*, *) 'line302'
    """
    _arespy_pkg.f90wrap_real_first_rrstep(nps=nps, nev=nev, veff=veff, x=x, d=d)

def init_uplow_real(nps, k, veff, v):
    """
    a, b, al = init_uplow_real(nps, k, veff, v)
    
    
    Defined at Chebyshev_fliter.fpp lines 279-335
    
    Parameters
    ----------
    nps : int
    k : int
    veff : float array
    v : float array
    
    Returns
    -------
    a : float
    b : float
    al : float
    
    """
    a, b, al = _arespy_pkg.f90wrap_init_uplow_real(nps=nps, k=k, veff=veff, v=v)
    return a, b, al

def real_first_filter(nps, nst, veff, x, eval):
    """
    real_first_filter(nps, nst, veff, x, eval)
    
    
    Defined at Chebyshev_fliter.fpp lines 338-419
    
    Parameters
    ----------
    nps : int
    nst : int
    veff : float array
    x : float array
    eval : float array
    
    -----------------
    """
    _arespy_pkg.f90wrap_real_first_filter(nps=nps, nst=nst, veff=veff, x=x, \
        eval=eval)

def rayleigh_quotient_real(nps, nst, veff, x, xhx):
    """
    rayleigh_quotient_real(nps, nst, veff, x, xhx)
    
    
    Defined at Chebyshev_fliter.fpp lines 427-446
    
    Parameters
    ----------
    nps : int
    nst : int
    veff : float array
    x : float array
    xhx : float array
    
    """
    _arespy_pkg.f90wrap_rayleigh_quotient_real(nps=nps, nst=nst, veff=veff, x=x, \
        xhx=xhx)

def cal_hx_real(nps, nst, veff, v, hv):
    """
    cal_hx_real(nps, nst, veff, v, hv)
    
    
    Defined at Chebyshev_fliter.fpp lines 449-463
    
    Parameters
    ----------
    nps : int
    nst : int
    veff : float array
    v : float array
    hv : float array
    
    """
    _arespy_pkg.f90wrap_cal_hx_real(nps=nps, nst=nst, veff=veff, v=v, hv=hv)

def estupb_real(nps, k, veff, vec):
    """
    b = estupb_real(nps, k, veff, vec)
    
    
    Defined at Chebyshev_fliter.fpp lines 466-521
    
    Parameters
    ----------
    nps : int
    k : int
    veff : float array
    vec : float array
    
    Returns
    -------
    b : float
    
    """
    b = _arespy_pkg.f90wrap_estupb_real(nps=nps, k=k, veff=veff, vec=vec)
    return b

def chebyshev_filter_real(nps, nst, veff, x, m, a, b):
    """
    chebyshev_filter_real(nps, nst, veff, x, m, a, b)
    
    
    Defined at Chebyshev_fliter.fpp lines 524-558
    
    Parameters
    ----------
    nps : int
    nst : int
    veff : float array
    x : float array
    m : int
    a : float
    b : float
    
    """
    _arespy_pkg.f90wrap_chebyshev_filter_real(nps=nps, nst=nst, veff=veff, x=x, m=m, \
        a=a, b=b)

def chebyshev_filter_scaled_real(nps, nst, veff, x, m, a, b, al):
    """
    chebyshev_filter_scaled_real(nps, nst, veff, x, m, a, b, al)
    
    
    Defined at Chebyshev_fliter.fpp lines 561-599
    
    Parameters
    ----------
    nps : int
    nst : int
    veff : float array
    x : float array
    m : int
    a : float
    b : float
    al : float
    
    """
    _arespy_pkg.f90wrap_chebyshev_filter_scaled_real(nps=nps, nst=nst, veff=veff, \
        x=x, m=m, a=a, b=b, al=al)

def grayleigh_ritz_real(nps, nev, veff, x, d):
    """
    grayleigh_ritz_real(nps, nev, veff, x, d)
    
    
    Defined at Chebyshev_fliter.fpp lines 602-651
    
    Parameters
    ----------
    nps : int
    nev : int
    veff : float array
    x : float array
    d : float array
    
    """
    _arespy_pkg.f90wrap_grayleigh_ritz_real(nps=nps, nev=nev, veff=veff, x=x, d=d)

def rayleigh_ritz_real(nps, sn, veff, x, d):
    """
    rayleigh_ritz_real(nps, sn, veff, x, d)
    
    
    Defined at Chebyshev_fliter.fpp lines 654-700
    
    Parameters
    ----------
    nps : int
    sn : int
    veff : float array
    x : float array
    d : float array
    
    """
    _arespy_pkg.f90wrap_rayleigh_ritz_real(nps=nps, sn=sn, veff=veff, x=x, d=d)

def cheby_filtering_grrr(nps, nev, veff, x, d):
    """
    cheby_filtering_grrr(nps, nev, veff, x, d)
    
    
    Defined at Chebyshev_fliter.fpp lines 703-746
    
    Parameters
    ----------
    nps : int
    nev : int
    veff : float array
    x : float array
    d : float array
    
    """
    _arespy_pkg.f90wrap_cheby_filtering_grrr(nps=nps, nev=nev, veff=veff, x=x, d=d)

def get_larged():
    """
    Element larged ftype=real(dp) pytype=float
    
    
    Defined at Chebyshev_fliter.fpp line 16
    
    """
    return _arespy_pkg.f90wrap_chebyshev_module__get__larged()

LARGED = get_larged()

def get_array_ad():
    """
    Element ad ftype=real(dp) pytype=float
    
    
    Defined at Chebyshev_fliter.fpp line 17
    
    """
    global ad
    array_ndim, array_type, array_shape, array_handle = \
        _arespy_pkg.f90wrap_chebyshev_module__array__ad(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        ad = _arrays[array_handle]
    else:
        ad = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_chebyshev_module__array__ad)
        _arrays[array_handle] = ad
    return ad

def set_array_ad(ad):
    ad[...] = ad


_array_initialisers = [get_array_ad]
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "chebyshev_module".')

for func in _dt_array_initialisers:
    func()

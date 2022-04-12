"""
Module mathsplines


Defined at MathSplines.fpp lines 5-836

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def spline_cubic_set(n, t, y, ibcbeg, ybcbeg, ibcend, ybcend, ypp):
    """
    spline_cubic_set(n, t, y, ibcbeg, ybcbeg, ibcend, ybcend, ypp)
    
    
    Defined at MathSplines.fpp lines 178-387
    
    Parameters
    ----------
    n : int
    t : float array
    y : float array
    ibcbeg : int
    ybcbeg : float
    ibcend : int
    ybcend : float
    ypp : float array
    
    """
    _AresMainPy_pkg.f90wrap_spline_cubic_set(n=n, t=t, y=y, ibcbeg=ibcbeg, \
        ybcbeg=ybcbeg, ibcend=ibcend, ybcend=ybcend, ypp=ypp)

def spline_cubic_val(n, t, y, ypp, tval):
    """
    yval, ypval, yppval = spline_cubic_val(n, t, y, ypp, tval)
    
    
    Defined at MathSplines.fpp lines 389-483
    
    Parameters
    ----------
    n : int
    t : float array
    y : float array
    ypp : float array
    tval : float
    
    Returns
    -------
    yval : float
    ypval : float
    yppval : float
    
    """
    yval, ypval, yppval = _AresMainPy_pkg.f90wrap_spline_cubic_val(n=n, t=t, y=y, \
        ypp=ypp, tval=tval)
    return yval, ypval, yppval

def rvec_bracket(n, x, xval):
    """
    left, right = rvec_bracket(n, x, xval)
    
    
    Defined at MathSplines.fpp lines 485-581
    
    Parameters
    ----------
    n : int
    x : float array
    xval : float
    
    Returns
    -------
    left : int
    right : int
    
    """
    left, right = _AresMainPy_pkg.f90wrap_rvec_bracket(n=n, x=x, xval=xval)
    return left, right

def s3_fs(a1, a2, a3, n, b, x):
    """
    s3_fs(a1, a2, a3, n, b, x)
    
    
    Defined at MathSplines.fpp lines 583-650
    
    Parameters
    ----------
    a1 : float array
    a2 : float array
    a3 : float array
    n : int
    b : float array
    x : float array
    
    """
    _AresMainPy_pkg.f90wrap_s3_fs(a1=a1, a2=a2, a3=a3, n=n, b=b, x=x)

def polynom(m, np, xa, ya, c, x):
    """
    polynom = polynom(m, np, xa, ya, c, x)
    
    
    Defined at MathSplines.fpp lines 658-836
    
    Parameters
    ----------
    m : int
    np : int
    xa : float array
    ya : float array
    c : float array
    x : float
    
    Returns
    -------
    polynom : float
    
    """
    polynom = _AresMainPy_pkg.f90wrap_polynom(m=m, np=np, xa=xa, ya=ya, c=c, x=x)
    return polynom


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "mathsplines".')

for func in _dt_array_initialisers:
    func()

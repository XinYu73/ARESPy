"""
Module getvlocalpseudopotential


Defined at IonLocalPotentialAssignment.fpp lines 5-940

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def calvlpp():
    """
    calvlpp()
    
    
    Defined at IonLocalPotentialAssignment.fpp lines 26-107
    
    
    ==========================================================
    """
    _AresMainPy_pkg.f90wrap_calvlpp()

def ionpotentialassignment(ity, zion, poscar, temp):
    """
    ionpotentialassignment(ity, zion, poscar, temp)
    
    
    Defined at IonLocalPotentialAssignment.fpp lines 110-186
    
    Parameters
    ----------
    ity : int
    zion : float
    poscar : float array
    temp : float array
    
    ========================================================================
    """
    _AresMainPy_pkg.f90wrap_ionpotentialassignment(ity=ity, zion=zion, \
        poscar=poscar, temp=temp)

def sphbess(l, x):
    """
    sphbess = sphbess(l, x)
    
    
    Defined at IonLocalPotentialAssignment.fpp lines 191-230
    
    Parameters
    ----------
    l : int
    x : float
    
    Returns
    -------
    sphbess : float
    
    """
    sphbess = _AresMainPy_pkg.f90wrap_sphbess(l=l, x=x)
    return sphbess

def fourbess_gr(g, fg, r, fr):
    """
    fourbess_gr(g, fg, r, fr)
    
    
    Defined at IonLocalPotentialAssignment.fpp lines 234-257
    
    Parameters
    ----------
    g : float array
    fg : float array
    r : float array
    fr : float array
    
    """
    _AresMainPy_pkg.f90wrap_fourbess_gr(g=g, fg=fg, r=r, fr=fr)

def dir2car_single(cry_coo, ort_coo, lat):
    """
    dir2car_single(cry_coo, ort_coo, lat)
    
    
    Defined at IonLocalPotentialAssignment.fpp lines 598-605
    
    Parameters
    ----------
    cry_coo : float array
    ort_coo : float array
    lat : float array
    
    """
    _AresMainPy_pkg.f90wrap_dir2car_single(cry_coo=cry_coo, ort_coo=ort_coo, \
        lat=lat)

def ionpotentialassignment_dg(ity, zion, poscar, temp):
    """
    ionpotentialassignment_dg(ity, zion, poscar, temp)
    
    
    Defined at IonLocalPotentialAssignment.fpp lines 800-866
    
    Parameters
    ----------
    ity : int
    zion : float
    poscar : float array
    temp : float array
    
    ========================================================================
    ==========================================================xlt test
    """
    _AresMainPy_pkg.f90wrap_ionpotentialassignment_dg(ity=ity, zion=zion, \
        poscar=poscar, temp=temp)

def set_vloc_dg(ity, xyz):
    """
    vloc = set_vloc_dg(ity, xyz)
    
    
    Defined at IonLocalPotentialAssignment.fpp lines 869-915
    
    Parameters
    ----------
    ity : int
    xyz : float array
    
    Returns
    -------
    vloc : float
    
    """
    vloc = _AresMainPy_pkg.f90wrap_set_vloc_dg(ity=ity, xyz=xyz)
    return vloc

def set_vcomp(n_inp, zion, rgauss, r_inp, v_inp, vcomp):
    """
    set_vcomp(n_inp, zion, rgauss, r_inp, v_inp, vcomp)
    
    
    Defined at IonLocalPotentialAssignment.fpp lines 918-938
    
    Parameters
    ----------
    n_inp : int
    zion : int
    rgauss : float
    r_inp : float array
    v_inp : float array
    vcomp : float array
    
    """
    _AresMainPy_pkg.f90wrap_set_vcomp(n_inp=n_inp, zion=zion, rgauss=rgauss, \
        r_inp=r_inp, v_inp=v_inp, vcomp=vcomp)

def get_radiusn():
    """
    Element radiusn ftype=integer(i4b) pytype=int
    
    
    Defined at IonLocalPotentialAssignment.fpp line 15
    
    """
    return _AresMainPy_pkg.f90wrap_getvlocalpseudopotential__get__radiusn()

def set_radiusn(radiusn):
    _AresMainPy_pkg.f90wrap_getvlocalpseudopotential__set__radiusn(radiusn)

def get_array_fr():
    """
    Element fr ftype=real(dp) pytype=float
    
    
    Defined at IonLocalPotentialAssignment.fpp line 16
    
    """
    global fr
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_getvlocalpseudopotential__array__fr(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        fr = _arrays[array_handle]
    else:
        fr = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_getvlocalpseudopotential__array__fr)
        _arrays[array_handle] = fr
    return fr

def set_array_fr(fr):
    fr[...] = fr

def get_array_wijk():
    """
    Element wijk ftype=real(dp) pytype=float
    
    
    Defined at IonLocalPotentialAssignment.fpp line 19
    
    """
    global wijk
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_getvlocalpseudopotential__array__wijk(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        wijk = _arrays[array_handle]
    else:
        wijk = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_getvlocalpseudopotential__array__wijk)
        _arrays[array_handle] = wijk
    return wijk

def set_array_wijk(wijk):
    wijk[...] = wijk

def get_array_drijk():
    """
    Element drijk ftype=real(dp) pytype=float
    
    
    Defined at IonLocalPotentialAssignment.fpp line 19
    
    """
    global drijk
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_getvlocalpseudopotential__array__drijk(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        drijk = _arrays[array_handle]
    else:
        drijk = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_getvlocalpseudopotential__array__drijk)
        _arrays[array_handle] = drijk
    return drijk

def set_array_drijk(drijk):
    drijk[...] = drijk

def get_ilft():
    """
    Element ilft ftype=integer(i4b) pytype=int
    
    
    Defined at IonLocalPotentialAssignment.fpp line 20
    
    """
    return _AresMainPy_pkg.f90wrap_getvlocalpseudopotential__get__ilft()

def set_ilft(ilft):
    _AresMainPy_pkg.f90wrap_getvlocalpseudopotential__set__ilft(ilft)

def get_irit():
    """
    Element irit ftype=integer(i4b) pytype=int
    
    
    Defined at IonLocalPotentialAssignment.fpp line 20
    
    """
    return _AresMainPy_pkg.f90wrap_getvlocalpseudopotential__get__irit()

def set_irit(irit):
    _AresMainPy_pkg.f90wrap_getvlocalpseudopotential__set__irit(irit)

def get_numdg():
    """
    Element numdg ftype=integer(i4b) pytype=int
    
    
    Defined at IonLocalPotentialAssignment.fpp line 20
    
    """
    return _AresMainPy_pkg.f90wrap_getvlocalpseudopotential__get__numdg()

def set_numdg(numdg):
    _AresMainPy_pkg.f90wrap_getvlocalpseudopotential__set__numdg(numdg)


_array_initialisers = [get_array_fr, get_array_wijk, get_array_drijk]
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "getvlocalpseudopotential".')

for func in _dt_array_initialisers:
    func()

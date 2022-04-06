"""
Module arpack_module


Defined at Arpack_module.fpp lines 5-192

"""
from __future__ import print_function, absolute_import, division
import _arespy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def real_diagh_arpack(n, veff, nev, evec, eval, resid_restart, nec, info, \
    maxmvs, tol):
    """
    real_diagh_arpack(n, veff, nev, evec, eval, resid_restart, nec, info, maxmvs, \
        tol)
    
    
    Defined at Arpack_module.fpp lines 39-191
    
    Parameters
    ----------
    n : int
    veff : float array
    nev : int
    evec : float array
    eval : float array
    resid_restart : float array
    nec : int
    info : int
    maxmvs : int
    tol : float
    
    ------
    """
    _arespy_pkg.f90wrap_real_diagh_arpack(n=n, veff=veff, nev=nev, evec=evec, \
        eval=eval, resid_restart=resid_restart, nec=nec, info=info, maxmvs=maxmvs, \
        tol=tol)

def get_maxn():
    """
    Element maxn ftype=integer(i4b) pytype=int
    
    
    Defined at Arpack_module.fpp line 13
    
    """
    return _arespy_pkg.f90wrap_arpack_module__get__maxn()

maxn = get_maxn()

def get_maxnev():
    """
    Element maxnev ftype=integer(i4b) pytype=int
    
    
    Defined at Arpack_module.fpp line 13
    
    """
    return _arespy_pkg.f90wrap_arpack_module__get__maxnev()

maxnev = get_maxnev()

def get_maxncv():
    """
    Element maxncv ftype=integer(i4b) pytype=int
    
    
    Defined at Arpack_module.fpp line 13
    
    """
    return _arespy_pkg.f90wrap_arpack_module__get__maxncv()

maxncv = get_maxncv()

def get_logfil():
    """
    Element logfil ftype=integer  pytype=int
    
    
    Defined at Arpack_module.fpp line 24
    
    """
    return _arespy_pkg.f90wrap_arpack_module__get__logfil()

def set_logfil(logfil):
    _arespy_pkg.f90wrap_arpack_module__set__logfil(logfil)

def get_ndigit():
    """
    Element ndigit ftype=integer  pytype=int
    
    
    Defined at Arpack_module.fpp line 24
    
    """
    return _arespy_pkg.f90wrap_arpack_module__get__ndigit()

def set_ndigit(ndigit):
    _arespy_pkg.f90wrap_arpack_module__set__ndigit(ndigit)

def get_mgetv0():
    """
    Element mgetv0 ftype=integer  pytype=int
    
    
    Defined at Arpack_module.fpp line 24
    
    """
    return _arespy_pkg.f90wrap_arpack_module__get__mgetv0()

def set_mgetv0(mgetv0):
    _arespy_pkg.f90wrap_arpack_module__set__mgetv0(mgetv0)

def get_msaupd():
    """
    Element msaupd ftype=integer  pytype=int
    
    
    Defined at Arpack_module.fpp line 24
    
    """
    return _arespy_pkg.f90wrap_arpack_module__get__msaupd()

def set_msaupd(msaupd):
    _arespy_pkg.f90wrap_arpack_module__set__msaupd(msaupd)

def get_msaup2():
    """
    Element msaup2 ftype=integer  pytype=int
    
    
    Defined at Arpack_module.fpp line 24
    
    """
    return _arespy_pkg.f90wrap_arpack_module__get__msaup2()

def set_msaup2(msaup2):
    _arespy_pkg.f90wrap_arpack_module__set__msaup2(msaup2)

def get_msaitr():
    """
    Element msaitr ftype=integer  pytype=int
    
    
    Defined at Arpack_module.fpp line 24
    
    """
    return _arespy_pkg.f90wrap_arpack_module__get__msaitr()

def set_msaitr(msaitr):
    _arespy_pkg.f90wrap_arpack_module__set__msaitr(msaitr)

def get_mseigt():
    """
    Element mseigt ftype=integer  pytype=int
    
    
    Defined at Arpack_module.fpp line 24
    
    """
    return _arespy_pkg.f90wrap_arpack_module__get__mseigt()

def set_mseigt(mseigt):
    _arespy_pkg.f90wrap_arpack_module__set__mseigt(mseigt)

def get_msapps():
    """
    Element msapps ftype=integer  pytype=int
    
    
    Defined at Arpack_module.fpp line 24
    
    """
    return _arespy_pkg.f90wrap_arpack_module__get__msapps()

def set_msapps(msapps):
    _arespy_pkg.f90wrap_arpack_module__set__msapps(msapps)

def get_msgets():
    """
    Element msgets ftype=integer  pytype=int
    
    
    Defined at Arpack_module.fpp line 24
    
    """
    return _arespy_pkg.f90wrap_arpack_module__get__msgets()

def set_msgets(msgets):
    _arespy_pkg.f90wrap_arpack_module__set__msgets(msgets)

def get_mseupd():
    """
    Element mseupd ftype=integer  pytype=int
    
    
    Defined at Arpack_module.fpp line 24
    
    """
    return _arespy_pkg.f90wrap_arpack_module__get__mseupd()

def set_mseupd(mseupd):
    _arespy_pkg.f90wrap_arpack_module__set__mseupd(mseupd)

def get_mnaupd():
    """
    Element mnaupd ftype=integer  pytype=int
    
    
    Defined at Arpack_module.fpp line 24
    
    """
    return _arespy_pkg.f90wrap_arpack_module__get__mnaupd()

def set_mnaupd(mnaupd):
    _arespy_pkg.f90wrap_arpack_module__set__mnaupd(mnaupd)

def get_mnaup2():
    """
    Element mnaup2 ftype=integer  pytype=int
    
    
    Defined at Arpack_module.fpp line 24
    
    """
    return _arespy_pkg.f90wrap_arpack_module__get__mnaup2()

def set_mnaup2(mnaup2):
    _arespy_pkg.f90wrap_arpack_module__set__mnaup2(mnaup2)

def get_mnaitr():
    """
    Element mnaitr ftype=integer  pytype=int
    
    
    Defined at Arpack_module.fpp line 24
    
    """
    return _arespy_pkg.f90wrap_arpack_module__get__mnaitr()

def set_mnaitr(mnaitr):
    _arespy_pkg.f90wrap_arpack_module__set__mnaitr(mnaitr)

def get_mneigh():
    """
    Element mneigh ftype=integer  pytype=int
    
    
    Defined at Arpack_module.fpp line 24
    
    """
    return _arespy_pkg.f90wrap_arpack_module__get__mneigh()

def set_mneigh(mneigh):
    _arespy_pkg.f90wrap_arpack_module__set__mneigh(mneigh)

def get_mnapps():
    """
    Element mnapps ftype=integer  pytype=int
    
    
    Defined at Arpack_module.fpp line 24
    
    """
    return _arespy_pkg.f90wrap_arpack_module__get__mnapps()

def set_mnapps(mnapps):
    _arespy_pkg.f90wrap_arpack_module__set__mnapps(mnapps)

def get_mngets():
    """
    Element mngets ftype=integer  pytype=int
    
    
    Defined at Arpack_module.fpp line 24
    
    """
    return _arespy_pkg.f90wrap_arpack_module__get__mngets()

def set_mngets(mngets):
    _arespy_pkg.f90wrap_arpack_module__set__mngets(mngets)

def get_mneupd():
    """
    Element mneupd ftype=integer  pytype=int
    
    
    Defined at Arpack_module.fpp line 24
    
    """
    return _arespy_pkg.f90wrap_arpack_module__get__mneupd()

def set_mneupd(mneupd):
    _arespy_pkg.f90wrap_arpack_module__set__mneupd(mneupd)

def get_mcaupd():
    """
    Element mcaupd ftype=integer  pytype=int
    
    
    Defined at Arpack_module.fpp line 24
    
    """
    return _arespy_pkg.f90wrap_arpack_module__get__mcaupd()

def set_mcaupd(mcaupd):
    _arespy_pkg.f90wrap_arpack_module__set__mcaupd(mcaupd)

def get_mcaup2():
    """
    Element mcaup2 ftype=integer  pytype=int
    
    
    Defined at Arpack_module.fpp line 24
    
    """
    return _arespy_pkg.f90wrap_arpack_module__get__mcaup2()

def set_mcaup2(mcaup2):
    _arespy_pkg.f90wrap_arpack_module__set__mcaup2(mcaup2)

def get_mcaitr():
    """
    Element mcaitr ftype=integer  pytype=int
    
    
    Defined at Arpack_module.fpp line 24
    
    """
    return _arespy_pkg.f90wrap_arpack_module__get__mcaitr()

def set_mcaitr(mcaitr):
    _arespy_pkg.f90wrap_arpack_module__set__mcaitr(mcaitr)

def get_mceigh():
    """
    Element mceigh ftype=integer  pytype=int
    
    
    Defined at Arpack_module.fpp line 24
    
    """
    return _arespy_pkg.f90wrap_arpack_module__get__mceigh()

def set_mceigh(mceigh):
    _arespy_pkg.f90wrap_arpack_module__set__mceigh(mceigh)

def get_mcapps():
    """
    Element mcapps ftype=integer  pytype=int
    
    
    Defined at Arpack_module.fpp line 24
    
    """
    return _arespy_pkg.f90wrap_arpack_module__get__mcapps()

def set_mcapps(mcapps):
    _arespy_pkg.f90wrap_arpack_module__set__mcapps(mcapps)

def get_mcgets():
    """
    Element mcgets ftype=integer  pytype=int
    
    
    Defined at Arpack_module.fpp line 24
    
    """
    return _arespy_pkg.f90wrap_arpack_module__get__mcgets()

def set_mcgets(mcgets):
    _arespy_pkg.f90wrap_arpack_module__set__mcgets(mcgets)

def get_mceupd():
    """
    Element mceupd ftype=integer  pytype=int
    
    
    Defined at Arpack_module.fpp line 24
    
    """
    return _arespy_pkg.f90wrap_arpack_module__get__mceupd()

def set_mceupd(mceupd):
    _arespy_pkg.f90wrap_arpack_module__set__mceupd(mceupd)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "arpack_module".')

for func in _dt_array_initialisers:
    func()

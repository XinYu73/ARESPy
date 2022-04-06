"""
Module lapack_module


Defined at Lapack_module.fpp lines 5-526

"""
from __future__ import print_function, absolute_import, division
import _arespy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def diagm(mat, evec, eval):
    """
    diagm(mat, evec, eval)
    
    
    Defined at Lapack_module.fpp lines 20-58
    
    Parameters
    ----------
    mat : complex array
    evec : complex array
    eval : float array
    
    """
    _arespy_pkg.f90wrap_diagm(mat=mat, evec=evec, eval=eval)

def generalizeeigen(dime, mata, matb, evec, eval):
    """
    generalizeeigen(dime, mata, matb, evec, eval)
    
    
    Defined at Lapack_module.fpp lines 67-103
    
    Parameters
    ----------
    dime : int
    mata : complex array
    matb : complex array
    evec : complex array
    eval : float array
    
    """
    _arespy_pkg.f90wrap_generalizeeigen(dime=dime, mata=mata, matb=matb, evec=evec, \
        eval=eval)

def orthnorm(mat):
    """
    orthnorm(mat)
    
    
    Defined at Lapack_module.fpp lines 108-146
    
    Parameters
    ----------
    mat : complex array
    
    """
    _arespy_pkg.f90wrap_orthnorm(mat=mat)

def norm_2(mat, k):
    """
    norm_2 = norm_2(mat, k)
    
    
    Defined at Lapack_module.fpp lines 149-163
    
    Parameters
    ----------
    mat : complex array
    k : int
    
    Returns
    -------
    norm_2 : float
    
    """
    norm_2 = _arespy_pkg.f90wrap_norm_2(mat=mat, k=k)
    return norm_2

def matmat(mata, matb, opa, opb, matc):
    """
    matmat(mata, matb, opa, opb, matc)
    
    
    Defined at Lapack_module.fpp lines 171-201
    
    Parameters
    ----------
    mata : complex array
    matb : complex array
    opa : str
    opb : str
    matc : complex array
    
    """
    _arespy_pkg.f90wrap_matmat(mata=mata, matb=matb, opa=opa, opb=opb, matc=matc)

def invmat(mat):
    """
    invmat(mat)
    
    
    Defined at Lapack_module.fpp lines 204-236
    
    Parameters
    ----------
    mat : complex array
    
    """
    _arespy_pkg.f90wrap_invmat(mat=mat)

def orthnorm_real(mat):
    """
    orthnorm_real(mat)
    
    
    Defined at Lapack_module.fpp lines 260-300
    
    Parameters
    ----------
    mat : float array
    
    """
    _arespy_pkg.f90wrap_orthnorm_real(mat=mat)

def diagm_real(mat, evec, eval):
    """
    diagm_real(mat, evec, eval)
    
    
    Defined at Lapack_module.fpp lines 305-346
    
    Parameters
    ----------
    mat : float array
    evec : float array
    eval : float array
    
    """
    _arespy_pkg.f90wrap_diagm_real(mat=mat, evec=evec, eval=eval)

def generalizeeigen_real(dime, mata, matb, evec, eval):
    """
    generalizeeigen_real(dime, mata, matb, evec, eval)
    
    
    Defined at Lapack_module.fpp lines 354-391
    
    Parameters
    ----------
    dime : int
    mata : float array
    matb : float array
    evec : float array
    eval : float array
    
    """
    _arespy_pkg.f90wrap_generalizeeigen_real(dime=dime, mata=mata, matb=matb, \
        evec=evec, eval=eval)

def diagmx_real(mat, dime, num, il, iu, evec, eval):
    """
    diagmx_real(mat, dime, num, il, iu, evec, eval)
    
    
    Defined at Lapack_module.fpp lines 394-434
    
    Parameters
    ----------
    mat : float array
    dime : int
    num : int
    il : int
    iu : int
    evec : float array
    eval : float array
    
    -------------------OUT PUT-----------------------------
    eigen-values
    """
    _arespy_pkg.f90wrap_diagmx_real(mat=mat, dime=dime, num=num, il=il, iu=iu, \
        evec=evec, eval=eval)

def matmat_real(mata, matb, opa, opb, matc):
    """
    matmat_real(mata, matb, opa, opb, matc)
    
    
    Defined at Lapack_module.fpp lines 437-468
    
    Parameters
    ----------
    mata : float array
    matb : float array
    opa : str
    opb : str
    matc : float array
    
    """
    _arespy_pkg.f90wrap_matmat_real(mata=mata, matb=matb, opa=opa, opb=opb, \
        matc=matc)

def cholesky_factor_real(mat):
    """
    cholesky_factor_real(mat)
    
    
    Defined at Lapack_module.fpp lines 474-490
    
    Parameters
    ----------
    mat : float array
    
    """
    _arespy_pkg.f90wrap_cholesky_factor_real(mat=mat)

def invmat_real(mat):
    """
    invmat_real(mat)
    
    
    Defined at Lapack_module.fpp lines 493-525
    
    Parameters
    ----------
    mat : float array
    
    """
    _arespy_pkg.f90wrap_invmat_real(mat=mat)

def get_mmax():
    """
    Element mmax ftype=integer(i4b) pytype=int
    
    
    Defined at Lapack_module.fpp line 14
    
    """
    return _arespy_pkg.f90wrap_lapack_module__get__mmax()

def set_mmax(mmax):
    _arespy_pkg.f90wrap_lapack_module__set__mmax(mmax)

def get_nmax():
    """
    Element nmax ftype=integer(i4b) pytype=int
    
    
    Defined at Lapack_module.fpp line 14
    
    """
    return _arespy_pkg.f90wrap_lapack_module__get__nmax()

def set_nmax(nmax):
    _arespy_pkg.f90wrap_lapack_module__set__nmax(nmax)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "lapack_module".')

for func in _dt_array_initialisers:
    func()

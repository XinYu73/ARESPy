"""
Module finite_module


Defined at Finite_module.fpp lines 5-592

"""
from __future__ import print_function, absolute_import, division
import _arespy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def destroy_finite():
    """
    destroy_finite()
    
    
    Defined at Finite_module.fpp lines 35-40
    
    
    """
    _arespy_pkg.f90wrap_destroy_finite()

def init_finite(h):
    """
    init_finite(h)
    
    
    Defined at Finite_module.fpp lines 43-109
    
    Parameters
    ----------
    h : float array
    
    """
    _arespy_pkg.f90wrap_init_finite(h=h)

def trans_mat_full(mat, factor, int_miu, lad_gap, err):
    """
    trans_mat_full(mat, factor, int_miu, lad_gap, err)
    
    
    Defined at Finite_module.fpp lines 112-320
    
    Parameters
    ----------
    mat : float array
    factor : float array
    int_miu : int array
    lad_gap : float array
    err : float array
    
    """
    _arespy_pkg.f90wrap_trans_mat_full(mat=mat, factor=factor, int_miu=int_miu, \
        lad_gap=lad_gap, err=err)

def cmplx_keop(uk, ik, ts):
    """
    cmplx_keop(uk, ik, ts)
    
    
    Defined at Finite_module.fpp lines 324-343
    
    Parameters
    ----------
    uk : complex array
    ik : int
    ts : complex array
    
    """
    _arespy_pkg.f90wrap_cmplx_keop(uk=uk, ik=ik, ts=ts)

def real_comm(ifun):
    """
    real_comm(ifun)
    
    
    Defined at Finite_module.fpp lines 350-378
    
    Parameters
    ----------
    ifun : float array
    
    """
    _arespy_pkg.f90wrap_real_comm(ifun=ifun)

def real_comm_clean():
    """
    real_comm_clean()
    
    
    Defined at Finite_module.fpp lines 381-383
    
    
    """
    _arespy_pkg.f90wrap_real_comm_clean()

def real_nabla1_3d(ifun, norder, derf):
    """
    real_nabla1_3d(ifun, norder, derf)
    
    
    Defined at Finite_module.fpp lines 386-434
    
    Parameters
    ----------
    ifun : float array
    norder : int
    derf : float array
    
    """
    _arespy_pkg.f90wrap_real_nabla1_3d(ifun=ifun, norder=norder, derf=derf)

def real_nabla2_3d(func, norder, ofun):
    """
    real_nabla2_3d(func, norder, ofun)
    
    
    Defined at Finite_module.fpp lines 437-504
    
    Parameters
    ----------
    func : float array
    norder : int
    ofun : float array
    
    """
    _arespy_pkg.f90wrap_real_nabla2_3d(func=func, norder=norder, ofun=ofun)

def real_nabla1_1d(ifun, derf, mgfun=None):
    """
    real_nabla1_1d(ifun, derf[, mgfun])
    
    
    Defined at Finite_module.fpp lines 508-548
    
    Parameters
    ----------
    ifun : float array
    derf : float array
    mgfun : float array
    
    """
    _arespy_pkg.f90wrap_real_nabla1_1d(ifun=ifun, derf=derf, mgfun=mgfun)

def real_nabla2_1d(ifun, ofun):
    """
    real_nabla2_1d(ifun, ofun)
    
    
    Defined at Finite_module.fpp lines 551-590
    
    Parameters
    ----------
    ifun : float array
    ofun : float array
    
    """
    _arespy_pkg.f90wrap_real_nabla2_1d(ifun=ifun, ofun=ofun)

def get_array_lapl():
    """
    Element lapl ftype=real(dp) pytype=float
    
    
    Defined at Finite_module.fpp line 19
    
    """
    global lapl
    array_ndim, array_type, array_shape, array_handle = \
        _arespy_pkg.f90wrap_finite_module__array__lapl(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        lapl = _arrays[array_handle]
    else:
        lapl = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_finite_module__array__lapl)
        _arrays[array_handle] = lapl
    return lapl

def set_array_lapl(lapl):
    lapl[...] = lapl

def get_array_grad():
    """
    Element grad ftype=real(dp) pytype=float
    
    
    Defined at Finite_module.fpp line 21
    
    """
    global grad
    array_ndim, array_type, array_shape, array_handle = \
        _arespy_pkg.f90wrap_finite_module__array__grad(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        grad = _arrays[array_handle]
    else:
        grad = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_finite_module__array__grad)
        _arrays[array_handle] = grad
    return grad

def set_array_grad(grad):
    grad[...] = grad

def get_array_tbmat():
    """
    Element tbmat ftype=real(dp) pytype=float
    
    
    Defined at Finite_module.fpp line 22
    
    """
    global tbmat
    array_ndim, array_type, array_shape, array_handle = \
        _arespy_pkg.f90wrap_finite_module__array__tbmat(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        tbmat = _arrays[array_handle]
    else:
        tbmat = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_finite_module__array__tbmat)
        _arrays[array_handle] = tbmat
    return tbmat

def set_array_tbmat(tbmat):
    tbmat[...] = tbmat

def get_array_lap_add():
    """
    Element lap_add ftype=integer(i4b) pytype=int
    
    
    Defined at Finite_module.fpp line 24
    
    """
    global lap_add
    array_ndim, array_type, array_shape, array_handle = \
        _arespy_pkg.f90wrap_finite_module__array__lap_add(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        lap_add = _arrays[array_handle]
    else:
        lap_add = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_finite_module__array__lap_add)
        _arrays[array_handle] = lap_add
    return lap_add

def set_array_lap_add(lap_add):
    lap_add[...] = lap_add

def get_array_cell_mu():
    """
    Element cell_mu ftype=integer  pytype=int
    
    
    Defined at Finite_module.fpp line 26
    
    """
    global cell_mu
    array_ndim, array_type, array_shape, array_handle = \
        _arespy_pkg.f90wrap_finite_module__array__cell_mu(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        cell_mu = _arrays[array_handle]
    else:
        cell_mu = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_finite_module__array__cell_mu)
        _arrays[array_handle] = cell_mu
    return cell_mu

def set_array_cell_mu(cell_mu):
    cell_mu[...] = cell_mu

def get_array_cell_factor():
    """
    Element cell_factor ftype=real(dp) pytype=float
    
    
    Defined at Finite_module.fpp line 27
    
    """
    global cell_factor
    array_ndim, array_type, array_shape, array_handle = \
        _arespy_pkg.f90wrap_finite_module__array__cell_factor(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        cell_factor = _arrays[array_handle]
    else:
        cell_factor = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_finite_module__array__cell_factor)
        _arrays[array_handle] = cell_factor
    return cell_factor

def set_array_cell_factor(cell_factor):
    cell_factor[...] = cell_factor

def get_array_wrap_real():
    """
    Element wrap_real ftype=real(dp) pytype=float
    
    
    Defined at Finite_module.fpp line 28
    
    """
    global wrap_real
    array_ndim, array_type, array_shape, array_handle = \
        _arespy_pkg.f90wrap_finite_module__array__wrap_real(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        wrap_real = _arrays[array_handle]
    else:
        wrap_real = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_finite_module__array__wrap_real)
        _arrays[array_handle] = wrap_real
    return wrap_real

def set_array_wrap_real(wrap_real):
    wrap_real[...] = wrap_real


_array_initialisers = [get_array_lapl, get_array_grad, get_array_tbmat, \
    get_array_lap_add, get_array_cell_mu, get_array_cell_factor, \
    get_array_wrap_real]
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "finite_module".')

for func in _dt_array_initialisers:
    func()

"""
Module finite_module


Defined at Finite_module.fpp lines 5-2216

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def init_finite(norder, h):
    """
    init_finite(norder, h)
    
    
    Defined at Finite_module.fpp lines 39-88
    
    Parameters
    ----------
    norder : int
    h : float array
    
    """
    _AresMainPy_pkg.f90wrap_init_finite(norder=norder, h=h)

def destroy_finite():
    """
    destroy_finite()
    
    
    Defined at Finite_module.fpp lines 91-103
    
    
    """
    _AresMainPy_pkg.f90wrap_destroy_finite()

def trans_mat_full(mat, factor, int_miu, lad_gap, err):
    """
    trans_mat_full(mat, factor, int_miu, lad_gap, err)
    
    
    Defined at Finite_module.fpp lines 106-309
    
    Parameters
    ----------
    mat : float array
    factor : float array
    int_miu : int array
    lad_gap : float array
    err : float array
    
    """
    _AresMainPy_pkg.f90wrap_trans_mat_full(mat=mat, factor=factor, int_miu=int_miu, \
        lad_gap=lad_gap, err=err)

def cmplx_nabla2(ifun, norder, ofun):
    """
    cmplx_nabla2(ifun, norder, ofun)
    
    
    Defined at Finite_module.fpp lines 312-439
    
    Parameters
    ----------
    ifun : complex array
    norder : int
    ofun : complex array
    
    """
    _AresMainPy_pkg.f90wrap_cmplx_nabla2(ifun=ifun, norder=norder, ofun=ofun)

def cmplx_nabla1(ifun, norder, derf):
    """
    cmplx_nabla1(ifun, norder, derf)
    
    
    Defined at Finite_module.fpp lines 442-634
    
    Parameters
    ----------
    ifun : complex array
    norder : int
    derf : complex array
    
    """
    _AresMainPy_pkg.f90wrap_cmplx_nabla1(ifun=ifun, norder=norder, derf=derf)

def nabla(ifun, norder, derf):
    """
    nabla(ifun, norder, derf)
    
    
    Defined at Finite_module.fpp lines 637-741
    
    Parameters
    ----------
    ifun : complex array
    norder : int
    derf : complex array
    
    """
    _AresMainPy_pkg.f90wrap_nabla(ifun=ifun, norder=norder, derf=derf)

def kskeperiod_op(psi, ik, kepsi):
    """
    kskeperiod_op(psi, ik, kepsi)
    
    
    Defined at Finite_module.fpp lines 743-866
    
    Parameters
    ----------
    psi : complex array
    ik : int
    kepsi : complex array
    
    """
    _AresMainPy_pkg.f90wrap_kskeperiod_op(psi=psi, ik=ik, kepsi=kepsi)

def ke1(psi, ts1):
    """
    ke1(psi, ts1)
    
    
    Defined at Finite_module.fpp lines 869-880
    
    Parameters
    ----------
    psi : complex array
    ts1 : complex array
    
    """
    _AresMainPy_pkg.f90wrap_ke1(psi=psi, ts1=ts1)

def ke2(psi, ik, ts2):
    """
    ke2(psi, ik, ts2)
    
    
    Defined at Finite_module.fpp lines 883-909
    
    Parameters
    ----------
    psi : complex array
    ik : int
    ts2 : complex array
    
    """
    _AresMainPy_pkg.f90wrap_ke2(psi=psi, ik=ik, ts2=ts2)

def realnabla2(ifun, norder, ofun):
    """
    realnabla2(ifun, norder, ofun)
    
    
    Defined at Finite_module.fpp lines 918-1008
    
    Parameters
    ----------
    ifun : float array
    norder : int
    ofun : float array
    
    """
    _AresMainPy_pkg.f90wrap_realnabla2(ifun=ifun, norder=norder, ofun=ofun)

def realke(psi, kepsi):
    """
    realke(psi, kepsi)
    
    
    Defined at Finite_module.fpp lines 1011-1023
    
    Parameters
    ----------
    psi : float array
    kepsi : float array
    
    """
    _AresMainPy_pkg.f90wrap_realke(psi=psi, kepsi=kepsi)

def iso_realnabla2(ifun, norder, ofun):
    """
    iso_realnabla2(ifun, norder, ofun)
    
    
    Defined at Finite_module.fpp lines 1033-1188
    
    Parameters
    ----------
    ifun : float array
    norder : int
    ofun : float array
    
    ==================================================
    """
    _AresMainPy_pkg.f90wrap_iso_realnabla2(ifun=ifun, norder=norder, ofun=ofun)

def iso_realke(psi, kepsi):
    """
    iso_realke(psi, kepsi)
    
    
    Defined at Finite_module.fpp lines 1243-1255
    
    Parameters
    ----------
    psi : float array
    kepsi : float array
    
    """
    _AresMainPy_pkg.f90wrap_iso_realke(psi=psi, kepsi=kepsi)

def kskeperiod_op_band(psi, ik, kepsi):
    """
    kskeperiod_op_band(psi, ik, kepsi)
    
    
    Defined at Finite_module.fpp lines 1258-1308
    
    Parameters
    ----------
    psi : complex array
    ik : int
    kepsi : complex array
    
    """
    _AresMainPy_pkg.f90wrap_kskeperiod_op_band(psi=psi, ik=ik, kepsi=kepsi)

def nabla2_np(ifun, norder, ofun):
    """
    nabla2_np(ifun, norder, ofun)
    
    
    Defined at Finite_module.fpp lines 1310-1431
    
    Parameters
    ----------
    ifun : complex array
    norder : int
    ofun : complex array
    
    """
    _AresMainPy_pkg.f90wrap_nabla2_np(ifun=ifun, norder=norder, ofun=ofun)

def nabla1_np(ifun, norder, derf):
    """
    nabla1_np(ifun, norder, derf)
    
    
    Defined at Finite_module.fpp lines 1434-1512
    
    Parameters
    ----------
    ifun : complex array
    norder : int
    derf : complex array
    
    """
    _AresMainPy_pkg.f90wrap_nabla1_np(ifun=ifun, norder=norder, derf=derf)

def ke1_np(psi, ts1):
    """
    ke1_np(psi, ts1)
    
    
    Defined at Finite_module.fpp lines 1515-1526
    
    Parameters
    ----------
    psi : complex array
    ts1 : complex array
    
    """
    _AresMainPy_pkg.f90wrap_ke1_np(psi=psi, ts1=ts1)

def ke2_np(psi, ik, ts2):
    """
    ke2_np(psi, ik, ts2)
    
    
    Defined at Finite_module.fpp lines 1529-1555
    
    Parameters
    ----------
    psi : complex array
    ik : int
    ts2 : complex array
    
    """
    _AresMainPy_pkg.f90wrap_ke2_np(psi=psi, ik=ik, ts2=ts2)

def kskeperiod_op_gamma(psi, ik, kepsi):
    """
    kskeperiod_op_gamma(psi, ik, kepsi)
    
    
    Defined at Finite_module.fpp lines 1557-1654
    
    Parameters
    ----------
    psi : float array
    ik : int
    kepsi : float array
    
    """
    _AresMainPy_pkg.f90wrap_kskeperiod_op_gamma(psi=psi, ik=ik, kepsi=kepsi)

def ke_gamma(psi, ts1):
    """
    ke_gamma(psi, ts1)
    
    
    Defined at Finite_module.fpp lines 1656-1668
    
    Parameters
    ----------
    psi : float array
    ts1 : float array
    
    """
    _AresMainPy_pkg.f90wrap_ke_gamma(psi=psi, ts1=ts1)

def real_nabla2(ifun, norder, ofun):
    """
    real_nabla2(ifun, norder, ofun)
    
    
    Defined at Finite_module.fpp lines 1670-1795
    
    Parameters
    ----------
    ifun : float array
    norder : int
    ofun : float array
    
    """
    _AresMainPy_pkg.f90wrap_real_nabla2(ifun=ifun, norder=norder, ofun=ofun)

def real_nabla2_1d(ifun, ofun):
    """
    real_nabla2_1d(ifun, ofun)
    
    
    Defined at Finite_module.fpp lines 1797-1928
    
    Parameters
    ----------
    ifun : float array
    ofun : float array
    
    """
    _AresMainPy_pkg.f90wrap_real_nabla2_1d(ifun=ifun, ofun=ofun)

def nabla_gamma(ifun, norder, derf):
    """
    nabla_gamma(ifun, norder, derf)
    
    
    Defined at Finite_module.fpp lines 1930-2022
    
    Parameters
    ----------
    ifun : float array
    norder : int
    derf : float array
    
    """
    _AresMainPy_pkg.f90wrap_nabla_gamma(ifun=ifun, norder=norder, derf=derf)

def real_nabla1(ifun, norder, derf):
    """
    real_nabla1(ifun, norder, derf)
    
    
    Defined at Finite_module.fpp lines 2024-2216
    
    Parameters
    ----------
    ifun : float array
    norder : int
    derf : float array
    
    """
    _AresMainPy_pkg.f90wrap_real_nabla1(ifun=ifun, norder=norder, derf=derf)

def get_array_lapl():
    """
    Element lapl ftype=real(dp) pytype=float
    
    
    Defined at Finite_module.fpp line 19
    
    """
    global lapl
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_finite_module__array__lapl(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        lapl = _arrays[array_handle]
    else:
        lapl = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_finite_module__array__lapl)
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
        _AresMainPy_pkg.f90wrap_finite_module__array__grad(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        grad = _arrays[array_handle]
    else:
        grad = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_finite_module__array__grad)
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
        _AresMainPy_pkg.f90wrap_finite_module__array__tbmat(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        tbmat = _arrays[array_handle]
    else:
        tbmat = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_finite_module__array__tbmat)
        _arrays[array_handle] = tbmat
    return tbmat

def set_array_tbmat(tbmat):
    tbmat[...] = tbmat

def get_array_lap_add():
    """
    Element lap_add ftype=integer(i4b) pytype=int
    
    
    Defined at Finite_module.fpp line 25
    
    """
    global lap_add
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_finite_module__array__lap_add(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        lap_add = _arrays[array_handle]
    else:
        lap_add = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_finite_module__array__lap_add)
        _arrays[array_handle] = lap_add
    return lap_add

def set_array_lap_add(lap_add):
    lap_add[...] = lap_add

def get_array_cell_mu():
    """
    Element cell_mu ftype=integer  pytype=int
    
    
    Defined at Finite_module.fpp line 27
    
    """
    global cell_mu
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_finite_module__array__cell_mu(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        cell_mu = _arrays[array_handle]
    else:
        cell_mu = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_finite_module__array__cell_mu)
        _arrays[array_handle] = cell_mu
    return cell_mu

def set_array_cell_mu(cell_mu):
    cell_mu[...] = cell_mu

def get_array_cell_factor():
    """
    Element cell_factor ftype=real(dp) pytype=float
    
    
    Defined at Finite_module.fpp line 28
    
    """
    global cell_factor
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_finite_module__array__cell_factor(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        cell_factor = _arrays[array_handle]
    else:
        cell_factor = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_finite_module__array__cell_factor)
        _arrays[array_handle] = cell_factor
    return cell_factor

def set_array_cell_factor(cell_factor):
    cell_factor[...] = cell_factor

def get_array_wrap_box():
    """
    Element wrap_box ftype=complex(dcp) pytype=complex
    
    
    Defined at Finite_module.fpp line 29
    
    """
    global wrap_box
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_finite_module__array__wrap_box(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        wrap_box = _arrays[array_handle]
    else:
        wrap_box = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_finite_module__array__wrap_box)
        _arrays[array_handle] = wrap_box
    return wrap_box

def set_array_wrap_box(wrap_box):
    wrap_box[...] = wrap_box

def get_array_fun_global():
    """
    Element fun_global ftype=complex(dcp) pytype=complex
    
    
    Defined at Finite_module.fpp line 29
    
    """
    global fun_global
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_finite_module__array__fun_global(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        fun_global = _arrays[array_handle]
    else:
        fun_global = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_finite_module__array__fun_global)
        _arrays[array_handle] = fun_global
    return fun_global

def set_array_fun_global(fun_global):
    fun_global[...] = fun_global

def get_array_fun_1d():
    """
    Element fun_1d ftype=complex(dcp) pytype=complex
    
    
    Defined at Finite_module.fpp line 30
    
    """
    global fun_1d
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_finite_module__array__fun_1d(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        fun_1d = _arrays[array_handle]
    else:
        fun_1d = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_finite_module__array__fun_1d)
        _arrays[array_handle] = fun_1d
    return fun_1d

def set_array_fun_1d(fun_1d):
    fun_1d[...] = fun_1d

def get_array_wrap_box1d():
    """
    Element wrap_box1d ftype=complex(dcp) pytype=complex
    
    
    Defined at Finite_module.fpp line 30
    
    """
    global wrap_box1d
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_finite_module__array__wrap_box1d(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        wrap_box1d = _arrays[array_handle]
    else:
        wrap_box1d = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_finite_module__array__wrap_box1d)
        _arrays[array_handle] = wrap_box1d
    return wrap_box1d

def set_array_wrap_box1d(wrap_box1d):
    wrap_box1d[...] = wrap_box1d

def get_array_wrap_box_real():
    """
    Element wrap_box_real ftype=real(dp) pytype=float
    
    
    Defined at Finite_module.fpp line 31
    
    """
    global wrap_box_real
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_finite_module__array__wrap_box_real(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        wrap_box_real = _arrays[array_handle]
    else:
        wrap_box_real = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_finite_module__array__wrap_box_real)
        _arrays[array_handle] = wrap_box_real
    return wrap_box_real

def set_array_wrap_box_real(wrap_box_real):
    wrap_box_real[...] = wrap_box_real

def get_array_fun_global_real():
    """
    Element fun_global_real ftype=real(dp) pytype=float
    
    
    Defined at Finite_module.fpp line 31
    
    """
    global fun_global_real
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_finite_module__array__fun_global_real(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        fun_global_real = _arrays[array_handle]
    else:
        fun_global_real = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_finite_module__array__fun_global_real)
        _arrays[array_handle] = fun_global_real
    return fun_global_real

def set_array_fun_global_real(fun_global_real):
    fun_global_real[...] = fun_global_real

def get_array_fun_1d_real():
    """
    Element fun_1d_real ftype=real(dp) pytype=float
    
    
    Defined at Finite_module.fpp line 32
    
    """
    global fun_1d_real
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_finite_module__array__fun_1d_real(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        fun_1d_real = _arrays[array_handle]
    else:
        fun_1d_real = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_finite_module__array__fun_1d_real)
        _arrays[array_handle] = fun_1d_real
    return fun_1d_real

def set_array_fun_1d_real(fun_1d_real):
    fun_1d_real[...] = fun_1d_real

def get_array_wrap_box1d_real():
    """
    Element wrap_box1d_real ftype=real(dp) pytype=float
    
    
    Defined at Finite_module.fpp line 32
    
    """
    global wrap_box1d_real
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_finite_module__array__wrap_box1d_real(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        wrap_box1d_real = _arrays[array_handle]
    else:
        wrap_box1d_real = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_finite_module__array__wrap_box1d_real)
        _arrays[array_handle] = wrap_box1d_real
    return wrap_box1d_real

def set_array_wrap_box1d_real(wrap_box1d_real):
    wrap_box1d_real[...] = wrap_box1d_real


_array_initialisers = [get_array_lapl, get_array_grad, get_array_tbmat, \
    get_array_lap_add, get_array_cell_mu, get_array_cell_factor, \
    get_array_wrap_box, get_array_fun_global, get_array_fun_1d, \
    get_array_wrap_box1d, get_array_wrap_box_real, get_array_fun_global_real, \
    get_array_fun_1d_real, get_array_wrap_box1d_real]
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "finite_module".')

for func in _dt_array_initialisers:
    func()

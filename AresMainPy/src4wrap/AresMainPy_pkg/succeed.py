"""
Module succeed


Defined at Succeed_module.fpp lines 5-487

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def init_succeed_rho_real(n_rho, n_r, n_s, nspin, dv):
    """
    init_succeed_rho_real(n_rho, n_r, n_s, nspin, dv)
    
    
    Defined at Succeed_module.fpp lines 31-70
    
    Parameters
    ----------
    n_rho : int
    n_r : int
    n_s : int
    nspin : int
    dv : float
    
    """
    _AresMainPy_pkg.f90wrap_init_succeed_rho_real(n_rho=n_rho, n_r=n_r, n_s=n_s, \
        nspin=nspin, dv=dv)

def init_succeed_rho_cmplx(n_rho, n_r, n_s, n_k, nspin, dv):
    """
    init_succeed_rho_cmplx(n_rho, n_r, n_s, n_k, nspin, dv)
    
    
    Defined at Succeed_module.fpp lines 72-113
    
    Parameters
    ----------
    n_rho : int
    n_r : int
    n_s : int
    n_k : int
    nspin : int
    dv : float
    
    """
    _AresMainPy_pkg.f90wrap_init_succeed_rho_cmplx(n_rho=n_rho, n_r=n_r, n_s=n_s, \
        n_k=n_k, nspin=nspin, dv=dv)

def destroy_succeed():
    """
    destroy_succeed()
    
    
    Defined at Succeed_module.fpp lines 115-124
    
    
    """
    _AresMainPy_pkg.f90wrap_destroy_succeed()

def store_rho(n_rho, nspin, rho, dvol):
    """
    store_rho(n_rho, nspin, rho, dvol)
    
    
    Defined at Succeed_module.fpp lines 126-140
    
    Parameters
    ----------
    n_rho : int
    nspin : int
    rho : float array
    dvol : float
    
    """
    _AresMainPy_pkg.f90wrap_store_rho(n_rho=n_rho, nspin=nspin, rho=rho, dvol=dvol)

def store_rho_at(n_rho, nspin, rho_in):
    """
    store_rho_at(n_rho, nspin, rho_in)
    
    
    Defined at Succeed_module.fpp lines 142-147
    
    Parameters
    ----------
    n_rho : int
    nspin : int
    rho_in : float array
    
    """
    _AresMainPy_pkg.f90wrap_store_rho_at(n_rho=n_rho, nspin=nspin, rho_in=rho_in)

def store_r(nr, r):
    """
    store_r(nr, r)
    
    
    Defined at Succeed_module.fpp lines 149-158
    
    Parameters
    ----------
    nr : int
    r : float array
    
    """
    _AresMainPy_pkg.f90wrap_store_r(nr=nr, r=r)

def get_rho(nr, r_new, nrho, nspin, rho_new, dvol):
    """
    get_rho(nr, r_new, nrho, nspin, rho_new, dvol)
    
    
    Defined at Succeed_module.fpp lines 186-267
    
    Parameters
    ----------
    nr : int
    r_new : float array
    nrho : int
    nspin : int
    rho_new : float array
    dvol : float
    
    """
    _AresMainPy_pkg.f90wrap_get_rho(nr=nr, r_new=r_new, nrho=nrho, nspin=nspin, \
        rho_new=rho_new, dvol=dvol)

def cal_trans_phase(nr, nspin, r_new, n1xy, n2xy, n3xy, ng1, ng2, ng3, gvec, \
    trans_phase):
    """
    cal_trans_phase(nr, nspin, r_new, n1xy, n2xy, n3xy, ng1, ng2, ng3, gvec, \
        trans_phase)
    
    
    Defined at Succeed_module.fpp lines 323-375
    
    Parameters
    ----------
    nr : int
    nspin : int
    r_new : float array
    n1xy : int
    n2xy : int
    n3xy : int
    ng1 : int
    ng2 : int
    ng3 : int
    gvec : float array
    trans_phase : complex array
    
    """
    _AresMainPy_pkg.f90wrap_cal_trans_phase(nr=nr, nspin=nspin, r_new=r_new, \
        n1xy=n1xy, n2xy=n2xy, n3xy=n3xy, ng1=ng1, ng2=ng2, ng3=ng3, gvec=gvec, \
        trans_phase=trans_phase)

def get_new_rho_psi(nr, r_new, nrho, n1xy, n2xy, n3xy, nspin, rho_new, n_s, \
    psi_new, gvec):
    """
    get_new_rho_psi(nr, r_new, nrho, n1xy, n2xy, n3xy, nspin, rho_new, n_s, psi_new, \
        gvec)
    
    
    Defined at Succeed_module.fpp lines 377-436
    
    Parameters
    ----------
    nr : int
    r_new : float array
    nrho : int
    n1xy : int
    n2xy : int
    n3xy : int
    nspin : int
    rho_new : float array
    n_s : int
    psi_new : float array
    gvec : float array
    
    """
    _AresMainPy_pkg.f90wrap_get_new_rho_psi(nr=nr, r_new=r_new, nrho=nrho, \
        n1xy=n1xy, n2xy=n2xy, n3xy=n3xy, nspin=nspin, rho_new=rho_new, n_s=n_s, \
        psi_new=psi_new, gvec=gvec)

def store_rho_fft_trans(n_rho, nspin, rho):
    """
    store_rho_fft_trans(n_rho, nspin, rho)
    
    
    Defined at Succeed_module.fpp lines 438-446
    
    Parameters
    ----------
    n_rho : int
    nspin : int
    rho : float array
    
    """
    _AresMainPy_pkg.f90wrap_store_rho_fft_trans(n_rho=n_rho, nspin=nspin, rho=rho)

def store_rho_at_fft_trans(n_rho, nspin, na, rho_in, rho_in2):
    """
    store_rho_at_fft_trans(n_rho, nspin, na, rho_in, rho_in2)
    
    
    Defined at Succeed_module.fpp lines 448-463
    
    Parameters
    ----------
    n_rho : int
    nspin : int
    na : int
    rho_in : float array
    rho_in2 : float array
    
    """
    _AresMainPy_pkg.f90wrap_store_rho_at_fft_trans(n_rho=n_rho, nspin=nspin, na=na, \
        rho_in=rho_in, rho_in2=rho_in2)

def store_r_fft_trans(nr, r):
    """
    store_r_fft_trans(nr, r)
    
    
    Defined at Succeed_module.fpp lines 465-473
    
    Parameters
    ----------
    nr : int
    r : float array
    
    """
    _AresMainPy_pkg.f90wrap_store_r_fft_trans(nr=nr, r=r)

def store_psi_fft_trans(n_rho, n_s, nspin, psi):
    """
    store_psi_fft_trans(n_rho, n_s, nspin, psi)
    
    
    Defined at Succeed_module.fpp lines 475-487
    
    Parameters
    ----------
    n_rho : int
    n_s : int
    nspin : int
    psi : float array
    
    """
    _AresMainPy_pkg.f90wrap_store_psi_fft_trans(n_rho=n_rho, n_s=n_s, nspin=nspin, \
        psi=psi)

def _store_psi_cmplx(n_rho, n_s, n_k, nspin, psi):
    """
    _store_psi_cmplx(n_rho, n_s, n_k, nspin, psi)
    
    
    Defined at Succeed_module.fpp lines 173-184
    
    Parameters
    ----------
    n_rho : int
    n_s : int
    n_k : int
    nspin : int
    psi : complex array
    
    """
    _AresMainPy_pkg.f90wrap_store_psi_cmplx(n_rho=n_rho, n_s=n_s, n_k=n_k, \
        nspin=nspin, psi=psi)

def _store_psi_real(n_rho, n_s, nspin, psi):
    """
    _store_psi_real(n_rho, n_s, nspin, psi)
    
    
    Defined at Succeed_module.fpp lines 160-171
    
    Parameters
    ----------
    n_rho : int
    n_s : int
    nspin : int
    psi : float array
    
    """
    _AresMainPy_pkg.f90wrap_store_psi_real(n_rho=n_rho, n_s=n_s, nspin=nspin, \
        psi=psi)

def store_psi(*args, **kwargs):
    """
    store_psi(*args, **kwargs)
    
    
    Defined at Succeed_module.fpp lines 22-24
    
    Overloaded interface containing the following procedures:
      _store_psi_cmplx
      _store_psi_real
    
    """
    for proc in [_store_psi_cmplx, _store_psi_real]:
        try:
            return proc(*args, **kwargs)
        except TypeError:
            continue

def _get_psi_cmplx(nrho, n_s, n_k, nspin, psi_new):
    """
    _get_psi_cmplx(nrho, n_s, n_k, nspin, psi_new)
    
    
    Defined at Succeed_module.fpp lines 292-313
    
    Parameters
    ----------
    nrho : int
    n_s : int
    n_k : int
    nspin : int
    psi_new : complex array
    
    """
    _AresMainPy_pkg.f90wrap_get_psi_cmplx(nrho=nrho, n_s=n_s, n_k=n_k, nspin=nspin, \
        psi_new=psi_new)

def _get_psi_real(nrho, n_s, nspin, psi_new):
    """
    _get_psi_real(nrho, n_s, nspin, psi_new)
    
    
    Defined at Succeed_module.fpp lines 269-290
    
    Parameters
    ----------
    nrho : int
    n_s : int
    nspin : int
    psi_new : float array
    
    """
    _AresMainPy_pkg.f90wrap_get_psi_real(nrho=nrho, n_s=n_s, nspin=nspin, \
        psi_new=psi_new)

def get_psi(*args, **kwargs):
    """
    get_psi(*args, **kwargs)
    
    
    Defined at Succeed_module.fpp lines 26-28
    
    Overloaded interface containing the following procedures:
      _get_psi_cmplx
      _get_psi_real
    
    """
    for proc in [_get_psi_cmplx, _get_psi_real]:
        try:
            return proc(*args, **kwargs)
        except TypeError:
            continue

def get_llastrho():
    """
    Element llastrho ftype=logical pytype=bool
    
    
    Defined at Succeed_module.fpp line 12
    
    """
    return _AresMainPy_pkg.f90wrap_succeed__get__llastrho()

def set_llastrho(llastrho):
    _AresMainPy_pkg.f90wrap_succeed__set__llastrho(llastrho)

def get_lsr():
    """
    Element lsr ftype=logical pytype=bool
    
    
    Defined at Succeed_module.fpp line 12
    
    """
    return _AresMainPy_pkg.f90wrap_succeed__get__lsr()

def set_lsr(lsr):
    _AresMainPy_pkg.f90wrap_succeed__set__lsr(lsr)

def get_lsrho():
    """
    Element lsrho ftype=logical pytype=bool
    
    
    Defined at Succeed_module.fpp line 12
    
    """
    return _AresMainPy_pkg.f90wrap_succeed__get__lsrho()

def set_lsrho(lsrho):
    _AresMainPy_pkg.f90wrap_succeed__set__lsrho(lsrho)

def get_lspsi():
    """
    Element lspsi ftype=logical pytype=bool
    
    
    Defined at Succeed_module.fpp line 12
    
    """
    return _AresMainPy_pkg.f90wrap_succeed__get__lspsi()

def set_lspsi(lspsi):
    _AresMainPy_pkg.f90wrap_succeed__set__lspsi(lspsi)

def get_lsrho_at():
    """
    Element lsrho_at ftype=logical pytype=bool
    
    
    Defined at Succeed_module.fpp line 12
    
    """
    return _AresMainPy_pkg.f90wrap_succeed__get__lsrho_at()

def set_lsrho_at(lsrho_at):
    _AresMainPy_pkg.f90wrap_succeed__set__lsrho_at(lsrho_at)

def get_array_rho1():
    """
    Element rho1 ftype=real(dp) pytype=float
    
    
    Defined at Succeed_module.fpp line 17
    
    """
    global rho1
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_succeed__array__rho1(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        rho1 = _arrays[array_handle]
    else:
        rho1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_succeed__array__rho1)
        _arrays[array_handle] = rho1
    return rho1

def set_array_rho1(rho1):
    rho1[...] = rho1

def get_array_rho2():
    """
    Element rho2 ftype=real(dp) pytype=float
    
    
    Defined at Succeed_module.fpp line 17
    
    """
    global rho2
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_succeed__array__rho2(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        rho2 = _arrays[array_handle]
    else:
        rho2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_succeed__array__rho2)
        _arrays[array_handle] = rho2
    return rho2

def set_array_rho2(rho2):
    rho2[...] = rho2

def get_array_rho3():
    """
    Element rho3 ftype=real(dp) pytype=float
    
    
    Defined at Succeed_module.fpp line 17
    
    """
    global rho3
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_succeed__array__rho3(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        rho3 = _arrays[array_handle]
    else:
        rho3 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_succeed__array__rho3)
        _arrays[array_handle] = rho3
    return rho3

def set_array_rho3(rho3):
    rho3[...] = rho3

def get_array_r1():
    """
    Element r1 ftype=real(dp) pytype=float
    
    
    Defined at Succeed_module.fpp line 17
    
    """
    global r1
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_succeed__array__r1(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        r1 = _arrays[array_handle]
    else:
        r1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_succeed__array__r1)
        _arrays[array_handle] = r1
    return r1

def set_array_r1(r1):
    r1[...] = r1

def get_array_r2():
    """
    Element r2 ftype=real(dp) pytype=float
    
    
    Defined at Succeed_module.fpp line 17
    
    """
    global r2
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_succeed__array__r2(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        r2 = _arrays[array_handle]
    else:
        r2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_succeed__array__r2)
        _arrays[array_handle] = r2
    return r2

def set_array_r2(r2):
    r2[...] = r2

def get_array_r3():
    """
    Element r3 ftype=real(dp) pytype=float
    
    
    Defined at Succeed_module.fpp line 17
    
    """
    global r3
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_succeed__array__r3(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        r3 = _arrays[array_handle]
    else:
        r3 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_succeed__array__r3)
        _arrays[array_handle] = r3
    return r3

def set_array_r3(r3):
    r3[...] = r3

def get_array_psi1():
    """
    Element psi1 ftype=real(dp) pytype=float
    
    
    Defined at Succeed_module.fpp line 17
    
    """
    global psi1
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_succeed__array__psi1(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        psi1 = _arrays[array_handle]
    else:
        psi1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_succeed__array__psi1)
        _arrays[array_handle] = psi1
    return psi1

def set_array_psi1(psi1):
    psi1[...] = psi1

def get_array_psi2():
    """
    Element psi2 ftype=real(dp) pytype=float
    
    
    Defined at Succeed_module.fpp line 17
    
    """
    global psi2
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_succeed__array__psi2(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        psi2 = _arrays[array_handle]
    else:
        psi2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_succeed__array__psi2)
        _arrays[array_handle] = psi2
    return psi2

def set_array_psi2(psi2):
    psi2[...] = psi2

def get_array_psi3():
    """
    Element psi3 ftype=real(dp) pytype=float
    
    
    Defined at Succeed_module.fpp line 17
    
    """
    global psi3
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_succeed__array__psi3(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        psi3 = _arrays[array_handle]
    else:
        psi3 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_succeed__array__psi3)
        _arrays[array_handle] = psi3
    return psi3

def set_array_psi3(psi3):
    psi3[...] = psi3

def get_array_rho_at():
    """
    Element rho_at ftype=real(dp) pytype=float
    
    
    Defined at Succeed_module.fpp line 17
    
    """
    global rho_at
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_succeed__array__rho_at(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        rho_at = _arrays[array_handle]
    else:
        rho_at = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_succeed__array__rho_at)
        _arrays[array_handle] = rho_at
    return rho_at

def set_array_rho_at(rho_at):
    rho_at[...] = rho_at

def get_array_rhoi():
    """
    Element rhoi ftype=real(dp) pytype=float
    
    
    Defined at Succeed_module.fpp line 17
    
    """
    global rhoi
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_succeed__array__rhoi(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        rhoi = _arrays[array_handle]
    else:
        rhoi = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_succeed__array__rhoi)
        _arrays[array_handle] = rhoi
    return rhoi

def set_array_rhoi(rhoi):
    rhoi[...] = rhoi

def get_array_rho_at1():
    """
    Element rho_at1 ftype=real(dp) pytype=float
    
    
    Defined at Succeed_module.fpp line 17
    
    """
    global rho_at1
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_succeed__array__rho_at1(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        rho_at1 = _arrays[array_handle]
    else:
        rho_at1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_succeed__array__rho_at1)
        _arrays[array_handle] = rho_at1
    return rho_at1

def set_array_rho_at1(rho_at1):
    rho_at1[...] = rho_at1

def get_array_rhoi1():
    """
    Element rhoi1 ftype=real(dp) pytype=float
    
    
    Defined at Succeed_module.fpp line 17
    
    """
    global rhoi1
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_succeed__array__rhoi1(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        rhoi1 = _arrays[array_handle]
    else:
        rhoi1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_succeed__array__rhoi1)
        _arrays[array_handle] = rhoi1
    return rhoi1

def set_array_rhoi1(rhoi1):
    rhoi1[...] = rhoi1

def get_array_psi1_c():
    """
    Element psi1_c ftype=complex(dcp) pytype=complex
    
    
    Defined at Succeed_module.fpp line 18
    
    """
    global psi1_c
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_succeed__array__psi1_c(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        psi1_c = _arrays[array_handle]
    else:
        psi1_c = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_succeed__array__psi1_c)
        _arrays[array_handle] = psi1_c
    return psi1_c

def set_array_psi1_c(psi1_c):
    psi1_c[...] = psi1_c

def get_array_psi2_c():
    """
    Element psi2_c ftype=complex(dcp) pytype=complex
    
    
    Defined at Succeed_module.fpp line 18
    
    """
    global psi2_c
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_succeed__array__psi2_c(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        psi2_c = _arrays[array_handle]
    else:
        psi2_c = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_succeed__array__psi2_c)
        _arrays[array_handle] = psi2_c
    return psi2_c

def set_array_psi2_c(psi2_c):
    psi2_c[...] = psi2_c

def get_array_psi3_c():
    """
    Element psi3_c ftype=complex(dcp) pytype=complex
    
    
    Defined at Succeed_module.fpp line 18
    
    """
    global psi3_c
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_succeed__array__psi3_c(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        psi3_c = _arrays[array_handle]
    else:
        psi3_c = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_succeed__array__psi3_c)
        _arrays[array_handle] = psi3_c
    return psi3_c

def set_array_psi3_c(psi3_c):
    psi3_c[...] = psi3_c

def get_counter1():
    """
    Element counter1 ftype=integer(i4b) pytype=int
    
    
    Defined at Succeed_module.fpp line 20
    
    """
    return _AresMainPy_pkg.f90wrap_succeed__get__counter1()

def set_counter1(counter1):
    _AresMainPy_pkg.f90wrap_succeed__set__counter1(counter1)

def get_counter2():
    """
    Element counter2 ftype=integer(i4b) pytype=int
    
    
    Defined at Succeed_module.fpp line 20
    
    """
    return _AresMainPy_pkg.f90wrap_succeed__get__counter2()

def set_counter2(counter2):
    _AresMainPy_pkg.f90wrap_succeed__set__counter2(counter2)

def get_counter3():
    """
    Element counter3 ftype=integer(i4b) pytype=int
    
    
    Defined at Succeed_module.fpp line 20
    
    """
    return _AresMainPy_pkg.f90wrap_succeed__get__counter3()

def set_counter3(counter3):
    _AresMainPy_pkg.f90wrap_succeed__set__counter3(counter3)

def get_alpha():
    """
    Element alpha ftype=real(dp) pytype=float
    
    
    Defined at Succeed_module.fpp line 21
    
    """
    return _AresMainPy_pkg.f90wrap_succeed__get__alpha()

def set_alpha(alpha):
    _AresMainPy_pkg.f90wrap_succeed__set__alpha(alpha)

def get_beta():
    """
    Element beta ftype=real(dp) pytype=float
    
    
    Defined at Succeed_module.fpp line 21
    
    """
    return _AresMainPy_pkg.f90wrap_succeed__get__beta()

def set_beta(beta):
    _AresMainPy_pkg.f90wrap_succeed__set__beta(beta)


_array_initialisers = [get_array_rho1, get_array_rho2, get_array_rho3, \
    get_array_r1, get_array_r2, get_array_r3, get_array_psi1, get_array_psi2, \
    get_array_psi3, get_array_rho_at, get_array_rhoi, get_array_rho_at1, \
    get_array_rhoi1, get_array_psi1_c, get_array_psi2_c, get_array_psi3_c]
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "succeed".')

for func in _dt_array_initialisers:
    func()

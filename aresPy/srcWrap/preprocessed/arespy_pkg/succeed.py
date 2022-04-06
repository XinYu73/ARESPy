"""
Module succeed


Defined at Succeed_module.fpp lines 5-396

"""
from __future__ import print_function, absolute_import, division
import _arespy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def init_succeed_rho(n_rho, n_r, n_s, nspin, dv):
    """
    init_succeed_rho(n_rho, n_r, n_s, nspin, dv)
    
    
    Defined at Succeed_module.fpp lines 23-62
    
    Parameters
    ----------
    n_rho : int
    n_r : int
    n_s : int
    nspin : int
    dv : float
    
    """
    _arespy_pkg.f90wrap_init_succeed_rho(n_rho=n_rho, n_r=n_r, n_s=n_s, nspin=nspin, \
        dv=dv)

def destory_succeed():
    """
    destory_succeed()
    
    
    Defined at Succeed_module.fpp lines 64-71
    
    
    """
    _arespy_pkg.f90wrap_destory_succeed()

def store_rho(n_rho, nspin, rho):
    """
    store_rho(n_rho, nspin, rho)
    
    
    Defined at Succeed_module.fpp lines 73-86
    
    Parameters
    ----------
    n_rho : int
    nspin : int
    rho : float array
    
    """
    _arespy_pkg.f90wrap_store_rho(n_rho=n_rho, nspin=nspin, rho=rho)

def store_rho_at(n_rho, nspin, rho_in):
    """
    store_rho_at(n_rho, nspin, rho_in)
    
    
    Defined at Succeed_module.fpp lines 88-93
    
    Parameters
    ----------
    n_rho : int
    nspin : int
    rho_in : float array
    
    """
    _arespy_pkg.f90wrap_store_rho_at(n_rho=n_rho, nspin=nspin, rho_in=rho_in)

def store_r(nr, r):
    """
    store_r(nr, r)
    
    
    Defined at Succeed_module.fpp lines 95-104
    
    Parameters
    ----------
    nr : int
    r : float array
    
    """
    _arespy_pkg.f90wrap_store_r(nr=nr, r=r)

def store_psi(n_rho, n_s, nspin, psi):
    """
    store_psi(n_rho, n_s, nspin, psi)
    
    
    Defined at Succeed_module.fpp lines 106-117
    
    Parameters
    ----------
    n_rho : int
    n_s : int
    nspin : int
    psi : float array
    
    """
    _arespy_pkg.f90wrap_store_psi(n_rho=n_rho, n_s=n_s, nspin=nspin, psi=psi)

def get_rho(nr, r_new, nrho, nspin, rho_new):
    """
    get_rho(nr, r_new, nrho, nspin, rho_new)
    
    
    Defined at Succeed_module.fpp lines 119-199
    
    Parameters
    ----------
    nr : int
    r_new : float array
    nrho : int
    nspin : int
    rho_new : float array
    
    """
    _arespy_pkg.f90wrap_get_rho(nr=nr, r_new=r_new, nrho=nrho, nspin=nspin, \
        rho_new=rho_new)

def get_psi(nrho, n_s, nspin, psi_new):
    """
    get_psi(nrho, n_s, nspin, psi_new)
    
    
    Defined at Succeed_module.fpp lines 201-222
    
    Parameters
    ----------
    nrho : int
    n_s : int
    nspin : int
    psi_new : float array
    
    """
    _arespy_pkg.f90wrap_get_psi(nrho=nrho, n_s=n_s, nspin=nspin, psi_new=psi_new)

def cal_trans_phase(nr, nspin, r_new, n1, n2, n3, ng1, ng2, ng3, gvec, \
    trans_phase):
    """
    cal_trans_phase(nr, nspin, r_new, n1, n2, n3, ng1, ng2, ng3, gvec, trans_phase)
    
    
    Defined at Succeed_module.fpp lines 232-284
    
    Parameters
    ----------
    nr : int
    nspin : int
    r_new : float array
    n1 : int
    n2 : int
    n3 : int
    ng1 : int
    ng2 : int
    ng3 : int
    gvec : float array
    trans_phase : complex array
    
    """
    _arespy_pkg.f90wrap_cal_trans_phase(nr=nr, nspin=nspin, r_new=r_new, n1=n1, \
        n2=n2, n3=n3, ng1=ng1, ng2=ng2, ng3=ng3, gvec=gvec, trans_phase=trans_phase)

def get_new_rho_psi(nr, r_new, nrho, n1, n2, n3, nspin, rho_new, n_s, psi_new, \
    gvec):
    """
    get_new_rho_psi(nr, r_new, nrho, n1, n2, n3, nspin, rho_new, n_s, psi_new, gvec)
    
    
    Defined at Succeed_module.fpp lines 286-345
    
    Parameters
    ----------
    nr : int
    r_new : float array
    nrho : int
    n1 : int
    n2 : int
    n3 : int
    nspin : int
    rho_new : float array
    n_s : int
    psi_new : float array
    gvec : float array
    
    """
    _arespy_pkg.f90wrap_get_new_rho_psi(nr=nr, r_new=r_new, nrho=nrho, n1=n1, n2=n2, \
        n3=n3, nspin=nspin, rho_new=rho_new, n_s=n_s, psi_new=psi_new, gvec=gvec)

def store_rho_fft_trans(n_rho, nspin, rho):
    """
    store_rho_fft_trans(n_rho, nspin, rho)
    
    
    Defined at Succeed_module.fpp lines 347-355
    
    Parameters
    ----------
    n_rho : int
    nspin : int
    rho : float array
    
    """
    _arespy_pkg.f90wrap_store_rho_fft_trans(n_rho=n_rho, nspin=nspin, rho=rho)

def store_rho_at_fft_trans(n_rho, nspin, na, rho_in, rho_in2):
    """
    store_rho_at_fft_trans(n_rho, nspin, na, rho_in, rho_in2)
    
    
    Defined at Succeed_module.fpp lines 357-372
    
    Parameters
    ----------
    n_rho : int
    nspin : int
    na : int
    rho_in : float array
    rho_in2 : float array
    
    """
    _arespy_pkg.f90wrap_store_rho_at_fft_trans(n_rho=n_rho, nspin=nspin, na=na, \
        rho_in=rho_in, rho_in2=rho_in2)

def store_r_fft_trans(nr, r):
    """
    store_r_fft_trans(nr, r)
    
    
    Defined at Succeed_module.fpp lines 374-382
    
    Parameters
    ----------
    nr : int
    r : float array
    
    """
    _arespy_pkg.f90wrap_store_r_fft_trans(nr=nr, r=r)

def store_psi_fft_trans(n_rho, n_s, nspin, psi):
    """
    store_psi_fft_trans(n_rho, n_s, nspin, psi)
    
    
    Defined at Succeed_module.fpp lines 384-396
    
    Parameters
    ----------
    n_rho : int
    n_s : int
    nspin : int
    psi : float array
    
    """
    _arespy_pkg.f90wrap_store_psi_fft_trans(n_rho=n_rho, n_s=n_s, nspin=nspin, \
        psi=psi)

def get_llastrho():
    """
    Element llastrho ftype=logical pytype=bool
    
    
    Defined at Succeed_module.fpp line 12
    
    """
    return _arespy_pkg.f90wrap_succeed__get__llastrho()

def set_llastrho(llastrho):
    _arespy_pkg.f90wrap_succeed__set__llastrho(llastrho)

def get_lsr():
    """
    Element lsr ftype=logical pytype=bool
    
    
    Defined at Succeed_module.fpp line 12
    
    """
    return _arespy_pkg.f90wrap_succeed__get__lsr()

def set_lsr(lsr):
    _arespy_pkg.f90wrap_succeed__set__lsr(lsr)

def get_lsrho():
    """
    Element lsrho ftype=logical pytype=bool
    
    
    Defined at Succeed_module.fpp line 12
    
    """
    return _arespy_pkg.f90wrap_succeed__get__lsrho()

def set_lsrho(lsrho):
    _arespy_pkg.f90wrap_succeed__set__lsrho(lsrho)

def get_lspsi():
    """
    Element lspsi ftype=logical pytype=bool
    
    
    Defined at Succeed_module.fpp line 12
    
    """
    return _arespy_pkg.f90wrap_succeed__get__lspsi()

def set_lspsi(lspsi):
    _arespy_pkg.f90wrap_succeed__set__lspsi(lspsi)

def get_lsrho_at():
    """
    Element lsrho_at ftype=logical pytype=bool
    
    
    Defined at Succeed_module.fpp line 12
    
    """
    return _arespy_pkg.f90wrap_succeed__get__lsrho_at()

def set_lsrho_at(lsrho_at):
    _arespy_pkg.f90wrap_succeed__set__lsrho_at(lsrho_at)

def get_array_rho1():
    """
    Element rho1 ftype=real(dp) pytype=float
    
    
    Defined at Succeed_module.fpp line 17
    
    """
    global rho1
    array_ndim, array_type, array_shape, array_handle = \
        _arespy_pkg.f90wrap_succeed__array__rho1(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        rho1 = _arrays[array_handle]
    else:
        rho1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_succeed__array__rho1)
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
        _arespy_pkg.f90wrap_succeed__array__rho2(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        rho2 = _arrays[array_handle]
    else:
        rho2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_succeed__array__rho2)
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
        _arespy_pkg.f90wrap_succeed__array__rho3(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        rho3 = _arrays[array_handle]
    else:
        rho3 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_succeed__array__rho3)
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
        _arespy_pkg.f90wrap_succeed__array__r1(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        r1 = _arrays[array_handle]
    else:
        r1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_succeed__array__r1)
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
        _arespy_pkg.f90wrap_succeed__array__r2(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        r2 = _arrays[array_handle]
    else:
        r2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_succeed__array__r2)
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
        _arespy_pkg.f90wrap_succeed__array__r3(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        r3 = _arrays[array_handle]
    else:
        r3 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_succeed__array__r3)
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
        _arespy_pkg.f90wrap_succeed__array__psi1(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        psi1 = _arrays[array_handle]
    else:
        psi1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_succeed__array__psi1)
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
        _arespy_pkg.f90wrap_succeed__array__psi2(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        psi2 = _arrays[array_handle]
    else:
        psi2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_succeed__array__psi2)
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
        _arespy_pkg.f90wrap_succeed__array__psi3(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        psi3 = _arrays[array_handle]
    else:
        psi3 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_succeed__array__psi3)
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
        _arespy_pkg.f90wrap_succeed__array__rho_at(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        rho_at = _arrays[array_handle]
    else:
        rho_at = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_succeed__array__rho_at)
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
        _arespy_pkg.f90wrap_succeed__array__rhoi(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        rhoi = _arrays[array_handle]
    else:
        rhoi = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_succeed__array__rhoi)
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
        _arespy_pkg.f90wrap_succeed__array__rho_at1(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        rho_at1 = _arrays[array_handle]
    else:
        rho_at1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_succeed__array__rho_at1)
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
        _arespy_pkg.f90wrap_succeed__array__rhoi1(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        rhoi1 = _arrays[array_handle]
    else:
        rhoi1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_succeed__array__rhoi1)
        _arrays[array_handle] = rhoi1
    return rhoi1

def set_array_rhoi1(rhoi1):
    rhoi1[...] = rhoi1

def get_dvol():
    """
    Element dvol ftype=real(dp) pytype=float
    
    
    Defined at Succeed_module.fpp line 19
    
    """
    return _arespy_pkg.f90wrap_succeed__get__dvol()

def set_dvol(dvol):
    _arespy_pkg.f90wrap_succeed__set__dvol(dvol)

def get_counter1():
    """
    Element counter1 ftype=integer(i4b) pytype=int
    
    
    Defined at Succeed_module.fpp line 20
    
    """
    return _arespy_pkg.f90wrap_succeed__get__counter1()

def set_counter1(counter1):
    _arespy_pkg.f90wrap_succeed__set__counter1(counter1)

def get_counter2():
    """
    Element counter2 ftype=integer(i4b) pytype=int
    
    
    Defined at Succeed_module.fpp line 20
    
    """
    return _arespy_pkg.f90wrap_succeed__get__counter2()

def set_counter2(counter2):
    _arespy_pkg.f90wrap_succeed__set__counter2(counter2)

def get_counter3():
    """
    Element counter3 ftype=integer(i4b) pytype=int
    
    
    Defined at Succeed_module.fpp line 20
    
    """
    return _arespy_pkg.f90wrap_succeed__get__counter3()

def set_counter3(counter3):
    _arespy_pkg.f90wrap_succeed__set__counter3(counter3)

def get_alpha():
    """
    Element alpha ftype=real(dp) pytype=float
    
    
    Defined at Succeed_module.fpp line 21
    
    """
    return _arespy_pkg.f90wrap_succeed__get__alpha()

def set_alpha(alpha):
    _arespy_pkg.f90wrap_succeed__set__alpha(alpha)

def get_beta():
    """
    Element beta ftype=real(dp) pytype=float
    
    
    Defined at Succeed_module.fpp line 21
    
    """
    return _arespy_pkg.f90wrap_succeed__get__beta()

def set_beta(beta):
    _arespy_pkg.f90wrap_succeed__set__beta(beta)


_array_initialisers = [get_array_rho1, get_array_rho2, get_array_rho3, \
    get_array_r1, get_array_r2, get_array_r3, get_array_psi1, get_array_psi2, \
    get_array_psi3, get_array_rho_at, get_array_rhoi, get_array_rho_at1, \
    get_array_rhoi1]
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "succeed".')

for func in _dt_array_initialisers:
    func()

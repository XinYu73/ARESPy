"""
Module constants


Defined at Constants.fpp lines 5-97

"""
from __future__ import print_function, absolute_import, division
import _arespy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def get_i4b():
    """
    Element i4b ftype=integer pytype=int
    
    
    Defined at Constants.fpp line 25
    
    """
    return _arespy_pkg.f90wrap_constants__get__i4b()

I4B = get_i4b()

def get_i2b():
    """
    Element i2b ftype=integer pytype=int
    
    
    Defined at Constants.fpp line 26
    
    """
    return _arespy_pkg.f90wrap_constants__get__i2b()

I2B = get_i2b()

def get_i1b():
    """
    Element i1b ftype=integer pytype=int
    
    
    Defined at Constants.fpp line 27
    
    """
    return _arespy_pkg.f90wrap_constants__get__i1b()

I1B = get_i1b()

def get_sp():
    """
    Element sp ftype=integer pytype=int
    
    
    Defined at Constants.fpp line 28
    
    """
    return _arespy_pkg.f90wrap_constants__get__sp()

SP = get_sp()

def get_dp():
    """
    Element dp ftype=integer pytype=int
    
    
    Defined at Constants.fpp line 29
    
    """
    return _arespy_pkg.f90wrap_constants__get__dp()

DP = get_dp()

def get_scp():
    """
    Element scp ftype=integer pytype=int
    
    
    Defined at Constants.fpp line 30
    
    """
    return _arespy_pkg.f90wrap_constants__get__scp()

SCP = get_scp()

def get_dcp():
    """
    Element dcp ftype=integer pytype=int
    
    
    Defined at Constants.fpp line 31
    
    """
    return _arespy_pkg.f90wrap_constants__get__dcp()

DCP = get_dcp()

def get_lgt():
    """
    Element lgt ftype=integer pytype=int
    
    
    Defined at Constants.fpp line 33
    
    """
    return _arespy_pkg.f90wrap_constants__get__lgt()

LGT = get_lgt()

def get_pi():
    """
    Element pi ftype=real(dp) pytype=float
    
    
    Defined at Constants.fpp line 36
    
    """
    return _arespy_pkg.f90wrap_constants__get__pi()

pi = get_pi()

def get_pio2():
    """
    Element pio2 ftype=real(dp) pytype=float
    
    
    Defined at Constants.fpp line 37
    
    """
    return _arespy_pkg.f90wrap_constants__get__pio2()

pio2 = get_pio2()

def get_twopi():
    """
    Element twopi ftype=real(dp) pytype=float
    
    
    Defined at Constants.fpp line 38
    
    """
    return _arespy_pkg.f90wrap_constants__get__twopi()

twopi = get_twopi()

def get_imag():
    """
    Element imag ftype=complex(dcp) pytype=complex
    
    
    Defined at Constants.fpp line 39
    
    """
    return _arespy_pkg.f90wrap_constants__get__imag()

imag = get_imag()

def get_inputunit():
    """
    Element inputunit ftype=integer(i4b) pytype=int
    
    
    Defined at Constants.fpp line 40
    
    """
    return _arespy_pkg.f90wrap_constants__get__inputunit()

inputunit = get_inputunit()

def get_errunit():
    """
    Element errunit ftype=integer(i4b) pytype=int
    
    
    Defined at Constants.fpp line 41
    
    """
    return _arespy_pkg.f90wrap_constants__get__errunit()

errunit = get_errunit()

def get_outputunit():
    """
    Element outputunit ftype=integer(i4b) pytype=int
    
    
    Defined at Constants.fpp line 42
    
    """
    return _arespy_pkg.f90wrap_constants__get__outputunit()

outputunit = get_outputunit()

def get_mdposunit():
    """
    Element mdposunit ftype=integer(i4b) pytype=int
    
    
    Defined at Constants.fpp line 43
    
    """
    return _arespy_pkg.f90wrap_constants__get__mdposunit()

mdposunit = get_mdposunit()

def get_rydberg():
    """
    Element rydberg ftype=real(dp) pytype=float
    
    
    Defined at Constants.fpp line 44
    
    """
    return _arespy_pkg.f90wrap_constants__get__rydberg()

rydberg = get_rydberg()

def get_bohr2ang():
    """
    Element bohr2ang ftype=real(dp) pytype=float
    
    
    Defined at Constants.fpp line 45
    
    """
    return _arespy_pkg.f90wrap_constants__get__bohr2ang()

bohr2ang = get_bohr2ang()

def get_ang2bohr():
    """
    Element ang2bohr ftype=real(dp) pytype=float
    
    
    Defined at Constants.fpp line 46
    
    """
    return _arespy_pkg.f90wrap_constants__get__ang2bohr()

ang2bohr = get_ang2bohr()

def get_hart2ev():
    """
    Element hart2ev ftype=real(dp) pytype=float
    
    
    Defined at Constants.fpp line 47
    
    """
    return _arespy_pkg.f90wrap_constants__get__hart2ev()

hart2ev = get_hart2ev()

def get_ev2hart():
    """
    Element ev2hart ftype=real(dp) pytype=float
    
    
    Defined at Constants.fpp line 48
    
    """
    return _arespy_pkg.f90wrap_constants__get__ev2hart()

ev2hart = get_ev2hart()

def get_scf_tol():
    """
    Element scf_tol ftype=real(dp) pytype=float
    
    
    Defined at Constants.fpp line 49
    
    """
    return _arespy_pkg.f90wrap_constants__get__scf_tol()

scf_tol = get_scf_tol()

def get_au2ev_force():
    """
    Element au2ev_force ftype=real(dp) pytype=float
    
    
    Defined at Constants.fpp line 50
    
    """
    return _arespy_pkg.f90wrap_constants__get__au2ev_force()

au2ev_force = get_au2ev_force()

def get_au2gpa():
    """
    Element au2gpa ftype=real(dp) pytype=float
    
    
    Defined at Constants.fpp line 51
    
    """
    return _arespy_pkg.f90wrap_constants__get__au2gpa()

au2gpa = get_au2gpa()

def get_golden():
    """
    Element golden ftype=real(dp) pytype=float
    
    
    Defined at Constants.fpp line 52
    
    """
    return _arespy_pkg.f90wrap_constants__get__golden()

golden = get_golden()

def get_vlight():
    """
    Element vlight ftype=real(dp) pytype=float
    
    
    Defined at Constants.fpp line 53
    
    """
    return _arespy_pkg.f90wrap_constants__get__vlight()

vlight = get_vlight()

def get_const_me_au():
    """
    Element const_me_au ftype=real(dp) pytype=float
    
    
    Defined at Constants.fpp line 57
    
    """
    return _arespy_pkg.f90wrap_constants__get__const_me_au()

CONST_ME_AU = get_const_me_au()

def get_const_e_au():
    """
    Element const_e_au ftype=real(dp) pytype=float
    
    
    Defined at Constants.fpp line 59
    
    """
    return _arespy_pkg.f90wrap_constants__get__const_e_au()

CONST_E_AU = get_const_e_au()

def get_const_eh_au():
    """
    Element const_eh_au ftype=real(dp) pytype=float
    
    
    Defined at Constants.fpp line 61
    
    """
    return _arespy_pkg.f90wrap_constants__get__const_eh_au()

CONST_EH_AU = get_const_eh_au()

def get_const_len_au():
    """
    Element const_len_au ftype=real(dp) pytype=float
    
    
    Defined at Constants.fpp line 63
    
    """
    return _arespy_pkg.f90wrap_constants__get__const_len_au()

CONST_LEN_AU = get_const_len_au()

def get_const_hbar_au():
    """
    Element const_hbar_au ftype=real(dp) pytype=float
    
    
    Defined at Constants.fpp line 65
    
    """
    return _arespy_pkg.f90wrap_constants__get__const_hbar_au()

CONST_HBAR_AU = get_const_hbar_au()

def get_const_ep_au():
    """
    Element const_ep_au ftype=real(dp) pytype=float
    
    
    Defined at Constants.fpp line 67
    
    """
    return _arespy_pkg.f90wrap_constants__get__const_ep_au()

CONST_EP_AU = get_const_ep_au()

def get_const_f_au():
    """
    Element const_f_au ftype=real(dp) pytype=float
    
    
    Defined at Constants.fpp line 69
    
    """
    return _arespy_pkg.f90wrap_constants__get__const_f_au()

CONST_F_AU = get_const_f_au()

def get_const_mt_au():
    """
    Element const_mt_au ftype=real(dp) pytype=float
    
    
    Defined at Constants.fpp line 71
    
    """
    return _arespy_pkg.f90wrap_constants__get__const_mt_au()

CONST_MT_AU = get_const_mt_au()

def get_const_time_au():
    """
    Element const_time_au ftype=real(dp) pytype=float
    
    
    Defined at Constants.fpp line 73
    
    """
    return _arespy_pkg.f90wrap_constants__get__const_time_au()

CONST_TIME_AU = get_const_time_au()

def get_const_i_au():
    """
    Element const_i_au ftype=real(dp) pytype=float
    
    
    Defined at Constants.fpp line 75
    
    """
    return _arespy_pkg.f90wrap_constants__get__const_i_au()

CONST_I_AU = get_const_i_au()

def get_const_temp_au():
    """
    Element const_temp_au ftype=real(dp) pytype=float
    
    
    Defined at Constants.fpp line 77
    
    """
    return _arespy_pkg.f90wrap_constants__get__const_temp_au()

CONST_TEMP_AU = get_const_temp_au()

def get_const_p_au():
    """
    Element const_p_au ftype=real(dp) pytype=float
    
    
    Defined at Constants.fpp line 79
    
    """
    return _arespy_pkg.f90wrap_constants__get__const_p_au()

CONST_P_AU = get_const_p_au()

def get_const_v_au():
    """
    Element const_v_au ftype=real(dp) pytype=float
    
    
    Defined at Constants.fpp line 81
    
    """
    return _arespy_pkg.f90wrap_constants__get__const_v_au()

CONST_V_AU = get_const_v_au()

def get_const_ke_au():
    """
    Element const_ke_au ftype=real(dp) pytype=float
    
    
    Defined at Constants.fpp line 83
    
    """
    return _arespy_pkg.f90wrap_constants__get__const_ke_au()

CONST_KE_AU = get_const_ke_au()

def get_const_mu_au():
    """
    Element const_mu_au ftype=real(dp) pytype=float
    
    
    Defined at Constants.fpp line 85
    
    """
    return _arespy_pkg.f90wrap_constants__get__const_mu_au()

CONST_MU_AU = get_const_mu_au()

def get_const_ma_au():
    """
    Element const_ma_au ftype=real(dp) pytype=float
    
    
    Defined at Constants.fpp line 87
    
    """
    return _arespy_pkg.f90wrap_constants__get__const_ma_au()

CONST_MA_AU = get_const_ma_au()

def get_const_kb_si():
    """
    Element const_kb_si ftype=real(dp) pytype=float
    
    
    Defined at Constants.fpp line 89
    
    """
    return _arespy_pkg.f90wrap_constants__get__const_kb_si()

CONST_KB_SI = get_const_kb_si()

def get_angs():
    """
    Element angs ftype=real(dp) pytype=float
    
    
    Defined at Constants.fpp line 92
    
    """
    return _arespy_pkg.f90wrap_constants__get__angs()

angs = get_angs()

def get_force2ev():
    """
    Element force2ev ftype=real(dp) pytype=float
    
    
    Defined at Constants.fpp line 93
    
    """
    return _arespy_pkg.f90wrap_constants__get__force2ev()

force2ev = get_force2ev()

def get_clen():
    """
    Element clen ftype=integer pytype=int
    
    
    Defined at Constants.fpp line 95
    
    """
    return _arespy_pkg.f90wrap_constants__get__clen()

CLen = get_clen()

def get_xtiny():
    """
    Element xtiny ftype=real(dp) pytype=float
    
    
    Defined at Constants.fpp line 96
    
    """
    return _arespy_pkg.f90wrap_constants__get__xtiny()

xtiny = get_xtiny()

def get_autogpa():
    """
    Element autogpa ftype=real(dp) pytype=float
    
    
    Defined at Constants.fpp line 97
    
    """
    return _arespy_pkg.f90wrap_constants__get__autogpa()

autogpa = get_autogpa()

def get_hartree2ev():
    """
    Element hartree2ev ftype=real(dp) pytype=float
    
    
    Defined at Constants.fpp line 98
    
    """
    return _arespy_pkg.f90wrap_constants__get__hartree2ev()

hartree2ev = get_hartree2ev()


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "constants".')

for func in _dt_array_initialisers:
    func()

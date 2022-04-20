"""
Module scf_module


Defined at Scf_module.fpp lines 5-1462

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def arpackscf(rhos, psi, eval):
    """
    arpackscf(rhos, psi, eval)
    
    
    Defined at Scf_module.fpp lines 21-43
    
    Parameters
    ----------
    rhos : float array
    psi : complex array
    eval : float array
    
    """
    _AresMainPy_pkg.f90wrap_arpackscf(rhos=rhos, psi=psi, eval=eval)

def solver_spin(rhos, psi, eval, nev, diagtol):
    """
    solver_spin(rhos, psi, eval, nev, diagtol)
    
    
    Defined at Scf_module.fpp lines 46-66
    
    Parameters
    ----------
    rhos : float array
    psi : complex array
    eval : float array
    nev : int
    diagtol : float
    
    """
    _AresMainPy_pkg.f90wrap_solver_spin(rhos=rhos, psi=psi, eval=eval, nev=nev, \
        diagtol=diagtol)

def ksolver(veff, kpsi, keval, nev, tol):
    """
    ksolver(veff, kpsi, keval, nev, tol)
    
    
    Defined at Scf_module.fpp lines 69-113
    
    Parameters
    ----------
    veff : float array
    kpsi : complex array
    keval : float array
    nev : int
    tol : float
    
    """
    _AresMainPy_pkg.f90wrap_ksolver(veff=veff, kpsi=kpsi, keval=keval, nev=nev, \
        tol=tol)

def smear_updaterho(nev, ne, psi, eval, wke_l, rhos):
    """
    fme, ets = smear_updaterho(nev, ne, psi, eval, wke_l, rhos)
    
    
    Defined at Scf_module.fpp lines 117-174
    
    Parameters
    ----------
    nev : int
    ne : float
    psi : complex array
    eval : float array
    wke_l : float array
    rhos : float array
    
    Returns
    -------
    fme : float
    ets : float
    
    """
    fme, ets = _AresMainPy_pkg.f90wrap_smear_updaterho(nev=nev, ne=ne, psi=psi, \
        eval=eval, wke_l=wke_l, rhos=rhos)
    return fme, ets

def chefsi(rhos, psi, eval):
    """
    chefsi(rhos, psi, eval)
    
    
    Defined at Scf_module.fpp lines 182-386
    
    Parameters
    ----------
    rhos : float array
    psi : complex array
    eval : float array
    
    ===============================================================
     WRITE(6,"(1X,A8,I4,1X,A6,E15.7,1X,A4,E15.7)") &
          &      '>CheFSI:',iter,'dTOTEN',dtoten,      &
          &  'dRHO',drho
     IF(drho<RTOL.AND.ABS(dtoten)<ETOL) EXIT
    ===============================================================
    """
    _AresMainPy_pkg.f90wrap_chefsi(rhos=rhos, psi=psi, eval=eval)

def filter_spin(veff, psi, eval):
    """
    filter_spin(veff, psi, eval)
    
    
    Defined at Scf_module.fpp lines 389-411
    
    Parameters
    ----------
    veff : float array
    psi : complex array
    eval : float array
    
    """
    _AresMainPy_pkg.f90wrap_filter_spin(veff=veff, psi=psi, eval=eval)

def arpackscf_r(rhos, psi, eval):
    """
    arpackscf_r(rhos, psi, eval)
    
    
    Defined at Scf_module.fpp lines 435-502
    
    Parameters
    ----------
    rhos : float array
    psi : float array
    eval : float array
    
    ===============================================================
    """
    _AresMainPy_pkg.f90wrap_arpackscf_r(rhos=rhos, psi=psi, eval=eval)

def solver_spin_r(rhos, psi, eval, nev, diagtol):
    """
    solver_spin_r(rhos, psi, eval, nev, diagtol)
    
    
    Defined at Scf_module.fpp lines 505-558
    
    Parameters
    ----------
    rhos : float array
    psi : float array
    eval : float array
    nev : int
    diagtol : float
    
    """
    _AresMainPy_pkg.f90wrap_solver_spin_r(rhos=rhos, psi=psi, eval=eval, nev=nev, \
        diagtol=diagtol)

def real_smear_updaterho(nev, ne, psi, eval, wke_l, rhos):
    """
    fme, ets = real_smear_updaterho(nev, ne, psi, eval, wke_l, rhos)
    
    
    Defined at Scf_module.fpp lines 561-606
    
    Parameters
    ----------
    nev : int
    ne : float
    psi : float array
    eval : float array
    wke_l : float array
    rhos : float array
    
    Returns
    -------
    fme : float
    ets : float
    
    """
    fme, ets = _AresMainPy_pkg.f90wrap_real_smear_updaterho(nev=nev, ne=ne, psi=psi, \
        eval=eval, wke_l=wke_l, rhos=rhos)
    return fme, ets

def bvk_chefsi(rhos, psi, eval):
    """
    bvk_chefsi(rhos, psi, eval)
    
    
    Defined at Scf_module.fpp lines 613-710
    
    Parameters
    ----------
    rhos : float array
    psi : float array
    eval : float array
    
    ================================================
    ##XLT test I/O FOR DETECT ISO
    open(unit=1118,file="RHO_ISO")
    READ(1118,*)
    READ(1118,*)rhoS
    close(1118)
    goto 1188
    ================================================
    initial mixer
    """
    _AresMainPy_pkg.f90wrap_bvk_chefsi(rhos=rhos, psi=psi, eval=eval)

def bvk_filter_spin(rhos, x, d):
    """
    bvk_filter_spin(rhos, x, d)
    
    
    Defined at Scf_module.fpp lines 713-732
    
    Parameters
    ----------
    rhos : float array
    x : float array
    d : float array
    
    """
    _AresMainPy_pkg.f90wrap_bvk_filter_spin(rhos=rhos, x=x, d=d)

def prr_filter_spin(nfs, nfe, rhos, x, efr, pbar):
    """
    fme, ets = prr_filter_spin(nfs, nfe, rhos, x, efr, pbar)
    
    
    Defined at Scf_module.fpp lines 739-827
    
    Parameters
    ----------
    nfs : int
    nfe : float
    rhos : float array
    x : float array
    efr : float array
    pbar : float array
    
    Returns
    -------
    fme : float
    ets : float
    
    """
    fme, ets = _AresMainPy_pkg.f90wrap_prr_filter_spin(nfs=nfs, nfe=nfe, rhos=rhos, \
        x=x, efr=efr, pbar=pbar)
    return fme, ets

def prr_chefsi(rhos, phi):
    """
    prr_chefsi(rhos, phi)
    
    
    Defined at Scf_module.fpp lines 830-908
    
    Parameters
    ----------
    rhos : float array
    phi : float array
    
    -------------------------------------------------------
    """
    _AresMainPy_pkg.f90wrap_prr_chefsi(rhos=rhos, phi=phi)

def iso_solver_spin_r(rhos, psi, eval, nev, diagtol):
    """
    iso_solver_spin_r(rhos, psi, eval, nev, diagtol)
    
    
    Defined at Scf_module.fpp lines 918-971
    
    Parameters
    ----------
    rhos : float array
    psi : float array
    eval : float array
    nev : int
    diagtol : float
    
    """
    _AresMainPy_pkg.f90wrap_iso_solver_spin_r(rhos=rhos, psi=psi, eval=eval, \
        nev=nev, diagtol=diagtol)

def iso_smear_updaterho(nev, ne, psi, eval, wke_l, rhos_out):
    """
    fme, ets = iso_smear_updaterho(nev, ne, psi, eval, wke_l, rhos_out)
    
    
    Defined at Scf_module.fpp lines 974-1023
    
    Parameters
    ----------
    nev : int
    ne : float
    psi : float array
    eval : float array
    wke_l : float array
    rhos_out : float array
    
    Returns
    -------
    fme : float
    ets : float
    
    """
    fme, ets = _AresMainPy_pkg.f90wrap_iso_smear_updaterho(nev=nev, ne=ne, psi=psi, \
        eval=eval, wke_l=wke_l, rhos_out=rhos_out)
    return fme, ets

def iso_chefsi(rhos, psi, eval):
    """
    iso_chefsi(rhos, psi, eval)
    
    
    Defined at Scf_module.fpp lines 1030-1204
    
    Parameters
    ----------
    rhos : float array
    psi : float array
    eval : float array
    
    ================================================
    xlt 2018-5-18 USED FOR TEST
     USE grid_module , ONLY : ISO_Rho2grid,rho_calc
    ================================================
    """
    _AresMainPy_pkg.f90wrap_iso_chefsi(rhos=rhos, psi=psi, eval=eval)

def iso_filter_spin(rhos, x, d):
    """
    iso_filter_spin(rhos, x, d)
    
    
    Defined at Scf_module.fpp lines 1207-1238
    
    Parameters
    ----------
    rhos : float array
    x : float array
    d : float array
    
    =============================================
    print *,'calculate effectial time -->',(tf-te)/10000.d0
    print *,'filter GeneralRealyRitz  time -->',(tg-tf)/10000.d0
    print *,'total chebyshev filter time -->',(tg-te)/10000.d0
    =============================================
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    """
    _AresMainPy_pkg.f90wrap_iso_filter_spin(rhos=rhos, x=x, d=d)

def chefsi_gamma(rhos, psi, eval):
    """
    chefsi_gamma(rhos, psi, eval)
    
    
    Defined at Scf_module.fpp lines 1241-1378
    
    Parameters
    ----------
    rhos : float array
    psi : float array
    eval : float array
    
    """
    _AresMainPy_pkg.f90wrap_chefsi_gamma(rhos=rhos, psi=psi, eval=eval)

def smear_updaterho_gamma(nev, ne, psi, eval, wke_l, rhos):
    """
    fme, ets = smear_updaterho_gamma(nev, ne, psi, eval, wke_l, rhos)
    
    
    Defined at Scf_module.fpp lines 1381-1438
    
    Parameters
    ----------
    nev : int
    ne : float
    psi : float array
    eval : float array
    wke_l : float array
    rhos : float array
    
    Returns
    -------
    fme : float
    ets : float
    
    """
    fme, ets = _AresMainPy_pkg.f90wrap_smear_updaterho_gamma(nev=nev, ne=ne, \
        psi=psi, eval=eval, wke_l=wke_l, rhos=rhos)
    return fme, ets

def filter_spin_gamma(veff, psi, eval):
    """
    filter_spin_gamma(veff, psi, eval)
    
    
    Defined at Scf_module.fpp lines 1441-1461
    
    Parameters
    ----------
    veff : float array
    psi : float array
    eval : float array
    
    """
    _AresMainPy_pkg.f90wrap_filter_spin_gamma(veff=veff, psi=psi, eval=eval)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "scf_module".')

for func in _dt_array_initialisers:
    func()

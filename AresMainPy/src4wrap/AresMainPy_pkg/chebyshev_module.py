"""
Module chebyshev_module


Defined at Chebyshev_fliter.fpp lines 10-2437

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def cheby_filter_rr(ik, veff, x, d):
    """
    cheby_filter_rr(ik, veff, x, d)
    
    
    Defined at Chebyshev_fliter.fpp lines 24-66
    
    Parameters
    ----------
    ik : int
    veff : float array
    x : complex array
    d : float array
    
    """
    _AresMainPy_pkg.f90wrap_cheby_filter_rr(ik=ik, veff=veff, x=x, d=d)

def grayleigh_ritz(ik, veff, x, d):
    """
    grayleigh_ritz(ik, veff, x, d)
    
    
    Defined at Chebyshev_fliter.fpp lines 69-127
    
    Parameters
    ----------
    ik : int
    veff : float array
    x : complex array
    d : float array
    
    """
    _AresMainPy_pkg.f90wrap_grayleigh_ritz(ik=ik, veff=veff, x=x, d=d)

def chebyshev_filter(ik, veff, x, m, a, b):
    """
    chebyshev_filter(ik, veff, x, m, a, b)
    
    
    Defined at Chebyshev_fliter.fpp lines 130-164
    
    Parameters
    ----------
    ik : int
    veff : float array
    x : complex array
    m : int
    a : float
    b : float
    
    """
    _AresMainPy_pkg.f90wrap_chebyshev_filter(ik=ik, veff=veff, x=x, m=m, a=a, b=b)

def chebyshev_filter_scaled(ik, veff, x, m, a, b, al):
    """
    chebyshev_filter_scaled(ik, veff, x, m, a, b, al)
    
    
    Defined at Chebyshev_fliter.fpp lines 167-211
    
    Parameters
    ----------
    ik : int
    veff : float array
    x : complex array
    m : int
    a : float
    b : float
    al : float
    
    """
    _AresMainPy_pkg.f90wrap_chebyshev_filter_scaled(ik=ik, veff=veff, x=x, m=m, a=a, \
        b=b, al=al)

def cal_hx(ik, veff, nst, v, hv):
    """
    cal_hx(ik, veff, nst, v, hv)
    
    
    Defined at Chebyshev_fliter.fpp lines 214-231
    
    Parameters
    ----------
    ik : int
    veff : float array
    nst : int
    v : complex array
    hv : complex array
    
    """
    _AresMainPy_pkg.f90wrap_cal_hx(ik=ik, veff=veff, nst=nst, v=v, hv=hv)

def rayleigh_quotient(ik, veff, nst, x, xhx):
    """
    rayleigh_quotient(ik, veff, nst, x, xhx)
    
    
    Defined at Chebyshev_fliter.fpp lines 234-281
    
    Parameters
    ----------
    ik : int
    veff : float array
    nst : int
    x : complex array
    xhx : complex array
    
    """
    _AresMainPy_pkg.f90wrap_rayleigh_quotient(ik=ik, veff=veff, nst=nst, x=x, \
        xhx=xhx)

def rayleigh_ritz(ik, veff, x, d):
    """
    rayleigh_ritz(ik, veff, x, d)
    
    
    Defined at Chebyshev_fliter.fpp lines 284-315
    
    Parameters
    ----------
    ik : int
    veff : float array
    x : complex array
    d : float array
    
    """
    _AresMainPy_pkg.f90wrap_rayleigh_ritz(ik=ik, veff=veff, x=x, d=d)

def estupb(k, ik, veff, vec):
    """
    b = estupb(k, ik, veff, vec)
    
    
    Defined at Chebyshev_fliter.fpp lines 318-377
    
    Parameters
    ----------
    k : int
    ik : int
    veff : float array
    vec : complex array
    
    Returns
    -------
    b : float
    
    """
    b = _AresMainPy_pkg.f90wrap_estupb(k=k, ik=ik, veff=veff, vec=vec)
    return b

def inituplow(k, ik, veff, v, eval):
    """
    a, b, al = inituplow(k, ik, veff, v, eval)
    
    
    Defined at Chebyshev_fliter.fpp lines 384-427
    
    Parameters
    ----------
    k : int
    ik : int
    veff : float array
    v : complex array
    eval : float array
    
    Returns
    -------
    a : float
    b : float
    al : float
    
    """
    a, b, al = _AresMainPy_pkg.f90wrap_inituplow(k=k, ik=ik, veff=veff, v=v, \
        eval=eval)
    return a, b, al

def first_scfstep_filter(ik, veff, x, eval):
    """
    first_scfstep_filter(ik, veff, x, eval)
    
    
    Defined at Chebyshev_fliter.fpp lines 430-484
    
    Parameters
    ----------
    ik : int
    veff : float array
    x : complex array
    eval : float array
    
    -----------------
    """
    _AresMainPy_pkg.f90wrap_first_scfstep_filter(ik=ik, veff=veff, x=x, eval=eval)

def randomfast(psi, nev, radius):
    """
    randomfast(psi, nev, radius)
    
    
    Defined at Chebyshev_fliter.fpp lines 487-511
    
    Parameters
    ----------
    psi : complex array
    nev : int
    radius : float
    
    """
    _AresMainPy_pkg.f90wrap_randomfast(psi=psi, nev=nev, radius=radius)

def first_chescf_random(rhos, nev, psi, eval):
    """
    first_chescf_random(rhos, nev, psi, eval)
    
    
    Defined at Chebyshev_fliter.fpp lines 514-537
    
    Parameters
    ----------
    rhos : float array
    nev : int
    psi : complex array
    eval : float array
    
    """
    _AresMainPy_pkg.f90wrap_first_chescf_random(rhos=rhos, nev=nev, psi=psi, \
        eval=eval)

def first_chescf(veff, psi_ran, nev, psi, eval):
    """
    first_chescf(veff, psi_ran, nev, psi, eval)
    
    
    Defined at Chebyshev_fliter.fpp lines 540-559
    
    Parameters
    ----------
    veff : float array
    psi_ran : complex array
    nev : int
    psi : complex array
    eval : float array
    
    """
    _AresMainPy_pkg.f90wrap_first_chescf(veff=veff, psi_ran=psi_ran, nev=nev, \
        psi=psi, eval=eval)

def cheby_init_sto(nmax, npw, initx_sto):
    """
    cheby_init_sto(nmax, npw, initx_sto)
    
    
    Defined at Chebyshev_fliter.fpp lines 565-655
    
    Parameters
    ----------
    nmax : int
    npw : int
    initx_sto : complex array
    
    """
    _AresMainPy_pkg.f90wrap_cheby_init_sto(nmax=nmax, npw=npw, initx_sto=initx_sto)

def first_chescf_sto(veff, nev, psi, eval):
    """
    first_chescf_sto(veff, nev, psi, eval)
    
    
    Defined at Chebyshev_fliter.fpp lines 658-720
    
    Parameters
    ----------
    veff : float array
    nev : int
    psi : complex array
    eval : float array
    
    """
    _AresMainPy_pkg.f90wrap_first_chescf_sto(veff=veff, nev=nev, psi=psi, eval=eval)

def first_subspace_stopw(ik, veff, nmax, nev, initx, x, d):
    """
    first_subspace_stopw(ik, veff, nmax, nev, initx, x, d)
    
    
    Defined at Chebyshev_fliter.fpp lines 723-823
    
    Parameters
    ----------
    ik : int
    veff : float array
    nmax : int
    nev : int
    initx : complex array
    x : complex array
    d : float array
    
    """
    _AresMainPy_pkg.f90wrap_first_subspace_stopw(ik=ik, veff=veff, nmax=nmax, \
        nev=nev, initx=initx, x=x, d=d)

def bvk_first_chescf_sto_rand(rhos, nev, psi, eval):
    """
    bvk_first_chescf_sto_rand(rhos, nev, psi, eval)
    
    
    Defined at Chebyshev_fliter.fpp lines 832-875
    
    Parameters
    ----------
    rhos : float array
    nev : int
    psi : float array
    eval : float array
    
    """
    _AresMainPy_pkg.f90wrap_bvk_first_chescf_sto_rand(rhos=rhos, nev=nev, psi=psi, \
        eval=eval)

def bvk_cheby_init_sto_rand(nmax, nrand, initx_sto):
    """
    bvk_cheby_init_sto_rand(nmax, nrand, initx_sto)
    
    
    Defined at Chebyshev_fliter.fpp lines 878-958
    
    Parameters
    ----------
    nmax : int
    nrand : int
    initx_sto : float array
    
    """
    _AresMainPy_pkg.f90wrap_bvk_cheby_init_sto_rand(nmax=nmax, nrand=nrand, \
        initx_sto=initx_sto)

def first_subspace_sto_rand(veff, nmax, nev, initx, x, d):
    """
    first_subspace_sto_rand(veff, nmax, nev, initx, x, d)
    
    
    Defined at Chebyshev_fliter.fpp lines 961-988
    
    Parameters
    ----------
    veff : float array
    nmax : int
    nev : int
    initx : float array
    x : float array
    d : float array
    
    """
    _AresMainPy_pkg.f90wrap_first_subspace_sto_rand(veff=veff, nmax=nmax, nev=nev, \
        initx=initx, x=x, d=d)

def rayleigh_quotient_real(veff, nst, x, xhx):
    """
    rayleigh_quotient_real(veff, nst, x, xhx)
    
    
    Defined at Chebyshev_fliter.fpp lines 996-1021
    
    Parameters
    ----------
    veff : float array
    nst : int
    x : float array
    xhx : float array
    
    """
    _AresMainPy_pkg.f90wrap_rayleigh_quotient_real(veff=veff, nst=nst, x=x, xhx=xhx)

def cal_hx_real(veff, nst, v, hv):
    """
    cal_hx_real(veff, nst, v, hv)
    
    
    Defined at Chebyshev_fliter.fpp lines 1024-1043
    
    Parameters
    ----------
    veff : float array
    nst : int
    v : float array
    hv : float array
    
    =========================================
    print *,"nst",nst
    print *,'cal HX time -->',(ts-tt)/10000.d0
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    """
    _AresMainPy_pkg.f90wrap_cal_hx_real(veff=veff, nst=nst, v=v, hv=hv)

def estupb_real(k, veff, vec):
    """
    b = estupb_real(k, veff, vec)
    
    
    Defined at Chebyshev_fliter.fpp lines 1046-1089
    
    Parameters
    ----------
    k : int
    veff : float array
    vec : float array
    
    Returns
    -------
    b : float
    
    """
    b = _AresMainPy_pkg.f90wrap_estupb_real(k=k, veff=veff, vec=vec)
    return b

def chebyshev_filter_scaled_real(veff, x, m, a, b, al):
    """
    chebyshev_filter_scaled_real(veff, x, m, a, b, al)
    
    
    Defined at Chebyshev_fliter.fpp lines 1092-1135
    
    Parameters
    ----------
    veff : float array
    x : float array
    m : int
    a : float
    b : float
    al : float
    
    ===========================
    print*,"{-------------------"
    ===========================
    """
    _AresMainPy_pkg.f90wrap_chebyshev_filter_scaled_real(veff=veff, x=x, m=m, a=a, \
        b=b, al=al)

def grayleigh_ritz_real(veff, x, d):
    """
    grayleigh_ritz_real(veff, x, d)
    
    
    Defined at Chebyshev_fliter.fpp lines 1138-1181
    
    Parameters
    ----------
    veff : float array
    x : float array
    d : float array
    
    """
    _AresMainPy_pkg.f90wrap_grayleigh_ritz_real(veff=veff, x=x, d=d)

def cheby_filtering_grrr(veff, x, d):
    """
    cheby_filtering_grrr(veff, x, d)
    
    
    Defined at Chebyshev_fliter.fpp lines 1184-1202
    
    Parameters
    ----------
    veff : float array
    x : float array
    d : float array
    
    """
    _AresMainPy_pkg.f90wrap_cheby_filtering_grrr(veff=veff, x=x, d=d)

def cheby_filtering_prrr(veff, x, nfs, cfr, efr):
    """
    cheby_filtering_prrr(veff, x, nfs, cfr, efr)
    
    
    Defined at Chebyshev_fliter.fpp lines 1210-1230
    
    Parameters
    ----------
    veff : float array
    x : float array
    nfs : int
    cfr : float array
    efr : float array
    
    """
    _AresMainPy_pkg.f90wrap_cheby_filtering_prrr(veff=veff, x=x, nfs=nfs, cfr=cfr, \
        efr=efr)

def partialrayleighritz(veff, x, nfs, cfr, efr):
    """
    partialrayleighritz(veff, x, nfs, cfr, efr)
    
    
    Defined at Chebyshev_fliter.fpp lines 1233-1250
    
    Parameters
    ----------
    veff : float array
    x : float array
    nfs : int
    cfr : float array
    efr : float array
    
    """
    _AresMainPy_pkg.f90wrap_partialrayleighritz(veff=veff, x=x, nfs=nfs, cfr=cfr, \
        efr=efr)

def prr_orthnorm(veff, x, nfs, cfr, efr):
    """
    prr_orthnorm(veff, x, nfs, cfr, efr)
    
    
    Defined at Chebyshev_fliter.fpp lines 1253-1271
    
    Parameters
    ----------
    veff : float array
    x : float array
    nfs : int
    cfr : float array
    efr : float array
    
    """
    _AresMainPy_pkg.f90wrap_prr_orthnorm(veff=veff, x=x, nfs=nfs, cfr=cfr, efr=efr)

def iso_first_chescf_sto_rand(rhos, nev, psi, eval):
    """
    iso_first_chescf_sto_rand(rhos, nev, psi, eval)
    
    
    Defined at Chebyshev_fliter.fpp lines 1279-1375
    
    Parameters
    ----------
    rhos : float array
    nev : int
    psi : float array
    eval : float array
    
    """
    _AresMainPy_pkg.f90wrap_iso_first_chescf_sto_rand(rhos=rhos, nev=nev, psi=psi, \
        eval=eval)

def iso_cheby_init_sto_rand(nmax, nrand, initx_sto):
    """
    iso_cheby_init_sto_rand(nmax, nrand, initx_sto)
    
    
    Defined at Chebyshev_fliter.fpp lines 1379-1457
    
    Parameters
    ----------
    nmax : int
    nrand : int
    initx_sto : float array
    
    """
    _AresMainPy_pkg.f90wrap_iso_cheby_init_sto_rand(nmax=nmax, nrand=nrand, \
        initx_sto=initx_sto)

def iso_cheby_init_rand(nmax, nrand, initx_sto):
    """
    iso_cheby_init_rand(nmax, nrand, initx_sto)
    
    
    Defined at Chebyshev_fliter.fpp lines 1460-1498
    
    Parameters
    ----------
    nmax : int
    nrand : int
    initx_sto : float array
    
    """
    _AresMainPy_pkg.f90wrap_iso_cheby_init_rand(nmax=nmax, nrand=nrand, \
        initx_sto=initx_sto)

def iso_first_subspace_sto_rand(veff_3d, nmax, nev, initx, x, d):
    """
    iso_first_subspace_sto_rand(veff_3d, nmax, nev, initx, x, d)
    
    
    Defined at Chebyshev_fliter.fpp lines 1501-1628
    
    Parameters
    ----------
    veff_3d : float array
    nmax : int
    nev : int
    initx : float array
    x : float array
    d : float array
    
    """
    _AresMainPy_pkg.f90wrap_iso_first_subspace_sto_rand(veff_3d=veff_3d, nmax=nmax, \
        nev=nev, initx=initx, x=x, d=d)

def rayleigh_quotient_iso(veff_3d, nst, x, xhx):
    """
    rayleigh_quotient_iso(veff_3d, nst, x, xhx)
    
    
    Defined at Chebyshev_fliter.fpp lines 1632-1657
    
    Parameters
    ----------
    veff_3d : float array
    nst : int
    x : float array
    xhx : float array
    
    """
    _AresMainPy_pkg.f90wrap_rayleigh_quotient_iso(veff_3d=veff_3d, nst=nst, x=x, \
        xhx=xhx)

def cal_hx_iso(veff, nst, v, hv):
    """
    cal_hx_iso(veff, nst, v, hv)
    
    
    Defined at Chebyshev_fliter.fpp lines 1661-1689
    
    Parameters
    ----------
    veff : float array
    nst : int
    v : float array
    hv : float array
    
    =========================================
    """
    _AresMainPy_pkg.f90wrap_cal_hx_iso(veff=veff, nst=nst, v=v, hv=hv)

def cheby_filtering_grriso(veff_3d, x_sphe, d):
    """
    cheby_filtering_grriso(veff_3d, x_sphe, d)
    
    
    Defined at Chebyshev_fliter.fpp lines 1694-1715
    
    Parameters
    ----------
    veff_3d : float array
    x_sphe : float array
    d : float array
    
    """
    _AresMainPy_pkg.f90wrap_cheby_filtering_grriso(veff_3d=veff_3d, x_sphe=x_sphe, \
        d=d)

def estupb_iso(k, veff_3d, vec):
    """
    b = estupb_iso(k, veff_3d, vec)
    
    
    Defined at Chebyshev_fliter.fpp lines 1720-1783
    
    Parameters
    ----------
    k : int
    veff_3d : float array
    vec : float array
    
    Returns
    -------
    b : float
    
    """
    b = _AresMainPy_pkg.f90wrap_estupb_iso(k=k, veff_3d=veff_3d, vec=vec)
    return b

def chebyshev_filter_scaled_iso(veff_3d, x, m, a, b, al):
    """
    chebyshev_filter_scaled_iso(veff_3d, x, m, a, b, al)
    
    
    Defined at Chebyshev_fliter.fpp lines 1788-1841
    
    Parameters
    ----------
    veff_3d : float array
    x : float array
    m : int
    a : float
    b : float
    al : float
    
    """
    _AresMainPy_pkg.f90wrap_chebyshev_filter_scaled_iso(veff_3d=veff_3d, x=x, m=m, \
        a=a, b=b, al=al)

def grayleigh_ritz_iso(veff_3d, x, d):
    """
    grayleigh_ritz_iso(veff_3d, x, d)
    
    
    Defined at Chebyshev_fliter.fpp lines 1845-1925
    
    Parameters
    ----------
    veff_3d : float array
    x : float array
    d : float array
    
    """
    _AresMainPy_pkg.f90wrap_grayleigh_ritz_iso(veff_3d=veff_3d, x=x, d=d)

def first_chescf_sto_gamma(veff, nev, psi, eval):
    """
    first_chescf_sto_gamma(veff, nev, psi, eval)
    
    
    Defined at Chebyshev_fliter.fpp lines 1929-1990
    
    Parameters
    ----------
    veff : float array
    nev : int
    psi : float array
    eval : float array
    
    """
    _AresMainPy_pkg.f90wrap_first_chescf_sto_gamma(veff=veff, nev=nev, psi=psi, \
        eval=eval)

def cheby_init_sto_gamma(nmax, npw, initx_sto):
    """
    cheby_init_sto_gamma(nmax, npw, initx_sto)
    
    
    Defined at Chebyshev_fliter.fpp lines 1993-2083
    
    Parameters
    ----------
    nmax : int
    npw : int
    initx_sto : float array
    
    """
    _AresMainPy_pkg.f90wrap_cheby_init_sto_gamma(nmax=nmax, npw=npw, \
        initx_sto=initx_sto)

def first_subspace_stopw_gamma(ik, veff, nmax, nev, initx, x, d):
    """
    first_subspace_stopw_gamma(ik, veff, nmax, nev, initx, x, d)
    
    
    Defined at Chebyshev_fliter.fpp lines 2086-2185
    
    Parameters
    ----------
    ik : int
    veff : float array
    nmax : int
    nev : int
    initx : float array
    x : float array
    d : float array
    
    """
    _AresMainPy_pkg.f90wrap_first_subspace_stopw_gamma(ik=ik, veff=veff, nmax=nmax, \
        nev=nev, initx=initx, x=x, d=d)

def rayleigh_quotient_gamma(ik, veff, nst, x, xhx):
    """
    rayleigh_quotient_gamma(ik, veff, nst, x, xhx)
    
    
    Defined at Chebyshev_fliter.fpp lines 2188-2213
    
    Parameters
    ----------
    ik : int
    veff : float array
    nst : int
    x : float array
    xhx : float array
    
    """
    _AresMainPy_pkg.f90wrap_rayleigh_quotient_gamma(ik=ik, veff=veff, nst=nst, x=x, \
        xhx=xhx)

def cal_hx_gamma(ik, veff, nst, v, hv):
    """
    cal_hx_gamma(ik, veff, nst, v, hv)
    
    
    Defined at Chebyshev_fliter.fpp lines 2216-2233
    
    Parameters
    ----------
    ik : int
    veff : float array
    nst : int
    v : float array
    hv : float array
    
    """
    _AresMainPy_pkg.f90wrap_cal_hx_gamma(ik=ik, veff=veff, nst=nst, v=v, hv=hv)

def cheby_filter_rr_gamma(ik, veff, x, d):
    """
    cheby_filter_rr_gamma(ik, veff, x, d)
    
    
    Defined at Chebyshev_fliter.fpp lines 2235-2269
    
    Parameters
    ----------
    ik : int
    veff : float array
    x : float array
    d : float array
    
    """
    _AresMainPy_pkg.f90wrap_cheby_filter_rr_gamma(ik=ik, veff=veff, x=x, d=d)

def estupb_gamma(k, ik, veff, vec):
    """
    b = estupb_gamma(k, ik, veff, vec)
    
    
    Defined at Chebyshev_fliter.fpp lines 2271-2330
    
    Parameters
    ----------
    k : int
    ik : int
    veff : float array
    vec : float array
    
    Returns
    -------
    b : float
    
    """
    b = _AresMainPy_pkg.f90wrap_estupb_gamma(k=k, ik=ik, veff=veff, vec=vec)
    return b

def chebyshev_filter_scaled_gamma(ik, veff, x, m, a, b, al):
    """
    chebyshev_filter_scaled_gamma(ik, veff, x, m, a, b, al)
    
    
    Defined at Chebyshev_fliter.fpp lines 2332-2376
    
    Parameters
    ----------
    ik : int
    veff : float array
    x : float array
    m : int
    a : float
    b : float
    al : float
    
    """
    _AresMainPy_pkg.f90wrap_chebyshev_filter_scaled_gamma(ik=ik, veff=veff, x=x, \
        m=m, a=a, b=b, al=al)

def grayleigh_ritz_gamma(ik, veff, x, d):
    """
    grayleigh_ritz_gamma(ik, veff, x, d)
    
    
    Defined at Chebyshev_fliter.fpp lines 2378-2436
    
    Parameters
    ----------
    ik : int
    veff : float array
    x : float array
    d : float array
    
    """
    _AresMainPy_pkg.f90wrap_grayleigh_ritz_gamma(ik=ik, veff=veff, x=x, d=d)

def get_iiii():
    """
    Element iiii ftype=integer(i4b) pytype=int
    
    
    Defined at Chebyshev_fliter.fpp line 16
    
    """
    return _AresMainPy_pkg.f90wrap_chebyshev_module__get__iiii()

def set_iiii(iiii):
    _AresMainPy_pkg.f90wrap_chebyshev_module__set__iiii(iiii)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "chebyshev_module".')

for func in _dt_array_initialisers:
    func()

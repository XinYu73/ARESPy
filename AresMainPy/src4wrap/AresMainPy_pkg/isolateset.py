"""
Module isolateset


Defined at Isolate_module.fpp lines 5-2205

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def isolatevcoulomb_a(rho_s, v_s):
    """
    isolatevcoulomb_a(rho_s, v_s)
    
    
    Defined at Isolate_module.fpp lines 25-93
    
    Parameters
    ----------
    rho_s : float array
    v_s : float array
    
    ===========================================================
    """
    _AresMainPy_pkg.f90wrap_isolatevcoulomb_a(rho_s=rho_s, v_s=v_s)

def evaluateboundrygrid_b(rhos):
    """
    evaluateboundrygrid_b(rhos)
    
    
    Defined at Isolate_module.fpp lines 97-122
    
    Parameters
    ----------
    rhos : float array
    
    """
    _AresMainPy_pkg.f90wrap_evaluateboundrygrid_b(rhos=rhos)

def endevaluate_b():
    """
    endevaluate_b()
    
    
    Defined at Isolate_module.fpp lines 126-137
    
    
    ------------
    """
    _AresMainPy_pkg.f90wrap_endevaluate_b()

def calculateboundrypotential_b(rho, lmax):
    """
    calculateboundrypotential_b(rho, lmax)
    
    
    Defined at Isolate_module.fpp lines 141-177
    
    Parameters
    ----------
    rho : float array
    lmax : int
    
    """
    _AresMainPy_pkg.f90wrap_calculateboundrypotential_b(rho=rho, lmax=lmax)

def laplaceequationslover_b(iter):
    """
    laplaceequationslover_b(iter)
    
    
    Defined at Isolate_module.fpp lines 181-249
    
    Parameters
    ----------
    iter : int
    
    ---------------------------------------------------------
    """
    _AresMainPy_pkg.f90wrap_laplaceequationslover_b(iter=iter)

def vcoulombassignment_b(iter, niter_poisson, phi):
    """
    vcoulombassignment_b(iter, niter_poisson, phi)
    
    
    Defined at Isolate_module.fpp lines 253-263
    
    Parameters
    ----------
    iter : int
    niter_poisson : int
    phi : float array
    
    -------------------------------------------------
    """
    _AresMainPy_pkg.f90wrap_vcoulombassignment_b(iter=iter, \
        niter_poisson=niter_poisson, phi=phi)

def cal_qlm(rho, lmax, q_l0, q_lm):
    """
    cal_qlm(rho, lmax, q_l0, q_lm)
    
    
    Defined at Isolate_module.fpp lines 267-310
    
    Parameters
    ----------
    rho : float array
    lmax : int
    q_l0 : float array
    q_lm : complex array
    
    """
    _AresMainPy_pkg.f90wrap_cal_qlm(rho=rho, lmax=lmax, q_l0=q_l0, q_lm=q_lm)

def car2spe(orig, x, y, z):
    """
    r, cost, sint, cosp, sinp = car2spe(orig, x, y, z)
    
    
    Defined at Isolate_module.fpp lines 314-345
    
    Parameters
    ----------
    orig : float array
    x : int
    y : int
    z : int
    
    Returns
    -------
    r : float
    cost : float
    sint : float
    cosp : float
    sinp : float
    
    """
    r, cost, sint, cosp, sinp = _AresMainPy_pkg.f90wrap_car2spe(orig=orig, x=x, y=y, \
        z=z)
    return r, cost, sint, cosp, sinp

def cal_plm(lmax, x, z, plm):
    """
    cal_plm(lmax, x, z, plm)
    
    
    Defined at Isolate_module.fpp lines 349-377
    
    Parameters
    ----------
    lmax : int
    x : float
    z : float
    plm : float array
    
    """
    _AresMainPy_pkg.f90wrap_cal_plm(lmax=lmax, x=x, z=z, plm=plm)

def cal1pot(iex, iey, iez, orig, lmax, q_l0, q_lm):
    """
    phi = cal1pot(iex, iey, iez, orig, lmax, q_l0, q_lm)
    
    
    Defined at Isolate_module.fpp lines 381-411
    
    Parameters
    ----------
    iex : int
    iey : int
    iez : int
    orig : float array
    lmax : int
    q_l0 : float array
    q_lm : complex array
    
    Returns
    -------
    phi : float
    
    """
    phi = _AresMainPy_pkg.f90wrap_cal1pot(iex=iex, iey=iey, iez=iez, orig=orig, \
        lmax=lmax, q_l0=q_l0, q_lm=q_lm)
    return phi

def apply_boundary(lmax, orig, q_l0, q_lm, bphi):
    """
    apply_boundary(lmax, orig, q_l0, q_lm, bphi)
    
    
    Defined at Isolate_module.fpp lines 415-446
    
    Parameters
    ----------
    lmax : int
    orig : float array
    q_l0 : float array
    q_lm : complex array
    bphi : float array
    
    """
    _AresMainPy_pkg.f90wrap_apply_boundary(lmax=lmax, orig=orig, q_l0=q_l0, \
        q_lm=q_lm, bphi=bphi)

def apply_boundary00(rho, bphi):
    """
    apply_boundary00(rho, bphi)
    
    
    Defined at Isolate_module.fpp lines 449-483
    
    Parameters
    ----------
    rho : float array
    bphi : float array
    
    """
    _AresMainPy_pkg.f90wrap_apply_boundary00(rho=rho, bphi=bphi)

def laplb(ford, bf, f):
    """
    laplb(ford, bf, f)
    
    
    Defined at Isolate_module.fpp lines 486-499
    
    Parameters
    ----------
    ford : int
    bf : float array
    f : float array
    
    """
    _AresMainPy_pkg.f90wrap_laplb(ford=ford, bf=bf, f=f)

def lapla(ford, dxyz, f, af):
    """
    lapla(ford, dxyz, f, af)
    
    
    Defined at Isolate_module.fpp lines 503-520
    
    Parameters
    ----------
    ford : int
    dxyz : float array
    f : float array
    af : float array
    
    """
    _AresMainPy_pkg.f90wrap_lapla(ford=ford, dxyz=dxyz, f=f, af=af)

def lapla_coe(ford, dxyz, coe):
    """
    lapla_coe(ford, dxyz, coe)
    
    
    Defined at Isolate_module.fpp lines 524-537
    
    Parameters
    ----------
    ford : int
    dxyz : float array
    coe : float array
    
    """
    _AresMainPy_pkg.f90wrap_lapla_coe(ford=ford, dxyz=dxyz, coe=coe)

def damped_jacobi_iterate(ford, dxyz, func, rhs):
    """
    damped_jacobi_iterate(ford, dxyz, func, rhs)
    
    
    Defined at Isolate_module.fpp lines 541-595
    
    Parameters
    ----------
    ford : int
    dxyz : float array
    func : float array
    rhs : float array
    
    """
    _AresMainPy_pkg.f90wrap_damped_jacobi_iterate(ford=ford, dxyz=dxyz, func=func, \
        rhs=rhs)

def gauss_seidel_iterate(ford, dxyz, func, rhs):
    """
    gauss_seidel_iterate(ford, dxyz, func, rhs)
    
    
    Defined at Isolate_module.fpp lines 600-656
    
    Parameters
    ----------
    ford : int
    dxyz : float array
    func : float array
    rhs : float array
    
    """
    _AresMainPy_pkg.f90wrap_gauss_seidel_iterate(ford=ford, dxyz=dxyz, func=func, \
        rhs=rhs)

def refine(var, nxyz, nxyz_fine, var_fine):
    """
    refine(var, nxyz, nxyz_fine, var_fine)
    
    
    Defined at Isolate_module.fpp lines 661-674
    
    Parameters
    ----------
    var : float array
    nxyz : int array
    nxyz_fine : int array
    var_fine : float array
    
    """
    _AresMainPy_pkg.f90wrap_refine(var=var, nxyz=nxyz, nxyz_fine=nxyz_fine, \
        var_fine=var_fine)

def interpolate(var1, nxyz1, nxyz2, var2):
    """
    interpolate(var1, nxyz1, nxyz2, var2)
    
    
    Defined at Isolate_module.fpp lines 679-741
    
    Parameters
    ----------
    var1 : float array
    nxyz1 : int array
    nxyz2 : int array
    var2 : float array
    
    -----------------------------------------------------------------------
    """
    _AresMainPy_pkg.f90wrap_interpolate(var1=var1, nxyz1=nxyz1, nxyz2=nxyz2, \
        var2=var2)

def smooth(f):
    """
    smooth(f)
    
    
    Defined at Isolate_module.fpp lines 746-797
    
    Parameters
    ----------
    f : float array
    
    """
    _AresMainPy_pkg.f90wrap_smooth(f=f)

def restrict(var, nxyz, nxyz_coarse, var_coarse):
    """
    restrict(var, nxyz, nxyz_coarse, var_coarse)
    
    
    Defined at Isolate_module.fpp lines 802-821
    
    Parameters
    ----------
    var : float array
    nxyz : int array
    nxyz_coarse : int array
    var_coarse : float array
    
    """
    _AresMainPy_pkg.f90wrap_restrict(var=var, nxyz=nxyz, nxyz_coarse=nxyz_coarse, \
        var_coarse=var_coarse)

def residual(ford, dxyz, f, rhs, res):
    """
    residual(ford, dxyz, f, rhs, res)
    
    
    Defined at Isolate_module.fpp lines 826-844
    
    Parameters
    ----------
    ford : int
    dxyz : float array
    f : float array
    rhs : float array
    res : float array
    
    """
    _AresMainPy_pkg.f90wrap_residual(ford=ford, dxyz=dxyz, f=f, rhs=rhs, res=res)

def v_cycle(ford, f, rhs, dxyz, nvc):
    """
    v_cycle(ford, f, rhs, dxyz, nvc)
    
    
    Defined at Isolate_module.fpp lines 849-965
    
    Parameters
    ----------
    ford : int
    f : float array
    rhs : float array
    dxyz : float array
    nvc : int
    
    """
    _AresMainPy_pkg.f90wrap_v_cycle(ford=ford, f=f, rhs=rhs, dxyz=dxyz, nvc=nvc)

def csix_coe(dxyz, coe):
    """
    csix_coe(dxyz, coe)
    
    
    Defined at Isolate_module.fpp lines 970-983
    
    Parameters
    ----------
    dxyz : float array
    coe : float array
    
    """
    _AresMainPy_pkg.f90wrap_csix_coe(dxyz=dxyz, coe=coe)

def lapla_csix(f, af, coe):
    """
    lapla_csix(f, af, coe)
    
    
    Defined at Isolate_module.fpp lines 987-1021
    
    Parameters
    ----------
    f : float array
    af : float array
    coe : float array
    
    """
    _AresMainPy_pkg.f90wrap_lapla_csix(f=f, af=af, coe=coe)

def laplb_csix(bff, f):
    """
    laplb_csix(bff, f)
    
    
    Defined at Isolate_module.fpp lines 1025-1066
    
    Parameters
    ----------
    bff : float array
    f : float array
    
    """
    _AresMainPy_pkg.f90wrap_laplb_csix(bff=bff, f=f)

def cfour_coe(gaps, coe):
    """
    cfour_coe(gaps, coe)
    
    
    Defined at Isolate_module.fpp lines 1070-1082
    
    Parameters
    ----------
    gaps : float array
    coe : float array
    
    """
    _AresMainPy_pkg.f90wrap_cfour_coe(gaps=gaps, coe=coe)

def lapla_cfour(u, au, coe):
    """
    lapla_cfour(u, au, coe)
    
    
    Defined at Isolate_module.fpp lines 1086-1106
    
    Parameters
    ----------
    u : float array
    au : float array
    coe : float array
    
    """
    _AresMainPy_pkg.f90wrap_lapla_cfour(u=u, au=au, coe=coe)

def laplb_cfour(rhs, nrhs):
    """
    laplb_cfour(rhs, nrhs)
    
    
    Defined at Isolate_module.fpp lines 1110-1132
    
    Parameters
    ----------
    rhs : float array
    nrhs : float array
    
    """
    _AresMainPy_pkg.f90wrap_laplb_cfour(rhs=rhs, nrhs=nrhs)

def calclm(lmax, clm):
    """
    calclm(lmax, clm)
    
    
    Defined at Isolate_module.fpp lines 1136-1145
    
    Parameters
    ----------
    lmax : int
    clm : float array
    
    """
    _AresMainPy_pkg.f90wrap_calclm(lmax=lmax, clm=clm)

def c(l, m):
    """
    c = c(l, m)
    
    
    Defined at Isolate_module.fpp lines 1149-1162
    
    Parameters
    ----------
    l : int
    m : int
    
    Returns
    -------
    c : float
    
    """
    c = _AresMainPy_pkg.f90wrap_c(l=l, m=m)
    return c

def rcs(n1, n2, n3, dr, cost, sint, cns, indx):
    """
    rcs(n1, n2, n3, dr, cost, sint, cns, indx)
    
    
    Defined at Isolate_module.fpp lines 1166-1195
    
    Parameters
    ----------
    n1 : int
    n2 : int
    n3 : int
    dr : float array
    cost : float array
    sint : float array
    cns : complex array
    indx : int array
    
    """
    _AresMainPy_pkg.f90wrap_rcs(n1=n1, n2=n2, n3=n3, dr=dr, cost=cost, sint=sint, \
        cns=cns, indx=indx)

def vhartree_fmm(rhos, vcoulomb):
    """
    vhartree_fmm(rhos, vcoulomb)
    
    
    Defined at Isolate_module.fpp lines 1199-1371
    
    Parameters
    ----------
    rhos : float array
    vcoulomb : float array
    
    ===============================================
    ##CALCULATE VHARTREE VIA FMM CALLED FROM LFMM3D
    ===============================================
    """
    _AresMainPy_pkg.f90wrap_vhartree_fmm(rhos=rhos, vcoulomb=vcoulomb)

def vhartree_direct(rhos, vcoulomb):
    """
    vhartree_direct(rhos, vcoulomb)
    
    
    Defined at Isolate_module.fpp lines 1374-1459
    
    Parameters
    ----------
    rhos : float array
    vcoulomb : float array
    
    ===============================================
    ##CALCULATE VHARTREE VIA FMM CALLED FROM LFMM3D
    ===============================================
    """
    _AresMainPy_pkg.f90wrap_vhartree_direct(rhos=rhos, vcoulomb=vcoulomb)

def cal_vsrcpot(size_src, src, pos_src, size_tar, pot_tar, pos_tar):
    """
    cal_vsrcpot(size_src, src, pos_src, size_tar, pot_tar, pos_tar)
    
    
    Defined at Isolate_module.fpp lines 1462-1492
    
    Parameters
    ----------
    size_src : int
    src : float array
    pos_src : float array
    size_tar : int
    pot_tar : float array
    pos_tar : float array
    
    """
    _AresMainPy_pkg.f90wrap_cal_vsrcpot(size_src=size_src, src=src, pos_src=pos_src, \
        size_tar=size_tar, pot_tar=pot_tar, pos_tar=pos_tar)

def vhart_iso(rho, vhart):
    """
    vhart_iso(rho, vhart)
    
    
    Defined at Isolate_module.fpp lines 1494-1524
    
    Parameters
    ----------
    rho : float array
    vhart : float array
    
    """
    _AresMainPy_pkg.f90wrap_vhart_iso(rho=rho, vhart=vhart)

def poisson_mcm(rho, pot):
    """
    poisson_mcm(rho, pot)
    
    
    Defined at Isolate_module.fpp lines 1526-1591
    
    Parameters
    ----------
    rho : float array
    pot : float array
    
    """
    _AresMainPy_pkg.f90wrap_poisson_mcm(rho=rho, pot=pot)

def mcm_calvcorr(q_l0, q_lm, v_corr):
    """
    mcm_calvcorr(q_l0, q_lm, v_corr)
    
    
    Defined at Isolate_module.fpp lines 1593-1740
    
    Parameters
    ----------
    q_l0 : float array
    q_lm : complex array
    v_corr : float array
    
    """
    _AresMainPy_pkg.f90wrap_mcm_calvcorr(q_l0=q_l0, q_lm=q_lm, v_corr=v_corr)

def mcm_cal_naux(q_l0, q_lm, rhoaux):
    """
    mcm_cal_naux(q_l0, q_lm, rhoaux)
    
    
    Defined at Isolate_module.fpp lines 1742-1801
    
    Parameters
    ----------
    q_l0 : float array
    q_lm : complex array
    rhoaux : float array
    
    """
    _AresMainPy_pkg.f90wrap_mcm_cal_naux(q_l0=q_l0, q_lm=q_lm, rhoaux=rhoaux)

def cal_ylm(mycost, mysint, e_phi, ylm):
    """
    cal_ylm(mycost, mysint, e_phi, ylm)
    
    
    Defined at Isolate_module.fpp lines 1804-1823
    
    Parameters
    ----------
    mycost : float
    mysint : float
    e_phi : complex
    ylm : complex array
    
    """
    _AresMainPy_pkg.f90wrap_cal_ylm(mycost=mycost, mysint=mysint, e_phi=e_phi, \
        ylm=ylm)

def rcs_rec(ng1, ng2, ng3, cost, sint, cns):
    """
    rcs_rec(ng1, ng2, ng3, cost, sint, cns)
    
    
    Defined at Isolate_module.fpp lines 1825-1858
    
    Parameters
    ----------
    ng1 : int
    ng2 : int
    ng3 : int
    cost : float array
    sint : float array
    cns : complex array
    
    """
    _AresMainPy_pkg.f90wrap_rcs_rec(ng1=ng1, ng2=ng2, ng3=ng3, cost=cost, sint=sint, \
        cns=cns)

def car2spe_rec(orig, x, y, z):
    """
    cost, sint, cosp, sinp = car2spe_rec(orig, x, y, z)
    
    
    Defined at Isolate_module.fpp lines 1860-1898
    
    Parameters
    ----------
    orig : float array
    x : float
    y : float
    z : float
    
    Returns
    -------
    cost : float
    sint : float
    cosp : float
    sinp : float
    
    """
    cost, sint, cosp, sinp = _AresMainPy_pkg.f90wrap_car2spe_rec(orig=orig, x=x, \
        y=y, z=z)
    return cost, sint, cosp, sinp

def integrade_i(l, x):
    """
    integrade_i = integrade_i(l, x)
    
    
    Defined at Isolate_module.fpp lines 1900-2130
    
    Parameters
    ----------
    l : int
    x : float
    
    Returns
    -------
    integrade_i : float
    
    """
    integrade_i = _AresMainPy_pkg.f90wrap_integrade_i(l=l, x=x)
    return integrade_i

def gmg_cg_hart(rhos, vcoulomb):
    """
    gmg_cg_hart(rhos, vcoulomb)
    
    
    Defined at Isolate_module.fpp lines 2132-2171
    
    Parameters
    ----------
    rhos : float array
    vcoulomb : float array
    
    ===========================================================
    """
    _AresMainPy_pkg.f90wrap_gmg_cg_hart(rhos=rhos, vcoulomb=vcoulomb)

def cutoff_method(rho, vh):
    """
    cutoff_method(rho, vh)
    
    
    Defined at Isolate_module.fpp lines 2174-2204
    
    Parameters
    ----------
    rho : float array
    vh : float array
    
    """
    _AresMainPy_pkg.f90wrap_cutoff_method(rho=rho, vh=vh)

def get_array_center():
    """
    Element center ftype=real(dp) pytype=float
    
    
    Defined at Isolate_module.fpp line 16
    
    """
    global center
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_isolateset__array__center(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        center = _arrays[array_handle]
    else:
        center = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_isolateset__array__center)
        _arrays[array_handle] = center
    return center

def set_array_center(center):
    center[...] = center

def get_array_celllengthn():
    """
    Element celllengthn ftype=integer(i4b) pytype=int
    
    
    Defined at Isolate_module.fpp line 17
    
    """
    global celllengthn
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_isolateset__array__celllengthn(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        celllengthn = _arrays[array_handle]
    else:
        celllengthn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_isolateset__array__celllengthn)
        _arrays[array_handle] = celllengthn
    return celllengthn

def set_array_celllengthn(celllengthn):
    celllengthn[...] = celllengthn

def get_thickness():
    """
    Element thickness ftype=integer(i4b) pytype=int
    
    
    Defined at Isolate_module.fpp line 17
    
    """
    return _AresMainPy_pkg.f90wrap_isolateset__get__thickness()

def set_thickness(thickness):
    _AresMainPy_pkg.f90wrap_isolateset__set__thickness(thickness)

def get_cellleft():
    """
    Element cellleft ftype=integer(i4b) pytype=int
    
    
    Defined at Isolate_module.fpp line 17
    
    """
    return _AresMainPy_pkg.f90wrap_isolateset__get__cellleft()

def set_cellleft(cellleft):
    _AresMainPy_pkg.f90wrap_isolateset__set__cellleft(cellleft)

def get_array_cellright():
    """
    Element cellright ftype=integer(i4b) pytype=int
    
    
    Defined at Isolate_module.fpp line 17
    
    """
    global cellright
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_isolateset__array__cellright(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        cellright = _arrays[array_handle]
    else:
        cellright = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_isolateset__array__cellright)
        _arrays[array_handle] = cellright
    return cellright

def set_array_cellright(cellright):
    cellright[...] = cellright

def get_array_bandright():
    """
    Element bandright ftype=integer(i4b) pytype=int
    
    
    Defined at Isolate_module.fpp line 17
    
    """
    global bandright
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_isolateset__array__bandright(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        bandright = _arrays[array_handle]
    else:
        bandright = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_isolateset__array__bandright)
        _arrays[array_handle] = bandright
    return bandright

def set_array_bandright(bandright):
    bandright[...] = bandright

def get_array_boundaryrhos():
    """
    Element boundaryrhos ftype=real(dp) pytype=float
    
    
    Defined at Isolate_module.fpp line 19
    
    """
    global boundaryrhos
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_isolateset__array__boundaryrhos(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        boundaryrhos = _arrays[array_handle]
    else:
        boundaryrhos = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_isolateset__array__boundaryrhos)
        _arrays[array_handle] = boundaryrhos
    return boundaryrhos

def set_array_boundaryrhos(boundaryrhos):
    boundaryrhos[...] = boundaryrhos

def get_array_boundaryvcoulomb():
    """
    Element boundaryvcoulomb ftype=real(dp) pytype=float
    
    
    Defined at Isolate_module.fpp line 19
    
    """
    global boundaryvcoulomb
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_isolateset__array__boundaryvcoulomb(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        boundaryvcoulomb = _arrays[array_handle]
    else:
        boundaryvcoulomb = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_isolateset__array__boundaryvcoulomb)
        _arrays[array_handle] = boundaryvcoulomb
    return boundaryvcoulomb

def set_array_boundaryvcoulomb(boundaryvcoulomb):
    boundaryvcoulomb[...] = boundaryvcoulomb

def get_array_ap():
    """
    Element ap ftype=real(dp) pytype=float
    
    
    Defined at Isolate_module.fpp line 19
    
    """
    global ap
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_isolateset__array__ap(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        ap = _arrays[array_handle]
    else:
        ap = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_isolateset__array__ap)
        _arrays[array_handle] = ap
    return ap

def set_array_ap(ap):
    ap[...] = ap

def get_array_aphi():
    """
    Element aphi ftype=real(dp) pytype=float
    
    
    Defined at Isolate_module.fpp line 19
    
    """
    global aphi
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_isolateset__array__aphi(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        aphi = _arrays[array_handle]
    else:
        aphi = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_isolateset__array__aphi)
        _arrays[array_handle] = aphi
    return aphi

def set_array_aphi(aphi):
    aphi[...] = aphi

def get_array_res():
    """
    Element res ftype=real(dp) pytype=float
    
    
    Defined at Isolate_module.fpp line 19
    
    """
    global res
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_isolateset__array__res(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        res = _arrays[array_handle]
    else:
        res = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_isolateset__array__res)
        _arrays[array_handle] = res
    return res

def set_array_res(res):
    res[...] = res

def get_array_res1():
    """
    Element res1 ftype=real(dp) pytype=float
    
    
    Defined at Isolate_module.fpp line 19
    
    """
    global res1
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_isolateset__array__res1(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        res1 = _arrays[array_handle]
    else:
        res1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_isolateset__array__res1)
        _arrays[array_handle] = res1
    return res1

def set_array_res1(res1):
    res1[...] = res1

def get_array_vcoulomb_old():
    """
    Element vcoulomb_old ftype=real(dp) pytype=float
    
    
    Defined at Isolate_module.fpp line 22
    
    """
    global vcoulomb_old
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_isolateset__array__vcoulomb_old(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        vcoulomb_old = _arrays[array_handle]
    else:
        vcoulomb_old = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_isolateset__array__vcoulomb_old)
        _arrays[array_handle] = vcoulomb_old
    return vcoulomb_old

def set_array_vcoulomb_old(vcoulomb_old):
    vcoulomb_old[...] = vcoulomb_old


_array_initialisers = [get_array_center, get_array_celllengthn, \
    get_array_cellright, get_array_bandright, get_array_boundaryrhos, \
    get_array_boundaryvcoulomb, get_array_ap, get_array_aphi, get_array_res, \
    get_array_res1, get_array_vcoulomb_old]
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "isolateset".')

for func in _dt_array_initialisers:
    func()

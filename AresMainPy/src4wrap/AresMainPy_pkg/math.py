"""
Module math


Defined at Math.fpp lines 5-3439

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def change_case(instr, str, fun):
    """
    change_case(instr, str, fun)
    
    
    Defined at Math.fpp lines 55-95
    
    Parameters
    ----------
    instr : str
    str : str
    fun : int
    
    """
    _AresMainPy_pkg.f90wrap_change_case(instr=instr, str=str, fun=fun)

def find_keywords(str, ch_mark, id_key, id_value):
    """
    find_keywords(str, ch_mark, id_key, id_value)
    
    
    Defined at Math.fpp lines 99-128
    
    Parameters
    ----------
    str : str
    ch_mark : str
    id_key : int
    id_value : int
    
    """
    _AresMainPy_pkg.f90wrap_find_keywords(str=str, ch_mark=ch_mark, id_key=id_key, \
        id_value=id_value)

def find_nword(str, ch_comma, nword):
    """
    find_nword(str, ch_comma, nword)
    
    
    Defined at Math.fpp lines 132-158
    
    Parameters
    ----------
    str : str
    ch_comma : str
    nword : int
    
    """
    _AresMainPy_pkg.f90wrap_find_nword(str=str, ch_comma=ch_comma, nword=nword)

def det(matrix):
    """
    det = det(matrix)
    
    
    Defined at Math.fpp lines 190-194
    
    Parameters
    ----------
    matrix : float array
    
    Returns
    -------
    det : float
    
    """
    det = _AresMainPy_pkg.f90wrap_det(matrix=matrix)
    return det

def inv_33(m):
    """
    inv_33 = inv_33(m)
    
    
    Defined at Math.fpp lines 198-228
    
    Parameters
    ----------
    m : float array
    
    Returns
    -------
    inv_33 : float array
    
    """
    inv_33 = _AresMainPy_pkg.f90wrap_inv_33(m=m)
    return inv_33

def lindg(eta, lambda_, mu):
    """
    lindg = lindg(eta, lambda_, mu)
    
    
    Defined at Math.fpp lines 232-278
    
    Parameters
    ----------
    eta : float
    lambda_ : float
    mu : float
    
    Returns
    -------
    lindg : float
    
    """
    lindg = _AresMainPy_pkg.f90wrap_lindg(eta=eta, lambda_=lambda_, mu=mu)
    return lindg

def int_to_char(int_bn):
    """
    int_to_char = int_to_char(int_bn)
    
    
    Defined at Math.fpp lines 282-303
    
    Parameters
    ----------
    int_bn : int
    
    Returns
    -------
    int_to_char : str
    
    -----------------------------------------------------------------------
    """
    int_to_char = _AresMainPy_pkg.f90wrap_int_to_char(int_bn=int_bn)
    return int_to_char

def lat2matrix(lat_para, lat_mat, flag):
    """
    lat2matrix(lat_para, lat_mat, flag)
    
    
    Defined at Math.fpp lines 307-350
    
    Parameters
    ----------
    lat_para : float array
    lat_mat : float array
    flag : int
    
    """
    _AresMainPy_pkg.f90wrap_lat2matrix(lat_para=lat_para, lat_mat=lat_mat, \
        flag=flag)

def one2three(id, n_dens, pos):
    """
    one2three(id, n_dens, pos)
    
    
    Defined at Math.fpp lines 354-376
    
    Parameters
    ----------
    id : int
    n_dens : int array
    pos : int array
    
    """
    _AresMainPy_pkg.f90wrap_one2three(id=id, n_dens=n_dens, pos=pos)

def init_random_seed():
    """
    init_random_seed()
    
    
    Defined at Math.fpp lines 470-534
    
    
    """
    _AresMainPy_pkg.f90wrap_init_random_seed()

def atom_mass(atom_name, mass):
    """
    atom_mass(atom_name, mass)
    
    
    Defined at Math.fpp lines 694-835
    
    Parameters
    ----------
    atom_name : str
    mass : float
    
    --------------------------------------------
    """
    _AresMainPy_pkg.f90wrap_atom_mass(atom_name=atom_name, mass=mass)

def newton_inter(n, x, y, m, tx, ty):
    """
    newton_inter(n, x, y, m, tx, ty)
    
    
    Defined at Math.fpp lines 839-862
    
    Parameters
    ----------
    n : int
    x : float array
    y : float array
    m : int
    tx : float array
    ty : float array
    
    ----------------------------------------------------
    """
    _AresMainPy_pkg.f90wrap_newton_inter(n=n, x=x, y=y, m=m, tx=tx, ty=ty)

def diag(n, ina, w, q):
    """
    diag(n, ina, w, q)
    
    
    Defined at Math.fpp lines 866-1000
    
    Parameters
    ----------
    n : int
    ina : float array
    w : float array
    q : float array
    
    """
    _AresMainPy_pkg.f90wrap_diag(n=n, ina=ina, w=w, q=q)

def boltzmann_distribution(rnull, width):
    """
    boltzmann_distribution = boltzmann_distribution(rnull, width)
    
    
    Defined at Math.fpp lines 1004-1019
    
    Parameters
    ----------
    rnull : float
    width : float
    
    Returns
    -------
    boltzmann_distribution : float
    
    """
    boltzmann_distribution = \
        _AresMainPy_pkg.f90wrap_boltzmann_distribution(rnull=rnull, width=width)
    return boltzmann_distribution

def dir2car(cry_coo, ort_coo, lat):
    """
    dir2car(cry_coo, ort_coo, lat)
    
    
    Defined at Math.fpp lines 1026-1036
    
    Parameters
    ----------
    cry_coo : float array
    ort_coo : float array
    lat : float array
    
    """
    _AresMainPy_pkg.f90wrap_dir2car(cry_coo=cry_coo, ort_coo=ort_coo, lat=lat)

def car2dir(ort_coo, cry_coo, lat):
    """
    car2dir(ort_coo, cry_coo, lat)
    
    
    Defined at Math.fpp lines 1042-1053
    
    Parameters
    ----------
    ort_coo : float array
    cry_coo : float array
    lat : float array
    
    """
    _AresMainPy_pkg.f90wrap_car2dir(ort_coo=ort_coo, cry_coo=cry_coo, lat=lat)

def thr2mat(n1, n2, n3, i, j, k):
    """
    dimnu = thr2mat(n1, n2, n3, i, j, k)
    
    
    Defined at Math.fpp lines 1056-1062
    
    Parameters
    ----------
    n1 : int
    n2 : int
    n3 : int
    i : int
    j : int
    k : int
    
    Returns
    -------
    dimnu : int
    
    """
    dimnu = _AresMainPy_pkg.f90wrap_thr2mat(n1=n1, n2=n2, n3=n3, i=i, j=j, k=k)
    return dimnu

def mat2thr(n1, n2, n3, i):
    """
    ix, iy, iz = mat2thr(n1, n2, n3, i)
    
    
    Defined at Math.fpp lines 1065-1096
    
    Parameters
    ----------
    n1 : int
    n2 : int
    n3 : int
    i : int
    
    Returns
    -------
    ix : int
    iy : int
    iz : int
    
    ------------------------
    """
    ix, iy, iz = _AresMainPy_pkg.f90wrap_mat2thr(n1=n1, n2=n2, n3=n3, i=i)
    return ix, iy, iz

def sopo(a, lda, n, b, ldb, m, w):
    """
    sopo(a, lda, n, b, ldb, m, w)
    
    
    Defined at Math.fpp lines 1099-1155
    
    Parameters
    ----------
    a : float array
    lda : int
    n : int
    b : float array
    ldb : int
    m : int
    w : float array
    
    """
    _AresMainPy_pkg.f90wrap_sopo(a=a, lda=lda, n=n, b=b, ldb=ldb, m=m, w=w)

def csort_eigen(nev, arr, brr):
    """
    csort_eigen(nev, arr, brr)
    
    
    Defined at Math.fpp lines 1158-1191
    
    Parameters
    ----------
    nev : int
    arr : float array
    brr : complex array
    
    """
    _AresMainPy_pkg.f90wrap_csort_eigen(nev=nev, arr=arr, brr=brr)

def rsort_eigen(nev, arr, brr):
    """
    rsort_eigen(nev, arr, brr)
    
    
    Defined at Math.fpp lines 1194-1227
    
    Parameters
    ----------
    nev : int
    arr : float array
    brr : float array
    
    """
    _AresMainPy_pkg.f90wrap_rsort_eigen(nev=nev, arr=arr, brr=brr)

def sort_eigval(n, arr):
    """
    sort_eigval(n, arr)
    
    
    Defined at Math.fpp lines 1230-1252
    
    Parameters
    ----------
    n : int
    arr : float array
    
    """
    _AresMainPy_pkg.f90wrap_sort_eigval(n=n, arr=arr)

def realint_sort(nev, arr, brr, crr=None):
    """
    realint_sort(nev, arr, brr[, crr])
    
    
    Defined at Math.fpp lines 1255-1293
    
    Parameters
    ----------
    nev : int
    arr : float array
    brr : int array
    crr : float array
    
    """
    _AresMainPy_pkg.f90wrap_realint_sort(nev=nev, arr=arr, brr=brr, crr=crr)

def pgfo(omat, n1, n2, n3, ix, iy, iz):
    """
    ex, ey, ez = pgfo(omat, n1, n2, n3, ix, iy, iz)
    
    
    Defined at Math.fpp lines 1297-1327
    
    Parameters
    ----------
    omat : float array
    n1 : int
    n2 : int
    n3 : int
    ix : int
    iy : int
    iz : int
    
    Returns
    -------
    ex : int
    ey : int
    ez : int
    
    """
    ex, ey, ez = _AresMainPy_pkg.f90wrap_pgfo(omat=omat, n1=n1, n2=n2, n3=n3, ix=ix, \
        iy=iy, iz=iz)
    return ex, ey, ez

def cubicsplineinterp(fun, ddfdx2, xmax, dx, x, zion=None):
    """
    cubicsplineinterp = cubicsplineinterp(fun, ddfdx2, xmax, dx, x[, zion])
    
    
    Defined at Math.fpp lines 1330-1406
    
    Parameters
    ----------
    fun : float array
    ddfdx2 : float array
    xmax : float
    dx : float
    x : float
    zion : float
    
    Returns
    -------
    cubicsplineinterp : float
    
    """
    cubicsplineinterp = _AresMainPy_pkg.f90wrap_cubicsplineinterp(fun=fun, \
        ddfdx2=ddfdx2, xmax=xmax, dx=dx, x=x, zion=zion)
    return cubicsplineinterp

def finite_factor(fnor, norder, coe):
    """
    finite_factor(fnor, norder, coe)
    
    
    Defined at Math.fpp lines 1409-1550
    
    Parameters
    ----------
    fnor : int
    norder : int
    coe : float array
    
    """
    _AresMainPy_pkg.f90wrap_finite_factor(fnor=fnor, norder=norder, coe=coe)

def finite_factor_new(fnor, norder, coe):
    """
    finite_factor_new(fnor, norder, coe)
    
    
    Defined at Math.fpp lines 1552-1608
    
    Parameters
    ----------
    fnor : int
    norder : int
    coe : float array
    
    """
    _AresMainPy_pkg.f90wrap_finite_factor_new(fnor=fnor, norder=norder, coe=coe)

def dfdr(np, h, f, df):
    """
    dfdr(np, h, f, df)
    
    
    Defined at Math.fpp lines 1612-1646
    
    Parameters
    ----------
    np : int
    h : float
    f : float array
    df : float array
    
    """
    _AresMainPy_pkg.f90wrap_dfdr(np=np, h=h, f=f, df=df)

def cubichermiteinterp(fun, dfdx, xmax, h, x):
    """
    cubichermiteinterp = cubichermiteinterp(fun, dfdx, xmax, h, x)
    
    
    Defined at Math.fpp lines 1649-1708
    
    Parameters
    ----------
    fun : float array
    dfdx : float array
    xmax : float
    h : float
    x : float
    
    Returns
    -------
    cubichermiteinterp : float
    
    """
    cubichermiteinterp = _AresMainPy_pkg.f90wrap_cubichermiteinterp(fun=fun, \
        dfdx=dfdx, xmax=xmax, h=h, x=x)
    return cubichermiteinterp

def simpleinterp(fun, xmax, h, x):
    """
    simpleinterp = simpleinterp(fun, xmax, h, x)
    
    
    Defined at Math.fpp lines 1711-1750
    
    Parameters
    ----------
    fun : float array
    xmax : float
    h : float
    x : float
    
    Returns
    -------
    simpleinterp : float
    
    """
    simpleinterp = _AresMainPy_pkg.f90wrap_simpleinterp(fun=fun, xmax=xmax, h=h, \
        x=x)
    return simpleinterp

def r_dylm(l, m, x, y, z, rmod, f):
    """
    r_dylm(l, m, x, y, z, rmod, f)
    
    
    Defined at Math.fpp lines 1753-1893
    
    Parameters
    ----------
    l : int
    m : int
    x : float
    y : float
    z : float
    rmod : float
    f : float array
    
    """
    _AresMainPy_pkg.f90wrap_r_dylm(l=l, m=m, x=x, y=y, z=z, rmod=rmod, f=f)

def atom_effcharge(atom_name, nquan, zeta):
    """
    lmax = atom_effcharge(atom_name, nquan, zeta)
    
    
    Defined at Math.fpp lines 2067-2161
    
    Parameters
    ----------
    atom_name : str
    nquan : int array
    zeta : float array
    
    Returns
    -------
    lmax : int
    
    --------------------------------------------
    """
    lmax = _AresMainPy_pkg.f90wrap_atom_effcharge(atom_name=atom_name, nquan=nquan, \
        zeta=zeta)
    return lmax

def atom_sto(p, l, m, zeta, r):
    """
    f = atom_sto(p, l, m, zeta, r)
    
    
    Defined at Math.fpp lines 2446-2647
    
    Parameters
    ----------
    p : int
    l : int
    m : int
    zeta : float
    r : float array
    
    Returns
    -------
    f : float
    
    """
    f = _AresMainPy_pkg.f90wrap_atom_sto(p=p, l=l, m=m, zeta=zeta, r=r)
    return f

def gammp(a, x):
    """
    gammp = gammp(a, x)
    
    
    Defined at Math.fpp lines 2805-2822
    
    Parameters
    ----------
    a : float
    x : float
    
    Returns
    -------
    gammp : float
    
    """
    gammp = _AresMainPy_pkg.f90wrap_gammp(a=a, x=x)
    return gammp

def gser(gamser, a, x, gln):
    """
    gser(gamser, a, x, gln)
    
    
    Defined at Math.fpp lines 2825-2852
    
    Parameters
    ----------
    gamser : float
    a : float
    x : float
    gln : float
    
    """
    _AresMainPy_pkg.f90wrap_gser(gamser=gamser, a=a, x=x, gln=gln)

def gcf(gammcf, a, x, gln):
    """
    gcf(gammcf, a, x, gln)
    
    
    Defined at Math.fpp lines 2855-2891
    
    Parameters
    ----------
    gammcf : float
    a : float
    x : float
    gln : float
    
    """
    _AresMainPy_pkg.f90wrap_gcf(gammcf=gammcf, a=a, x=x, gln=gln)

def gammln(xx):
    """
    gammln = gammln(xx)
    
    
    Defined at Math.fpp lines 2894-2914
    
    Parameters
    ----------
    xx : float
    
    Returns
    -------
    gammln : float
    
    """
    gammln = _AresMainPy_pkg.f90wrap_gammln(xx=xx)
    return gammln

def integral(l, dx, x):
    """
    integral = integral(l, dx, x)
    
    
    Defined at Math.fpp lines 2917-2930
    
    Parameters
    ----------
    l : int
    dx : float
    x : float
    
    Returns
    -------
    integral : float
    
    """
    integral = _AresMainPy_pkg.f90wrap_integral(l=l, dx=dx, x=x)
    return integral

def myfun(l, t):
    """
    myfun = myfun(l, t)
    
    
    Defined at Math.fpp lines 2933-2938
    
    Parameters
    ----------
    l : int
    t : float
    
    Returns
    -------
    myfun : float
    
    """
    myfun = _AresMainPy_pkg.f90wrap_myfun(l=l, t=t)
    return myfun

def plgndr(l, m, x):
    """
    plgndr = plgndr(l, m, x)
    
    
    Defined at Math.fpp lines 2941-2974
    
    Parameters
    ----------
    l : int
    m : int
    x : float
    
    Returns
    -------
    plgndr : float
    
    """
    plgndr = _AresMainPy_pkg.f90wrap_plgndr(l=l, m=m, x=x)
    return plgndr

def lagrange_interpolation_coe(npoint, scatter_x, coe):
    """
    lagrange_interpolation_coe(npoint, scatter_x, coe)
    
    
    Defined at Math.fpp lines 2977-2994
    
    Parameters
    ----------
    npoint : int
    scatter_x : float array
    coe : float array
    
    """
    _AresMainPy_pkg.f90wrap_lagrange_interpolation_coe(npoint=npoint, \
        scatter_x=scatter_x, coe=coe)

def lagrange_interpolation_x(npoint, x_sample, x_in):
    """
    lagrange_interpolation_x = lagrange_interpolation_x(npoint, x_sample, x_in)
    
    
    Defined at Math.fpp lines 2996-3011
    
    Parameters
    ----------
    npoint : int
    x_sample : float array
    x_in : float
    
    Returns
    -------
    lagrange_interpolation_x : float array
    
    """
    lagrange_interpolation_x = \
        _AresMainPy_pkg.f90wrap_lagrange_interpolation_x(npoint=npoint, \
        x_sample=x_sample, x_in=x_in)
    return lagrange_interpolation_x

def interpolation_test():
    """
    interpolation_test()
    
    
    Defined at Math.fpp lines 3013-3037
    
    
    """
    _AresMainPy_pkg.f90wrap_interpolation_test()

def direct_productlm(nll, nml, index_ll, index_ml, mat_in, mat_out):
    """
    direct_productlm(nll, nml, index_ll, index_ml, mat_in, mat_out)
    
    
    Defined at Math.fpp lines 3039-3116
    
    Parameters
    ----------
    nll : int
    nml : int
    index_ll : int array
    index_ml : int array
    mat_in : float array
    mat_out : float array
    
    """
    _AresMainPy_pkg.f90wrap_direct_productlm(nll=nll, nml=nml, index_ll=index_ll, \
        index_ml=index_ml, mat_in=mat_in, mat_out=mat_out)

def fourier_1d(nr, rr, rab, vr, ll, nql, yp, vql, vt):
    """
    fourier_1d(nr, rr, rab, vr, ll, nql, yp, vql, vt)
    
    
    Defined at Math.fpp lines 3161-3219
    
    Parameters
    ----------
    nr : int
    rr : float array
    rab : float array
    vr : float array
    ll : int
    nql : int
    yp : float array
    vql : float array
    vt : float
    
    """
    _AresMainPy_pkg.f90wrap_fourier_1d(nr=nr, rr=rr, rab=rab, vr=vr, ll=ll, nql=nql, \
        yp=yp, vql=vql, vt=vt)

def invfourier_1d(g, fg, ll, r, fr):
    """
    invfourier_1d(g, fg, ll, r, fr)
    
    
    Defined at Math.fpp lines 3222-3253
    
    Parameters
    ----------
    g : float array
    fg : float array
    ll : int
    r : float array
    fr : float array
    
    """
    _AresMainPy_pkg.f90wrap_invfourier_1d(g=g, fg=fg, ll=ll, r=r, fr=fr)

def integ_new(rab, y):
    """
    f = integ_new(rab, y)
    
    
    Defined at Math.fpp lines 3256-3276
    
    Parameters
    ----------
    rab : float array
    y : float array
    
    Returns
    -------
    f : float
    
    """
    f = _AresMainPy_pkg.f90wrap_integ_new(rab=rab, y=y)
    return f

def interp(np, f, r, rnorm, z=None):
    """
    interp = interp(np, f, r, rnorm[, z])
    
    
    Defined at Math.fpp lines 3279-3333
    
    Parameters
    ----------
    np : int
    f : float array
    r : float array
    rnorm : float
    z : float
    
    Returns
    -------
    interp : float
    
    """
    interp = _AresMainPy_pkg.f90wrap_interp(np=np, f=f, r=r, rnorm=rnorm, z=z)
    return interp

def getfileunit():
    """
    getfileunit = getfileunit()
    
    
    Defined at Math.fpp lines 3408-3421
    
    
    Returns
    -------
    getfileunit : int
    
    """
    getfileunit = _AresMainPy_pkg.f90wrap_getfileunit()
    return getfileunit

def kahan_sum(n, array):
    """
    kahan_sum = kahan_sum(n, array)
    
    
    Defined at Math.fpp lines 3424-3439
    
    Parameters
    ----------
    n : int
    array : float array
    
    Returns
    -------
    kahan_sum : float
    
    """
    kahan_sum = _AresMainPy_pkg.f90wrap_kahan_sum(n=n, array=array)
    return kahan_sum

def _gasdev_s_sp(inmu=None, insigma=None):
    """
    harvest = _gasdev_s_sp([inmu, insigma])
    
    
    Defined at Math.fpp lines 538-567
    
    Parameters
    ----------
    inmu : float
    insigma : float
    
    Returns
    -------
    harvest : float
    
    """
    harvest = _AresMainPy_pkg.f90wrap_gasdev_s_sp(inmu=inmu, insigma=insigma)
    return harvest

def _gasdev_s_dp(inmu=None, insigma=None):
    """
    harvest = _gasdev_s_dp([inmu, insigma])
    
    
    Defined at Math.fpp lines 571-600
    
    Parameters
    ----------
    inmu : float
    insigma : float
    
    Returns
    -------
    harvest : float
    
    """
    harvest = _AresMainPy_pkg.f90wrap_gasdev_s_dp(inmu=inmu, insigma=insigma)
    return harvest

def _gasdev_v_sp(harvest, inmu=None, insigma=None):
    """
    _gasdev_v_sp(harvest[, inmu, insigma])
    
    
    Defined at Math.fpp lines 604-645
    
    Parameters
    ----------
    harvest : float array
    inmu : float
    insigma : float
    
    """
    _AresMainPy_pkg.f90wrap_gasdev_v_sp(harvest=harvest, inmu=inmu, insigma=insigma)

def _gasdev_v_dp(harvest, inmu=None, insigma=None):
    """
    _gasdev_v_dp(harvest[, inmu, insigma])
    
    
    Defined at Math.fpp lines 649-690
    
    Parameters
    ----------
    harvest : float array
    inmu : float
    insigma : float
    
    """
    _AresMainPy_pkg.f90wrap_gasdev_v_dp(harvest=harvest, inmu=inmu, insigma=insigma)

def gasdev(*args, **kwargs):
    """
    gasdev(*args, **kwargs)
    
    
    Defined at Math.fpp lines 21-22
    
    Overloaded interface containing the following procedures:
      _gasdev_s_sp
      _gasdev_s_dp
      _gasdev_v_sp
      _gasdev_v_dp
    
    """
    for proc in [_gasdev_s_sp, _gasdev_s_dp, _gasdev_v_sp, _gasdev_v_dp]:
        try:
            return proc(*args, **kwargs)
        except TypeError:
            continue

def _norm_real(a):
    """
    norm_real = _norm_real(a)
    
    
    Defined at Math.fpp lines 162-164
    
    Parameters
    ----------
    a : float array
    
    Returns
    -------
    norm_real : float
    
    """
    norm_real = _AresMainPy_pkg.f90wrap_norm_real(a=a)
    return norm_real

def _norm_complex(a):
    """
    norm_complex = _norm_complex(a)
    
    
    Defined at Math.fpp lines 168-170
    
    Parameters
    ----------
    a : complex array
    
    Returns
    -------
    norm_complex : float
    
    """
    norm_complex = _AresMainPy_pkg.f90wrap_norm_complex(a=a)
    return norm_complex

def norm(*args, **kwargs):
    """
    norm(*args, **kwargs)
    
    
    Defined at Math.fpp lines 26-27
    
    Overloaded interface containing the following procedures:
      _norm_real
      _norm_complex
    
    """
    for proc in [_norm_real, _norm_complex]:
        try:
            return proc(*args, **kwargs)
        except TypeError:
            continue

def _cross_real(a, b):
    """
    cross_real = _cross_real(a, b)
    
    
    Defined at Math.fpp lines 182-186
    
    Parameters
    ----------
    a : float array
    b : float array
    
    Returns
    -------
    cross_real : float array
    
    """
    cross_real = _AresMainPy_pkg.f90wrap_cross_real(a=a, b=b)
    return cross_real

def _cross_complex(a, b):
    """
    cross_complex = _cross_complex(a, b)
    
    
    Defined at Math.fpp lines 174-178
    
    Parameters
    ----------
    a : complex array
    b : complex array
    
    Returns
    -------
    cross_complex : complex array
    
    """
    cross_complex = _AresMainPy_pkg.f90wrap_cross_complex(a=a, b=b)
    return cross_complex

def cross(*args, **kwargs):
    """
    cross(*args, **kwargs)
    
    
    Defined at Math.fpp lines 31-32
    
    Overloaded interface containing the following procedures:
      _cross_real
      _cross_complex
    
    """
    for proc in [_cross_real, _cross_complex]:
        try:
            return proc(*args, **kwargs)
        except TypeError:
            continue

def _integer_index(array, n, id):
    """
    _integer_index(array, n, id)
    
    
    Defined at Math.fpp lines 425-466
    
    Parameters
    ----------
    array : int array
    n : int
    id : int array
    
    """
    _AresMainPy_pkg.f90wrap_integer_index(array=array, n=n, id=id)

def _real_index(array, n, id):
    """
    _real_index(array, n, id)
    
    
    Defined at Math.fpp lines 380-421
    
    Parameters
    ----------
    array : float array
    n : int
    id : int array
    
    """
    _AresMainPy_pkg.f90wrap_real_index(array=array, n=n, id=id)

def sort_id(*args, **kwargs):
    """
    sort_id(*args, **kwargs)
    
    
    Defined at Math.fpp lines 36-37
    
    Overloaded interface containing the following procedures:
      _integer_index
      _real_index
    
    """
    for proc in [_integer_index, _real_index]:
        try:
            return proc(*args, **kwargs)
        except TypeError:
            continue

def _three2one_d_cplx(amat, bmat):
    """
    _three2one_d_cplx(amat, bmat)
    
    
    Defined at Math.fpp lines 3337-3351
    
    Parameters
    ----------
    amat : complex array
    bmat : complex array
    
    -------------------------------------------------------------
    allocate(bmat(grid%n))
    """
    _AresMainPy_pkg.f90wrap_three2one_d_cplx(amat=amat, bmat=bmat)

def _three2one_d_real(amat, bmat):
    """
    _three2one_d_real(amat, bmat)
    
    
    Defined at Math.fpp lines 3355-3369
    
    Parameters
    ----------
    amat : float array
    bmat : float array
    
    -------------------------------------------------------------
    allocate(bmat(grid%n))
    """
    _AresMainPy_pkg.f90wrap_three2one_d_real(amat=amat, bmat=bmat)

def three2one_dim(*args, **kwargs):
    """
    three2one_dim(*args, **kwargs)
    
    
    Defined at Math.fpp lines 42-44
    
    Overloaded interface containing the following procedures:
      _three2one_d_cplx
      _three2one_d_real
    
    """
    for proc in [_three2one_d_cplx, _three2one_d_real]:
        try:
            return proc(*args, **kwargs)
        except TypeError:
            continue

def _one2three_d_cplx(amat, bmat):
    """
    _one2three_d_cplx(amat, bmat)
    
    
    Defined at Math.fpp lines 3373-3387
    
    Parameters
    ----------
    amat : complex array
    bmat : complex array
    
    -------------------------------------------------------------
    allocate(bmat(grid%n1,grid%n2,grid%n3))
    """
    _AresMainPy_pkg.f90wrap_one2three_d_cplx(amat=amat, bmat=bmat)

def _one2three_d_real(amat, bmat):
    """
    _one2three_d_real(amat, bmat)
    
    
    Defined at Math.fpp lines 3391-3405
    
    Parameters
    ----------
    amat : float array
    bmat : float array
    
    -------------------------------------------------------------
    allocate(bmat(grid%n1,grid%n2,grid%n3))
    """
    _AresMainPy_pkg.f90wrap_one2three_d_real(amat=amat, bmat=bmat)

def one2three_dim(*args, **kwargs):
    """
    one2three_dim(*args, **kwargs)
    
    
    Defined at Math.fpp lines 46-47
    
    Overloaded interface containing the following procedures:
      _one2three_d_cplx
      _one2three_d_real
    
    """
    for proc in [_one2three_d_cplx, _one2three_d_real]:
        try:
            return proc(*args, **kwargs)
        except TypeError:
            continue


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "math".')

for func in _dt_array_initialisers:
    func()

"""
Module math


Defined at Math.fpp lines 5-2413

"""
from __future__ import print_function, absolute_import, division
import _arespy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def change_case(instr, str, fun):
    """
    change_case(instr, str, fun)
    
    
    Defined at Math.fpp lines 42-82
    
    Parameters
    ----------
    instr : str
    str : str
    fun : int
    
    """
    _arespy_pkg.f90wrap_change_case(instr=instr, str=str, fun=fun)

def find_keywords(str, ch_mark, id_key, id_value):
    """
    find_keywords(str, ch_mark, id_key, id_value)
    
    
    Defined at Math.fpp lines 86-115
    
    Parameters
    ----------
    str : str
    ch_mark : str
    id_key : int
    id_value : int
    
    """
    _arespy_pkg.f90wrap_find_keywords(str=str, ch_mark=ch_mark, id_key=id_key, \
        id_value=id_value)

def find_nword(str, ch_comma, nword):
    """
    find_nword(str, ch_comma, nword)
    
    
    Defined at Math.fpp lines 119-145
    
    Parameters
    ----------
    str : str
    ch_comma : str
    nword : int
    
    """
    _arespy_pkg.f90wrap_find_nword(str=str, ch_comma=ch_comma, nword=nword)

def det(matrix):
    """
    det = det(matrix)
    
    
    Defined at Math.fpp lines 177-181
    
    Parameters
    ----------
    matrix : float array
    
    Returns
    -------
    det : float
    
    """
    det = _arespy_pkg.f90wrap_det(matrix=matrix)
    return det

def inv_33(m):
    """
    inv_33 = inv_33(m)
    
    
    Defined at Math.fpp lines 185-215
    
    Parameters
    ----------
    m : float array
    
    Returns
    -------
    inv_33 : float array
    
    """
    inv_33 = _arespy_pkg.f90wrap_inv_33(m=m)
    return inv_33

def lindg(eta, lambda_, mu):
    """
    lindg = lindg(eta, lambda_, mu)
    
    
    Defined at Math.fpp lines 219-265
    
    Parameters
    ----------
    eta : float
    lambda_ : float
    mu : float
    
    Returns
    -------
    lindg : float
    
    """
    lindg = _arespy_pkg.f90wrap_lindg(eta=eta, lambda_=lambda_, mu=mu)
    return lindg

def int_to_char(int_bn):
    """
    int_to_char = int_to_char(int_bn)
    
    
    Defined at Math.fpp lines 269-290
    
    Parameters
    ----------
    int_bn : int
    
    Returns
    -------
    int_to_char : str
    
    -----------------------------------------------------------------------
    """
    int_to_char = _arespy_pkg.f90wrap_int_to_char(int_bn=int_bn)
    return int_to_char

def lat2matrix(lat_para, lat_mat, flag):
    """
    lat2matrix(lat_para, lat_mat, flag)
    
    
    Defined at Math.fpp lines 294-337
    
    Parameters
    ----------
    lat_para : float array
    lat_mat : float array
    flag : int
    
    """
    _arespy_pkg.f90wrap_lat2matrix(lat_para=lat_para, lat_mat=lat_mat, flag=flag)

def one2three(id, n_dens, pos):
    """
    one2three(id, n_dens, pos)
    
    
    Defined at Math.fpp lines 341-363
    
    Parameters
    ----------
    id : int
    n_dens : int array
    pos : int array
    
    """
    _arespy_pkg.f90wrap_one2three(id=id, n_dens=n_dens, pos=pos)

def init_random_seed():
    """
    init_random_seed()
    
    
    Defined at Math.fpp lines 457-521
    
    
    """
    _arespy_pkg.f90wrap_init_random_seed()

def atom_mass(atom_name, mass):
    """
    atom_mass(atom_name, mass)
    
    
    Defined at Math.fpp lines 529-670
    
    Parameters
    ----------
    atom_name : str
    mass : float
    
    --------------------------------------------
    """
    _arespy_pkg.f90wrap_atom_mass(atom_name=atom_name, mass=mass)

def newton_inter(n, x, y, m, tx, ty):
    """
    newton_inter(n, x, y, m, tx, ty)
    
    
    Defined at Math.fpp lines 674-697
    
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
    _arespy_pkg.f90wrap_newton_inter(n=n, x=x, y=y, m=m, tx=tx, ty=ty)

def diag(n, ina, w, q):
    """
    diag(n, ina, w, q)
    
    
    Defined at Math.fpp lines 701-835
    
    Parameters
    ----------
    n : int
    ina : float array
    w : float array
    q : float array
    
    """
    _arespy_pkg.f90wrap_diag(n=n, ina=ina, w=w, q=q)

def boltzmann_distribution(rnull, width):
    """
    boltzmann_distribution = boltzmann_distribution(rnull, width)
    
    
    Defined at Math.fpp lines 839-854
    
    Parameters
    ----------
    rnull : float
    width : float
    
    Returns
    -------
    boltzmann_distribution : float
    
    """
    boltzmann_distribution = _arespy_pkg.f90wrap_boltzmann_distribution(rnull=rnull, \
        width=width)
    return boltzmann_distribution

def dir2car(cry_coo, ort_coo, lat):
    """
    dir2car(cry_coo, ort_coo, lat)
    
    
    Defined at Math.fpp lines 861-871
    
    Parameters
    ----------
    cry_coo : float array
    ort_coo : float array
    lat : float array
    
    """
    _arespy_pkg.f90wrap_dir2car(cry_coo=cry_coo, ort_coo=ort_coo, lat=lat)

def car2dir(ort_coo, cry_coo, lat):
    """
    car2dir(ort_coo, cry_coo, lat)
    
    
    Defined at Math.fpp lines 877-888
    
    Parameters
    ----------
    ort_coo : float array
    cry_coo : float array
    lat : float array
    
    """
    _arespy_pkg.f90wrap_car2dir(ort_coo=ort_coo, cry_coo=cry_coo, lat=lat)

def thr2mat(n1, n2, n3, i, j, k):
    """
    dimnu = thr2mat(n1, n2, n3, i, j, k)
    
    
    Defined at Math.fpp lines 891-897
    
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
    dimnu = _arespy_pkg.f90wrap_thr2mat(n1=n1, n2=n2, n3=n3, i=i, j=j, k=k)
    return dimnu

def mat2thr(n1, n2, n3, i, offset=None):
    """
    ix, iy, iz = mat2thr(n1, n2, n3, i[, offset])
    
    
    Defined at Math.fpp lines 900-938
    
    Parameters
    ----------
    n1 : int
    n2 : int
    n3 : int
    i : int
    offset : int array
    
    Returns
    -------
    ix : int
    iy : int
    iz : int
    
    ------------------------
    """
    ix, iy, iz = _arespy_pkg.f90wrap_mat2thr(n1=n1, n2=n2, n3=n3, i=i, \
        offset=offset)
    return ix, iy, iz

def sopo(a, lda, n, b, ldb, m, w):
    """
    sopo(a, lda, n, b, ldb, m, w)
    
    
    Defined at Math.fpp lines 941-996
    
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
    _arespy_pkg.f90wrap_sopo(a=a, lda=lda, n=n, b=b, ldb=ldb, m=m, w=w)

def csort_eigen(nev, arr, brr):
    """
    csort_eigen(nev, arr, brr)
    
    
    Defined at Math.fpp lines 999-1037
    
    Parameters
    ----------
    nev : int
    arr : float array
    brr : complex array
    
    """
    _arespy_pkg.f90wrap_csort_eigen(nev=nev, arr=arr, brr=brr)

def rsort_eigen(nev, arr, brr):
    """
    rsort_eigen(nev, arr, brr)
    
    
    Defined at Math.fpp lines 1040-1074
    
    Parameters
    ----------
    nev : int
    arr : float array
    brr : float array
    
    """
    _arespy_pkg.f90wrap_rsort_eigen(nev=nev, arr=arr, brr=brr)

def realint_sort(nev, arr, brr, crr=None):
    """
    realint_sort(nev, arr, brr[, crr])
    
    
    Defined at Math.fpp lines 1077-1116
    
    Parameters
    ----------
    nev : int
    arr : float array
    brr : int array
    crr : float array
    
    """
    _arespy_pkg.f90wrap_realint_sort(nev=nev, arr=arr, brr=brr, crr=crr)

def sort_eigval(n, arr):
    """
    sort_eigval(n, arr)
    
    
    Defined at Math.fpp lines 1119-1141
    
    Parameters
    ----------
    n : int
    arr : float array
    
    """
    _arespy_pkg.f90wrap_sort_eigval(n=n, arr=arr)

def pgfo(omat, n1, n2, n3, ix, iy, iz):
    """
    ex, ey, ez = pgfo(omat, n1, n2, n3, ix, iy, iz)
    
    
    Defined at Math.fpp lines 1144-1174
    
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
    ex, ey, ez = _arespy_pkg.f90wrap_pgfo(omat=omat, n1=n1, n2=n2, n3=n3, ix=ix, \
        iy=iy, iz=iz)
    return ex, ey, ez

def cubicsplineinterp(fun, ddfdx2, xmax, dx, x, zion=None):
    """
    cubicsplineinterp = cubicsplineinterp(fun, ddfdx2, xmax, dx, x[, zion])
    
    
    Defined at Math.fpp lines 1177-1253
    
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
    cubicsplineinterp = _arespy_pkg.f90wrap_cubicsplineinterp(fun=fun, \
        ddfdx2=ddfdx2, xmax=xmax, dx=dx, x=x, zion=zion)
    return cubicsplineinterp

def finite_factor(fnor, norder, coe):
    """
    finite_factor(fnor, norder, coe)
    
    
    Defined at Math.fpp lines 1256-1397
    
    Parameters
    ----------
    fnor : int
    norder : int
    coe : float array
    
    """
    _arespy_pkg.f90wrap_finite_factor(fnor=fnor, norder=norder, coe=coe)

def finite_factor_new(fnor, norder, coe):
    """
    finite_factor_new(fnor, norder, coe)
    
    
    Defined at Math.fpp lines 1400-1457
    
    Parameters
    ----------
    fnor : int
    norder : int
    coe : float array
    
    """
    _arespy_pkg.f90wrap_finite_factor_new(fnor=fnor, norder=norder, coe=coe)

def dfdr(np, h, f, df):
    """
    dfdr(np, h, f, df)
    
    
    Defined at Math.fpp lines 1460-1494
    
    Parameters
    ----------
    np : int
    h : float
    f : float array
    df : float array
    
    """
    _arespy_pkg.f90wrap_dfdr(np=np, h=h, f=f, df=df)

def cubichermiteinterp(fun, dfdx, xmax, h, x):
    """
    cubichermiteinterp = cubichermiteinterp(fun, dfdx, xmax, h, x)
    
    
    Defined at Math.fpp lines 1497-1556
    
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
    cubichermiteinterp = _arespy_pkg.f90wrap_cubichermiteinterp(fun=fun, dfdx=dfdx, \
        xmax=xmax, h=h, x=x)
    return cubichermiteinterp

def simpleinterp(fun, xmax, h, x):
    """
    simpleinterp = simpleinterp(fun, xmax, h, x)
    
    
    Defined at Math.fpp lines 1559-1598
    
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
    simpleinterp = _arespy_pkg.f90wrap_simpleinterp(fun=fun, xmax=xmax, h=h, x=x)
    return simpleinterp

def r_dylm(l, m, x, y, z, rmod, f):
    """
    r_dylm(l, m, x, y, z, rmod, f)
    
    
    Defined at Math.fpp lines 1601-1741
    
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
    _arespy_pkg.f90wrap_r_dylm(l=l, m=m, x=x, y=y, z=z, rmod=rmod, f=f)

def rcsntable(n1, n2, n3, dn, h, srmax, rcs, cns, sphindx):
    """
    nspt = rcsntable(n1, n2, n3, dn, h, srmax, rcs, cns, sphindx)
    
    
    Defined at Math.fpp lines 1744-1803
    
    Parameters
    ----------
    n1 : int
    n2 : int
    n3 : int
    dn : int
    h : float
    srmax : float
    rcs : float array
    cns : complex array
    sphindx : bool array
    
    Returns
    -------
    nspt : int
    
    """
    nspt = _arespy_pkg.f90wrap_rcsntable(n1=n1, n2=n2, n3=n3, dn=dn, h=h, \
        srmax=srmax, rcs=rcs, cns=cns, sphindx=sphindx)
    return nspt

def rcsntable_atoms(n1, n2, n3, dn, na, h, poscar, srmax, atomr, rvec, rcs, cns, \
    sphindx):
    """
    nspt = rcsntable_atoms(n1, n2, n3, dn, na, h, poscar, srmax, atomr, rvec, rcs, \
        cns, sphindx)
    
    
    Defined at Math.fpp lines 1806-1885
    
    Parameters
    ----------
    n1 : int
    n2 : int
    n3 : int
    dn : int
    na : int
    h : float
    poscar : float array
    srmax : float
    atomr : float
    rvec : float array
    rcs : float array
    cns : complex array
    sphindx : bool array
    
    Returns
    -------
    nspt : int
    
    """
    nspt = _arespy_pkg.f90wrap_rcsntable_atoms(n1=n1, n2=n2, n3=n3, dn=dn, na=na, \
        h=h, poscar=poscar, srmax=srmax, atomr=atomr, rvec=rvec, rcs=rcs, cns=cns, \
        sphindx=sphindx)
    return nspt

def car2spe(orig, h, x, y, z):
    """
    r, cost, sint, cosp, sinp = car2spe(orig, h, x, y, z)
    
    
    Defined at Math.fpp lines 1888-1919
    
    Parameters
    ----------
    orig : float array
    h : float
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
    r, cost, sint, cosp, sinp = _arespy_pkg.f90wrap_car2spe(orig=orig, h=h, x=x, \
        y=y, z=z)
    return r, cost, sint, cosp, sinp

def calclm(lmax, clm):
    """
    calclm(lmax, clm)
    
    
    Defined at Math.fpp lines 1922-1935
    
    Parameters
    ----------
    lmax : int
    clm : float array
    
    """
    _arespy_pkg.f90wrap_calclm(lmax=lmax, clm=clm)

def c(l, m):
    """
    c = c(l, m)
    
    
    Defined at Math.fpp lines 1938-1951
    
    Parameters
    ----------
    l : int
    m : int
    
    Returns
    -------
    c : float
    
    """
    c = _arespy_pkg.f90wrap_c(l=l, m=m)
    return c

def cal_plm(lmax, x, sx, plm):
    """
    cal_plm(lmax, x, sx, plm)
    
    
    Defined at Math.fpp lines 1954-1981
    
    Parameters
    ----------
    lmax : int
    x : float
    sx : float
    plm : float array
    
    """
    _arespy_pkg.f90wrap_cal_plm(lmax=lmax, x=x, sx=sx, plm=plm)

def interp(np, f, r, rnorm, z=None):
    """
    interp = interp(np, f, r, rnorm[, z])
    
    
    Defined at Math.fpp lines 1984-2037
    
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
    interp = _arespy_pkg.f90wrap_interp(np=np, f=f, r=r, rnorm=rnorm, z=z)
    return interp

def interp_dnf(n, np, f, r, rnorm, z=None):
    """
    interp_dnf = interp_dnf(n, np, f, r, rnorm[, z])
    
    
    Defined at Math.fpp lines 2040-2101
    
    Parameters
    ----------
    n : int
    np : int
    f : float array
    r : float array
    rnorm : float
    z : float
    
    Returns
    -------
    interp_dnf : float
    
    """
    interp_dnf = _arespy_pkg.f90wrap_interp_dnf(n=n, np=np, f=f, r=r, rnorm=rnorm, \
        z=z)
    return interp_dnf

def polynom(m, np, xa, ya, c, x):
    """
    polynom = polynom(m, np, xa, ya, c, x)
    
    
    Defined at Math.fpp lines 2104-2284
    
    Parameters
    ----------
    m : int
    np : int
    xa : float array
    ya : float array
    c : float array
    x : float
    
    Returns
    -------
    polynom : float
    
    """
    polynom = _arespy_pkg.f90wrap_polynom(m=m, np=np, xa=xa, ya=ya, c=c, x=x)
    return polynom

def fourier_1d(nr, rr, rab, vr, ll, nql, yp, vql, vt):
    """
    fourier_1d(nr, rr, rab, vr, ll, nql, yp, vql, vt)
    
    
    Defined at Math.fpp lines 2287-2345
    
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
    _arespy_pkg.f90wrap_fourier_1d(nr=nr, rr=rr, rab=rab, vr=vr, ll=ll, nql=nql, \
        yp=yp, vql=vql, vt=vt)

def integ_new(rab, y):
    """
    f = integ_new(rab, y)
    
    
    Defined at Math.fpp lines 2348-2368
    
    Parameters
    ----------
    rab : float array
    y : float array
    
    Returns
    -------
    f : float
    
    """
    f = _arespy_pkg.f90wrap_integ_new(rab=rab, y=y)
    return f

def sphbess(l, x):
    """
    sphbess = sphbess(l, x)
    
    
    Defined at Math.fpp lines 2371-2410
    
    Parameters
    ----------
    l : int
    x : float
    
    Returns
    -------
    sphbess : float
    
    """
    sphbess = _arespy_pkg.f90wrap_sphbess(l=l, x=x)
    return sphbess

def _norm_real(a):
    """
    norm_real = _norm_real(a)
    
    
    Defined at Math.fpp lines 149-151
    
    Parameters
    ----------
    a : float array
    
    Returns
    -------
    norm_real : float
    
    """
    norm_real = _arespy_pkg.f90wrap_norm_real(a=a)
    return norm_real

def _norm_complex(a):
    """
    norm_complex = _norm_complex(a)
    
    
    Defined at Math.fpp lines 155-157
    
    Parameters
    ----------
    a : complex array
    
    Returns
    -------
    norm_complex : float
    
    """
    norm_complex = _arespy_pkg.f90wrap_norm_complex(a=a)
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
    
    
    Defined at Math.fpp lines 169-173
    
    Parameters
    ----------
    a : float array
    b : float array
    
    Returns
    -------
    cross_real : float array
    
    """
    cross_real = _arespy_pkg.f90wrap_cross_real(a=a, b=b)
    return cross_real

def _cross_complex(a, b):
    """
    cross_complex = _cross_complex(a, b)
    
    
    Defined at Math.fpp lines 161-165
    
    Parameters
    ----------
    a : complex array
    b : complex array
    
    Returns
    -------
    cross_complex : complex array
    
    """
    cross_complex = _arespy_pkg.f90wrap_cross_complex(a=a, b=b)
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
    
    
    Defined at Math.fpp lines 412-453
    
    Parameters
    ----------
    array : int array
    n : int
    id : int array
    
    """
    _arespy_pkg.f90wrap_integer_index(array=array, n=n, id=id)

def _real_index(array, n, id):
    """
    _real_index(array, n, id)
    
    
    Defined at Math.fpp lines 367-408
    
    Parameters
    ----------
    array : float array
    n : int
    id : int array
    
    """
    _arespy_pkg.f90wrap_real_index(array=array, n=n, id=id)

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


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "math".')

for func in _dt_array_initialisers:
    func()

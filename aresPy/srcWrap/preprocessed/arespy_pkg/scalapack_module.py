"""
Module scalapack_module


Defined at ScaLapack_module.fpp lines 5-1152

"""
from __future__ import print_function, absolute_import, division
import _arespy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def init_scala():
    """
    init_scala()
    
    
    Defined at ScaLapack_module.fpp lines 29-49
    
    
    """
    _arespy_pkg.f90wrap_init_scala()

def sl_orthnorm(amat, m, n, mb, nb):
    """
    sl_orthnorm(amat, m, n, mb, nb)
    
    
    Defined at ScaLapack_module.fpp lines 61-80
    
    Parameters
    ----------
    amat : complex array
    m : int
    n : int
    mb : int
    nb : int
    
    """
    _arespy_pkg.f90wrap_sl_orthnorm(amat=amat, m=m, n=n, mb=mb, nb=nb)

def sl_orthnorm_real(amat, m, n, mb, nb):
    """
    sl_orthnorm_real(amat, m, n, mb, nb)
    
    
    Defined at ScaLapack_module.fpp lines 83-118
    
    Parameters
    ----------
    amat : float array
    m : int
    n : int
    mb : int
    nb : int
    
    """
    _arespy_pkg.f90wrap_sl_orthnorm_real(amat=amat, m=m, n=n, mb=mb, nb=nb)

def sl_matmat(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb, \
    cmbin=None, cnbin=None):
    """
    sl_matmat(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb[, \
        cmbin, cnbin])
    
    
    Defined at ScaLapack_module.fpp lines 121-177
    
    Parameters
    ----------
    opa : str
    opb : str
    amat : complex array
    bmat : complex array
    cmat : complex array
    am : int
    an : int
    bm : int
    bn : int
    amb : int
    anb : int
    bmb : int
    bnb : int
    cmbin : int
    cnbin : int
    
    ---------------------------------------------------------------------
    """
    _arespy_pkg.f90wrap_sl_matmat(opa=opa, opb=opb, amat=amat, bmat=bmat, cmat=cmat, \
        am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, cmbin=cmbin, \
        cnbin=cnbin)

def sl_matmat_cmplx_cn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, \
    bmb, bnb, cmb, cnb):
    """
    sl_matmat_cmplx_cn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, \
        bnb, cmb, cnb)
    
    
    Defined at ScaLapack_module.fpp lines 180-289
    
    Parameters
    ----------
    opa : str
    opb : str
    amat : complex array
    bmat : complex array
    cmat : complex array
    am : int
    an : int
    bm : int
    bn : int
    amb : int
    anb : int
    bmb : int
    bnb : int
    cmb : int
    cnb : int
    
    ---------------------------------------------------------------------
    """
    _arespy_pkg.f90wrap_sl_matmat_cmplx_cn(opa=opa, opb=opb, amat=amat, bmat=bmat, \
        cmat=cmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, \
        cmb=cmb, cnb=cnb)

def sl_matmat_cmplx_nn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, \
    bmb, bnb, cmb, cnb):
    """
    sl_matmat_cmplx_nn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, \
        bnb, cmb, cnb)
    
    
    Defined at ScaLapack_module.fpp lines 292-351
    
    Parameters
    ----------
    opa : str
    opb : str
    amat : complex array
    bmat : complex array
    cmat : complex array
    am : int
    an : int
    bm : int
    bn : int
    amb : int
    anb : int
    bmb : int
    bnb : int
    cmb : int
    cnb : int
    
    ---------------------------------------------------------------------
    > cmat_local
    """
    _arespy_pkg.f90wrap_sl_matmat_cmplx_nn(opa=opa, opb=opb, amat=amat, bmat=bmat, \
        cmat=cmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, \
        cmb=cmb, cnb=cnb)

def sl_matmat_sub(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, \
    bnb, cmbin=None, cnbin=None):
    """
    sl_matmat_sub(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb[, \
        cmbin, cnbin])
    
    
    Defined at ScaLapack_module.fpp lines 354-439
    
    Parameters
    ----------
    opa : str
    opb : str
    amat : complex array
    bmat : complex array
    cmat : complex array
    am : int
    an : int
    bm : int
    bn : int
    amb : int
    anb : int
    bmb : int
    bnb : int
    cmbin : int
    cnbin : int
    
    ---------------------------------------------------------------------
    blacs_contxt = parallel%comm
    nprow = parallel%numprocs
    npcol = 1
    CALL BLACS_GET( -1, 0, blacs_contxt )
    CALL BLACS_GRIDINIT( blacs_contxt, 'Row-major', NPROW, NPCOL )
    CALL BLACS_GRIDINFO( blacs_contxt, NPROW, NPCOL, MYROW, MYCOL )
    call blacs_pinfo(iam,nprocs)
    am = grid%tn
    bm = grid%tn
    if( .not. init_called) then
     call init_scala_sub()
    end if
    """
    _arespy_pkg.f90wrap_sl_matmat_sub(opa=opa, opb=opb, amat=amat, bmat=bmat, \
        cmat=cmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, \
        cmbin=cmbin, cnbin=cnbin)

def sl_matmat_sub_real(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, \
    bmb, bnb, cmbin=None, cnbin=None):
    """
    sl_matmat_sub_real(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, \
        bnb[, cmbin, cnbin])
    
    
    Defined at ScaLapack_module.fpp lines 442-527
    
    Parameters
    ----------
    opa : str
    opb : str
    amat : float array
    bmat : float array
    cmat : float array
    am : int
    an : int
    bm : int
    bn : int
    amb : int
    anb : int
    bmb : int
    bnb : int
    cmbin : int
    cnbin : int
    
    ---------------------------------------------------------------------
    blacs_contxt = parallel%comm
    nprow = parallel%numprocs
    npcol = 1
    CALL BLACS_GET( -1, 0, blacs_contxt )
    CALL BLACS_GRIDINIT( blacs_contxt, 'Row-major', NPROW, NPCOL )
    CALL BLACS_GRIDINFO( blacs_contxt, NPROW, NPCOL, MYROW, MYCOL )
    call blacs_pinfo(iam,nprocs)
    am = grid%tn
    bm = grid%tn
    if( .not. init_called) then
     call init_scala_sub()
    end if
    """
    _arespy_pkg.f90wrap_sl_matmat_sub_real(opa=opa, opb=opb, amat=amat, bmat=bmat, \
        cmat=cmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, \
        cmbin=cmbin, cnbin=cnbin)

def sl_matmat_real(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, \
    bnb, cmbin=None, cnbin=None):
    """
    sl_matmat_real(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb[, \
        cmbin, cnbin])
    
    
    Defined at ScaLapack_module.fpp lines 530-627
    
    Parameters
    ----------
    opa : str
    opb : str
    amat : float array
    bmat : float array
    cmat : float array
    am : int
    an : int
    bm : int
    bn : int
    amb : int
    anb : int
    bmb : int
    bnb : int
    cmbin : int
    cnbin : int
    
    ---------------------------------------------------------------------
    """
    _arespy_pkg.f90wrap_sl_matmat_real(opa=opa, opb=opb, amat=amat, bmat=bmat, \
        cmat=cmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, \
        cmbin=cmbin, cnbin=cnbin)

def sl_generalizeeigen(dime, amat, bmat, am, an, bm, bn, amb, anb, bmb, bnb, \
    evec, cm, cn, eval, cmbin=None, cnbin=None):
    """
    sl_generalizeeigen(dime, amat, bmat, am, an, bm, bn, amb, anb, bmb, bnb, evec, \
        cm, cn, eval[, cmbin, cnbin])
    
    
    Defined at ScaLapack_module.fpp lines 630-713
    
    Parameters
    ----------
    dime : int
    amat : complex array
    bmat : complex array
    am : int
    an : int
    bm : int
    bn : int
    amb : int
    anb : int
    bmb : int
    bnb : int
    evec : complex array
    cm : int
    cn : int
    eval : float array
    cmbin : int
    cnbin : int
    
    """
    _arespy_pkg.f90wrap_sl_generalizeeigen(dime=dime, amat=amat, bmat=bmat, am=am, \
        an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, evec=evec, cm=cm, \
        cn=cn, eval=eval, cmbin=cmbin, cnbin=cnbin)

def sl_generalizeeigen_real(dime, amat, bmat, am, an, bm, bn, amb, anb, bmb, \
    bnb, evec, cm, cn, eval, cmbin=None, cnbin=None):
    """
    sl_generalizeeigen_real(dime, amat, bmat, am, an, bm, bn, amb, anb, bmb, bnb, \
        evec, cm, cn, eval[, cmbin, cnbin])
    
    
    Defined at ScaLapack_module.fpp lines 716-861
    
    Parameters
    ----------
    dime : int
    amat : float array
    bmat : float array
    am : int
    an : int
    bm : int
    bn : int
    amb : int
    anb : int
    bmb : int
    bnb : int
    evec : float array
    cm : int
    cn : int
    eval : float array
    cmbin : int
    cnbin : int
    
    """
    _arespy_pkg.f90wrap_sl_generalizeeigen_real(dime=dime, amat=amat, bmat=bmat, \
        am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, evec=evec, \
        cm=cm, cn=cn, eval=eval, cmbin=cmbin, cnbin=cnbin)

def sl_diagm_real(dime, amat, am, an, bm, bn, amb, anb, bmb, bnb, evec, cm, cn, \
    eval, cmb, cnb):
    """
    sl_diagm_real(dime, amat, am, an, bm, bn, amb, anb, bmb, bnb, evec, cm, cn, \
        eval, cmb, cnb)
    
    
    Defined at ScaLapack_module.fpp lines 864-947
    
    Parameters
    ----------
    dime : int
    amat : float array
    am : int
    an : int
    bm : int
    bn : int
    amb : int
    anb : int
    bmb : int
    bnb : int
    evec : float array
    cm : int
    cn : int
    eval : float array
    cmb : int
    cnb : int
    
    """
    _arespy_pkg.f90wrap_sl_diagm_real(dime=dime, amat=amat, am=am, an=an, bm=bm, \
        bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, evec=evec, cm=cm, cn=cn, \
        eval=eval, cmb=cmb, cnb=cnb)

def twod_map_set(nstates, nrow, ncol, twod_map):
    """
    twod_map_set(nstates, nrow, ncol, twod_map)
    
    
    Defined at ScaLapack_module.fpp lines 951-969
    
    Parameters
    ----------
    nstates : int
    nrow : int
    ncol : int
    twod_map : int array
    
    """
    _arespy_pkg.f90wrap_twod_map_set(nstates=nstates, nrow=nrow, ncol=ncol, \
        twod_map=twod_map)

def sl_matmat_real_tn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, \
    bnb, cmb, cnb):
    """
    sl_matmat_real_tn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, \
        bnb, cmb, cnb)
    
    
    Defined at ScaLapack_module.fpp lines 971-1089
    
    Parameters
    ----------
    opa : str
    opb : str
    amat : float array
    bmat : float array
    cmat : float array
    am : int
    an : int
    bm : int
    bn : int
    amb : int
    anb : int
    bmb : int
    bnb : int
    cmb : int
    cnb : int
    
    ---------------------------------------------------------------------
    """
    _arespy_pkg.f90wrap_sl_matmat_real_tn(opa=opa, opb=opb, amat=amat, bmat=bmat, \
        cmat=cmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, \
        cmb=cmb, cnb=cnb)

def sl_matmat_real_nn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, \
    bnb, cmb, cnb):
    """
    sl_matmat_real_nn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, \
        bnb, cmb, cnb)
    
    
    Defined at ScaLapack_module.fpp lines 1092-1151
    
    Parameters
    ----------
    opa : str
    opb : str
    amat : float array
    bmat : float array
    cmat : float array
    am : int
    an : int
    bm : int
    bn : int
    amb : int
    anb : int
    bmb : int
    bnb : int
    cmb : int
    cnb : int
    
    ---------------------------------------------------------------------
    > cmat_local
    """
    _arespy_pkg.f90wrap_sl_matmat_real_nn(opa=opa, opb=opb, amat=amat, bmat=bmat, \
        cmat=cmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, \
        cmb=cmb, cnb=cnb)

def get_blacs_contxt():
    """
    Element blacs_contxt ftype=integer(i4b) pytype=int
    
    
    Defined at ScaLapack_module.fpp line 15
    
    """
    return _arespy_pkg.f90wrap_scalapack_module__get__blacs_contxt()

def set_blacs_contxt(blacs_contxt):
    _arespy_pkg.f90wrap_scalapack_module__set__blacs_contxt(blacs_contxt)

def get_dlen():
    """
    Element dlen ftype=integer(i4b) pytype=int
    
    
    Defined at ScaLapack_module.fpp line 16
    
    """
    return _arespy_pkg.f90wrap_scalapack_module__get__dlen()

DLEN = get_dlen()

def get_myrow():
    """
    Element myrow ftype=integer(i4b) pytype=int
    
    
    Defined at ScaLapack_module.fpp line 17
    
    """
    return _arespy_pkg.f90wrap_scalapack_module__get__myrow()

def set_myrow(myrow):
    _arespy_pkg.f90wrap_scalapack_module__set__myrow(myrow)

def get_mycol():
    """
    Element mycol ftype=integer(i4b) pytype=int
    
    
    Defined at ScaLapack_module.fpp line 17
    
    """
    return _arespy_pkg.f90wrap_scalapack_module__get__mycol()

def set_mycol(mycol):
    _arespy_pkg.f90wrap_scalapack_module__set__mycol(mycol)

def get_npcol():
    """
    Element npcol ftype=integer(i4b) pytype=int
    
    
    Defined at ScaLapack_module.fpp line 18
    
    """
    return _arespy_pkg.f90wrap_scalapack_module__get__npcol()

def set_npcol(npcol):
    _arespy_pkg.f90wrap_scalapack_module__set__npcol(npcol)

def get_nprow():
    """
    Element nprow ftype=integer(i4b) pytype=int
    
    
    Defined at ScaLapack_module.fpp line 18
    
    """
    return _arespy_pkg.f90wrap_scalapack_module__get__nprow()

def set_nprow(nprow):
    _arespy_pkg.f90wrap_scalapack_module__set__nprow(nprow)

def get_my_blacs_id():
    """
    Element my_blacs_id ftype=integer(i4b) pytype=int
    
    
    Defined at ScaLapack_module.fpp line 19
    
    """
    return _arespy_pkg.f90wrap_scalapack_module__get__my_blacs_id()

def set_my_blacs_id(my_blacs_id):
    _arespy_pkg.f90wrap_scalapack_module__set__my_blacs_id(my_blacs_id)

def get_np():
    """
    Element np ftype=integer(i4b) pytype=int
    
    
    Defined at ScaLapack_module.fpp line 20
    
    """
    return _arespy_pkg.f90wrap_scalapack_module__get__np()

def set_np(np):
    _arespy_pkg.f90wrap_scalapack_module__set__np(np)

def get_nq():
    """
    Element nq ftype=integer(i4b) pytype=int
    
    
    Defined at ScaLapack_module.fpp line 20
    
    """
    return _arespy_pkg.f90wrap_scalapack_module__get__nq()

def set_nq(nq):
    _arespy_pkg.f90wrap_scalapack_module__set__nq(nq)

def get_iam():
    """
    Element iam ftype=integer(i4b) pytype=int
    
    
    Defined at ScaLapack_module.fpp line 21
    
    """
    return _arespy_pkg.f90wrap_scalapack_module__get__iam()

def set_iam(iam):
    _arespy_pkg.f90wrap_scalapack_module__set__iam(iam)

def get_nprocs():
    """
    Element nprocs ftype=integer(i4b) pytype=int
    
    
    Defined at ScaLapack_module.fpp line 21
    
    """
    return _arespy_pkg.f90wrap_scalapack_module__get__nprocs()

def set_nprocs(nprocs):
    _arespy_pkg.f90wrap_scalapack_module__set__nprocs(nprocs)

def get_init_called():
    """
    Element init_called ftype=logical pytype=bool
    
    
    Defined at ScaLapack_module.fpp line 22
    
    """
    return _arespy_pkg.f90wrap_scalapack_module__get__init_called()

def set_init_called(init_called):
    _arespy_pkg.f90wrap_scalapack_module__set__init_called(init_called)

def get_info():
    """
    Element info ftype=integer(i4b) pytype=int
    
    
    Defined at ScaLapack_module.fpp line 24
    
    """
    return _arespy_pkg.f90wrap_scalapack_module__get__info()

def set_info(info):
    _arespy_pkg.f90wrap_scalapack_module__set__info(info)

def get_l_useless():
    """
    Element l_useless ftype=logical pytype=bool
    
    
    Defined at ScaLapack_module.fpp line 25
    
    """
    return _arespy_pkg.f90wrap_scalapack_module__get__l_useless()

def set_l_useless(l_useless):
    _arespy_pkg.f90wrap_scalapack_module__set__l_useless(l_useless)

def get_array_twod_map():
    """
    Element twod_map ftype=integer(i4b) pytype=int
    
    
    Defined at ScaLapack_module.fpp line 26
    
    """
    global twod_map
    array_ndim, array_type, array_shape, array_handle = \
        _arespy_pkg.f90wrap_scalapack_module__array__twod_map(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        twod_map = _arrays[array_handle]
    else:
        twod_map = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_scalapack_module__array__twod_map)
        _arrays[array_handle] = twod_map
    return twod_map

def set_array_twod_map(twod_map):
    twod_map[...] = twod_map


_array_initialisers = [get_array_twod_map]
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "scalapack_module".')

for func in _dt_array_initialisers:
    func()

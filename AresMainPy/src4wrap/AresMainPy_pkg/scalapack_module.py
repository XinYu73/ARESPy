"""
Module scalapack_module


Defined at Scala.fpp lines 5-1149

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def init_scala():
    """
    init_scala()
    
    
    Defined at Scala.fpp lines 29-47
    
    
    """
    _AresMainPy_pkg.f90wrap_init_scala()

def sl_orthnorm(amat, m, n, mb, nb):
    """
    sl_orthnorm(amat, m, n, mb, nb)
    
    
    Defined at Scala.fpp lines 59-78
    
    Parameters
    ----------
    amat : complex array
    m : int
    n : int
    mb : int
    nb : int
    
    """
    _AresMainPy_pkg.f90wrap_sl_orthnorm(amat=amat, m=m, n=n, mb=mb, nb=nb)

def sl_orthnorm_real(amat, m, n, mb, nb):
    """
    sl_orthnorm_real(amat, m, n, mb, nb)
    
    
    Defined at Scala.fpp lines 81-117
    
    Parameters
    ----------
    amat : float array
    m : int
    n : int
    mb : int
    nb : int
    
    """
    _AresMainPy_pkg.f90wrap_sl_orthnorm_real(amat=amat, m=m, n=n, mb=mb, nb=nb)

def sl_matmat(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb, \
    cmbin=None, cnbin=None):
    """
    sl_matmat(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb[, \
        cmbin, cnbin])
    
    
    Defined at Scala.fpp lines 120-176
    
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
    _AresMainPy_pkg.f90wrap_sl_matmat(opa=opa, opb=opb, amat=amat, bmat=bmat, \
        cmat=cmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, \
        cmbin=cmbin, cnbin=cnbin)

def sl_matmat_gridcn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, \
    bnb, cmb, cnb):
    """
    sl_matmat_gridcn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb, \
        cmb, cnb)
    
    
    Defined at Scala.fpp lines 179-288
    
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
    _AresMainPy_pkg.f90wrap_sl_matmat_gridcn(opa=opa, opb=opb, amat=amat, bmat=bmat, \
        cmat=cmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, \
        cmb=cmb, cnb=cnb)

def sl_matmat_gridnn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, \
    bnb, cmb, cnb):
    """
    sl_matmat_gridnn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb, \
        cmb, cnb)
    
    
    Defined at Scala.fpp lines 291-350
    
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
    _AresMainPy_pkg.f90wrap_sl_matmat_gridnn(opa=opa, opb=opb, amat=amat, bmat=bmat, \
        cmat=cmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, \
        cmb=cmb, cnb=cnb)

def sl_matmat_sub(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, \
    bnb, cmbin=None, cnbin=None):
    """
    sl_matmat_sub(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb[, \
        cmbin, cnbin])
    
    
    Defined at Scala.fpp lines 353-438
    
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
    _AresMainPy_pkg.f90wrap_sl_matmat_sub(opa=opa, opb=opb, amat=amat, bmat=bmat, \
        cmat=cmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, \
        cmbin=cmbin, cnbin=cnbin)

def sl_matmat_sub_real(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, \
    bmb, bnb, cmbin=None, cnbin=None):
    """
    sl_matmat_sub_real(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, \
        bnb[, cmbin, cnbin])
    
    
    Defined at Scala.fpp lines 441-526
    
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
    _AresMainPy_pkg.f90wrap_sl_matmat_sub_real(opa=opa, opb=opb, amat=amat, \
        bmat=bmat, cmat=cmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, \
        bnb=bnb, cmbin=cmbin, cnbin=cnbin)

def sl_matmat_real(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, \
    bnb, cmbin=None, cnbin=None):
    """
    sl_matmat_real(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb[, \
        cmbin, cnbin])
    
    
    Defined at Scala.fpp lines 529-626
    
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
    _AresMainPy_pkg.f90wrap_sl_matmat_real(opa=opa, opb=opb, amat=amat, bmat=bmat, \
        cmat=cmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, \
        cmbin=cmbin, cnbin=cnbin)

def sl_generalizeeigen(dime, amat, bmat, am, an, bm, bn, amb, anb, bmb, bnb, \
    evec, cm, cn, eval, cmbin=None, cnbin=None):
    """
    sl_generalizeeigen(dime, amat, bmat, am, an, bm, bn, amb, anb, bmb, bnb, evec, \
        cm, cn, eval[, cmbin, cnbin])
    
    
    Defined at Scala.fpp lines 629-712
    
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
    _AresMainPy_pkg.f90wrap_sl_generalizeeigen(dime=dime, amat=amat, bmat=bmat, \
        am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, evec=evec, \
        cm=cm, cn=cn, eval=eval, cmbin=cmbin, cnbin=cnbin)

def sl_generalizeeigen_real2(dime, amat, bmat, am, an, bm, bn, amb, anb, bmb, \
    bnb, evec, cm, cn, eval, cmbin=None, cnbin=None):
    """
    sl_generalizeeigen_real2(dime, amat, bmat, am, an, bm, bn, amb, anb, bmb, bnb, \
        evec, cm, cn, eval[, cmbin, cnbin])
    
    
    Defined at Scala.fpp lines 715-801
    
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
    _AresMainPy_pkg.f90wrap_sl_generalizeeigen_real2(dime=dime, amat=amat, \
        bmat=bmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, \
        evec=evec, cm=cm, cn=cn, eval=eval, cmbin=cmbin, cnbin=cnbin)

def sl_generalizeeigen_real(dime, amat, bmat, am, an, bm, bn, amb, anb, bmb, \
    bnb, evec, cm, cn, eval, cmbin=None, cnbin=None):
    """
    sl_generalizeeigen_real(dime, amat, bmat, am, an, bm, bn, amb, anb, bmb, bnb, \
        evec, cm, cn, eval[, cmbin, cnbin])
    
    
    Defined at Scala.fpp lines 804-953
    
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
    _AresMainPy_pkg.f90wrap_sl_generalizeeigen_real(dime=dime, amat=amat, bmat=bmat, \
        am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, evec=evec, \
        cm=cm, cn=cn, eval=eval, cmbin=cmbin, cnbin=cnbin)

def twod_map_set(nstates, nrow, ncol, twod_map):
    """
    twod_map_set(nstates, nrow, ncol, twod_map)
    
    
    Defined at Scala.fpp lines 957-975
    
    Parameters
    ----------
    nstates : int
    nrow : int
    ncol : int
    twod_map : int array
    
    """
    _AresMainPy_pkg.f90wrap_twod_map_set(nstates=nstates, nrow=nrow, ncol=ncol, \
        twod_map=twod_map)

def sl_matmat_realtn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, \
    bnb, cmb, cnb):
    """
    sl_matmat_realtn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb, \
        cmb, cnb)
    
    
    Defined at Scala.fpp lines 977-1086
    
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
    _AresMainPy_pkg.f90wrap_sl_matmat_realtn(opa=opa, opb=opb, amat=amat, bmat=bmat, \
        cmat=cmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, \
        cmb=cmb, cnb=cnb)

def sl_matmat_realnn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, \
    bnb, cmb, cnb):
    """
    sl_matmat_realnn(opa, opb, amat, bmat, cmat, am, an, bm, bn, amb, anb, bmb, bnb, \
        cmb, cnb)
    
    
    Defined at Scala.fpp lines 1089-1148
    
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
    _AresMainPy_pkg.f90wrap_sl_matmat_realnn(opa=opa, opb=opb, amat=amat, bmat=bmat, \
        cmat=cmat, am=am, an=an, bm=bm, bn=bn, amb=amb, anb=anb, bmb=bmb, bnb=bnb, \
        cmb=cmb, cnb=cnb)

def get_blacs_contxt():
    """
    Element blacs_contxt ftype=integer(i4b) pytype=int
    
    
    Defined at Scala.fpp line 15
    
    """
    return _AresMainPy_pkg.f90wrap_scalapack_module__get__blacs_contxt()

def set_blacs_contxt(blacs_contxt):
    _AresMainPy_pkg.f90wrap_scalapack_module__set__blacs_contxt(blacs_contxt)

def get_dlen():
    """
    Element dlen ftype=integer(i4b) pytype=int
    
    
    Defined at Scala.fpp line 16
    
    """
    return _AresMainPy_pkg.f90wrap_scalapack_module__get__dlen()

DLEN = get_dlen()

def get_myrow():
    """
    Element myrow ftype=integer(i4b) pytype=int
    
    
    Defined at Scala.fpp line 17
    
    """
    return _AresMainPy_pkg.f90wrap_scalapack_module__get__myrow()

def set_myrow(myrow):
    _AresMainPy_pkg.f90wrap_scalapack_module__set__myrow(myrow)

def get_mycol():
    """
    Element mycol ftype=integer(i4b) pytype=int
    
    
    Defined at Scala.fpp line 17
    
    """
    return _AresMainPy_pkg.f90wrap_scalapack_module__get__mycol()

def set_mycol(mycol):
    _AresMainPy_pkg.f90wrap_scalapack_module__set__mycol(mycol)

def get_npcol():
    """
    Element npcol ftype=integer(i4b) pytype=int
    
    
    Defined at Scala.fpp line 18
    
    """
    return _AresMainPy_pkg.f90wrap_scalapack_module__get__npcol()

def set_npcol(npcol):
    _AresMainPy_pkg.f90wrap_scalapack_module__set__npcol(npcol)

def get_nprow():
    """
    Element nprow ftype=integer(i4b) pytype=int
    
    
    Defined at Scala.fpp line 18
    
    """
    return _AresMainPy_pkg.f90wrap_scalapack_module__get__nprow()

def set_nprow(nprow):
    _AresMainPy_pkg.f90wrap_scalapack_module__set__nprow(nprow)

def get_my_blacs_id():
    """
    Element my_blacs_id ftype=integer(i4b) pytype=int
    
    
    Defined at Scala.fpp line 19
    
    """
    return _AresMainPy_pkg.f90wrap_scalapack_module__get__my_blacs_id()

def set_my_blacs_id(my_blacs_id):
    _AresMainPy_pkg.f90wrap_scalapack_module__set__my_blacs_id(my_blacs_id)

def get_np():
    """
    Element np ftype=integer(i4b) pytype=int
    
    
    Defined at Scala.fpp line 20
    
    """
    return _AresMainPy_pkg.f90wrap_scalapack_module__get__np()

def set_np(np):
    _AresMainPy_pkg.f90wrap_scalapack_module__set__np(np)

def get_nq():
    """
    Element nq ftype=integer(i4b) pytype=int
    
    
    Defined at Scala.fpp line 20
    
    """
    return _AresMainPy_pkg.f90wrap_scalapack_module__get__nq()

def set_nq(nq):
    _AresMainPy_pkg.f90wrap_scalapack_module__set__nq(nq)

def get_iam():
    """
    Element iam ftype=integer(i4b) pytype=int
    
    
    Defined at Scala.fpp line 21
    
    """
    return _AresMainPy_pkg.f90wrap_scalapack_module__get__iam()

def set_iam(iam):
    _AresMainPy_pkg.f90wrap_scalapack_module__set__iam(iam)

def get_nprocs():
    """
    Element nprocs ftype=integer(i4b) pytype=int
    
    
    Defined at Scala.fpp line 21
    
    """
    return _AresMainPy_pkg.f90wrap_scalapack_module__get__nprocs()

def set_nprocs(nprocs):
    _AresMainPy_pkg.f90wrap_scalapack_module__set__nprocs(nprocs)

def get_init_called():
    """
    Element init_called ftype=logical pytype=bool
    
    
    Defined at Scala.fpp line 22
    
    """
    return _AresMainPy_pkg.f90wrap_scalapack_module__get__init_called()

def set_init_called(init_called):
    _AresMainPy_pkg.f90wrap_scalapack_module__set__init_called(init_called)

def get_info():
    """
    Element info ftype=integer(i4b) pytype=int
    
    
    Defined at Scala.fpp line 24
    
    """
    return _AresMainPy_pkg.f90wrap_scalapack_module__get__info()

def set_info(info):
    _AresMainPy_pkg.f90wrap_scalapack_module__set__info(info)

def get_l_useless():
    """
    Element l_useless ftype=logical pytype=bool
    
    
    Defined at Scala.fpp line 25
    
    """
    return _AresMainPy_pkg.f90wrap_scalapack_module__get__l_useless()

def set_l_useless(l_useless):
    _AresMainPy_pkg.f90wrap_scalapack_module__set__l_useless(l_useless)

def get_array_twod_map():
    """
    Element twod_map ftype=integer(i4b) pytype=int
    
    
    Defined at Scala.fpp line 26
    
    """
    global twod_map
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_scalapack_module__array__twod_map(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        twod_map = _arrays[array_handle]
    else:
        twod_map = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_scalapack_module__array__twod_map)
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

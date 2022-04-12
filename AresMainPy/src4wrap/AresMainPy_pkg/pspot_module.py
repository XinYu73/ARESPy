"""
Module pspot_module


Defined at Psp_module.fpp lines 10-48

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("AresMainPy_pkg.pspot")
class pspot(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=pspot)
    
    
    Defined at Psp_module.fpp lines 14-42
    
    """
    def __init__(self, handle=None):
        """
        self = Pspot()
        
        
        Defined at Psp_module.fpp lines 14-42
        
        
        Returns
        -------
        this : Pspot
        	Object to be constructed
        
        
        Automatically generated constructor for pspot
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _AresMainPy_pkg.f90wrap_pspot_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Pspot
        
        
        Defined at Psp_module.fpp lines 14-42
        
        Parameters
        ----------
        this : Pspot
        	Object to be destructed
        
        
        Automatically generated destructor for pspot
        """
        if self._alloc:
            _AresMainPy_pkg.f90wrap_pspot_finalise(this=self._handle)
    
    @property
    def zion(self):
        """
        Element zion ftype=real(dp) pytype=float
        
        
        Defined at Psp_module.fpp line 15
        
        """
        return _AresMainPy_pkg.f90wrap_pspot__get__zion(self._handle)
    
    @zion.setter
    def zion(self, zion):
        _AresMainPy_pkg.f90wrap_pspot__set__zion(self._handle, zion)
    
    @property
    def numps(self):
        """
        Element numps ftype=integer(i4b) pytype=int
        
        
        Defined at Psp_module.fpp line 16
        
        """
        return _AresMainPy_pkg.f90wrap_pspot__get__numps(self._handle)
    
    @numps.setter
    def numps(self, numps):
        _AresMainPy_pkg.f90wrap_pspot__set__numps(self._handle, numps)
    
    @property
    def qnumps(self):
        """
        Element qnumps ftype=integer(i4b) pytype=int
        
        
        Defined at Psp_module.fpp line 17
        
        """
        return _AresMainPy_pkg.f90wrap_pspot__get__qnumps(self._handle)
    
    @qnumps.setter
    def qnumps(self, qnumps):
        _AresMainPy_pkg.f90wrap_pspot__set__qnumps(self._handle, qnumps)
    
    @property
    def numps_den(self):
        """
        Element numps_den ftype=integer(i4b) pytype=int
        
        
        Defined at Psp_module.fpp line 19
        
        """
        return _AresMainPy_pkg.f90wrap_pspot__get__numps_den(self._handle)
    
    @numps_den.setter
    def numps_den(self, numps_den):
        _AresMainPy_pkg.f90wrap_pspot__set__numps_den(self._handle, numps_den)
    
    @property
    def qmax(self):
        """
        Element qmax ftype=real(dp) pytype=float
        
        
        Defined at Psp_module.fpp line 21
        
        """
        return _AresMainPy_pkg.f90wrap_pspot__get__qmax(self._handle)
    
    @qmax.setter
    def qmax(self, qmax):
        _AresMainPy_pkg.f90wrap_pspot__set__qmax(self._handle, qmax)
    
    @property
    def qspacing(self):
        """
        Element qspacing ftype=real(dp) pytype=float
        
        
        Defined at Psp_module.fpp line 22
        
        """
        return _AresMainPy_pkg.f90wrap_pspot__get__qspacing(self._handle)
    
    @qspacing.setter
    def qspacing(self, qspacing):
        _AresMainPy_pkg.f90wrap_pspot__set__qspacing(self._handle, qspacing)
    
    @property
    def qmesh(self):
        """
        Element qmesh ftype=real(dp) pytype=float
        
        
        Defined at Psp_module.fpp line 23
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_pspot__array__qmesh(self._handle)
        if array_handle in self._arrays:
            qmesh = self._arrays[array_handle]
        else:
            qmesh = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_pspot__array__qmesh)
            self._arrays[array_handle] = qmesh
        return qmesh
    
    @qmesh.setter
    def qmesh(self, qmesh):
        self.qmesh[...] = qmesh
    
    @property
    def vlocq(self):
        """
        Element vlocq ftype=real(dp) pytype=float
        
        
        Defined at Psp_module.fpp line 24
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_pspot__array__vlocq(self._handle)
        if array_handle in self._arrays:
            vlocq = self._arrays[array_handle]
        else:
            vlocq = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_pspot__array__vlocq)
            self._arrays[array_handle] = vlocq
        return vlocq
    
    @vlocq.setter
    def vlocq(self, vlocq):
        self.vlocq[...] = vlocq
    
    @property
    def vlocqs(self):
        """
        Element vlocqs ftype=real(dp) pytype=float
        
        
        Defined at Psp_module.fpp line 25
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_pspot__array__vlocqs(self._handle)
        if array_handle in self._arrays:
            vlocqs = self._arrays[array_handle]
        else:
            vlocqs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_pspot__array__vlocqs)
            self._arrays[array_handle] = vlocqs
        return vlocqs
    
    @vlocqs.setter
    def vlocqs(self, vlocqs):
        self.vlocqs[...] = vlocqs
    
    @property
    def ddvl_dq2(self):
        """
        Element ddvl_dq2 ftype=real(dp) pytype=float
        
        
        Defined at Psp_module.fpp line 26
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_pspot__array__ddvl_dq2(self._handle)
        if array_handle in self._arrays:
            ddvl_dq2 = self._arrays[array_handle]
        else:
            ddvl_dq2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_pspot__array__ddvl_dq2)
            self._arrays[array_handle] = ddvl_dq2
        return ddvl_dq2
    
    @ddvl_dq2.setter
    def ddvl_dq2(self, ddvl_dq2):
        self.ddvl_dq2[...] = ddvl_dq2
    
    @property
    def nproj(self):
        """
        Element nproj ftype=integer(i4b) pytype=int
        
        
        Defined at Psp_module.fpp line 28
        
        """
        return _AresMainPy_pkg.f90wrap_pspot__get__nproj(self._handle)
    
    @nproj.setter
    def nproj(self, nproj):
        _AresMainPy_pkg.f90wrap_pspot__set__nproj(self._handle, nproj)
    
    @property
    def proj_l(self):
        """
        Element proj_l ftype=integer(i4b) pytype=int
        
        
        Defined at Psp_module.fpp line 29
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_pspot__array__proj_l(self._handle)
        if array_handle in self._arrays:
            proj_l = self._arrays[array_handle]
        else:
            proj_l = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_pspot__array__proj_l)
            self._arrays[array_handle] = proj_l
        return proj_l
    
    @proj_l.setter
    def proj_l(self, proj_l):
        self.proj_l[...] = proj_l
    
    @property
    def proj_m(self):
        """
        Element proj_m ftype=integer(i4b) pytype=int
        
        
        Defined at Psp_module.fpp line 30
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_pspot__array__proj_m(self._handle)
        if array_handle in self._arrays:
            proj_m = self._arrays[array_handle]
        else:
            proj_m = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_pspot__array__proj_m)
            self._arrays[array_handle] = proj_m
        return proj_m
    
    @proj_m.setter
    def proj_m(self, proj_m):
        self.proj_m[...] = proj_m
    
    @property
    def rcut(self):
        """
        Element rcut ftype=real(dp) pytype=float
        
        
        Defined at Psp_module.fpp line 31
        
        """
        return _AresMainPy_pkg.f90wrap_pspot__get__rcut(self._handle)
    
    @rcut.setter
    def rcut(self, rcut):
        _AresMainPy_pkg.f90wrap_pspot__set__rcut(self._handle, rcut)
    
    @property
    def rmax(self):
        """
        Element rmax ftype=real(dp) pytype=float
        
        
        Defined at Psp_module.fpp line 32
        
        """
        return _AresMainPy_pkg.f90wrap_pspot__get__rmax(self._handle)
    
    @rmax.setter
    def rmax(self, rmax):
        _AresMainPy_pkg.f90wrap_pspot__set__rmax(self._handle, rmax)
    
    @property
    def rspacing(self):
        """
        Element rspacing ftype=real(dp) pytype=float
        
        
        Defined at Psp_module.fpp line 33
        
        """
        return _AresMainPy_pkg.f90wrap_pspot__get__rspacing(self._handle)
    
    @rspacing.setter
    def rspacing(self, rspacing):
        _AresMainPy_pkg.f90wrap_pspot__set__rspacing(self._handle, rspacing)
    
    @property
    def d0(self):
        """
        Element d0 ftype=real(dp) pytype=float
        
        
        Defined at Psp_module.fpp line 34
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_pspot__array__d0(self._handle)
        if array_handle in self._arrays:
            d0 = self._arrays[array_handle]
        else:
            d0 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_pspot__array__d0)
            self._arrays[array_handle] = d0
        return d0
    
    @d0.setter
    def d0(self, d0):
        self.d0[...] = d0
    
    @property
    def beta_r(self):
        """
        Element beta_r ftype=real(dp) pytype=float
        
        
        Defined at Psp_module.fpp line 35
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_pspot__array__beta_r(self._handle)
        if array_handle in self._arrays:
            beta_r = self._arrays[array_handle]
        else:
            beta_r = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_pspot__array__beta_r)
            self._arrays[array_handle] = beta_r
        return beta_r
    
    @beta_r.setter
    def beta_r(self, beta_r):
        self.beta_r[...] = beta_r
    
    @property
    def dbeta_dr(self):
        """
        Element dbeta_dr ftype=real(dp) pytype=float
        
        
        Defined at Psp_module.fpp line 36
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_pspot__array__dbeta_dr(self._handle)
        if array_handle in self._arrays:
            dbeta_dr = self._arrays[array_handle]
        else:
            dbeta_dr = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_pspot__array__dbeta_dr)
            self._arrays[array_handle] = dbeta_dr
        return dbeta_dr
    
    @dbeta_dr.setter
    def dbeta_dr(self, dbeta_dr):
        self.dbeta_dr[...] = dbeta_dr
    
    @property
    def denr(self):
        """
        Element denr ftype=real(dp) pytype=float
        
        
        Defined at Psp_module.fpp line 38
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_pspot__array__denr(self._handle)
        if array_handle in self._arrays:
            denr = self._arrays[array_handle]
        else:
            denr = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_pspot__array__denr)
            self._arrays[array_handle] = denr
        return denr
    
    @denr.setter
    def denr(self, denr):
        self.denr[...] = denr
    
    @property
    def ddden_dr2(self):
        """
        Element ddden_dr2 ftype=real(dp) pytype=float
        
        
        Defined at Psp_module.fpp line 39
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_pspot__array__ddden_dr2(self._handle)
        if array_handle in self._arrays:
            ddden_dr2 = self._arrays[array_handle]
        else:
            ddden_dr2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_pspot__array__ddden_dr2)
            self._arrays[array_handle] = ddden_dr2
        return ddden_dr2
    
    @ddden_dr2.setter
    def ddden_dr2(self, ddden_dr2):
        self.ddden_dr2[...] = ddden_dr2
    
    @property
    def r_real(self):
        """
        Element r_real ftype=real(dp) pytype=float
        
        
        Defined at Psp_module.fpp line 41
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_pspot__array__r_real(self._handle)
        if array_handle in self._arrays:
            r_real = self._arrays[array_handle]
        else:
            r_real = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_pspot__array__r_real)
            self._arrays[array_handle] = r_real
        return r_real
    
    @r_real.setter
    def r_real(self, r_real):
        self.r_real[...] = r_real
    
    @property
    def v_loc(self):
        """
        Element v_loc ftype=real(dp) pytype=float
        
        
        Defined at Psp_module.fpp line 42
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_pspot__array__v_loc(self._handle)
        if array_handle in self._arrays:
            v_loc = self._arrays[array_handle]
        else:
            v_loc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_pspot__array__v_loc)
            self._arrays[array_handle] = v_loc
        return v_loc
    
    @v_loc.setter
    def v_loc(self, v_loc):
        self.v_loc[...] = v_loc
    
    def __str__(self):
        ret = ['<pspot>{\n']
        ret.append('    zion : ')
        ret.append(repr(self.zion))
        ret.append(',\n    numps : ')
        ret.append(repr(self.numps))
        ret.append(',\n    qnumps : ')
        ret.append(repr(self.qnumps))
        ret.append(',\n    numps_den : ')
        ret.append(repr(self.numps_den))
        ret.append(',\n    qmax : ')
        ret.append(repr(self.qmax))
        ret.append(',\n    qspacing : ')
        ret.append(repr(self.qspacing))
        ret.append(',\n    qmesh : ')
        ret.append(repr(self.qmesh))
        ret.append(',\n    vlocq : ')
        ret.append(repr(self.vlocq))
        ret.append(',\n    vlocqs : ')
        ret.append(repr(self.vlocqs))
        ret.append(',\n    ddvl_dq2 : ')
        ret.append(repr(self.ddvl_dq2))
        ret.append(',\n    nproj : ')
        ret.append(repr(self.nproj))
        ret.append(',\n    proj_l : ')
        ret.append(repr(self.proj_l))
        ret.append(',\n    proj_m : ')
        ret.append(repr(self.proj_m))
        ret.append(',\n    rcut : ')
        ret.append(repr(self.rcut))
        ret.append(',\n    rmax : ')
        ret.append(repr(self.rmax))
        ret.append(',\n    rspacing : ')
        ret.append(repr(self.rspacing))
        ret.append(',\n    d0 : ')
        ret.append(repr(self.d0))
        ret.append(',\n    beta_r : ')
        ret.append(repr(self.beta_r))
        ret.append(',\n    dbeta_dr : ')
        ret.append(repr(self.dbeta_dr))
        ret.append(',\n    denr : ')
        ret.append(repr(self.denr))
        ret.append(',\n    ddden_dr2 : ')
        ret.append(repr(self.ddden_dr2))
        ret.append(',\n    r_real : ')
        ret.append(repr(self.r_real))
        ret.append(',\n    v_loc : ')
        ret.append(repr(self.v_loc))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def get_max_nproj():
    """
    Element max_nproj ftype=integer(i4b) pytype=int
    
    
    Defined at Psp_module.fpp line 47
    
    """
    return _AresMainPy_pkg.f90wrap_pspot_module__get__max_nproj()

def set_max_nproj(max_nproj):
    _AresMainPy_pkg.f90wrap_pspot_module__set__max_nproj(max_nproj)

def get_max_rcut():
    """
    Element max_rcut ftype=real(dp) pytype=float
    
    
    Defined at Psp_module.fpp line 48
    
    """
    return _AresMainPy_pkg.f90wrap_pspot_module__get__max_rcut()

def set_max_rcut(max_rcut):
    _AresMainPy_pkg.f90wrap_pspot_module__set__max_rcut(max_rcut)

def get_array_tknots():
    """
    Element tknots ftype=real(dp) pytype=float
    
    
    Defined at Psp_module.fpp line 49
    
    """
    global tknots
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_pspot_module__array__tknots(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        tknots = _arrays[array_handle]
    else:
        tknots = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_pspot_module__array__tknots)
        _arrays[array_handle] = tknots
    return tknots

def set_array_tknots(tknots):
    tknots[...] = tknots


_array_initialisers = [get_array_tknots]
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "pspot_module".')

for func in _dt_array_initialisers:
    func()

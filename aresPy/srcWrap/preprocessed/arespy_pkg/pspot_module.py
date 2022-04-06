"""
Module pspot_module


Defined at Struct_module.fpp lines 105-1051

"""
from __future__ import print_function, absolute_import, division
import _arespy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("arespy_pkg.pspot")
class pspot(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=pspot)
    
    
    Defined at Struct_module.fpp lines 109-149
    
    """
    def __init__(self, handle=None):
        """
        self = Pspot()
        
        
        Defined at Struct_module.fpp lines 109-149
        
        
        Returns
        -------
        this : Pspot
        	Object to be constructed
        
        
        Automatically generated constructor for pspot
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _arespy_pkg.f90wrap_pspot_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Pspot
        
        
        Defined at Struct_module.fpp lines 109-149
        
        Parameters
        ----------
        this : Pspot
        	Object to be destructed
        
        
        Automatically generated destructor for pspot
        """
        if self._alloc:
            _arespy_pkg.f90wrap_pspot_finalise(this=self._handle)
    
    @property
    def zion(self):
        """
        Element zion ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 110
        
        """
        return _arespy_pkg.f90wrap_pspot__get__zion(self._handle)
    
    @zion.setter
    def zion(self, zion):
        _arespy_pkg.f90wrap_pspot__set__zion(self._handle, zion)
    
    @property
    def numps(self):
        """
        Element numps ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.fpp line 111
        
        """
        return _arespy_pkg.f90wrap_pspot__get__numps(self._handle)
    
    @numps.setter
    def numps(self, numps):
        _arespy_pkg.f90wrap_pspot__set__numps(self._handle, numps)
    
    @property
    def qnumps(self):
        """
        Element qnumps ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.fpp line 112
        
        """
        return _arespy_pkg.f90wrap_pspot__get__qnumps(self._handle)
    
    @qnumps.setter
    def qnumps(self, qnumps):
        _arespy_pkg.f90wrap_pspot__set__qnumps(self._handle, qnumps)
    
    @property
    def elename(self):
        """
        Element elename ftype=character(len=3) pytype=str
        
        
        Defined at Struct_module.fpp line 113
        
        """
        return _arespy_pkg.f90wrap_pspot__get__elename(self._handle)
    
    @elename.setter
    def elename(self, elename):
        _arespy_pkg.f90wrap_pspot__set__elename(self._handle, elename)
    
    @property
    def qmax(self):
        """
        Element qmax ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 115
        
        """
        return _arespy_pkg.f90wrap_pspot__get__qmax(self._handle)
    
    @qmax.setter
    def qmax(self, qmax):
        _arespy_pkg.f90wrap_pspot__set__qmax(self._handle, qmax)
    
    @property
    def qspacing(self):
        """
        Element qspacing ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 116
        
        """
        return _arespy_pkg.f90wrap_pspot__get__qspacing(self._handle)
    
    @qspacing.setter
    def qspacing(self, qspacing):
        _arespy_pkg.f90wrap_pspot__set__qspacing(self._handle, qspacing)
    
    @property
    def qmesh(self):
        """
        Element qmesh ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 117
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_pspot__array__qmesh(self._handle)
        if array_handle in self._arrays:
            qmesh = self._arrays[array_handle]
        else:
            qmesh = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_pspot__array__qmesh)
            self._arrays[array_handle] = qmesh
        return qmesh
    
    @qmesh.setter
    def qmesh(self, qmesh):
        self.qmesh[...] = qmesh
    
    @property
    def vlocq(self):
        """
        Element vlocq ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 118
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_pspot__array__vlocq(self._handle)
        if array_handle in self._arrays:
            vlocq = self._arrays[array_handle]
        else:
            vlocq = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_pspot__array__vlocq)
            self._arrays[array_handle] = vlocq
        return vlocq
    
    @vlocq.setter
    def vlocq(self, vlocq):
        self.vlocq[...] = vlocq
    
    @property
    def vlocqs(self):
        """
        Element vlocqs ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 119
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_pspot__array__vlocqs(self._handle)
        if array_handle in self._arrays:
            vlocqs = self._arrays[array_handle]
        else:
            vlocqs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_pspot__array__vlocqs)
            self._arrays[array_handle] = vlocqs
        return vlocqs
    
    @vlocqs.setter
    def vlocqs(self, vlocqs):
        self.vlocqs[...] = vlocqs
    
    @property
    def ddvl_dq2(self):
        """
        Element ddvl_dq2 ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 120
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_pspot__array__ddvl_dq2(self._handle)
        if array_handle in self._arrays:
            ddvl_dq2 = self._arrays[array_handle]
        else:
            ddvl_dq2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_pspot__array__ddvl_dq2)
            self._arrays[array_handle] = ddvl_dq2
        return ddvl_dq2
    
    @ddvl_dq2.setter
    def ddvl_dq2(self, ddvl_dq2):
        self.ddvl_dq2[...] = ddvl_dq2
    
    @property
    def r(self):
        """
        Element r ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 122
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_pspot__array__r(self._handle)
        if array_handle in self._arrays:
            r = self._arrays[array_handle]
        else:
            r = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_pspot__array__r)
            self._arrays[array_handle] = r
        return r
    
    @r.setter
    def r(self, r):
        self.r[...] = r
    
    @property
    def rab(self):
        """
        Element rab ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 122
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_pspot__array__rab(self._handle)
        if array_handle in self._arrays:
            rab = self._arrays[array_handle]
        else:
            rab = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_pspot__array__rab)
            self._arrays[array_handle] = rab
        return rab
    
    @rab.setter
    def rab(self, rab):
        self.rab[...] = rab
    
    @property
    def vlocr(self):
        """
        Element vlocr ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 123
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_pspot__array__vlocr(self._handle)
        if array_handle in self._arrays:
            vlocr = self._arrays[array_handle]
        else:
            vlocr = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_pspot__array__vlocr)
            self._arrays[array_handle] = vlocr
        return vlocr
    
    @vlocr.setter
    def vlocr(self, vlocr):
        self.vlocr[...] = vlocr
    
    @property
    def nproj(self):
        """
        Element nproj ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.fpp line 125
        
        """
        return _arespy_pkg.f90wrap_pspot__get__nproj(self._handle)
    
    @nproj.setter
    def nproj(self, nproj):
        _arespy_pkg.f90wrap_pspot__set__nproj(self._handle, nproj)
    
    @property
    def proj_l(self):
        """
        Element proj_l ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.fpp line 126
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_pspot__array__proj_l(self._handle)
        if array_handle in self._arrays:
            proj_l = self._arrays[array_handle]
        else:
            proj_l = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_pspot__array__proj_l)
            self._arrays[array_handle] = proj_l
        return proj_l
    
    @proj_l.setter
    def proj_l(self, proj_l):
        self.proj_l[...] = proj_l
    
    @property
    def proj_m(self):
        """
        Element proj_m ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.fpp line 127
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_pspot__array__proj_m(self._handle)
        if array_handle in self._arrays:
            proj_m = self._arrays[array_handle]
        else:
            proj_m = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_pspot__array__proj_m)
            self._arrays[array_handle] = proj_m
        return proj_m
    
    @proj_m.setter
    def proj_m(self, proj_m):
        self.proj_m[...] = proj_m
    
    @property
    def indx(self):
        """
        Element indx ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.fpp line 127
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_pspot__array__indx(self._handle)
        if array_handle in self._arrays:
            indx = self._arrays[array_handle]
        else:
            indx = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_pspot__array__indx)
            self._arrays[array_handle] = indx
        return indx
    
    @indx.setter
    def indx(self, indx):
        self.indx[...] = indx
    
    @property
    def rcut(self):
        """
        Element rcut ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 128
        
        """
        return _arespy_pkg.f90wrap_pspot__get__rcut(self._handle)
    
    @rcut.setter
    def rcut(self, rcut):
        _arespy_pkg.f90wrap_pspot__set__rcut(self._handle, rcut)
    
    @property
    def rmax(self):
        """
        Element rmax ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 129
        
        """
        return _arespy_pkg.f90wrap_pspot__get__rmax(self._handle)
    
    @rmax.setter
    def rmax(self, rmax):
        _arespy_pkg.f90wrap_pspot__set__rmax(self._handle, rmax)
    
    @property
    def rspacing(self):
        """
        Element rspacing ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 130
        
        """
        return _arespy_pkg.f90wrap_pspot__get__rspacing(self._handle)
    
    @rspacing.setter
    def rspacing(self, rspacing):
        _arespy_pkg.f90wrap_pspot__set__rspacing(self._handle, rspacing)
    
    @property
    def dij(self):
        """
        Element dij ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 131
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_pspot__array__dij(self._handle)
        if array_handle in self._arrays:
            dij = self._arrays[array_handle]
        else:
            dij = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_pspot__array__dij)
            self._arrays[array_handle] = dij
        return dij
    
    @dij.setter
    def dij(self, dij):
        self.dij[...] = dij
    
    @property
    def beta_r(self):
        """
        Element beta_r ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 132
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_pspot__array__beta_r(self._handle)
        if array_handle in self._arrays:
            beta_r = self._arrays[array_handle]
        else:
            beta_r = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_pspot__array__beta_r)
            self._arrays[array_handle] = beta_r
        return beta_r
    
    @beta_r.setter
    def beta_r(self, beta_r):
        self.beta_r[...] = beta_r
    
    @property
    def dbeta_dr(self):
        """
        Element dbeta_dr ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 133
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_pspot__array__dbeta_dr(self._handle)
        if array_handle in self._arrays:
            dbeta_dr = self._arrays[array_handle]
        else:
            dbeta_dr = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_pspot__array__dbeta_dr)
            self._arrays[array_handle] = dbeta_dr
        return dbeta_dr
    
    @dbeta_dr.setter
    def dbeta_dr(self, dbeta_dr):
        self.dbeta_dr[...] = dbeta_dr
    
    @property
    def lden(self):
        """
        Element lden ftype=logical pytype=bool
        
        
        Defined at Struct_module.fpp line 135
        
        """
        return _arespy_pkg.f90wrap_pspot__get__lden(self._handle)
    
    @lden.setter
    def lden(self, lden):
        _arespy_pkg.f90wrap_pspot__set__lden(self._handle, lden)
    
    @property
    def denr(self):
        """
        Element denr ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 136
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_pspot__array__denr(self._handle)
        if array_handle in self._arrays:
            denr = self._arrays[array_handle]
        else:
            denr = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_pspot__array__denr)
            self._arrays[array_handle] = denr
        return denr
    
    @denr.setter
    def denr(self, denr):
        self.denr[...] = denr
    
    @property
    def dden_dr(self):
        """
        Element dden_dr ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 137
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_pspot__array__dden_dr(self._handle)
        if array_handle in self._arrays:
            dden_dr = self._arrays[array_handle]
        else:
            dden_dr = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_pspot__array__dden_dr)
            self._arrays[array_handle] = dden_dr
        return dden_dr
    
    @dden_dr.setter
    def dden_dr(self, dden_dr):
        self.dden_dr[...] = dden_dr
    
    @property
    def lcore(self):
        """
        Element lcore ftype=logical pytype=bool
        
        
        Defined at Struct_module.fpp line 139
        
        """
        return _arespy_pkg.f90wrap_pspot__get__lcore(self._handle)
    
    @lcore.setter
    def lcore(self, lcore):
        _arespy_pkg.f90wrap_pspot__set__lcore(self._handle, lcore)
    
    @property
    def denc(self):
        """
        Element denc ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 141
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_pspot__array__denc(self._handle)
        if array_handle in self._arrays:
            denc = self._arrays[array_handle]
        else:
            denc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_pspot__array__denc)
            self._arrays[array_handle] = denc
        return denc
    
    @denc.setter
    def denc(self, denc):
        self.denc[...] = denc
    
    @property
    def ddenc_dr(self):
        """
        Element ddenc_dr ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 141
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_pspot__array__ddenc_dr(self._handle)
        if array_handle in self._arrays:
            ddenc_dr = self._arrays[array_handle]
        else:
            ddenc_dr = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_pspot__array__ddenc_dr)
            self._arrays[array_handle] = ddenc_dr
        return ddenc_dr
    
    @ddenc_dr.setter
    def ddenc_dr(self, ddenc_dr):
        self.ddenc_dr[...] = ddenc_dr
    
    @property
    def rnoverlap(self):
        """
        Element rnoverlap ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 143
        
        """
        return _arespy_pkg.f90wrap_pspot__get__rnoverlap(self._handle)
    
    @rnoverlap.setter
    def rnoverlap(self, rnoverlap):
        _arespy_pkg.f90wrap_pspot__set__rnoverlap(self._handle, rnoverlap)
    
    @property
    def eps(self):
        """
        Element eps ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 145
        
        """
        return _arespy_pkg.f90wrap_pspot__get__eps(self._handle)
    
    @eps.setter
    def eps(self, eps):
        _arespy_pkg.f90wrap_pspot__set__eps(self._handle, eps)
    
    @property
    def eae(self):
        """
        Element eae ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 145
        
        """
        return _arespy_pkg.f90wrap_pspot__get__eae(self._handle)
    
    @eae.setter
    def eae(self, eae):
        _arespy_pkg.f90wrap_pspot__set__eae(self._handle, eae)
    
    @property
    def nwfa(self):
        """
        Element nwfa ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.fpp line 147
        
        """
        return _arespy_pkg.f90wrap_pspot__get__nwfa(self._handle)
    
    @nwfa.setter
    def nwfa(self, nwfa):
        _arespy_pkg.f90wrap_pspot__set__nwfa(self._handle, nwfa)
    
    @property
    def wfal(self):
        """
        Element wfal ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.fpp line 148
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_pspot__array__wfal(self._handle)
        if array_handle in self._arrays:
            wfal = self._arrays[array_handle]
        else:
            wfal = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_pspot__array__wfal)
            self._arrays[array_handle] = wfal
        return wfal
    
    @wfal.setter
    def wfal(self, wfal):
        self.wfal[...] = wfal
    
    @property
    def wfar(self):
        """
        Element wfar ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 149
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_pspot__array__wfar(self._handle)
        if array_handle in self._arrays:
            wfar = self._arrays[array_handle]
        else:
            wfar = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_pspot__array__wfar)
            self._arrays[array_handle] = wfar
        return wfar
    
    @wfar.setter
    def wfar(self, wfar):
        self.wfar[...] = wfar
    
    def __str__(self):
        ret = ['<pspot>{\n']
        ret.append('    zion : ')
        ret.append(repr(self.zion))
        ret.append(',\n    numps : ')
        ret.append(repr(self.numps))
        ret.append(',\n    qnumps : ')
        ret.append(repr(self.qnumps))
        ret.append(',\n    elename : ')
        ret.append(repr(self.elename))
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
        ret.append(',\n    r : ')
        ret.append(repr(self.r))
        ret.append(',\n    rab : ')
        ret.append(repr(self.rab))
        ret.append(',\n    vlocr : ')
        ret.append(repr(self.vlocr))
        ret.append(',\n    nproj : ')
        ret.append(repr(self.nproj))
        ret.append(',\n    proj_l : ')
        ret.append(repr(self.proj_l))
        ret.append(',\n    proj_m : ')
        ret.append(repr(self.proj_m))
        ret.append(',\n    indx : ')
        ret.append(repr(self.indx))
        ret.append(',\n    rcut : ')
        ret.append(repr(self.rcut))
        ret.append(',\n    rmax : ')
        ret.append(repr(self.rmax))
        ret.append(',\n    rspacing : ')
        ret.append(repr(self.rspacing))
        ret.append(',\n    dij : ')
        ret.append(repr(self.dij))
        ret.append(',\n    beta_r : ')
        ret.append(repr(self.beta_r))
        ret.append(',\n    dbeta_dr : ')
        ret.append(repr(self.dbeta_dr))
        ret.append(',\n    lden : ')
        ret.append(repr(self.lden))
        ret.append(',\n    denr : ')
        ret.append(repr(self.denr))
        ret.append(',\n    dden_dr : ')
        ret.append(repr(self.dden_dr))
        ret.append(',\n    lcore : ')
        ret.append(repr(self.lcore))
        ret.append(',\n    denc : ')
        ret.append(repr(self.denc))
        ret.append(',\n    ddenc_dr : ')
        ret.append(repr(self.ddenc_dr))
        ret.append(',\n    rnoverlap : ')
        ret.append(repr(self.rnoverlap))
        ret.append(',\n    eps : ')
        ret.append(repr(self.eps))
        ret.append(',\n    eae : ')
        ret.append(repr(self.eae))
        ret.append(',\n    nwfa : ')
        ret.append(repr(self.nwfa))
        ret.append(',\n    wfal : ')
        ret.append(repr(self.wfal))
        ret.append(',\n    wfar : ')
        ret.append(repr(self.wfar))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("arespy_pkg.attribute")
class attribute(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=attribute)
    
    
    Defined at Struct_module.fpp lines 162-163
    
    """
    def __init__(self, handle=None):
        """
        self = Attribute()
        
        
        Defined at Struct_module.fpp lines 162-163
        
        
        Returns
        -------
        this : Attribute
        	Object to be constructed
        
        
        Automatically generated constructor for attribute
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _arespy_pkg.f90wrap_attribute_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Attribute
        
        
        Defined at Struct_module.fpp lines 162-163
        
        Parameters
        ----------
        this : Attribute
        	Object to be destructed
        
        
        Automatically generated destructor for attribute
        """
        if self._alloc:
            _arespy_pkg.f90wrap_attribute_finalise(this=self._handle)
    
    @property
    def value(self):
        """
        Element value ftype=character(len=150) pytype=str
        
        
        Defined at Struct_module.fpp line 163
        
        """
        return _arespy_pkg.f90wrap_attribute__get__value(self._handle)
    
    @value.setter
    def value(self, value):
        _arespy_pkg.f90wrap_attribute__set__value(self._handle, value)
    
    def __str__(self):
        ret = ['<attribute>{\n']
        ret.append('    value : ')
        ret.append(repr(self.value))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def read_pspot(nty, filenames):
    """
    read_pspot(nty, filenames)
    
    
    Defined at Struct_module.fpp lines 169-358
    
    Parameters
    ----------
    nty : int
    filenames : str array
    
    """
    _arespy_pkg.f90wrap_read_pspot(nty=nty, filenames=filenames)

def read_psupf_atom(ity, filename, ps):
    """
    read_psupf_atom(ity, filename, ps)
    
    
    Defined at Struct_module.fpp lines 361-578
    
    Parameters
    ----------
    ity : int
    filename : str
    ps : Pspot
    
    -------->Search for Local potential
    """
    _arespy_pkg.f90wrap_read_psupf_atom(ity=ity, filename=filename, ps=ps._handle)

def scan_head(file_unit, title, start_old):
    """
    scan_head(file_unit, title, start_old)
    
    
    Defined at Struct_module.fpp lines 581-639
    
    Parameters
    ----------
    file_unit : int
    title : str
    start_old : bool
    
    """
    _arespy_pkg.f90wrap_scan_head(file_unit=file_unit, title=title, \
        start_old=start_old)

def scan_tail(file_unit, title):
    """
    scan_tail(file_unit, title)
    
    
    Defined at Struct_module.fpp lines 642-665
    
    Parameters
    ----------
    file_unit : int
    title : str
    
    """
    _arespy_pkg.f90wrap_scan_tail(file_unit=file_unit, title=title)

def read_pseudo_header():
    """
    ps, nproj = read_pseudo_header()
    
    
    Defined at Struct_module.fpp lines 668-783
    
    
    Returns
    -------
    ps : Pspot
    nproj : int
    
    """
    ps, nproj = _arespy_pkg.f90wrap_read_pseudo_header()
    ps = f90wrap.runtime.lookup_class("arespy_pkg.pspot").from_handle(ps, \
        alloc=True)
    return ps, nproj

def exist_in(string1, string2):
    """
    exist_in = exist_in(string1, string2)
    
    
    Defined at Struct_module.fpp lines 906-929
    
    Parameters
    ----------
    string1 : str
    string2 : str
    
    Returns
    -------
    exist_in : bool
    
    ====== ====== ======
    """
    exist_in = _arespy_pkg.f90wrap_exist_in(string1=string1, string2=string2)
    return exist_in

def exist_ibegin(string1, string2):
    """
    exist_ibegin = exist_ibegin(string1, string2)
    
    
    Defined at Struct_module.fpp lines 932-955
    
    Parameters
    ----------
    string1 : str
    string2 : str
    
    Returns
    -------
    exist_ibegin : int
    
    ====== ====== ======
    """
    exist_ibegin = _arespy_pkg.f90wrap_exist_ibegin(string1=string1, \
        string2=string2)
    return exist_ibegin

def read_pseudo_nonlocal(unit_upf, nl, beta_r, d0, rcut, proj_l):
    """
    read_pseudo_nonlocal(unit_upf, nl, beta_r, d0, rcut, proj_l)
    
    
    Defined at Struct_module.fpp lines 958-990
    
    Parameters
    ----------
    unit_upf : int
    nl : int
    beta_r : float array
    d0 : float array
    rcut : float
    proj_l : int array
    
    """
    _arespy_pkg.f90wrap_read_pseudo_nonlocal(unit_upf=unit_upf, nl=nl, \
        beta_r=beta_r, d0=d0, rcut=rcut, proj_l=proj_l)

def read_pseudo_pswfc(unit_upf, nwfc, nps, wfcl, wfcr):
    """
    read_pseudo_pswfc(unit_upf, nwfc, nps, wfcl, wfcr)
    
    
    Defined at Struct_module.fpp lines 993-1018
    
    Parameters
    ----------
    unit_upf : int
    nwfc : int
    nps : int
    wfcl : int array
    wfcr : float array
    
    ---------------------------------------------------------------------
    """
    _arespy_pkg.f90wrap_read_pseudo_pswfc(unit_upf=unit_upf, nwfc=nwfc, nps=nps, \
        wfcl=wfcl, wfcr=wfcr)

def aep_generator(nz):
    """
    ps = aep_generator(nz)
    
    
    Defined at Struct_module.fpp lines 1021-1048
    
    Parameters
    ----------
    nz : int
    
    Returns
    -------
    ps : Pspot
    
    """
    ps = _arespy_pkg.f90wrap_aep_generator(nz=nz)
    ps = f90wrap.runtime.lookup_class("arespy_pkg.pspot").from_handle(ps, \
        alloc=True)
    return ps

def _get_value_int(char_in, char_find, variable, find_flag):
    """
    _get_value_int(char_in, char_find, variable, find_flag)
    
    
    Defined at Struct_module.fpp lines 786-813
    
    Parameters
    ----------
    char_in : str
    char_find : str
    variable : int
    find_flag : bool
    
    """
    _arespy_pkg.f90wrap_get_value_int(char_in=char_in, char_find=char_find, \
        variable=variable, find_flag=find_flag)

def _get_value_real(char_in, char_find, variable, find_flag):
    """
    _get_value_real(char_in, char_find, variable, find_flag)
    
    
    Defined at Struct_module.fpp lines 816-843
    
    Parameters
    ----------
    char_in : str
    char_find : str
    variable : float
    find_flag : bool
    
    """
    _arespy_pkg.f90wrap_get_value_real(char_in=char_in, char_find=char_find, \
        variable=variable, find_flag=find_flag)

def _get_value_char(char_in, char_find, variable, find_flag):
    """
    _get_value_char(char_in, char_find, variable, find_flag)
    
    
    Defined at Struct_module.fpp lines 846-873
    
    Parameters
    ----------
    char_in : str
    char_find : str
    variable : str
    find_flag : bool
    
    """
    _arespy_pkg.f90wrap_get_value_char(char_in=char_in, char_find=char_find, \
        variable=variable, find_flag=find_flag)

def _get_value_logic(char_in, char_find, variable, find_flag):
    """
    _get_value_logic(char_in, char_find, variable, find_flag)
    
    
    Defined at Struct_module.fpp lines 876-903
    
    Parameters
    ----------
    char_in : str
    char_find : str
    variable : bool
    find_flag : bool
    
    """
    _arespy_pkg.f90wrap_get_value_logic(char_in=char_in, char_find=char_find, \
        variable=variable, find_flag=find_flag)

def get_value(*args, **kwargs):
    """
    get_value(*args, **kwargs)
    
    
    Defined at Struct_module.fpp lines 158-160
    
    Overloaded interface containing the following procedures:
      _get_value_int
      _get_value_real
      _get_value_char
      _get_value_logic
    
    """
    for proc in [_get_value_int, _get_value_real, _get_value_char, \
        _get_value_logic]:
        try:
            return proc(*args, **kwargs)
        except TypeError:
            continue

def get_max_nproj():
    """
    Element max_nproj ftype=integer(i4b) pytype=int
    
    
    Defined at Struct_module.fpp line 154
    
    """
    return _arespy_pkg.f90wrap_pspot_module__get__max_nproj()

def set_max_nproj(max_nproj):
    _arespy_pkg.f90wrap_pspot_module__set__max_nproj(max_nproj)

def get_max_nwfa():
    """
    Element max_nwfa ftype=integer(i4b) pytype=int
    
    
    Defined at Struct_module.fpp line 154
    
    """
    return _arespy_pkg.f90wrap_pspot_module__get__max_nwfa()

def set_max_nwfa(max_nwfa):
    _arespy_pkg.f90wrap_pspot_module__set__max_nwfa(max_nwfa)

def get_max_rcut():
    """
    Element max_rcut ftype=real(dp) pytype=float
    
    
    Defined at Struct_module.fpp line 155
    
    """
    return _arespy_pkg.f90wrap_pspot_module__get__max_rcut()

def set_max_rcut(max_rcut):
    _arespy_pkg.f90wrap_pspot_module__set__max_rcut(max_rcut)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "pspot_module".')

for func in _dt_array_initialisers:
    func()

"""
Module nlpot_module


Defined at Nonlocalpot_module.fpp lines 5-1371

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("AresMainPy_pkg.nol_type")
class nol_type(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=nol_type)
    
    
    Defined at Nonlocalpot_module.fpp lines 14-27
    
    """
    def __init__(self, handle=None):
        """
        self = Nol_Type()
        
        
        Defined at Nonlocalpot_module.fpp lines 14-27
        
        
        Returns
        -------
        this : Nol_Type
        	Object to be constructed
        
        
        Automatically generated constructor for nol_type
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _AresMainPy_pkg.f90wrap_nol_type_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Nol_Type
        
        
        Defined at Nonlocalpot_module.fpp lines 14-27
        
        Parameters
        ----------
        this : Nol_Type
        	Object to be destructed
        
        
        Automatically generated destructor for nol_type
        """
        if self._alloc:
            _AresMainPy_pkg.f90wrap_nol_type_finalise(this=self._handle)
    
    @property
    def npts(self):
        """
        Element npts ftype=integer(i4b) pytype=int
        
        
        Defined at Nonlocalpot_module.fpp line 15
        
        """
        return _AresMainPy_pkg.f90wrap_nol_type__get__npts(self._handle)
    
    @npts.setter
    def npts(self, npts):
        _AresMainPy_pkg.f90wrap_nol_type__set__npts(self._handle, npts)
    
    @property
    def id(self):
        """
        Element id ftype=integer(i4b) pytype=int
        
        
        Defined at Nonlocalpot_module.fpp line 16
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_nol_type__array__id(self._handle)
        if array_handle in self._arrays:
            id = self._arrays[array_handle]
        else:
            id = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_nol_type__array__id)
            self._arrays[array_handle] = id
        return id
    
    @id.setter
    def id(self, id):
        self.id[...] = id
    
    @property
    def rrvec(self):
        """
        Element rrvec ftype=real(dp) pytype=float
        
        
        Defined at Nonlocalpot_module.fpp line 17
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_nol_type__array__rrvec(self._handle)
        if array_handle in self._arrays:
            rrvec = self._arrays[array_handle]
        else:
            rrvec = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_nol_type__array__rrvec)
            self._arrays[array_handle] = rrvec
        return rrvec
    
    @rrvec.setter
    def rrvec(self, rrvec):
        self.rrvec[...] = rrvec
    
    @property
    def proj0(self):
        """
        Element proj0 ftype=real(dp) pytype=float
        
        
        Defined at Nonlocalpot_module.fpp line 18
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_nol_type__array__proj0(self._handle)
        if array_handle in self._arrays:
            proj0 = self._arrays[array_handle]
        else:
            proj0 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_nol_type__array__proj0)
            self._arrays[array_handle] = proj0
        return proj0
    
    @proj0.setter
    def proj0(self, proj0):
        self.proj0[...] = proj0
    
    @property
    def proj(self):
        """
        Element proj ftype=real(dp) pytype=float
        
        
        Defined at Nonlocalpot_module.fpp line 19
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_nol_type__array__proj(self._handle)
        if array_handle in self._arrays:
            proj = self._arrays[array_handle]
        else:
            proj = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_nol_type__array__proj)
            self._arrays[array_handle] = proj
        return proj
    
    @proj.setter
    def proj(self, proj):
        self.proj[...] = proj
    
    @property
    def proj_phs(self):
        """
        Element proj_phs ftype=complex(dp) pytype=complex
        
        
        Defined at Nonlocalpot_module.fpp line 20
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_nol_type__array__proj_phs(self._handle)
        if array_handle in self._arrays:
            proj_phs = self._arrays[array_handle]
        else:
            proj_phs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_nol_type__array__proj_phs)
            self._arrays[array_handle] = proj_phs
        return proj_phs
    
    @proj_phs.setter
    def proj_phs(self, proj_phs):
        self.proj_phs[...] = proj_phs
    
    @property
    def proj0_dg(self):
        """
        Element proj0_dg ftype=real(dp) pytype=float
        
        
        Defined at Nonlocalpot_module.fpp line 22
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_nol_type__array__proj0_dg(self._handle)
        if array_handle in self._arrays:
            proj0_dg = self._arrays[array_handle]
        else:
            proj0_dg = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_nol_type__array__proj0_dg)
            self._arrays[array_handle] = proj0_dg
        return proj0_dg
    
    @proj0_dg.setter
    def proj0_dg(self, proj0_dg):
        self.proj0_dg[...] = proj0_dg
    
    @property
    def proj_dg(self):
        """
        Element proj_dg ftype=real(dp) pytype=float
        
        
        Defined at Nonlocalpot_module.fpp line 23
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_nol_type__array__proj_dg(self._handle)
        if array_handle in self._arrays:
            proj_dg = self._arrays[array_handle]
        else:
            proj_dg = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_nol_type__array__proj_dg)
            self._arrays[array_handle] = proj_dg
        return proj_dg
    
    @proj_dg.setter
    def proj_dg(self, proj_dg):
        self.proj_dg[...] = proj_dg
    
    @property
    def proj_phs_dg(self):
        """
        Element proj_phs_dg ftype=complex(dp) pytype=complex
        
        
        Defined at Nonlocalpot_module.fpp line 24
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_nol_type__array__proj_phs_dg(self._handle)
        if array_handle in self._arrays:
            proj_phs_dg = self._arrays[array_handle]
        else:
            proj_phs_dg = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_nol_type__array__proj_phs_dg)
            self._arrays[array_handle] = proj_phs_dg
        return proj_phs_dg
    
    @proj_phs_dg.setter
    def proj_phs_dg(self, proj_phs_dg):
        self.proj_phs_dg[...] = proj_phs_dg
    
    @property
    def id_iso(self):
        """
        Element id_iso ftype=integer(i4b) pytype=int
        
        
        Defined at Nonlocalpot_module.fpp line 26
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_nol_type__array__id_iso(self._handle)
        if array_handle in self._arrays:
            id_iso = self._arrays[array_handle]
        else:
            id_iso = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_nol_type__array__id_iso)
            self._arrays[array_handle] = id_iso
        return id_iso
    
    @id_iso.setter
    def id_iso(self, id_iso):
        self.id_iso[...] = id_iso
    
    @property
    def rrvec_iso(self):
        """
        Element rrvec_iso ftype=real(dp) pytype=float
        
        
        Defined at Nonlocalpot_module.fpp line 27
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_nol_type__array__rrvec_iso(self._handle)
        if array_handle in self._arrays:
            rrvec_iso = self._arrays[array_handle]
        else:
            rrvec_iso = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_nol_type__array__rrvec_iso)
            self._arrays[array_handle] = rrvec_iso
        return rrvec_iso
    
    @rrvec_iso.setter
    def rrvec_iso(self, rrvec_iso):
        self.rrvec_iso[...] = rrvec_iso
    
    def __str__(self):
        ret = ['<nol_type>{\n']
        ret.append('    npts : ')
        ret.append(repr(self.npts))
        ret.append(',\n    id : ')
        ret.append(repr(self.id))
        ret.append(',\n    rrvec : ')
        ret.append(repr(self.rrvec))
        ret.append(',\n    proj0 : ')
        ret.append(repr(self.proj0))
        ret.append(',\n    proj : ')
        ret.append(repr(self.proj))
        ret.append(',\n    proj_phs : ')
        ret.append(repr(self.proj_phs))
        ret.append(',\n    proj0_dg : ')
        ret.append(repr(self.proj0_dg))
        ret.append(',\n    proj_dg : ')
        ret.append(repr(self.proj_dg))
        ret.append(',\n    proj_phs_dg : ')
        ret.append(repr(self.proj_phs_dg))
        ret.append(',\n    id_iso : ')
        ret.append(repr(self.id_iso))
        ret.append(',\n    rrvec_iso : ')
        ret.append(repr(self.rrvec_iso))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def initialize_nlpot():
    """
    initialize_nlpot()
    
    
    Defined at Nonlocalpot_module.fpp lines 37-44
    
    
    """
    _AresMainPy_pkg.f90wrap_initialize_nlpot()

def initialize_nlpot_per():
    """
    initialize_nlpot_per()
    
    
    Defined at Nonlocalpot_module.fpp lines 47-191
    
    
    ===================store the data================
    """
    _AresMainPy_pkg.f90wrap_initialize_nlpot_per()

def initialize_nlpot_iso():
    """
    initialize_nlpot_iso()
    
    
    Defined at Nonlocalpot_module.fpp lines 194-326
    
    
    ===================store the data================
    """
    _AresMainPy_pkg.f90wrap_initialize_nlpot_iso()

def destroy_nlpot():
    """
    destroy_nlpot()
    
    
    Defined at Nonlocalpot_module.fpp lines 329-385
    
    
    """
    _AresMainPy_pkg.f90wrap_destroy_nlpot()

def set_beta_real():
    """
    set_beta_real()
    
    
    Defined at Nonlocalpot_module.fpp lines 388-451
    
    
    ======================
    ##PLOT
    """
    _AresMainPy_pkg.f90wrap_set_beta_real()

def nlp_beta_interp_r(ity, ia, beta_init):
    """
    nlp_beta_interp_r(ity, ia, beta_init)
    
    
    Defined at Nonlocalpot_module.fpp lines 454-516
    
    Parameters
    ----------
    ity : int
    ia : int
    beta_init : float array
    
    """
    _AresMainPy_pkg.f90wrap_nlp_beta_interp_r(ity=ity, ia=ia, beta_init=beta_init)

def nlp_beta_ylm_r(ity, ia, beta_ylm):
    """
    nlp_beta_ylm_r(ity, ia, beta_ylm)
    
    
    Defined at Nonlocalpot_module.fpp lines 519-561
    
    Parameters
    ----------
    ity : int
    ia : int
    beta_ylm : float array
    
    """
    _AresMainPy_pkg.f90wrap_nlp_beta_ylm_r(ity=ity, ia=ia, beta_ylm=beta_ylm)

def nlp_beta_phase_r(ity, ia, beta_ylm, beta_phase):
    """
    nlp_beta_phase_r(ity, ia, beta_ylm, beta_phase)
    
    
    Defined at Nonlocalpot_module.fpp lines 564-591
    
    Parameters
    ----------
    ity : int
    ia : int
    beta_ylm : float array
    beta_phase : complex array
    
    """
    _AresMainPy_pkg.f90wrap_nlp_beta_phase_r(ity=ity, ia=ia, beta_ylm=beta_ylm, \
        beta_phase=beta_phase)

def apply_ylm(l, m, fac, x, y, z, f):
    """
    apply_ylm(l, m, fac, x, y, z, f)
    
    
    Defined at Nonlocalpot_module.fpp lines 594-670
    
    Parameters
    ----------
    l : int
    m : int
    fac : float
    x : float
    y : float
    z : float
    f : float
    
    """
    _AresMainPy_pkg.f90wrap_apply_ylm(l=l, m=m, fac=fac, x=x, y=y, z=z, f=f)

def nlp_wijk_dg(ndg, inpol, numdg, lft, rit, wijk, drijk):
    """
    nlp_wijk_dg(ndg, inpol, numdg, lft, rit, wijk, drijk)
    
    
    Defined at Nonlocalpot_module.fpp lines 673-786
    
    Parameters
    ----------
    ndg : int
    inpol : int
    numdg : int
    lft : int
    rit : int
    wijk : float array
    drijk : float array
    
    """
    _AresMainPy_pkg.f90wrap_nlp_wijk_dg(ndg=ndg, inpol=inpol, numdg=numdg, lft=lft, \
        rit=rit, wijk=wijk, drijk=drijk)

def nlp_wijk_dg_old(ndg, n_near, numdg, lft, rit, wijk, drijk):
    """
    nlp_wijk_dg_old(ndg, n_near, numdg, lft, rit, wijk, drijk)
    
    
    Defined at Nonlocalpot_module.fpp lines 788-891
    
    Parameters
    ----------
    ndg : int
    n_near : int
    numdg : int
    lft : int
    rit : int
    wijk : float array
    drijk : float array
    
    """
    _AresMainPy_pkg.f90wrap_nlp_wijk_dg_old(ndg=ndg, n_near=n_near, numdg=numdg, \
        lft=lft, rit=rit, wijk=wijk, drijk=drijk)

def nlp_interp_betalm_dg(ity, ia, ip, numdg, lft, rit, wijk, drijk, beta0_ijk, \
    betalm_ijk):
    """
    nlp_interp_betalm_dg(ity, ia, ip, numdg, lft, rit, wijk, drijk, beta0_ijk, \
        betalm_ijk)
    
    
    Defined at Nonlocalpot_module.fpp lines 894-995
    
    Parameters
    ----------
    ity : int
    ia : int
    ip : int
    numdg : int
    lft : int
    rit : int
    wijk : float array
    drijk : float array
    beta0_ijk : float array
    betalm_ijk : float array
    
    =============================
    ##FOR ISO
    """
    _AresMainPy_pkg.f90wrap_nlp_interp_betalm_dg(ity=ity, ia=ia, ip=ip, numdg=numdg, \
        lft=lft, rit=rit, wijk=wijk, drijk=drijk, beta0_ijk=beta0_ijk, \
        betalm_ijk=betalm_ijk)

def nlp_interp_betalm_dg_iso(ity, ia, ip, numdg, lft, rit, wijk, drijk, \
    beta0_ijk, betalm_ijk):
    """
    nlp_interp_betalm_dg_iso(ity, ia, ip, numdg, lft, rit, wijk, drijk, beta0_ijk, \
        betalm_ijk)
    
    
    Defined at Nonlocalpot_module.fpp lines 998-1104
    
    Parameters
    ----------
    ity : int
    ia : int
    ip : int
    numdg : int
    lft : int
    rit : int
    wijk : float array
    drijk : float array
    beta0_ijk : float array
    betalm_ijk : float array
    
    =============================
    ##FOR ISO
    """
    _AresMainPy_pkg.f90wrap_nlp_interp_betalm_dg_iso(ity=ity, ia=ia, ip=ip, \
        numdg=numdg, lft=lft, rit=rit, wijk=wijk, drijk=drijk, beta0_ijk=beta0_ijk, \
        betalm_ijk=betalm_ijk)

def nlp_beta_ylm_r_dg(ity, ia, numdg, lft, rit, wijk, drijk, beta0, beta_ylm):
    """
    nlp_beta_ylm_r_dg(ity, ia, numdg, lft, rit, wijk, drijk, beta0, beta_ylm)
    
    
    Defined at Nonlocalpot_module.fpp lines 1107-1122
    
    Parameters
    ----------
    ity : int
    ia : int
    numdg : int
    lft : int
    rit : int
    wijk : float array
    drijk : float array
    beta0 : float array
    beta_ylm : float array
    
    """
    _AresMainPy_pkg.f90wrap_nlp_beta_ylm_r_dg(ity=ity, ia=ia, numdg=numdg, lft=lft, \
        rit=rit, wijk=wijk, drijk=drijk, beta0=beta0, beta_ylm=beta_ylm)

def nlp_beta_ylm_r_dg_iso(ity, ia, numdg, lft, rit, wijk, drijk, beta0, \
    beta_ylm):
    """
    nlp_beta_ylm_r_dg_iso(ity, ia, numdg, lft, rit, wijk, drijk, beta0, beta_ylm)
    
    
    Defined at Nonlocalpot_module.fpp lines 1125-1161
    
    Parameters
    ----------
    ity : int
    ia : int
    numdg : int
    lft : int
    rit : int
    wijk : float array
    drijk : float array
    beta0 : float array
    beta_ylm : float array
    
    """
    _AresMainPy_pkg.f90wrap_nlp_beta_ylm_r_dg_iso(ity=ity, ia=ia, numdg=numdg, \
        lft=lft, rit=rit, wijk=wijk, drijk=drijk, beta0=beta0, beta_ylm=beta_ylm)

def set_beta_real_dg():
    """
    set_beta_real_dg()
    
    
    Defined at Nonlocalpot_module.fpp lines 1164-1220
    
    
    """
    _AresMainPy_pkg.f90wrap_set_beta_real_dg()

def initialize_nlpot_band():
    """
    initialize_nlpot_band()
    
    
    Defined at Nonlocalpot_module.fpp lines 1222-1370
    
    
    ===================store the data================
    """
    _AresMainPy_pkg.f90wrap_initialize_nlpot_band()

def get_max_nlnpts():
    """
    Element max_nlnpts ftype=integer(i4b) pytype=int
    
    
    Defined at Nonlocalpot_module.fpp line 32
    
    """
    return _AresMainPy_pkg.f90wrap_nlpot_module__get__max_nlnpts()

def set_max_nlnpts(max_nlnpts):
    _AresMainPy_pkg.f90wrap_nlpot_module__set__max_nlnpts(max_nlnpts)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "nlpot_module".')

for func in _dt_array_initialisers:
    func()

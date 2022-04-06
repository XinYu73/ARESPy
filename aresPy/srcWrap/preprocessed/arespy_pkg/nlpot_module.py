"""
Module nlpot_module


Defined at Nonlocalpot_module.fpp lines 5-483

"""
from __future__ import print_function, absolute_import, division
import _arespy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("arespy_pkg.nol_type")
class nol_type(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=nol_type)
    
    
    Defined at Nonlocalpot_module.fpp lines 14-24
    
    """
    def __init__(self, handle=None):
        """
        self = Nol_Type()
        
        
        Defined at Nonlocalpot_module.fpp lines 14-24
        
        
        Returns
        -------
        this : Nol_Type
        	Object to be constructed
        
        
        Automatically generated constructor for nol_type
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _arespy_pkg.f90wrap_nol_type_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Nol_Type
        
        
        Defined at Nonlocalpot_module.fpp lines 14-24
        
        Parameters
        ----------
        this : Nol_Type
        	Object to be destructed
        
        
        Automatically generated destructor for nol_type
        """
        if self._alloc:
            _arespy_pkg.f90wrap_nol_type_finalise(this=self._handle)
    
    @property
    def npts(self):
        """
        Element npts ftype=integer(i4b) pytype=int
        
        
        Defined at Nonlocalpot_module.fpp line 15
        
        """
        return _arespy_pkg.f90wrap_nol_type__get__npts(self._handle)
    
    @npts.setter
    def npts(self, npts):
        _arespy_pkg.f90wrap_nol_type__set__npts(self._handle, npts)
    
    @property
    def id(self):
        """
        Element id ftype=integer(i4b) pytype=int
        
        
        Defined at Nonlocalpot_module.fpp line 16
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_nol_type__array__id(self._handle)
        if array_handle in self._arrays:
            id = self._arrays[array_handle]
        else:
            id = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_nol_type__array__id)
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
            _arespy_pkg.f90wrap_nol_type__array__rrvec(self._handle)
        if array_handle in self._arrays:
            rrvec = self._arrays[array_handle]
        else:
            rrvec = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_nol_type__array__rrvec)
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
            _arespy_pkg.f90wrap_nol_type__array__proj0(self._handle)
        if array_handle in self._arrays:
            proj0 = self._arrays[array_handle]
        else:
            proj0 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_nol_type__array__proj0)
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
            _arespy_pkg.f90wrap_nol_type__array__proj(self._handle)
        if array_handle in self._arrays:
            proj = self._arrays[array_handle]
        else:
            proj = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_nol_type__array__proj)
            self._arrays[array_handle] = proj
        return proj
    
    @proj.setter
    def proj(self, proj):
        self.proj[...] = proj
    
    @property
    def proj_phs(self):
        """
        Element proj_phs ftype=complex(dcp) pytype=complex
        
        
        Defined at Nonlocalpot_module.fpp line 20
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_nol_type__array__proj_phs(self._handle)
        if array_handle in self._arrays:
            proj_phs = self._arrays[array_handle]
        else:
            proj_phs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_nol_type__array__proj_phs)
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
            _arespy_pkg.f90wrap_nol_type__array__proj0_dg(self._handle)
        if array_handle in self._arrays:
            proj0_dg = self._arrays[array_handle]
        else:
            proj0_dg = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_nol_type__array__proj0_dg)
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
            _arespy_pkg.f90wrap_nol_type__array__proj_dg(self._handle)
        if array_handle in self._arrays:
            proj_dg = self._arrays[array_handle]
        else:
            proj_dg = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_nol_type__array__proj_dg)
            self._arrays[array_handle] = proj_dg
        return proj_dg
    
    @proj_dg.setter
    def proj_dg(self, proj_dg):
        self.proj_dg[...] = proj_dg
    
    @property
    def proj_phs_dg(self):
        """
        Element proj_phs_dg ftype=complex(dcp) pytype=complex
        
        
        Defined at Nonlocalpot_module.fpp line 24
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_nol_type__array__proj_phs_dg(self._handle)
        if array_handle in self._arrays:
            proj_phs_dg = self._arrays[array_handle]
        else:
            proj_phs_dg = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_nol_type__array__proj_phs_dg)
            self._arrays[array_handle] = proj_phs_dg
        return proj_phs_dg
    
    @proj_phs_dg.setter
    def proj_phs_dg(self, proj_phs_dg):
        self.proj_phs_dg[...] = proj_phs_dg
    
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
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def initialize_nlpot():
    """
    initialize_nlpot()
    
    
    Defined at Nonlocalpot_module.fpp lines 32-178
    
    
    ===================store the data================
    """
    _arespy_pkg.f90wrap_initialize_nlpot()

def destroy_nlpot():
    """
    destroy_nlpot()
    
    
    Defined at Nonlocalpot_module.fpp lines 181-200
    
    
    """
    _arespy_pkg.f90wrap_destroy_nlpot()

def set_beta_real():
    """
    set_beta_real()
    
    
    Defined at Nonlocalpot_module.fpp lines 203-234
    
    
    ============apply Ylm in real space===============
    cycle the all species
    """
    _arespy_pkg.f90wrap_set_beta_real()

def nlp_beta_interp_r(ity, ia, beta_init):
    """
    nlp_beta_interp_r(ity, ia, beta_init)
    
    
    Defined at Nonlocalpot_module.fpp lines 237-268
    
    Parameters
    ----------
    ity : int
    ia : int
    beta_init : float array
    
    """
    _arespy_pkg.f90wrap_nlp_beta_interp_r(ity=ity, ia=ia, beta_init=beta_init)

def nlp_beta_ylm_r(ity, ia, beta_ylm):
    """
    nlp_beta_ylm_r(ity, ia, beta_ylm)
    
    
    Defined at Nonlocalpot_module.fpp lines 271-308
    
    Parameters
    ----------
    ity : int
    ia : int
    beta_ylm : float array
    
    """
    _arespy_pkg.f90wrap_nlp_beta_ylm_r(ity=ity, ia=ia, beta_ylm=beta_ylm)

def nlp_beta_phase_r(ity, ia, beta_ylm, beta_phase):
    """
    nlp_beta_phase_r(ity, ia, beta_ylm, beta_phase)
    
    
    Defined at Nonlocalpot_module.fpp lines 311-347
    
    Parameters
    ----------
    ity : int
    ia : int
    beta_ylm : float array
    beta_phase : complex array
    
    """
    _arespy_pkg.f90wrap_nlp_beta_phase_r(ity=ity, ia=ia, beta_ylm=beta_ylm, \
        beta_phase=beta_phase)

def apply_ylm(l, m, fac, x, y, z, f):
    """
    apply_ylm(l, m, fac, x, y, z, f)
    
    
    Defined at Nonlocalpot_module.fpp lines 350-426
    
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
    _arespy_pkg.f90wrap_apply_ylm(l=l, m=m, fac=fac, x=x, y=y, z=z, f=f)

def nlp_init_partialcore(rhoc):
    """
    nlp_init_partialcore(rhoc)
    
    
    Defined at Nonlocalpot_module.fpp lines 429-482
    
    Parameters
    ----------
    rhoc : float array
    
    """
    _arespy_pkg.f90wrap_nlp_init_partialcore(rhoc=rhoc)

def get_max_nlnpts():
    """
    Element max_nlnpts ftype=integer(i4b) pytype=int
    
    
    Defined at Nonlocalpot_module.fpp line 29
    
    """
    return _arespy_pkg.f90wrap_nlpot_module__get__max_nlnpts()

def set_max_nlnpts(max_nlnpts):
    _arespy_pkg.f90wrap_nlpot_module__set__max_nlnpts(max_nlnpts)


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

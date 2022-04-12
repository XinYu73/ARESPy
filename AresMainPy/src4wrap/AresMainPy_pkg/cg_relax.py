"""
Module cg_relax


Defined at cg_module.fpp lines 5-289

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("AresMainPy_pkg.lattice")
class lattice(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=lattice)
    
    
    Defined at cg_module.fpp lines 16-20
    
    """
    def __init__(self, handle=None):
        """
        self = Lattice()
        
        
        Defined at cg_module.fpp lines 16-20
        
        
        Returns
        -------
        this : Lattice
        	Object to be constructed
        
        
        Automatically generated constructor for lattice
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _AresMainPy_pkg.f90wrap_lattice_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Lattice
        
        
        Defined at cg_module.fpp lines 16-20
        
        Parameters
        ----------
        this : Lattice
        	Object to be destructed
        
        
        Automatically generated destructor for lattice
        """
        if self._alloc:
            _AresMainPy_pkg.f90wrap_lattice_finalise(this=self._handle)
    
    @property
    def a(self):
        """
        Element a ftype=real(dp) pytype=float
        
        
        Defined at cg_module.fpp line 17
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_lattice__array__a(self._handle)
        if array_handle in self._arrays:
            a = self._arrays[array_handle]
        else:
            a = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_lattice__array__a)
            self._arrays[array_handle] = a
        return a
    
    @a.setter
    def a(self, a):
        self.a[...] = a
    
    @property
    def b(self):
        """
        Element b ftype=real(dp) pytype=float
        
        
        Defined at cg_module.fpp line 18
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_lattice__array__b(self._handle)
        if array_handle in self._arrays:
            b = self._arrays[array_handle]
        else:
            b = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_lattice__array__b)
            self._arrays[array_handle] = b
        return b
    
    @b.setter
    def b(self, b):
        self.b[...] = b
    
    @property
    def anorm(self):
        """
        Element anorm ftype=real(dp) pytype=float
        
        
        Defined at cg_module.fpp line 19
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_lattice__array__anorm(self._handle)
        if array_handle in self._arrays:
            anorm = self._arrays[array_handle]
        else:
            anorm = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_lattice__array__anorm)
            self._arrays[array_handle] = anorm
        return anorm
    
    @anorm.setter
    def anorm(self, anorm):
        self.anorm[...] = anorm
    
    @property
    def bnorm(self):
        """
        Element bnorm ftype=real(dp) pytype=float
        
        
        Defined at cg_module.fpp line 19
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_lattice__array__bnorm(self._handle)
        if array_handle in self._arrays:
            bnorm = self._arrays[array_handle]
        else:
            bnorm = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_lattice__array__bnorm)
            self._arrays[array_handle] = bnorm
        return bnorm
    
    @bnorm.setter
    def bnorm(self, bnorm):
        self.bnorm[...] = bnorm
    
    @property
    def omega(self):
        """
        Element omega ftype=real(dp) pytype=float
        
        
        Defined at cg_module.fpp line 20
        
        """
        return _AresMainPy_pkg.f90wrap_lattice__get__omega(self._handle)
    
    @omega.setter
    def omega(self, omega):
        _AresMainPy_pkg.f90wrap_lattice__set__omega(self._handle, omega)
    
    def __str__(self):
        ret = ['<lattice>{\n']
        ret.append('    a : ')
        ret.append(repr(self.a))
        ret.append(',\n    b : ')
        ret.append(repr(self.b))
        ret.append(',\n    anorm : ')
        ret.append(repr(self.anorm))
        ret.append(',\n    bnorm : ')
        ret.append(repr(self.bnorm))
        ret.append(',\n    omega : ')
        ret.append(repr(self.omega))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def cg_relax_vasp_interface(nstep, lfopt, lopt, lexceed):
    """
    cg_relax_vasp_interface(nstep, lfopt, lopt, lexceed)
    
    
    Defined at cg_module.fpp lines 23-196
    
    Parameters
    ----------
    nstep : int
    lfopt : bool
    lopt : bool
    lexceed : bool
    
    ----------------------------------------------------------------
    """
    _AresMainPy_pkg.f90wrap_cg_relax_vasp_interface(nstep=nstep, lfopt=lfopt, \
        lopt=lopt, lexceed=lexceed)

def check_opt(na, flag):
    """
    check_opt(na, flag)
    
    
    Defined at cg_module.fpp lines 199-236
    
    Parameters
    ----------
    na : int
    flag : bool
    
    """
    _AresMainPy_pkg.f90wrap_check_opt(na=na, flag=flag)

def lattic(self):
    """
    lattic(self)
    
    
    Defined at cg_module.fpp lines 239-259
    
    Parameters
    ----------
    mylatt : Lattice
    
    """
    _AresMainPy_pkg.f90wrap_lattic(mylatt=self._handle)

def expro(h, u1, u2):
    """
    expro(h, u1, u2)
    
    
    Defined at cg_module.fpp lines 262-267
    
    Parameters
    ----------
    h : float array
    u1 : float array
    u2 : float array
    
    """
    _AresMainPy_pkg.f90wrap_expro(h=h, u1=u1, u2=u2)

def check_distance(l, pos, lwrong):
    """
    check_distance(l, pos, lwrong)
    
    
    Defined at cg_module.fpp lines 270-287
    
    Parameters
    ----------
    l : float
    pos : float array
    lwrong : bool
    
    -----------------------------------------------
    """
    _AresMainPy_pkg.f90wrap_check_distance(l=l, pos=pos, lwrong=lwrong)

def get_maxforce():
    """
    Element maxforce ftype=real(dp) pytype=float
    
    
    Defined at cg_module.fpp line 13
    
    """
    return _AresMainPy_pkg.f90wrap_cg_relax__get__maxforce()

def set_maxforce(maxforce):
    _AresMainPy_pkg.f90wrap_cg_relax__set__maxforce(maxforce)

def get_maxstress():
    """
    Element maxstress ftype=real(dp) pytype=float
    
    
    Defined at cg_module.fpp line 13
    
    """
    return _AresMainPy_pkg.f90wrap_cg_relax__get__maxstress()

def set_maxstress(maxstress):
    _AresMainPy_pkg.f90wrap_cg_relax__set__maxstress(maxstress)

def get_array_tstress():
    """
    Element tstress ftype=real(dp) pytype=float
    
    
    Defined at cg_module.fpp line 14
    
    """
    global tstress
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_cg_relax__array__tstress(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        tstress = _arrays[array_handle]
    else:
        tstress = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_cg_relax__array__tstress)
        _arrays[array_handle] = tstress
    return tstress

def set_array_tstress(tstress):
    tstress[...] = tstress

def get_array_latmark():
    """
    Element latmark ftype=real(dp) pytype=float
    
    
    Defined at cg_module.fpp line 15
    
    """
    global latmark
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_cg_relax__array__latmark(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        latmark = _arrays[array_handle]
    else:
        latmark = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_cg_relax__array__latmark)
        _arrays[array_handle] = latmark
    return latmark

def set_array_latmark(latmark):
    latmark[...] = latmark


_array_initialisers = [get_array_tstress, get_array_latmark]
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "cg_relax".')

for func in _dt_array_initialisers:
    func()

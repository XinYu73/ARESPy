"""
Module dyn_module


Defined at dyn_module.fpp lines 5-44

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("AresMainPy_pkg.dynamics")
class dynamics(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=dynamics)
    
    
    Defined at dyn_module.fpp lines 8-23
    
    """
    def __init__(self, handle=None):
        """
        self = Dynamics()
        
        
        Defined at dyn_module.fpp lines 8-23
        
        
        Returns
        -------
        this : Dynamics
        	Object to be constructed
        
        
        Automatically generated constructor for dynamics
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _AresMainPy_pkg.f90wrap_dynamics_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Dynamics
        
        
        Defined at dyn_module.fpp lines 8-23
        
        Parameters
        ----------
        this : Dynamics
        	Object to be destructed
        
        
        Automatically generated destructor for dynamics
        """
        if self._alloc:
            _AresMainPy_pkg.f90wrap_dynamics_finalise(this=self._handle)
    
    @property
    def lstop(self):
        """
        Element lstop ftype=logical pytype=bool
        
        
        Defined at dyn_module.fpp line 9
        
        """
        return _AresMainPy_pkg.f90wrap_dynamics__get__lstop(self._handle)
    
    @lstop.setter
    def lstop(self, lstop):
        _AresMainPy_pkg.f90wrap_dynamics__set__lstop(self._handle, lstop)
    
    @property
    def posion(self):
        """
        Element posion ftype=real(dp) pytype=float
        
        
        Defined at dyn_module.fpp line 10
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_dynamics__array__posion(self._handle)
        if array_handle in self._arrays:
            posion = self._arrays[array_handle]
        else:
            posion = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_dynamics__array__posion)
            self._arrays[array_handle] = posion
        return posion
    
    @posion.setter
    def posion(self, posion):
        self.posion[...] = posion
    
    @property
    def posioc(self):
        """
        Element posioc ftype=real(dp) pytype=float
        
        
        Defined at dyn_module.fpp line 11
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_dynamics__array__posioc(self._handle)
        if array_handle in self._arrays:
            posioc = self._arrays[array_handle]
        else:
            posioc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_dynamics__array__posioc)
            self._arrays[array_handle] = posioc
        return posioc
    
    @posioc.setter
    def posioc(self, posioc):
        self.posioc[...] = posioc
    
    @property
    def d2(self):
        """
        Element d2 ftype=real(dp) pytype=float
        
        
        Defined at dyn_module.fpp line 12
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_dynamics__array__d2(self._handle)
        if array_handle in self._arrays:
            d2 = self._arrays[array_handle]
        else:
            d2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_dynamics__array__d2)
            self._arrays[array_handle] = d2
        return d2
    
    @d2.setter
    def d2(self, d2):
        self.d2[...] = d2
    
    @property
    def d2c(self):
        """
        Element d2c ftype=real(dp) pytype=float
        
        
        Defined at dyn_module.fpp line 13
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_dynamics__array__d2c(self._handle)
        if array_handle in self._arrays:
            d2c = self._arrays[array_handle]
        else:
            d2c = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_dynamics__array__d2c)
            self._arrays[array_handle] = d2c
        return d2c
    
    @d2c.setter
    def d2c(self, d2c):
        self.d2c[...] = d2c
    
    @property
    def d3(self):
        """
        Element d3 ftype=real(dp) pytype=float
        
        
        Defined at dyn_module.fpp line 14
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_dynamics__array__d3(self._handle)
        if array_handle in self._arrays:
            d3 = self._arrays[array_handle]
        else:
            d3 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_dynamics__array__d3)
            self._arrays[array_handle] = d3
        return d3
    
    @d3.setter
    def d3(self, d3):
        self.d3[...] = d3
    
    @property
    def a(self):
        """
        Element a ftype=real(dp) pytype=float
        
        
        Defined at dyn_module.fpp line 15
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_dynamics__array__a(self._handle)
        if array_handle in self._arrays:
            a = self._arrays[array_handle]
        else:
            a = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_dynamics__array__a)
            self._arrays[array_handle] = a
        return a
    
    @a.setter
    def a(self, a):
        self.a[...] = a
    
    @property
    def b(self):
        """
        Element b ftype=real(dp) pytype=float
        
        
        Defined at dyn_module.fpp line 16
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_dynamics__array__b(self._handle)
        if array_handle in self._arrays:
            b = self._arrays[array_handle]
        else:
            b = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_dynamics__array__b)
            self._arrays[array_handle] = b
        return b
    
    @b.setter
    def b(self, b):
        self.b[...] = b
    
    @property
    def ac(self):
        """
        Element ac ftype=real(dp) pytype=float
        
        
        Defined at dyn_module.fpp line 17
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_dynamics__array__ac(self._handle)
        if array_handle in self._arrays:
            ac = self._arrays[array_handle]
        else:
            ac = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_dynamics__array__ac)
            self._arrays[array_handle] = ac
        return ac
    
    @ac.setter
    def ac(self, ac):
        self.ac[...] = ac
    
    @property
    def potim(self):
        """
        Element potim ftype=real(dp) pytype=float
        
        
        Defined at dyn_module.fpp line 18
        
        """
        return _AresMainPy_pkg.f90wrap_dynamics__get__potim(self._handle)
    
    @potim.setter
    def potim(self, potim):
        _AresMainPy_pkg.f90wrap_dynamics__set__potim(self._handle, potim)
    
    @property
    def ediffg(self):
        """
        Element ediffg ftype=real(dp) pytype=float
        
        
        Defined at dyn_module.fpp line 19
        
        """
        return _AresMainPy_pkg.f90wrap_dynamics__get__ediffg(self._handle)
    
    @ediffg.setter
    def ediffg(self, ediffg):
        _AresMainPy_pkg.f90wrap_dynamics__set__ediffg(self._handle, ediffg)
    
    @property
    def pstress(self):
        """
        Element pstress ftype=real(dp) pytype=float
        
        
        Defined at dyn_module.fpp line 20
        
        """
        return _AresMainPy_pkg.f90wrap_dynamics__get__pstress(self._handle)
    
    @pstress.setter
    def pstress(self, pstress):
        _AresMainPy_pkg.f90wrap_dynamics__set__pstress(self._handle, pstress)
    
    @property
    def ibrion(self):
        """
        Element ibrion ftype=integer(i4b) pytype=int
        
        
        Defined at dyn_module.fpp line 21
        
        """
        return _AresMainPy_pkg.f90wrap_dynamics__get__ibrion(self._handle)
    
    @ibrion.setter
    def ibrion(self, ibrion):
        _AresMainPy_pkg.f90wrap_dynamics__set__ibrion(self._handle, ibrion)
    
    @property
    def isif(self):
        """
        Element isif ftype=integer(i4b) pytype=int
        
        
        Defined at dyn_module.fpp line 22
        
        """
        return _AresMainPy_pkg.f90wrap_dynamics__get__isif(self._handle)
    
    @isif.setter
    def isif(self, isif):
        _AresMainPy_pkg.f90wrap_dynamics__set__isif(self._handle, isif)
    
    @property
    def nfree(self):
        """
        Element nfree ftype=integer(i4b) pytype=int
        
        
        Defined at dyn_module.fpp line 23
        
        """
        return _AresMainPy_pkg.f90wrap_dynamics__get__nfree(self._handle)
    
    @nfree.setter
    def nfree(self, nfree):
        _AresMainPy_pkg.f90wrap_dynamics__set__nfree(self._handle, nfree)
    
    def __str__(self):
        ret = ['<dynamics>{\n']
        ret.append('    lstop : ')
        ret.append(repr(self.lstop))
        ret.append(',\n    posion : ')
        ret.append(repr(self.posion))
        ret.append(',\n    posioc : ')
        ret.append(repr(self.posioc))
        ret.append(',\n    d2 : ')
        ret.append(repr(self.d2))
        ret.append(',\n    d2c : ')
        ret.append(repr(self.d2c))
        ret.append(',\n    d3 : ')
        ret.append(repr(self.d3))
        ret.append(',\n    a : ')
        ret.append(repr(self.a))
        ret.append(',\n    b : ')
        ret.append(repr(self.b))
        ret.append(',\n    ac : ')
        ret.append(repr(self.ac))
        ret.append(',\n    potim : ')
        ret.append(repr(self.potim))
        ret.append(',\n    ediffg : ')
        ret.append(repr(self.ediffg))
        ret.append(',\n    pstress : ')
        ret.append(repr(self.pstress))
        ret.append(',\n    ibrion : ')
        ret.append(repr(self.ibrion))
        ret.append(',\n    isif : ')
        ret.append(repr(self.isif))
        ret.append(',\n    nfree : ')
        ret.append(repr(self.nfree))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def create_dyn(na, dyn):
    """
    create_dyn(na, dyn)
    
    
    Defined at dyn_module.fpp lines 27-35
    
    Parameters
    ----------
    na : int
    dyn : Dynamics
    
    """
    _AresMainPy_pkg.f90wrap_create_dyn(na=na, dyn=dyn._handle)

def destroy_dyn(self):
    """
    destroy_dyn(self)
    
    
    Defined at dyn_module.fpp lines 37-44
    
    Parameters
    ----------
    dyn : Dynamics
    
    """
    _AresMainPy_pkg.f90wrap_destroy_dyn(dyn=self._handle)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "dyn_module".')

for func in _dt_array_initialisers:
    func()

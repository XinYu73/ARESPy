"""
Module datatypes_allocatable


Defined at datatypes.fpp lines 6-27

"""
from __future__ import print_function, absolute_import, division
import _ExampleDerivedTypes_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("ExampleDerivedTypes_pkg.alloc_arrays")
class alloc_arrays(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=alloc_arrays)
    
    
    Defined at datatypes.fpp lines 10-14
    
    """
    def __init__(self, handle=None):
        """
        self = Alloc_Arrays()
        
        
        Defined at datatypes.fpp lines 10-14
        
        
        Returns
        -------
        this : Alloc_Arrays
        	Object to be constructed
        
        
        Automatically generated constructor for alloc_arrays
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _ExampleDerivedTypes_pkg.f90wrap_alloc_arrays_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Alloc_Arrays
        
        
        Defined at datatypes.fpp lines 10-14
        
        Parameters
        ----------
        this : Alloc_Arrays
        	Object to be destructed
        
        
        Automatically generated destructor for alloc_arrays
        """
        if self._alloc:
            _ExampleDerivedTypes_pkg.f90wrap_alloc_arrays_finalise(this=self._handle)
    
    @property
    def chi(self):
        """
        Element chi ftype=real(idp) pytype=float
        
        
        Defined at datatypes.fpp line 11
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _ExampleDerivedTypes_pkg.f90wrap_alloc_arrays__array__chi(self._handle)
        if array_handle in self._arrays:
            chi = self._arrays[array_handle]
        else:
            chi = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _ExampleDerivedTypes_pkg.f90wrap_alloc_arrays__array__chi)
            self._arrays[array_handle] = chi
        return chi
    
    @chi.setter
    def chi(self, chi):
        self.chi[...] = chi
    
    @property
    def psi(self):
        """
        Element psi ftype=real(idp) pytype=float
        
        
        Defined at datatypes.fpp line 12
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _ExampleDerivedTypes_pkg.f90wrap_alloc_arrays__array__psi(self._handle)
        if array_handle in self._arrays:
            psi = self._arrays[array_handle]
        else:
            psi = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _ExampleDerivedTypes_pkg.f90wrap_alloc_arrays__array__psi)
            self._arrays[array_handle] = psi
        return psi
    
    @psi.setter
    def psi(self, psi):
        self.psi[...] = psi
    
    @property
    def chi_shape(self):
        """
        Element chi_shape ftype=integer(4) pytype=int
        
        
        Defined at datatypes.fpp line 13
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _ExampleDerivedTypes_pkg.f90wrap_alloc_arrays__array__chi_shape(self._handle)
        if array_handle in self._arrays:
            chi_shape = self._arrays[array_handle]
        else:
            chi_shape = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _ExampleDerivedTypes_pkg.f90wrap_alloc_arrays__array__chi_shape)
            self._arrays[array_handle] = chi_shape
        return chi_shape
    
    @chi_shape.setter
    def chi_shape(self, chi_shape):
        self.chi_shape[...] = chi_shape
    
    @property
    def psi_shape(self):
        """
        Element psi_shape ftype=integer(4) pytype=int
        
        
        Defined at datatypes.fpp line 14
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _ExampleDerivedTypes_pkg.f90wrap_alloc_arrays__array__psi_shape(self._handle)
        if array_handle in self._arrays:
            psi_shape = self._arrays[array_handle]
        else:
            psi_shape = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _ExampleDerivedTypes_pkg.f90wrap_alloc_arrays__array__psi_shape)
            self._arrays[array_handle] = psi_shape
        return psi_shape
    
    @psi_shape.setter
    def psi_shape(self, psi_shape):
        self.psi_shape[...] = psi_shape
    
    def __str__(self):
        ret = ['<alloc_arrays>{\n']
        ret.append('    chi : ')
        ret.append(repr(self.chi))
        ret.append(',\n    psi : ')
        ret.append(repr(self.psi))
        ret.append(',\n    chi_shape : ')
        ret.append(repr(self.chi_shape))
        ret.append(',\n    psi_shape : ')
        ret.append(repr(self.psi_shape))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def init_alloc_arrays(self, m, n):
    """
    init_alloc_arrays(self, m, n)
    
    
    Defined at datatypes.fpp lines 17-21
    
    Parameters
    ----------
    dertype : Alloc_Arrays
    m : int
    n : int
    
    """
    _ExampleDerivedTypes_pkg.f90wrap_init_alloc_arrays(dertype=self._handle, m=m, \
        n=n)

def destroy_alloc_arrays(self):
    """
    destroy_alloc_arrays(self)
    
    
    Defined at datatypes.fpp lines 23-26
    
    Parameters
    ----------
    dertype : Alloc_Arrays
    
    """
    _ExampleDerivedTypes_pkg.f90wrap_destroy_alloc_arrays(dertype=self._handle)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "datatypes_allocatable".')

for func in _dt_array_initialisers:
    func()

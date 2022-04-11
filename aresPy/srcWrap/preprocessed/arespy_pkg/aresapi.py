"""
Module aresapi


Defined at aresAPI.fpp lines 5-55

"""
from __future__ import print_function, absolute_import, division
import _arespy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("arespy_pkg.aresOut")
class aresOut(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=aresout)
    
    
    Defined at aresAPI.fpp lines 9-14
    
    """
    def __init__(self, handle=None):
        """
        self = Aresout()
        
        
        Defined at aresAPI.fpp lines 9-14
        
        
        Returns
        -------
        this : Aresout
        	Object to be constructed
        
        
        Automatically generated constructor for aresout
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _arespy_pkg.f90wrap_aresout_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Aresout
        
        
        Defined at aresAPI.fpp lines 9-14
        
        Parameters
        ----------
        this : Aresout
        	Object to be destructed
        
        
        Automatically generated destructor for aresout
        """
        if self._alloc:
            _arespy_pkg.f90wrap_aresout_finalise(this=self._handle)
    
    @property
    def forces(self):
        """
        Element forces ftype=real(dp) pytype=float
        
        
        Defined at aresAPI.fpp line 10
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_aresout__array__forces(self._handle)
        if array_handle in self._arrays:
            forces = self._arrays[array_handle]
        else:
            forces = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_aresout__array__forces)
            self._arrays[array_handle] = forces
        return forces
    
    @forces.setter
    def forces(self, forces):
        self.forces[...] = forces
    
    @property
    def stress(self):
        """
        Element stress ftype=real(dp) pytype=float
        
        
        Defined at aresAPI.fpp line 11
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_aresout__array__stress(self._handle)
        if array_handle in self._arrays:
            stress = self._arrays[array_handle]
        else:
            stress = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_aresout__array__stress)
            self._arrays[array_handle] = stress
        return stress
    
    @stress.setter
    def stress(self, stress):
        self.stress[...] = stress
    
    @property
    def poscar(self):
        """
        Element poscar ftype=real(dp) pytype=float
        
        
        Defined at aresAPI.fpp line 12
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_aresout__array__poscar(self._handle)
        if array_handle in self._arrays:
            poscar = self._arrays[array_handle]
        else:
            poscar = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_aresout__array__poscar)
            self._arrays[array_handle] = poscar
        return poscar
    
    @poscar.setter
    def poscar(self, poscar):
        self.poscar[...] = poscar
    
    @property
    def pos(self):
        """
        Element pos ftype=real(dp) pytype=float
        
        
        Defined at aresAPI.fpp line 13
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_aresout__array__pos(self._handle)
        if array_handle in self._arrays:
            pos = self._arrays[array_handle]
        else:
            pos = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_aresout__array__pos)
            self._arrays[array_handle] = pos
        return pos
    
    @pos.setter
    def pos(self, pos):
        self.pos[...] = pos
    
    @property
    def chargerho(self):
        """
        Element chargerho ftype=real(dp) pytype=float
        
        
        Defined at aresAPI.fpp line 14
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_aresout__array__chargerho(self._handle)
        if array_handle in self._arrays:
            chargerho = self._arrays[array_handle]
        else:
            chargerho = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_aresout__array__chargerho)
            self._arrays[array_handle] = chargerho
        return chargerho
    
    @chargerho.setter
    def chargerho(self, chargerho):
        self.chargerho[...] = chargerho
    
    def __str__(self):
        ret = ['<aresout>{\n']
        ret.append('    forces : ')
        ret.append(repr(self.forces))
        ret.append(',\n    stress : ')
        ret.append(repr(self.stress))
        ret.append(',\n    poscar : ')
        ret.append(repr(self.poscar))
        ret.append(',\n    pos : ')
        ret.append(repr(self.pos))
        ret.append(',\n    chargerho : ')
        ret.append(repr(self.chargerho))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def init_alloc_arrays(self, nnatom):
    """
    init_alloc_arrays(self, nnatom)
    
    
    Defined at aresAPI.fpp lines 17-47
    
    Parameters
    ----------
    dertype : Aresout
    nnatom : int
    
    """
    _arespy_pkg.f90wrap_init_alloc_arrays(dertype=self._handle, nnatom=nnatom)

def destroy_alloc_arrays(self):
    """
    destroy_alloc_arrays(self)
    
    
    Defined at aresAPI.fpp lines 49-55
    
    Parameters
    ----------
    dertype : Aresout
    
    """
    _arespy_pkg.f90wrap_destroy_alloc_arrays(dertype=self._handle)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "aresapi".')

for func in _dt_array_initialisers:
    func()

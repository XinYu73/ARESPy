"""
Module relax_module


Defined at Relax_module.fpp lines 10-305

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("AresMainPy_pkg.relax_type")
class relax_type(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=relax_type)
    
    
    Defined at Relax_module.fpp lines 17-29
    
    """
    def __init__(self, handle=None):
        """
        self = Relax_Type()
        
        
        Defined at Relax_module.fpp lines 17-29
        
        
        Returns
        -------
        this : Relax_Type
        	Object to be constructed
        
        
        Automatically generated constructor for relax_type
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _AresMainPy_pkg.f90wrap_relax_type_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Relax_Type
        
        
        Defined at Relax_module.fpp lines 17-29
        
        Parameters
        ----------
        this : Relax_Type
        	Object to be destructed
        
        
        Automatically generated destructor for relax_type
        """
        if self._alloc:
            _AresMainPy_pkg.f90wrap_relax_type_finalise(this=self._handle)
    
    @property
    def timsmax(self):
        """
        Element timsmax ftype=real(dp) pytype=float
        
        
        Defined at Relax_module.fpp line 21
        
        """
        return _AresMainPy_pkg.f90wrap_relax_type__get__timsmax(self._handle)
    
    @timsmax.setter
    def timsmax(self, timsmax):
        _AresMainPy_pkg.f90wrap_relax_type__set__timsmax(self._handle, timsmax)
    
    @property
    def tims(self):
        """
        Element tims ftype=real(dp) pytype=float
        
        
        Defined at Relax_module.fpp line 21
        
        """
        return _AresMainPy_pkg.f90wrap_relax_type__get__tims(self._handle)
    
    @tims.setter
    def tims(self, tims):
        _AresMainPy_pkg.f90wrap_relax_type__set__tims(self._handle, tims)
    
    @property
    def alpha(self):
        """
        Element alpha ftype=real(dp) pytype=float
        
        
        Defined at Relax_module.fpp line 21
        
        """
        return _AresMainPy_pkg.f90wrap_relax_type__get__alpha(self._handle)
    
    @alpha.setter
    def alpha(self, alpha):
        _AresMainPy_pkg.f90wrap_relax_type__set__alpha(self._handle, alpha)
    
    @property
    def faction(self):
        """
        Element faction ftype=real(dp) pytype=float
        
        
        Defined at Relax_module.fpp line 24
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_relax_type__array__faction(self._handle)
        if array_handle in self._arrays:
            faction = self._arrays[array_handle]
        else:
            faction = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_relax_type__array__faction)
            self._arrays[array_handle] = faction
        return faction
    
    @faction.setter
    def faction(self, faction):
        self.faction[...] = faction
    
    @property
    def fiond(self):
        """
        Element fiond ftype=real(dp) pytype=float
        
        
        Defined at Relax_module.fpp line 24
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_relax_type__array__fiond(self._handle)
        if array_handle in self._arrays:
            fiond = self._arrays[array_handle]
        else:
            fiond = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_relax_type__array__fiond)
            self._arrays[array_handle] = fiond
        return fiond
    
    @fiond.setter
    def fiond(self, fiond):
        self.fiond[...] = fiond
    
    @property
    def fceld(self):
        """
        Element fceld ftype=real(dp) pytype=float
        
        
        Defined at Relax_module.fpp line 24
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_relax_type__array__fceld(self._handle)
        if array_handle in self._arrays:
            fceld = self._arrays[array_handle]
        else:
            fceld = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_relax_type__array__fceld)
            self._arrays[array_handle] = fceld
        return fceld
    
    @fceld.setter
    def fceld(self, fceld):
        self.fceld[...] = fceld
    
    @property
    def velion(self):
        """
        Element velion ftype=real(dp) pytype=float
        
        
        Defined at Relax_module.fpp line 25
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_relax_type__array__velion(self._handle)
        if array_handle in self._arrays:
            velion = self._arrays[array_handle]
        else:
            velion = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_relax_type__array__velion)
            self._arrays[array_handle] = velion
        return velion
    
    @velion.setter
    def velion(self, velion):
        self.velion[...] = velion
    
    @property
    def velcel(self):
        """
        Element velcel ftype=real(dp) pytype=float
        
        
        Defined at Relax_module.fpp line 26
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_relax_type__array__velcel(self._handle)
        if array_handle in self._arrays:
            velcel = self._arrays[array_handle]
        else:
            velcel = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_relax_type__array__velcel)
            self._arrays[array_handle] = velcel
        return velcel
    
    @velcel.setter
    def velcel(self, velcel):
        self.velcel[...] = velcel
    
    @property
    def factcel(self):
        """
        Element factcel ftype=real(dp) pytype=float
        
        
        Defined at Relax_module.fpp line 27
        
        """
        return _AresMainPy_pkg.f90wrap_relax_type__get__factcel(self._handle)
    
    @factcel.setter
    def factcel(self, factcel):
        _AresMainPy_pkg.f90wrap_relax_type__set__factcel(self._handle, factcel)
    
    @property
    def neg(self):
        """
        Element neg ftype=integer(i4b) pytype=int
        
        
        Defined at Relax_module.fpp line 28
        
        """
        return _AresMainPy_pkg.f90wrap_relax_type__get__neg(self._handle)
    
    @neg.setter
    def neg(self, neg):
        _AresMainPy_pkg.f90wrap_relax_type__set__neg(self._handle, neg)
    
    @property
    def lneg(self):
        """
        Element lneg ftype=logical pytype=bool
        
        
        Defined at Relax_module.fpp line 29
        
        """
        return _AresMainPy_pkg.f90wrap_relax_type__get__lneg(self._handle)
    
    @lneg.setter
    def lneg(self, lneg):
        _AresMainPy_pkg.f90wrap_relax_type__set__lneg(self._handle, lneg)
    
    def __str__(self):
        ret = ['<relax_type>{\n']
        ret.append('    timsmax : ')
        ret.append(repr(self.timsmax))
        ret.append(',\n    tims : ')
        ret.append(repr(self.tims))
        ret.append(',\n    alpha : ')
        ret.append(repr(self.alpha))
        ret.append(',\n    faction : ')
        ret.append(repr(self.faction))
        ret.append(',\n    fiond : ')
        ret.append(repr(self.fiond))
        ret.append(',\n    fceld : ')
        ret.append(repr(self.fceld))
        ret.append(',\n    velion : ')
        ret.append(repr(self.velion))
        ret.append(',\n    velcel : ')
        ret.append(repr(self.velcel))
        ret.append(',\n    factcel : ')
        ret.append(repr(self.factcel))
        ret.append(',\n    neg : ')
        ret.append(repr(self.neg))
        ret.append(',\n    lneg : ')
        ret.append(repr(self.lneg))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def initialize_relax():
    """
    initialize_relax()
    
    
    Defined at Relax_module.fpp lines 36-80
    
    
    """
    _AresMainPy_pkg.f90wrap_initialize_relax()

def destroy_relax():
    """
    destroy_relax()
    
    
    Defined at Relax_module.fpp lines 83-107
    
    
    """
    _AresMainPy_pkg.f90wrap_destroy_relax()

def relaxer(nstep):
    """
    relaxer(nstep)
    
    
    Defined at Relax_module.fpp lines 110-136
    
    Parameters
    ----------
    nstep : int
    
    """
    _AresMainPy_pkg.f90wrap_relaxer(nstep=nstep)

def fire_relax(lat, pos, fion, fcel):
    """
    fire_relax(lat, pos, fion, fcel)
    
    
    Defined at Relax_module.fpp lines 139-261
    
    Parameters
    ----------
    lat : float array
    pos : float array
    fion : float array
    fcel : float array
    
    ========================================================
    """
    _AresMainPy_pkg.f90wrap_fire_relax(lat=lat, pos=pos, fion=fion, fcel=fcel)

def force_on_cell(stress, lat, fcel):
    """
    force_on_cell(stress, lat, fcel)
    
    
    Defined at Relax_module.fpp lines 264-281
    
    Parameters
    ----------
    stress : float array
    lat : float array
    fcel : float array
    
    """
    _AresMainPy_pkg.f90wrap_force_on_cell(stress=stress, lat=lat, fcel=fcel)

def fix_direction(stress):
    """
    fix_direction(stress)
    
    
    Defined at Relax_module.fpp lines 284-304
    
    Parameters
    ----------
    stress : float array
    
    """
    _AresMainPy_pkg.f90wrap_fix_direction(stress=stress)

def get_pstress():
    """
    Element pstress ftype=real(dp) pytype=float
    
    
    Defined at Relax_module.fpp line 14
    
    """
    return _AresMainPy_pkg.f90wrap_relax_module__get__pstress()

def set_pstress(pstress):
    _AresMainPy_pkg.f90wrap_relax_module__set__pstress(pstress)

def get_ldone():
    """
    Element ldone ftype=logical pytype=bool
    
    
    Defined at Relax_module.fpp line 16
    
    """
    return _AresMainPy_pkg.f90wrap_relax_module__get__ldone()

def set_ldone(ldone):
    _AresMainPy_pkg.f90wrap_relax_module__set__ldone(ldone)

def get_lfirst():
    """
    Element lfirst ftype=logical pytype=bool
    
    
    Defined at Relax_module.fpp line 16
    
    """
    return _AresMainPy_pkg.f90wrap_relax_module__get__lfirst()

def set_lfirst(lfirst):
    _AresMainPy_pkg.f90wrap_relax_module__set__lfirst(lfirst)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "relax_module".')

for func in _dt_array_initialisers:
    func()

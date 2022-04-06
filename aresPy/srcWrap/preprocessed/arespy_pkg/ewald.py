"""
Module ewald


Defined at Ewald.fpp lines 5-826

"""
from __future__ import print_function, absolute_import, division
import _arespy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("arespy_pkg.ion")
class ion(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=ion)
    
    
    Defined at Ewald.fpp lines 36-38
    
    """
    def __init__(self, handle=None):
        """
        self = Ion()
        
        
        Defined at Ewald.fpp lines 36-38
        
        
        Returns
        -------
        this : Ion
        	Object to be constructed
        
        
        Automatically generated constructor for ion
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _arespy_pkg.f90wrap_ion_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Ion
        
        
        Defined at Ewald.fpp lines 36-38
        
        Parameters
        ----------
        this : Ion
        	Object to be destructed
        
        
        Automatically generated destructor for ion
        """
        if self._alloc:
            _arespy_pkg.f90wrap_ion_finalise(this=self._handle)
    
    @property
    def charge(self):
        """
        Element charge ftype=real(dp) pytype=float
        
        
        Defined at Ewald.fpp line 37
        
        """
        return _arespy_pkg.f90wrap_ion__get__charge(self._handle)
    
    @charge.setter
    def charge(self, charge):
        _arespy_pkg.f90wrap_ion__set__charge(self._handle, charge)
    
    @property
    def fcd(self):
        """
        Element fcd ftype=real(dp) pytype=float
        
        
        Defined at Ewald.fpp line 38
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_ion__array__fcd(self._handle)
        if array_handle in self._arrays:
            fcd = self._arrays[array_handle]
        else:
            fcd = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_ion__array__fcd)
            self._arrays[array_handle] = fcd
        return fcd
    
    @fcd.setter
    def fcd(self, fcd):
        self.fcd[...] = fcd
    
    def __str__(self):
        ret = ['<ion>{\n']
        ret.append('    charge : ')
        ret.append(repr(self.charge))
        ret.append(',\n    fcd : ')
        ret.append(repr(self.fcd))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def ewald_energy(latticev, ionpositions, iontpid, ioncharges):
    """
    ewald_energy = ewald_energy(latticev, ionpositions, iontpid, ioncharges)
    
    
    Defined at Ewald.fpp lines 58-76
    
    Parameters
    ----------
    latticev : float array
    ionpositions : float array
    iontpid : int array
    ioncharges : float array
    
    Returns
    -------
    ewald_energy : float
    
    """
    ewald_energy = _arespy_pkg.f90wrap_ewald_energy(latticev=latticev, \
        ionpositions=ionpositions, iontpid=iontpid, ioncharges=ioncharges)
    return ewald_energy

def iso_ewald_energy(latticev, ionpositions, iontpid, ioncharges):
    """
    iso_ewald_energy = iso_ewald_energy(latticev, ionpositions, iontpid, ioncharges)
    
    
    Defined at Ewald.fpp lines 80-108
    
    Parameters
    ----------
    latticev : float array
    ionpositions : float array
    iontpid : int array
    ioncharges : float array
    
    Returns
    -------
    iso_ewald_energy : float
    
    """
    iso_ewald_energy = _arespy_pkg.f90wrap_iso_ewald_energy(latticev=latticev, \
        ionpositions=ionpositions, iontpid=iontpid, ioncharges=ioncharges)
    return iso_ewald_energy

def iso_ewald_forces(latticev, ionpositions, iontpid, ioncharges):
    """
    iso_ewald_forces = iso_ewald_forces(latticev, ionpositions, iontpid, ioncharges)
    
    
    Defined at Ewald.fpp lines 112-193
    
    Parameters
    ----------
    latticev : float array
    ionpositions : float array
    iontpid : int array
    ioncharges : float array
    
    Returns
    -------
    iso_ewald_forces : float array
    
    ==========================================
    """
    iso_ewald_forces = _arespy_pkg.f90wrap_iso_ewald_forces(latticev=latticev, \
        ionpositions=ionpositions, iontpid=iontpid, ioncharges=ioncharges)
    return iso_ewald_forces

def ewald_forces(latticev, ionpositions, iontpid, ioncharges):
    """
    ewald_forces = ewald_forces(latticev, ionpositions, iontpid, ioncharges)
    
    
    Defined at Ewald.fpp lines 196-219
    
    Parameters
    ----------
    latticev : float array
    ionpositions : float array
    iontpid : int array
    ioncharges : float array
    
    Returns
    -------
    ewald_forces : float array
    
    """
    ewald_forces = _arespy_pkg.f90wrap_ewald_forces(latticev=latticev, \
        ionpositions=ionpositions, iontpid=iontpid, ioncharges=ioncharges)
    return ewald_forces

def ewald_stress(latticev, ionpositions, iontpid, ioncharges):
    """
    ewald_stress = ewald_stress(latticev, ionpositions, iontpid, ioncharges)
    
    
    Defined at Ewald.fpp lines 222-240
    
    Parameters
    ----------
    latticev : float array
    ionpositions : float array
    iontpid : int array
    ioncharges : float array
    
    Returns
    -------
    ewald_stress : float array
    
    """
    ewald_stress = _arespy_pkg.f90wrap_ewald_stress(latticev=latticev, \
        ionpositions=ionpositions, iontpid=iontpid, ioncharges=ioncharges)
    return ewald_stress

def vectorlength(vc):
    """
    vectorlength = vectorlength(vc)
    
    
    Defined at Ewald.fpp lines 715-716
    
    Parameters
    ----------
    vc : float array
    
    Returns
    -------
    vectorlength : float
    
    """
    vectorlength = _arespy_pkg.f90wrap_vectorlength(vc=vc)
    return vectorlength

def recipvector(lat):
    """
    recipvector = recipvector(lat)
    
    
    Defined at Ewald.fpp lines 719-724
    
    Parameters
    ----------
    lat : float array
    
    Returns
    -------
    recipvector : float array
    
    """
    recipvector = _arespy_pkg.f90wrap_recipvector(lat=lat)
    return recipvector

def volume(lat):
    """
    volume = volume(lat)
    
    
    Defined at Ewald.fpp lines 727-729
    
    Parameters
    ----------
    lat : float array
    
    Returns
    -------
    volume : float
    
    """
    volume = _arespy_pkg.f90wrap_volume(lat=lat)
    return volume

def crossp(va, vb):
    """
    crossp = crossp(va, vb)
    
    
    Defined at Ewald.fpp lines 732-736
    
    Parameters
    ----------
    va : float array
    vb : float array
    
    Returns
    -------
    crossp : float array
    
    """
    crossp = _arespy_pkg.f90wrap_crossp(va=va, vb=vb)
    return crossp

def erfc(x):
    """
    erfc = erfc(x)
    
    
    Defined at Ewald.fpp lines 739-825
    
    Parameters
    ----------
    x : float
    
    Returns
    -------
    erfc : float
    
    """
    erfc = _arespy_pkg.f90wrap_erfc(x=x)
    return erfc

def get_bohr():
    """
    Element bohr ftype=real(dp) pytype=float
    
    
    Defined at Ewald.fpp line 40
    
    """
    return _arespy_pkg.f90wrap_ewald__get__bohr()

bohr = get_bohr()

def get_hartreetoev():
    """
    Element hartreetoev ftype=real(dp) pytype=float
    
    
    Defined at Ewald.fpp line 40
    
    """
    return _arespy_pkg.f90wrap_ewald__get__hartreetoev()

hartreetoev = get_hartreetoev()


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "ewald".')

for func in _dt_array_initialisers:
    func()

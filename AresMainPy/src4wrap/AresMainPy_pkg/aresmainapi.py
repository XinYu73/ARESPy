"""
Module aresmainapi


Defined at AresMainAPI.fpp lines 5-77

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("AresMainPy_pkg.aresOut")
class aresOut(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=aresout)
    
    
    Defined at AresMainAPI.fpp lines 10-21
    
    """
    def __init__(self, handle=None):
        """
        self = Aresout()
        
        
        Defined at AresMainAPI.fpp lines 10-21
        
        
        Returns
        -------
        this : Aresout
        	Object to be constructed
        
        
        Automatically generated constructor for aresout
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _AresMainPy_pkg.f90wrap_aresout_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Aresout
        
        
        Defined at AresMainAPI.fpp lines 10-21
        
        Parameters
        ----------
        this : Aresout
        	Object to be destructed
        
        
        Automatically generated destructor for aresout
        """
        if self._alloc:
            _AresMainPy_pkg.f90wrap_aresout_finalise(this=self._handle)
    
    @property
    def forces(self):
        """
        Element forces ftype=real(dp) pytype=float
        
        
        Defined at AresMainAPI.fpp line 11
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_aresout__array__forces(self._handle)
        if array_handle in self._arrays:
            forces = self._arrays[array_handle]
        else:
            forces = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_aresout__array__forces)
            self._arrays[array_handle] = forces
        return forces
    
    @forces.setter
    def forces(self, forces):
        self.forces[...] = forces
    
    @property
    def stress(self):
        """
        Element stress ftype=real(dp) pytype=float
        
        
        Defined at AresMainAPI.fpp line 12
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_aresout__array__stress(self._handle)
        if array_handle in self._arrays:
            stress = self._arrays[array_handle]
        else:
            stress = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_aresout__array__stress)
            self._arrays[array_handle] = stress
        return stress
    
    @stress.setter
    def stress(self, stress):
        self.stress[...] = stress
    
    @property
    def poscar(self):
        """
        Element poscar ftype=real(dp) pytype=float
        
        
        Defined at AresMainAPI.fpp line 13
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_aresout__array__poscar(self._handle)
        if array_handle in self._arrays:
            poscar = self._arrays[array_handle]
        else:
            poscar = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_aresout__array__poscar)
            self._arrays[array_handle] = poscar
        return poscar
    
    @poscar.setter
    def poscar(self, poscar):
        self.poscar[...] = poscar
    
    @property
    def pos(self):
        """
        Element pos ftype=real(dp) pytype=float
        
        
        Defined at AresMainAPI.fpp line 14
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_aresout__array__pos(self._handle)
        if array_handle in self._arrays:
            pos = self._arrays[array_handle]
        else:
            pos = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_aresout__array__pos)
            self._arrays[array_handle] = pos
        return pos
    
    @pos.setter
    def pos(self, pos):
        self.pos[...] = pos
    
    @property
    def chargerho(self):
        """
        Element chargerho ftype=real(dp) pytype=float
        
        
        Defined at AresMainAPI.fpp line 15
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_aresout__array__chargerho(self._handle)
        if array_handle in self._arrays:
            chargerho = self._arrays[array_handle]
        else:
            chargerho = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_aresout__array__chargerho)
            self._arrays[array_handle] = chargerho
        return chargerho
    
    @chargerho.setter
    def chargerho(self, chargerho):
        self.chargerho[...] = chargerho
    
    @property
    def apilat_mat(self):
        """
        Element apilat_mat ftype=real(dp) pytype=float
        
        
        Defined at AresMainAPI.fpp line 16
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_aresout__array__apilat_mat(self._handle)
        if array_handle in self._arrays:
            apilat_mat = self._arrays[array_handle]
        else:
            apilat_mat = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_aresout__array__apilat_mat)
            self._arrays[array_handle] = apilat_mat
        return apilat_mat
    
    @apilat_mat.setter
    def apilat_mat(self, apilat_mat):
        self.apilat_mat[...] = apilat_mat
    
    @property
    def apilat_para(self):
        """
        Element apilat_para ftype=real(dp) pytype=float
        
        
        Defined at AresMainAPI.fpp line 17
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_aresout__array__apilat_para(self._handle)
        if array_handle in self._arrays:
            apilat_para = self._arrays[array_handle]
        else:
            apilat_para = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_aresout__array__apilat_para)
            self._arrays[array_handle] = apilat_para
        return apilat_para
    
    @apilat_para.setter
    def apilat_para(self, apilat_para):
        self.apilat_para[...] = apilat_para
    
    @property
    def comm(self):
        """
        Element comm ftype=integer(i4b) pytype=int
        
        
        Defined at AresMainAPI.fpp line 18
        
        """
        return _AresMainPy_pkg.f90wrap_aresout__get__comm(self._handle)
    
    @comm.setter
    def comm(self, comm):
        _AresMainPy_pkg.f90wrap_aresout__set__comm(self._handle, comm)
    
    @property
    def myid(self):
        """
        Element myid ftype=integer(i4b) pytype=int
        
        
        Defined at AresMainAPI.fpp line 19
        
        """
        return _AresMainPy_pkg.f90wrap_aresout__get__myid(self._handle)
    
    @myid.setter
    def myid(self, myid):
        _AresMainPy_pkg.f90wrap_aresout__set__myid(self._handle, myid)
    
    @property
    def numprocs(self):
        """
        Element numprocs ftype=integer(i4b) pytype=int
        
        
        Defined at AresMainAPI.fpp line 20
        
        """
        return _AresMainPy_pkg.f90wrap_aresout__get__numprocs(self._handle)
    
    @numprocs.setter
    def numprocs(self, numprocs):
        _AresMainPy_pkg.f90wrap_aresout__set__numprocs(self._handle, numprocs)
    
    @property
    def rootid(self):
        """
        Element rootid ftype=integer(i4b) pytype=int
        
        
        Defined at AresMainAPI.fpp line 21
        
        """
        return _AresMainPy_pkg.f90wrap_aresout__get__rootid(self._handle)
    
    @rootid.setter
    def rootid(self, rootid):
        _AresMainPy_pkg.f90wrap_aresout__set__rootid(self._handle, rootid)
    
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
        ret.append(',\n    apilat_mat : ')
        ret.append(repr(self.apilat_mat))
        ret.append(',\n    apilat_para : ')
        ret.append(repr(self.apilat_para))
        ret.append(',\n    comm : ')
        ret.append(repr(self.comm))
        ret.append(',\n    myid : ')
        ret.append(repr(self.myid))
        ret.append(',\n    numprocs : ')
        ret.append(repr(self.numprocs))
        ret.append(',\n    rootid : ')
        ret.append(repr(self.rootid))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def init_alloc_arrays(self, nnatom):
    """
    init_alloc_arrays(self, nnatom)
    
    
    Defined at AresMainAPI.fpp lines 24-36
    
    Parameters
    ----------
    dertype : Aresout
    nnatom : int
    
    """
    _AresMainPy_pkg.f90wrap_init_alloc_arrays(dertype=self._handle, nnatom=nnatom)

def assignment(self):
    """
    assignment(self)
    
    
    Defined at AresMainAPI.fpp lines 39-48
    
    Parameters
    ----------
    dertype : Aresout
    
    """
    _AresMainPy_pkg.f90wrap_assignment(dertype=self._handle)

def destroy_alloc_arrays(self):
    """
    destroy_alloc_arrays(self)
    
    
    Defined at AresMainAPI.fpp lines 50-55
    
    Parameters
    ----------
    dertype : Aresout
    
    """
    _AresMainPy_pkg.f90wrap_destroy_alloc_arrays(dertype=self._handle)

def updateions(pos, lattice=None, ikind=None):
    """
    updateions(pos[, lattice, ikind])
    
    
    Defined at AresMainAPI.fpp lines 58-77
    
    Parameters
    ----------
    pos : float array
    lattice : float array
    ikind : int
    
    """
    _AresMainPy_pkg.f90wrap_updateions(pos=pos, lattice=lattice, ikind=ikind)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "aresmainapi".')

for func in _dt_array_initialisers:
    func()

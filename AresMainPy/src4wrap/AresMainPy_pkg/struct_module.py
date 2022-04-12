"""
Module struct_module


Defined at Struct_module.fpp lines 5-106

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("AresMainPy_pkg.struct_type")
class struct_type(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=struct_type)
    
    
    Defined at Struct_module.fpp lines 13-30
    
    """
    def __init__(self, handle=None):
        """
        self = Struct_Type()
        
        
        Defined at Struct_module.fpp lines 13-30
        
        
        Returns
        -------
        this : Struct_Type
        	Object to be constructed
        
        
        Automatically generated constructor for struct_type
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _AresMainPy_pkg.f90wrap_struct_type_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Struct_Type
        
        
        Defined at Struct_module.fpp lines 13-30
        
        Parameters
        ----------
        this : Struct_Type
        	Object to be destructed
        
        
        Automatically generated destructor for struct_type
        """
        if self._alloc:
            _AresMainPy_pkg.f90wrap_struct_type_finalise(this=self._handle)
    
    @property
    def zion(self):
        """
        Element zion ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 14
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_struct_type__array__zion(self._handle)
        if array_handle in self._arrays:
            zion = self._arrays[array_handle]
        else:
            zion = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_struct_type__array__zion)
            self._arrays[array_handle] = zion
        return zion
    
    @zion.setter
    def zion(self, zion):
        self.zion[...] = zion
    
    @property
    def nati(self):
        """
        Element nati ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.fpp line 15
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_struct_type__array__nati(self._handle)
        if array_handle in self._arrays:
            nati = self._arrays[array_handle]
        else:
            nati = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_struct_type__array__nati)
            self._arrays[array_handle] = nati
        return nati
    
    @nati.setter
    def nati(self, nati):
        self.nati[...] = nati
    
    @property
    def eleid(self):
        """
        Element eleid ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.fpp line 16
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_struct_type__array__eleid(self._handle)
        if array_handle in self._arrays:
            eleid = self._arrays[array_handle]
        else:
            eleid = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_struct_type__array__eleid)
            self._arrays[array_handle] = eleid
        return eleid
    
    @eleid.setter
    def eleid(self, eleid):
        self.eleid[...] = eleid
    
    @property
    def pos(self):
        """
        Element pos ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 17
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_struct_type__array__pos(self._handle)
        if array_handle in self._arrays:
            pos = self._arrays[array_handle]
        else:
            pos = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_struct_type__array__pos)
            self._arrays[array_handle] = pos
        return pos
    
    @pos.setter
    def pos(self, pos):
        self.pos[...] = pos
    
    @property
    def poscar(self):
        """
        Element poscar ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 18
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_struct_type__array__poscar(self._handle)
        if array_handle in self._arrays:
            poscar = self._arrays[array_handle]
        else:
            poscar = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_struct_type__array__poscar)
            self._arrays[array_handle] = poscar
        return poscar
    
    @poscar.setter
    def poscar(self, poscar):
        self.poscar[...] = poscar
    
    @property
    def stress(self):
        """
        Element stress ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 19
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_struct_type__array__stress(self._handle)
        if array_handle in self._arrays:
            stress = self._arrays[array_handle]
        else:
            stress = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_struct_type__array__stress)
            self._arrays[array_handle] = stress
        return stress
    
    @stress.setter
    def stress(self, stress):
        self.stress[...] = stress
    
    @property
    def forces(self):
        """
        Element forces ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 20
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_struct_type__array__forces(self._handle)
        if array_handle in self._arrays:
            forces = self._arrays[array_handle]
        else:
            forces = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_struct_type__array__forces)
            self._arrays[array_handle] = forces
        return forces
    
    @forces.setter
    def forces(self, forces):
        self.forces[...] = forces
    
    @property
    def mass(self):
        """
        Element mass ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 21
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_struct_type__array__mass(self._handle)
        if array_handle in self._arrays:
            mass = self._arrays[array_handle]
        else:
            mass = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_struct_type__array__mass)
            self._arrays[array_handle] = mass
        return mass
    
    @mass.setter
    def mass(self, mass):
        self.mass[...] = mass
    
    @property
    def zeta(self):
        """
        Element zeta ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 23
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_struct_type__array__zeta(self._handle)
        if array_handle in self._arrays:
            zeta = self._arrays[array_handle]
        else:
            zeta = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_struct_type__array__zeta)
            self._arrays[array_handle] = zeta
        return zeta
    
    @zeta.setter
    def zeta(self, zeta):
        self.zeta[...] = zeta
    
    @property
    def prinq(self):
        """
        Element prinq ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.fpp line 24
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_struct_type__array__prinq(self._handle)
        if array_handle in self._arrays:
            prinq = self._arrays[array_handle]
        else:
            prinq = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_struct_type__array__prinq)
            self._arrays[array_handle] = prinq
        return prinq
    
    @prinq.setter
    def prinq(self, prinq):
        self.prinq[...] = prinq
    
    @property
    def lmax(self):
        """
        Element lmax ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.fpp line 25
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_struct_type__array__lmax(self._handle)
        if array_handle in self._arrays:
            lmax = self._arrays[array_handle]
        else:
            lmax = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_struct_type__array__lmax)
            self._arrays[array_handle] = lmax
        return lmax
    
    @lmax.setter
    def lmax(self, lmax):
        self.lmax[...] = lmax
    
    @property
    def elements(self):
        """
        Element elements ftype=character(len=3) pytype=str
        
        
        Defined at Struct_module.fpp line 26
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_struct_type__array__elements(self._handle)
        if array_handle in self._arrays:
            elements = self._arrays[array_handle]
        else:
            elements = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_struct_type__array__elements)
            self._arrays[array_handle] = elements
        return elements
    
    @elements.setter
    def elements(self, elements):
        self.elements[...] = elements
    
    @property
    def coeff(self):
        """
        Element coeff ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 28
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_struct_type__array__coeff(self._handle)
        if array_handle in self._arrays:
            coeff = self._arrays[array_handle]
        else:
            coeff = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_struct_type__array__coeff)
            self._arrays[array_handle] = coeff
        return coeff
    
    @coeff.setter
    def coeff(self, coeff):
        self.coeff[...] = coeff
    
    @property
    def occupy(self):
        """
        Element occupy ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 29
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_struct_type__array__occupy(self._handle)
        if array_handle in self._arrays:
            occupy = self._arrays[array_handle]
        else:
            occupy = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_struct_type__array__occupy)
            self._arrays[array_handle] = occupy
        return occupy
    
    @occupy.setter
    def occupy(self, occupy):
        self.occupy[...] = occupy
    
    @property
    def noccupy(self):
        """
        Element noccupy ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.fpp line 30
        
        """
        return _AresMainPy_pkg.f90wrap_struct_type__get__noccupy(self._handle)
    
    @noccupy.setter
    def noccupy(self, noccupy):
        _AresMainPy_pkg.f90wrap_struct_type__set__noccupy(self._handle, noccupy)
    
    def __str__(self):
        ret = ['<struct_type>{\n']
        ret.append('    zion : ')
        ret.append(repr(self.zion))
        ret.append(',\n    nati : ')
        ret.append(repr(self.nati))
        ret.append(',\n    eleid : ')
        ret.append(repr(self.eleid))
        ret.append(',\n    pos : ')
        ret.append(repr(self.pos))
        ret.append(',\n    poscar : ')
        ret.append(repr(self.poscar))
        ret.append(',\n    stress : ')
        ret.append(repr(self.stress))
        ret.append(',\n    forces : ')
        ret.append(repr(self.forces))
        ret.append(',\n    mass : ')
        ret.append(repr(self.mass))
        ret.append(',\n    zeta : ')
        ret.append(repr(self.zeta))
        ret.append(',\n    prinq : ')
        ret.append(repr(self.prinq))
        ret.append(',\n    lmax : ')
        ret.append(repr(self.lmax))
        ret.append(',\n    elements : ')
        ret.append(repr(self.elements))
        ret.append(',\n    coeff : ')
        ret.append(repr(self.coeff))
        ret.append(',\n    occupy : ')
        ret.append(repr(self.occupy))
        ret.append(',\n    noccupy : ')
        ret.append(repr(self.noccupy))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def creat_struct(numtyp, numatom):
    """
    creat_struct(numtyp, numatom)
    
    
    Defined at Struct_module.fpp lines 47-76
    
    Parameters
    ----------
    numtyp : int
    numatom : int
    
    """
    _AresMainPy_pkg.f90wrap_creat_struct(numtyp=numtyp, numatom=numatom)

def destroy_struct():
    """
    destroy_struct()
    
    
    Defined at Struct_module.fpp lines 79-104
    
    
    """
    _AresMainPy_pkg.f90wrap_destroy_struct()

def get_natom():
    """
    Element natom ftype=integer(i4b) pytype=int
    
    
    Defined at Struct_module.fpp line 33
    
    """
    return _AresMainPy_pkg.f90wrap_struct_module__get__natom()

def set_natom(natom):
    _AresMainPy_pkg.f90wrap_struct_module__set__natom(natom)

def get_naty():
    """
    Element naty ftype=integer(i4b) pytype=int
    
    
    Defined at Struct_module.fpp line 34
    
    """
    return _AresMainPy_pkg.f90wrap_struct_module__get__naty()

def set_naty(naty):
    _AresMainPy_pkg.f90wrap_struct_module__set__naty(naty)

def get_ncharge():
    """
    Element ncharge ftype=integer(i4b) pytype=int
    
    
    Defined at Struct_module.fpp line 35
    
    """
    return _AresMainPy_pkg.f90wrap_struct_module__get__ncharge()

def set_ncharge(ncharge):
    _AresMainPy_pkg.f90wrap_struct_module__set__ncharge(ncharge)

def get_charge_ave():
    """
    Element charge_ave ftype=real(dp) pytype=float
    
    
    Defined at Struct_module.fpp line 36
    
    """
    return _AresMainPy_pkg.f90wrap_struct_module__get__charge_ave()

def set_charge_ave(charge_ave):
    _AresMainPy_pkg.f90wrap_struct_module__set__charge_ave(charge_ave)

def get_volume():
    """
    Element volume ftype=real(dp) pytype=float
    
    
    Defined at Struct_module.fpp line 37
    
    """
    return _AresMainPy_pkg.f90wrap_struct_module__get__volume()

def set_volume(volume):
    _AresMainPy_pkg.f90wrap_struct_module__set__volume(volume)

def get_array_lat_mat():
    """
    Element lat_mat ftype=real(dp) pytype=float
    
    
    Defined at Struct_module.fpp line 38
    
    """
    global lat_mat
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_struct_module__array__lat_mat(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        lat_mat = _arrays[array_handle]
    else:
        lat_mat = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_struct_module__array__lat_mat)
        _arrays[array_handle] = lat_mat
    return lat_mat

def set_array_lat_mat(lat_mat):
    lat_mat[...] = lat_mat

def get_array_lat_para():
    """
    Element lat_para ftype=real(dp) pytype=float
    
    
    Defined at Struct_module.fpp line 39
    
    """
    global lat_para
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_struct_module__array__lat_para(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        lat_para = _arrays[array_handle]
    else:
        lat_para = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_struct_module__array__lat_para)
        _arrays[array_handle] = lat_para
    return lat_para

def set_array_lat_para(lat_para):
    lat_para[...] = lat_para

def get_array_recip_lat():
    """
    Element recip_lat ftype=real(dp) pytype=float
    
    
    Defined at Struct_module.fpp line 40
    
    """
    global recip_lat
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_struct_module__array__recip_lat(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        recip_lat = _arrays[array_handle]
    else:
        recip_lat = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_struct_module__array__recip_lat)
        _arrays[array_handle] = recip_lat
    return recip_lat

def set_array_recip_lat(recip_lat):
    recip_lat[...] = recip_lat

def get_array_reclat_para():
    """
    Element reclat_para ftype=real(dp) pytype=float
    
    
    Defined at Struct_module.fpp line 41
    
    """
    global reclat_para
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_struct_module__array__reclat_para(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        reclat_para = _arrays[array_handle]
    else:
        reclat_para = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_struct_module__array__reclat_para)
        _arrays[array_handle] = reclat_para
    return reclat_para

def set_array_reclat_para(reclat_para):
    reclat_para[...] = reclat_para

def get_array_energy():
    """
    Element energy ftype=real(dp) pytype=float
    
    
    Defined at Struct_module.fpp line 42
    
    """
    global energy
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_struct_module__array__energy(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        energy = _arrays[array_handle]
    else:
        energy = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_struct_module__array__energy)
        _arrays[array_handle] = energy
    return energy

def set_array_energy(energy):
    energy[...] = energy


_array_initialisers = [get_array_lat_mat, get_array_lat_para, \
    get_array_recip_lat, get_array_reclat_para, get_array_energy]
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "struct_module".')

for func in _dt_array_initialisers:
    func()

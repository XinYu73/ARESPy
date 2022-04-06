"""
Module struct_module


Defined at Struct_module.fpp lines 11-102

"""
from __future__ import print_function, absolute_import, division
import _arespy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("arespy_pkg.struct_type")
class struct_type(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=struct_type)
    
    
    Defined at Struct_module.fpp lines 19-32
    
    """
    def __init__(self, handle=None):
        """
        self = Struct_Type()
        
        
        Defined at Struct_module.fpp lines 19-32
        
        
        Returns
        -------
        this : Struct_Type
        	Object to be constructed
        
        
        Automatically generated constructor for struct_type
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _arespy_pkg.f90wrap_struct_type_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Struct_Type
        
        
        Defined at Struct_module.fpp lines 19-32
        
        Parameters
        ----------
        this : Struct_Type
        	Object to be destructed
        
        
        Automatically generated destructor for struct_type
        """
        if self._alloc:
            _arespy_pkg.f90wrap_struct_type_finalise(this=self._handle)
    
    @property
    def zion(self):
        """
        Element zion ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 20
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_struct_type__array__zion(self._handle)
        if array_handle in self._arrays:
            zion = self._arrays[array_handle]
        else:
            zion = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_struct_type__array__zion)
            self._arrays[array_handle] = zion
        return zion
    
    @zion.setter
    def zion(self, zion):
        self.zion[...] = zion
    
    @property
    def nati(self):
        """
        Element nati ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.fpp line 21
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_struct_type__array__nati(self._handle)
        if array_handle in self._arrays:
            nati = self._arrays[array_handle]
        else:
            nati = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_struct_type__array__nati)
            self._arrays[array_handle] = nati
        return nati
    
    @nati.setter
    def nati(self, nati):
        self.nati[...] = nati
    
    @property
    def eleid(self):
        """
        Element eleid ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.fpp line 22
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_struct_type__array__eleid(self._handle)
        if array_handle in self._arrays:
            eleid = self._arrays[array_handle]
        else:
            eleid = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_struct_type__array__eleid)
            self._arrays[array_handle] = eleid
        return eleid
    
    @eleid.setter
    def eleid(self, eleid):
        self.eleid[...] = eleid
    
    @property
    def pos(self):
        """
        Element pos ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 23
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_struct_type__array__pos(self._handle)
        if array_handle in self._arrays:
            pos = self._arrays[array_handle]
        else:
            pos = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_struct_type__array__pos)
            self._arrays[array_handle] = pos
        return pos
    
    @pos.setter
    def pos(self, pos):
        self.pos[...] = pos
    
    @property
    def poscar(self):
        """
        Element poscar ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 24
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_struct_type__array__poscar(self._handle)
        if array_handle in self._arrays:
            poscar = self._arrays[array_handle]
        else:
            poscar = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_struct_type__array__poscar)
            self._arrays[array_handle] = poscar
        return poscar
    
    @poscar.setter
    def poscar(self, poscar):
        self.poscar[...] = poscar
    
    @property
    def stress(self):
        """
        Element stress ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 25
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_struct_type__array__stress(self._handle)
        if array_handle in self._arrays:
            stress = self._arrays[array_handle]
        else:
            stress = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_struct_type__array__stress)
            self._arrays[array_handle] = stress
        return stress
    
    @stress.setter
    def stress(self, stress):
        self.stress[...] = stress
    
    @property
    def forces(self):
        """
        Element forces ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 26
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_struct_type__array__forces(self._handle)
        if array_handle in self._arrays:
            forces = self._arrays[array_handle]
        else:
            forces = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_struct_type__array__forces)
            self._arrays[array_handle] = forces
        return forces
    
    @forces.setter
    def forces(self, forces):
        self.forces[...] = forces
    
    @property
    def mass(self):
        """
        Element mass ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 27
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_struct_type__array__mass(self._handle)
        if array_handle in self._arrays:
            mass = self._arrays[array_handle]
        else:
            mass = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_struct_type__array__mass)
            self._arrays[array_handle] = mass
        return mass
    
    @mass.setter
    def mass(self, mass):
        self.mass[...] = mass
    
    @property
    def zeta(self):
        """
        Element zeta ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 29
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_struct_type__array__zeta(self._handle)
        if array_handle in self._arrays:
            zeta = self._arrays[array_handle]
        else:
            zeta = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_struct_type__array__zeta)
            self._arrays[array_handle] = zeta
        return zeta
    
    @zeta.setter
    def zeta(self, zeta):
        self.zeta[...] = zeta
    
    @property
    def prinq(self):
        """
        Element prinq ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.fpp line 30
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_struct_type__array__prinq(self._handle)
        if array_handle in self._arrays:
            prinq = self._arrays[array_handle]
        else:
            prinq = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_struct_type__array__prinq)
            self._arrays[array_handle] = prinq
        return prinq
    
    @prinq.setter
    def prinq(self, prinq):
        self.prinq[...] = prinq
    
    @property
    def lmax(self):
        """
        Element lmax ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.fpp line 31
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_struct_type__array__lmax(self._handle)
        if array_handle in self._arrays:
            lmax = self._arrays[array_handle]
        else:
            lmax = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_struct_type__array__lmax)
            self._arrays[array_handle] = lmax
        return lmax
    
    @lmax.setter
    def lmax(self, lmax):
        self.lmax[...] = lmax
    
    @property
    def elements(self):
        """
        Element elements ftype=character(len=3) pytype=str
        
        
        Defined at Struct_module.fpp line 32
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_struct_type__array__elements(self._handle)
        if array_handle in self._arrays:
            elements = self._arrays[array_handle]
        else:
            elements = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_struct_type__array__elements)
            self._arrays[array_handle] = elements
        return elements
    
    @elements.setter
    def elements(self, elements):
        self.elements[...] = elements
    
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
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def creat_struct(numtyp, numatom):
    """
    creat_struct(numtyp, numatom)
    
    
    Defined at Struct_module.fpp lines 57-81
    
    Parameters
    ----------
    numtyp : int
    numatom : int
    
    """
    _arespy_pkg.f90wrap_creat_struct(numtyp=numtyp, numatom=numatom)

def destroy_struct():
    """
    destroy_struct()
    
    
    Defined at Struct_module.fpp lines 84-100
    
    
    """
    _arespy_pkg.f90wrap_destroy_struct()

def get_natom():
    """
    Element natom ftype=integer(i4b) pytype=int
    
    
    Defined at Struct_module.fpp line 35
    
    """
    return _arespy_pkg.f90wrap_struct_module__get__natom()

def set_natom(natom):
    _arespy_pkg.f90wrap_struct_module__set__natom(natom)

def get_nzion():
    """
    Element nzion ftype=integer(i4b) pytype=int
    
    
    Defined at Struct_module.fpp line 35
    
    """
    return _arespy_pkg.f90wrap_struct_module__get__nzion()

def set_nzion(nzion):
    _arespy_pkg.f90wrap_struct_module__set__nzion(nzion)

def get_naty():
    """
    Element naty ftype=integer(i4b) pytype=int
    
    
    Defined at Struct_module.fpp line 36
    
    """
    return _arespy_pkg.f90wrap_struct_module__get__naty()

def set_naty(naty):
    _arespy_pkg.f90wrap_struct_module__set__naty(naty)

def get_ncharge():
    """
    Element ncharge ftype=real(dp) pytype=float
    
    
    Defined at Struct_module.fpp line 37
    
    """
    return _arespy_pkg.f90wrap_struct_module__get__ncharge()

def set_ncharge(ncharge):
    _arespy_pkg.f90wrap_struct_module__set__ncharge(ncharge)

def get_charge_ave():
    """
    Element charge_ave ftype=real(dp) pytype=float
    
    
    Defined at Struct_module.fpp line 38
    
    """
    return _arespy_pkg.f90wrap_struct_module__get__charge_ave()

def set_charge_ave(charge_ave):
    _arespy_pkg.f90wrap_struct_module__set__charge_ave(charge_ave)

def get_volume():
    """
    Element volume ftype=real(dp) pytype=float
    
    
    Defined at Struct_module.fpp line 39
    
    """
    return _arespy_pkg.f90wrap_struct_module__get__volume()

def set_volume(volume):
    _arespy_pkg.f90wrap_struct_module__set__volume(volume)

def get_volsp():
    """
    Element volsp ftype=real(dp) pytype=float
    
    
    Defined at Struct_module.fpp line 39
    
    """
    return _arespy_pkg.f90wrap_struct_module__get__volsp()

def set_volsp(volsp):
    _arespy_pkg.f90wrap_struct_module__set__volsp(volsp)

def get_array_lat_mat():
    """
    Element lat_mat ftype=real(dp) pytype=float
    
    
    Defined at Struct_module.fpp line 40
    
    """
    global lat_mat
    array_ndim, array_type, array_shape, array_handle = \
        _arespy_pkg.f90wrap_struct_module__array__lat_mat(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        lat_mat = _arrays[array_handle]
    else:
        lat_mat = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_struct_module__array__lat_mat)
        _arrays[array_handle] = lat_mat
    return lat_mat

def set_array_lat_mat(lat_mat):
    lat_mat[...] = lat_mat

def get_array_lat_para():
    """
    Element lat_para ftype=real(dp) pytype=float
    
    
    Defined at Struct_module.fpp line 41
    
    """
    global lat_para
    array_ndim, array_type, array_shape, array_handle = \
        _arespy_pkg.f90wrap_struct_module__array__lat_para(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        lat_para = _arrays[array_handle]
    else:
        lat_para = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_struct_module__array__lat_para)
        _arrays[array_handle] = lat_para
    return lat_para

def set_array_lat_para(lat_para):
    lat_para[...] = lat_para

def get_array_recip_lat():
    """
    Element recip_lat ftype=real(dp) pytype=float
    
    
    Defined at Struct_module.fpp line 42
    
    """
    global recip_lat
    array_ndim, array_type, array_shape, array_handle = \
        _arespy_pkg.f90wrap_struct_module__array__recip_lat(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        recip_lat = _arrays[array_handle]
    else:
        recip_lat = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_struct_module__array__recip_lat)
        _arrays[array_handle] = recip_lat
    return recip_lat

def set_array_recip_lat(recip_lat):
    recip_lat[...] = recip_lat

def get_array_reclat_para():
    """
    Element reclat_para ftype=real(dp) pytype=float
    
    
    Defined at Struct_module.fpp line 43
    
    """
    global reclat_para
    array_ndim, array_type, array_shape, array_handle = \
        _arespy_pkg.f90wrap_struct_module__array__reclat_para(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        reclat_para = _arrays[array_handle]
    else:
        reclat_para = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_struct_module__array__reclat_para)
        _arrays[array_handle] = reclat_para
    return reclat_para

def set_array_reclat_para(reclat_para):
    reclat_para[...] = reclat_para

def get_eionion():
    """
    Element eionion ftype=real(dp) pytype=float
    
    
    Defined at Struct_module.fpp line 44
    
    """
    return _arespy_pkg.f90wrap_struct_module__get__eionion()

def set_eionion(eionion):
    _arespy_pkg.f90wrap_struct_module__set__eionion(eionion)

def get_eshift_ps():
    """
    Element eshift_ps ftype=real(dp) pytype=float
    
    
    Defined at Struct_module.fpp line 45
    
    """
    return _arespy_pkg.f90wrap_struct_module__get__eshift_ps()

def set_eshift_ps(eshift_ps):
    _arespy_pkg.f90wrap_struct_module__set__eshift_ps(eshift_ps)

def get_eshift_tot():
    """
    Element eshift_tot ftype=real(dp) pytype=float
    
    
    Defined at Struct_module.fpp line 45
    
    """
    return _arespy_pkg.f90wrap_struct_module__get__eshift_tot()

def set_eshift_tot(eshift_tot):
    _arespy_pkg.f90wrap_struct_module__set__eshift_tot(eshift_tot)

def get_array_opsym():
    """
    Element opsym ftype=real(dp) pytype=float
    
    
    Defined at Struct_module.fpp line 47
    
    """
    global opsym
    array_ndim, array_type, array_shape, array_handle = \
        _arespy_pkg.f90wrap_struct_module__array__opsym(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        opsym = _arrays[array_handle]
    else:
        opsym = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_struct_module__array__opsym)
        _arrays[array_handle] = opsym
    return opsym

def set_array_opsym(opsym):
    opsym[...] = opsym

def get_array_otrans():
    """
    Element otrans ftype=real(dp) pytype=float
    
    
    Defined at Struct_module.fpp line 48
    
    """
    global otrans
    array_ndim, array_type, array_shape, array_handle = \
        _arespy_pkg.f90wrap_struct_module__array__otrans(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        otrans = _arrays[array_handle]
    else:
        otrans = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_struct_module__array__otrans)
        _arrays[array_handle] = otrans
    return otrans

def set_array_otrans(otrans):
    otrans[...] = otrans

def get_nsym():
    """
    Element nsym ftype=integer(i4b) pytype=int
    
    
    Defined at Struct_module.fpp line 49
    
    """
    return _arespy_pkg.f90wrap_struct_module__get__nsym()

def set_nsym(nsym):
    _arespy_pkg.f90wrap_struct_module__set__nsym(nsym)

def get_num_t():
    """
    Element num_t ftype=integer(i4b) pytype=int
    
    
    Defined at Struct_module.fpp line 50
    
    """
    return _arespy_pkg.f90wrap_struct_module__get__num_t()

def set_num_t(num_t):
    _arespy_pkg.f90wrap_struct_module__set__num_t(num_t)

def get_array_c_i():
    """
    Element c_i ftype=integer(i4b) pytype=int
    
    
    Defined at Struct_module.fpp line 51
    
    """
    global c_i
    array_ndim, array_type, array_shape, array_handle = \
        _arespy_pkg.f90wrap_struct_module__array__c_i(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        c_i = _arrays[array_handle]
    else:
        c_i = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_struct_module__array__c_i)
        _arrays[array_handle] = c_i
    return c_i

def set_array_c_i(c_i):
    c_i[...] = c_i

def get_array_odet():
    """
    Element odet ftype=integer(i4b) pytype=int
    
    
    Defined at Struct_module.fpp line 52
    
    """
    global odet
    array_ndim, array_type, array_shape, array_handle = \
        _arespy_pkg.f90wrap_struct_module__array__odet(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        odet = _arrays[array_handle]
    else:
        odet = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_struct_module__array__odet)
        _arrays[array_handle] = odet
    return odet

def set_array_odet(odet):
    odet[...] = odet


_array_initialisers = [get_array_lat_mat, get_array_lat_para, \
    get_array_recip_lat, get_array_reclat_para, get_array_opsym, \
    get_array_otrans, get_array_c_i, get_array_odet]
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "struct_module".')

for func in _dt_array_initialisers:
    func()

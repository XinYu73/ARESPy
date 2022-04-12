"""
Module read_module


Defined at Read_module.fpp lines 5-2613

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("AresMainPy_pkg.attribute")
class attribute(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=attribute)
    
    
    Defined at Read_module.fpp lines 19-20
    
    """
    def __init__(self, handle=None):
        """
        self = Attribute()
        
        
        Defined at Read_module.fpp lines 19-20
        
        
        Returns
        -------
        this : Attribute
        	Object to be constructed
        
        
        Automatically generated constructor for attribute
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _AresMainPy_pkg.f90wrap_attribute_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Attribute
        
        
        Defined at Read_module.fpp lines 19-20
        
        Parameters
        ----------
        this : Attribute
        	Object to be destructed
        
        
        Automatically generated destructor for attribute
        """
        if self._alloc:
            _AresMainPy_pkg.f90wrap_attribute_finalise(this=self._handle)
    
    @property
    def value(self):
        """
        Element value ftype=character(len=120) pytype=str
        
        
        Defined at Read_module.fpp line 20
        
        """
        return _AresMainPy_pkg.f90wrap_attribute__get__value(self._handle)
    
    @value.setter
    def value(self, value):
        _AresMainPy_pkg.f90wrap_attribute__set__value(self._handle, value)
    
    def __str__(self):
        ret = ['<attribute>{\n']
        ret.append('    value : ')
        ret.append(repr(self.value))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def read_file(infile):
    """
    read_file(infile)
    
    
    Defined at Read_module.fpp lines 26-703
    
    Parameters
    ----------
    infile : str
    
    -----------------------------------------------------------
    """
    _AresMainPy_pkg.f90wrap_read_file(infile=infile)

def read_pos(nty, filename):
    """
    read_pos(nty, filename)
    
    
    Defined at Read_module.fpp lines 706-828
    
    Parameters
    ----------
    nty : int
    filename : str
    
    """
    _AresMainPy_pkg.f90wrap_read_pos(nty=nty, filename=filename)

def resetlattice():
    """
    resetlattice()
    
    
    Defined at Read_module.fpp lines 831-851
    
    
    """
    _AresMainPy_pkg.f90wrap_resetlattice()

def read_pspot_atom(ity, filename):
    """
    ps = read_pspot_atom(ity, filename)
    
    
    Defined at Read_module.fpp lines 854-1172
    
    Parameters
    ----------
    ity : int
    filename : str
    
    Returns
    -------
    ps : Pspot
    
    ======================read comment=========================
    """
    ps = _AresMainPy_pkg.f90wrap_read_pspot_atom(ity=ity, filename=filename)
    ps = f90wrap.runtime.lookup_class("AresMainPy_pkg.pspot").from_handle(ps, \
        alloc=True)
    return ps

def read_pspot(nty, filenames):
    """
    read_pspot(nty, filenames)
    
    
    Defined at Read_module.fpp lines 1175-1456
    
    Parameters
    ----------
    nty : int
    filenames : str array
    
    """
    _AresMainPy_pkg.f90wrap_read_pspot(nty=nty, filenames=filenames)

def read_realpot_atom(ity, filename):
    """
    ps = read_realpot_atom(ity, filename)
    
    
    Defined at Read_module.fpp lines 1459-1804
    
    Parameters
    ----------
    ity : int
    filename : str
    
    Returns
    -------
    ps : Pspot
    
    ======================read comment=========================
    """
    ps = _AresMainPy_pkg.f90wrap_read_realpot_atom(ity=ity, filename=filename)
    ps = f90wrap.runtime.lookup_class("AresMainPy_pkg.pspot").from_handle(ps, \
        alloc=True)
    return ps

def exist_in(string1, string2):
    """
    exist_in = exist_in(string1, string2)
    
    
    Defined at Read_module.fpp lines 1808-1831
    
    Parameters
    ----------
    string1 : str
    string2 : str
    
    Returns
    -------
    exist_in : bool
    
    ====== ====== ======
    """
    exist_in = _AresMainPy_pkg.f90wrap_exist_in(string1=string1, string2=string2)
    return exist_in

def exist_ibegin(string1, string2):
    """
    exist_ibegin = exist_ibegin(string1, string2)
    
    
    Defined at Read_module.fpp lines 1833-1856
    
    Parameters
    ----------
    string1 : str
    string2 : str
    
    Returns
    -------
    exist_ibegin : int
    
    ====== ====== ======
    """
    exist_ibegin = _AresMainPy_pkg.f90wrap_exist_ibegin(string1=string1, \
        string2=string2)
    return exist_ibegin

def scan_head(file_unit, title, start_old):
    """
    scan_head(file_unit, title, start_old)
    
    
    Defined at Read_module.fpp lines 1859-1916
    
    Parameters
    ----------
    file_unit : int
    title : str
    start_old : bool
    
    ====== ====== ======
     print*,"read"//title
    """
    _AresMainPy_pkg.f90wrap_scan_head(file_unit=file_unit, title=title, \
        start_old=start_old)

def scan_tail(file_unit, title):
    """
    scan_tail(file_unit, title)
    
    
    Defined at Read_module.fpp lines 1919-1942
    
    Parameters
    ----------
    file_unit : int
    title : str
    
    ====== ====== ======
    """
    _AresMainPy_pkg.f90wrap_scan_tail(file_unit=file_unit, title=title)

def read_upf(ity, filename, ps):
    """
    read_upf(ity, filename, ps)
    
    
    Defined at Read_module.fpp lines 1945-2224
    
    Parameters
    ----------
    ity : int
    filename : str
    ps : Pspot
    
    -------->Search for Nonlocal potential
    """
    _AresMainPy_pkg.f90wrap_read_upf(ity=ity, filename=filename, ps=ps._handle)

def read_pseudo_header(zion, mesh_size, nproj):
    """
    read_pseudo_header(zion, mesh_size, nproj)
    
    
    Defined at Read_module.fpp lines 2226-2289
    
    Parameters
    ----------
    zion : float
    mesh_size : int
    nproj : int
    
    """
    _AresMainPy_pkg.f90wrap_read_pseudo_header(zion=zion, mesh_size=mesh_size, \
        nproj=nproj)

def read_pseudo_nonlocal(unit_upf, nl, beta_r, d0, rcut, proj_l):
    """
    read_pseudo_nonlocal(unit_upf, nl, beta_r, d0, rcut, proj_l)
    
    
    Defined at Read_module.fpp lines 2291-2323
    
    Parameters
    ----------
    unit_upf : int
    nl : int
    beta_r : float array
    d0 : float array
    rcut : float
    proj_l : int array
    
    """
    _AresMainPy_pkg.f90wrap_read_pseudo_nonlocal(unit_upf=unit_upf, nl=nl, \
        beta_r=beta_r, d0=d0, rcut=rcut, proj_l=proj_l)

def get_mo_coefficient():
    """
    get_mo_coefficient()
    
    
    Defined at Read_module.fpp lines 2386-2493
    
    
    """
    _AresMainPy_pkg.f90wrap_get_mo_coefficient()

def parse_headline(str, atom_sign, atom_id, atom_orbital):
    """
    parse_headline(str, atom_sign, atom_id, atom_orbital)
    
    
    Defined at Read_module.fpp lines 2496-2547
    
    Parameters
    ----------
    str : str
    atom_sign : str array
    atom_id : int array
    atom_orbital : str array
    
    """
    _AresMainPy_pkg.f90wrap_parse_headline(str=str, atom_sign=atom_sign, \
        atom_id=atom_id, atom_orbital=atom_orbital)

def sign2lm(atom_orbital, l, m):
    """
    sign2lm(atom_orbital, l, m)
    
    
    Defined at Read_module.fpp lines 2550-2572
    
    Parameters
    ----------
    atom_orbital : str array
    l : int array
    m : int array
    
    """
    _AresMainPy_pkg.f90wrap_sign2lm(atom_orbital=atom_orbital, l=l, m=m)

def destroy_moinit():
    """
    destroy_moinit()
    
    
    Defined at Read_module.fpp lines 2575-2580
    
    
    """
    _AresMainPy_pkg.f90wrap_destroy_moinit()

def out_concar(filename):
    """
    out_concar(filename)
    
    
    Defined at Read_module.fpp lines 2583-2601
    
    Parameters
    ----------
    filename : str
    
    """
    _AresMainPy_pkg.f90wrap_out_concar(filename=filename)

def read_chgcar():
    """
    read_chgcar()
    
    
    Defined at Read_module.fpp lines 2603-2612
    
    
    """
    _AresMainPy_pkg.f90wrap_read_chgcar()

def _get_value_int(char_in, char_find, variable, find_flag):
    """
    _get_value_int(char_in, char_find, variable, find_flag)
    
    
    Defined at Read_module.fpp lines 2326-2353
    
    Parameters
    ----------
    char_in : str
    char_find : str
    variable : int
    find_flag : bool
    
    """
    _AresMainPy_pkg.f90wrap_get_value_int(char_in=char_in, char_find=char_find, \
        variable=variable, find_flag=find_flag)

def _get_value_real(char_in, char_find, variable, find_flag):
    """
    _get_value_real(char_in, char_find, variable, find_flag)
    
    
    Defined at Read_module.fpp lines 2355-2382
    
    Parameters
    ----------
    char_in : str
    char_find : str
    variable : float
    find_flag : bool
    
    """
    _AresMainPy_pkg.f90wrap_get_value_real(char_in=char_in, char_find=char_find, \
        variable=variable, find_flag=find_flag)

def get_value(*args, **kwargs):
    """
    get_value(*args, **kwargs)
    
    
    Defined at Read_module.fpp lines 16-17
    
    Overloaded interface containing the following procedures:
      _get_value_int
      _get_value_real
    
    """
    for proc in [_get_value_int, _get_value_real]:
        try:
            return proc(*args, **kwargs)
        except TypeError:
            continue


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "read_module".')

for func in _dt_array_initialisers:
    func()

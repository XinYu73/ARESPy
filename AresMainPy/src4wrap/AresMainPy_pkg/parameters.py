"""
Module parameters


Defined at Parameters.fpp lines 5-354

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("AresMainPy_pkg.IO_INFO_TYPE")
class IO_INFO_TYPE(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=io_info_type)
    
    
    Defined at Parameters.fpp lines 125-133
    
    """
    def __init__(self, handle=None):
        """
        self = Io_Info_Type()
        
        
        Defined at Parameters.fpp lines 125-133
        
        
        Returns
        -------
        this : Io_Info_Type
        	Object to be constructed
        
        
        Automatically generated constructor for io_info_type
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _AresMainPy_pkg.f90wrap_io_info_type_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Io_Info_Type
        
        
        Defined at Parameters.fpp lines 125-133
        
        Parameters
        ----------
        this : Io_Info_Type
        	Object to be destructed
        
        
        Automatically generated destructor for io_info_type
        """
        if self._alloc:
            _AresMainPy_pkg.f90wrap_io_info_type_finalise(this=self._handle)
    
    @property
    def iu0(self):
        """
        Element iu0 ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 126
        
        """
        return _AresMainPy_pkg.f90wrap_io_info_type__get__iu0(self._handle)
    
    @iu0.setter
    def iu0(self, iu0):
        _AresMainPy_pkg.f90wrap_io_info_type__set__iu0(self._handle, iu0)
    
    @property
    def iu1(self):
        """
        Element iu1 ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 127
        
        """
        return _AresMainPy_pkg.f90wrap_io_info_type__get__iu1(self._handle)
    
    @iu1.setter
    def iu1(self, iu1):
        _AresMainPy_pkg.f90wrap_io_info_type__set__iu1(self._handle, iu1)
    
    @property
    def iu2(self):
        """
        Element iu2 ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 128
        
        """
        return _AresMainPy_pkg.f90wrap_io_info_type__get__iu2(self._handle)
    
    @iu2.setter
    def iu2(self, iu2):
        _AresMainPy_pkg.f90wrap_io_info_type__set__iu2(self._handle, iu2)
    
    @property
    def iu3(self):
        """
        Element iu3 ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 129
        
        """
        return _AresMainPy_pkg.f90wrap_io_info_type__get__iu3(self._handle)
    
    @iu3.setter
    def iu3(self, iu3):
        _AresMainPy_pkg.f90wrap_io_info_type__set__iu3(self._handle, iu3)
    
    @property
    def iu4(self):
        """
        Element iu4 ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 130
        
        """
        return _AresMainPy_pkg.f90wrap_io_info_type__get__iu4(self._handle)
    
    @iu4.setter
    def iu4(self, iu4):
        _AresMainPy_pkg.f90wrap_io_info_type__set__iu4(self._handle, iu4)
    
    @property
    def iu5(self):
        """
        Element iu5 ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 131
        
        """
        return _AresMainPy_pkg.f90wrap_io_info_type__get__iu5(self._handle)
    
    @iu5.setter
    def iu5(self, iu5):
        _AresMainPy_pkg.f90wrap_io_info_type__set__iu5(self._handle, iu5)
    
    @property
    def iu6(self):
        """
        Element iu6 ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 132
        
        """
        return _AresMainPy_pkg.f90wrap_io_info_type__get__iu6(self._handle)
    
    @iu6.setter
    def iu6(self, iu6):
        _AresMainPy_pkg.f90wrap_io_info_type__set__iu6(self._handle, iu6)
    
    @property
    def iu8(self):
        """
        Element iu8 ftype=integer(i4b) pytype=int
        
        
        Defined at Parameters.fpp line 133
        
        """
        return _AresMainPy_pkg.f90wrap_io_info_type__get__iu8(self._handle)
    
    @iu8.setter
    def iu8(self, iu8):
        _AresMainPy_pkg.f90wrap_io_info_type__set__iu8(self._handle, iu8)
    
    def __str__(self):
        ret = ['<io_info_type>{\n']
        ret.append('    iu0 : ')
        ret.append(repr(self.iu0))
        ret.append(',\n    iu1 : ')
        ret.append(repr(self.iu1))
        ret.append(',\n    iu2 : ')
        ret.append(repr(self.iu2))
        ret.append(',\n    iu3 : ')
        ret.append(repr(self.iu3))
        ret.append(',\n    iu4 : ')
        ret.append(repr(self.iu4))
        ret.append(',\n    iu5 : ')
        ret.append(repr(self.iu5))
        ret.append(',\n    iu6 : ')
        ret.append(repr(self.iu6))
        ret.append(',\n    iu8 : ')
        ret.append(repr(self.iu8))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def generate_infile(infile_name):
    """
    generate_infile(infile_name)
    
    
    Defined at Parameters.fpp lines 216-354
    
    Parameters
    ----------
    infile_name : str
    
    """
    _AresMainPy_pkg.f90wrap_generate_infile(infile_name=infile_name)

def get_ikedf():
    """
    Element ikedf ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 10
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__ikedf()

def set_ikedf(ikedf):
    _AresMainPy_pkg.f90wrap_parameters__set__ikedf(ikedf)

def get_ixcdf():
    """
    Element ixcdf ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 11
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__ixcdf()

def set_ixcdf(ixcdf):
    _AresMainPy_pkg.f90wrap_parameters__set__ixcdf(ixcdf)

def get_finite_order():
    """
    Element finite_order ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 12
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__finite_order()

def set_finite_order(finite_order):
    _AresMainPy_pkg.f90wrap_parameters__set__finite_order(finite_order)

def get_ntype():
    """
    Element ntype ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 13
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__ntype()

def set_ntype(ntype):
    _AresMainPy_pkg.f90wrap_parameters__set__ntype(ntype)

def get_naddstates():
    """
    Element naddstates ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 22
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__naddstates()

def set_naddstates(naddstates):
    _AresMainPy_pkg.f90wrap_parameters__set__naddstates(naddstates)

def get_istart():
    """
    Element istart ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 22
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__istart()

def set_istart(istart):
    _AresMainPy_pkg.f90wrap_parameters__set__istart(istart)

def get_idiag():
    """
    Element idiag ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 22
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__idiag()

def set_idiag(idiag):
    _AresMainPy_pkg.f90wrap_parameters__set__idiag(idiag)

def get_chem():
    """
    Element chem ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 22
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__chem()

def set_chem(chem):
    _AresMainPy_pkg.f90wrap_parameters__set__chem(chem)

def get_chem0():
    """
    Element chem0 ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 22
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__chem0()

def set_chem0(chem0):
    _AresMainPy_pkg.f90wrap_parameters__set__chem0(chem0)

def get_nstates():
    """
    Element nstates ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 22
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__nstates()

def set_nstates(nstates):
    _AresMainPy_pkg.f90wrap_parameters__set__nstates(nstates)

def get_nstates_global():
    """
    Element nstates_global ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 22
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__nstates_global()

def set_nstates_global(nstates_global):
    _AresMainPy_pkg.f90wrap_parameters__set__nstates_global(nstates_global)

def get_nprr():
    """
    Element nprr ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 22
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__nprr()

def set_nprr(nprr):
    _AresMainPy_pkg.f90wrap_parameters__set__nprr(nprr)

def get_nssp():
    """
    Element nssp ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 22
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__nssp()

def set_nssp(nssp):
    _AresMainPy_pkg.f90wrap_parameters__set__nssp(nssp)

def get_kspacing():
    """
    Element kspacing ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 24
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__kspacing()

def set_kspacing(kspacing):
    _AresMainPy_pkg.f90wrap_parameters__set__kspacing(kspacing)

def get_array_kshift():
    """
    Element kshift ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 25
    
    """
    global kshift
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_parameters__array__kshift(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        kshift = _arrays[array_handle]
    else:
        kshift = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_parameters__array__kshift)
        _arrays[array_handle] = kshift
    return kshift

def set_array_kshift(kshift):
    kshift[...] = kshift

def get_init_gap():
    """
    Element init_gap ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 27
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__init_gap()

def set_init_gap(init_gap):
    _AresMainPy_pkg.f90wrap_parameters__set__init_gap(init_gap)

def get_ecut():
    """
    Element ecut ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 27
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__ecut()

def set_ecut(ecut):
    _AresMainPy_pkg.f90wrap_parameters__set__ecut(ecut)

def get_lkp():
    """
    Element lkp ftype=logical pytype=bool
    
    
    Defined at Parameters.fpp line 29
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__lkp()

def set_lkp(lkp):
    _AresMainPy_pkg.f90wrap_parameters__set__lkp(lkp)

def get_pp_identifer():
    """
    Element pp_identifer ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 30
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__pp_identifer()

def set_pp_identifer(pp_identifer):
    _AresMainPy_pkg.f90wrap_parameters__set__pp_identifer(pp_identifer)

def get_lpbc():
    """
    Element lpbc ftype=logical pytype=bool
    
    
    Defined at Parameters.fpp line 33
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__lpbc()

def set_lpbc(lpbc):
    _AresMainPy_pkg.f90wrap_parameters__set__lpbc(lpbc)

def get_calforce():
    """
    Element calforce ftype=logical pytype=bool
    
    
    Defined at Parameters.fpp line 34
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__calforce()

def set_calforce(calforce):
    _AresMainPy_pkg.f90wrap_parameters__set__calforce(calforce)

def get_hartree_method():
    """
    Element hartree_method ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 41
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__hartree_method()

def set_hartree_method(hartree_method):
    _AresMainPy_pkg.f90wrap_parameters__set__hartree_method(hartree_method)

def get_lcell():
    """
    Element lcell ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 42
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__lcell()

def set_lcell(lcell):
    _AresMainPy_pkg.f90wrap_parameters__set__lcell(lcell)

def get_lcellbohr():
    """
    Element lcellbohr ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 42
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__lcellbohr()

def set_lcellbohr(lcellbohr):
    _AresMainPy_pkg.f90wrap_parameters__set__lcellbohr(lcellbohr)

def get_isormax():
    """
    Element isormax ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 43
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__isormax()

def set_isormax(isormax):
    _AresMainPy_pkg.f90wrap_parameters__set__isormax(isormax)

def get_isormaxbohr():
    """
    Element isormaxbohr ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 43
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__isormaxbohr()

def set_isormaxbohr(isormaxbohr):
    _AresMainPy_pkg.f90wrap_parameters__set__isormaxbohr(isormaxbohr)

def get_radiusmax():
    """
    Element radiusmax ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 44
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__radiusmax()

def set_radiusmax(radiusmax):
    _AresMainPy_pkg.f90wrap_parameters__set__radiusmax(radiusmax)

def get_nvc():
    """
    Element nvc ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 45
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__nvc()

def set_nvc(nvc):
    _AresMainPy_pkg.f90wrap_parameters__set__nvc(nvc)

def get_isonorder():
    """
    Element isonorder ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 45
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__isonorder()

def set_isonorder(isonorder):
    _AresMainPy_pkg.f90wrap_parameters__set__isonorder(isonorder)

def get_tolcg():
    """
    Element tolcg ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 46
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__tolcg()

def set_tolcg(tolcg):
    _AresMainPy_pkg.f90wrap_parameters__set__tolcg(tolcg)

def get_nfcd():
    """
    Element nfcd ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 47
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__nfcd()

def set_nfcd(nfcd):
    _AresMainPy_pkg.f90wrap_parameters__set__nfcd(nfcd)

def get_isolmax():
    """
    Element isolmax ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 48
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__isolmax()

def set_isolmax(isolmax):
    _AresMainPy_pkg.f90wrap_parameters__set__isolmax(isolmax)

def get_iprec_fmm():
    """
    Element iprec_fmm ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 49
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__iprec_fmm()

def set_iprec_fmm(iprec_fmm):
    _AresMainPy_pkg.f90wrap_parameters__set__iprec_fmm(iprec_fmm)

def get_addcharge():
    """
    Element addcharge ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 50
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__addcharge()

def set_addcharge(addcharge):
    _AresMainPy_pkg.f90wrap_parameters__set__addcharge(addcharge)

def get_cell_shape():
    """
    Element cell_shape ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 51
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__cell_shape()

def set_cell_shape(cell_shape):
    _AresMainPy_pkg.f90wrap_parameters__set__cell_shape(cell_shape)

def get_cell_thick():
    """
    Element cell_thick ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 52
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__cell_thick()

def set_cell_thick(cell_thick):
    _AresMainPy_pkg.f90wrap_parameters__set__cell_thick(cell_thick)

def get_lpbc2iso():
    """
    Element lpbc2iso ftype=logical pytype=bool
    
    
    Defined at Parameters.fpp line 53
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__lpbc2iso()

def set_lpbc2iso(lpbc2iso):
    _AresMainPy_pkg.f90wrap_parameters__set__lpbc2iso(lpbc2iso)

def get_lradius_auto():
    """
    Element lradius_auto ftype=logical pytype=bool
    
    
    Defined at Parameters.fpp line 54
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__lradius_auto()

def set_lradius_auto(lradius_auto):
    _AresMainPy_pkg.f90wrap_parameters__set__lradius_auto(lradius_auto)

def get_block_mbnb():
    """
    Element block_mbnb ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 56
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__block_mbnb()

def set_block_mbnb(block_mbnb):
    _AresMainPy_pkg.f90wrap_parameters__set__block_mbnb(block_mbnb)

def get_debug_out():
    """
    Element debug_out ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 57
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__debug_out()

def set_debug_out(debug_out):
    _AresMainPy_pkg.f90wrap_parameters__set__debug_out(debug_out)

def get_wexict():
    """
    Element wexict ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 59
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__wexict()

def set_wexict(wexict):
    _AresMainPy_pkg.f90wrap_parameters__set__wexict(wexict)

def get_lbvk():
    """
    Element lbvk ftype=logical pytype=bool
    
    
    Defined at Parameters.fpp line 66
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__lbvk()

def set_lbvk(lbvk):
    _AresMainPy_pkg.f90wrap_parameters__set__lbvk(lbvk)

def get_lfirst():
    """
    Element lfirst ftype=logical pytype=bool
    
    
    Defined at Parameters.fpp line 66
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__lfirst()

def set_lfirst(lfirst):
    _AresMainPy_pkg.f90wrap_parameters__set__lfirst(lfirst)

def get_linrho():
    """
    Element linrho ftype=logical pytype=bool
    
    
    Defined at Parameters.fpp line 66
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__linrho()

def set_linrho(linrho):
    _AresMainPy_pkg.f90wrap_parameters__set__linrho(linrho)

def get_lradrho():
    """
    Element lradrho ftype=logical pytype=bool
    
    
    Defined at Parameters.fpp line 66
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__lradrho()

def set_lradrho(lradrho):
    _AresMainPy_pkg.f90wrap_parameters__set__lradrho(lradrho)

def get_lrrorthnorm():
    """
    Element lrrorthnorm ftype=logical pytype=bool
    
    
    Defined at Parameters.fpp line 66
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__lrrorthnorm()

def set_lrrorthnorm(lrrorthnorm):
    _AresMainPy_pkg.f90wrap_parameters__set__lrrorthnorm(lrrorthnorm)

def get_lrandom():
    """
    Element lrandom ftype=logical pytype=bool
    
    
    Defined at Parameters.fpp line 66
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__lrandom()

def set_lrandom(lrandom):
    _AresMainPy_pkg.f90wrap_parameters__set__lrandom(lrandom)

def get_system_name():
    """
    Element system_name ftype=character(30) pytype=str
    
    
    Defined at Parameters.fpp line 68
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__system_name()

def set_system_name(system_name):
    _AresMainPy_pkg.f90wrap_parameters__set__system_name(system_name)

def get_cellfile_name():
    """
    Element cellfile_name ftype=character(30) pytype=str
    
    
    Defined at Parameters.fpp line 69
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__cellfile_name()

def set_cellfile_name(cellfile_name):
    _AresMainPy_pkg.f90wrap_parameters__set__cellfile_name(cellfile_name)

def get_outfile():
    """
    Element outfile ftype=character(30) pytype=str
    
    
    Defined at Parameters.fpp line 70
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__outfile()

def set_outfile(outfile):
    _AresMainPy_pkg.f90wrap_parameters__set__outfile(outfile)

def get_array_ppfile_name():
    """
    Element ppfile_name ftype=character(30) pytype=str
    
    
    Defined at Parameters.fpp line 71
    
    """
    global ppfile_name
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_parameters__array__ppfile_name(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        ppfile_name = _arrays[array_handle]
    else:
        ppfile_name = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_parameters__array__ppfile_name)
        _arrays[array_handle] = ppfile_name
    return ppfile_name

def set_array_ppfile_name(ppfile_name):
    ppfile_name[...] = ppfile_name

def get_array_elements():
    """
    Element elements ftype=character(30) pytype=str
    
    
    Defined at Parameters.fpp line 72
    
    """
    global elements
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_parameters__array__elements(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        elements = _arrays[array_handle]
    else:
        elements = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_parameters__array__elements)
        _arrays[array_handle] = elements
    return elements

def set_array_elements(elements):
    elements[...] = elements

def get_nspin():
    """
    Element nspin ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 74
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__nspin()

def set_nspin(nspin):
    _AresMainPy_pkg.f90wrap_parameters__set__nspin(nspin)

def get_nsmear():
    """
    Element nsmear ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 76
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__nsmear()

def set_nsmear(nsmear):
    _AresMainPy_pkg.f90wrap_parameters__set__nsmear(nsmear)

def get_wsmear():
    """
    Element wsmear ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 77
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__wsmear()

def set_wsmear(wsmear):
    _AresMainPy_pkg.f90wrap_parameters__set__wsmear(wsmear)

def get_imixer():
    """
    Element imixer ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 80
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__imixer()

def set_imixer(imixer):
    _AresMainPy_pkg.f90wrap_parameters__set__imixer(imixer)

def get_nmiter():
    """
    Element nmiter ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 80
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__nmiter()

def set_nmiter(nmiter):
    _AresMainPy_pkg.f90wrap_parameters__set__nmiter(nmiter)

def get_nsmix():
    """
    Element nsmix ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 83
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__nsmix()

def set_nsmix(nsmix):
    _AresMainPy_pkg.f90wrap_parameters__set__nsmix(nsmix)

def get_nhmix():
    """
    Element nhmix ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 83
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__nhmix()

def set_nhmix(nhmix):
    _AresMainPy_pkg.f90wrap_parameters__set__nhmix(nhmix)

def get_nhmin():
    """
    Element nhmin ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 83
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__nhmin()

def set_nhmin(nhmin):
    _AresMainPy_pkg.f90wrap_parameters__set__nhmin(nhmin)

def get_malpha():
    """
    Element malpha ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 87
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__malpha()

def set_malpha(malpha):
    _AresMainPy_pkg.f90wrap_parameters__set__malpha(malpha)

def get_mbeta():
    """
    Element mbeta ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 87
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__mbeta()

def set_mbeta(mbeta):
    _AresMainPy_pkg.f90wrap_parameters__set__mbeta(mbeta)

def get_amix():
    """
    Element amix ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 87
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__amix()

def set_amix(amix):
    _AresMainPy_pkg.f90wrap_parameters__set__amix(amix)

def get_bmix():
    """
    Element bmix ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 87
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__bmix()

def set_bmix(bmix):
    _AresMainPy_pkg.f90wrap_parameters__set__bmix(bmix)

def get_array_resta():
    """
    Element resta ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 88
    
    """
    global resta
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_parameters__array__resta(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        resta = _arrays[array_handle]
    else:
        resta = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_parameters__array__resta)
        _arrays[array_handle] = resta
    return resta

def set_array_resta(resta):
    resta[...] = resta

def get_w0am():
    """
    Element w0am ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 90
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__w0am()

def set_w0am(w0am):
    _AresMainPy_pkg.f90wrap_parameters__set__w0am(w0am)

def get_rtol():
    """
    Element rtol ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 92
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__rtol()

def set_rtol(rtol):
    _AresMainPy_pkg.f90wrap_parameters__set__rtol(rtol)

def get_etol():
    """
    Element etol ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 92
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__etol()

def set_etol(etol):
    _AresMainPy_pkg.f90wrap_parameters__set__etol(etol)

def get_lsub():
    """
    Element lsub ftype=logical pytype=bool
    
    
    Defined at Parameters.fpp line 94
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__lsub()

def set_lsub(lsub):
    _AresMainPy_pkg.f90wrap_parameters__set__lsub(lsub)

def get_lfat():
    """
    Element lfat ftype=logical pytype=bool
    
    
    Defined at Parameters.fpp line 95
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__lfat()

def set_lfat(lfat):
    _AresMainPy_pkg.f90wrap_parameters__set__lfat(lfat)

def get_lone():
    """
    Element lone ftype=logical pytype=bool
    
    
    Defined at Parameters.fpp line 96
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__lone()

def set_lone(lone):
    _AresMainPy_pkg.f90wrap_parameters__set__lone(lone)

def get_nsub():
    """
    Element nsub ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 97
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__nsub()

def set_nsub(nsub):
    _AresMainPy_pkg.f90wrap_parameters__set__nsub(nsub)

def get_nfat():
    """
    Element nfat ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 98
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__nfat()

def set_nfat(nfat):
    _AresMainPy_pkg.f90wrap_parameters__set__nfat(nfat)

def get_array_tfvw():
    """
    Element tfvw ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 99
    
    """
    global tfvw
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_parameters__array__tfvw(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        tfvw = _arrays[array_handle]
    else:
        tfvw = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_parameters__array__tfvw)
        _arrays[array_handle] = tfvw
    return tfvw

def set_array_tfvw(tfvw):
    tfvw[...] = tfvw

def get_lofdft():
    """
    Element lofdft ftype=logical pytype=bool
    
    
    Defined at Parameters.fpp line 101
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__lofdft()

def set_lofdft(lofdft):
    _AresMainPy_pkg.f90wrap_parameters__set__lofdft(lofdft)

def get_lband():
    """
    Element lband ftype=logical pytype=bool
    
    
    Defined at Parameters.fpp line 103
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__lband()

def set_lband(lband):
    _AresMainPy_pkg.f90wrap_parameters__set__lband(lband)

def get_iopm():
    """
    Element iopm ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 106
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__iopm()

def set_iopm(iopm):
    _AresMainPy_pkg.f90wrap_parameters__set__iopm(iopm)

def get_igoal():
    """
    Element igoal ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 106
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__igoal()

def set_igoal(igoal):
    _AresMainPy_pkg.f90wrap_parameters__set__igoal(igoal)

def get_press():
    """
    Element press ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 110
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__press()

def set_press(press):
    _AresMainPy_pkg.f90wrap_parameters__set__press(press)

def get_tolf():
    """
    Element tolf ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 110
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__tolf()

def set_tolf(tolf):
    _AresMainPy_pkg.f90wrap_parameters__set__tolf(tolf)

def get_tolp():
    """
    Element tolp ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 110
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__tolp()

def set_tolp(tolp):
    _AresMainPy_pkg.f90wrap_parameters__set__tolp(tolp)

def get_times():
    """
    Element times ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 110
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__times()

def set_times(times):
    _AresMainPy_pkg.f90wrap_parameters__set__times(times)

def get_isp():
    """
    Element isp ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 111
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__isp()

def set_isp(isp):
    _AresMainPy_pkg.f90wrap_parameters__set__isp(isp)

def get_ldg():
    """
    Element ldg ftype=logical pytype=bool
    
    
    Defined at Parameters.fpp line 113
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__ldg()

def set_ldg(ldg):
    _AresMainPy_pkg.f90wrap_parameters__set__ldg(ldg)

def get_ndg():
    """
    Element ndg ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 114
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__ndg()

def set_ndg(ndg):
    _AresMainPy_pkg.f90wrap_parameters__set__ndg(ndg)

def get_n_near():
    """
    Element n_near ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 115
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__n_near()

def set_n_near(n_near):
    _AresMainPy_pkg.f90wrap_parameters__set__n_near(n_near)

def get_inpol():
    """
    Element inpol ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 116
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__inpol()

def set_inpol(inpol):
    _AresMainPy_pkg.f90wrap_parameters__set__inpol(inpol)

def get_idinit():
    """
    Element idinit ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 118
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__idinit()

def set_idinit(idinit):
    _AresMainPy_pkg.f90wrap_parameters__set__idinit(idinit)

def get_mo_file():
    """
    Element mo_file ftype=character(30) pytype=str
    
    
    Defined at Parameters.fpp line 119
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__mo_file()

def set_mo_file(mo_file):
    _AresMainPy_pkg.f90wrap_parameters__set__mo_file(mo_file)

def get_potim():
    """
    Element potim ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 121
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__potim()

def set_potim(potim):
    _AresMainPy_pkg.f90wrap_parameters__set__potim(potim)

def get_ediffg():
    """
    Element ediffg ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 122
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__ediffg()

def set_ediffg(ediffg):
    _AresMainPy_pkg.f90wrap_parameters__set__ediffg(ediffg)

def get_step_fixrho():
    """
    Element step_fixrho ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 123
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__step_fixrho()

def set_step_fixrho(step_fixrho):
    _AresMainPy_pkg.f90wrap_parameters__set__step_fixrho(step_fixrho)

def get_array_fix_xyz():
    """
    Element fix_xyz ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 136
    
    """
    global fix_xyz
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_parameters__array__fix_xyz(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        fix_xyz = _arrays[array_handle]
    else:
        fix_xyz = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_parameters__array__fix_xyz)
        _arrays[array_handle] = fix_xyz
    return fix_xyz

def set_array_fix_xyz(fix_xyz):
    fix_xyz[...] = fix_xyz

def get_pstress():
    """
    Element pstress ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 147
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__pstress()

def set_pstress(pstress):
    _AresMainPy_pkg.f90wrap_parameters__set__pstress(pstress)

def get_ramax():
    """
    Element ramax ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 148
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__ramax()

def set_ramax(ramax):
    _AresMainPy_pkg.f90wrap_parameters__set__ramax(ramax)

def get_grad_order():
    """
    Element grad_order ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 150
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__grad_order()

def set_grad_order(grad_order):
    _AresMainPy_pkg.f90wrap_parameters__set__grad_order(grad_order)

def get_maxnpts():
    """
    Element maxnpts ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 151
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__maxnpts()

def set_maxnpts(maxnpts):
    _AresMainPy_pkg.f90wrap_parameters__set__maxnpts(maxnpts)

def get_nevshift():
    """
    Element nevshift ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 152
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__nevshift()

def set_nevshift(nevshift):
    _AresMainPy_pkg.f90wrap_parameters__set__nevshift(nevshift)

def get_array_gridn():
    """
    Element gridn ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 153
    
    """
    global gridn
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_parameters__array__gridn(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        gridn = _arrays[array_handle]
    else:
        gridn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_parameters__array__gridn)
        _arrays[array_handle] = gridn
    return gridn

def set_array_gridn(gridn):
    gridn[...] = gridn

def get_lgamma():
    """
    Element lgamma ftype=logical pytype=bool
    
    
    Defined at Parameters.fpp line 154
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__lgamma()

def set_lgamma(lgamma):
    _AresMainPy_pkg.f90wrap_parameters__set__lgamma(lgamma)

def get_lcore_val():
    """
    Element lcore_val ftype=logical pytype=bool
    
    
    Defined at Parameters.fpp line 155
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__lcore_val()

def set_lcore_val(lcore_val):
    _AresMainPy_pkg.f90wrap_parameters__set__lcore_val(lcore_val)

def get_array_iounits():
    """
    Element iounits ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 156
    
    """
    global iounits
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_parameters__array__iounits(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        iounits = _arrays[array_handle]
    else:
        iounits = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_parameters__array__iounits)
        _arrays[array_handle] = iounits
    return iounits

def set_array_iounits(iounits):
    iounits[...] = iounits

def get_isym():
    """
    Element isym ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 157
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__isym()

def set_isym(isym):
    _AresMainPy_pkg.f90wrap_parameters__set__isym(isym)

def get_lfinite_full():
    """
    Element lfinite_full ftype=logical pytype=bool
    
    
    Defined at Parameters.fpp line 158
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__lfinite_full()

def set_lfinite_full(lfinite_full):
    _AresMainPy_pkg.f90wrap_parameters__set__lfinite_full(lfinite_full)

def get_lforce():
    """
    Element lforce ftype=logical pytype=bool
    
    
    Defined at Parameters.fpp line 159
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__lforce()

def set_lforce(lforce):
    _AresMainPy_pkg.f90wrap_parameters__set__lforce(lforce)

def get_lstress():
    """
    Element lstress ftype=logical pytype=bool
    
    
    Defined at Parameters.fpp line 160
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__lstress()

def set_lstress(lstress):
    _AresMainPy_pkg.f90wrap_parameters__set__lstress(lstress)

def get_temper_elec():
    """
    Element temper_elec ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 161
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__temper_elec()

def set_temper_elec(temper_elec):
    _AresMainPy_pkg.f90wrap_parameters__set__temper_elec(temper_elec)

def get_lwstyle():
    """
    Element lwstyle ftype=character(len=20) pytype=str
    
    
    Defined at Parameters.fpp line 162
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__lwstyle()

def set_lwstyle(lwstyle):
    _AresMainPy_pkg.f90wrap_parameters__set__lwstyle(lwstyle)

def get_lopt():
    """
    Element lopt ftype=logical pytype=bool
    
    
    Defined at Parameters.fpp line 163
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__lopt()

def set_lopt(lopt):
    _AresMainPy_pkg.f90wrap_parameters__set__lopt(lopt)

def get_ibron():
    """
    Element ibron ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 164
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__ibron()

def set_ibron(ibron):
    _AresMainPy_pkg.f90wrap_parameters__set__ibron(ibron)

def get_maxsave():
    """
    Element maxsave ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 165
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__maxsave()

def set_maxsave(maxsave):
    _AresMainPy_pkg.f90wrap_parameters__set__maxsave(maxsave)

def get_igamma():
    """
    Element igamma ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 166
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__igamma()

def set_igamma(igamma):
    _AresMainPy_pkg.f90wrap_parameters__set__igamma(igamma)

def get_array_kgrid():
    """
    Element kgrid ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 167
    
    """
    global kgrid
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_parameters__array__kgrid(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        kgrid = _arrays[array_handle]
    else:
        kgrid = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_parameters__array__kgrid)
        _arrays[array_handle] = kgrid
    return kgrid

def set_array_kgrid(kgrid):
    kgrid[...] = kgrid

def get_lmd():
    """
    Element lmd ftype=logical pytype=bool
    
    
    Defined at Parameters.fpp line 169
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__lmd()

def set_lmd(lmd):
    _AresMainPy_pkg.f90wrap_parameters__set__lmd(lmd)

def get_nhis():
    """
    Element nhis ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 171
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__nhis()

def set_nhis(nhis):
    _AresMainPy_pkg.f90wrap_parameters__set__nhis(nhis)

def get_rdfr():
    """
    Element rdfr ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 172
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__rdfr()

def set_rdfr(rdfr):
    _AresMainPy_pkg.f90wrap_parameters__set__rdfr(rdfr)

def get_thermostat():
    """
    Element thermostat ftype=character(len=30) pytype=str
    
    
    Defined at Parameters.fpp line 173
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__thermostat()

def set_thermostat(thermostat):
    _AresMainPy_pkg.f90wrap_parameters__set__thermostat(thermostat)

def get_integrator():
    """
    Element integrator ftype=character(len=30) pytype=str
    
    
    Defined at Parameters.fpp line 174
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__integrator()

def set_integrator(integrator):
    _AresMainPy_pkg.f90wrap_parameters__set__integrator(integrator)

def get_sfreq():
    """
    Element sfreq ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 176
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__sfreq()

def set_sfreq(sfreq):
    _AresMainPy_pkg.f90wrap_parameters__set__sfreq(sfreq)

def get_cfreq():
    """
    Element cfreq ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 177
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__cfreq()

def set_cfreq(cfreq):
    _AresMainPy_pkg.f90wrap_parameters__set__cfreq(cfreq)

def get_temperature():
    """
    Element temperature ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 178
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__temperature()

def set_temperature(temperature):
    _AresMainPy_pkg.f90wrap_parameters__set__temperature(temperature)

def get_mste():
    """
    Element mste ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 179
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__mste()

def set_mste(mste):
    _AresMainPy_pkg.f90wrap_parameters__set__mste(mste)

def get_delt():
    """
    Element delt ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 180
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__delt()

def set_delt(delt):
    _AresMainPy_pkg.f90wrap_parameters__set__delt(delt)

def get_relaxt():
    """
    Element relaxt ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 181
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__relaxt()

def set_relaxt(relaxt):
    _AresMainPy_pkg.f90wrap_parameters__set__relaxt(relaxt)

def get_ensemble():
    """
    Element ensemble ftype=character(10) pytype=str
    
    
    Defined at Parameters.fpp line 182
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__ensemble()

def set_ensemble(ensemble):
    _AresMainPy_pkg.f90wrap_parameters__set__ensemble(ensemble)

def get_dof():
    """
    Element dof ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 183
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__dof()

def set_dof(dof):
    _AresMainPy_pkg.f90wrap_parameters__set__dof(dof)

def get_pext():
    """
    Element pext ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 184
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__pext()

def set_pext(pext):
    _AresMainPy_pkg.f90wrap_parameters__set__pext(pext)

def get_iresmd():
    """
    Element iresmd ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 185
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__iresmd()

def set_iresmd(iresmd):
    _AresMainPy_pkg.f90wrap_parameters__set__iresmd(iresmd)

def get_nresn():
    """
    Element nresn ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 188
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__nresn()

def set_nresn(nresn):
    _AresMainPy_pkg.f90wrap_parameters__set__nresn(nresn)

def get_nyosh():
    """
    Element nyosh ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 189
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__nyosh()

def set_nyosh(nyosh):
    _AresMainPy_pkg.f90wrap_parameters__set__nyosh(nyosh)

def get_nnhc():
    """
    Element nnhc ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 190
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__nnhc()

def set_nnhc(nnhc):
    _AresMainPy_pkg.f90wrap_parameters__set__nnhc(nnhc)

def get_array_wdti2():
    """
    Element wdti2 ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 191
    
    """
    global wdti2
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_parameters__array__wdti2(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        wdti2 = _arrays[array_handle]
    else:
        wdti2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_parameters__array__wdti2)
        _arrays[array_handle] = wdti2
    return wdti2

def set_array_wdti2(wdti2):
    wdti2[...] = wdti2

def get_array_wdti4():
    """
    Element wdti4 ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 191
    
    """
    global wdti4
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_parameters__array__wdti4(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        wdti4 = _arrays[array_handle]
    else:
        wdti4 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_parameters__array__wdti4)
        _arrays[array_handle] = wdti4
    return wdti4

def set_array_wdti4(wdti4):
    wdti4[...] = wdti4

def get_array_wdti8():
    """
    Element wdti8 ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 191
    
    """
    global wdti8
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_parameters__array__wdti8(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        wdti8 = _arrays[array_handle]
    else:
        wdti8 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_parameters__array__wdti8)
        _arrays[array_handle] = wdti8
    return wdti8

def set_array_wdti8(wdti8):
    wdti8[...] = wdti8

def get_array_qmass():
    """
    Element qmass ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 192
    
    """
    global qmass
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_parameters__array__qmass(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        qmass = _arrays[array_handle]
    else:
        qmass = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_parameters__array__qmass)
        _arrays[array_handle] = qmass
    return qmass

def set_array_qmass(qmass):
    qmass[...] = qmass

def get_array_syin_coeff():
    """
    Element syin_coeff ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 193
    
    """
    global syin_coeff
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_parameters__array__syin_coeff(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        syin_coeff = _arrays[array_handle]
    else:
        syin_coeff = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_parameters__array__syin_coeff)
        _arrays[array_handle] = syin_coeff
    return syin_coeff

def set_array_syin_coeff(syin_coeff):
    syin_coeff[...] = syin_coeff

def get_bmass():
    """
    Element bmass ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 194
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__bmass()

def set_bmass(bmass):
    _AresMainPy_pkg.f90wrap_parameters__set__bmass(bmass)

def get_fthermo():
    """
    Element fthermo ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 195
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__fthermo()

def set_fthermo(fthermo):
    _AresMainPy_pkg.f90wrap_parameters__set__fthermo(fthermo)

def get_fbaro():
    """
    Element fbaro ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 196
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__fbaro()

def set_fbaro(fbaro):
    _AresMainPy_pkg.f90wrap_parameters__set__fbaro(fbaro)

def get_fthermown():
    """
    Element fthermown ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 197
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__fthermown()

def set_fthermown(fthermown):
    _AresMainPy_pkg.f90wrap_parameters__set__fthermown(fthermown)

def get_fbarown():
    """
    Element fbarown ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 198
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__fbarown()

def set_fbarown(fbarown):
    _AresMainPy_pkg.f90wrap_parameters__set__fbarown(fbarown)

def get_lbin():
    """
    Element lbin ftype=logical pytype=bool
    
    
    Defined at Parameters.fpp line 200
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__lbin()

def set_lbin(lbin):
    _AresMainPy_pkg.f90wrap_parameters__set__lbin(lbin)

def get_array_p_flag():
    """
    Element p_flag ftype=integer(i4b) pytype=int
    
    
    Defined at Parameters.fpp line 202
    
    """
    global p_flag
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_parameters__array__p_flag(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        p_flag = _arrays[array_handle]
    else:
        p_flag = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_parameters__array__p_flag)
        _arrays[array_handle] = p_flag
    return p_flag

def set_array_p_flag(p_flag):
    p_flag[...] = p_flag

def get_array_p_start():
    """
    Element p_start ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 203
    
    """
    global p_start
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_parameters__array__p_start(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        p_start = _arrays[array_handle]
    else:
        p_start = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_parameters__array__p_start)
        _arrays[array_handle] = p_start
    return p_start

def set_array_p_start(p_start):
    p_start[...] = p_start

def get_array_p_stop():
    """
    Element p_stop ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 203
    
    """
    global p_stop
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_parameters__array__p_stop(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        p_stop = _arrays[array_handle]
    else:
        p_stop = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_parameters__array__p_stop)
        _arrays[array_handle] = p_stop
    return p_stop

def set_array_p_stop(p_stop):
    p_stop[...] = p_stop

def get_array_p_freq():
    """
    Element p_freq ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 203
    
    """
    global p_freq
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_parameters__array__p_freq(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        p_freq = _arrays[array_handle]
    else:
        p_freq = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_parameters__array__p_freq)
        _arrays[array_handle] = p_freq
    return p_freq

def set_array_p_freq(p_freq):
    p_freq[...] = p_freq

def get_array_p_target():
    """
    Element p_target ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 203
    
    """
    global p_target
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_parameters__array__p_target(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        p_target = _arrays[array_handle]
    else:
        p_target = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_parameters__array__p_target)
        _arrays[array_handle] = p_target
    return p_target

def set_array_p_target(p_target):
    p_target[...] = p_target

def get_t_start():
    """
    Element t_start ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 204
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__t_start()

def set_t_start(t_start):
    _AresMainPy_pkg.f90wrap_parameters__set__t_start(t_start)

def get_t_stop():
    """
    Element t_stop ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 204
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__t_stop()

def set_t_stop(t_stop):
    _AresMainPy_pkg.f90wrap_parameters__set__t_stop(t_stop)

def get_t_freq():
    """
    Element t_freq ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 204
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__t_freq()

def set_t_freq(t_freq):
    _AresMainPy_pkg.f90wrap_parameters__set__t_freq(t_freq)

def get_t_target():
    """
    Element t_target ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 204
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__t_target()

def set_t_target(t_target):
    _AresMainPy_pkg.f90wrap_parameters__set__t_target(t_target)

def get_press_control():
    """
    Element press_control ftype=character(len=clen) pytype=str
    
    
    Defined at Parameters.fpp line 205
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__press_control()

def set_press_control(press_control):
    _AresMainPy_pkg.f90wrap_parameters__set__press_control(press_control)

def get_pcontrol():
    """
    Element pcontrol ftype=character(len=clen) pytype=str
    
    
    Defined at Parameters.fpp line 206
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__pcontrol()

def set_pcontrol(pcontrol):
    _AresMainPy_pkg.f90wrap_parameters__set__pcontrol(pcontrol)

def get_pmass():
    """
    Element pmass ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 207
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__pmass()

def set_pmass(pmass):
    _AresMainPy_pkg.f90wrap_parameters__set__pmass(pmass)

def get_pdrag():
    """
    Element pdrag ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 208
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__pdrag()

def set_pdrag(pdrag):
    _AresMainPy_pkg.f90wrap_parameters__set__pdrag(pdrag)

def get_tdrag():
    """
    Element tdrag ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 209
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__tdrag()

def set_tdrag(tdrag):
    _AresMainPy_pkg.f90wrap_parameters__set__tdrag(tdrag)

def get_erate():
    """
    Element erate ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 210
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__erate()

def set_erate(erate):
    _AresMainPy_pkg.f90wrap_parameters__set__erate(erate)

def get_mstrain():
    """
    Element mstrain ftype=real(dp) pytype=float
    
    
    Defined at Parameters.fpp line 211
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__mstrain()

def set_mstrain(mstrain):
    _AresMainPy_pkg.f90wrap_parameters__set__mstrain(mstrain)

def get_sdir():
    """
    Element sdir ftype=character(len=clen) pytype=str
    
    
    Defined at Parameters.fpp line 212
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__sdir()

def set_sdir(sdir):
    _AresMainPy_pkg.f90wrap_parameters__set__sdir(sdir)

def get_lke_pot():
    """
    Element lke_pot ftype=logical pytype=bool
    
    
    Defined at Parameters.fpp line 213
    
    """
    return _AresMainPy_pkg.f90wrap_parameters__get__lke_pot()

def set_lke_pot(lke_pot):
    _AresMainPy_pkg.f90wrap_parameters__set__lke_pot(lke_pot)


_array_initialisers = [get_array_kshift, get_array_ppfile_name, \
    get_array_elements, get_array_resta, get_array_tfvw, get_array_fix_xyz, \
    get_array_gridn, get_array_iounits, get_array_kgrid, get_array_wdti2, \
    get_array_wdti4, get_array_wdti8, get_array_qmass, get_array_syin_coeff, \
    get_array_p_flag, get_array_p_start, get_array_p_stop, get_array_p_freq, \
    get_array_p_target]
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "parameters".')

for func in _dt_array_initialisers:
    func()

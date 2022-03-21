from __future__ import print_function, absolute_import, division
import _bader
import f90wrap.runtime
import logging

class Kind_Mod(f90wrap.runtime.FortranModule):
    """
    Module kind_mod
    
    
    Defined at kind_mod.f90 lines 5-12
    
    """
    @property
    def q1(self):
        """
        Element q1 ftype=integer pytype=int
        
        
        Defined at kind_mod.f90 line 9
        
        """
        return _bader.f90wrap_kind_mod__get__q1()
    
    @property
    def q2(self):
        """
        Element q2 ftype=integer pytype=int
        
        
        Defined at kind_mod.f90 line 10
        
        """
        return _bader.f90wrap_kind_mod__get__q2()
    
    @property
    def pi(self):
        """
        Element pi ftype=real(q2) pytype=float
        
        
        Defined at kind_mod.f90 line 11
        
        """
        return _bader.f90wrap_kind_mod__get__pi()
    
    def __str__(self):
        ret = ['<kind_mod>{\n']
        ret.append('    q1 : ')
        ret.append(repr(self.q1))
        ret.append(',\n    q2 : ')
        ret.append(repr(self.q2))
        ret.append(',\n    pi : ')
        ret.append(repr(self.pi))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

kind_mod = Kind_Mod()

class Matrix_Mod(f90wrap.runtime.FortranModule):
    """
    Module matrix_mod
    
    
    Defined at matrix_mod.f90 lines 5-155
    
    """
    @staticmethod
    def inverse_(a):
        """
        inverse_ = inverse_(a)
        
        
        Defined at matrix_mod.f90 lines 16-26
        
        Parameters
        ----------
        a : float array
        
        Returns
        -------
        inverse_ : float array
        
        """
        inverse_ = _bader.f90wrap_inverse_(a=a)
        return inverse_
    
    @staticmethod
    def adjoint(a):
        """
        adjoint = adjoint(a)
        
        
        Defined at matrix_mod.f90 lines 31-43
        
        Parameters
        ----------
        a : float array
        
        Returns
        -------
        adjoint : float array
        
        """
        adjoint = _bader.f90wrap_adjoint(a=a)
        return adjoint
    
    @staticmethod
    def matrix_volume(h):
        """
        matrix_volume = matrix_volume(h)
        
        
        Defined at matrix_mod.f90 lines 48-55
        
        Parameters
        ----------
        h : float array
        
        Returns
        -------
        matrix_volume : float
        
        """
        matrix_volume = _bader.f90wrap_matrix_volume(h=h)
        return matrix_volume
    
    @staticmethod
    def triple_product(a, b, c):
        """
        triple_product = triple_product(a, b, c)
        
        
        Defined at matrix_mod.f90 lines 60-66
        
        Parameters
        ----------
        a : float array
        b : float array
        c : float array
        
        Returns
        -------
        triple_product : float
        
        """
        triple_product = _bader.f90wrap_triple_product(a=a, b=b, c=c)
        return triple_product
    
    @staticmethod
    def cross_product(a, b):
        """
        cross_product = cross_product(a, b)
        
        
        Defined at matrix_mod.f90 lines 71-77
        
        Parameters
        ----------
        a : float array
        b : float array
        
        Returns
        -------
        cross_product : float array
        
        """
        cross_product = _bader.f90wrap_cross_product(a=a, b=b)
        return cross_product
    
    @staticmethod
    def determinant(a):
        """
        determinant = determinant(a)
        
        
        Defined at matrix_mod.f90 lines 82-91
        
        Parameters
        ----------
        a : float array
        
        Returns
        -------
        determinant : float
        
        """
        determinant = _bader.f90wrap_determinant(a=a)
        return determinant
    
    @staticmethod
    def eigenvectors(yita2, idm, dm, s1, s2, v1, v2, v3):
        """
        eigenvectors(yita2, idm, dm, s1, s2, v1, v2, v3)
        
        
        Defined at matrix_mod.f90 lines 98-114
        
        Parameters
        ----------
        yita2 : float
        idm : float array
        dm : float array
        s1 : float array
        s2 : float array
        v1 : float array
        v2 : float array
        v3 : float array
        
        """
        _bader.f90wrap_eigenvectors(yita2=yita2, idm=idm, dm=dm, s1=s1, s2=s2, v1=v1, \
            v2=v2, v3=v3)
    
    @staticmethod
    def v2area(v1, v2):
        """
        v2area = v2area(v1, v2)
        
        
        Defined at matrix_mod.f90 lines 119-125
        
        Parameters
        ----------
        v1 : float array
        v2 : float array
        
        Returns
        -------
        v2area : float
        
        """
        v2area = _bader.f90wrap_v2area(v1=v1, v2=v2)
        return v2area
    
    @staticmethod
    def normalizelatticevectors(v1):
        """
        normalizelatticevectors = normalizelatticevectors(v1)
        
        
        Defined at matrix_mod.f90 lines 131-154
        
        Parameters
        ----------
        v1 : float array
        
        Returns
        -------
        normalizelatticevectors : float array
        
        """
        normalizelatticevectors = _bader.f90wrap_normalizelatticevectors(v1=v1)
        return normalizelatticevectors
    
    _dt_array_initialisers = []
    

matrix_mod = Matrix_Mod()

class Ions_Mod(f90wrap.runtime.FortranModule):
    """
    Module ions_mod
    
    
    Defined at ions_mod.f90 lines 5-19
    
    """
    @f90wrap.runtime.register_class("bader.ions_obj")
    class ions_obj(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=ions_obj)
        
        
        Defined at ions_mod.f90 lines 8-16
        
        """
        def __init__(self, handle=None):
            """
            self = Ions_Obj()
            
            
            Defined at ions_mod.f90 lines 8-16
            
            
            Returns
            -------
            this : Ions_Obj
            	Object to be constructed
            
            
            Automatically generated constructor for ions_obj
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _bader.f90wrap_ions_obj_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Ions_Obj
            
            
            Defined at ions_mod.f90 lines 8-16
            
            Parameters
            ----------
            this : Ions_Obj
            	Object to be destructed
            
            
            Automatically generated destructor for ions_obj
            """
            if self._alloc:
                _bader.f90wrap_ions_obj_finalise(this=self._handle)
        
        @property
        def r_car(self):
            """
            Element r_car ftype=real(q2) pytype=float
            
            
            Defined at ions_mod.f90 line 9
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bader.f90wrap_ions_obj__array__r_car(self._handle)
            if array_handle in self._arrays:
                r_car = self._arrays[array_handle]
            else:
                r_car = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bader.f90wrap_ions_obj__array__r_car)
                self._arrays[array_handle] = r_car
            return r_car
        
        @r_car.setter
        def r_car(self, r_car):
            self.r_car[...] = r_car
        
        @property
        def r_dir(self):
            """
            Element r_dir ftype=real(q2) pytype=float
            
            
            Defined at ions_mod.f90 line 9
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bader.f90wrap_ions_obj__array__r_dir(self._handle)
            if array_handle in self._arrays:
                r_dir = self._arrays[array_handle]
            else:
                r_dir = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bader.f90wrap_ions_obj__array__r_dir)
                self._arrays[array_handle] = r_dir
            return r_dir
        
        @r_dir.setter
        def r_dir(self, r_dir):
            self.r_dir[...] = r_dir
        
        @property
        def r_lat(self):
            """
            Element r_lat ftype=real(q2) pytype=float
            
            
            Defined at ions_mod.f90 line 9
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bader.f90wrap_ions_obj__array__r_lat(self._handle)
            if array_handle in self._arrays:
                r_lat = self._arrays[array_handle]
            else:
                r_lat = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bader.f90wrap_ions_obj__array__r_lat)
                self._arrays[array_handle] = r_lat
            return r_lat
        
        @r_lat.setter
        def r_lat(self, r_lat):
            self.r_lat[...] = r_lat
        
        @property
        def ion_chg(self):
            """
            Element ion_chg ftype=real(q2) pytype=float
            
            
            Defined at ions_mod.f90 line 10
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bader.f90wrap_ions_obj__array__ion_chg(self._handle)
            if array_handle in self._arrays:
                ion_chg = self._arrays[array_handle]
            else:
                ion_chg = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bader.f90wrap_ions_obj__array__ion_chg)
                self._arrays[array_handle] = ion_chg
            return ion_chg
        
        @ion_chg.setter
        def ion_chg(self, ion_chg):
            self.ion_chg[...] = ion_chg
        
        @property
        def lattice(self):
            """
            Element lattice ftype=real(q2) pytype=float
            
            
            Defined at ions_mod.f90 line 11
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bader.f90wrap_ions_obj__array__lattice(self._handle)
            if array_handle in self._arrays:
                lattice = self._arrays[array_handle]
            else:
                lattice = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bader.f90wrap_ions_obj__array__lattice)
                self._arrays[array_handle] = lattice
            return lattice
        
        @lattice.setter
        def lattice(self, lattice):
            self.lattice[...] = lattice
        
        @property
        def dir2car(self):
            """
            Element dir2car ftype=real(q2) pytype=float
            
            
            Defined at ions_mod.f90 line 11
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bader.f90wrap_ions_obj__array__dir2car(self._handle)
            if array_handle in self._arrays:
                dir2car = self._arrays[array_handle]
            else:
                dir2car = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bader.f90wrap_ions_obj__array__dir2car)
                self._arrays[array_handle] = dir2car
            return dir2car
        
        @dir2car.setter
        def dir2car(self, dir2car):
            self.dir2car[...] = dir2car
        
        @property
        def car2dir(self):
            """
            Element car2dir ftype=real(q2) pytype=float
            
            
            Defined at ions_mod.f90 line 11
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bader.f90wrap_ions_obj__array__car2dir(self._handle)
            if array_handle in self._arrays:
                car2dir = self._arrays[array_handle]
            else:
                car2dir = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bader.f90wrap_ions_obj__array__car2dir)
                self._arrays[array_handle] = car2dir
            return car2dir
        
        @car2dir.setter
        def car2dir(self, car2dir):
            self.car2dir[...] = car2dir
        
        @property
        def num_ion(self):
            """
            Element num_ion ftype=integer pytype=int
            
            
            Defined at ions_mod.f90 line 12
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bader.f90wrap_ions_obj__array__num_ion(self._handle)
            if array_handle in self._arrays:
                num_ion = self._arrays[array_handle]
            else:
                num_ion = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bader.f90wrap_ions_obj__array__num_ion)
                self._arrays[array_handle] = num_ion
            return num_ion
        
        @num_ion.setter
        def num_ion(self, num_ion):
            self.num_ion[...] = num_ion
        
        @property
        def atomic_num(self):
            """
            Element atomic_num ftype=integer pytype=int
            
            
            Defined at ions_mod.f90 line 13
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bader.f90wrap_ions_obj__array__atomic_num(self._handle)
            if array_handle in self._arrays:
                atomic_num = self._arrays[array_handle]
            else:
                atomic_num = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bader.f90wrap_ions_obj__array__atomic_num)
                self._arrays[array_handle] = atomic_num
            return atomic_num
        
        @atomic_num.setter
        def atomic_num(self, atomic_num):
            self.atomic_num[...] = atomic_num
        
        @property
        def niontypes(self):
            """
            Element niontypes ftype=integer  pytype=int
            
            
            Defined at ions_mod.f90 line 15
            
            """
            return _bader.f90wrap_ions_obj__get__niontypes(self._handle)
        
        @niontypes.setter
        def niontypes(self, niontypes):
            _bader.f90wrap_ions_obj__set__niontypes(self._handle, niontypes)
        
        @property
        def nions(self):
            """
            Element nions ftype=integer  pytype=int
            
            
            Defined at ions_mod.f90 line 15
            
            """
            return _bader.f90wrap_ions_obj__get__nions(self._handle)
        
        @nions.setter
        def nions(self, nions):
            _bader.f90wrap_ions_obj__set__nions(self._handle, nions)
        
        @property
        def scalefactor(self):
            """
            Element scalefactor ftype=real(q2) pytype=float
            
            
            Defined at ions_mod.f90 line 16
            
            """
            return _bader.f90wrap_ions_obj__get__scalefactor(self._handle)
        
        @scalefactor.setter
        def scalefactor(self, scalefactor):
            _bader.f90wrap_ions_obj__set__scalefactor(self._handle, scalefactor)
        
        def __str__(self):
            ret = ['<ions_obj>{\n']
            ret.append('    r_car : ')
            ret.append(repr(self.r_car))
            ret.append(',\n    r_dir : ')
            ret.append(repr(self.r_dir))
            ret.append(',\n    r_lat : ')
            ret.append(repr(self.r_lat))
            ret.append(',\n    ion_chg : ')
            ret.append(repr(self.ion_chg))
            ret.append(',\n    lattice : ')
            ret.append(repr(self.lattice))
            ret.append(',\n    dir2car : ')
            ret.append(repr(self.dir2car))
            ret.append(',\n    car2dir : ')
            ret.append(repr(self.car2dir))
            ret.append(',\n    num_ion : ')
            ret.append(repr(self.num_ion))
            ret.append(',\n    atomic_num : ')
            ret.append(repr(self.atomic_num))
            ret.append(',\n    niontypes : ')
            ret.append(repr(self.niontypes))
            ret.append(',\n    nions : ')
            ret.append(repr(self.nions))
            ret.append(',\n    scalefactor : ')
            ret.append(repr(self.scalefactor))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    _dt_array_initialisers = []
    

ions_mod = Ions_Mod()

class Options_Mod(f90wrap.runtime.FortranModule):
    """
    Module options_mod
    
    
    Defined at options_mod.f90 lines 5-632
    
    """
    @f90wrap.runtime.register_class("bader.options_obj")
    class options_obj(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=options_obj)
        
        
        Defined at options_mod.f90 lines 8-33
        
        """
        def __init__(self, handle=None):
            """
            self = Options_Obj()
            
            
            Defined at options_mod.f90 lines 8-33
            
            
            Returns
            -------
            this : Options_Obj
            	Object to be constructed
            
            
            Automatically generated constructor for options_obj
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _bader.f90wrap_options_obj_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Options_Obj
            
            
            Defined at options_mod.f90 lines 8-33
            
            Parameters
            ----------
            this : Options_Obj
            	Object to be destructed
            
            
            Automatically generated destructor for options_obj
            """
            if self._alloc:
                _bader.f90wrap_options_obj_finalise(this=self._handle)
        
        @property
        def chargefile(self):
            """
            Element chargefile ftype=character(len=128) pytype=str
            
            
            Defined at options_mod.f90 line 9
            
            """
            return _bader.f90wrap_options_obj__get__chargefile(self._handle)
        
        @chargefile.setter
        def chargefile(self, chargefile):
            _bader.f90wrap_options_obj__set__chargefile(self._handle, chargefile)
        
        @property
        def refchgfile(self):
            """
            Element refchgfile ftype=character(len=128) pytype=str
            
            
            Defined at options_mod.f90 line 9
            
            """
            return _bader.f90wrap_options_obj__get__refchgfile(self._handle)
        
        @refchgfile.setter
        def refchgfile(self, refchgfile):
            _bader.f90wrap_options_obj__set__refchgfile(self._handle, refchgfile)
        
        @property
        def badertol(self):
            """
            Element badertol ftype=real(q2) pytype=float
            
            
            Defined at options_mod.f90 line 10
            
            """
            return _bader.f90wrap_options_obj__get__badertol(self._handle)
        
        @badertol.setter
        def badertol(self, badertol):
            _bader.f90wrap_options_obj__set__badertol(self._handle, badertol)
        
        @property
        def stepsize(self):
            """
            Element stepsize ftype=real(q2) pytype=float
            
            
            Defined at options_mod.f90 line 10
            
            """
            return _bader.f90wrap_options_obj__get__stepsize(self._handle)
        
        @stepsize.setter
        def stepsize(self, stepsize):
            _bader.f90wrap_options_obj__set__stepsize(self._handle, stepsize)
        
        @property
        def vacval(self):
            """
            Element vacval ftype=real(q2) pytype=float
            
            
            Defined at options_mod.f90 line 10
            
            """
            return _bader.f90wrap_options_obj__get__vacval(self._handle)
        
        @vacval.setter
        def vacval(self, vacval):
            _bader.f90wrap_options_obj__set__vacval(self._handle, vacval)
        
        @property
        def out_opt(self):
            """
            Element out_opt ftype=integer  pytype=int
            
            
            Defined at options_mod.f90 line 11
            
            """
            return _bader.f90wrap_options_obj__get__out_opt(self._handle)
        
        @out_opt.setter
        def out_opt(self, out_opt):
            _bader.f90wrap_options_obj__set__out_opt(self._handle, out_opt)
        
        @property
        def out_auto(self):
            """
            Element out_auto ftype=integer  pytype=int
            
            
            Defined at options_mod.f90 line 11
            
            """
            return _bader.f90wrap_options_obj__get__out_auto(self._handle)
        
        @out_auto.setter
        def out_auto(self, out_auto):
            _bader.f90wrap_options_obj__set__out_auto(self._handle, out_auto)
        
        @property
        def out_cube(self):
            """
            Element out_cube ftype=integer  pytype=int
            
            
            Defined at options_mod.f90 line 11
            
            """
            return _bader.f90wrap_options_obj__get__out_cube(self._handle)
        
        @out_cube.setter
        def out_cube(self, out_cube):
            _bader.f90wrap_options_obj__set__out_cube(self._handle, out_cube)
        
        @property
        def out_chgcar4(self):
            """
            Element out_chgcar4 ftype=integer  pytype=int
            
            
            Defined at options_mod.f90 line 11
            
            """
            return _bader.f90wrap_options_obj__get__out_chgcar4(self._handle)
        
        @out_chgcar4.setter
        def out_chgcar4(self, out_chgcar4):
            _bader.f90wrap_options_obj__set__out_chgcar4(self._handle, out_chgcar4)
        
        @property
        def out_chgcar5(self):
            """
            Element out_chgcar5 ftype=integer  pytype=int
            
            
            Defined at options_mod.f90 line 11
            
            """
            return _bader.f90wrap_options_obj__get__out_chgcar5(self._handle)
        
        @out_chgcar5.setter
        def out_chgcar5(self, out_chgcar5):
            _bader.f90wrap_options_obj__set__out_chgcar5(self._handle, out_chgcar5)
        
        @property
        def in_opt(self):
            """
            Element in_opt ftype=integer  pytype=int
            
            
            Defined at options_mod.f90 line 12
            
            """
            return _bader.f90wrap_options_obj__get__in_opt(self._handle)
        
        @in_opt.setter
        def in_opt(self, in_opt):
            _bader.f90wrap_options_obj__set__in_opt(self._handle, in_opt)
        
        @property
        def in_auto(self):
            """
            Element in_auto ftype=integer  pytype=int
            
            
            Defined at options_mod.f90 line 12
            
            """
            return _bader.f90wrap_options_obj__get__in_auto(self._handle)
        
        @in_auto.setter
        def in_auto(self, in_auto):
            _bader.f90wrap_options_obj__set__in_auto(self._handle, in_auto)
        
        @property
        def in_cube(self):
            """
            Element in_cube ftype=integer  pytype=int
            
            
            Defined at options_mod.f90 line 12
            
            """
            return _bader.f90wrap_options_obj__get__in_cube(self._handle)
        
        @in_cube.setter
        def in_cube(self, in_cube):
            _bader.f90wrap_options_obj__set__in_cube(self._handle, in_cube)
        
        @property
        def in_chgcar(self):
            """
            Element in_chgcar ftype=integer  pytype=int
            
            
            Defined at options_mod.f90 line 12
            
            """
            return _bader.f90wrap_options_obj__get__in_chgcar(self._handle)
        
        @in_chgcar.setter
        def in_chgcar(self, in_chgcar):
            _bader.f90wrap_options_obj__set__in_chgcar(self._handle, in_chgcar)
        
        @property
        def in_chgcar4(self):
            """
            Element in_chgcar4 ftype=integer  pytype=int
            
            
            Defined at options_mod.f90 line 12
            
            """
            return _bader.f90wrap_options_obj__get__in_chgcar4(self._handle)
        
        @in_chgcar4.setter
        def in_chgcar4(self, in_chgcar4):
            _bader.f90wrap_options_obj__set__in_chgcar4(self._handle, in_chgcar4)
        
        @property
        def in_chgcar5(self):
            """
            Element in_chgcar5 ftype=integer  pytype=int
            
            
            Defined at options_mod.f90 line 12
            
            """
            return _bader.f90wrap_options_obj__get__in_chgcar5(self._handle)
        
        @in_chgcar5.setter
        def in_chgcar5(self, in_chgcar5):
            _bader.f90wrap_options_obj__set__in_chgcar5(self._handle, in_chgcar5)
        
        @property
        def bader_opt(self):
            """
            Element bader_opt ftype=integer  pytype=int
            
            
            Defined at options_mod.f90 line 13
            
            """
            return _bader.f90wrap_options_obj__get__bader_opt(self._handle)
        
        @bader_opt.setter
        def bader_opt(self, bader_opt):
            _bader.f90wrap_options_obj__set__bader_opt(self._handle, bader_opt)
        
        @property
        def bader_offgrid(self):
            """
            Element bader_offgrid ftype=integer  pytype=int
            
            
            Defined at options_mod.f90 line 13
            
            """
            return _bader.f90wrap_options_obj__get__bader_offgrid(self._handle)
        
        @bader_offgrid.setter
        def bader_offgrid(self, bader_offgrid):
            _bader.f90wrap_options_obj__set__bader_offgrid(self._handle, bader_offgrid)
        
        @property
        def bader_ongrid(self):
            """
            Element bader_ongrid ftype=integer  pytype=int
            
            
            Defined at options_mod.f90 line 13
            
            """
            return _bader.f90wrap_options_obj__get__bader_ongrid(self._handle)
        
        @bader_ongrid.setter
        def bader_ongrid(self, bader_ongrid):
            _bader.f90wrap_options_obj__set__bader_ongrid(self._handle, bader_ongrid)
        
        @property
        def bader_neargrid(self):
            """
            Element bader_neargrid ftype=integer  pytype=int
            
            
            Defined at options_mod.f90 line 13
            
            """
            return _bader.f90wrap_options_obj__get__bader_neargrid(self._handle)
        
        @bader_neargrid.setter
        def bader_neargrid(self, bader_neargrid):
            _bader.f90wrap_options_obj__set__bader_neargrid(self._handle, bader_neargrid)
        
        @property
        def bader_weight(self):
            """
            Element bader_weight ftype=integer  pytype=int
            
            
            Defined at options_mod.f90 line 13
            
            """
            return _bader.f90wrap_options_obj__get__bader_weight(self._handle)
        
        @bader_weight.setter
        def bader_weight(self, bader_weight):
            _bader.f90wrap_options_obj__set__bader_weight(self._handle, bader_weight)
        
        @property
        def quit_opt(self):
            """
            Element quit_opt ftype=integer  pytype=int
            
            
            Defined at options_mod.f90 line 14
            
            """
            return _bader.f90wrap_options_obj__get__quit_opt(self._handle)
        
        @quit_opt.setter
        def quit_opt(self, quit_opt):
            _bader.f90wrap_options_obj__set__quit_opt(self._handle, quit_opt)
        
        @property
        def quit_max(self):
            """
            Element quit_max ftype=integer  pytype=int
            
            
            Defined at options_mod.f90 line 14
            
            """
            return _bader.f90wrap_options_obj__get__quit_max(self._handle)
        
        @quit_max.setter
        def quit_max(self, quit_max):
            _bader.f90wrap_options_obj__set__quit_max(self._handle, quit_max)
        
        @property
        def quit_known(self):
            """
            Element quit_known ftype=integer  pytype=int
            
            
            Defined at options_mod.f90 line 14
            
            """
            return _bader.f90wrap_options_obj__get__quit_known(self._handle)
        
        @quit_known.setter
        def quit_known(self, quit_known):
            _bader.f90wrap_options_obj__set__quit_known(self._handle, quit_known)
        
        @property
        def refine_edge_itrs(self):
            """
            Element refine_edge_itrs ftype=integer  pytype=int
            
            
            Defined at options_mod.f90 line 15
            
            """
            return _bader.f90wrap_options_obj__get__refine_edge_itrs(self._handle)
        
        @refine_edge_itrs.setter
        def refine_edge_itrs(self, refine_edge_itrs):
            _bader.f90wrap_options_obj__set__refine_edge_itrs(self._handle, \
                refine_edge_itrs)
        
        @property
        def selanum(self):
            """
            Element selanum ftype=integer  pytype=int
            
            
            Defined at options_mod.f90 line 18
            
            """
            return _bader.f90wrap_options_obj__get__selanum(self._handle)
        
        @selanum.setter
        def selanum(self, selanum):
            _bader.f90wrap_options_obj__set__selanum(self._handle, selanum)
        
        @property
        def selbnum(self):
            """
            Element selbnum ftype=integer  pytype=int
            
            
            Defined at options_mod.f90 line 18
            
            """
            return _bader.f90wrap_options_obj__get__selbnum(self._handle)
        
        @selbnum.setter
        def selbnum(self, selbnum):
            _bader.f90wrap_options_obj__set__selbnum(self._handle, selbnum)
        
        @property
        def sumanum(self):
            """
            Element sumanum ftype=integer  pytype=int
            
            
            Defined at options_mod.f90 line 18
            
            """
            return _bader.f90wrap_options_obj__get__sumanum(self._handle)
        
        @sumanum.setter
        def sumanum(self, sumanum):
            _bader.f90wrap_options_obj__set__sumanum(self._handle, sumanum)
        
        @property
        def sumbnum(self):
            """
            Element sumbnum ftype=integer  pytype=int
            
            
            Defined at options_mod.f90 line 18
            
            """
            return _bader.f90wrap_options_obj__get__sumbnum(self._handle)
        
        @sumbnum.setter
        def sumbnum(self, sumbnum):
            _bader.f90wrap_options_obj__set__sumbnum(self._handle, sumbnum)
        
        @property
        def selavol(self):
            """
            Element selavol ftype=integer pytype=int
            
            
            Defined at options_mod.f90 line 19
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bader.f90wrap_options_obj__array__selavol(self._handle)
            if array_handle in self._arrays:
                selavol = self._arrays[array_handle]
            else:
                selavol = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bader.f90wrap_options_obj__array__selavol)
                self._arrays[array_handle] = selavol
            return selavol
        
        @selavol.setter
        def selavol(self, selavol):
            self.selavol[...] = selavol
        
        @property
        def selbvol(self):
            """
            Element selbvol ftype=integer pytype=int
            
            
            Defined at options_mod.f90 line 19
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bader.f90wrap_options_obj__array__selbvol(self._handle)
            if array_handle in self._arrays:
                selbvol = self._arrays[array_handle]
            else:
                selbvol = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bader.f90wrap_options_obj__array__selbvol)
                self._arrays[array_handle] = selbvol
            return selbvol
        
        @selbvol.setter
        def selbvol(self, selbvol):
            self.selbvol[...] = selbvol
        
        @property
        def sumavol(self):
            """
            Element sumavol ftype=integer pytype=int
            
            
            Defined at options_mod.f90 line 19
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bader.f90wrap_options_obj__array__sumavol(self._handle)
            if array_handle in self._arrays:
                sumavol = self._arrays[array_handle]
            else:
                sumavol = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bader.f90wrap_options_obj__array__sumavol)
                self._arrays[array_handle] = sumavol
            return sumavol
        
        @sumavol.setter
        def sumavol(self, sumavol):
            self.sumavol[...] = sumavol
        
        @property
        def sumbvol(self):
            """
            Element sumbvol ftype=integer pytype=int
            
            
            Defined at options_mod.f90 line 19
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bader.f90wrap_options_obj__array__sumbvol(self._handle)
            if array_handle in self._arrays:
                sumbvol = self._arrays[array_handle]
            else:
                sumbvol = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bader.f90wrap_options_obj__array__sumbvol)
                self._arrays[array_handle] = sumbvol
            return sumbvol
        
        @sumbvol.setter
        def sumbvol(self, sumbvol):
            self.sumbvol[...] = sumbvol
        
        @property
        def vac_flag(self):
            """
            Element vac_flag ftype=logical pytype=bool
            
            
            Defined at options_mod.f90 line 20
            
            """
            return _bader.f90wrap_options_obj__get__vac_flag(self._handle)
        
        @vac_flag.setter
        def vac_flag(self, vac_flag):
            _bader.f90wrap_options_obj__set__vac_flag(self._handle, vac_flag)
        
        @property
        def bader_flag(self):
            """
            Element bader_flag ftype=logical pytype=bool
            
            
            Defined at options_mod.f90 line 21
            
            """
            return _bader.f90wrap_options_obj__get__bader_flag(self._handle)
        
        @bader_flag.setter
        def bader_flag(self, bader_flag):
            _bader.f90wrap_options_obj__set__bader_flag(self._handle, bader_flag)
        
        @property
        def voronoi_flag(self):
            """
            Element voronoi_flag ftype=logical pytype=bool
            
            
            Defined at options_mod.f90 line 21
            
            """
            return _bader.f90wrap_options_obj__get__voronoi_flag(self._handle)
        
        @voronoi_flag.setter
        def voronoi_flag(self, voronoi_flag):
            _bader.f90wrap_options_obj__set__voronoi_flag(self._handle, voronoi_flag)
        
        @property
        def dipole_flag(self):
            """
            Element dipole_flag ftype=logical pytype=bool
            
            
            Defined at options_mod.f90 line 21
            
            """
            return _bader.f90wrap_options_obj__get__dipole_flag(self._handle)
        
        @dipole_flag.setter
        def dipole_flag(self, dipole_flag):
            _bader.f90wrap_options_obj__set__dipole_flag(self._handle, dipole_flag)
        
        @property
        def ldos_flag(self):
            """
            Element ldos_flag ftype=logical pytype=bool
            
            
            Defined at options_mod.f90 line 21
            
            """
            return _bader.f90wrap_options_obj__get__ldos_flag(self._handle)
        
        @ldos_flag.setter
        def ldos_flag(self, ldos_flag):
            _bader.f90wrap_options_obj__set__ldos_flag(self._handle, ldos_flag)
        
        @property
        def print_all_bader(self):
            """
            Element print_all_bader ftype=logical pytype=bool
            
            
            Defined at options_mod.f90 line 22
            
            """
            return _bader.f90wrap_options_obj__get__print_all_bader(self._handle)
        
        @print_all_bader.setter
        def print_all_bader(self, print_all_bader):
            _bader.f90wrap_options_obj__set__print_all_bader(self._handle, print_all_bader)
        
        @property
        def print_all_atom(self):
            """
            Element print_all_atom ftype=logical pytype=bool
            
            
            Defined at options_mod.f90 line 22
            
            """
            return _bader.f90wrap_options_obj__get__print_all_atom(self._handle)
        
        @print_all_atom.setter
        def print_all_atom(self, print_all_atom):
            _bader.f90wrap_options_obj__set__print_all_atom(self._handle, print_all_atom)
        
        @property
        def print_sel_bader(self):
            """
            Element print_sel_bader ftype=logical pytype=bool
            
            
            Defined at options_mod.f90 line 23
            
            """
            return _bader.f90wrap_options_obj__get__print_sel_bader(self._handle)
        
        @print_sel_bader.setter
        def print_sel_bader(self, print_sel_bader):
            _bader.f90wrap_options_obj__set__print_sel_bader(self._handle, print_sel_bader)
        
        @property
        def print_sel_atom(self):
            """
            Element print_sel_atom ftype=logical pytype=bool
            
            
            Defined at options_mod.f90 line 23
            
            """
            return _bader.f90wrap_options_obj__get__print_sel_atom(self._handle)
        
        @print_sel_atom.setter
        def print_sel_atom(self, print_sel_atom):
            _bader.f90wrap_options_obj__set__print_sel_atom(self._handle, print_sel_atom)
        
        @property
        def print_sum_bader(self):
            """
            Element print_sum_bader ftype=logical pytype=bool
            
            
            Defined at options_mod.f90 line 24
            
            """
            return _bader.f90wrap_options_obj__get__print_sum_bader(self._handle)
        
        @print_sum_bader.setter
        def print_sum_bader(self, print_sum_bader):
            _bader.f90wrap_options_obj__set__print_sum_bader(self._handle, print_sum_bader)
        
        @property
        def print_sum_atom(self):
            """
            Element print_sum_atom ftype=logical pytype=bool
            
            
            Defined at options_mod.f90 line 24
            
            """
            return _bader.f90wrap_options_obj__get__print_sum_atom(self._handle)
        
        @print_sum_atom.setter
        def print_sum_atom(self, print_sum_atom):
            _bader.f90wrap_options_obj__set__print_sum_atom(self._handle, print_sum_atom)
        
        @property
        def print_bader_index(self):
            """
            Element print_bader_index ftype=logical pytype=bool
            
            
            Defined at options_mod.f90 line 25
            
            """
            return _bader.f90wrap_options_obj__get__print_bader_index(self._handle)
        
        @print_bader_index.setter
        def print_bader_index(self, print_bader_index):
            _bader.f90wrap_options_obj__set__print_bader_index(self._handle, \
                print_bader_index)
        
        @property
        def print_atom_index(self):
            """
            Element print_atom_index ftype=logical pytype=bool
            
            
            Defined at options_mod.f90 line 25
            
            """
            return _bader.f90wrap_options_obj__get__print_atom_index(self._handle)
        
        @print_atom_index.setter
        def print_atom_index(self, print_atom_index):
            _bader.f90wrap_options_obj__set__print_atom_index(self._handle, \
                print_atom_index)
        
        @property
        def verbose_flag(self):
            """
            Element verbose_flag ftype=logical pytype=bool
            
            
            Defined at options_mod.f90 line 26
            
            """
            return _bader.f90wrap_options_obj__get__verbose_flag(self._handle)
        
        @verbose_flag.setter
        def verbose_flag(self, verbose_flag):
            _bader.f90wrap_options_obj__set__verbose_flag(self._handle, verbose_flag)
        
        @property
        def ref_flag(self):
            """
            Element ref_flag ftype=logical pytype=bool
            
            
            Defined at options_mod.f90 line 26
            
            """
            return _bader.f90wrap_options_obj__get__ref_flag(self._handle)
        
        @ref_flag.setter
        def ref_flag(self, ref_flag):
            _bader.f90wrap_options_obj__set__ref_flag(self._handle, ref_flag)
        
        @property
        def find_critpoints_flag(self):
            """
            Element find_critpoints_flag ftype=logical pytype=bool
            
            
            Defined at options_mod.f90 line 26
            
            """
            return _bader.f90wrap_options_obj__get__find_critpoints_flag(self._handle)
        
        @find_critpoints_flag.setter
        def find_critpoints_flag(self, find_critpoints_flag):
            _bader.f90wrap_options_obj__set__find_critpoints_flag(self._handle, \
                find_critpoints_flag)
        
        @property
        def leastsquare_flag(self):
            """
            Element leastsquare_flag ftype=logical pytype=bool
            
            
            Defined at options_mod.f90 line 27
            
            """
            return _bader.f90wrap_options_obj__get__leastsquare_flag(self._handle)
        
        @leastsquare_flag.setter
        def leastsquare_flag(self, leastsquare_flag):
            _bader.f90wrap_options_obj__set__leastsquare_flag(self._handle, \
                leastsquare_flag)
        
        @property
        def print_surfaces_atoms(self):
            """
            Element print_surfaces_atoms ftype=logical pytype=bool
            
            
            Defined at options_mod.f90 line 28
            
            """
            return _bader.f90wrap_options_obj__get__print_surfaces_atoms(self._handle)
        
        @print_surfaces_atoms.setter
        def print_surfaces_atoms(self, print_surfaces_atoms):
            _bader.f90wrap_options_obj__set__print_surfaces_atoms(self._handle, \
                print_surfaces_atoms)
        
        @property
        def knob_tem(self):
            """
            Element knob_tem ftype=real(q2) pytype=float
            
            
            Defined at options_mod.f90 line 29
            
            """
            return _bader.f90wrap_options_obj__get__knob_tem(self._handle)
        
        @knob_tem.setter
        def knob_tem(self, knob_tem):
            _bader.f90wrap_options_obj__set__knob_tem(self._handle, knob_tem)
        
        @property
        def knob_distance(self):
            """
            Element knob_distance ftype=real(q2) pytype=float
            
            
            Defined at options_mod.f90 line 30
            
            """
            return _bader.f90wrap_options_obj__get__knob_distance(self._handle)
        
        @knob_distance.setter
        def knob_distance(self, knob_distance):
            _bader.f90wrap_options_obj__set__knob_distance(self._handle, knob_distance)
        
        @property
        def knob_newtonr(self):
            """
            Element knob_newtonr ftype=real(q2) pytype=float
            
            
            Defined at options_mod.f90 line 31
            
            """
            return _bader.f90wrap_options_obj__get__knob_newtonr(self._handle)
        
        @knob_newtonr.setter
        def knob_newtonr(self, knob_newtonr):
            _bader.f90wrap_options_obj__set__knob_newtonr(self._handle, knob_newtonr)
        
        @property
        def ismolecule(self):
            """
            Element ismolecule ftype=logical pytype=bool
            
            
            Defined at options_mod.f90 line 32
            
            """
            return _bader.f90wrap_options_obj__get__ismolecule(self._handle)
        
        @ismolecule.setter
        def ismolecule(self, ismolecule):
            _bader.f90wrap_options_obj__set__ismolecule(self._handle, ismolecule)
        
        @property
        def iscrystal(self):
            """
            Element iscrystal ftype=logical pytype=bool
            
            
            Defined at options_mod.f90 line 33
            
            """
            return _bader.f90wrap_options_obj__get__iscrystal(self._handle)
        
        @iscrystal.setter
        def iscrystal(self, iscrystal):
            _bader.f90wrap_options_obj__set__iscrystal(self._handle, iscrystal)
        
        def __str__(self):
            ret = ['<options_obj>{\n']
            ret.append('    chargefile : ')
            ret.append(repr(self.chargefile))
            ret.append(',\n    refchgfile : ')
            ret.append(repr(self.refchgfile))
            ret.append(',\n    badertol : ')
            ret.append(repr(self.badertol))
            ret.append(',\n    stepsize : ')
            ret.append(repr(self.stepsize))
            ret.append(',\n    vacval : ')
            ret.append(repr(self.vacval))
            ret.append(',\n    out_opt : ')
            ret.append(repr(self.out_opt))
            ret.append(',\n    out_auto : ')
            ret.append(repr(self.out_auto))
            ret.append(',\n    out_cube : ')
            ret.append(repr(self.out_cube))
            ret.append(',\n    out_chgcar4 : ')
            ret.append(repr(self.out_chgcar4))
            ret.append(',\n    out_chgcar5 : ')
            ret.append(repr(self.out_chgcar5))
            ret.append(',\n    in_opt : ')
            ret.append(repr(self.in_opt))
            ret.append(',\n    in_auto : ')
            ret.append(repr(self.in_auto))
            ret.append(',\n    in_cube : ')
            ret.append(repr(self.in_cube))
            ret.append(',\n    in_chgcar : ')
            ret.append(repr(self.in_chgcar))
            ret.append(',\n    in_chgcar4 : ')
            ret.append(repr(self.in_chgcar4))
            ret.append(',\n    in_chgcar5 : ')
            ret.append(repr(self.in_chgcar5))
            ret.append(',\n    bader_opt : ')
            ret.append(repr(self.bader_opt))
            ret.append(',\n    bader_offgrid : ')
            ret.append(repr(self.bader_offgrid))
            ret.append(',\n    bader_ongrid : ')
            ret.append(repr(self.bader_ongrid))
            ret.append(',\n    bader_neargrid : ')
            ret.append(repr(self.bader_neargrid))
            ret.append(',\n    bader_weight : ')
            ret.append(repr(self.bader_weight))
            ret.append(',\n    quit_opt : ')
            ret.append(repr(self.quit_opt))
            ret.append(',\n    quit_max : ')
            ret.append(repr(self.quit_max))
            ret.append(',\n    quit_known : ')
            ret.append(repr(self.quit_known))
            ret.append(',\n    refine_edge_itrs : ')
            ret.append(repr(self.refine_edge_itrs))
            ret.append(',\n    selanum : ')
            ret.append(repr(self.selanum))
            ret.append(',\n    selbnum : ')
            ret.append(repr(self.selbnum))
            ret.append(',\n    sumanum : ')
            ret.append(repr(self.sumanum))
            ret.append(',\n    sumbnum : ')
            ret.append(repr(self.sumbnum))
            ret.append(',\n    selavol : ')
            ret.append(repr(self.selavol))
            ret.append(',\n    selbvol : ')
            ret.append(repr(self.selbvol))
            ret.append(',\n    sumavol : ')
            ret.append(repr(self.sumavol))
            ret.append(',\n    sumbvol : ')
            ret.append(repr(self.sumbvol))
            ret.append(',\n    vac_flag : ')
            ret.append(repr(self.vac_flag))
            ret.append(',\n    bader_flag : ')
            ret.append(repr(self.bader_flag))
            ret.append(',\n    voronoi_flag : ')
            ret.append(repr(self.voronoi_flag))
            ret.append(',\n    dipole_flag : ')
            ret.append(repr(self.dipole_flag))
            ret.append(',\n    ldos_flag : ')
            ret.append(repr(self.ldos_flag))
            ret.append(',\n    print_all_bader : ')
            ret.append(repr(self.print_all_bader))
            ret.append(',\n    print_all_atom : ')
            ret.append(repr(self.print_all_atom))
            ret.append(',\n    print_sel_bader : ')
            ret.append(repr(self.print_sel_bader))
            ret.append(',\n    print_sel_atom : ')
            ret.append(repr(self.print_sel_atom))
            ret.append(',\n    print_sum_bader : ')
            ret.append(repr(self.print_sum_bader))
            ret.append(',\n    print_sum_atom : ')
            ret.append(repr(self.print_sum_atom))
            ret.append(',\n    print_bader_index : ')
            ret.append(repr(self.print_bader_index))
            ret.append(',\n    print_atom_index : ')
            ret.append(repr(self.print_atom_index))
            ret.append(',\n    verbose_flag : ')
            ret.append(repr(self.verbose_flag))
            ret.append(',\n    ref_flag : ')
            ret.append(repr(self.ref_flag))
            ret.append(',\n    find_critpoints_flag : ')
            ret.append(repr(self.find_critpoints_flag))
            ret.append(',\n    leastsquare_flag : ')
            ret.append(repr(self.leastsquare_flag))
            ret.append(',\n    print_surfaces_atoms : ')
            ret.append(repr(self.print_surfaces_atoms))
            ret.append(',\n    knob_tem : ')
            ret.append(repr(self.knob_tem))
            ret.append(',\n    knob_distance : ')
            ret.append(repr(self.knob_distance))
            ret.append(',\n    knob_newtonr : ')
            ret.append(repr(self.knob_newtonr))
            ret.append(',\n    ismolecule : ')
            ret.append(repr(self.ismolecule))
            ret.append(',\n    iscrystal : ')
            ret.append(repr(self.iscrystal))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def get_options(self):
        """
        get_options(self)
        
        
        Defined at options_mod.f90 lines 41-517
        
        Parameters
        ----------
        opts : Options_Obj
        
        """
        _bader.f90wrap_get_options(opts=self._handle)
    
    _dt_array_initialisers = []
    

options_mod = Options_Mod()

class Charge_Mod(f90wrap.runtime.FortranModule):
    """
    Module charge_mod
    
    
    Defined at charge_mod.f90 lines 5-524
    
    """
    @f90wrap.runtime.register_class("bader.charge_obj")
    class charge_obj(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=charge_obj)
        
        
        Defined at charge_mod.f90 lines 10-17
        
        """
        def __init__(self, handle=None):
            """
            self = Charge_Obj()
            
            
            Defined at charge_mod.f90 lines 10-17
            
            
            Returns
            -------
            this : Charge_Obj
            	Object to be constructed
            
            
            Automatically generated constructor for charge_obj
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _bader.f90wrap_charge_obj_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Charge_Obj
            
            
            Defined at charge_mod.f90 lines 10-17
            
            Parameters
            ----------
            this : Charge_Obj
            	Object to be destructed
            
            
            Automatically generated destructor for charge_obj
            """
            if self._alloc:
                _bader.f90wrap_charge_obj_finalise(this=self._handle)
        
        @property
        def rho(self):
            """
            Element rho ftype=real(q2) pytype=float
            
            
            Defined at charge_mod.f90 line 11
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bader.f90wrap_charge_obj__array__rho(self._handle)
            if array_handle in self._arrays:
                rho = self._arrays[array_handle]
            else:
                rho = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bader.f90wrap_charge_obj__array__rho)
                self._arrays[array_handle] = rho
            return rho
        
        @rho.setter
        def rho(self, rho):
            self.rho[...] = rho
        
        @property
        def lat2car(self):
            """
            Element lat2car ftype=real(q2) pytype=float
            
            
            Defined at charge_mod.f90 line 12
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bader.f90wrap_charge_obj__array__lat2car(self._handle)
            if array_handle in self._arrays:
                lat2car = self._arrays[array_handle]
            else:
                lat2car = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bader.f90wrap_charge_obj__array__lat2car)
                self._arrays[array_handle] = lat2car
            return lat2car
        
        @lat2car.setter
        def lat2car(self, lat2car):
            self.lat2car[...] = lat2car
        
        @property
        def car2lat(self):
            """
            Element car2lat ftype=real(q2) pytype=float
            
            
            Defined at charge_mod.f90 line 12
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bader.f90wrap_charge_obj__array__car2lat(self._handle)
            if array_handle in self._arrays:
                car2lat = self._arrays[array_handle]
            else:
                car2lat = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bader.f90wrap_charge_obj__array__car2lat)
                self._arrays[array_handle] = car2lat
            return car2lat
        
        @car2lat.setter
        def car2lat(self, car2lat):
            self.car2lat[...] = car2lat
        
        @property
        def lat_dist(self):
            """
            Element lat_dist ftype=real(q2) pytype=float
            
            
            Defined at charge_mod.f90 line 13
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bader.f90wrap_charge_obj__array__lat_dist(self._handle)
            if array_handle in self._arrays:
                lat_dist = self._arrays[array_handle]
            else:
                lat_dist = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bader.f90wrap_charge_obj__array__lat_dist)
                self._arrays[array_handle] = lat_dist
            return lat_dist
        
        @lat_dist.setter
        def lat_dist(self, lat_dist):
            self.lat_dist[...] = lat_dist
        
        @property
        def lat_i_dist(self):
            """
            Element lat_i_dist ftype=real(q2) pytype=float
            
            
            Defined at charge_mod.f90 line 13
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bader.f90wrap_charge_obj__array__lat_i_dist(self._handle)
            if array_handle in self._arrays:
                lat_i_dist = self._arrays[array_handle]
            else:
                lat_i_dist = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bader.f90wrap_charge_obj__array__lat_i_dist)
                self._arrays[array_handle] = lat_i_dist
            return lat_i_dist
        
        @lat_i_dist.setter
        def lat_i_dist(self, lat_i_dist):
            self.lat_i_dist[...] = lat_i_dist
        
        @property
        def org_lat(self):
            """
            Element org_lat ftype=real(q2) pytype=float
            
            
            Defined at charge_mod.f90 line 14
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bader.f90wrap_charge_obj__array__org_lat(self._handle)
            if array_handle in self._arrays:
                org_lat = self._arrays[array_handle]
            else:
                org_lat = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bader.f90wrap_charge_obj__array__org_lat)
                self._arrays[array_handle] = org_lat
            return org_lat
        
        @org_lat.setter
        def org_lat(self, org_lat):
            self.org_lat[...] = org_lat
        
        @property
        def org_car(self):
            """
            Element org_car ftype=real(q2) pytype=float
            
            
            Defined at charge_mod.f90 line 14
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bader.f90wrap_charge_obj__array__org_car(self._handle)
            if array_handle in self._arrays:
                org_car = self._arrays[array_handle]
            else:
                org_car = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bader.f90wrap_charge_obj__array__org_car)
                self._arrays[array_handle] = org_car
            return org_car
        
        @org_car.setter
        def org_car(self, org_car):
            self.org_car[...] = org_car
        
        @property
        def i_npts(self):
            """
            Element i_npts ftype=real(q2) pytype=float
            
            
            Defined at charge_mod.f90 line 15
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bader.f90wrap_charge_obj__array__i_npts(self._handle)
            if array_handle in self._arrays:
                i_npts = self._arrays[array_handle]
            else:
                i_npts = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bader.f90wrap_charge_obj__array__i_npts)
                self._arrays[array_handle] = i_npts
            return i_npts
        
        @i_npts.setter
        def i_npts(self, i_npts):
            self.i_npts[...] = i_npts
        
        @property
        def npts(self):
            """
            Element npts ftype=integer pytype=int
            
            
            Defined at charge_mod.f90 line 16
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bader.f90wrap_charge_obj__array__npts(self._handle)
            if array_handle in self._arrays:
                npts = self._arrays[array_handle]
            else:
                npts = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bader.f90wrap_charge_obj__array__npts)
                self._arrays[array_handle] = npts
            return npts
        
        @npts.setter
        def npts(self, npts):
            self.npts[...] = npts
        
        @property
        def nrho(self):
            """
            Element nrho ftype=integer  pytype=int
            
            
            Defined at charge_mod.f90 line 17
            
            """
            return _bader.f90wrap_charge_obj__get__nrho(self._handle)
        
        @nrho.setter
        def nrho(self, nrho):
            _bader.f90wrap_charge_obj__set__nrho(self._handle, nrho)
        
        def __str__(self):
            ret = ['<charge_obj>{\n']
            ret.append('    rho : ')
            ret.append(repr(self.rho))
            ret.append(',\n    lat2car : ')
            ret.append(repr(self.lat2car))
            ret.append(',\n    car2lat : ')
            ret.append(repr(self.car2lat))
            ret.append(',\n    lat_dist : ')
            ret.append(repr(self.lat_dist))
            ret.append(',\n    lat_i_dist : ')
            ret.append(repr(self.lat_i_dist))
            ret.append(',\n    org_lat : ')
            ret.append(repr(self.org_lat))
            ret.append(',\n    org_car : ')
            ret.append(repr(self.org_car))
            ret.append(',\n    i_npts : ')
            ret.append(repr(self.i_npts))
            ret.append(',\n    npts : ')
            ret.append(repr(self.npts))
            ret.append(',\n    nrho : ')
            ret.append(repr(self.nrho))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def rho_val(self, p1, p2, p3):
        """
        rho_val = rho_val(self, p1, p2, p3)
        
        
        Defined at charge_mod.f90 lines 54-72
        
        Parameters
        ----------
        chg : Charge_Obj
        p1 : int
        p2 : int
        p3 : int
        
        Returns
        -------
        rho_val : float
        
        """
        rho_val = _bader.f90wrap_rho_val(chg=self._handle, p1=p1, p2=p2, p3=p3)
        return rho_val
    
    @staticmethod
    def rho_grad(self, r):
        """
        rho_grad_lat, rho = rho_grad(self, r)
        
        
        Defined at charge_mod.f90 lines 77-131
        
        Parameters
        ----------
        chg : Charge_Obj
        r : float array
        
        Returns
        -------
        rho_grad_lat : float array
        rho : float
        
        """
        rho_grad_lat, rho = _bader.f90wrap_rho_grad(chg=self._handle, r=r)
        return rho_grad_lat, rho
    
    @staticmethod
    def rho_grad_dir(self, p):
        """
        rho_grad_dir = rho_grad_dir(self, p)
        
        
        Defined at charge_mod.f90 lines 137-168
        
        Parameters
        ----------
        chg : Charge_Obj
        p : int array
        
        Returns
        -------
        rho_grad_dir : float array
        
        """
        rho_grad_dir = _bader.f90wrap_rho_grad_dir(chg=self._handle, p=p)
        return rho_grad_dir
    
    @staticmethod
    def pbc_r_lat(r_lat, pmax):
        """
        pbc_r_lat(r_lat, pmax)
        
        
        Defined at charge_mod.f90 lines 173-187
        
        Parameters
        ----------
        r_lat : float array
        pmax : int array
        
        """
        _bader.f90wrap_pbc_r_lat(r_lat=r_lat, pmax=pmax)
    
    @staticmethod
    def pbc(p, pmax):
        """
        pbc(p, pmax)
        
        
        Defined at charge_mod.f90 lines 192-206
        
        Parameters
        ----------
        p : int array
        pmax : int array
        
        """
        _bader.f90wrap_pbc(p=p, pmax=pmax)
    
    @staticmethod
    def dpbc(dr, nf, nf_2):
        """
        dpbc(dr, nf, nf_2)
        
        
        Defined at charge_mod.f90 lines 211-225
        
        Parameters
        ----------
        dr : float array
        nf : float array
        nf_2 : float array
        
        """
        _bader.f90wrap_dpbc(dr=dr, nf=nf, nf_2=nf_2)
    
    @staticmethod
    def dpbc_dir_org(dr):
        """
        dpbc_dir_org(dr)
        
        
        Defined at charge_mod.f90 lines 230-243
        
        Parameters
        ----------
        dr : float array
        
        """
        _bader.f90wrap_dpbc_dir_org(dr=dr)
    
    @staticmethod
    def dpbc_dir(self, dr_dir):
        """
        dpbc_dir(self, dr_dir)
        
        
        Defined at charge_mod.f90 lines 248-259
        
        Parameters
        ----------
        ions : Ions_Obj
        dr_dir : float array
        
        """
        _bader.f90wrap_dpbc_dir(ions=self._handle, dr_dir=dr_dir)
    
    @staticmethod
    def dpbc_car(self, dr_car):
        """
        dpbc_car(self, dr_car)
        
        
        Defined at charge_mod.f90 lines 264-292
        
        Parameters
        ----------
        ions : Ions_Obj
        dr_car : float array
        
        """
        _bader.f90wrap_dpbc_car(ions=self._handle, dr_car=dr_car)
    
    @staticmethod
    def to_lat(self, r):
        """
        to_lat = to_lat(self, r)
        
        
        Defined at charge_mod.f90 lines 297-328
        
        Parameters
        ----------
        chg : Charge_Obj
        r : float array
        
        Returns
        -------
        to_lat : int array
        
        """
        to_lat = _bader.f90wrap_to_lat(chg=self._handle, r=r)
        return to_lat
    
    @staticmethod
    def is_max(self, p):
        """
        is_max = is_max(self, p)
        
        
        Defined at charge_mod.f90 lines 333-356
        
        Parameters
        ----------
        chg : Charge_Obj
        p : int array
        
        Returns
        -------
        is_max : bool
        
        """
        is_max = _bader.f90wrap_is_max(chg=self._handle, p=p)
        return is_max
    
    @staticmethod
    def is_max_ongrid(self, p):
        """
        is_max_ongrid = is_max_ongrid(self, p)
        
        
        Defined at charge_mod.f90 lines 361-379
        
        Parameters
        ----------
        chg : Charge_Obj
        p : int array
        
        Returns
        -------
        is_max_ongrid : bool
        
        """
        is_max_ongrid = _bader.f90wrap_is_max_ongrid(chg=self._handle, p=p)
        return is_max_ongrid
    
    @staticmethod
    def lat2dir(self, p):
        """
        lat2dir = lat2dir(self, p)
        
        
        Defined at charge_mod.f90 lines 384-391
        
        Parameters
        ----------
        chg : Charge_Obj
        p : float array
        
        Returns
        -------
        lat2dir : float array
        
        """
        lat2dir = _bader.f90wrap_lat2dir(chg=self._handle, p=p)
        return lat2dir
    
    @staticmethod
    def lat2car(self, p):
        """
        lat2car = lat2car(self, p)
        
        
        Defined at charge_mod.f90 lines 396-404
        
        Parameters
        ----------
        chg : Charge_Obj
        p : float array
        
        Returns
        -------
        lat2car : float array
        
        """
        lat2car = _bader.f90wrap_lat2car(chg=self._handle, p=p)
        return lat2car
    
    @staticmethod
    def dir2lat(self, p):
        """
        dir2lat = dir2lat(self, p)
        
        
        Defined at charge_mod.f90 lines 409-417
        
        Parameters
        ----------
        chg : Charge_Obj
        p : float array
        
        Returns
        -------
        dir2lat : float array
        
        """
        dir2lat = _bader.f90wrap_dir2lat(chg=self._handle, p=p)
        return dir2lat
    
    @staticmethod
    def car2lat(self, p):
        """
        car2lat = car2lat(self, p)
        
        
        Defined at charge_mod.f90 lines 422-430
        
        Parameters
        ----------
        chg : Charge_Obj
        p : float array
        
        Returns
        -------
        car2lat : float array
        
        """
        car2lat = _bader.f90wrap_car2lat(chg=self._handle, p=p)
        return car2lat
    
    _dt_array_initialisers = []
    

charge_mod = Charge_Mod()

class Chgcar_Mod(f90wrap.runtime.FortranModule):
    """
    Module chgcar_mod
    
    
    Defined at chgcar_mod.f90 lines 5-142
    
    """
    @staticmethod
    def read_charge_chgcar(self, chg, chargefile, opts):
        """
        read_charge_chgcar(self, chg, chargefile, opts)
        
        
        Defined at chgcar_mod.f90 lines 18-113
        
        Parameters
        ----------
        ions : Ions_Obj
        chg : Charge_Obj
        chargefile : str
        opts : Options_Obj
        
        """
        _bader.f90wrap_read_charge_chgcar(ions=self._handle, chg=chg._handle, \
            chargefile=chargefile, opts=opts._handle)
    
    @staticmethod
    def write_charge_chgcar(self, chg, chargefile, opts):
        """
        write_charge_chgcar(self, chg, chargefile, opts)
        
        
        Defined at chgcar_mod.f90 lines 118-141
        
        Parameters
        ----------
        ions : Ions_Obj
        chg : Charge_Obj
        chargefile : str
        opts : Options_Obj
        
        """
        _bader.f90wrap_write_charge_chgcar(ions=self._handle, chg=chg._handle, \
            chargefile=chargefile, opts=opts._handle)
    
    _dt_array_initialisers = []
    

chgcar_mod = Chgcar_Mod()

class Cube_Mod(f90wrap.runtime.FortranModule):
    """
    Module cube_mod
    
    
    Defined at cube_mod.f90 lines 5-151
    
    """
    @staticmethod
    def read_charge_cube(self, chg, chargefile):
        """
        read_charge_cube(self, chg, chargefile)
        
        
        Defined at cube_mod.f90 lines 18-109
        
        Parameters
        ----------
        ions : Ions_Obj
        chg : Charge_Obj
        chargefile : str
        
        """
        _bader.f90wrap_read_charge_cube(ions=self._handle, chg=chg._handle, \
            chargefile=chargefile)
    
    @staticmethod
    def write_charge_cube(self, chg, chargefile):
        """
        write_charge_cube(self, chg, chargefile)
        
        
        Defined at cube_mod.f90 lines 114-150
        
        Parameters
        ----------
        ions : Ions_Obj
        chg : Charge_Obj
        chargefile : str
        
        """
        _bader.f90wrap_write_charge_cube(ions=self._handle, chg=chg._handle, \
            chargefile=chargefile)
    
    _dt_array_initialisers = []
    

cube_mod = Cube_Mod()

class Io_Mod(f90wrap.runtime.FortranModule):
    """
    Module io_mod
    
    
    Defined at io_mod.f90 lines 5-101
    
    """
    @staticmethod
    def read_charge(self, chg, opts):
        """
        read_charge(self, chg, opts)
        
        
        Defined at io_mod.f90 lines 21-60
        
        Parameters
        ----------
        ions : Ions_Obj
        chg : Charge_Obj
        opts : Options_Obj
        
        """
        _bader.f90wrap_read_charge(ions=self._handle, chg=chg._handle, \
            opts=opts._handle)
    
    @staticmethod
    def read_charge_ref(self, chg, opts):
        """
        read_charge_ref(self, chg, opts)
        
        
        Defined at io_mod.f90 lines 66-81
        
        Parameters
        ----------
        ions : Ions_Obj
        chg : Charge_Obj
        opts : Options_Obj
        
        """
        _bader.f90wrap_read_charge_ref(ions=self._handle, chg=chg._handle, \
            opts=opts._handle)
    
    @staticmethod
    def write_charge(self, chg, opts, chargefile):
        """
        write_charge(self, chg, opts, chargefile)
        
        
        Defined at io_mod.f90 lines 86-100
        
        Parameters
        ----------
        ions : Ions_Obj
        chg : Charge_Obj
        opts : Options_Obj
        chargefile : str
        
        """
        _bader.f90wrap_write_charge(ions=self._handle, chg=chg._handle, \
            opts=opts._handle, chargefile=chargefile)
    
    _dt_array_initialisers = []
    

io_mod = Io_Mod()

class Bader_Mod(f90wrap.runtime.FortranModule):
    """
    Module bader_mod
    
    
    Defined at bader_mod.f90 lines 5-1661
    
    """
    @f90wrap.runtime.register_class("bader.bader_obj")
    class bader_obj(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=bader_obj)
        
        
        Defined at bader_mod.f90 lines 23-31
        
        """
        def __init__(self, handle=None):
            """
            self = Bader_Obj()
            
            
            Defined at bader_mod.f90 lines 23-31
            
            
            Returns
            -------
            this : Bader_Obj
            	Object to be constructed
            
            
            Automatically generated constructor for bader_obj
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _bader.f90wrap_bader_obj_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Bader_Obj
            
            
            Defined at bader_mod.f90 lines 23-31
            
            Parameters
            ----------
            this : Bader_Obj
            	Object to be destructed
            
            
            Automatically generated destructor for bader_obj
            """
            if self._alloc:
                _bader.f90wrap_bader_obj_finalise(this=self._handle)
        
        @property
        def volpos_lat(self):
            """
            Element volpos_lat ftype=real(q2) pytype=float
            
            
            Defined at bader_mod.f90 line 24
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bader.f90wrap_bader_obj__array__volpos_lat(self._handle)
            if array_handle in self._arrays:
                volpos_lat = self._arrays[array_handle]
            else:
                volpos_lat = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bader.f90wrap_bader_obj__array__volpos_lat)
                self._arrays[array_handle] = volpos_lat
            return volpos_lat
        
        @volpos_lat.setter
        def volpos_lat(self, volpos_lat):
            self.volpos_lat[...] = volpos_lat
        
        @property
        def volpos_car(self):
            """
            Element volpos_car ftype=real(q2) pytype=float
            
            
            Defined at bader_mod.f90 line 24
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bader.f90wrap_bader_obj__array__volpos_car(self._handle)
            if array_handle in self._arrays:
                volpos_car = self._arrays[array_handle]
            else:
                volpos_car = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bader.f90wrap_bader_obj__array__volpos_car)
                self._arrays[array_handle] = volpos_car
            return volpos_car
        
        @volpos_car.setter
        def volpos_car(self, volpos_car):
            self.volpos_car[...] = volpos_car
        
        @property
        def volpos_dir(self):
            """
            Element volpos_dir ftype=real(q2) pytype=float
            
            
            Defined at bader_mod.f90 line 24
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bader.f90wrap_bader_obj__array__volpos_dir(self._handle)
            if array_handle in self._arrays:
                volpos_dir = self._arrays[array_handle]
            else:
                volpos_dir = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bader.f90wrap_bader_obj__array__volpos_dir)
                self._arrays[array_handle] = volpos_dir
            return volpos_dir
        
        @volpos_dir.setter
        def volpos_dir(self, volpos_dir):
            self.volpos_dir[...] = volpos_dir
        
        @property
        def volchg(self):
            """
            Element volchg ftype=real(q2) pytype=float
            
            
            Defined at bader_mod.f90 line 25
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bader.f90wrap_bader_obj__array__volchg(self._handle)
            if array_handle in self._arrays:
                volchg = self._arrays[array_handle]
            else:
                volchg = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bader.f90wrap_bader_obj__array__volchg)
                self._arrays[array_handle] = volchg
            return volchg
        
        @volchg.setter
        def volchg(self, volchg):
            self.volchg[...] = volchg
        
        @property
        def iondist(self):
            """
            Element iondist ftype=real(q2) pytype=float
            
            
            Defined at bader_mod.f90 line 25
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bader.f90wrap_bader_obj__array__iondist(self._handle)
            if array_handle in self._arrays:
                iondist = self._arrays[array_handle]
            else:
                iondist = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bader.f90wrap_bader_obj__array__iondist)
                self._arrays[array_handle] = iondist
            return iondist
        
        @iondist.setter
        def iondist(self, iondist):
            self.iondist[...] = iondist
        
        @property
        def ionchg(self):
            """
            Element ionchg ftype=real(q2) pytype=float
            
            
            Defined at bader_mod.f90 line 25
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bader.f90wrap_bader_obj__array__ionchg(self._handle)
            if array_handle in self._arrays:
                ionchg = self._arrays[array_handle]
            else:
                ionchg = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bader.f90wrap_bader_obj__array__ionchg)
                self._arrays[array_handle] = ionchg
            return ionchg
        
        @ionchg.setter
        def ionchg(self, ionchg):
            self.ionchg[...] = ionchg
        
        @property
        def minsurfdist(self):
            """
            Element minsurfdist ftype=real(q2) pytype=float
            
            
            Defined at bader_mod.f90 line 25
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bader.f90wrap_bader_obj__array__minsurfdist(self._handle)
            if array_handle in self._arrays:
                minsurfdist = self._arrays[array_handle]
            else:
                minsurfdist = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bader.f90wrap_bader_obj__array__minsurfdist)
                self._arrays[array_handle] = minsurfdist
            return minsurfdist
        
        @minsurfdist.setter
        def minsurfdist(self, minsurfdist):
            self.minsurfdist[...] = minsurfdist
        
        @property
        def ionvol(self):
            """
            Element ionvol ftype=real(q2) pytype=float
            
            
            Defined at bader_mod.f90 line 25
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bader.f90wrap_bader_obj__array__ionvol(self._handle)
            if array_handle in self._arrays:
                ionvol = self._arrays[array_handle]
            else:
                ionvol = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bader.f90wrap_bader_obj__array__ionvol)
                self._arrays[array_handle] = ionvol
            return ionvol
        
        @ionvol.setter
        def ionvol(self, ionvol):
            self.ionvol[...] = ionvol
        
        @property
        def volnum(self):
            """
            Element volnum ftype=integer pytype=int
            
            
            Defined at bader_mod.f90 line 26
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bader.f90wrap_bader_obj__array__volnum(self._handle)
            if array_handle in self._arrays:
                volnum = self._arrays[array_handle]
            else:
                volnum = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bader.f90wrap_bader_obj__array__volnum)
                self._arrays[array_handle] = volnum
            return volnum
        
        @volnum.setter
        def volnum(self, volnum):
            self.volnum[...] = volnum
        
        @property
        def known(self):
            """
            Element known ftype=integer pytype=int
            
            
            Defined at bader_mod.f90 line 26
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bader.f90wrap_bader_obj__array__known(self._handle)
            if array_handle in self._arrays:
                known = self._arrays[array_handle]
            else:
                known = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bader.f90wrap_bader_obj__array__known)
                self._arrays[array_handle] = known
            return known
        
        @known.setter
        def known(self, known):
            self.known[...] = known
        
        @property
        def path(self):
            """
            Element path ftype=integer pytype=int
            
            
            Defined at bader_mod.f90 line 27
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bader.f90wrap_bader_obj__array__path(self._handle)
            if array_handle in self._arrays:
                path = self._arrays[array_handle]
            else:
                path = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bader.f90wrap_bader_obj__array__path)
                self._arrays[array_handle] = path
            return path
        
        @path.setter
        def path(self, path):
            self.path[...] = path
        
        @property
        def nnion(self):
            """
            Element nnion ftype=integer pytype=int
            
            
            Defined at bader_mod.f90 line 28
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bader.f90wrap_bader_obj__array__nnion(self._handle)
            if array_handle in self._arrays:
                nnion = self._arrays[array_handle]
            else:
                nnion = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bader.f90wrap_bader_obj__array__nnion)
                self._arrays[array_handle] = nnion
            return nnion
        
        @nnion.setter
        def nnion(self, nnion):
            self.nnion[...] = nnion
        
        @property
        def stepsize(self):
            """
            Element stepsize ftype=real(q2) pytype=float
            
            
            Defined at bader_mod.f90 line 29
            
            """
            return _bader.f90wrap_bader_obj__get__stepsize(self._handle)
        
        @stepsize.setter
        def stepsize(self, stepsize):
            _bader.f90wrap_bader_obj__set__stepsize(self._handle, stepsize)
        
        @property
        def tol(self):
            """
            Element tol ftype=real(q2) pytype=float
            
            
            Defined at bader_mod.f90 line 29
            
            """
            return _bader.f90wrap_bader_obj__get__tol(self._handle)
        
        @tol.setter
        def tol(self, tol):
            _bader.f90wrap_bader_obj__set__tol(self._handle, tol)
        
        @property
        def vacchg(self):
            """
            Element vacchg ftype=real(q2) pytype=float
            
            
            Defined at bader_mod.f90 line 30
            
            """
            return _bader.f90wrap_bader_obj__get__vacchg(self._handle)
        
        @vacchg.setter
        def vacchg(self, vacchg):
            _bader.f90wrap_bader_obj__set__vacchg(self._handle, vacchg)
        
        @property
        def vacvol(self):
            """
            Element vacvol ftype=real(q2) pytype=float
            
            
            Defined at bader_mod.f90 line 30
            
            """
            return _bader.f90wrap_bader_obj__get__vacvol(self._handle)
        
        @vacvol.setter
        def vacvol(self, vacvol):
            _bader.f90wrap_bader_obj__set__vacvol(self._handle, vacvol)
        
        @property
        def nvols(self):
            """
            Element nvols ftype=integer  pytype=int
            
            
            Defined at bader_mod.f90 line 31
            
            """
            return _bader.f90wrap_bader_obj__get__nvols(self._handle)
        
        @nvols.setter
        def nvols(self, nvols):
            _bader.f90wrap_bader_obj__set__nvols(self._handle, nvols)
        
        @property
        def pnum(self):
            """
            Element pnum ftype=integer  pytype=int
            
            
            Defined at bader_mod.f90 line 31
            
            """
            return _bader.f90wrap_bader_obj__get__pnum(self._handle)
        
        @pnum.setter
        def pnum(self, pnum):
            _bader.f90wrap_bader_obj__set__pnum(self._handle, pnum)
        
        @property
        def bnum(self):
            """
            Element bnum ftype=integer  pytype=int
            
            
            Defined at bader_mod.f90 line 31
            
            """
            return _bader.f90wrap_bader_obj__get__bnum(self._handle)
        
        @bnum.setter
        def bnum(self, bnum):
            _bader.f90wrap_bader_obj__set__bnum(self._handle, bnum)
        
        @property
        def pdim(self):
            """
            Element pdim ftype=integer  pytype=int
            
            
            Defined at bader_mod.f90 line 31
            
            """
            return _bader.f90wrap_bader_obj__get__pdim(self._handle)
        
        @pdim.setter
        def pdim(self, pdim):
            _bader.f90wrap_bader_obj__set__pdim(self._handle, pdim)
        
        @property
        def bdim(self):
            """
            Element bdim ftype=integer  pytype=int
            
            
            Defined at bader_mod.f90 line 31
            
            """
            return _bader.f90wrap_bader_obj__get__bdim(self._handle)
        
        @bdim.setter
        def bdim(self, bdim):
            _bader.f90wrap_bader_obj__set__bdim(self._handle, bdim)
        
        @property
        def refine_edge_itrs(self):
            """
            Element refine_edge_itrs ftype=integer  pytype=int
            
            
            Defined at bader_mod.f90 line 31
            
            """
            return _bader.f90wrap_bader_obj__get__refine_edge_itrs(self._handle)
        
        @refine_edge_itrs.setter
        def refine_edge_itrs(self, refine_edge_itrs):
            _bader.f90wrap_bader_obj__set__refine_edge_itrs(self._handle, refine_edge_itrs)
        
        def __str__(self):
            ret = ['<bader_obj>{\n']
            ret.append('    volpos_lat : ')
            ret.append(repr(self.volpos_lat))
            ret.append(',\n    volpos_car : ')
            ret.append(repr(self.volpos_car))
            ret.append(',\n    volpos_dir : ')
            ret.append(repr(self.volpos_dir))
            ret.append(',\n    volchg : ')
            ret.append(repr(self.volchg))
            ret.append(',\n    iondist : ')
            ret.append(repr(self.iondist))
            ret.append(',\n    ionchg : ')
            ret.append(repr(self.ionchg))
            ret.append(',\n    minsurfdist : ')
            ret.append(repr(self.minsurfdist))
            ret.append(',\n    ionvol : ')
            ret.append(repr(self.ionvol))
            ret.append(',\n    volnum : ')
            ret.append(repr(self.volnum))
            ret.append(',\n    known : ')
            ret.append(repr(self.known))
            ret.append(',\n    path : ')
            ret.append(repr(self.path))
            ret.append(',\n    nnion : ')
            ret.append(repr(self.nnion))
            ret.append(',\n    stepsize : ')
            ret.append(repr(self.stepsize))
            ret.append(',\n    tol : ')
            ret.append(repr(self.tol))
            ret.append(',\n    vacchg : ')
            ret.append(repr(self.vacchg))
            ret.append(',\n    vacvol : ')
            ret.append(repr(self.vacvol))
            ret.append(',\n    nvols : ')
            ret.append(repr(self.nvols))
            ret.append(',\n    pnum : ')
            ret.append(repr(self.pnum))
            ret.append(',\n    bnum : ')
            ret.append(repr(self.bnum))
            ret.append(',\n    pdim : ')
            ret.append(repr(self.pdim))
            ret.append(',\n    bdim : ')
            ret.append(repr(self.bdim))
            ret.append(',\n    refine_edge_itrs : ')
            ret.append(repr(self.refine_edge_itrs))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def bader_calc(self, ions, chgval, opts):
        """
        bader_calc(self, ions, chgval, opts)
        
        
        Defined at bader_mod.f90 lines 47-233
        
        Parameters
        ----------
        bdr : Bader_Obj
        ions : Ions_Obj
        chgval : Charge_Obj
        opts : Options_Obj
        
        """
        _bader.f90wrap_bader_calc(bdr=self._handle, ions=ions._handle, \
            chgval=chgval._handle, opts=opts._handle)
    
    @staticmethod
    def assign_chg2atom(self, ions, chg):
        """
        assign_chg2atom(self, ions, chg)
        
        
        Defined at bader_mod.f90 lines 563-593
        
        Parameters
        ----------
        bdr : Bader_Obj
        ions : Ions_Obj
        chg : Charge_Obj
        
        """
        _bader.f90wrap_assign_chg2atom(bdr=self._handle, ions=ions._handle, \
            chg=chg._handle)
    
    @staticmethod
    def bader_mindist(self, ions, chg):
        """
        bader_mindist(self, ions, chg)
        
        
        Defined at bader_mod.f90 lines 599-649
        
        Parameters
        ----------
        bdr : Bader_Obj
        ions : Ions_Obj
        chg : Charge_Obj
        
        """
        _bader.f90wrap_bader_mindist(bdr=self._handle, ions=ions._handle, \
            chg=chg._handle)
    
    @staticmethod
    def write_all_bader(self, opts, ions, chg):
        """
        write_all_bader(self, opts, ions, chg)
        
        
        Defined at bader_mod.f90 lines 654-694
        
        Parameters
        ----------
        bdr : Bader_Obj
        opts : Options_Obj
        ions : Ions_Obj
        chg : Charge_Obj
        
        """
        _bader.f90wrap_write_all_bader(bdr=self._handle, opts=opts._handle, \
            ions=ions._handle, chg=chg._handle)
    
    @staticmethod
    def write_all_atom(self, opts, ions, chg):
        """
        write_all_atom(self, opts, ions, chg)
        
        
        Defined at bader_mod.f90 lines 701-755
        
        Parameters
        ----------
        bdr : Bader_Obj
        opts : Options_Obj
        ions : Ions_Obj
        chg : Charge_Obj
        
        """
        _bader.f90wrap_write_all_atom(bdr=self._handle, opts=opts._handle, \
            ions=ions._handle, chg=chg._handle)
    
    @staticmethod
    def write_sel_atom(self, opts, ions, chg):
        """
        write_sel_atom(self, opts, ions, chg)
        
        
        Defined at bader_mod.f90 lines 760-812
        
        Parameters
        ----------
        bdr : Bader_Obj
        opts : Options_Obj
        ions : Ions_Obj
        chg : Charge_Obj
        
        """
        _bader.f90wrap_write_sel_atom(bdr=self._handle, opts=opts._handle, \
            ions=ions._handle, chg=chg._handle)
    
    @staticmethod
    def write_sum_atom(self, opts, ions, chg):
        """
        write_sum_atom(self, opts, ions, chg)
        
        
        Defined at bader_mod.f90 lines 818-874
        
        Parameters
        ----------
        bdr : Bader_Obj
        opts : Options_Obj
        ions : Ions_Obj
        chg : Charge_Obj
        
        """
        _bader.f90wrap_write_sum_atom(bdr=self._handle, opts=opts._handle, \
            ions=ions._handle, chg=chg._handle)
    
    @staticmethod
    def write_sel_bader(self, opts, ions, chg):
        """
        write_sel_bader(self, opts, ions, chg)
        
        
        Defined at bader_mod.f90 lines 880-929
        
        Parameters
        ----------
        bdr : Bader_Obj
        opts : Options_Obj
        ions : Ions_Obj
        chg : Charge_Obj
        
        """
        _bader.f90wrap_write_sel_bader(bdr=self._handle, opts=opts._handle, \
            ions=ions._handle, chg=chg._handle)
    
    @staticmethod
    def write_sum_bader(self, opts, ions, chg):
        """
        write_sum_bader(self, opts, ions, chg)
        
        
        Defined at bader_mod.f90 lines 935-983
        
        Parameters
        ----------
        bdr : Bader_Obj
        opts : Options_Obj
        ions : Ions_Obj
        chg : Charge_Obj
        
        """
        _bader.f90wrap_write_sum_bader(bdr=self._handle, opts=opts._handle, \
            ions=ions._handle, chg=chg._handle)
    
    @staticmethod
    def write_bader_index(self, opts, ions, chg):
        """
        write_bader_index(self, opts, ions, chg)
        
        
        Defined at bader_mod.f90 lines 989-1023
        
        Parameters
        ----------
        bdr : Bader_Obj
        opts : Options_Obj
        ions : Ions_Obj
        chg : Charge_Obj
        
        """
        _bader.f90wrap_write_bader_index(bdr=self._handle, opts=opts._handle, \
            ions=ions._handle, chg=chg._handle)
    
    @staticmethod
    def write_atom_index(self, opts, ions, chg):
        """
        write_atom_index(self, opts, ions, chg)
        
        
        Defined at bader_mod.f90 lines 1029-1066
        
        Parameters
        ----------
        bdr : Bader_Obj
        opts : Options_Obj
        ions : Ions_Obj
        chg : Charge_Obj
        
        """
        _bader.f90wrap_write_atom_index(bdr=self._handle, opts=opts._handle, \
            ions=ions._handle, chg=chg._handle)
    
    @staticmethod
    def cal_atomic_vol(self, ions, chg):
        """
        cal_atomic_vol(self, ions, chg)
        
        
        Defined at bader_mod.f90 lines 1146-1170
        
        Parameters
        ----------
        bdr : Bader_Obj
        ions : Ions_Obj
        chg : Charge_Obj
        
        """
        _bader.f90wrap_cal_atomic_vol(bdr=self._handle, ions=ions._handle, \
            chg=chg._handle)
    
    @staticmethod
    def bader_output(self, ions, chg):
        """
        bader_output(self, ions, chg)
        
        
        Defined at bader_mod.f90 lines 1180-1252
        
        Parameters
        ----------
        bdr : Bader_Obj
        ions : Ions_Obj
        chg : Charge_Obj
        
        """
        _bader.f90wrap_bader_output(bdr=self._handle, ions=ions._handle, \
            chg=chg._handle)
    
    @staticmethod
    def reallocate_volpos(self, newsize):
        """
        reallocate_volpos(self, newsize)
        
        
        Defined at bader_mod.f90 lines 1566-1580
        
        Parameters
        ----------
        bdr : Bader_Obj
        newsize : int
        
        """
        _bader.f90wrap_reallocate_volpos(bdr=self._handle, newsize=newsize)
    
    @staticmethod
    def bader_check_partitioning(self, chgval, scell_dir_out):
        """
        bader_check_partitioning(self, chgval, scell_dir_out)
        
        
        Defined at bader_mod.f90 lines 1602-1661
        
        Parameters
        ----------
        bdr : Bader_Obj
        chgval : Charge_Obj
        scell_dir_out : bool array
        
        """
        _bader.f90wrap_bader_check_partitioning(bdr=self._handle, chgval=chgval._handle, \
            scell_dir_out=scell_dir_out)
    
    _dt_array_initialisers = []
    

bader_mod = Bader_Mod()

class Voronoi_Mod(f90wrap.runtime.FortranModule):
    """
    Module voronoi_mod
    
    
    Defined at voronoi_mod.f90 lines 5-97
    
    """
    @f90wrap.runtime.register_class("bader.voronoi_obj")
    class voronoi_obj(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=voronoi_obj)
        
        
        Defined at voronoi_mod.f90 lines 14-15
        
        """
        def __init__(self, handle=None):
            """
            self = Voronoi_Obj()
            
            
            Defined at voronoi_mod.f90 lines 14-15
            
            
            Returns
            -------
            this : Voronoi_Obj
            	Object to be constructed
            
            
            Automatically generated constructor for voronoi_obj
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _bader.f90wrap_voronoi_obj_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Voronoi_Obj
            
            
            Defined at voronoi_mod.f90 lines 14-15
            
            Parameters
            ----------
            this : Voronoi_Obj
            	Object to be destructed
            
            
            Automatically generated destructor for voronoi_obj
            """
            if self._alloc:
                _bader.f90wrap_voronoi_obj_finalise(this=self._handle)
        
        @property
        def vorchg(self):
            """
            Element vorchg ftype=real(q2) pytype=float
            
            
            Defined at voronoi_mod.f90 line 15
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _bader.f90wrap_voronoi_obj__array__vorchg(self._handle)
            if array_handle in self._arrays:
                vorchg = self._arrays[array_handle]
            else:
                vorchg = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _bader.f90wrap_voronoi_obj__array__vorchg)
                self._arrays[array_handle] = vorchg
            return vorchg
        
        @vorchg.setter
        def vorchg(self, vorchg):
            self.vorchg[...] = vorchg
        
        def __str__(self):
            ret = ['<voronoi_obj>{\n']
            ret.append('    vorchg : ')
            ret.append(repr(self.vorchg))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def voronoi(self, ions, chg):
        """
        voronoi(self, ions, chg)
        
        
        Defined at voronoi_mod.f90 lines 26-96
        
        Parameters
        ----------
        vor : Voronoi_Obj
        ions : Ions_Obj
        chg : Charge_Obj
        
        """
        _bader.f90wrap_voronoi(vor=self._handle, ions=ions._handle, chg=chg._handle)
    
    _dt_array_initialisers = []
    

voronoi_mod = Voronoi_Mod()

class Multipole_Mod(f90wrap.runtime.FortranModule):
    """
    Module multipole_mod
    
    
    Defined at multipole_mod.f90 lines 6-146
    
    """
    @staticmethod
    def multipole_calc(self, ions, chgval, opts):
        """
        multipole_calc(self, ions, chgval, opts)
        
        
        Defined at multipole_mod.f90 lines 19-146
        
        Parameters
        ----------
        bdr : Bader_Obj
        ions : Ions_Obj
        chgval : Charge_Obj
        opts : Options_Obj
        
        """
        _bader.f90wrap_multipole_calc(bdr=self._handle, ions=ions._handle, \
            chgval=chgval._handle, opts=opts._handle)
    
    _dt_array_initialisers = []
    

multipole_mod = Multipole_Mod()


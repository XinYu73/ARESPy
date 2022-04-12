"""
Module grid_module


Defined at Grid_module.fpp lines 5-1248

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("AresMainPy_pkg.grid_type")
class grid_type(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=grid_type)
    
    
    Defined at Grid_module.fpp lines 14-23
    
    """
    def __init__(self, handle=None):
        """
        self = Grid_Type()
        
        
        Defined at Grid_module.fpp lines 14-23
        
        
        Returns
        -------
        this : Grid_Type
        	Object to be constructed
        
        
        Automatically generated constructor for grid_type
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _AresMainPy_pkg.f90wrap_grid_type_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Grid_Type
        
        
        Defined at Grid_module.fpp lines 14-23
        
        Parameters
        ----------
        this : Grid_Type
        	Object to be destructed
        
        
        Automatically generated destructor for grid_type
        """
        if self._alloc:
            _AresMainPy_pkg.f90wrap_grid_type_finalise(this=self._handle)
    
    @property
    def rhos(self):
        """
        Element rhos ftype=real(dp) pytype=float
        
        
        Defined at Grid_module.fpp line 16
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_grid_type__array__rhos(self._handle)
        if array_handle in self._arrays:
            rhos = self._arrays[array_handle]
        else:
            rhos = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_grid_type__array__rhos)
            self._arrays[array_handle] = rhos
        return rhos
    
    @rhos.setter
    def rhos(self, rhos):
        self.rhos[...] = rhos
    
    @property
    def vlpp(self):
        """
        Element vlpp ftype=real(dp) pytype=float
        
        
        Defined at Grid_module.fpp line 17
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_grid_type__array__vlpp(self._handle)
        if array_handle in self._arrays:
            vlpp = self._arrays[array_handle]
        else:
            vlpp = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_grid_type__array__vlpp)
            self._arrays[array_handle] = vlpp
        return vlpp
    
    @vlpp.setter
    def vlpp(self, vlpp):
        self.vlpp[...] = vlpp
    
    @property
    def eval(self):
        """
        Element eval ftype=real(dp) pytype=float
        
        
        Defined at Grid_module.fpp line 18
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_grid_type__array__eval(self._handle)
        if array_handle in self._arrays:
            eval = self._arrays[array_handle]
        else:
            eval = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_grid_type__array__eval)
            self._arrays[array_handle] = eval
        return eval
    
    @eval.setter
    def eval(self, eval):
        self.eval[...] = eval
    
    @property
    def gvec(self):
        """
        Element gvec ftype=real(dp) pytype=float
        
        
        Defined at Grid_module.fpp line 20
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_grid_type__array__gvec(self._handle)
        if array_handle in self._arrays:
            gvec = self._arrays[array_handle]
        else:
            gvec = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_grid_type__array__gvec)
            self._arrays[array_handle] = gvec
        return gvec
    
    @gvec.setter
    def gvec(self, gvec):
        self.gvec[...] = gvec
    
    @property
    def gmask(self):
        """
        Element gmask ftype=logical pytype=bool
        
        
        Defined at Grid_module.fpp line 21
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_grid_type__array__gmask(self._handle)
        if array_handle in self._arrays:
            gmask = self._arrays[array_handle]
        else:
            gmask = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_grid_type__array__gmask)
            self._arrays[array_handle] = gmask
        return gmask
    
    @gmask.setter
    def gmask(self, gmask):
        self.gmask[...] = gmask
    
    @property
    def rvec(self):
        """
        Element rvec ftype=real(dp) pytype=float
        
        
        Defined at Grid_module.fpp line 23
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_grid_type__array__rvec(self._handle)
        if array_handle in self._arrays:
            rvec = self._arrays[array_handle]
        else:
            rvec = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_grid_type__array__rvec)
            self._arrays[array_handle] = rvec
        return rvec
    
    @rvec.setter
    def rvec(self, rvec):
        self.rvec[...] = rvec
    
    def __str__(self):
        ret = ['<grid_type>{\n']
        ret.append('    rhos : ')
        ret.append(repr(self.rhos))
        ret.append(',\n    vlpp : ')
        ret.append(repr(self.vlpp))
        ret.append(',\n    eval : ')
        ret.append(repr(self.eval))
        ret.append(',\n    gvec : ')
        ret.append(repr(self.gvec))
        ret.append(',\n    gmask : ')
        ret.append(repr(self.gmask))
        ret.append(',\n    rvec : ')
        ret.append(repr(self.rvec))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("AresMainPy_pkg.kgrid_type")
class kgrid_type(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=kgrid_type)
    
    
    Defined at Grid_module.fpp lines 26-30
    
    """
    def __init__(self, handle=None):
        """
        self = Kgrid_Type()
        
        
        Defined at Grid_module.fpp lines 26-30
        
        
        Returns
        -------
        this : Kgrid_Type
        	Object to be constructed
        
        
        Automatically generated constructor for kgrid_type
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _AresMainPy_pkg.f90wrap_kgrid_type_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Kgrid_Type
        
        
        Defined at Grid_module.fpp lines 26-30
        
        Parameters
        ----------
        this : Kgrid_Type
        	Object to be destructed
        
        
        Automatically generated destructor for kgrid_type
        """
        if self._alloc:
            _AresMainPy_pkg.f90wrap_kgrid_type_finalise(this=self._handle)
    
    @property
    def vec(self):
        """
        Element vec ftype=real(dp) pytype=float
        
        
        Defined at Grid_module.fpp line 28
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_kgrid_type__array__vec(self._handle)
        if array_handle in self._arrays:
            vec = self._arrays[array_handle]
        else:
            vec = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_kgrid_type__array__vec)
            self._arrays[array_handle] = vec
        return vec
    
    @vec.setter
    def vec(self, vec):
        self.vec[...] = vec
    
    @property
    def vcar(self):
        """
        Element vcar ftype=real(dp) pytype=float
        
        
        Defined at Grid_module.fpp line 29
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_kgrid_type__array__vcar(self._handle)
        if array_handle in self._arrays:
            vcar = self._arrays[array_handle]
        else:
            vcar = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_kgrid_type__array__vcar)
            self._arrays[array_handle] = vcar
        return vcar
    
    @vcar.setter
    def vcar(self, vcar):
        self.vcar[...] = vcar
    
    @property
    def wk(self):
        """
        Element wk ftype=real(dp) pytype=float
        
        
        Defined at Grid_module.fpp line 30
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_kgrid_type__array__wk(self._handle)
        if array_handle in self._arrays:
            wk = self._arrays[array_handle]
        else:
            wk = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_kgrid_type__array__wk)
            self._arrays[array_handle] = wk
        return wk
    
    @wk.setter
    def wk(self, wk):
        self.wk[...] = wk
    
    def __str__(self):
        ret = ['<kgrid_type>{\n']
        ret.append('    vec : ')
        ret.append(repr(self.vec))
        ret.append(',\n    vcar : ')
        ret.append(repr(self.vcar))
        ret.append(',\n    wk : ')
        ret.append(repr(self.wk))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("AresMainPy_pkg.eigen_type")
class eigen_type(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=eigen_type)
    
    
    Defined at Grid_module.fpp lines 33-37
    
    """
    def __init__(self, handle=None):
        """
        self = Eigen_Type()
        
        
        Defined at Grid_module.fpp lines 33-37
        
        
        Returns
        -------
        this : Eigen_Type
        	Object to be constructed
        
        
        Automatically generated constructor for eigen_type
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _AresMainPy_pkg.f90wrap_eigen_type_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Eigen_Type
        
        
        Defined at Grid_module.fpp lines 33-37
        
        Parameters
        ----------
        this : Eigen_Type
        	Object to be destructed
        
        
        Automatically generated destructor for eigen_type
        """
        if self._alloc:
            _AresMainPy_pkg.f90wrap_eigen_type_finalise(this=self._handle)
    
    @property
    def wvf(self):
        """
        Element wvf ftype=complex(dp) pytype=complex
        
        
        Defined at Grid_module.fpp line 35
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_eigen_type__array__wvf(self._handle)
        if array_handle in self._arrays:
            wvf = self._arrays[array_handle]
        else:
            wvf = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_eigen_type__array__wvf)
            self._arrays[array_handle] = wvf
        return wvf
    
    @wvf.setter
    def wvf(self, wvf):
        self.wvf[...] = wvf
    
    @property
    def val(self):
        """
        Element val ftype=real(dp) pytype=float
        
        
        Defined at Grid_module.fpp line 37
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_eigen_type__array__val(self._handle)
        if array_handle in self._arrays:
            val = self._arrays[array_handle]
        else:
            val = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_eigen_type__array__val)
            self._arrays[array_handle] = val
        return val
    
    @val.setter
    def val(self, val):
        self.val[...] = val
    
    def __str__(self):
        ret = ['<eigen_type>{\n']
        ret.append('    wvf : ')
        ret.append(repr(self.wvf))
        ret.append(',\n    val : ')
        ret.append(repr(self.val))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("AresMainPy_pkg.eigen_type_r")
class eigen_type_r(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=eigen_type_r)
    
    
    Defined at Grid_module.fpp lines 40-44
    
    """
    def __init__(self, handle=None):
        """
        self = Eigen_Type_R()
        
        
        Defined at Grid_module.fpp lines 40-44
        
        
        Returns
        -------
        this : Eigen_Type_R
        	Object to be constructed
        
        
        Automatically generated constructor for eigen_type_r
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _AresMainPy_pkg.f90wrap_eigen_type_r_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Eigen_Type_R
        
        
        Defined at Grid_module.fpp lines 40-44
        
        Parameters
        ----------
        this : Eigen_Type_R
        	Object to be destructed
        
        
        Automatically generated destructor for eigen_type_r
        """
        if self._alloc:
            _AresMainPy_pkg.f90wrap_eigen_type_r_finalise(this=self._handle)
    
    @property
    def wvf(self):
        """
        Element wvf ftype=real(dp) pytype=float
        
        
        Defined at Grid_module.fpp line 42
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_eigen_type_r__array__wvf(self._handle)
        if array_handle in self._arrays:
            wvf = self._arrays[array_handle]
        else:
            wvf = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_eigen_type_r__array__wvf)
            self._arrays[array_handle] = wvf
        return wvf
    
    @wvf.setter
    def wvf(self, wvf):
        self.wvf[...] = wvf
    
    @property
    def val(self):
        """
        Element val ftype=real(dp) pytype=float
        
        
        Defined at Grid_module.fpp line 44
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_eigen_type_r__array__val(self._handle)
        if array_handle in self._arrays:
            val = self._arrays[array_handle]
        else:
            val = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_eigen_type_r__array__val)
            self._arrays[array_handle] = val
        return val
    
    @val.setter
    def val(self, val):
        self.val[...] = val
    
    def __str__(self):
        ret = ['<eigen_type_r>{\n']
        ret.append('    wvf : ')
        ret.append(repr(self.wvf))
        ret.append(',\n    val : ')
        ret.append(repr(self.val))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("AresMainPy_pkg.charge_sphere")
class charge_sphere(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=charge_sphere)
    
    
    Defined at Grid_module.fpp lines 48-57
    
    """
    def __init__(self, handle=None):
        """
        self = Charge_Sphere()
        
        
        Defined at Grid_module.fpp lines 48-57
        
        
        Returns
        -------
        this : Charge_Sphere
        	Object to be constructed
        
        
        Automatically generated constructor for charge_sphere
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _AresMainPy_pkg.f90wrap_charge_sphere_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Charge_Sphere
        
        
        Defined at Grid_module.fpp lines 48-57
        
        Parameters
        ----------
        this : Charge_Sphere
        	Object to be destructed
        
        
        Automatically generated destructor for charge_sphere
        """
        if self._alloc:
            _AresMainPy_pkg.f90wrap_charge_sphere_finalise(this=self._handle)
    
    @property
    def onedlength(self):
        """
        Element onedlength ftype=integer(i4b) pytype=int
        
        
        Defined at Grid_module.fpp line 49
        
        """
        return _AresMainPy_pkg.f90wrap_charge_sphere__get__onedlength(self._handle)
    
    @onedlength.setter
    def onedlength(self, onedlength):
        _AresMainPy_pkg.f90wrap_charge_sphere__set__onedlength(self._handle, onedlength)
    
    @property
    def volume(self):
        """
        Element volume ftype=real(dp) pytype=float
        
        
        Defined at Grid_module.fpp line 50
        
        """
        return _AresMainPy_pkg.f90wrap_charge_sphere__get__volume(self._handle)
    
    @volume.setter
    def volume(self, volume):
        _AresMainPy_pkg.f90wrap_charge_sphere__set__volume(self._handle, volume)
    
    @property
    def x(self):
        """
        Element x ftype=integer(i4b) pytype=int
        
        
        Defined at Grid_module.fpp line 51
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_charge_sphere__array__x(self._handle)
        if array_handle in self._arrays:
            x = self._arrays[array_handle]
        else:
            x = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_charge_sphere__array__x)
            self._arrays[array_handle] = x
        return x
    
    @x.setter
    def x(self, x):
        self.x[...] = x
    
    @property
    def y(self):
        """
        Element y ftype=integer(i4b) pytype=int
        
        
        Defined at Grid_module.fpp line 51
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_charge_sphere__array__y(self._handle)
        if array_handle in self._arrays:
            y = self._arrays[array_handle]
        else:
            y = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_charge_sphere__array__y)
            self._arrays[array_handle] = y
        return y
    
    @y.setter
    def y(self, y):
        self.y[...] = y
    
    @property
    def z(self):
        """
        Element z ftype=integer(i4b) pytype=int
        
        
        Defined at Grid_module.fpp line 51
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_charge_sphere__array__z(self._handle)
        if array_handle in self._arrays:
            z = self._arrays[array_handle]
        else:
            z = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_charge_sphere__array__z)
            self._arrays[array_handle] = z
        return z
    
    @z.setter
    def z(self, z):
        self.z[...] = z
    
    @property
    def n(self):
        """
        Element n ftype=integer(i4b) pytype=int
        
        
        Defined at Grid_module.fpp line 51
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_charge_sphere__array__n(self._handle)
        if array_handle in self._arrays:
            n = self._arrays[array_handle]
        else:
            n = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_charge_sphere__array__n)
            self._arrays[array_handle] = n
        return n
    
    @n.setter
    def n(self, n):
        self.n[...] = n
    
    @property
    def onedsphere(self):
        """
        Element onedsphere ftype=real(dp) pytype=float
        
        
        Defined at Grid_module.fpp line 52
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_charge_sphere__array__onedsphere(self._handle)
        if array_handle in self._arrays:
            onedsphere = self._arrays[array_handle]
        else:
            onedsphere = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_charge_sphere__array__onedsphere)
            self._arrays[array_handle] = onedsphere
        return onedsphere
    
    @onedsphere.setter
    def onedsphere(self, onedsphere):
        self.onedsphere[...] = onedsphere
    
    @property
    def onedveff(self):
        """
        Element onedveff ftype=real(dp) pytype=float
        
        
        Defined at Grid_module.fpp line 53
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_charge_sphere__array__onedveff(self._handle)
        if array_handle in self._arrays:
            onedveff = self._arrays[array_handle]
        else:
            onedveff = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_charge_sphere__array__onedveff)
            self._arrays[array_handle] = onedveff
        return onedveff
    
    @onedveff.setter
    def onedveff(self, onedveff):
        self.onedveff[...] = onedveff
    
    @property
    def onedwvf(self):
        """
        Element onedwvf ftype=real(dp) pytype=float
        
        
        Defined at Grid_module.fpp line 54
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_charge_sphere__array__onedwvf(self._handle)
        if array_handle in self._arrays:
            onedwvf = self._arrays[array_handle]
        else:
            onedwvf = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_charge_sphere__array__onedwvf)
            self._arrays[array_handle] = onedwvf
        return onedwvf
    
    @onedwvf.setter
    def onedwvf(self, onedwvf):
        self.onedwvf[...] = onedwvf
    
    @property
    def onedval(self):
        """
        Element onedval ftype=real(dp) pytype=float
        
        
        Defined at Grid_module.fpp line 55
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_charge_sphere__array__onedval(self._handle)
        if array_handle in self._arrays:
            onedval = self._arrays[array_handle]
        else:
            onedval = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_charge_sphere__array__onedval)
            self._arrays[array_handle] = onedval
        return onedval
    
    @onedval.setter
    def onedval(self, onedval):
        self.onedval[...] = onedval
    
    @property
    def initmo(self):
        """
        Element initmo ftype=real(dp) pytype=float
        
        
        Defined at Grid_module.fpp line 56
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_charge_sphere__array__initmo(self._handle)
        if array_handle in self._arrays:
            initmo = self._arrays[array_handle]
        else:
            initmo = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_charge_sphere__array__initmo)
            self._arrays[array_handle] = initmo
        return initmo
    
    @initmo.setter
    def initmo(self, initmo):
        self.initmo[...] = initmo
    
    @property
    def vlpp(self):
        """
        Element vlpp ftype=real(dp) pytype=float
        
        
        Defined at Grid_module.fpp line 57
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_charge_sphere__array__vlpp(self._handle)
        if array_handle in self._arrays:
            vlpp = self._arrays[array_handle]
        else:
            vlpp = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_charge_sphere__array__vlpp)
            self._arrays[array_handle] = vlpp
        return vlpp
    
    @vlpp.setter
    def vlpp(self, vlpp):
        self.vlpp[...] = vlpp
    
    def __str__(self):
        ret = ['<charge_sphere>{\n']
        ret.append('    onedlength : ')
        ret.append(repr(self.onedlength))
        ret.append(',\n    volume : ')
        ret.append(repr(self.volume))
        ret.append(',\n    x : ')
        ret.append(repr(self.x))
        ret.append(',\n    y : ')
        ret.append(repr(self.y))
        ret.append(',\n    z : ')
        ret.append(repr(self.z))
        ret.append(',\n    n : ')
        ret.append(repr(self.n))
        ret.append(',\n    onedsphere : ')
        ret.append(repr(self.onedsphere))
        ret.append(',\n    onedveff : ')
        ret.append(repr(self.onedveff))
        ret.append(',\n    onedwvf : ')
        ret.append(repr(self.onedwvf))
        ret.append(',\n    onedval : ')
        ret.append(repr(self.onedval))
        ret.append(',\n    initmo : ')
        ret.append(repr(self.initmo))
        ret.append(',\n    vlpp : ')
        ret.append(repr(self.vlpp))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("AresMainPy_pkg.grid_diff_map_type")
class grid_diff_map_type(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=grid_diff_map_type)
    
    
    Defined at Grid_module.fpp lines 59-65
    
    """
    def __init__(self, handle=None):
        """
        self = Grid_Diff_Map_Type()
        
        
        Defined at Grid_module.fpp lines 59-65
        
        
        Returns
        -------
        this : Grid_Diff_Map_Type
        	Object to be constructed
        
        
        Automatically generated constructor for grid_diff_map_type
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _AresMainPy_pkg.f90wrap_grid_diff_map_type_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Grid_Diff_Map_Type
        
        
        Defined at Grid_module.fpp lines 59-65
        
        Parameters
        ----------
        this : Grid_Diff_Map_Type
        	Object to be destructed
        
        
        Automatically generated destructor for grid_diff_map_type
        """
        if self._alloc:
            _AresMainPy_pkg.f90wrap_grid_diff_map_type_finalise(this=self._handle)
    
    @property
    def nz_map(self):
        """
        Element nz_map ftype=integer(i4b) pytype=int
        
        
        Defined at Grid_module.fpp line 60
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_grid_diff_map_type__array__nz_map(self._handle)
        if array_handle in self._arrays:
            nz_map = self._arrays[array_handle]
        else:
            nz_map = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_grid_diff_map_type__array__nz_map)
            self._arrays[array_handle] = nz_map
        return nz_map
    
    @nz_map.setter
    def nz_map(self, nz_map):
        self.nz_map[...] = nz_map
    
    @property
    def mycomm_cores(self):
        """
        Element mycomm_cores ftype=integer(i4b) pytype=int
        
        
        Defined at Grid_module.fpp line 61
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_grid_diff_map_type__array__mycomm_cores(self._handle)
        if array_handle in self._arrays:
            mycomm_cores = self._arrays[array_handle]
        else:
            mycomm_cores = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_grid_diff_map_type__array__mycomm_cores)
            self._arrays[array_handle] = mycomm_cores
        return mycomm_cores
    
    @mycomm_cores.setter
    def mycomm_cores(self, mycomm_cores):
        self.mycomm_cores[...] = mycomm_cores
    
    @property
    def mycomm_size(self):
        """
        Element mycomm_size ftype=integer(i4b) pytype=int
        
        
        Defined at Grid_module.fpp line 62
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_grid_diff_map_type__array__mycomm_size(self._handle)
        if array_handle in self._arrays:
            mycomm_size = self._arrays[array_handle]
        else:
            mycomm_size = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_grid_diff_map_type__array__mycomm_size)
            self._arrays[array_handle] = mycomm_size
        return mycomm_size
    
    @mycomm_size.setter
    def mycomm_size(self, mycomm_size):
        self.mycomm_size[...] = mycomm_size
    
    @property
    def mysend_size(self):
        """
        Element mysend_size ftype=integer(i4b) pytype=int
        
        
        Defined at Grid_module.fpp line 63
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_grid_diff_map_type__array__mysend_size(self._handle)
        if array_handle in self._arrays:
            mysend_size = self._arrays[array_handle]
        else:
            mysend_size = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_grid_diff_map_type__array__mysend_size)
            self._arrays[array_handle] = mysend_size
        return mysend_size
    
    @mysend_size.setter
    def mysend_size(self, mysend_size):
        self.mysend_size[...] = mysend_size
    
    @property
    def local_map(self):
        """
        Element local_map ftype=integer(i4b) pytype=int
        
        
        Defined at Grid_module.fpp line 64
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_grid_diff_map_type__array__local_map(self._handle)
        if array_handle in self._arrays:
            local_map = self._arrays[array_handle]
        else:
            local_map = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_grid_diff_map_type__array__local_map)
            self._arrays[array_handle] = local_map
        return local_map
    
    @local_map.setter
    def local_map(self, local_map):
        self.local_map[...] = local_map
    
    @property
    def boundary(self):
        """
        Element boundary ftype=integer(i4b) pytype=int
        
        
        Defined at Grid_module.fpp line 65
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_grid_diff_map_type__array__boundary(self._handle)
        if array_handle in self._arrays:
            boundary = self._arrays[array_handle]
        else:
            boundary = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_grid_diff_map_type__array__boundary)
            self._arrays[array_handle] = boundary
        return boundary
    
    @boundary.setter
    def boundary(self, boundary):
        self.boundary[...] = boundary
    
    def __str__(self):
        ret = ['<grid_diff_map_type>{\n']
        ret.append('    nz_map : ')
        ret.append(repr(self.nz_map))
        ret.append(',\n    mycomm_cores : ')
        ret.append(repr(self.mycomm_cores))
        ret.append(',\n    mycomm_size : ')
        ret.append(repr(self.mycomm_size))
        ret.append(',\n    mysend_size : ')
        ret.append(repr(self.mysend_size))
        ret.append(',\n    local_map : ')
        ret.append(repr(self.local_map))
        ret.append(',\n    boundary : ')
        ret.append(repr(self.boundary))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def build_rgrid():
    """
    build_rgrid()
    
    
    Defined at Grid_module.fpp lines 102-165
    
    
    =============================================================
    isolate system and odd grid
     IF(.NOT.Lpbc.AND.MOD(n1,2)==0 )n1=n1+1
     IF(.NOT.Lpbc.AND.MOD(n2,2)==0 )n2=n2+1
     IF(.NOT.Lpbc.AND.MOD(n3,2)==0 )n3=n3+1
    =============================================================
    grid size
    """
    _AresMainPy_pkg.f90wrap_build_rgrid()

def build_rgrid_iso():
    """
    build_rgrid_iso()
    
    
    Defined at Grid_module.fpp lines 168-219
    
    
    =============================================================
    isolate system and odd grid
    """
    _AresMainPy_pkg.f90wrap_build_rgrid_iso()

def destroy_rgrid():
    """
    destroy_rgrid()
    
    
    Defined at Grid_module.fpp lines 222-251
    
    
    """
    _AresMainPy_pkg.f90wrap_destroy_rgrid()

def build_kgrid():
    """
    build_kgrid()
    
    
    Defined at Grid_module.fpp lines 254-353
    
    
    """
    _AresMainPy_pkg.f90wrap_build_kgrid()

def destroy_kpt():
    """
    destroy_kpt()
    
    
    Defined at Grid_module.fpp lines 356-372
    
    
    """
    _AresMainPy_pkg.f90wrap_destroy_kpt()

def build_eigen():
    """
    build_eigen()
    
    
    Defined at Grid_module.fpp lines 375-419
    
    
    """
    _AresMainPy_pkg.f90wrap_build_eigen()

def destroy_eigen():
    """
    destroy_eigen()
    
    
    Defined at Grid_module.fpp lines 422-444
    
    
    """
    _AresMainPy_pkg.f90wrap_destroy_eigen()

def fillqtable():
    """
    fillqtable()
    
    
    Defined at Grid_module.fpp lines 447-505
    
    
    """
    _AresMainPy_pkg.f90wrap_fillqtable()

def fillrtable():
    """
    fillrtable()
    
    
    Defined at Grid_module.fpp lines 508-589
    
    
    """
    _AresMainPy_pkg.f90wrap_fillrtable()

def fillrtable_iso():
    """
    fillrtable_iso()
    
    
    Defined at Grid_module.fpp lines 592-620
    
    
    """
    _AresMainPy_pkg.f90wrap_fillrtable_iso()

def build_iso_sphere_grid():
    """
    build_iso_sphere_grid()
    
    
    Defined at Grid_module.fpp lines 623-749
    
    
    ======================================
    ##XLT OBTAIN GRID INFORMATION
    print*,"orig",orig
    print*,"struct%poscar",struct%poscar
    ======================================
    """
    _AresMainPy_pkg.f90wrap_build_iso_sphere_grid()

def matchposcar_iso():
    """
    matchposcar_iso()
    
    
    Defined at Grid_module.fpp lines 752-766
    
    
    =========================================================
    ##move atoms in line with the grid coordinate system
    """
    _AresMainPy_pkg.f90wrap_matchposcar_iso()

def set_car2spe(orig, x, y, z):
    """
    r, cost, sint, cosp, sinp = set_car2spe(orig, x, y, z)
    
    
    Defined at Grid_module.fpp lines 771-802
    
    Parameters
    ----------
    orig : float array
    x : int
    y : int
    z : int
    
    Returns
    -------
    r : float
    cost : float
    sint : float
    cosp : float
    sinp : float
    
    """
    r, cost, sint, cosp, sinp = _AresMainPy_pkg.f90wrap_set_car2spe(orig=orig, x=x, \
        y=y, z=z)
    return r, cost, sint, cosp, sinp

def destroy_iso_sphere_grid():
    """
    destroy_iso_sphere_grid()
    
    
    Defined at Grid_module.fpp lines 806-857
    
    
    """
    _AresMainPy_pkg.f90wrap_destroy_iso_sphere_grid()

def iso_vsphere(vgrid, vsphere):
    """
    iso_vsphere(vgrid, vsphere)
    
    
    Defined at Grid_module.fpp lines 861-872
    
    Parameters
    ----------
    vgrid : float array
    vsphere : float array
    
    """
    _AresMainPy_pkg.f90wrap_iso_vsphere(vgrid=vgrid, vsphere=vsphere)

def iso_rho2grid(rhosphere, rhogrid):
    """
    iso_rho2grid(rhosphere, rhogrid)
    
    
    Defined at Grid_module.fpp lines 876-888
    
    Parameters
    ----------
    rhosphere : float array
    rhogrid : float array
    
    """
    _AresMainPy_pkg.f90wrap_iso_rho2grid(rhosphere=rhosphere, rhogrid=rhogrid)

def parallel_s2g(p, thrq):
    """
    parallel_s2g(p, thrq)
    
    
    Defined at Grid_module.fpp lines 891-914
    
    Parameters
    ----------
    p : float array
    thrq : float array
    
    """
    _AresMainPy_pkg.f90wrap_parallel_s2g(p=p, thrq=thrq)

def parallel_g2s(p, thrq):
    """
    parallel_g2s(p, thrq)
    
    
    Defined at Grid_module.fpp lines 917-932
    
    Parameters
    ----------
    p : float array
    thrq : float array
    
    """
    _AresMainPy_pkg.f90wrap_parallel_g2s(p=p, thrq=thrq)

def get_z_range(my_z_range, delta_z, ngrid_z, cell_thick, cell_shape):
    """
    get_z_range(my_z_range, delta_z, ngrid_z, cell_thick, cell_shape)
    
    
    Defined at Grid_module.fpp lines 935-953
    
    Parameters
    ----------
    my_z_range : int array
    delta_z : float
    ngrid_z : int
    cell_thick : float
    cell_shape : int
    
    """
    _AresMainPy_pkg.f90wrap_get_z_range(my_z_range=my_z_range, delta_z=delta_z, \
        ngrid_z=ngrid_z, cell_thick=cell_thick, cell_shape=cell_shape)

def confirm_iso_radius():
    """
    confirm_iso_radius()
    
    
    Defined at Grid_module.fpp lines 956-1006
    
    
    ========================================
    > find MAX(r)
    """
    _AresMainPy_pkg.f90wrap_confirm_iso_radius()

def reshape_center():
    """
    reshape_center()
    
    
    Defined at Grid_module.fpp lines 1009-1034
    
    
    -----------------------------------------------
    > reset the direct pos
    """
    _AresMainPy_pkg.f90wrap_reshape_center()

def reset_poscar(new_gap, new_ni):
    """
    reset_poscar(new_gap, new_ni)
    
    
    Defined at Grid_module.fpp lines 1038-1071
    
    Parameters
    ----------
    new_gap : float
    new_ni : int
    
    -----------------------------------------------
    > reset the direct pos
    """
    _AresMainPy_pkg.f90wrap_reset_poscar(new_gap=new_gap, new_ni=new_ni)

def build_parallel_cubic_grid():
    """
    build_parallel_cubic_grid()
    
    
    Defined at Grid_module.fpp lines 1074-1111
    
    
    """
    _AresMainPy_pkg.f90wrap_build_parallel_cubic_grid()

def rho_trans1d(rho3d, rho1d):
    """
    rho_trans1d(rho3d, rho1d)
    
    
    Defined at Grid_module.fpp lines 1113-1157
    
    Parameters
    ----------
    rho3d : float array
    rho1d : float array
    
    """
    _AresMainPy_pkg.f90wrap_rho_trans1d(rho3d=rho3d, rho1d=rho1d)

def rho_trans3d(rho1d, rho3d):
    """
    rho_trans3d(rho1d, rho3d)
    
    
    Defined at Grid_module.fpp lines 1159-1203
    
    Parameters
    ----------
    rho1d : float array
    rho3d : float array
    
    """
    _AresMainPy_pkg.f90wrap_rho_trans3d(rho1d=rho1d, rho3d=rho3d)

def _fft_sph_r2c(array_r):
    """
    array_c = _fft_sph_r2c(array_r)
    
    
    Defined at Grid_module.fpp lines 1206-1224
    
    Parameters
    ----------
    array_r : float array
    
    Returns
    -------
    array_c : complex array
    
    """
    array_c = _AresMainPy_pkg.f90wrap_fft_sph_r2c(array_r=array_r)
    return array_c

def _fft_sph_c2r(array_c):
    """
    array_r = _fft_sph_c2r(array_c)
    
    
    Defined at Grid_module.fpp lines 1227-1246
    
    Parameters
    ----------
    array_c : complex array
    
    Returns
    -------
    array_r : float array
    
    """
    array_r = _AresMainPy_pkg.f90wrap_fft_sph_c2r(array_c=array_c)
    return array_r

def fft_sph(*args, **kwargs):
    """
    fft_sph(*args, **kwargs)
    
    
    Defined at Grid_module.fpp lines 96-98
    
    Overloaded interface containing the following procedures:
      _fft_sph_r2c
      _fft_sph_c2r
    
    """
    for proc in [_fft_sph_r2c, _fft_sph_c2r]:
        try:
            return proc(*args, **kwargs)
        except TypeError:
            continue

def get_n1():
    """
    Element n1 ftype=integer(i4b) pytype=int
    
    
    Defined at Grid_module.fpp line 71
    
    """
    return _AresMainPy_pkg.f90wrap_grid_module__get__n1()

def set_n1(n1):
    _AresMainPy_pkg.f90wrap_grid_module__set__n1(n1)

def get_n2():
    """
    Element n2 ftype=integer(i4b) pytype=int
    
    
    Defined at Grid_module.fpp line 71
    
    """
    return _AresMainPy_pkg.f90wrap_grid_module__get__n2()

def set_n2(n2):
    _AresMainPy_pkg.f90wrap_grid_module__set__n2(n2)

def get_n3():
    """
    Element n3 ftype=integer(i4b) pytype=int
    
    
    Defined at Grid_module.fpp line 71
    
    """
    return _AresMainPy_pkg.f90wrap_grid_module__get__n3()

def set_n3(n3):
    _AresMainPy_pkg.f90wrap_grid_module__set__n3(n3)

def get_n():
    """
    Element n ftype=integer(i4b) pytype=int
    
    
    Defined at Grid_module.fpp line 71
    
    """
    return _AresMainPy_pkg.f90wrap_grid_module__get__n()

def set_n(n):
    _AresMainPy_pkg.f90wrap_grid_module__set__n(n)

def get_nsn():
    """
    Element nsn ftype=integer(i4b) pytype=int
    
    
    Defined at Grid_module.fpp line 71
    
    """
    return _AresMainPy_pkg.f90wrap_grid_module__get__nsn()

def set_nsn(nsn):
    _AresMainPy_pkg.f90wrap_grid_module__set__nsn(nsn)

def get_global_n1():
    """
    Element global_n1 ftype=integer(i4b) pytype=int
    
    
    Defined at Grid_module.fpp line 72
    
    """
    return _AresMainPy_pkg.f90wrap_grid_module__get__global_n1()

def set_global_n1(global_n1):
    _AresMainPy_pkg.f90wrap_grid_module__set__global_n1(global_n1)

def get_global_n2():
    """
    Element global_n2 ftype=integer(i4b) pytype=int
    
    
    Defined at Grid_module.fpp line 72
    
    """
    return _AresMainPy_pkg.f90wrap_grid_module__get__global_n2()

def set_global_n2(global_n2):
    _AresMainPy_pkg.f90wrap_grid_module__set__global_n2(global_n2)

def get_global_n3():
    """
    Element global_n3 ftype=integer(i4b) pytype=int
    
    
    Defined at Grid_module.fpp line 72
    
    """
    return _AresMainPy_pkg.f90wrap_grid_module__get__global_n3()

def set_global_n3(global_n3):
    _AresMainPy_pkg.f90wrap_grid_module__set__global_n3(global_n3)

def get_global_n():
    """
    Element global_n ftype=integer(i4b) pytype=int
    
    
    Defined at Grid_module.fpp line 72
    
    """
    return _AresMainPy_pkg.f90wrap_grid_module__get__global_n()

def set_global_n(global_n):
    _AresMainPy_pkg.f90wrap_grid_module__set__global_n(global_n)

def get_ng1():
    """
    Element ng1 ftype=integer(i4b) pytype=int
    
    
    Defined at Grid_module.fpp line 73
    
    """
    return _AresMainPy_pkg.f90wrap_grid_module__get__ng1()

def set_ng1(ng1):
    _AresMainPy_pkg.f90wrap_grid_module__set__ng1(ng1)

def get_ng2():
    """
    Element ng2 ftype=integer(i4b) pytype=int
    
    
    Defined at Grid_module.fpp line 73
    
    """
    return _AresMainPy_pkg.f90wrap_grid_module__get__ng2()

def set_ng2(ng2):
    _AresMainPy_pkg.f90wrap_grid_module__set__ng2(ng2)

def get_ng3():
    """
    Element ng3 ftype=integer(i4b) pytype=int
    
    
    Defined at Grid_module.fpp line 73
    
    """
    return _AresMainPy_pkg.f90wrap_grid_module__get__ng3()

def set_ng3(ng3):
    _AresMainPy_pkg.f90wrap_grid_module__set__ng3(ng3)

def get_ng():
    """
    Element ng ftype=integer(i4b) pytype=int
    
    
    Defined at Grid_module.fpp line 73
    
    """
    return _AresMainPy_pkg.f90wrap_grid_module__get__ng()

def set_ng(ng):
    _AresMainPy_pkg.f90wrap_grid_module__set__ng(ng)

def get_array_gap():
    """
    Element gap ftype=real(dp) pytype=float
    
    
    Defined at Grid_module.fpp line 74
    
    """
    global gap
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_grid_module__array__gap(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        gap = _arrays[array_handle]
    else:
        gap = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_grid_module__array__gap)
        _arrays[array_handle] = gap
    return gap

def set_array_gap(gap):
    gap[...] = gap

def get_dvol():
    """
    Element dvol ftype=real(dp) pytype=float
    
    
    Defined at Grid_module.fpp line 75
    
    """
    return _AresMainPy_pkg.f90wrap_grid_module__get__dvol()

def set_dvol(dvol):
    _AresMainPy_pkg.f90wrap_grid_module__set__dvol(dvol)

def get_nk1():
    """
    Element nk1 ftype=integer(i4b) pytype=int
    
    
    Defined at Grid_module.fpp line 77
    
    """
    return _AresMainPy_pkg.f90wrap_grid_module__get__nk1()

def set_nk1(nk1):
    _AresMainPy_pkg.f90wrap_grid_module__set__nk1(nk1)

def get_nk2():
    """
    Element nk2 ftype=integer(i4b) pytype=int
    
    
    Defined at Grid_module.fpp line 77
    
    """
    return _AresMainPy_pkg.f90wrap_grid_module__get__nk2()

def set_nk2(nk2):
    _AresMainPy_pkg.f90wrap_grid_module__set__nk2(nk2)

def get_nk3():
    """
    Element nk3 ftype=integer(i4b) pytype=int
    
    
    Defined at Grid_module.fpp line 77
    
    """
    return _AresMainPy_pkg.f90wrap_grid_module__get__nk3()

def set_nk3(nk3):
    _AresMainPy_pkg.f90wrap_grid_module__set__nk3(nk3)

def get_nk():
    """
    Element nk ftype=integer(i4b) pytype=int
    
    
    Defined at Grid_module.fpp line 77
    
    """
    return _AresMainPy_pkg.f90wrap_grid_module__get__nk()

def set_nk(nk):
    _AresMainPy_pkg.f90wrap_grid_module__set__nk(nk)

def get_array_kdispl():
    """
    Element kdispl ftype=real(dp) pytype=float
    
    
    Defined at Grid_module.fpp line 78
    
    """
    global kdispl
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_grid_module__array__kdispl(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        kdispl = _arrays[array_handle]
    else:
        kdispl = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_grid_module__array__kdispl)
        _arrays[array_handle] = kdispl
    return kdispl

def set_array_kdispl(kdispl):
    kdispl[...] = kdispl

def get_array_dr():
    """
    Element dr ftype=real(dp) pytype=float
    
    
    Defined at Grid_module.fpp line 88
    
    """
    global dr
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_grid_module__array__dr(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        dr = _arrays[array_handle]
    else:
        dr = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_grid_module__array__dr)
        _arrays[array_handle] = dr
    return dr

def set_array_dr(dr):
    dr[...] = dr

def get_array_cost():
    """
    Element cost ftype=real(dp) pytype=float
    
    
    Defined at Grid_module.fpp line 89
    
    """
    global cost
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_grid_module__array__cost(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        cost = _arrays[array_handle]
    else:
        cost = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_grid_module__array__cost)
        _arrays[array_handle] = cost
    return cost

def set_array_cost(cost):
    cost[...] = cost

def get_array_sint():
    """
    Element sint ftype=real(dp) pytype=float
    
    
    Defined at Grid_module.fpp line 90
    
    """
    global sint
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_grid_module__array__sint(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        sint = _arrays[array_handle]
    else:
        sint = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_grid_module__array__sint)
        _arrays[array_handle] = sint
    return sint

def set_array_sint(sint):
    sint[...] = sint

def get_array_cns():
    """
    Element cns ftype=complex(dcp) pytype=complex
    
    
    Defined at Grid_module.fpp line 91
    
    """
    global cns
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_grid_module__array__cns(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        cns = _arrays[array_handle]
    else:
        cns = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_grid_module__array__cns)
        _arrays[array_handle] = cns
    return cns

def set_array_cns(cns):
    cns[...] = cns

def get_array_indx():
    """
    Element indx ftype=integer(i4b) pytype=int
    
    
    Defined at Grid_module.fpp line 92
    
    """
    global indx
    array_ndim, array_type, array_shape, array_handle = \
        _AresMainPy_pkg.f90wrap_grid_module__array__indx(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        indx = _arrays[array_handle]
    else:
        indx = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _AresMainPy_pkg.f90wrap_grid_module__array__indx)
        _arrays[array_handle] = indx
    return indx

def set_array_indx(indx):
    indx[...] = indx


_array_initialisers = [get_array_gap, get_array_kdispl, get_array_dr, \
    get_array_cost, get_array_sint, get_array_cns, get_array_indx]
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "grid_module".')

for func in _dt_array_initialisers:
    func()

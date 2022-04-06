"""
Module grid_module


Defined at Struct_module.fpp lines 1054-2105

"""
from __future__ import print_function, absolute_import, division
import _arespy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("arespy_pkg.grid_type")
class grid_type(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=grid_type)
    
    
    Defined at Struct_module.fpp lines 1063-1082
    
    """
    def __init__(self, handle=None):
        """
        self = Grid_Type()
        
        
        Defined at Struct_module.fpp lines 1063-1082
        
        
        Returns
        -------
        this : Grid_Type
        	Object to be constructed
        
        
        Automatically generated constructor for grid_type
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _arespy_pkg.f90wrap_grid_type_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Grid_Type
        
        
        Defined at Struct_module.fpp lines 1063-1082
        
        Parameters
        ----------
        this : Grid_Type
        	Object to be destructed
        
        
        Automatically generated destructor for grid_type
        """
        if self._alloc:
            _arespy_pkg.f90wrap_grid_type_finalise(this=self._handle)
    
    @property
    def rhos(self):
        """
        Element rhos ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 1065
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_grid_type__array__rhos(self._handle)
        if array_handle in self._arrays:
            rhos = self._arrays[array_handle]
        else:
            rhos = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_grid_type__array__rhos)
            self._arrays[array_handle] = rhos
        return rhos
    
    @rhos.setter
    def rhos(self, rhos):
        self.rhos[...] = rhos
    
    @property
    def rho(self):
        """
        Element rho ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 1066
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_grid_type__array__rho(self._handle)
        if array_handle in self._arrays:
            rho = self._arrays[array_handle]
        else:
            rho = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_grid_type__array__rho)
            self._arrays[array_handle] = rho
        return rho
    
    @rho.setter
    def rho(self, rho):
        self.rho[...] = rho
    
    @property
    def vxcs(self):
        """
        Element vxcs ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 1067
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_grid_type__array__vxcs(self._handle)
        if array_handle in self._arrays:
            vxcs = self._arrays[array_handle]
        else:
            vxcs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_grid_type__array__vxcs)
            self._arrays[array_handle] = vxcs
        return vxcs
    
    @vxcs.setter
    def vxcs(self, vxcs):
        self.vxcs[...] = vxcs
    
    @property
    def vhxcd(self):
        """
        Element vhxcd ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 1068
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_grid_type__array__vhxcd(self._handle)
        if array_handle in self._arrays:
            vhxcd = self._arrays[array_handle]
        else:
            vhxcd = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_grid_type__array__vhxcd)
            self._arrays[array_handle] = vhxcd
        return vhxcd
    
    @vhxcd.setter
    def vhxcd(self, vhxcd):
        self.vhxcd[...] = vhxcd
    
    @property
    def vlpp(self):
        """
        Element vlpp ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 1069
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_grid_type__array__vlpp(self._handle)
        if array_handle in self._arrays:
            vlpp = self._arrays[array_handle]
        else:
            vlpp = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_grid_type__array__vlpp)
            self._arrays[array_handle] = vlpp
        return vlpp
    
    @vlpp.setter
    def vlpp(self, vlpp):
        self.vlpp[...] = vlpp
    
    @property
    def vh(self):
        """
        Element vh ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 1070
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_grid_type__array__vh(self._handle)
        if array_handle in self._arrays:
            vh = self._arrays[array_handle]
        else:
            vh = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_grid_type__array__vh)
            self._arrays[array_handle] = vh
        return vh
    
    @vh.setter
    def vh(self, vh):
        self.vh[...] = vh
    
    @property
    def eval(self):
        """
        Element eval ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 1071
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_grid_type__array__eval(self._handle)
        if array_handle in self._arrays:
            eval = self._arrays[array_handle]
        else:
            eval = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_grid_type__array__eval)
            self._arrays[array_handle] = eval
        return eval
    
    @eval.setter
    def eval(self, eval):
        self.eval[...] = eval
    
    @property
    def rhoc(self):
        """
        Element rhoc ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 1073
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_grid_type__array__rhoc(self._handle)
        if array_handle in self._arrays:
            rhoc = self._arrays[array_handle]
        else:
            rhoc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_grid_type__array__rhoc)
            self._arrays[array_handle] = rhoc
        return rhoc
    
    @rhoc.setter
    def rhoc(self, rhoc):
        self.rhoc[...] = rhoc
    
    @property
    def gvec(self):
        """
        Element gvec ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 1075
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_grid_type__array__gvec(self._handle)
        if array_handle in self._arrays:
            gvec = self._arrays[array_handle]
        else:
            gvec = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_grid_type__array__gvec)
            self._arrays[array_handle] = gvec
        return gvec
    
    @gvec.setter
    def gvec(self, gvec):
        self.gvec[...] = gvec
    
    @property
    def gmask(self):
        """
        Element gmask ftype=logical pytype=bool
        
        
        Defined at Struct_module.fpp line 1076
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_grid_type__array__gmask(self._handle)
        if array_handle in self._arrays:
            gmask = self._arrays[array_handle]
        else:
            gmask = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_grid_type__array__gmask)
            self._arrays[array_handle] = gmask
        return gmask
    
    @gmask.setter
    def gmask(self, gmask):
        self.gmask[...] = gmask
    
    @property
    def rvec(self):
        """
        Element rvec ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 1078
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_grid_type__array__rvec(self._handle)
        if array_handle in self._arrays:
            rvec = self._arrays[array_handle]
        else:
            rvec = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_grid_type__array__rvec)
            self._arrays[array_handle] = rvec
        return rvec
    
    @rvec.setter
    def rvec(self, rvec):
        self.rvec[...] = rvec
    
    @property
    def lsp(self):
        """
        Element lsp ftype=logical pytype=bool
        
        
        Defined at Struct_module.fpp line 1080
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_grid_type__array__lsp(self._handle)
        if array_handle in self._arrays:
            lsp = self._arrays[array_handle]
        else:
            lsp = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_grid_type__array__lsp)
            self._arrays[array_handle] = lsp
        return lsp
    
    @lsp.setter
    def lsp(self, lsp):
        self.lsp[...] = lsp
    
    @property
    def onedlength(self):
        """
        Element onedlength ftype=integer(i4b) pytype=int
        
        
        Defined at Struct_module.fpp line 1082
        
        """
        return _arespy_pkg.f90wrap_grid_type__get__onedlength(self._handle)
    
    @onedlength.setter
    def onedlength(self, onedlength):
        _arespy_pkg.f90wrap_grid_type__set__onedlength(self._handle, onedlength)
    
    def __str__(self):
        ret = ['<grid_type>{\n']
        ret.append('    rhos : ')
        ret.append(repr(self.rhos))
        ret.append(',\n    rho : ')
        ret.append(repr(self.rho))
        ret.append(',\n    vxcs : ')
        ret.append(repr(self.vxcs))
        ret.append(',\n    vhxcd : ')
        ret.append(repr(self.vhxcd))
        ret.append(',\n    vlpp : ')
        ret.append(repr(self.vlpp))
        ret.append(',\n    vh : ')
        ret.append(repr(self.vh))
        ret.append(',\n    eval : ')
        ret.append(repr(self.eval))
        ret.append(',\n    rhoc : ')
        ret.append(repr(self.rhoc))
        ret.append(',\n    gvec : ')
        ret.append(repr(self.gvec))
        ret.append(',\n    gmask : ')
        ret.append(repr(self.gmask))
        ret.append(',\n    rvec : ')
        ret.append(repr(self.rvec))
        ret.append(',\n    lsp : ')
        ret.append(repr(self.lsp))
        ret.append(',\n    onedlength : ')
        ret.append(repr(self.onedlength))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("arespy_pkg.kgrid_type")
class kgrid_type(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=kgrid_type)
    
    
    Defined at Struct_module.fpp lines 1085-1089
    
    """
    def __init__(self, handle=None):
        """
        self = Kgrid_Type()
        
        
        Defined at Struct_module.fpp lines 1085-1089
        
        
        Returns
        -------
        this : Kgrid_Type
        	Object to be constructed
        
        
        Automatically generated constructor for kgrid_type
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _arespy_pkg.f90wrap_kgrid_type_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Kgrid_Type
        
        
        Defined at Struct_module.fpp lines 1085-1089
        
        Parameters
        ----------
        this : Kgrid_Type
        	Object to be destructed
        
        
        Automatically generated destructor for kgrid_type
        """
        if self._alloc:
            _arespy_pkg.f90wrap_kgrid_type_finalise(this=self._handle)
    
    @property
    def vec(self):
        """
        Element vec ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 1087
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_kgrid_type__array__vec(self._handle)
        if array_handle in self._arrays:
            vec = self._arrays[array_handle]
        else:
            vec = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_kgrid_type__array__vec)
            self._arrays[array_handle] = vec
        return vec
    
    @vec.setter
    def vec(self, vec):
        self.vec[...] = vec
    
    @property
    def vcar(self):
        """
        Element vcar ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 1088
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_kgrid_type__array__vcar(self._handle)
        if array_handle in self._arrays:
            vcar = self._arrays[array_handle]
        else:
            vcar = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_kgrid_type__array__vcar)
            self._arrays[array_handle] = vcar
        return vcar
    
    @vcar.setter
    def vcar(self, vcar):
        self.vcar[...] = vcar
    
    @property
    def wk(self):
        """
        Element wk ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 1089
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_kgrid_type__array__wk(self._handle)
        if array_handle in self._arrays:
            wk = self._arrays[array_handle]
        else:
            wk = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_kgrid_type__array__wk)
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
    

@f90wrap.runtime.register_class("arespy_pkg.eigen_type")
class eigen_type(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=eigen_type)
    
    
    Defined at Struct_module.fpp lines 1092-1097
    
    """
    def __init__(self, handle=None):
        """
        self = Eigen_Type()
        
        
        Defined at Struct_module.fpp lines 1092-1097
        
        
        Returns
        -------
        this : Eigen_Type
        	Object to be constructed
        
        
        Automatically generated constructor for eigen_type
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _arespy_pkg.f90wrap_eigen_type_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Eigen_Type
        
        
        Defined at Struct_module.fpp lines 1092-1097
        
        Parameters
        ----------
        this : Eigen_Type
        	Object to be destructed
        
        
        Automatically generated destructor for eigen_type
        """
        if self._alloc:
            _arespy_pkg.f90wrap_eigen_type_finalise(this=self._handle)
    
    @property
    def wvf(self):
        """
        Element wvf ftype=complex(dcp) pytype=complex
        
        
        Defined at Struct_module.fpp line 1094
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_eigen_type__array__wvf(self._handle)
        if array_handle in self._arrays:
            wvf = self._arrays[array_handle]
        else:
            wvf = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_eigen_type__array__wvf)
            self._arrays[array_handle] = wvf
        return wvf
    
    @wvf.setter
    def wvf(self, wvf):
        self.wvf[...] = wvf
    
    @property
    def val(self):
        """
        Element val ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 1096
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_eigen_type__array__val(self._handle)
        if array_handle in self._arrays:
            val = self._arrays[array_handle]
        else:
            val = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_eigen_type__array__val)
            self._arrays[array_handle] = val
        return val
    
    @val.setter
    def val(self, val):
        self.val[...] = val
    
    @property
    def wvfg(self):
        """
        Element wvfg ftype=real(dp) pytype=float
        
        
        Defined at Struct_module.fpp line 1097
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_eigen_type__array__wvfg(self._handle)
        if array_handle in self._arrays:
            wvfg = self._arrays[array_handle]
        else:
            wvfg = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_eigen_type__array__wvfg)
            self._arrays[array_handle] = wvfg
        return wvfg
    
    @wvfg.setter
    def wvfg(self, wvfg):
        self.wvfg[...] = wvfg
    
    def __str__(self):
        ret = ['<eigen_type>{\n']
        ret.append('    wvf : ')
        ret.append(repr(self.wvf))
        ret.append(',\n    val : ')
        ret.append(repr(self.val))
        ret.append(',\n    wvfg : ')
        ret.append(repr(self.wvfg))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def build_rgrid():
    """
    build_rgrid()
    
    
    Defined at Struct_module.fpp lines 1121-1182
    
    
    """
    _arespy_pkg.f90wrap_build_rgrid()

def destroy_rgrid():
    """
    destroy_rgrid()
    
    
    Defined at Struct_module.fpp lines 1185-1202
    
    
    """
    _arespy_pkg.f90wrap_destroy_rgrid()

def build_kgrid():
    """
    build_kgrid()
    
    
    Defined at Struct_module.fpp lines 1205-1237
    
    
    """
    _arespy_pkg.f90wrap_build_kgrid()

def destroy_kpt():
    """
    destroy_kpt()
    
    
    Defined at Struct_module.fpp lines 1240-1247
    
    
    """
    _arespy_pkg.f90wrap_destroy_kpt()

def build_eigen():
    """
    build_eigen()
    
    
    Defined at Struct_module.fpp lines 1250-1269
    
    
    """
    _arespy_pkg.f90wrap_build_eigen()

def destroy_eigen():
    """
    destroy_eigen()
    
    
    Defined at Struct_module.fpp lines 1272-1278
    
    
    """
    _arespy_pkg.f90wrap_destroy_eigen()

def fillqtable():
    """
    fillqtable()
    
    
    Defined at Struct_module.fpp lines 1281-1339
    
    
    """
    _arespy_pkg.f90wrap_fillqtable()

def fillrtable():
    """
    fillrtable()
    
    
    Defined at Struct_module.fpp lines 1342-1372
    
    
    """
    _arespy_pkg.f90wrap_fillrtable()

def symm_kgrid():
    """
    symm_kgrid()
    
    
    Defined at Struct_module.fpp lines 1375-1549
    
    
    -----------------------
    """
    _arespy_pkg.f90wrap_symm_kgrid()

def symm_density(rho):
    """
    symm_density(rho)
    
    
    Defined at Struct_module.fpp lines 1552-1629
    
    Parameters
    ----------
    rho : float array
    
    """
    _arespy_pkg.f90wrap_symm_density(rho=rho)

def isymm_apply(isy, nsize, vin, vout):
    """
    lsymm = isymm_apply(isy, nsize, vin, vout)
    
    
    Defined at Struct_module.fpp lines 1632-1691
    
    Parameters
    ----------
    isy : int
    nsize : int array
    vin : int array
    vout : int array
    
    Returns
    -------
    lsymm : bool
    
    """
    lsymm = _arespy_pkg.f90wrap_isymm_apply(isy=isy, nsize=nsize, vin=vin, \
        vout=vout)
    return lsymm

def build_parallel_3d_grid():
    """
    build_parallel_3d_grid()
    
    
    Defined at Struct_module.fpp lines 1694-1714
    
    
    """
    _arespy_pkg.f90wrap_build_parallel_3d_grid()

def build_parallel_sph_grid():
    """
    build_parallel_sph_grid()
    
    
    Defined at Struct_module.fpp lines 1717-1730
    
    
    """
    _arespy_pkg.f90wrap_build_parallel_sph_grid()

def grid_split_sph(ngrid, ncore, comm, id, grid_range, recvcounts, displs, \
    gridrange_sum, n1, n2, n3, n):
    """
    grid_split_sph(ngrid, ncore, comm, id, grid_range, recvcounts, displs, \
        gridrange_sum, n1, n2, n3, n)
    
    
    Defined at Struct_module.fpp lines 1733-1832
    
    Parameters
    ----------
    ngrid : int
    ncore : int
    comm : int
    id : int
    grid_range : int array
    recvcounts : int array
    displs : int array
    gridrange_sum : int array
    n1 : int
    n2 : int
    n3 : int
    n : int
    
    """
    _arespy_pkg.f90wrap_grid_split_sph(ngrid=ngrid, ncore=ncore, comm=comm, id=id, \
        grid_range=grid_range, recvcounts=recvcounts, displs=displs, \
        gridrange_sum=gridrange_sum, n1=n1, n2=n2, n3=n3, n=n)

def sphere_region(n1, n2, n3, lsphere, nz_map, nsphere):
    """
    sphere_region(n1, n2, n3, lsphere, nz_map, nsphere)
    
    
    Defined at Struct_module.fpp lines 1835-1924
    
    Parameters
    ----------
    n1 : int
    n2 : int
    n3 : int
    lsphere : bool array
    nz_map : int array
    nsphere : int
    
    """
    _arespy_pkg.f90wrap_sphere_region(n1=n1, n2=n2, n3=n3, lsphere=lsphere, \
        nz_map=nz_map, nsphere=nsphere)

def sphere2cubic(nps, f1d, f3d, rfill=None):
    """
    sphere2cubic(nps, f1d, f3d[, rfill])
    
    
    Defined at Struct_module.fpp lines 1927-1963
    
    Parameters
    ----------
    nps : int
    f1d : float array
    f3d : float array
    rfill : float
    
    """
    _arespy_pkg.f90wrap_sphere2cubic(nps=nps, f1d=f1d, f3d=f3d, rfill=rfill)

def cubic2sphere(nps, f3d, f1d):
    """
    cubic2sphere(nps, f3d, f1d)
    
    
    Defined at Struct_module.fpp lines 1966-1985
    
    Parameters
    ----------
    nps : int
    f3d : float array
    f1d : float array
    
    """
    _arespy_pkg.f90wrap_cubic2sphere(nps=nps, f3d=f3d, f1d=f1d)

def cubic2sphere_fft(nps, f3d, f1d, shiftn, shiftz):
    """
    cubic2sphere_fft(nps, f3d, f1d, shiftn, shiftz)
    
    
    Defined at Struct_module.fpp lines 1988-2006
    
    Parameters
    ----------
    nps : int
    f3d : float array
    f1d : float array
    shiftn : int
    shiftz : int
    
    """
    _arespy_pkg.f90wrap_cubic2sphere_fft(nps=nps, f3d=f3d, f1d=f1d, shiftn=shiftn, \
        shiftz=shiftz)

def sphere2cubic_fft(nps, f1d, f3d, shiftn, shiftz):
    """
    sphere2cubic_fft(nps, f1d, f3d, shiftn, shiftz)
    
    
    Defined at Struct_module.fpp lines 2009-2027
    
    Parameters
    ----------
    nps : int
    f1d : float array
    f3d : float array
    shiftn : int
    shiftz : int
    
    """
    _arespy_pkg.f90wrap_sphere2cubic_fft(nps=nps, f1d=f1d, f3d=f3d, shiftn=shiftn, \
        shiftz=shiftz)

def sumrhos(nps, rhos, rho):
    """
    sumrhos(nps, rhos, rho)
    
    
    Defined at Struct_module.fpp lines 2087-2104
    
    Parameters
    ----------
    nps : int
    rhos : float array
    rho : float array
    
    """
    _arespy_pkg.f90wrap_sumrhos(nps=nps, rhos=rhos, rho=rho)

def _fft_sph_r2c(array_r):
    """
    array_c = _fft_sph_r2c(array_r)
    
    
    Defined at Struct_module.fpp lines 2030-2056
    
    Parameters
    ----------
    array_r : float array
    
    Returns
    -------
    array_c : complex array
    
    """
    array_c = _arespy_pkg.f90wrap_fft_sph_r2c(array_r=array_r)
    return array_c

def _fft_sph_c2r(array_c):
    """
    array_r = _fft_sph_c2r(array_c)
    
    
    Defined at Struct_module.fpp lines 2059-2083
    
    Parameters
    ----------
    array_c : complex array
    
    Returns
    -------
    array_r : float array
    
    """
    array_r = _arespy_pkg.f90wrap_fft_sph_c2r(array_c=array_c)
    return array_r

def fft_sph(*args, **kwargs):
    """
    fft_sph(*args, **kwargs)
    
    
    Defined at Struct_module.fpp lines 1115-1117
    
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
    
    
    Defined at Struct_module.fpp line 1100
    
    """
    return _arespy_pkg.f90wrap_grid_module__get__n1()

def set_n1(n1):
    _arespy_pkg.f90wrap_grid_module__set__n1(n1)

def get_n2():
    """
    Element n2 ftype=integer(i4b) pytype=int
    
    
    Defined at Struct_module.fpp line 1100
    
    """
    return _arespy_pkg.f90wrap_grid_module__get__n2()

def set_n2(n2):
    _arespy_pkg.f90wrap_grid_module__set__n2(n2)

def get_n3():
    """
    Element n3 ftype=integer(i4b) pytype=int
    
    
    Defined at Struct_module.fpp line 1100
    
    """
    return _arespy_pkg.f90wrap_grid_module__get__n3()

def set_n3(n3):
    _arespy_pkg.f90wrap_grid_module__set__n3(n3)

def get_n():
    """
    Element n ftype=integer(i4b) pytype=int
    
    
    Defined at Struct_module.fpp line 1100
    
    """
    return _arespy_pkg.f90wrap_grid_module__get__n()

def set_n(n):
    _arespy_pkg.f90wrap_grid_module__set__n(n)

def get_ng1():
    """
    Element ng1 ftype=integer(i4b) pytype=int
    
    
    Defined at Struct_module.fpp line 1101
    
    """
    return _arespy_pkg.f90wrap_grid_module__get__ng1()

def set_ng1(ng1):
    _arespy_pkg.f90wrap_grid_module__set__ng1(ng1)

def get_ng2():
    """
    Element ng2 ftype=integer(i4b) pytype=int
    
    
    Defined at Struct_module.fpp line 1101
    
    """
    return _arespy_pkg.f90wrap_grid_module__get__ng2()

def set_ng2(ng2):
    _arespy_pkg.f90wrap_grid_module__set__ng2(ng2)

def get_ng3():
    """
    Element ng3 ftype=integer(i4b) pytype=int
    
    
    Defined at Struct_module.fpp line 1101
    
    """
    return _arespy_pkg.f90wrap_grid_module__get__ng3()

def set_ng3(ng3):
    _arespy_pkg.f90wrap_grid_module__set__ng3(ng3)

def get_ng():
    """
    Element ng ftype=integer(i4b) pytype=int
    
    
    Defined at Struct_module.fpp line 1101
    
    """
    return _arespy_pkg.f90wrap_grid_module__get__ng()

def set_ng(ng):
    _arespy_pkg.f90wrap_grid_module__set__ng(ng)

def get_array_gap():
    """
    Element gap ftype=real(dp) pytype=float
    
    
    Defined at Struct_module.fpp line 1102
    
    """
    global gap
    array_ndim, array_type, array_shape, array_handle = \
        _arespy_pkg.f90wrap_grid_module__array__gap(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        gap = _arrays[array_handle]
    else:
        gap = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_grid_module__array__gap)
        _arrays[array_handle] = gap
    return gap

def set_array_gap(gap):
    gap[...] = gap

def get_dvol():
    """
    Element dvol ftype=real(dp) pytype=float
    
    
    Defined at Struct_module.fpp line 1103
    
    """
    return _arespy_pkg.f90wrap_grid_module__get__dvol()

def set_dvol(dvol):
    _arespy_pkg.f90wrap_grid_module__set__dvol(dvol)

def get_nk1():
    """
    Element nk1 ftype=integer(i4b) pytype=int
    
    
    Defined at Struct_module.fpp line 1105
    
    """
    return _arespy_pkg.f90wrap_grid_module__get__nk1()

def set_nk1(nk1):
    _arespy_pkg.f90wrap_grid_module__set__nk1(nk1)

def get_nk2():
    """
    Element nk2 ftype=integer(i4b) pytype=int
    
    
    Defined at Struct_module.fpp line 1105
    
    """
    return _arespy_pkg.f90wrap_grid_module__get__nk2()

def set_nk2(nk2):
    _arespy_pkg.f90wrap_grid_module__set__nk2(nk2)

def get_nk3():
    """
    Element nk3 ftype=integer(i4b) pytype=int
    
    
    Defined at Struct_module.fpp line 1105
    
    """
    return _arespy_pkg.f90wrap_grid_module__get__nk3()

def set_nk3(nk3):
    _arespy_pkg.f90wrap_grid_module__set__nk3(nk3)

def get_nk():
    """
    Element nk ftype=integer(i4b) pytype=int
    
    
    Defined at Struct_module.fpp line 1105
    
    """
    return _arespy_pkg.f90wrap_grid_module__get__nk()

def set_nk(nk):
    _arespy_pkg.f90wrap_grid_module__set__nk(nk)

def get_array_kdispl():
    """
    Element kdispl ftype=real(dp) pytype=float
    
    
    Defined at Struct_module.fpp line 1106
    
    """
    global kdispl
    array_ndim, array_type, array_shape, array_handle = \
        _arespy_pkg.f90wrap_grid_module__array__kdispl(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        kdispl = _arrays[array_handle]
    else:
        kdispl = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_grid_module__array__kdispl)
        _arrays[array_handle] = kdispl
    return kdispl

def set_array_kdispl(kdispl):
    kdispl[...] = kdispl

def get_global_n1():
    """
    Element global_n1 ftype=integer(i4b) pytype=int
    
    
    Defined at Struct_module.fpp line 1114
    
    """
    return _arespy_pkg.f90wrap_grid_module__get__global_n1()

def set_global_n1(global_n1):
    _arespy_pkg.f90wrap_grid_module__set__global_n1(global_n1)

def get_global_n2():
    """
    Element global_n2 ftype=integer(i4b) pytype=int
    
    
    Defined at Struct_module.fpp line 1114
    
    """
    return _arespy_pkg.f90wrap_grid_module__get__global_n2()

def set_global_n2(global_n2):
    _arespy_pkg.f90wrap_grid_module__set__global_n2(global_n2)

def get_global_n3():
    """
    Element global_n3 ftype=integer(i4b) pytype=int
    
    
    Defined at Struct_module.fpp line 1114
    
    """
    return _arespy_pkg.f90wrap_grid_module__get__global_n3()

def set_global_n3(global_n3):
    _arespy_pkg.f90wrap_grid_module__set__global_n3(global_n3)

def get_global_n():
    """
    Element global_n ftype=integer(i4b) pytype=int
    
    
    Defined at Struct_module.fpp line 1114
    
    """
    return _arespy_pkg.f90wrap_grid_module__get__global_n()

def set_global_n(global_n):
    _arespy_pkg.f90wrap_grid_module__set__global_n(global_n)


_array_initialisers = [get_array_gap, get_array_kdispl]
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "grid_module".')

for func in _dt_array_initialisers:
    func()

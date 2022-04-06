"""
Module smpi_math_module


Defined at Smpi_math_module.fpp lines 5-2105

"""
from __future__ import print_function, absolute_import, division
import _arespy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("arespy_pkg.parallel_type")
class parallel_type(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=parallel_type)
    
    
    Defined at Smpi_math_module.fpp lines 10-30
    
    """
    def __init__(self, handle=None):
        """
        self = Parallel_Type()
        
        
        Defined at Smpi_math_module.fpp lines 10-30
        
        
        Returns
        -------
        this : Parallel_Type
        	Object to be constructed
        
        
        Automatically generated constructor for parallel_type
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _arespy_pkg.f90wrap_parallel_type_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Parallel_Type
        
        
        Defined at Smpi_math_module.fpp lines 10-30
        
        Parameters
        ----------
        this : Parallel_Type
        	Object to be destructed
        
        
        Automatically generated destructor for parallel_type
        """
        if self._alloc:
            _arespy_pkg.f90wrap_parallel_type_finalise(this=self._handle)
    
    @property
    def comm(self):
        """
        Element comm ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 11
        
        """
        return _arespy_pkg.f90wrap_parallel_type__get__comm(self._handle)
    
    @comm.setter
    def comm(self, comm):
        _arespy_pkg.f90wrap_parallel_type__set__comm(self._handle, comm)
    
    @property
    def myid(self):
        """
        Element myid ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 12
        
        """
        return _arespy_pkg.f90wrap_parallel_type__get__myid(self._handle)
    
    @myid.setter
    def myid(self, myid):
        _arespy_pkg.f90wrap_parallel_type__set__myid(self._handle, myid)
    
    @property
    def numprocs(self):
        """
        Element numprocs ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 13
        
        """
        return _arespy_pkg.f90wrap_parallel_type__get__numprocs(self._handle)
    
    @numprocs.setter
    def numprocs(self, numprocs):
        _arespy_pkg.f90wrap_parallel_type__set__numprocs(self._handle, numprocs)
    
    @property
    def rootid(self):
        """
        Element rootid ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 14
        
        """
        return _arespy_pkg.f90wrap_parallel_type__get__rootid(self._handle)
    
    @rootid.setter
    def rootid(self, rootid):
        _arespy_pkg.f90wrap_parallel_type__set__rootid(self._handle, rootid)
    
    @property
    def isroot(self):
        """
        Element isroot ftype=logical pytype=bool
        
        
        Defined at Smpi_math_module.fpp line 15
        
        """
        return _arespy_pkg.f90wrap_parallel_type__get__isroot(self._handle)
    
    @isroot.setter
    def isroot(self, isroot):
        _arespy_pkg.f90wrap_parallel_type__set__isroot(self._handle, isroot)
    
    @property
    def nstate_proc(self):
        """
        Element nstate_proc ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 17
        
        """
        return _arespy_pkg.f90wrap_parallel_type__get__nstate_proc(self._handle)
    
    @nstate_proc.setter
    def nstate_proc(self, nstate_proc):
        _arespy_pkg.f90wrap_parallel_type__set__nstate_proc(self._handle, nstate_proc)
    
    @property
    def sub2sum(self):
        """
        Element sub2sum ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 19
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_parallel_type__array__sub2sum(self._handle)
        if array_handle in self._arrays:
            sub2sum = self._arrays[array_handle]
        else:
            sub2sum = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_parallel_type__array__sub2sum)
            self._arrays[array_handle] = sub2sum
        return sub2sum
    
    @sub2sum.setter
    def sub2sum(self, sub2sum):
        self.sub2sum[...] = sub2sum
    
    @property
    def mygrid_range(self):
        """
        Element mygrid_range ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 20
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_parallel_type__array__mygrid_range(self._handle)
        if array_handle in self._arrays:
            mygrid_range = self._arrays[array_handle]
        else:
            mygrid_range = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_parallel_type__array__mygrid_range)
            self._arrays[array_handle] = mygrid_range
        return mygrid_range
    
    @mygrid_range.setter
    def mygrid_range(self, mygrid_range):
        self.mygrid_range[...] = mygrid_range
    
    @property
    def recvcounts(self):
        """
        Element recvcounts ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 21
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_parallel_type__array__recvcounts(self._handle)
        if array_handle in self._arrays:
            recvcounts = self._arrays[array_handle]
        else:
            recvcounts = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_parallel_type__array__recvcounts)
            self._arrays[array_handle] = recvcounts
        return recvcounts
    
    @recvcounts.setter
    def recvcounts(self, recvcounts):
        self.recvcounts[...] = recvcounts
    
    @property
    def displs(self):
        """
        Element displs ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 22
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_parallel_type__array__displs(self._handle)
        if array_handle in self._arrays:
            displs = self._arrays[array_handle]
        else:
            displs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_parallel_type__array__displs)
            self._arrays[array_handle] = displs
        return displs
    
    @displs.setter
    def displs(self, displs):
        self.displs[...] = displs
    
    @property
    def global_gridrange(self):
        """
        Element global_gridrange ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 23
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_parallel_type__array__global_gridrange(self._handle)
        if array_handle in self._arrays:
            global_gridrange = self._arrays[array_handle]
        else:
            global_gridrange = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_parallel_type__array__global_gridrange)
            self._arrays[array_handle] = global_gridrange
        return global_gridrange
    
    @global_gridrange.setter
    def global_gridrange(self, global_gridrange):
        self.global_gridrange[...] = global_gridrange
    
    @property
    def comm2d(self):
        """
        Element comm2d ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 25
        
        """
        return _arespy_pkg.f90wrap_parallel_type__get__comm2d(self._handle)
    
    @comm2d.setter
    def comm2d(self, comm2d):
        _arespy_pkg.f90wrap_parallel_type__set__comm2d(self._handle, comm2d)
    
    @property
    def commx(self):
        """
        Element commx ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 25
        
        """
        return _arespy_pkg.f90wrap_parallel_type__get__commx(self._handle)
    
    @commx.setter
    def commx(self, commx):
        _arespy_pkg.f90wrap_parallel_type__set__commx(self._handle, commx)
    
    @property
    def commy(self):
        """
        Element commy ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 25
        
        """
        return _arespy_pkg.f90wrap_parallel_type__get__commy(self._handle)
    
    @commy.setter
    def commy(self, commy):
        _arespy_pkg.f90wrap_parallel_type__set__commy(self._handle, commy)
    
    @property
    def rankx(self):
        """
        Element rankx ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 25
        
        """
        return _arespy_pkg.f90wrap_parallel_type__get__rankx(self._handle)
    
    @rankx.setter
    def rankx(self, rankx):
        _arespy_pkg.f90wrap_parallel_type__set__rankx(self._handle, rankx)
    
    @property
    def ranky(self):
        """
        Element ranky ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 25
        
        """
        return _arespy_pkg.f90wrap_parallel_type__get__ranky(self._handle)
    
    @ranky.setter
    def ranky(self, ranky):
        _arespy_pkg.f90wrap_parallel_type__set__ranky(self._handle, ranky)
    
    @property
    def periods(self):
        """
        Element periods ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 25
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_parallel_type__array__periods(self._handle)
        if array_handle in self._arrays:
            periods = self._arrays[array_handle]
        else:
            periods = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_parallel_type__array__periods)
            self._arrays[array_handle] = periods
        return periods
    
    @periods.setter
    def periods(self, periods):
        self.periods[...] = periods
    
    @property
    def reorder(self):
        """
        Element reorder ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 25
        
        """
        return _arespy_pkg.f90wrap_parallel_type__get__reorder(self._handle)
    
    @reorder.setter
    def reorder(self, reorder):
        _arespy_pkg.f90wrap_parallel_type__set__reorder(self._handle, reorder)
    
    @property
    def remainx(self):
        """
        Element remainx ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 25
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_parallel_type__array__remainx(self._handle)
        if array_handle in self._arrays:
            remainx = self._arrays[array_handle]
        else:
            remainx = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_parallel_type__array__remainx)
            self._arrays[array_handle] = remainx
        return remainx
    
    @remainx.setter
    def remainx(self, remainx):
        self.remainx[...] = remainx
    
    @property
    def remainy(self):
        """
        Element remainy ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 25
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_parallel_type__array__remainy(self._handle)
        if array_handle in self._arrays:
            remainy = self._arrays[array_handle]
        else:
            remainy = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_parallel_type__array__remainy)
            self._arrays[array_handle] = remainy
        return remainy
    
    @remainy.setter
    def remainy(self, remainy):
        self.remainy[...] = remainy
    
    @property
    def ndims(self):
        """
        Element ndims ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 25
        
        """
        return _arespy_pkg.f90wrap_parallel_type__get__ndims(self._handle)
    
    @ndims.setter
    def ndims(self, ndims):
        _arespy_pkg.f90wrap_parallel_type__set__ndims(self._handle, ndims)
    
    @property
    def dims(self):
        """
        Element dims ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 25
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_parallel_type__array__dims(self._handle)
        if array_handle in self._arrays:
            dims = self._arrays[array_handle]
        else:
            dims = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_parallel_type__array__dims)
            self._arrays[array_handle] = dims
        return dims
    
    @dims.setter
    def dims(self, dims):
        self.dims[...] = dims
    
    @property
    def commfft(self):
        """
        Element commfft ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 27
        
        """
        return _arespy_pkg.f90wrap_parallel_type__get__commfft(self._handle)
    
    @commfft.setter
    def commfft(self, commfft):
        _arespy_pkg.f90wrap_parallel_type__set__commfft(self._handle, commfft)
    
    @property
    def local_z(self):
        """
        Element local_z ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 27
        
        """
        return _arespy_pkg.f90wrap_parallel_type__get__local_z(self._handle)
    
    @local_z.setter
    def local_z(self, local_z):
        _arespy_pkg.f90wrap_parallel_type__set__local_z(self._handle, local_z)
    
    @property
    def local_z_start(self):
        """
        Element local_z_start ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 27
        
        """
        return _arespy_pkg.f90wrap_parallel_type__get__local_z_start(self._handle)
    
    @local_z_start.setter
    def local_z_start(self, local_z_start):
        _arespy_pkg.f90wrap_parallel_type__set__local_z_start(self._handle, \
            local_z_start)
    
    @property
    def fft_grid_range(self):
        """
        Element fft_grid_range ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 28
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_parallel_type__array__fft_grid_range(self._handle)
        if array_handle in self._arrays:
            fft_grid_range = self._arrays[array_handle]
        else:
            fft_grid_range = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_parallel_type__array__fft_grid_range)
            self._arrays[array_handle] = fft_grid_range
        return fft_grid_range
    
    @fft_grid_range.setter
    def fft_grid_range(self, fft_grid_range):
        self.fft_grid_range[...] = fft_grid_range
    
    @property
    def fft_rcount(self):
        """
        Element fft_rcount ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 29
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_parallel_type__array__fft_rcount(self._handle)
        if array_handle in self._arrays:
            fft_rcount = self._arrays[array_handle]
        else:
            fft_rcount = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_parallel_type__array__fft_rcount)
            self._arrays[array_handle] = fft_rcount
        return fft_rcount
    
    @fft_rcount.setter
    def fft_rcount(self, fft_rcount):
        self.fft_rcount[...] = fft_rcount
    
    @property
    def fft_rdispls(self):
        """
        Element fft_rdispls ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 29
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_parallel_type__array__fft_rdispls(self._handle)
        if array_handle in self._arrays:
            fft_rdispls = self._arrays[array_handle]
        else:
            fft_rdispls = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_parallel_type__array__fft_rdispls)
            self._arrays[array_handle] = fft_rdispls
        return fft_rdispls
    
    @fft_rdispls.setter
    def fft_rdispls(self, fft_rdispls):
        self.fft_rdispls[...] = fft_rdispls
    
    @property
    def fft_scount(self):
        """
        Element fft_scount ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 30
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_parallel_type__array__fft_scount(self._handle)
        if array_handle in self._arrays:
            fft_scount = self._arrays[array_handle]
        else:
            fft_scount = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_parallel_type__array__fft_scount)
            self._arrays[array_handle] = fft_scount
        return fft_scount
    
    @fft_scount.setter
    def fft_scount(self, fft_scount):
        self.fft_scount[...] = fft_scount
    
    @property
    def fft_sdispls(self):
        """
        Element fft_sdispls ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 30
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_parallel_type__array__fft_sdispls(self._handle)
        if array_handle in self._arrays:
            fft_sdispls = self._arrays[array_handle]
        else:
            fft_sdispls = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_parallel_type__array__fft_sdispls)
            self._arrays[array_handle] = fft_sdispls
        return fft_sdispls
    
    @fft_sdispls.setter
    def fft_sdispls(self, fft_sdispls):
        self.fft_sdispls[...] = fft_sdispls
    
    def __str__(self):
        ret = ['<parallel_type>{\n']
        ret.append('    comm : ')
        ret.append(repr(self.comm))
        ret.append(',\n    myid : ')
        ret.append(repr(self.myid))
        ret.append(',\n    numprocs : ')
        ret.append(repr(self.numprocs))
        ret.append(',\n    rootid : ')
        ret.append(repr(self.rootid))
        ret.append(',\n    isroot : ')
        ret.append(repr(self.isroot))
        ret.append(',\n    nstate_proc : ')
        ret.append(repr(self.nstate_proc))
        ret.append(',\n    sub2sum : ')
        ret.append(repr(self.sub2sum))
        ret.append(',\n    mygrid_range : ')
        ret.append(repr(self.mygrid_range))
        ret.append(',\n    recvcounts : ')
        ret.append(repr(self.recvcounts))
        ret.append(',\n    displs : ')
        ret.append(repr(self.displs))
        ret.append(',\n    global_gridrange : ')
        ret.append(repr(self.global_gridrange))
        ret.append(',\n    comm2d : ')
        ret.append(repr(self.comm2d))
        ret.append(',\n    commx : ')
        ret.append(repr(self.commx))
        ret.append(',\n    commy : ')
        ret.append(repr(self.commy))
        ret.append(',\n    rankx : ')
        ret.append(repr(self.rankx))
        ret.append(',\n    ranky : ')
        ret.append(repr(self.ranky))
        ret.append(',\n    periods : ')
        ret.append(repr(self.periods))
        ret.append(',\n    reorder : ')
        ret.append(repr(self.reorder))
        ret.append(',\n    remainx : ')
        ret.append(repr(self.remainx))
        ret.append(',\n    remainy : ')
        ret.append(repr(self.remainy))
        ret.append(',\n    ndims : ')
        ret.append(repr(self.ndims))
        ret.append(',\n    dims : ')
        ret.append(repr(self.dims))
        ret.append(',\n    commfft : ')
        ret.append(repr(self.commfft))
        ret.append(',\n    local_z : ')
        ret.append(repr(self.local_z))
        ret.append(',\n    local_z_start : ')
        ret.append(repr(self.local_z_start))
        ret.append(',\n    fft_grid_range : ')
        ret.append(repr(self.fft_grid_range))
        ret.append(',\n    fft_rcount : ')
        ret.append(repr(self.fft_rcount))
        ret.append(',\n    fft_rdispls : ')
        ret.append(repr(self.fft_rdispls))
        ret.append(',\n    fft_scount : ')
        ret.append(repr(self.fft_scount))
        ret.append(',\n    fft_sdispls : ')
        ret.append(repr(self.fft_sdispls))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("arespy_pkg.smpi_root_type")
class smpi_root_type(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=smpi_root_type)
    
    
    Defined at Smpi_math_module.fpp lines 34-35
    
    """
    def __init__(self, handle=None):
        """
        self = Smpi_Root_Type()
        
        
        Defined at Smpi_math_module.fpp lines 34-35
        
        
        Returns
        -------
        this : Smpi_Root_Type
        	Object to be constructed
        
        
        Automatically generated constructor for smpi_root_type
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _arespy_pkg.f90wrap_smpi_root_type_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Smpi_Root_Type
        
        
        Defined at Smpi_math_module.fpp lines 34-35
        
        Parameters
        ----------
        this : Smpi_Root_Type
        	Object to be destructed
        
        
        Automatically generated destructor for smpi_root_type
        """
        if self._alloc:
            _arespy_pkg.f90wrap_smpi_root_type_finalise(this=self._handle)
    
    @property
    def natom_group(self):
        """
        Element natom_group ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 35
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_smpi_root_type__array__natom_group(self._handle)
        if array_handle in self._arrays:
            natom_group = self._arrays[array_handle]
        else:
            natom_group = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_smpi_root_type__array__natom_group)
            self._arrays[array_handle] = natom_group
        return natom_group
    
    @natom_group.setter
    def natom_group(self, natom_group):
        self.natom_group[...] = natom_group
    
    def __str__(self):
        ret = ['<smpi_root_type>{\n']
        ret.append('    natom_group : ')
        ret.append(repr(self.natom_group))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("arespy_pkg.smpi_comm_type")
class smpi_comm_type(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=smpi_comm_type)
    
    
    Defined at Smpi_math_module.fpp lines 39-41
    
    """
    def __init__(self, handle=None):
        """
        self = Smpi_Comm_Type()
        
        
        Defined at Smpi_math_module.fpp lines 39-41
        
        
        Returns
        -------
        this : Smpi_Comm_Type
        	Object to be constructed
        
        
        Automatically generated constructor for smpi_comm_type
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _arespy_pkg.f90wrap_smpi_comm_type_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Smpi_Comm_Type
        
        
        Defined at Smpi_math_module.fpp lines 39-41
        
        Parameters
        ----------
        this : Smpi_Comm_Type
        	Object to be destructed
        
        
        Automatically generated destructor for smpi_comm_type
        """
        if self._alloc:
            _arespy_pkg.f90wrap_smpi_comm_type_finalise(this=self._handle)
    
    @property
    def atoms(self):
        """
        Element atoms ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 40
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_smpi_comm_type__array__atoms(self._handle)
        if array_handle in self._arrays:
            atoms = self._arrays[array_handle]
        else:
            atoms = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_smpi_comm_type__array__atoms)
            self._arrays[array_handle] = atoms
        return atoms
    
    @atoms.setter
    def atoms(self, atoms):
        self.atoms[...] = atoms
    
    @property
    def displs(self):
        """
        Element displs ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 41
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_smpi_comm_type__array__displs(self._handle)
        if array_handle in self._arrays:
            displs = self._arrays[array_handle]
        else:
            displs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_smpi_comm_type__array__displs)
            self._arrays[array_handle] = displs
        return displs
    
    @displs.setter
    def displs(self, displs):
        self.displs[...] = displs
    
    def __str__(self):
        ret = ['<smpi_comm_type>{\n']
        ret.append('    atoms : ')
        ret.append(repr(self.atoms))
        ret.append(',\n    displs : ')
        ret.append(repr(self.displs))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("arespy_pkg.time_type")
class time_type(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=time_type)
    
    
    Defined at Smpi_math_module.fpp lines 45-51
    
    """
    def __init__(self, handle=None):
        """
        self = Time_Type()
        
        
        Defined at Smpi_math_module.fpp lines 45-51
        
        
        Returns
        -------
        this : Time_Type
        	Object to be constructed
        
        
        Automatically generated constructor for time_type
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _arespy_pkg.f90wrap_time_type_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Time_Type
        
        
        Defined at Smpi_math_module.fpp lines 45-51
        
        Parameters
        ----------
        this : Time_Type
        	Object to be destructed
        
        
        Automatically generated destructor for time_type
        """
        if self._alloc:
            _arespy_pkg.f90wrap_time_type_finalise(this=self._handle)
    
    @property
    def label(self):
        """
        Element label ftype=character(len=100) pytype=str
        
        
        Defined at Smpi_math_module.fpp line 46
        
        """
        return _arespy_pkg.f90wrap_time_type__get__label(self._handle)
    
    @label.setter
    def label(self, label):
        _arespy_pkg.f90wrap_time_type__set__label(self._handle, label)
    
    @property
    def tic(self):
        """
        Element tic ftype=real(dp) pytype=float
        
        
        Defined at Smpi_math_module.fpp line 47
        
        """
        return _arespy_pkg.f90wrap_time_type__get__tic(self._handle)
    
    @tic.setter
    def tic(self, tic):
        _arespy_pkg.f90wrap_time_type__set__tic(self._handle, tic)
    
    @property
    def toc(self):
        """
        Element toc ftype=real(dp) pytype=float
        
        
        Defined at Smpi_math_module.fpp line 48
        
        """
        return _arespy_pkg.f90wrap_time_type__get__toc(self._handle)
    
    @toc.setter
    def toc(self, toc):
        _arespy_pkg.f90wrap_time_type__set__toc(self._handle, toc)
    
    @property
    def total(self):
        """
        Element total ftype=real(dp) pytype=float
        
        
        Defined at Smpi_math_module.fpp line 49
        
        """
        return _arespy_pkg.f90wrap_time_type__get__total(self._handle)
    
    @total.setter
    def total(self, total):
        _arespy_pkg.f90wrap_time_type__set__total(self._handle, total)
    
    @property
    def sum_total(self):
        """
        Element sum_total ftype=real(dp) pytype=float
        
        
        Defined at Smpi_math_module.fpp line 50
        
        """
        return _arespy_pkg.f90wrap_time_type__get__sum_total(self._handle)
    
    @sum_total.setter
    def sum_total(self, sum_total):
        _arespy_pkg.f90wrap_time_type__set__sum_total(self._handle, sum_total)
    
    @property
    def num(self):
        """
        Element num ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 51
        
        """
        return _arespy_pkg.f90wrap_time_type__get__num(self._handle)
    
    @num.setter
    def num(self, num):
        _arespy_pkg.f90wrap_time_type__set__num(self._handle, num)
    
    def __str__(self):
        ret = ['<time_type>{\n']
        ret.append('    label : ')
        ret.append(repr(self.label))
        ret.append(',\n    tic : ')
        ret.append(repr(self.tic))
        ret.append(',\n    toc : ')
        ret.append(repr(self.toc))
        ret.append(',\n    total : ')
        ret.append(repr(self.total))
        ret.append(',\n    sum_total : ')
        ret.append(repr(self.sum_total))
        ret.append(',\n    num : ')
        ret.append(repr(self.num))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("arespy_pkg.mem_type")
class mem_type(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=mem_type)
    
    
    Defined at Smpi_math_module.fpp lines 53-57
    
    """
    def __init__(self, handle=None):
        """
        self = Mem_Type()
        
        
        Defined at Smpi_math_module.fpp lines 53-57
        
        
        Returns
        -------
        this : Mem_Type
        	Object to be constructed
        
        
        Automatically generated constructor for mem_type
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _arespy_pkg.f90wrap_mem_type_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Mem_Type
        
        
        Defined at Smpi_math_module.fpp lines 53-57
        
        Parameters
        ----------
        this : Mem_Type
        	Object to be destructed
        
        
        Automatically generated destructor for mem_type
        """
        if self._alloc:
            _arespy_pkg.f90wrap_mem_type_finalise(this=self._handle)
    
    @property
    def label(self):
        """
        Element label ftype=character(len=100) pytype=str
        
        
        Defined at Smpi_math_module.fpp line 54
        
        """
        return _arespy_pkg.f90wrap_mem_type__get__label(self._handle)
    
    @label.setter
    def label(self, label):
        _arespy_pkg.f90wrap_mem_type__set__label(self._handle, label)
    
    @property
    def memic(self):
        """
        Element memic ftype=real(dp) pytype=float
        
        
        Defined at Smpi_math_module.fpp line 55
        
        """
        return _arespy_pkg.f90wrap_mem_type__get__memic(self._handle)
    
    @memic.setter
    def memic(self, memic):
        _arespy_pkg.f90wrap_mem_type__set__memic(self._handle, memic)
    
    @property
    def total(self):
        """
        Element total ftype=real(dp) pytype=float
        
        
        Defined at Smpi_math_module.fpp line 56
        
        """
        return _arespy_pkg.f90wrap_mem_type__get__total(self._handle)
    
    @total.setter
    def total(self, total):
        _arespy_pkg.f90wrap_mem_type__set__total(self._handle, total)
    
    @property
    def num(self):
        """
        Element num ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 57
        
        """
        return _arespy_pkg.f90wrap_mem_type__get__num(self._handle)
    
    @num.setter
    def num(self, num):
        _arespy_pkg.f90wrap_mem_type__set__num(self._handle, num)
    
    def __str__(self):
        ret = ['<mem_type>{\n']
        ret.append('    label : ')
        ret.append(repr(self.label))
        ret.append(',\n    memic : ')
        ret.append(repr(self.memic))
        ret.append(',\n    total : ')
        ret.append(repr(self.total))
        ret.append(',\n    num : ')
        ret.append(repr(self.num))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("arespy_pkg.grid_diff_map_type")
class grid_diff_map_type(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=grid_diff_map_type)
    
    
    Defined at Smpi_math_module.fpp lines 64-74
    
    """
    def __init__(self, handle=None):
        """
        self = Grid_Diff_Map_Type()
        
        
        Defined at Smpi_math_module.fpp lines 64-74
        
        
        Returns
        -------
        this : Grid_Diff_Map_Type
        	Object to be constructed
        
        
        Automatically generated constructor for grid_diff_map_type
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _arespy_pkg.f90wrap_grid_diff_map_type_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Grid_Diff_Map_Type
        
        
        Defined at Smpi_math_module.fpp lines 64-74
        
        Parameters
        ----------
        this : Grid_Diff_Map_Type
        	Object to be destructed
        
        
        Automatically generated destructor for grid_diff_map_type
        """
        if self._alloc:
            _arespy_pkg.f90wrap_grid_diff_map_type_finalise(this=self._handle)
    
    @property
    def nz_map(self):
        """
        Element nz_map ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 65
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_grid_diff_map_type__array__nz_map(self._handle)
        if array_handle in self._arrays:
            nz_map = self._arrays[array_handle]
        else:
            nz_map = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_grid_diff_map_type__array__nz_map)
            self._arrays[array_handle] = nz_map
        return nz_map
    
    @nz_map.setter
    def nz_map(self, nz_map):
        self.nz_map[...] = nz_map
    
    @property
    def mycomm_cores(self):
        """
        Element mycomm_cores ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 66
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_grid_diff_map_type__array__mycomm_cores(self._handle)
        if array_handle in self._arrays:
            mycomm_cores = self._arrays[array_handle]
        else:
            mycomm_cores = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_grid_diff_map_type__array__mycomm_cores)
            self._arrays[array_handle] = mycomm_cores
        return mycomm_cores
    
    @mycomm_cores.setter
    def mycomm_cores(self, mycomm_cores):
        self.mycomm_cores[...] = mycomm_cores
    
    @property
    def mycomm_size(self):
        """
        Element mycomm_size ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 67
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_grid_diff_map_type__array__mycomm_size(self._handle)
        if array_handle in self._arrays:
            mycomm_size = self._arrays[array_handle]
        else:
            mycomm_size = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_grid_diff_map_type__array__mycomm_size)
            self._arrays[array_handle] = mycomm_size
        return mycomm_size
    
    @mycomm_size.setter
    def mycomm_size(self, mycomm_size):
        self.mycomm_size[...] = mycomm_size
    
    @property
    def mysend_size(self):
        """
        Element mysend_size ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 68
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_grid_diff_map_type__array__mysend_size(self._handle)
        if array_handle in self._arrays:
            mysend_size = self._arrays[array_handle]
        else:
            mysend_size = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_grid_diff_map_type__array__mysend_size)
            self._arrays[array_handle] = mysend_size
        return mysend_size
    
    @mysend_size.setter
    def mysend_size(self, mysend_size):
        self.mysend_size[...] = mysend_size
    
    @property
    def local_map(self):
        """
        Element local_map ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 69
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_grid_diff_map_type__array__local_map(self._handle)
        if array_handle in self._arrays:
            local_map = self._arrays[array_handle]
        else:
            local_map = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_grid_diff_map_type__array__local_map)
            self._arrays[array_handle] = local_map
        return local_map
    
    @local_map.setter
    def local_map(self, local_map):
        self.local_map[...] = local_map
    
    @property
    def local_map1d(self):
        """
        Element local_map1d ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 70
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_grid_diff_map_type__array__local_map1d(self._handle)
        if array_handle in self._arrays:
            local_map1d = self._arrays[array_handle]
        else:
            local_map1d = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_grid_diff_map_type__array__local_map1d)
            self._arrays[array_handle] = local_map1d
        return local_map1d
    
    @local_map1d.setter
    def local_map1d(self, local_map1d):
        self.local_map1d[...] = local_map1d
    
    @property
    def boundary(self):
        """
        Element boundary ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 71
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_grid_diff_map_type__array__boundary(self._handle)
        if array_handle in self._arrays:
            boundary = self._arrays[array_handle]
        else:
            boundary = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_grid_diff_map_type__array__boundary)
            self._arrays[array_handle] = boundary
        return boundary
    
    @boundary.setter
    def boundary(self, boundary):
        self.boundary[...] = boundary
    
    @property
    def boundary1d(self):
        """
        Element boundary1d ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 72
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_grid_diff_map_type__array__boundary1d(self._handle)
        if array_handle in self._arrays:
            boundary1d = self._arrays[array_handle]
        else:
            boundary1d = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_grid_diff_map_type__array__boundary1d)
            self._arrays[array_handle] = boundary1d
        return boundary1d
    
    @boundary1d.setter
    def boundary1d(self, boundary1d):
        self.boundary1d[...] = boundary1d
    
    @property
    def rcount(self):
        """
        Element rcount ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 73
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_grid_diff_map_type__array__rcount(self._handle)
        if array_handle in self._arrays:
            rcount = self._arrays[array_handle]
        else:
            rcount = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_grid_diff_map_type__array__rcount)
            self._arrays[array_handle] = rcount
        return rcount
    
    @rcount.setter
    def rcount(self, rcount):
        self.rcount[...] = rcount
    
    @property
    def rdispls(self):
        """
        Element rdispls ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 73
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_grid_diff_map_type__array__rdispls(self._handle)
        if array_handle in self._arrays:
            rdispls = self._arrays[array_handle]
        else:
            rdispls = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_grid_diff_map_type__array__rdispls)
            self._arrays[array_handle] = rdispls
        return rdispls
    
    @rdispls.setter
    def rdispls(self, rdispls):
        self.rdispls[...] = rdispls
    
    @property
    def scount(self):
        """
        Element scount ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 74
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_grid_diff_map_type__array__scount(self._handle)
        if array_handle in self._arrays:
            scount = self._arrays[array_handle]
        else:
            scount = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_grid_diff_map_type__array__scount)
            self._arrays[array_handle] = scount
        return scount
    
    @scount.setter
    def scount(self, scount):
        self.scount[...] = scount
    
    @property
    def sdispls(self):
        """
        Element sdispls ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 74
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_grid_diff_map_type__array__sdispls(self._handle)
        if array_handle in self._arrays:
            sdispls = self._arrays[array_handle]
        else:
            sdispls = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_grid_diff_map_type__array__sdispls)
            self._arrays[array_handle] = sdispls
        return sdispls
    
    @sdispls.setter
    def sdispls(self, sdispls):
        self.sdispls[...] = sdispls
    
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
        ret.append(',\n    local_map1d : ')
        ret.append(repr(self.local_map1d))
        ret.append(',\n    boundary : ')
        ret.append(repr(self.boundary))
        ret.append(',\n    boundary1d : ')
        ret.append(repr(self.boundary1d))
        ret.append(',\n    rcount : ')
        ret.append(repr(self.rcount))
        ret.append(',\n    rdispls : ')
        ret.append(repr(self.rdispls))
        ret.append(',\n    scount : ')
        ret.append(repr(self.scount))
        ret.append(',\n    sdispls : ')
        ret.append(repr(self.sdispls))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("arespy_pkg.sphere_type")
class sphere_type(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=sphere_type)
    
    
    Defined at Smpi_math_module.fpp lines 76-78
    
    """
    def __init__(self, handle=None):
        """
        self = Sphere_Type()
        
        
        Defined at Smpi_math_module.fpp lines 76-78
        
        
        Returns
        -------
        this : Sphere_Type
        	Object to be constructed
        
        
        Automatically generated constructor for sphere_type
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _arespy_pkg.f90wrap_sphere_type_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Sphere_Type
        
        
        Defined at Smpi_math_module.fpp lines 76-78
        
        Parameters
        ----------
        this : Sphere_Type
        	Object to be destructed
        
        
        Automatically generated destructor for sphere_type
        """
        if self._alloc:
            _arespy_pkg.f90wrap_sphere_type_finalise(this=self._handle)
    
    @property
    def length(self):
        """
        Element length ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 77
        
        """
        return _arespy_pkg.f90wrap_sphere_type__get__length(self._handle)
    
    @length.setter
    def length(self, length):
        _arespy_pkg.f90wrap_sphere_type__set__length(self._handle, length)
    
    @property
    def map3d(self):
        """
        Element map3d ftype=integer(i4b) pytype=int
        
        
        Defined at Smpi_math_module.fpp line 78
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_sphere_type__array__map3d(self._handle)
        if array_handle in self._arrays:
            map3d = self._arrays[array_handle]
        else:
            map3d = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_sphere_type__array__map3d)
            self._arrays[array_handle] = map3d
        return map3d
    
    @map3d.setter
    def map3d(self, map3d):
        self.map3d[...] = map3d
    
    def __str__(self):
        ret = ['<sphere_type>{\n']
        ret.append('    length : ')
        ret.append(repr(self.length))
        ret.append(',\n    map3d : ')
        ret.append(repr(self.map3d))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def smpi_init():
    """
    smpi_init()
    
    
    Defined at Smpi_math_module.fpp lines 128-139
    
    
    """
    _arespy_pkg.f90wrap_smpi_init()

def smpi_init_pbc():
    """
    smpi_init_pbc()
    
    
    Defined at Smpi_math_module.fpp lines 142-188
    
    
    """
    _arespy_pkg.f90wrap_smpi_init_pbc()

def smpi_exit():
    """
    smpi_exit()
    
    
    Defined at Smpi_math_module.fpp lines 192-196
    
    
    """
    _arespy_pkg.f90wrap_smpi_exit()

def smpi_stop(message):
    """
    smpi_stop(message)
    
    
    Defined at Smpi_math_module.fpp lines 199-204
    
    Parameters
    ----------
    message : str
    
    """
    _arespy_pkg.f90wrap_smpi_stop(message=message)

def smpi_stop_info(message):
    """
    smpi_stop_info(message)
    
    
    Defined at Smpi_math_module.fpp lines 207-212
    
    Parameters
    ----------
    message : str
    
    """
    _arespy_pkg.f90wrap_smpi_stop_info(message=message)

def nstates_split(m, np):
    """
    nstates_split(m, np)
    
    
    Defined at Smpi_math_module.fpp lines 216-235
    
    Parameters
    ----------
    m : int
    np : int
    
    -------------------------------------
     ms = mod(m ,np)
     if(ms == 0 ) then
       m = m /np
     else
    if( ms >= (parallel%myid + 1) ) then
      m = m/parallel%numprocs+1
    else
      m = m / parallel%numprocs
    end if
    if(parallel%myid == 0 ) then
       if(parallel%coords(1) == 0 ) then
         m = m/np + ms
       else
         m = m / np
       end if
     end if
    """
    _arespy_pkg.f90wrap_nstates_split(m=m, np=np)

def nstates_split_2(m, np):
    """
    nstates_split_2(m, np)
    
    
    Defined at Smpi_math_module.fpp lines 239-258
    
    Parameters
    ----------
    m : int
    np : int
    
    -------------------------------------
     ms = mod(m ,np)
     if(ms == 0 ) then
       m = m /np
     else
    if( ms >= (parallel%myid + 1) ) then
      m = m/parallel%numprocs+1
    else
      m = m / parallel%numprocs
    end if
    if(parallel%myid == 0 ) then
       if(parallel%coords(2) == 0 ) then
         m = m/np + ms
       else
         m = m / np
       end if
     end if
    """
    _arespy_pkg.f90wrap_nstates_split_2(m=m, np=np)

def smpi_reduce_sum_real(amat, na, ramat=None):
    """
    smpi_reduce_sum_real(amat, na[, ramat])
    
    
    Defined at Smpi_math_module.fpp lines 513-522
    
    Parameters
    ----------
    amat : float array
    na : int
    ramat : float array
    
    """
    _arespy_pkg.f90wrap_smpi_reduce_sum_real(amat=amat, na=na, ramat=ramat)

def start_time(inlabel, flag, tic=None):
    """
    start_time(inlabel, flag[, tic])
    
    
    Defined at Smpi_math_module.fpp lines 554-574
    
    Parameters
    ----------
    inlabel : str
    flag : bool
    tic : float
    
    """
    _arespy_pkg.f90wrap_start_time(inlabel=inlabel, flag=flag, tic=tic)

def end_time(inlabel, flag, toc=None):
    """
    end_time(inlabel, flag[, toc])
    
    
    Defined at Smpi_math_module.fpp lines 578-597
    
    Parameters
    ----------
    inlabel : str
    flag : bool
    toc : float
    
    """
    _arespy_pkg.f90wrap_end_time(inlabel=inlabel, flag=flag, toc=toc)

def write_time(inlabel, flag):
    """
    write_time(inlabel, flag)
    
    
    Defined at Smpi_math_module.fpp lines 601-613
    
    Parameters
    ----------
    inlabel : str
    flag : bool
    
    """
    _arespy_pkg.f90wrap_write_time(inlabel=inlabel, flag=flag)

def write_sum_time(inlabel, flag):
    """
    write_sum_time(inlabel, flag)
    
    
    Defined at Smpi_math_module.fpp lines 617-629
    
    Parameters
    ----------
    inlabel : str
    flag : bool
    
    """
    _arespy_pkg.f90wrap_write_sum_time(inlabel=inlabel, flag=flag)

def print_time(inlabel, t):
    """
    print_time(inlabel, t)
    
    
    Defined at Smpi_math_module.fpp lines 633-644
    
    Parameters
    ----------
    inlabel : str
    t : float
    
    """
    _arespy_pkg.f90wrap_print_time(inlabel=inlabel, t=t)

def states_split(nev):
    """
    states_split(nev)
    
    
    Defined at Smpi_math_module.fpp lines 684-704
    
    Parameters
    ----------
    nev : int
    
    """
    _arespy_pkg.f90wrap_states_split(nev=nev)

def array_split(nev):
    """
    array_split(nev)
    
    
    Defined at Smpi_math_module.fpp lines 707-755
    
    Parameters
    ----------
    nev : int
    
    """
    _arespy_pkg.f90wrap_array_split(nev=nev)

def grid_split(ngrid, ncore, comm, id, grid_range, recvcounts, displs, \
    gridrange_sum=None, n1=None, n2=None, n3=None, n=None):
    """
    grid_split(ngrid, ncore, comm, id, grid_range, recvcounts, displs[, \
        gridrange_sum, n1, n2, n3, n])
    
    
    Defined at Smpi_math_module.fpp lines 758-813
    
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
    _arespy_pkg.f90wrap_grid_split(ngrid=ngrid, ncore=ncore, comm=comm, id=id, \
        grid_range=grid_range, recvcounts=recvcounts, displs=displs, \
        gridrange_sum=gridrange_sum, n1=n1, n2=n2, n3=n3, n=n)

def atom_split(mysize, natom, atom_index):
    """
    atom_split(mysize, natom, atom_index)
    
    
    Defined at Smpi_math_module.fpp lines 816-852
    
    Parameters
    ----------
    mysize : int
    natom : int
    atom_index : int array
    
    """
    _arespy_pkg.f90wrap_atom_split(mysize=mysize, natom=natom, \
        atom_index=atom_index)

def grid_sphere_init(n1, n2, n3, norder):
    """
    grid_sphere_init(n1, n2, n3, norder)
    
    
    Defined at Smpi_math_module.fpp lines 854-971
    
    Parameters
    ----------
    n1 : int
    n2 : int
    n3 : int
    norder : int
    
    """
    _arespy_pkg.f90wrap_grid_sphere_init(n1=n1, n2=n2, n3=n3, norder=norder)

def set_wrap_grid_iso(myrho, wrap_box):
    """
    set_wrap_grid_iso(myrho, wrap_box)
    
    
    Defined at Smpi_math_module.fpp lines 973-1044
    
    Parameters
    ----------
    myrho : float array
    wrap_box : float array
    
    """
    _arespy_pkg.f90wrap_set_wrap_grid_iso(myrho=myrho, wrap_box=wrap_box)

def destroy_diff_map():
    """
    destroy_diff_map()
    
    
    Defined at Smpi_math_module.fpp lines 1046-1090
    
    
    """
    _arespy_pkg.f90wrap_destroy_diff_map()

def smpi_diff_init_sph(n1, n2, n3, n, norder, cell_mu, lsp):
    """
    smpi_diff_init_sph(n1, n2, n3, n, norder, cell_mu, lsp)
    
    
    Defined at Smpi_math_module.fpp lines 1092-1464
    
    Parameters
    ----------
    n1 : int
    n2 : int
    n3 : int
    n : int
    norder : int
    cell_mu : int array
    lsp : bool array
    
    """
    _arespy_pkg.f90wrap_smpi_diff_init_sph(n1=n1, n2=n2, n3=n3, n=n, norder=norder, \
        cell_mu=cell_mu, lsp=lsp)

def smpi_diff_init(n1, n2, n3, n, norder, cell_mu):
    """
    smpi_diff_init(n1, n2, n3, n, norder, cell_mu)
    
    
    Defined at Smpi_math_module.fpp lines 1466-1798
    
    Parameters
    ----------
    n1 : int
    n2 : int
    n3 : int
    n : int
    norder : int
    cell_mu : int array
    
    """
    _arespy_pkg.f90wrap_smpi_diff_init(n1=n1, n2=n2, n3=n3, n=n, norder=norder, \
        cell_mu=cell_mu)

def set_wrap_sph_pbc_ata_real(myrho, wrap_box1d):
    """
    set_wrap_sph_pbc_ata_real(myrho, wrap_box1d)
    
    
    Defined at Smpi_math_module.fpp lines 1800-1819
    
    Parameters
    ----------
    myrho : float array
    wrap_box1d : float array
    
    """
    _arespy_pkg.f90wrap_set_wrap_sph_pbc_ata_real(myrho=myrho, \
        wrap_box1d=wrap_box1d)

def set_fft_alltoallv(fft_grid_range_temp):
    """
    set_fft_alltoallv(fft_grid_range_temp)
    
    
    Defined at Smpi_math_module.fpp lines 2046-2089
    
    Parameters
    ----------
    fft_grid_range_temp : int array
    
    """
    _arespy_pkg.f90wrap_set_fft_alltoallv(fft_grid_range_temp=fft_grid_range_temp)

def destroy_fft_alltoallv():
    """
    destroy_fft_alltoallv()
    
    
    Defined at Smpi_math_module.fpp lines 2091-2104
    
    
    """
    _arespy_pkg.f90wrap_destroy_fft_alltoallv()

def _sum_real_1d(amat):
    """
    totals = _sum_real_1d(amat)
    
    
    Defined at Smpi_math_module.fpp lines 262-278
    
    Parameters
    ----------
    amat : float array
    
    Returns
    -------
    totals : float
    
    """
    totals = _arespy_pkg.f90wrap_sum_real_1d(amat=amat)
    return totals

def _sum_real_2d(amat, bmat):
    """
    totals = _sum_real_2d(amat, bmat)
    
    
    Defined at Smpi_math_module.fpp lines 302-319
    
    Parameters
    ----------
    amat : float array
    bmat : float array
    
    Returns
    -------
    totals : float
    
    """
    totals = _arespy_pkg.f90wrap_sum_real_2d(amat=amat, bmat=bmat)
    return totals

def _sum_real_3d(amat, bmat, cmat):
    """
    totals = _sum_real_3d(amat, bmat, cmat)
    
    
    Defined at Smpi_math_module.fpp lines 344-361
    
    Parameters
    ----------
    amat : float array
    bmat : float array
    cmat : float array
    
    Returns
    -------
    totals : float
    
    """
    totals = _arespy_pkg.f90wrap_sum_real_3d(amat=amat, bmat=bmat, cmat=cmat)
    return totals

def _sum_cplx_1d(amat):
    """
    totals = _sum_cplx_1d(amat)
    
    
    Defined at Smpi_math_module.fpp lines 282-298
    
    Parameters
    ----------
    amat : complex array
    
    Returns
    -------
    totals : complex
    
    """
    totals = _arespy_pkg.f90wrap_sum_cplx_1d(amat=amat)
    return totals

def _sum_cplx_2d(amat, bmat):
    """
    totals = _sum_cplx_2d(amat, bmat)
    
    
    Defined at Smpi_math_module.fpp lines 323-340
    
    Parameters
    ----------
    amat : complex array
    bmat : complex array
    
    Returns
    -------
    totals : complex
    
    """
    totals = _arespy_pkg.f90wrap_sum_cplx_2d(amat=amat, bmat=bmat)
    return totals

def _sum_cplx_3d(amat, bmat, cmat):
    """
    totals = _sum_cplx_3d(amat, bmat, cmat)
    
    
    Defined at Smpi_math_module.fpp lines 365-382
    
    Parameters
    ----------
    amat : complex array
    bmat : complex array
    cmat : complex array
    
    Returns
    -------
    totals : complex
    
    """
    totals = _arespy_pkg.f90wrap_sum_cplx_3d(amat=amat, bmat=bmat, cmat=cmat)
    return totals

def sompsum(*args, **kwargs):
    """
    sompsum(*args, **kwargs)
    
    
    Defined at Smpi_math_module.fpp lines 87-93
    
    Overloaded interface containing the following procedures:
      _sum_real_1d
      _sum_real_2d
      _sum_real_3d
      _sum_cplx_1d
      _sum_cplx_2d
      _sum_cplx_3d
    
    """
    for proc in [_sum_real_1d, _sum_real_2d, _sum_real_3d, _sum_cplx_1d, \
        _sum_cplx_2d, _sum_cplx_3d]:
        try:
            return proc(*args, **kwargs)
        except TypeError:
            continue

def _smpi_sum_int_1s(x):
    """
    sumx = _smpi_sum_int_1s(x)
    
    
    Defined at Smpi_math_module.fpp lines 428-430
    
    Parameters
    ----------
    x : int
    
    Returns
    -------
    sumx : int
    
    """
    sumx = _arespy_pkg.f90wrap_smpi_sum_int_1s(x=x)
    return sumx

def _smpi_sum_cplx_1s(x):
    """
    sumx = _smpi_sum_cplx_1s(x)
    
    
    Defined at Smpi_math_module.fpp lines 434-436
    
    Parameters
    ----------
    x : complex
    
    Returns
    -------
    sumx : complex
    
    """
    sumx = _arespy_pkg.f90wrap_smpi_sum_cplx_1s(x=x)
    return sumx

def _smpi_sum_real_1s(x):
    """
    sumx = _smpi_sum_real_1s(x)
    
    
    Defined at Smpi_math_module.fpp lines 440-442
    
    Parameters
    ----------
    x : float
    
    Returns
    -------
    sumx : float
    
    """
    sumx = _arespy_pkg.f90wrap_smpi_sum_real_1s(x=x)
    return sumx

def _smpi_sum_real_1d(amat):
    """
    suma = _smpi_sum_real_1d(amat)
    
    
    Defined at Smpi_math_module.fpp lines 446-449
    
    Parameters
    ----------
    amat : float array
    
    Returns
    -------
    suma : float
    
    """
    suma = _arespy_pkg.f90wrap_smpi_sum_real_1d(amat=amat)
    return suma

def _smpi_sum_real_2d(amat, bmat):
    """
    suma = _smpi_sum_real_2d(amat, bmat)
    
    
    Defined at Smpi_math_module.fpp lines 453-457
    
    Parameters
    ----------
    amat : float array
    bmat : float array
    
    Returns
    -------
    suma : float
    
    """
    suma = _arespy_pkg.f90wrap_smpi_sum_real_2d(amat=amat, bmat=bmat)
    return suma

def _smpi_sum_real_3d(amat, bmat, cmat):
    """
    suma = _smpi_sum_real_3d(amat, bmat, cmat)
    
    
    Defined at Smpi_math_module.fpp lines 461-465
    
    Parameters
    ----------
    amat : float array
    bmat : float array
    cmat : float array
    
    Returns
    -------
    suma : float
    
    """
    suma = _arespy_pkg.f90wrap_smpi_sum_real_3d(amat=amat, bmat=bmat, cmat=cmat)
    return suma

def smpisum(*args, **kwargs):
    """
    smpisum(*args, **kwargs)
    
    
    Defined at Smpi_math_module.fpp lines 95-101
    
    Overloaded interface containing the following procedures:
      _smpi_sum_int_1s
      _smpi_sum_cplx_1s
      _smpi_sum_real_1s
      _smpi_sum_real_1d
      _smpi_sum_real_2d
      _smpi_sum_real_3d
    
    """
    for proc in [_smpi_sum_int_1s, _smpi_sum_cplx_1s, _smpi_sum_real_1s, \
        _smpi_sum_real_1d, _smpi_sum_real_2d, _smpi_sum_real_3d]:
        try:
            return proc(*args, **kwargs)
        except TypeError:
            continue

def _smpi_sum_mem_1d(munit, amat):
    """
    summem = _smpi_sum_mem_1d(munit, amat)
    
    
    Defined at Smpi_math_module.fpp lines 648-656
    
    Parameters
    ----------
    munit : str
    amat : float array
    
    Returns
    -------
    summem : float
    
    """
    summem = _arespy_pkg.f90wrap_smpi_sum_mem_1d(munit=munit, amat=amat)
    return summem

def _smpi_sum_mem_2d(munit, amat):
    """
    summem = _smpi_sum_mem_2d(munit, amat)
    
    
    Defined at Smpi_math_module.fpp lines 660-668
    
    Parameters
    ----------
    munit : str
    amat : float array
    
    Returns
    -------
    summem : float
    
    """
    summem = _arespy_pkg.f90wrap_smpi_sum_mem_2d(munit=munit, amat=amat)
    return summem

def _smpi_sum_mem_3d(munit, amat):
    """
    summem = _smpi_sum_mem_3d(munit, amat)
    
    
    Defined at Smpi_math_module.fpp lines 672-680
    
    Parameters
    ----------
    munit : str
    amat : float array
    
    Returns
    -------
    summem : float
    
    """
    summem = _arespy_pkg.f90wrap_smpi_sum_mem_3d(munit=munit, amat=amat)
    return summem

def smpisummem(*args, **kwargs):
    """
    smpisummem(*args, **kwargs)
    
    
    Defined at Smpi_math_module.fpp lines 103-106
    
    Overloaded interface containing the following procedures:
      _smpi_sum_mem_1d
      _smpi_sum_mem_2d
      _smpi_sum_mem_3d
    
    """
    for proc in [_smpi_sum_mem_1d, _smpi_sum_mem_2d, _smpi_sum_mem_3d]:
        try:
            return proc(*args, **kwargs)
        except TypeError:
            continue

def _smpi_reduce_sum_real_1d(amat, ramat=None):
    """
    _smpi_reduce_sum_real_1d(amat[, ramat])
    
    
    Defined at Smpi_math_module.fpp lines 499-509
    
    Parameters
    ----------
    amat : float array
    ramat : float array
    
    """
    _arespy_pkg.f90wrap_smpi_reduce_sum_real_1d(amat=amat, ramat=ramat)

def _smpi_reduce_sum_int_1d(amat, ramat=None):
    """
    _smpi_reduce_sum_int_1d(amat[, ramat])
    
    
    Defined at Smpi_math_module.fpp lines 485-495
    
    Parameters
    ----------
    amat : int array
    ramat : int array
    
    """
    _arespy_pkg.f90wrap_smpi_reduce_sum_int_1d(amat=amat, ramat=ramat)

def _smpi_reduce_sum_cplx_1d(amat, ramat=None):
    """
    _smpi_reduce_sum_cplx_1d(amat[, ramat])
    
    
    Defined at Smpi_math_module.fpp lines 526-536
    
    Parameters
    ----------
    amat : complex array
    ramat : complex array
    
    """
    _arespy_pkg.f90wrap_smpi_reduce_sum_cplx_1d(amat=amat, ramat=ramat)

def _smpi_reduce_sum_real_2d(amat, ramat=None):
    """
    _smpi_reduce_sum_real_2d(amat[, ramat])
    
    
    Defined at Smpi_math_module.fpp lines 540-550
    
    Parameters
    ----------
    amat : float array
    ramat : float array
    
    """
    _arespy_pkg.f90wrap_smpi_reduce_sum_real_2d(amat=amat, ramat=ramat)

def smpireducesum(*args, **kwargs):
    """
    smpireducesum(*args, **kwargs)
    
    
    Defined at Smpi_math_module.fpp lines 108-112
    
    Overloaded interface containing the following procedures:
      _smpi_reduce_sum_real_1d
      _smpi_reduce_sum_int_1d
      _smpi_reduce_sum_cplx_1d
      _smpi_reduce_sum_real_2d
    
    """
    for proc in [_smpi_reduce_sum_real_1d, _smpi_reduce_sum_int_1d, \
        _smpi_reduce_sum_cplx_1d, _smpi_reduce_sum_real_2d]:
        try:
            return proc(*args, **kwargs)
        except TypeError:
            continue

def _sum_pow_int(amat, pow):
    """
    totals = _sum_pow_int(amat, pow)
    
    
    Defined at Smpi_math_module.fpp lines 407-424
    
    Parameters
    ----------
    amat : float array
    pow : int
    
    Returns
    -------
    totals : float
    
    """
    totals = _arespy_pkg.f90wrap_sum_pow_int(amat=amat, pow=pow)
    return totals

def _sum_pow_real(amat, pow):
    """
    totals = _sum_pow_real(amat, pow)
    
    
    Defined at Smpi_math_module.fpp lines 386-403
    
    Parameters
    ----------
    amat : float array
    pow : float
    
    Returns
    -------
    totals : float
    
    """
    totals = _arespy_pkg.f90wrap_sum_pow_real(amat=amat, pow=pow)
    return totals

def sompsumpow(*args, **kwargs):
    """
    sompsumpow(*args, **kwargs)
    
    
    Defined at Smpi_math_module.fpp lines 114-116
    
    Overloaded interface containing the following procedures:
      _sum_pow_int
      _sum_pow_real
    
    """
    for proc in [_sum_pow_int, _sum_pow_real]:
        try:
            return proc(*args, **kwargs)
        except TypeError:
            continue

def _smpi_sum_pow_int(amat, pow):
    """
    suma = _smpi_sum_pow_int(amat, pow)
    
    
    Defined at Smpi_math_module.fpp lines 477-481
    
    Parameters
    ----------
    amat : float array
    pow : int
    
    Returns
    -------
    suma : float
    
    """
    suma = _arespy_pkg.f90wrap_smpi_sum_pow_int(amat=amat, pow=pow)
    return suma

def _smpi_sum_pow_real(amat, pow):
    """
    suma = _smpi_sum_pow_real(amat, pow)
    
    
    Defined at Smpi_math_module.fpp lines 469-473
    
    Parameters
    ----------
    amat : float array
    pow : float
    
    Returns
    -------
    suma : float
    
    """
    suma = _arespy_pkg.f90wrap_smpi_sum_pow_real(amat=amat, pow=pow)
    return suma

def smpisumpow(*args, **kwargs):
    """
    smpisumpow(*args, **kwargs)
    
    
    Defined at Smpi_math_module.fpp lines 118-119
    
    Overloaded interface containing the following procedures:
      _smpi_sum_pow_int
      _smpi_sum_pow_real
    
    """
    for proc in [_smpi_sum_pow_int, _smpi_sum_pow_real]:
        try:
            return proc(*args, **kwargs)
        except TypeError:
            continue

def _set_wrap_grid_pbc_ata_cmplx(myrho, wrap_box1d, global_n, global_n1, \
    global_n2):
    """
    _set_wrap_grid_pbc_ata_cmplx(myrho, wrap_box1d, global_n, global_n1, global_n2)
    
    
    Defined at Smpi_math_module.fpp lines 1821-1968
    
    Parameters
    ----------
    myrho : complex array
    wrap_box1d : complex array
    global_n : int
    global_n1 : int
    global_n2 : int
    
    """
    _arespy_pkg.f90wrap_set_wrap_grid_pbc_ata_cmplx(myrho=myrho, \
        wrap_box1d=wrap_box1d, global_n=global_n, global_n1=global_n1, \
        global_n2=global_n2)

def _set_wrap_grid_pbc_ata_real(myrho, wrap_box1d, global_n, global_n1, \
    global_n2):
    """
    _set_wrap_grid_pbc_ata_real(myrho, wrap_box1d, global_n, global_n1, global_n2)
    
    
    Defined at Smpi_math_module.fpp lines 1972-2044
    
    Parameters
    ----------
    myrho : float array
    wrap_box1d : float array
    global_n : int
    global_n1 : int
    global_n2 : int
    
    """
    _arespy_pkg.f90wrap_set_wrap_grid_pbc_ata_real(myrho=myrho, \
        wrap_box1d=wrap_box1d, global_n=global_n, global_n1=global_n1, \
        global_n2=global_n2)

def set_wrap_grid_pbc_ata(*args, **kwargs):
    """
    set_wrap_grid_pbc_ata(*args, **kwargs)
    
    
    Defined at Smpi_math_module.fpp lines 122-124
    
    Overloaded interface containing the following procedures:
      _set_wrap_grid_pbc_ata_cmplx
      _set_wrap_grid_pbc_ata_real
    
    """
    for proc in [_set_wrap_grid_pbc_ata_cmplx, _set_wrap_grid_pbc_ata_real]:
        try:
            return proc(*args, **kwargs)
        except TypeError:
            continue

def get_rtic():
    """
    Element rtic ftype=real(dp) pytype=float
    
    
    Defined at Smpi_math_module.fpp line 62
    
    """
    return _arespy_pkg.f90wrap_smpi_math_module__get__rtic()

def set_rtic(rtic):
    _arespy_pkg.f90wrap_smpi_math_module__set__rtic(rtic)

def get_rtoc():
    """
    Element rtoc ftype=real(dp) pytype=float
    
    
    Defined at Smpi_math_module.fpp line 62
    
    """
    return _arespy_pkg.f90wrap_smpi_math_module__get__rtoc()

def set_rtoc(rtoc):
    _arespy_pkg.f90wrap_smpi_math_module__set__rtoc(rtoc)

def get_mpinfo():
    """
    Element mpinfo ftype=integer(i4b) pytype=int
    
    
    Defined at Smpi_math_module.fpp line 83
    
    """
    return _arespy_pkg.f90wrap_smpi_math_module__get__mpinfo()

def set_mpinfo(mpinfo):
    _arespy_pkg.f90wrap_smpi_math_module__set__mpinfo(mpinfo)

def get_array_smpi_status():
    """
    Element smpi_status ftype=integer(i4b) pytype=int
    
    
    Defined at Smpi_math_module.fpp line 84
    
    """
    global smpi_status
    array_ndim, array_type, array_shape, array_handle = \
        _arespy_pkg.f90wrap_smpi_math_module__array__smpi_status(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        smpi_status = _arrays[array_handle]
    else:
        smpi_status = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _arespy_pkg.f90wrap_smpi_math_module__array__smpi_status)
        _arrays[array_handle] = smpi_status
    return smpi_status

def set_array_smpi_status(smpi_status):
    smpi_status[...] = smpi_status

def get_lall_grid():
    """
    Element lall_grid ftype=logical pytype=bool
    
    
    Defined at Smpi_math_module.fpp line 85
    
    """
    return _arespy_pkg.f90wrap_smpi_math_module__get__lall_grid()

def set_lall_grid(lall_grid):
    _arespy_pkg.f90wrap_smpi_math_module__set__lall_grid(lall_grid)


_array_initialisers = [get_array_smpi_status]
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "smpi_math_module".')

for func in _dt_array_initialisers:
    func()

"""
Module mixer_module


Defined at Mixer_module.fpp lines 5-651

"""
from __future__ import print_function, absolute_import, division
import _arespy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("arespy_pkg.mixer_data")
class mixer_data(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=mixer_data)
    
    
    Defined at Mixer_module.fpp lines 18-27
    
    """
    def __init__(self, handle=None):
        """
        self = Mixer_Data()
        
        
        Defined at Mixer_module.fpp lines 18-27
        
        
        Returns
        -------
        this : Mixer_Data
        	Object to be constructed
        
        
        Automatically generated constructor for mixer_data
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _arespy_pkg.f90wrap_mixer_data_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Mixer_Data
        
        
        Defined at Mixer_module.fpp lines 18-27
        
        Parameters
        ----------
        this : Mixer_Data
        	Object to be destructed
        
        
        Automatically generated destructor for mixer_data
        """
        if self._alloc:
            _arespy_pkg.f90wrap_mixer_data_finalise(this=self._handle)
    
    @property
    def dxl(self):
        """
        Element dxl ftype=real(dp) pytype=float
        
        
        Defined at Mixer_module.fpp line 22
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_mixer_data__array__dxl(self._handle)
        if array_handle in self._arrays:
            dxl = self._arrays[array_handle]
        else:
            dxl = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_mixer_data__array__dxl)
            self._arrays[array_handle] = dxl
        return dxl
    
    @dxl.setter
    def dxl(self, dxl):
        self.dxl[...] = dxl
    
    @property
    def dfl(self):
        """
        Element dfl ftype=real(dp) pytype=float
        
        
        Defined at Mixer_module.fpp line 22
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_mixer_data__array__dfl(self._handle)
        if array_handle in self._arrays:
            dfl = self._arrays[array_handle]
        else:
            dfl = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_mixer_data__array__dfl)
            self._arrays[array_handle] = dfl
        return dfl
    
    @dfl.setter
    def dfl(self, dfl):
        self.dfl[...] = dfl
    
    @property
    def voma(self):
        """
        Element voma ftype=real(dp) pytype=float
        
        
        Defined at Mixer_module.fpp line 22
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_mixer_data__array__voma(self._handle)
        if array_handle in self._arrays:
            voma = self._arrays[array_handle]
        else:
            voma = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_mixer_data__array__voma)
            self._arrays[array_handle] = voma
        return voma
    
    @voma.setter
    def voma(self, voma):
        self.voma[...] = voma
    
    @property
    def kerker(self):
        """
        Element kerker ftype=real(dp) pytype=float
        
        
        Defined at Mixer_module.fpp line 24
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_mixer_data__array__kerker(self._handle)
        if array_handle in self._arrays:
            kerker = self._arrays[array_handle]
        else:
            kerker = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_mixer_data__array__kerker)
            self._arrays[array_handle] = kerker
        return kerker
    
    @kerker.setter
    def kerker(self, kerker):
        self.kerker[...] = kerker
    
    @property
    def dxgl(self):
        """
        Element dxgl ftype=complex(dcp) pytype=complex
        
        
        Defined at Mixer_module.fpp line 27
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_mixer_data__array__dxgl(self._handle)
        if array_handle in self._arrays:
            dxgl = self._arrays[array_handle]
        else:
            dxgl = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_mixer_data__array__dxgl)
            self._arrays[array_handle] = dxgl
        return dxgl
    
    @dxgl.setter
    def dxgl(self, dxgl):
        self.dxgl[...] = dxgl
    
    @property
    def drgl(self):
        """
        Element drgl ftype=complex(dcp) pytype=complex
        
        
        Defined at Mixer_module.fpp line 27
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _arespy_pkg.f90wrap_mixer_data__array__drgl(self._handle)
        if array_handle in self._arrays:
            drgl = self._arrays[array_handle]
        else:
            drgl = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _arespy_pkg.f90wrap_mixer_data__array__drgl)
            self._arrays[array_handle] = drgl
        return drgl
    
    @drgl.setter
    def drgl(self, drgl):
        self.drgl[...] = drgl
    
    def __str__(self):
        ret = ['<mixer_data>{\n']
        ret.append('    dxl : ')
        ret.append(repr(self.dxl))
        ret.append(',\n    dfl : ')
        ret.append(repr(self.dfl))
        ret.append(',\n    voma : ')
        ret.append(repr(self.voma))
        ret.append(',\n    kerker : ')
        ret.append(repr(self.kerker))
        ret.append(',\n    dxgl : ')
        ret.append(repr(self.dxgl))
        ret.append(',\n    drgl : ')
        ret.append(repr(self.drgl))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def init_mixer(nps):
    """
    init_mixer(nps)
    
    
    Defined at Mixer_module.fpp lines 35-75
    
    Parameters
    ----------
    nps : int
    
    """
    _arespy_pkg.f90wrap_init_mixer(nps=nps)

def destroy_mixer():
    """
    destroy_mixer()
    
    
    Defined at Mixer_module.fpp lines 78-89
    
    
    """
    _arespy_pkg.f90wrap_destroy_mixer()

def mixing(iter, xout, xin):
    """
    res = mixing(iter, xout, xin)
    
    
    Defined at Mixer_module.fpp lines 92-165
    
    Parameters
    ----------
    iter : int
    xout : float array
    xin : float array
    
    Returns
    -------
    res : float
    
    """
    res = _arespy_pkg.f90wrap_mixing(iter=iter, xout=xout, xin=xin)
    return res

def anderson_mixing(iiter, xin, xout):
    """
    err = anderson_mixing(iiter, xin, xout)
    
    
    Defined at Mixer_module.fpp lines 171-239
    
    Parameters
    ----------
    iiter : int
    xin : float array
    xout : float array
    
    Returns
    -------
    err : float
    
    """
    err = _arespy_pkg.f90wrap_anderson_mixing(iiter=iiter, xin=xin, xout=xout)
    return err

def om1c(nam, nuh, sp, dfp, voma):
    """
    om1c(nam, nuh, sp, dfp, voma)
    
    
    Defined at Mixer_module.fpp lines 242-288
    
    Parameters
    ----------
    nam : int
    nuh : int
    sp : float
    dfp : float
    voma : float
    
    """
    _arespy_pkg.f90wrap_om1c(nam=nam, nuh=nuh, sp=sp, dfp=dfp, voma=voma)

def amst(beta, w0, nam, nuh, dxp, dfp, sp, xl, fl, voma, xn):
    """
    amst(beta, w0, nam, nuh, dxp, dfp, sp, xl, fl, voma, xn)
    
    
    Defined at Mixer_module.fpp lines 291-394
    
    Parameters
    ----------
    beta : float
    w0 : float
    nam : int
    nuh : int
    dxp : float
    dfp : float
    sp : float
    xl : float
    fl : float
    voma : float
    xn : float
    
    """
    _arespy_pkg.f90wrap_amst(beta=beta, w0=w0, nam=nam, nuh=nuh, dxp=dxp, dfp=dfp, \
        sp=sp, xl=xl, fl=fl, voma=voma, xn=xn)

def init_kerker():
    """
    init_kerker()
    
    
    Defined at Mixer_module.fpp lines 397-413
    
    
    """
    _arespy_pkg.f90wrap_init_kerker()

def rpulayk_mixing(iter, rlg, xing):
    """
    rpulayk_mixing(iter, rlg, xing)
    
    
    Defined at Mixer_module.fpp lines 419-481
    
    Parameters
    ----------
    iter : int
    rlg : complex array
    xing : complex array
    
    ----------------------
    store RL
    """
    _arespy_pkg.f90wrap_rpulayk_mixing(iter=iter, rlg=rlg, xing=xing)

def rpulayk_mix(beta, w0, dime, nh, dxl, drl, xl, rl, xn):
    """
    rpulayk_mix(beta, w0, dime, nh, dxl, drl, xl, rl, xn)
    
    
    Defined at Mixer_module.fpp lines 484-542
    
    Parameters
    ----------
    beta : float
    w0 : float
    dime : int
    nh : int
    dxl : complex array
    drl : complex array
    xl : complex array
    rl : complex array
    xn : complex array
    
    """
    _arespy_pkg.f90wrap_rpulayk_mix(beta=beta, w0=w0, dime=dime, nh=nh, dxl=dxl, \
        drl=drl, xl=xl, rl=rl, xn=xn)

def resta_mixing(iter, rlg, xing):
    """
    resta_mixing(iter, rlg, xing)
    
    
    Defined at Mixer_module.fpp lines 547-611
    
    Parameters
    ----------
    iter : int
    rlg : complex array
    xing : complex array
    
    ----------------------
    store RL
    """
    _arespy_pkg.f90wrap_resta_mixing(iter=iter, rlg=rlg, xing=xing)

def init_resta():
    """
    init_resta()
    
    
    Defined at Mixer_module.fpp lines 614-648
    
    
    """
    _arespy_pkg.f90wrap_init_resta()

def get_nam():
    """
    Element nam ftype=integer(i4b) pytype=int
    
    
    Defined at Mixer_module.fpp line 30
    
    """
    return _arespy_pkg.f90wrap_mixer_module__get__nam()

def set_nam(nam):
    _arespy_pkg.f90wrap_mixer_module__set__nam(nam)

def get_nuh():
    """
    Element nuh ftype=integer(i4b) pytype=int
    
    
    Defined at Mixer_module.fpp line 30
    
    """
    return _arespy_pkg.f90wrap_mixer_module__get__nuh()

def set_nuh(nuh):
    _arespy_pkg.f90wrap_mixer_module__set__nuh(nuh)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "mixer_module".')

for func in _dt_array_initialisers:
    func()

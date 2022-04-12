"""
Module mixer_module


Defined at Mixer_module.fpp lines 5-1148

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("AresMainPy_pkg.mixer_data")
class mixer_data(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=mixer_data)
    
    
    Defined at Mixer_module.fpp lines 16-37
    
    """
    def __init__(self, handle=None):
        """
        self = Mixer_Data()
        
        
        Defined at Mixer_module.fpp lines 16-37
        
        
        Returns
        -------
        this : Mixer_Data
        	Object to be constructed
        
        
        Automatically generated constructor for mixer_data
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _AresMainPy_pkg.f90wrap_mixer_data_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Mixer_Data
        
        
        Defined at Mixer_module.fpp lines 16-37
        
        Parameters
        ----------
        this : Mixer_Data
        	Object to be destructed
        
        
        Automatically generated destructor for mixer_data
        """
        if self._alloc:
            _AresMainPy_pkg.f90wrap_mixer_data_finalise(this=self._handle)
    
    @property
    def nhmix(self):
        """
        Element nhmix ftype=integer  pytype=int
        
        
        Defined at Mixer_module.fpp line 21
        
        """
        return _AresMainPy_pkg.f90wrap_mixer_data__get__nhmix(self._handle)
    
    @nhmix.setter
    def nhmix(self, nhmix):
        _AresMainPy_pkg.f90wrap_mixer_data__set__nhmix(self._handle, nhmix)
    
    @property
    def nhmin(self):
        """
        Element nhmin ftype=integer  pytype=int
        
        
        Defined at Mixer_module.fpp line 21
        
        """
        return _AresMainPy_pkg.f90wrap_mixer_data__get__nhmin(self._handle)
    
    @nhmin.setter
    def nhmin(self, nhmin):
        _AresMainPy_pkg.f90wrap_mixer_data__set__nhmin(self._handle, nhmin)
    
    @property
    def nam(self):
        """
        Element nam ftype=integer  pytype=int
        
        
        Defined at Mixer_module.fpp line 21
        
        """
        return _AresMainPy_pkg.f90wrap_mixer_data__get__nam(self._handle)
    
    @nam.setter
    def nam(self, nam):
        _AresMainPy_pkg.f90wrap_mixer_data__set__nam(self._handle, nam)
    
    @property
    def nitra(self):
        """
        Element nitra ftype=integer  pytype=int
        
        
        Defined at Mixer_module.fpp line 21
        
        """
        return _AresMainPy_pkg.f90wrap_mixer_data__get__nitra(self._handle)
    
    @nitra.setter
    def nitra(self, nitra):
        _AresMainPy_pkg.f90wrap_mixer_data__set__nitra(self._handle, nitra)
    
    @property
    def alpha(self):
        """
        Element alpha ftype=real(dp) pytype=float
        
        
        Defined at Mixer_module.fpp line 25
        
        """
        return _AresMainPy_pkg.f90wrap_mixer_data__get__alpha(self._handle)
    
    @alpha.setter
    def alpha(self, alpha):
        _AresMainPy_pkg.f90wrap_mixer_data__set__alpha(self._handle, alpha)
    
    @property
    def beta(self):
        """
        Element beta ftype=real(dp) pytype=float
        
        
        Defined at Mixer_module.fpp line 25
        
        """
        return _AresMainPy_pkg.f90wrap_mixer_data__get__beta(self._handle)
    
    @beta.setter
    def beta(self, beta):
        _AresMainPy_pkg.f90wrap_mixer_data__set__beta(self._handle, beta)
    
    @property
    def w0(self):
        """
        Element w0 ftype=real(dp) pytype=float
        
        
        Defined at Mixer_module.fpp line 25
        
        """
        return _AresMainPy_pkg.f90wrap_mixer_data__get__w0(self._handle)
    
    @w0.setter
    def w0(self, w0):
        _AresMainPy_pkg.f90wrap_mixer_data__set__w0(self._handle, w0)
    
    @property
    def dxl(self):
        """
        Element dxl ftype=real(dp) pytype=float
        
        
        Defined at Mixer_module.fpp line 29
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_mixer_data__array__dxl(self._handle)
        if array_handle in self._arrays:
            dxl = self._arrays[array_handle]
        else:
            dxl = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_mixer_data__array__dxl)
            self._arrays[array_handle] = dxl
        return dxl
    
    @dxl.setter
    def dxl(self, dxl):
        self.dxl[...] = dxl
    
    @property
    def dfl(self):
        """
        Element dfl ftype=real(dp) pytype=float
        
        
        Defined at Mixer_module.fpp line 29
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_mixer_data__array__dfl(self._handle)
        if array_handle in self._arrays:
            dfl = self._arrays[array_handle]
        else:
            dfl = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_mixer_data__array__dfl)
            self._arrays[array_handle] = dfl
        return dfl
    
    @dfl.setter
    def dfl(self, dfl):
        self.dfl[...] = dfl
    
    @property
    def voma(self):
        """
        Element voma ftype=real(dp) pytype=float
        
        
        Defined at Mixer_module.fpp line 29
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_mixer_data__array__voma(self._handle)
        if array_handle in self._arrays:
            voma = self._arrays[array_handle]
        else:
            voma = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_mixer_data__array__voma)
            self._arrays[array_handle] = voma
        return voma
    
    @voma.setter
    def voma(self, voma):
        self.voma[...] = voma
    
    @property
    def sp(self):
        """
        Element sp ftype=real(dp) pytype=float
        
        
        Defined at Mixer_module.fpp line 30
        
        """
        return _AresMainPy_pkg.f90wrap_mixer_data__get__sp(self._handle)
    
    @sp.setter
    def sp(self, sp):
        _AresMainPy_pkg.f90wrap_mixer_data__set__sp(self._handle, sp)
    
    @property
    def kerker(self):
        """
        Element kerker ftype=real(dp) pytype=float
        
        
        Defined at Mixer_module.fpp line 32
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_mixer_data__array__kerker(self._handle)
        if array_handle in self._arrays:
            kerker = self._arrays[array_handle]
        else:
            kerker = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_mixer_data__array__kerker)
            self._arrays[array_handle] = kerker
        return kerker
    
    @kerker.setter
    def kerker(self, kerker):
        self.kerker[...] = kerker
    
    @property
    def dxgl(self):
        """
        Element dxgl ftype=complex(dp) pytype=complex
        
        
        Defined at Mixer_module.fpp line 35
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_mixer_data__array__dxgl(self._handle)
        if array_handle in self._arrays:
            dxgl = self._arrays[array_handle]
        else:
            dxgl = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_mixer_data__array__dxgl)
            self._arrays[array_handle] = dxgl
        return dxgl
    
    @dxgl.setter
    def dxgl(self, dxgl):
        self.dxgl[...] = dxgl
    
    @property
    def dfgl(self):
        """
        Element dfgl ftype=complex(dp) pytype=complex
        
        
        Defined at Mixer_module.fpp line 35
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_mixer_data__array__dfgl(self._handle)
        if array_handle in self._arrays:
            dfgl = self._arrays[array_handle]
        else:
            dfgl = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_mixer_data__array__dfgl)
            self._arrays[array_handle] = dfgl
        return dfgl
    
    @dfgl.setter
    def dfgl(self, dfgl):
        self.dfgl[...] = dfgl
    
    @property
    def bmix(self):
        """
        Element bmix ftype=real(dp) pytype=float
        
        
        Defined at Mixer_module.fpp line 37
        
        """
        return _AresMainPy_pkg.f90wrap_mixer_data__get__bmix(self._handle)
    
    @bmix.setter
    def bmix(self, bmix):
        _AresMainPy_pkg.f90wrap_mixer_data__set__bmix(self._handle, bmix)
    
    @property
    def amix(self):
        """
        Element amix ftype=real(dp) pytype=float
        
        
        Defined at Mixer_module.fpp line 37
        
        """
        return _AresMainPy_pkg.f90wrap_mixer_data__get__amix(self._handle)
    
    @amix.setter
    def amix(self, amix):
        _AresMainPy_pkg.f90wrap_mixer_data__set__amix(self._handle, amix)
    
    def __str__(self):
        ret = ['<mixer_data>{\n']
        ret.append('    nhmix : ')
        ret.append(repr(self.nhmix))
        ret.append(',\n    nhmin : ')
        ret.append(repr(self.nhmin))
        ret.append(',\n    nam : ')
        ret.append(repr(self.nam))
        ret.append(',\n    nitra : ')
        ret.append(repr(self.nitra))
        ret.append(',\n    alpha : ')
        ret.append(repr(self.alpha))
        ret.append(',\n    beta : ')
        ret.append(repr(self.beta))
        ret.append(',\n    w0 : ')
        ret.append(repr(self.w0))
        ret.append(',\n    dxl : ')
        ret.append(repr(self.dxl))
        ret.append(',\n    dfl : ')
        ret.append(repr(self.dfl))
        ret.append(',\n    voma : ')
        ret.append(repr(self.voma))
        ret.append(',\n    sp : ')
        ret.append(repr(self.sp))
        ret.append(',\n    kerker : ')
        ret.append(repr(self.kerker))
        ret.append(',\n    dxgl : ')
        ret.append(repr(self.dxgl))
        ret.append(',\n    dfgl : ')
        ret.append(repr(self.dfgl))
        ret.append(',\n    bmix : ')
        ret.append(repr(self.bmix))
        ret.append(',\n    amix : ')
        ret.append(repr(self.amix))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def init_mixer_data():
    """
    init_mixer_data()
    
    
    Defined at Mixer_module.fpp lines 44-51
    
    
    """
    _AresMainPy_pkg.f90wrap_init_mixer_data()

def init_mixer_data_per():
    """
    init_mixer_data_per()
    
    
    Defined at Mixer_module.fpp lines 54-129
    
    
    """
    _AresMainPy_pkg.f90wrap_init_mixer_data_per()

def init_mixer_data_iso():
    """
    init_mixer_data_iso()
    
    
    Defined at Mixer_module.fpp lines 133-200
    
    
    """
    _AresMainPy_pkg.f90wrap_init_mixer_data_iso()

def destroy_mixer():
    """
    destroy_mixer()
    
    
    Defined at Mixer_module.fpp lines 204-234
    
    
    """
    _AresMainPy_pkg.f90wrap_destroy_mixer()

def mixing(iter, xout, xin):
    """
    res = mixing(iter, xout, xin)
    
    
    Defined at Mixer_module.fpp lines 238-321
    
    Parameters
    ----------
    iter : int
    xout : float array
    xin : float array
    
    Returns
    -------
    res : float
    
    """
    res = _AresMainPy_pkg.f90wrap_mixing(iter=iter, xout=xout, xin=xin)
    return res

def mixing_iso(iter, xout, xin):
    """
    res = mixing_iso(iter, xout, xin)
    
    
    Defined at Mixer_module.fpp lines 325-382
    
    Parameters
    ----------
    iter : int
    xout : float array
    xin : float array
    
    Returns
    -------
    res : float
    
    """
    res = _AresMainPy_pkg.f90wrap_mixing_iso(iter=iter, xout=xout, xin=xin)
    return res

def anderson_mixing(iiter, xout, xin):
    """
    err = anderson_mixing(iiter, xout, xin)
    
    
    Defined at Mixer_module.fpp lines 389-512
    
    Parameters
    ----------
    iiter : int
    xout : float array
    xin : float array
    
    Returns
    -------
    err : float
    
    """
    err = _AresMainPy_pkg.f90wrap_anderson_mixing(iiter=iiter, xout=xout, xin=xin)
    return err

def om1c(nam, nuh, sp, dfp, voma):
    """
    om1c(nam, nuh, sp, dfp, voma)
    
    
    Defined at Mixer_module.fpp lines 516-571
    
    Parameters
    ----------
    nam : int
    nuh : int
    sp : float
    dfp : float array
    voma : float array
    
    """
    _AresMainPy_pkg.f90wrap_om1c(nam=nam, nuh=nuh, sp=sp, dfp=dfp, voma=voma)

def amst(beta, w0, nam, nuh, dxp, dfp, sp, xl, fl, voma, xn):
    """
    amst(beta, w0, nam, nuh, dxp, dfp, sp, xl, fl, voma, xn)
    
    
    Defined at Mixer_module.fpp lines 575-685
    
    Parameters
    ----------
    beta : float
    w0 : float
    nam : int
    nuh : int
    dxp : float array
    dfp : float array
    sp : float
    xl : float array
    fl : float array
    voma : float array
    xn : float array
    
    """
    _AresMainPy_pkg.f90wrap_amst(beta=beta, w0=w0, nam=nam, nuh=nuh, dxp=dxp, \
        dfp=dfp, sp=sp, xl=xl, fl=fl, voma=voma, xn=xn)

def rpulay_mixing(iter, xout, xin):
    """
    err = rpulay_mixing(iter, xout, xin)
    
    
    Defined at Mixer_module.fpp lines 693-757
    
    Parameters
    ----------
    iter : int
    xout : float array
    xin : float array
    
    Returns
    -------
    err : float
    
    ----------------------
    store RL
    """
    err = _AresMainPy_pkg.f90wrap_rpulay_mixing(iter=iter, xout=xout, xin=xin)
    return err

def rpulay_mix(beta, w0, dime, nh, dxl, drl, xl, rl, xn):
    """
    rpulay_mix(beta, w0, dime, nh, dxl, drl, xl, rl, xn)
    
    
    Defined at Mixer_module.fpp lines 761-799
    
    Parameters
    ----------
    beta : float
    w0 : float
    dime : int
    nh : int
    dxl : float array
    drl : float array
    xl : float array
    rl : float array
    xn : float array
    
    """
    _AresMainPy_pkg.f90wrap_rpulay_mix(beta=beta, w0=w0, dime=dime, nh=nh, dxl=dxl, \
        drl=drl, xl=xl, rl=rl, xn=xn)

def init_kerker():
    """
    init_kerker()
    
    
    Defined at Mixer_module.fpp lines 892-909
    
    
    """
    _AresMainPy_pkg.f90wrap_init_kerker()

def rpulayk_mixing(iter, rlg, xing):
    """
    rpulayk_mixing(iter, rlg, xing)
    
    
    Defined at Mixer_module.fpp lines 917-978
    
    Parameters
    ----------
    iter : int
    rlg : complex array
    xing : complex array
    
    ----------------------
    store RL
    """
    _AresMainPy_pkg.f90wrap_rpulayk_mixing(iter=iter, rlg=rlg, xing=xing)

def rpulayk_mix(beta, w0, dime, nh, dxl, drl, xl, rl, xn):
    """
    rpulayk_mix(beta, w0, dime, nh, dxl, drl, xl, rl, xn)
    
    
    Defined at Mixer_module.fpp lines 982-1037
    
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
    _AresMainPy_pkg.f90wrap_rpulayk_mix(beta=beta, w0=w0, dime=dime, nh=nh, dxl=dxl, \
        drl=drl, xl=xl, rl=rl, xn=xn)

def resta_mixing(iter, rlg, xing):
    """
    resta_mixing(iter, rlg, xing)
    
    
    Defined at Mixer_module.fpp lines 1043-1107
    
    Parameters
    ----------
    iter : int
    rlg : complex array
    xing : complex array
    
    ----------------------
    store RL
    """
    _AresMainPy_pkg.f90wrap_resta_mixing(iter=iter, rlg=rlg, xing=xing)

def init_resta():
    """
    init_resta()
    
    
    Defined at Mixer_module.fpp lines 1110-1145
    
    
    """
    _AresMainPy_pkg.f90wrap_init_resta()


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

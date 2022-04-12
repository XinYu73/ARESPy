"""
Module fourier


Defined at Fourier.fpp lines 5-497

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def planfft(dimx, dimy, dimz):
    """
    planfft(dimx, dimy, dimz)
    
    
    Defined at Fourier.fpp lines 73-147
    
    Parameters
    ----------
    dimx : int
    dimy : int
    dimz : int
    
    ------------------------------------------------------------------------------
     DESCRIPTION:
       This is the initialization procedure that first gets the system name as is
       called as an argument to OFDFT, and turns it into the various input file
       names.  Then, it calls all the programs necessary to set variables to
       default values, then reads the geometry file to get all the variables sets
       to the correct values.
     GLOBAL/MODULE VARIABLES CHANGED:
       realRA, cplxRA, offset
     CONDITIONS AND ASSUMPTIONS:
     FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
     REFERENCES:
    ------------------------------------------------------------------------------
     REVISION LOG:
       11/20/2003  File created.  (VLL)
    ------------------------------------------------------------------------------
    """
    _AresMainPy_pkg.f90wrap_planfft(dimx=dimx, dimy=dimy, dimz=dimz)

def planfst(dimx, dimy, dimz):
    """
    planfst(dimx, dimy, dimz)
    
    
    Defined at Fourier.fpp lines 149-181
    
    Parameters
    ----------
    dimx : int
    dimy : int
    dimz : int
    
    ------------------------------------------------------------------------------
     DESCRIPTION:
       This is the same as PlanFFT, except for the Fast Sine Transform.
     CONDITIONS AND ASSUMPTIONS:
     FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
     REFERENCES:
    ------------------------------------------------------------------------------
     REVISION LOG:
       11/20/2003  File created.  (VLL)
    ------------------------------------------------------------------------------
    """
    _AresMainPy_pkg.f90wrap_planfst(dimx=dimx, dimy=dimy, dimz=dimz)

def getfftdims():
    """
    dimx, dimy, dimz = getfftdims()
    
    
    Defined at Fourier.fpp lines 183-210
    
    
    Returns
    -------
    dimx : int
    dimy : int
    dimz : int
    
    ------------------------------------------------------------------------------
     DESCRIPTION:
       Gets the dimensions of the FFT(real-space part)
     GLOBAL/MODULE VARIABLES CHANGED:
     CONDITIONS AND ASSUMPTIONS:
     FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
     REFERENCES:
    ------------------------------------------------------------------------------
     REVISION LOG:
       4/25/2006  Added(GSH)
    ------------------------------------------------------------------------------
    """
    dimx, dimy, dimz = _AresMainPy_pkg.f90wrap_getfftdims()
    return dimx, dimy, dimz

def getfftcomplexdims():
    """
    dimx, dimy, dimz = getfftcomplexdims()
    
    
    Defined at Fourier.fpp lines 212-239
    
    
    Returns
    -------
    dimx : int
    dimy : int
    dimz : int
    
    ------------------------------------------------------------------------------
     DESCRIPTION:
       Gets the dimensions of the FFT(reciprocal space part)
     GLOBAL/MODULE VARIABLES CHANGED:
     CONDITIONS AND ASSUMPTIONS:
     FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
     REFERENCES:
    ------------------------------------------------------------------------------
     REVISION LOG:
       4/26/2006 Added(GSH)
    ------------------------------------------------------------------------------
    """
    dimx, dimy, dimz = _AresMainPy_pkg.f90wrap_getfftcomplexdims()
    return dimx, dimy, dimz

def forwardfst(array):
    """
    forwardfst(array)
    
    
    Defined at Fourier.fpp lines 410-435
    
    Parameters
    ----------
    array : float array
    
    ------------------------------------------------------------------------------
     DESCRIPTION:
     GLOBAL/MODULE VARIABLES CHANGED:
     CONDITIONS AND ASSUMPTIONS:
     FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
     REFERENCES:
    ------------------------------------------------------------------------------
     REVISION LOG:
    ------------------------------------------------------------------------------
    """
    _AresMainPy_pkg.f90wrap_forwardfst(array=array)

def backfst(array):
    """
    backfst(array)
    
    
    Defined at Fourier.fpp lines 437-460
    
    Parameters
    ----------
    array : float array
    
    ------------------------------------------------------------------------------
     DESCRIPTION:
     GLOBAL/MODULE VARIABLES CHANGED:
     CONDITIONS AND ASSUMPTIONS:
     FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
     REFERENCES:
    ------------------------------------------------------------------------------
     REVISION LOG:
    ------------------------------------------------------------------------------
    """
    _AresMainPy_pkg.f90wrap_backfst(array=array)

def cleanfft():
    """
    cleanfft()
    
    
    Defined at Fourier.fpp lines 462-497
    
    
    ------------------------------------------------------------------------------
     DESCRIPTION:
       This subroutine is called at the end of the run to free the memory
       associated with the plan.
     GLOBAL/MODULE VARIABLES CHANGED:
       realRA, cplxRA
     CONDITIONS AND ASSUMPTIONS:
     FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
     REFERENCES:
    ------------------------------------------------------------------------------
     REVISION LOG:
       11/20/2003  File created.  (VLL)
    ------------------------------------------------------------------------------
    """
    _AresMainPy_pkg.f90wrap_cleanfft()

def _forwardfft_4d(array):
    """
    transform = _forwardfft_4d(array)
    
    
    Defined at Fourier.fpp lines 241-284
    
    Parameters
    ----------
    array : float array
    
    Returns
    -------
    transform : complex array
    
    ------------------------------------------------------------------------------
     DESCRIPTION:
       This function is not called directly from the OFDFT code. Use the FFT
       interface instead. It performs the transformation of a real 4-dimensional
       array into its complex 4-dimensional transform. The first dimension is
       halved.
     GLOBAL/MODULE VARIABLES CHANGED:
     CONDITIONS AND ASSUMPTIONS:
     FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
     REFERENCES:
    ------------------------------------------------------------------------------
     REVISION LOG:
       11/20/2003  File created.  (VLL)
    ------------------------------------------------------------------------------
    """
    transform = _AresMainPy_pkg.f90wrap_forwardfft_4d(array=array)
    return transform

def _backfft_4d(array):
    """
    transform = _backfft_4d(array)
    
    
    Defined at Fourier.fpp lines 286-327
    
    Parameters
    ----------
    array : complex array
    
    Returns
    -------
    transform : float array
    
    ------------------------------------------------------------------------------
     DESCRIPTION:
       This function is not called directly from the OFDFT code, but rather
       through the FFT interface. It performs the reverse Fourier transform of
       a complex function over the half-box in reciprocal space back to real
       space. It acts on 4-dimensional arrays, the fourth dimension being spin.
     GLOBAL/MODULE VARIABLES CHANGED:
     CONDITIONS AND ASSUMPTIONS:
     FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
     REFERENCES:
    ------------------------------------------------------------------------------
     REVISION LOG:
       11/20/2003    File created.  (VLL)
    ------------------------------------------------------------------------------
    """
    transform = _AresMainPy_pkg.f90wrap_backfft_4d(array=array)
    return transform

def _forwardfft_3d(array):
    """
    transform = _forwardfft_3d(array)
    
    
    Defined at Fourier.fpp lines 329-367
    
    Parameters
    ----------
    array : float array
    
    Returns
    -------
    transform : complex array
    
    ------------------------------------------------------------------------------
     DESCRIPTION:
       This function is not called directly from the OFDFT code. Use the FFT
       interface instead. It performs the transformation of a real 4-dimensional
       array into its complex 4-dimensional transform. The first dimension is
       halved.
     GLOBAL/MODULE VARIABLES CHANGED:
     CONDITIONS AND ASSUMPTIONS:
     FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
     REFERENCES:
    ------------------------------------------------------------------------------
     REVISION LOG:
       11/20/2003  File created.  (VLL)
    ------------------------------------------------------------------------------
    """
    transform = _AresMainPy_pkg.f90wrap_forwardfft_3d(array=array)
    return transform

def _backfft_3d(array):
    """
    transform = _backfft_3d(array)
    
    
    Defined at Fourier.fpp lines 369-408
    
    Parameters
    ----------
    array : complex array
    
    Returns
    -------
    transform : float array
    
    ------------------------------------------------------------------------------
     DESCRIPTION:
       This function is not called directly from the OFDFT code, but rather
       through the FFT interface. It performs the reverse Fourier transform of a
       complex function over the half-box in reciprocal space back to real
       space. It acts on 3-dimensional arrays.
     GLOBAL/MODULE VARIABLES CHANGED:
     CONDITIONS AND ASSUMPTIONS:
     FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
     REFERENCES:
    ------------------------------------------------------------------------------
     REVISION LOG:
       11/20/2003  File created.  (VLL)
    ------------------------------------------------------------------------------
    """
    transform = _AresMainPy_pkg.f90wrap_backfft_3d(array=array)
    return transform

def fft(*args, **kwargs):
    """
    fft(*args, **kwargs)
    
    
    Defined at Fourier.fpp lines 66-70
    
    Overloaded interface containing the following procedures:
      _forwardfft_4d
      _backfft_4d
      _forwardfft_3d
      _backfft_3d
    
    """
    for proc in [_forwardfft_4d, _backfft_4d, _forwardfft_3d, _backfft_3d]:
        try:
            return proc(*args, **kwargs)
        except TypeError:
            continue

def get_offset():
    """
    Element offset ftype=integer(i4b) pytype=int
    
    
    Defined at Fourier.fpp line 61
    
    """
    return _AresMainPy_pkg.f90wrap_fourier__get__offset()

def set_offset(offset):
    _AresMainPy_pkg.f90wrap_fourier__set__offset(offset)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "fourier".')

for func in _dt_array_initialisers:
    func()

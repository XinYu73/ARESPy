"""
Module array_io


Defined at array_io.fpp lines 5-96

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("AresMainPy_pkg.out_label")
class out_label(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=out_label)
    
    
    Defined at array_io.fpp lines 8-10
    
    """
    def __init__(self, handle=None):
        """
        self = Out_Label()
        
        
        Defined at array_io.fpp lines 8-10
        
        
        Returns
        -------
        this : Out_Label
        	Object to be constructed
        
        
        Automatically generated constructor for out_label
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _AresMainPy_pkg.f90wrap_out_label_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Out_Label
        
        
        Defined at array_io.fpp lines 8-10
        
        Parameters
        ----------
        this : Out_Label
        	Object to be destructed
        
        
        Automatically generated destructor for out_label
        """
        if self._alloc:
            _AresMainPy_pkg.f90wrap_out_label_finalise(this=self._handle)
    
    @property
    def label(self):
        """
        Element label ftype=character(len=100) pytype=str
        
        
        Defined at array_io.fpp line 9
        
        """
        return _AresMainPy_pkg.f90wrap_out_label__get__label(self._handle)
    
    @label.setter
    def label(self, label):
        _AresMainPy_pkg.f90wrap_out_label__set__label(self._handle, label)
    
    @property
    def num(self):
        """
        Element num ftype=integer(i4b) pytype=int
        
        
        Defined at array_io.fpp line 10
        
        """
        return _AresMainPy_pkg.f90wrap_out_label__get__num(self._handle)
    
    @num.setter
    def num(self, num):
        _AresMainPy_pkg.f90wrap_out_label__set__num(self._handle, num)
    
    def __str__(self):
        ret = ['<out_label>{\n']
        ret.append('    label : ')
        ret.append(repr(self.label))
        ret.append(',\n    num : ')
        ret.append(repr(self.num))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def init_outfile():
    """
    init_outfile()
    
    
    Defined at array_io.fpp lines 21-27
    
    
    """
    _AresMainPy_pkg.f90wrap_init_outfile()

def _output_r(size_bn, array, name):
    """
    _output_r(size_bn, array, name)
    
    
    Defined at array_io.fpp lines 29-68
    
    Parameters
    ----------
    size_bn : int
    array : float array
    name : str
    
    """
    _AresMainPy_pkg.f90wrap_output_r(size_bn=size_bn, array=array, name=name)

def _output_i(size_bn, array, name):
    """
    _output_i(size_bn, array, name)
    
    
    Defined at array_io.fpp lines 79-86
    
    Parameters
    ----------
    size_bn : int
    array : int array
    name : str
    
    """
    _AresMainPy_pkg.f90wrap_output_i(size_bn=size_bn, array=array, name=name)

def output(*args, **kwargs):
    """
    output(*args, **kwargs)
    
    
    Defined at array_io.fpp lines 14-15
    
    Overloaded interface containing the following procedures:
      _output_r
      _output_i
    
    """
    for proc in [_output_r, _output_i]:
        try:
            return proc(*args, **kwargs)
        except TypeError:
            continue

def _input_r(size_bn, array, name):
    """
    _input_r(size_bn, array, name)
    
    
    Defined at array_io.fpp lines 70-77
    
    Parameters
    ----------
    size_bn : int
    array : float array
    name : str
    
    """
    _AresMainPy_pkg.f90wrap_input_r(size_bn=size_bn, array=array, name=name)

def _input_i(size_bn, array, name):
    """
    _input_i(size_bn, array, name)
    
    
    Defined at array_io.fpp lines 88-95
    
    Parameters
    ----------
    size_bn : int
    array : int array
    name : str
    
    """
    _AresMainPy_pkg.f90wrap_input_i(size_bn=size_bn, array=array, name=name)

def input(*args, **kwargs):
    """
    input(*args, **kwargs)
    
    
    Defined at array_io.fpp lines 17-18
    
    Overloaded interface containing the following procedures:
      _input_r
      _input_i
    
    """
    for proc in [_input_r, _input_i]:
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
    logging.debug('unallocated array(s) detected on import of module "array_io".')

for func in _dt_array_initialisers:
    func()

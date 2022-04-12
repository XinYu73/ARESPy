"""
Module m_time_evaluate


Defined at m_time_evaluate.fpp lines 5-202

"""
from __future__ import print_function, absolute_import, division
import _AresMainPy_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("AresMainPy_pkg.time_record")
class time_record(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=time_record)
    
    
    Defined at m_time_evaluate.fpp lines 15-18
    
    """
    def __init__(self, handle=None):
        """
        self = Time_Record()
        
        
        Defined at m_time_evaluate.fpp lines 15-18
        
        
        Returns
        -------
        this : Time_Record
        	Object to be constructed
        
        
        Automatically generated constructor for time_record
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _AresMainPy_pkg.f90wrap_time_record_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Time_Record
        
        
        Defined at m_time_evaluate.fpp lines 15-18
        
        Parameters
        ----------
        this : Time_Record
        	Object to be destructed
        
        
        Automatically generated destructor for time_record
        """
        if self._alloc:
            _AresMainPy_pkg.f90wrap_time_record_finalise(this=self._handle)
    
    @property
    def str(self):
        """
        Element str ftype=character(len=100) pytype=str
        
        
        Defined at m_time_evaluate.fpp line 16
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_time_record__array__str(self._handle)
        if array_handle in self._arrays:
            str = self._arrays[array_handle]
        else:
            str = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_time_record__array__str)
            self._arrays[array_handle] = str
        return str
    
    @str.setter
    def str(self, str):
        self.str[...] = str
    
    @property
    def t1(self):
        """
        Element t1 ftype=integer(i4b) pytype=int
        
        
        Defined at m_time_evaluate.fpp line 17
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_time_record__array__t1(self._handle)
        if array_handle in self._arrays:
            t1 = self._arrays[array_handle]
        else:
            t1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_time_record__array__t1)
            self._arrays[array_handle] = t1
        return t1
    
    @t1.setter
    def t1(self, t1):
        self.t1[...] = t1
    
    @property
    def t2(self):
        """
        Element t2 ftype=integer(i4b) pytype=int
        
        
        Defined at m_time_evaluate.fpp line 17
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _AresMainPy_pkg.f90wrap_time_record__array__t2(self._handle)
        if array_handle in self._arrays:
            t2 = self._arrays[array_handle]
        else:
            t2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _AresMainPy_pkg.f90wrap_time_record__array__t2)
            self._arrays[array_handle] = t2
        return t2
    
    @t2.setter
    def t2(self, t2):
        self.t2[...] = t2
    
    @property
    def length(self):
        """
        Element length ftype=integer(i4b) pytype=int
        
        
        Defined at m_time_evaluate.fpp line 18
        
        """
        return _AresMainPy_pkg.f90wrap_time_record__get__length(self._handle)
    
    @length.setter
    def length(self, length):
        _AresMainPy_pkg.f90wrap_time_record__set__length(self._handle, length)
    
    def __str__(self):
        ret = ['<time_record>{\n']
        ret.append('    str : ')
        ret.append(repr(self.str))
        ret.append(',\n    t1 : ')
        ret.append(repr(self.t1))
        ret.append(',\n    t2 : ')
        ret.append(repr(self.t2))
        ret.append(',\n    length : ')
        ret.append(repr(self.length))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def record(str, t):
    """
    record(str, t)
    
    
    Defined at m_time_evaluate.fpp lines 25-50
    
    Parameters
    ----------
    str : str
    t : int
    
    """
    _AresMainPy_pkg.f90wrap_record(str=str, t=t)

def reset_data(myflag, str, t):
    """
    reset_data(myflag, str, t)
    
    
    Defined at m_time_evaluate.fpp lines 52-81
    
    Parameters
    ----------
    myflag : bool
    str : str
    t : int
    
    """
    _AresMainPy_pkg.f90wrap_reset_data(myflag=myflag, str=str, t=t)

def init_recrod():
    """
    init_recrod()
    
    
    Defined at m_time_evaluate.fpp lines 83-88
    
    
    """
    _AresMainPy_pkg.f90wrap_init_recrod()

def time_start(str, flag):
    """
    time_start(str, flag)
    
    
    Defined at m_time_evaluate.fpp lines 90-100
    
    Parameters
    ----------
    str : str
    flag : bool
    
    """
    _AresMainPy_pkg.f90wrap_time_start(str=str, flag=flag)

def time_end(str, flag):
    """
    time_end(str, flag)
    
    
    Defined at m_time_evaluate.fpp lines 102-112
    
    Parameters
    ----------
    str : str
    flag : bool
    
    """
    _AresMainPy_pkg.f90wrap_time_end(str=str, flag=flag)

def time_output(str, flag):
    """
    time_output(str, flag)
    
    
    Defined at m_time_evaluate.fpp lines 114-129
    
    Parameters
    ----------
    str : str
    flag : bool
    
    """
    _AresMainPy_pkg.f90wrap_time_output(str=str, flag=flag)

def memory_sum(label, memory_size):
    """
    memory_sum(label, memory_size)
    
    
    Defined at m_time_evaluate.fpp lines 131-159
    
    Parameters
    ----------
    label : str
    memory_size : float
    
    """
    _AresMainPy_pkg.f90wrap_memory_sum(label=label, memory_size=memory_size)

def memory_free(label, memory_size):
    """
    memory_free(label, memory_size)
    
    
    Defined at m_time_evaluate.fpp lines 161-183
    
    Parameters
    ----------
    label : str
    memory_size : float
    
    """
    _AresMainPy_pkg.f90wrap_memory_free(label=label, memory_size=memory_size)

def get_total_memory_consum():
    """
    Element total_memory_consum ftype=real(dp) pytype=float
    
    
    Defined at m_time_evaluate.fpp line 21
    
    """
    return _AresMainPy_pkg.f90wrap_m_time_evaluate__get__total_memory_consum()

def set_total_memory_consum(total_memory_consum):
    _AresMainPy_pkg.f90wrap_m_time_evaluate__set__total_memory_consum(total_memory_consum)

def get_filename():
    """
    Element filename ftype=character(len=20) pytype=str
    
    
    Defined at m_time_evaluate.fpp line 22
    
    """
    return _AresMainPy_pkg.f90wrap_m_time_evaluate__get__filename()

def set_filename(filename):
    _AresMainPy_pkg.f90wrap_m_time_evaluate__set__filename(filename)

def get_file_unit():
    """
    Element file_unit ftype=integer(i4b) pytype=int
    
    
    Defined at m_time_evaluate.fpp line 23
    
    """
    return _AresMainPy_pkg.f90wrap_m_time_evaluate__get__file_unit()

def set_file_unit(file_unit):
    _AresMainPy_pkg.f90wrap_m_time_evaluate__set__file_unit(file_unit)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "m_time_evaluate".')

for func in _dt_array_initialisers:
    func()
